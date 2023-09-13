#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <chrono>
#include <numbers>
#include <random>
#include <mpi.h>

#include <gyronimo/version.hh>
#include <gyronimo/core/codata.hh>
#include <gyronimo/core/contraction.hh>
#include <gyronimo/fields/equilibrium_vmec.hh>
#include <gyronimo/interpolators/cubic_gsl.hh>

#include <gyronimo/dynamics/guiding_centre.hh>
#include <gyronimo/dynamics/odeint_adapter.hh>
#include <gyronimo/dynamics/classical_boris.hh>
#include <gyronimo/dynamics/curvilinear_boris.hh>

#include <boost/numeric/odeint/stepper/runge_kutta4.hpp>

#include <gsl/gsl_errno.h>

#include "argh.h"

void print_help() {
	std::cout << "lost-fraction, powered by ::gyronimo::v"
		<< gyronimo::version_major << "." << gyronimo::version_minor << ".\n";
	std::cout << "usage: lost-fraction.exe [options] vmec_netcdf_file\n";
	std::cout <<
		"reads a vmec output file, prints the required orbit to stdout.\n";
	std::cout << "options:\n";
	std::cout << "  -stepper=val   Name of the stepper algorithm:\n";
	std::cout << "      + guiding_centre (default)\n";
	std::cout << "      + classical_boris\n";
	std::cout << "      + curvilinear_boris\n";
	std::cout << "  -lref=val      Reference length   (in SI, default 1).\n";
	std::cout << "  -vref=val      Reference velocity (in SI, default 1).\n";
	std::cout << "  -mass=val      Particle mass   (in m_proton, default 1).\n";
	std::cout << "  -charge=val    Particle charge (in q_proton, default 1).\n";
	std::cout << "  -flux=val      Flux surface label (default random).\n";
	std::cout << "  -zeta=val      Toroidal angle (in rad, default random).\n";
	std::cout << "  -theta=val     Poloidal angle (in rad, default random).\n";
	std::cout << "  -energy=val    Energy (in eV, default 3.5e6 eV).\n";
	std::cout << "  -pitch=val     Pitch angle (between [-1,1], default random).\n";
	std::cout << "  -gyrophase=val Gyrophase angle (default random).\n";
	std::cout << "  -tfinal=val    Time limit (in Tref, default 1).\n";
	std::cout << "  -nsamples=val  Number of time samples (default 512).\n";
	std::cout << "  -ngyrons=val   Number of particles (default 1).\n";
	std::cout << "  -seed=val      Random generator seed (default 0 sets CPU time).\n";
	std::exit(0);
}

using doubles = std::vector<double>;

class particle {
public:
	virtual ~particle() {};
	virtual void do_step(double time, double dt) = 0;
	virtual gyronimo::IR3 position() = 0;
	virtual double energy_kinetic(double time) = 0;
	virtual double pitch(double time) = 0;
};

template<typename stepperclass>
class particle_fo : public particle {
public:
	particle_fo(const stepperclass *stepper, 
			typename stepperclass::state init) 
			: stepper_(stepper), s_(init) {};
	virtual ~particle_fo() override {};

	virtual void do_step(double time, double dt) override {
		s_ = stepper_->do_step(s_, time, dt);
		return;
	};
	virtual gyronimo::IR3 position() override {
		return stepper_->get_position(s_);
	};
    virtual double energy_kinetic(double time) override {
		return stepper_->energy_kinetic(s_);
	};
    virtual double pitch(double time) override {
		double Ekin = stepper_->energy_kinetic(s_);
		double Eperp = stepper_->energy_perpendicular(s_, time);
        return std::sqrt(1.0-Eperp/Ekin);
	};

private:
	const stepperclass *stepper_;
	typename stepperclass::state s_;
};

class particle_gc : public particle {
public:
	particle_gc(double Lref, double Vref, double qom, gyronimo::IR3field_c1 *Bfield, 
			gyronimo::IR3 Xgc, double mu_tilde, 
			double energy_tilde, gyronimo::guiding_centre::vpp_sign sign) 
			: gc_(Lref, Vref, qom, mu_tilde, Bfield), sys_(&gc_) {
		s_ = gc_.generate_state(Xgc, energy_tilde, sign);
	};
	virtual ~particle_gc() override {};

	void do_step(double time, double dt) {
		boost::numeric::odeint::runge_kutta4<gyronimo::guiding_centre::state> rk4;
		rk4.do_step(sys_, s_, time, dt);
		return;
	};
	virtual gyronimo::IR3 position() override {
		return gc_.get_position(s_);
	};
	virtual double energy_kinetic(double time) override {
		return gc_.energy_parallel(s_) + gc_.energy_perpendicular(s_, time);
	};
    virtual double pitch(double time) override {
        double Epar = gc_.energy_parallel(s_);
        double Eperp = gc_.energy_perpendicular(s_, time);
        return std::sqrt(1.0-Eperp/(Epar+Eperp));
	};

private:
	gyronimo::guiding_centre gc_;
	gyronimo::odeint_adapter<gyronimo::guiding_centre> sys_;
	gyronimo::guiding_centre::state s_;
};

class particle_factory {
public:
	particle_factory() {};
	virtual ~particle_factory() {};
	virtual particle* generate_particle(size_t particle_index) = 0;
};

template<typename stepperclass>
class particle_fo_factory : public particle_factory {
public:
	particle_fo_factory(stepperclass *stepper, const gyronimo::morphism *morph,
		doubles& vflux, doubles& vzeta, doubles& vtheta, 
		double energy_tilde, doubles& vpitch, doubles& vgyrophase) 
		: particle_factory(), stepper_(stepper), Bfield_(stepper->magnetic_field()), 
		morph_(morph), vflux_(vflux), vzeta_(vzeta), vtheta_(vtheta), 
		vel_tilde_(std::sqrt(energy_tilde)), vpitch_(vpitch), vgyrophase_(vgyrophase) {};
	virtual ~particle_fo_factory() override {};
	virtual particle_fo<stepperclass>* generate_particle(size_t particle_index) override {

		gyronimo::IR3 x0 = {vflux_[particle_index], vzeta_[particle_index], vtheta_[particle_index]};

		gyronimo::IR3 Bversor_con = Bfield_->contravariant_versor(x0, 0.0);
		gyronimo::IR3 Bversor = morph_->from_contravariant(Bversor_con, x0);

		double sin_phi = std::sin(vgyrophase_[particle_index]);
		double cos_phi = std::cos(vgyrophase_[particle_index]);
		double bx = Bversor[gyronimo::IR3::u];
		double by = Bversor[gyronimo::IR3::v];
		double bz = Bversor[gyronimo::IR3::w];
		double r = std::sqrt(bx*bx+by*by);
		double ca = bz, sa = r, cb = by/r, sb = bx/r;
		gyronimo::IR3 vperp_versor = {
			sin_phi*cb - cos_phi*ca*sb,
			-sin_phi*sb - cos_phi*ca*cb,
			cos_phi*sa
		};

		double vpar_tilde = vel_tilde_ * vpitch_[particle_index];
		double vperp_tilde = vel_tilde_ * std::sqrt(1-vpitch_[particle_index]*vpitch_[particle_index]);
		double sign = std::copysign(1.0, stepper_->qom());

		gyronimo::IR3 v0 = vpar_tilde*Bversor + (sign*vperp_tilde)*vperp_versor;

		typename stepperclass::state init = stepper_->generate_state(x0, v0);
		return new particle_fo<stepperclass>(stepper_, init);
	};
private:
	const stepperclass *stepper_;
	const gyronimo::IR3field *Bfield_;
	const gyronimo::morphism *morph_;
	doubles& vflux_;
	doubles& vzeta_;
	doubles& vtheta_;
	double vel_tilde_;
	doubles& vpitch_;
	doubles& vgyrophase_;
};

class particle_gc_factory : public particle_factory {
public:
	particle_gc_factory(double Lref, double Vref, double qom, 
		gyronimo::IR3field_c1 *Bfield, doubles& vflux, doubles& vzeta, 
			doubles& vtheta, double energy_tilde, doubles& vpitch) 
		: particle_factory(), Lref_(Lref), Vref_(Vref), qom_(qom), Bfield_(Bfield), 
		vflux_(vflux), vzeta_(vzeta), vtheta_(vtheta), energy_tilde_(energy_tilde), vpitch_(vpitch) {};
	virtual ~particle_gc_factory() override {};
	virtual particle_gc* generate_particle(size_t particle_index) override {
		gyronimo::IR3 x0 = {vflux_[particle_index], vzeta_[particle_index], vtheta_[particle_index]};
		double lambda = (1 - vpitch_[particle_index]*vpitch_[particle_index]) / Bfield_->magnitude(x0, 0.0);
		double mu_tilde = lambda * energy_tilde_;
		gyronimo::guiding_centre::vpp_sign sign = vpitch_[particle_index] > 0 ?
			gyronimo::guiding_centre::plus : 
			gyronimo::guiding_centre::minus;
		return new particle_gc(Lref_, Vref_, qom_, Bfield_, x0, mu_tilde, energy_tilde_, sign);
	};
private:
	double Lref_, Vref_, qom_;
	gyronimo::IR3field_c1 *Bfield_;
	doubles& vflux_;
	doubles& vzeta_;
	doubles& vtheta_;
	double energy_tilde_;
	doubles& vpitch_;
};

void new_handler(const char *reason, const char *file, int line, int gsl_errno) {
	if(gsl_errno == GSL_EDOM) {
		throw std::domain_error("Outside domain of interpolation.");
	}
	std::cout << "gsl: " << file << ":" << line << ": ERROR: " << reason << "\n";
	std::cout << "Custom GSL error handler invoked.\nAborted\n";
	std::exit(1);
	return;
}



int main(int argc, char *argv[]) {

	auto command_line = argh::parser(argv);
	if (command_line[{"h", "help"}]) print_help();
	if (!command_line(1)) {  // the 1st non-option argument is the mapping file.
		std::cerr << "lost-fraction: no vmec mapping file provided; -h for help.\n";
		std::exit(0);
	}

	int rank, size;
	/* Initialize MPI */
	MPI_Init(&argc, &argv);
	/* Get total number of procs. */
	MPI_Comm_size( MPI_COMM_WORLD, &size );
	/* Get my rank */
	MPI_Comm_rank( MPI_COMM_WORLD, &rank );

	gsl_set_error_handler_off();
	gsl_set_error_handler(&new_handler);

	/* Initialize helena equilibrium */
	gyronimo::parser_vmec vmap(command_line[1]);
	gyronimo::cubic_gsl_factory ifactory;
	gyronimo::morphism_vmec m(&vmap, &ifactory);
	gyronimo::metric_vmec g(&m, &ifactory);
	gyronimo::equilibrium_vmec veq(&g, &ifactory);
    double infp = 1.0 / vmap.nfp();

	std::function<bool(gyronimo::IR3)> escape_condition = 
		[](gyronimo::IR3 pos) {
			if(pos[gyronimo::IR3::u] > 0.99) return true;
			else return false;
		};

	/* Reads parameters from the command line: */
	size_t ngyrons_total;  command_line("ngyrons",  1) >> ngyrons_total;
	size_t ngyrons = std::ceil(ngyrons_total/((double)size));

	size_t seed; command_line("seed", 0) >> seed;
	if(seed == 0) {
		auto time = std::chrono::system_clock::now();
		auto duration = std::chrono::duration_cast<std::chrono::hours>(time.time_since_epoch());
		seed = duration.count() + rank;
	} else seed += rank;
	std::mt19937 rng(seed);
	std::uniform_real_distribution uniform_01(0.0, 1.0);
	std::uniform_real_distribution uniform_angle(0.0, 2*std::numbers::pi);

	std::string steppername; command_line("stepper", "guiding_centre") >> steppername;
	double Lref;   command_line("lref",   1.0) >> Lref;   // SI units.
	double Vref;   command_line("vref",   1.0) >> Vref;   // SI units.
	double mass;   command_line("mass",   1.0) >> mass;   // m_proton units.
	double charge; command_line("charge", 1.0) >> charge; // q_electron units.
	double energy; command_line("energy", 3.5e6) >> energy; // energy in eV.

	doubles vflux(ngyrons);
	doubles vzeta(ngyrons);
	doubles vtheta(ngyrons);
	doubles vpitch(ngyrons);
	doubles vgyrophase(ngyrons);

	double flux; std::string flux_str; command_line("flux", "random") >> flux_str;
	if(flux_str != "random") {
		command_line("flux", -1.0) >> flux; // normalized units.
		if(flux >= 0.0 && flux <= 1.0) {
			for(size_t i = 0; i < ngyrons; ++i) vflux[i] = flux;
		} else {
			std::cerr << "Invalid value for flux: " << flux << std::endl;
			MPI_Finalize();
			return 0;
		}
	} else {
		for(size_t i = 0; i < ngyrons; ++i) vflux[i] = uniform_01(rng);
	}

	double zeta; std::string zeta_str; command_line("zeta", "random") >> zeta_str;
	if(zeta_str != "random") {
		command_line("zeta", 0.0) >> zeta; // angle units.
		for(size_t i = 0; i < ngyrons; ++i) vzeta[i] = zeta;
	} else for(size_t i = 0; i < ngyrons; ++i) vzeta[i] = uniform_angle(rng) * infp;

	double theta; std::string theta_str; command_line("theta", "random") >> theta_str;
	if(theta_str != "random") {
		command_line("theta", 0.0) >> theta; // angle units.
		for(size_t i = 0; i < ngyrons; ++i) vtheta[i] = theta;
	} else for(size_t i = 0; i < ngyrons; ++i) vtheta[i] = uniform_angle(rng);

	double pitch; std::string pitch_str; command_line("pitch", "random") >> pitch_str;
	double lambda; std::string lambda_str; command_line("lambda", "random") >> lambda_str;
	if(pitch_str != "random") {
		command_line("pitch", 0.5) >> pitch;
		for(size_t i = 0; i < ngyrons; ++i) vpitch[i] = pitch;
	} else if (lambda_str != "random") {
		command_line("lambda", 0.5) >> lambda;
		double vpp_sign = std::copysign(1.0, lambda);
		for(size_t i = 0; i < ngyrons; ++i) vpitch[i] = vpp_sign*std::sqrt(1 - 
			std::abs(lambda)*veq.magnitude({vflux[i], vzeta[i], vtheta[i]}, 0.0));
	} else {
		for(size_t i = 0; i < ngyrons; ++i) // vpitch[i] = std::cos(0.5 * uniform_angle(rng));
            vpitch[i] = 2.0 * uniform_01(rng)-1.0;
	}

	double gyrophase; std::string phase_str; command_line("gyrophase", "random") >> phase_str;
	if(phase_str != "random") {
		command_line("gyrophase", 0.0) >> gyrophase; // angle units.
		for(size_t i = 0; i < ngyrons; ++i) vgyrophase[i] = gyrophase;
	} else for(size_t i = 0; i < ngyrons; ++i) vgyrophase[i] = uniform_angle(rng);

	double Tfinal; command_line("tfinal", 1.0) >> Tfinal; // in Tref units
	size_t nsamples; command_line("nsamples", 512) >> nsamples;
	double dt = Tfinal / nsamples;

	/* Computes normalisation constants: */
	double qom = charge / mass;
	double Tref = Lref / Vref;
	double Uref = 0.5*gyronimo::codata::m_proton*mass*Vref*Vref/gyronimo::codata::e; // Uref in eV
	double energy_tilde = energy / Uref;

	std::unique_ptr<particle_factory> pfactory = nullptr;
	std::unique_ptr<gyronimo::classical_boris>  boris_class = nullptr;
	std::unique_ptr<gyronimo::curvilinear_boris> boris_curv = nullptr;

	if(steppername == "guiding_centre") {

		pfactory.reset(new particle_gc_factory(Lref, Vref, qom, &veq, vflux, vzeta, vtheta, energy_tilde, vpitch));

	} else if(steppername == "classical_boris") {

		boris_class.reset(new gyronimo::classical_boris(Lref, Vref, qom, &veq, nullptr));
		pfactory.reset(new particle_fo_factory<gyronimo::classical_boris>(boris_class.get(), 
			&m, vflux, vzeta, vtheta, energy_tilde, vpitch, vgyrophase));

	} else if(steppername == "curvilinear_boris") {

		boris_curv.reset(new gyronimo::curvilinear_boris(Lref, Vref, qom, &veq, nullptr));
		pfactory.reset(new particle_fo_factory<gyronimo::curvilinear_boris>(boris_curv.get(), 
			&m, vflux, vzeta, vtheta, energy_tilde, vpitch, vgyrophase));

	} else {
		std::cerr << "lost-fraction: stepper name not recognized; -h for help.\n";
		MPI_Finalize();
		std::exit(1);
	}

	/* Write results to output */
	std::string outfile; command_line("out", "lfrac") >> outfile;
	std::string filename = outfile + std::to_string(rank) + ".dat";
	std::ofstream file(filename);
	if(!file.is_open()) {
		std::cerr << "Could not open file '" << filename << "'" << std::endl;
		// MPI_Finalize();
		return 0;
	}

	/* Print output header */
	file << "# lost-fraction, powered by ::gyronimo:: v"
		<< gyronimo::version_major << "." << gyronimo::version_minor << ".\n";
	file << "# args: ";
	for(int i = 1; i < argc; i++) file << argv[i] << " ";
	file << std::endl;
	file <<"# l_ref = " << Lref << " [m];";
	file << " v_ref = " << Vref << " [m/s];";
	file << " t_ref = " << Tref << " [s];";
	file << " u_ref = " << Uref << " [eV];\n";
	file << "# vars: t_escape s_i zeta_i theta_i energy_i pitch_i Bmag_i jac_i s_f zeta_f theta_f energy_f pitch_f Bmag_f \n";
    file.precision(16);
    file.setf(std::ios::scientific);

	/* LAUNCH PARTICLES */
    for(size_t i = 0; i < ngyrons; ++i) {
        std::unique_ptr<particle> p = nullptr;
        p.reset(pfactory->generate_particle(i));
        gyronimo::IR3 xcurrent = p->position();
        double escape_time = -1.0;
        double time = 0;
        for(size_t j = 0; j <= nsamples; ++j) {
            time = j * dt;
            xcurrent = p->position();
            if(escape_condition(xcurrent)) {
                escape_time = time;
                break;
            }
            try { p->do_step(time, dt); }
            catch(std::domain_error &e) {
                escape_time = time;
                break; 
            }
        }

        // print to file
        double Bmag_i = veq.magnitude({vflux[i], vzeta[i], vtheta[i]}, 0.0);
        double Bmag_f = veq.magnitude(xcurrent, time);
        file << escape_time << ' '
            << vflux[i] << ' ' << vzeta[i] << ' ' << vtheta[i] << ' '
            << energy_tilde << ' ' << vpitch[i] << ' ' << Bmag_i << ' '
            << m.jacobian({vflux[i], vzeta[i], vtheta[i]}) << ' '
            << xcurrent[gyronimo::IR3::u] << ' '
            << xcurrent[gyronimo::IR3::v] << ' '
            << xcurrent[gyronimo::IR3::w] << ' '
            << p->energy_kinetic(time) << ' '
            << p->pitch(time) << ' ' << Bmag_f << '\n';
    }

	file.close();
	
	/* Close MPI */
	MPI_Finalize();

	return 0;
}