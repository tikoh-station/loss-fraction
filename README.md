# ::lost-fraction::
============

*- Compute the lost-fraction of particles in stellarators -*

*App based on ::gyronimo:: object-oriented library on GitHub*

Licensing and terms of use:
---------------------------

`lost-fraction` is free software: you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by the
Free Software Foundation, either version 3 of the License, or (at your
option) any later version.

`lost-fraction` is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
Public License for more details.

You should have received a copy of the GNU General Public License along
with this program. If not, see <https://www.gnu.org/licenses/>.

Compiling and installing, a minimal and very optimistic howto:
--------------------------------------------------------------

The app requires an mpi compiler supporting the **c++20**
standard (e.g., gcc-10.1.0 or later). 

#### Basic build procedure:

1. Move to the build folder,
2. Update makefile with paths to the required libraries,
3. Run `make all` to compile the executable;

External dependencies and other nightmares:
-------------------------------------------

+ [gyronimo](https://github.com/prodrigs/gyronimo.git), an 
  object-oriented library for gyromotion applications in plasma physics;
+ [open-mpi](https://www.open-mpi.org/faq/?category=mpi-apps), a 
  high-performance message passing library;
+ [boost](https://www.boost.org), a broad spectrum c++ library;
+ [GSL](https://www.gnu.org/software/gsl), the GNU scientific library;
+ [netcdf-cxx4](https://github.com/Unidata/netcdf-cxx4), a c++ extension
  to NetCDF-4;