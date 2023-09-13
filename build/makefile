
# VARIABLES

MPICC := mpic++
CFLAGS := -O3 -Wall -std=c++20 -fopenmp

EXE := lost-fraction.exe
OBJ := $(patsubst ../src/%.cc, %.o, $(shell find ../src -maxdepth 1 -name "*.cc"))
TEST := $(shell find ../src -maxdepth 1 -name "*.cc")

# Path to 'src' and 'gyronimo' library headers
INCLUDE := -I../src -I/home/master/boris/gyrofork/gyronimo
# Path to gyronimo compiled library file
LIBGYRONIMO := -L/home/master/boris/gyrofork/build -lgyronimo
# Path to gyronimo compiled library file
GYRONIMOPATH := -Wl,-rpath /home/master/boris/gyrofork/build
# Path to GSL and BOOST
LIBGSL := -L/usr/local/lib -lgsl -lgslcblas -lm -ltbb


#	TARGETS

#	displays the different commands
.PHONY: options
options:
	@echo "\n--- OPTIONS ---\n"
	@echo "\"make options\"   - displays the different commands"
	@echo "\"make variables\" - dumps all variables and their values"
	@echo "\"make all\"       - compiles all programs in the directory"
	@echo "\"make %.exe\"     - compiles the file into the executable"
	@echo "\"make %.o\"       - produces object file"
	@echo "\"make touch\"     - touches everything"
	@echo "\"make clean\"     - removes executables, object files and libraries"
	@echo " "

#	dumps all variables and their values
.PHONY: variables
variables:
	@echo "\n--- MPICC ---"
	@echo " $(MPICC)"
	@echo "\n--- CFLAGS ---"
	@echo " $(CFLAGS)"
	@echo "\n--- EXE ---"
	@echo " $(EXE)"
	@echo "\n--- OBJ ---"
	@echo " $(OBJ)"
	@echo "\n--- INCLUDE ---"
	@echo " $(INCLUDE)"
	@echo "\n--- LIBGYRONIMO ---"
	@echo " $(LIBGYRONIMO)"
	@echo "\n--- GYRONIMOPATH ---"
	@echo " $(GYRONIMOPATH)"
	@echo "\n--- LIBGSL ---"
	@echo " $(LIBGSL)"
	@echo "\n--- TEST ---"
	@echo " $(TEST)"

# Compiles all programs
.PHONY: all
all: $(EXE)
	@echo "Compiled everything successfully."

# Rule to compile any executable with gyronimo
.PRECIOUS: %.exe
%.exe: ../%.cc $(OBJ)
	@$(MPICC) $(CFLAGS) -o $@ $< $(OBJ) $(INCLUDE) $(LIBGYRONIMO) $(GYRONIMOPATH)  $(LIBGSL)
	@echo "Program [ $@ ] ready for execution."

# Rule to compile any object file
.PRECIOUS: %.o
%.o: ../src/%.cc
	@$(CC) $(CFLAGS) -c $^ -o $@ $(INCLUDE) $(LIBGYRONIMO)  $(LIBGSL)
	@echo "Compiled [ $@ ]."

#	touches everything
.PHONY: touch
touch:
	@touch ../src/* ../*.cc
	@echo "Everything touched"

# Delete all .exe and .o files
.PHONY: clean
clean: 
	@rm -rfv *.exe *.o