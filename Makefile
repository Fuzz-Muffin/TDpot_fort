OBJ = arnie.o\
			io.o\
			potentials.o\
			force.o\
			main.o

#FC = gfortran
FC = mpifort
GOAL = TDpot

FFLAGS = -O3
#FFLAGS = -fbacktrace -Wall -Wextra -pedantic -fcheck=all -Og

LIBS = fortran_stdlib

# VSC
#LIB_DIR = /Users/filip/bin/fortran_stdlib/lib
#MOD_DIR = /Users/filip/bin/fortran_stdlib/include/fortran_stdlib/GNU-13.2.0

# FILIP
#LIB_DIR = /home/fs70998/vukovicf/.local/lib64
#MOD_DIR = /home/fs70998/vukovicf/.local/include/fortran_stdlib/GNU-10.2.0/

# LUKAS
LIB_DIR = /home/lukas/flib/lib
MOD_DIR = /home/lukas/flib/include/fortran_stdlib/GNU-12.3.0/

# final build
$(GOAL): $(OBJ)
	$(FC) $(OBJ) $(FFLAGS) -J $(MOD_DIR) -L $(LIB_DIR) -l $(LIBS) -o $(GOAL)

%.o: %.f90
	$(FC) $(FFLAGS) -J $(MOD_DIR) -L $(LIB_DIR) -l $(LIBS) -c $< -o $@

force.o: potentials.o

clean:;\
	rm -rf *.o *.mod $(GOAL)
