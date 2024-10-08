OBJ = mod_types.o\
			arnie.o\
			io.o\
			potentials.o\
			force.o\
			main.o

## SET COMPILER
#FC = gfortran
#FC = mpifort
FC = gfortran -I/usr/local/Cellar/open-mpi/5.0.3_1/include -Wl,-flat_namespace -Wl,-commons,use_dylibs -I/usr/local/Cellar/open-mpi/5.0.3_1/lib -L/usr/local/Cellar/open-mpi/5.0.3_1/lib -lmpi_usempif08 -lmpi_usempi_ignore_tkr -lmpi_mpifh -lmpi

## EXE FILE NAME 
GOAL = TDpot

## COMPILING FLAGS
FFLAGS = -O3
#FFLAGS = -fbacktrace -Wall -Wextra -pedantic -fcheck=all -Og

##========##
## IGNORE ##
##========##
# FILIP
#LIB = -L/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk/usr/lib
#LIB_DIR = /Users/filip/bin/fortran_stdlib/lib
#MOD_DIR = /Users/filip/bin/fortran_stdlib/include/fortran_stdlib/GNU-13.2.0

# FILIP VSC
#LIB_DIR = /home/fs70998/vukovicf/.local/lib64
#MOD_DIR = /home/fs70998/vukovicf/.local/include/fortran_stdlib/GNU-10.2.0/

# LUKAS
#LIB_DIR = /home/lukas/flib/lib
#MOD_DIR = /home/lukas/flib/include/fortran_stdlib/GNU-12.3.0/

# LUKAS VSC
#LIB_DIR = /home/fs71431/essletzbich/.local/lib64
#MOD_DIR = /home/fs71431/essletzbich/.local/include/fortran_stdlib/GNU-10.2.0/
##========##

# final build
$(GOAL): $(OBJ)
	$(FC) $(OBJ) $(FFLAGS) -o $(GOAL) $(LIB)
#	$(FC) $(OBJ) $(FFLAGS) -J $(MOD_DIR) -L $(LIB_DIR) -l $(LIBS) -o $(GOAL)

%.o: %.f90
	$(FC) $(FFLAGS) -c $< -o $@

force.o: potentials.o

clean:;\
	rm -rf *.o *.mod $(GOAL)
