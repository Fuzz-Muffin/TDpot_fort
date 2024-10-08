# TDPot

## Requirements:

* Fortran compiler (fortran 90 and above), openmpi (and the required wrapper to compile)

## Setup

Clone this repository. 
Modify the Makefile so that the your fortran (mpifort) compilers are used.
If all is in order, go ahead and run

```bash
make clean && make

```

to compile the source code. Should be quick and easy :)

## Usage

Switch to your work folder. You need an ***input file*** named

* indat.dat

to run TDPot. Here is an example of what it have to look like, including all needed parameters:

```
SLG_Xe20_50keV            ! jobname
hollow-krc                ! potential
Xe                        ! ion
54 129 50 20              ! Z, mass, KE, qin
1                         ! ff
0.0 0.0 1 2               ! fwhm_qout, sigma_therm, frozen_par, alpha_max
-30 30 0.01, 0.0001 1     ! ion_zi, ion_zf, dx_step, acc, nions
0.0105 2.8 0.8 15.0       ! gam_a, gam_b, gam_c, gam_cut
0 0                       ! logmode, print_xyz_files
SLG_tdpot.xyz             ! target filename
```

You also need a ***target file*** in your work folder to which the ***input file*** is referencing in its last line.
In this example, the ***target file*** is named

* SLG_tdpot.xyz

and is an 'xyz' file of sorts. looks like this:

```
60
12.30000 12.78000 3.35000
0.61488  2.48500  1.67500  6  12
0.61488  3.90500  1.67500  6  12
...
```

The first line is the number of particles. The second line is the size of the simulation cell (x, y and z).
All other lines are representing properties of the particles. The first three columns are the x, y and z position. The fourth column is the particle type and the fifth column is the mass.

!NOTE! this version of TDPot assumes the target to be periodic in x and y! Please be careful when providing target files. If you want to model ions incident on a non-periodic cell, some tricks will have to be used. Please email me if unsure.

To run the code, we need to used the mpirun:

```bash
mpirun -n 1 /Path/to/exe_file/.../TDPot

```

in your work folder, using the correct path to the TDPot executable that we just built. The output file(s) will be written in that folder.

To use several cores in parallel, in this case four, use:

```bash
mpirun -n 4 /Path/to/exe_file/.../TDPot

```

Happy flying :)

Filip V. and Lukas E. 
