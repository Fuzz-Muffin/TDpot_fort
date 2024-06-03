# TDPot

## Installation

* Install Fortran (https://fortran-lang.org/learn/os_setup/install_gfortran/)

* Install mpifort (https://cloud-gc.readthedocs.io/en/latest/chapter04_developer-guide/install-basic.html#mpi-library)

* Install the Fortran Standard Library (https://github.com/fortran-lang/stdlib)

## Setup

Clone this repository

Use

```bash
make clean

```

and then

```bash
make

```

to compile the source code.

## Usage

Switch to your work folder. You need an ***input file*** named

* indat.dat

to run TDPot. Here is an example of what it have to look like, including all needed parameters:

SLG_Xe20_50keV            ! jobname
hollow-krc                ! potential
Xe                        ! ion
54 129 50 20              ! Z, mass, KE, qin
1                         ! ff
0.0 0.0 1 2               ! fwhm_qout, sigma_therm, frozen_par, alpha_max
-30 30 0.01, 0.0001 1     ! ion_zi, ion_zf, dx_step, acc, nions
0.0105 2.8 0.8 15.0       ! gam_a, gam_b, gam_c, gam_cut
1 0                       ! logmode, print_xyz_files
SLG_tdpot.xyz             ! target filename

You also need a ***target file*** in your work folder to which the ***input file*** is referencing in its last line.
In this example, the ***target file*** has to be named

* SLG_tdpot.xyz

and looks like this:

60
12.30000 12.78000 3.35000
0.61488  2.48500  1.67500  6  12
0.61488  3.90500  1.67500  6  12
...

The first line is the number of particles. The second line is the size of the simulation cell.
All other lines are representing properties of the particles. The first three columns are the x, y and z position. The fourth column is the particle type and the fifth column is the mass.

To start TDPot, use

```bash
mpirun ~/.../TDPot

```

in your work folder. The output file(s) will be written in that folder.

To use several cores, in this case four, use:

```bash
mpirun -n 4 ~/.../TDPot

```
