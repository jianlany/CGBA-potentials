This folder contain the cpp source code and header file of the BA potential implementation.

## Installation
For quick install on Linux system just run ./get_lammps.sh, the script will download the LAMMPS, unzip and compile it if all the required libraries are ready.

The binary to be called is ```./lammps-stable_5Jun2019/src/lmp_ubuntu.```

If need to build manually, the procedures are:
1. Download LAMMPS from https://www.lammps.org/
2. Create soft links of cpp and header file in the LAMMPS `src/` directory
3. Go to the LAMMPS `src/` directory and include the Molecule package with `make yes-molecule` (and install other needed packages)
4. Build LAMMPS with `make -j 4 ubuntu`

More detailed installation guide for customized installzation and required compiler of LAMMPS can be found in LAMMPS documentation at https://docs.lammps.org/Build.html

## Usage

In the LAMMPS input script include
```
 angle_style cgangle bicubic [ba table resolution nl*nq]
 angle_coeff [angle type] [keyword]
 bond_style zero
 bond_coeff *
 ```


## Example LAMMPS script
```
 units           real
 atom_style      angle
 boundary        p p p

 # Forcefield parameters
 bond_style     zero
 angle_style    cgangle bicubic 1002001
 pair_style     table linear 16000

 read_data      PE_50chains_160beads-6_resized.lammps nocoeff

 bond_coeff     *
 angle_coeff    1 potentials/bond_angle.table.EEE.37 EEE
 pair_coeff     1 1 potentials/pair.table.EE.37 EE
 ```