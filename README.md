# CGBA-potentials
A repository for bond angle potentials for publication

The repository contains 12 CG models trained with different KDE bandwidth, each model has a pair potential and a bond_angle (BA) potential, tabulated in individual files. 
The file name contains the information what potential it is and what bandwidth it was trained with, e.g., pair.table.EE.w1.2.1x means this is the pair potential for the model trained with ![(w_r^*,&space;2w_l^*,&space;w_\theta^*)](https://latex.codecogs.com/svg.image?(w_r^*,&space;2w_l^*,&space;w_\theta^*)).
Potentials from different models should not be used together.

To build lammps for the CG potentials, go to lammps_src/ and follow the guideline in lammps_src/README.md.

