#!/bin/bash

if [ ! -f stable_5Jun2019.zip ]; then
  wget https://github.com/lammps/lammps/archive/stable_5Jun2019.zip
fi

if [ ! -d lammps-stable_5Jun2019/src ]; then
  unzip stable_5Jun2019.zip lammps-stable_5Jun2019/src/* 
fi

# Symlink source files.
cd lammps-stable_5Jun2019/src
ln -s ../../ba-src/angle_cgangle.cpp .
ln -s ../../ba-src/angle_cgangle.h .

make yes-molecule
make -j8 ubuntu

