#!/bin/bash

if [ ! -f stable_23Jun2022.zip ]; then
  wget https://github.com/lammps/lammps/archive/refs/tags/stable_23Jun2022.zip
fi

if [ ! -d lammps-stable_23Jun2022/src ]; then
  unzip stable_23Jun2022.zip lammps-stable_23Jun2022/src/*
fi

# Symlink source files.
cd lammps-stable_23Jun2022/src
ln -s ../../ba-src/angle_cgangle.cpp .
ln -s ../../ba-src/angle_cgangle.h .

make yes-molecule
make -j8 ubuntu

