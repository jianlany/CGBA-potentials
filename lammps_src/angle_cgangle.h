/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef ANGLE_CLASS

AngleStyle(cgangle,CgAngle)

#else

#ifndef LMP_ANGLE_CGANGLE_H
#define LMP_ANGLE_CGANGLE_H

#include "stdio.h"
#include "angle.h"

namespace LAMMPS_NS {

class CgAngle : public Angle {
 public:
  CgAngle(class LAMMPS *);
  virtual ~CgAngle();
  virtual void compute(int, int);
  void settings(int, char **);
  void coeff(int, char **);
  double equilibrium_angle(int);
  void write_restart(FILE *);
  void read_restart(FILE *);
  double single(int, int, int, int);
  void gather_num_bonded();

 protected:
  // Table style can be LINEAR or BICUBIC.
  int tabstyle;
  // tablength is currently ignored (could be used to resample table later).
  int tablength;
  double *theta0;

  struct Table {
    // Number of values in the potential tables.
    int ninput;
    // Number of values for bond length and bond angle.
    int nl, nq;
    double theta0;
    // Bond length and bond angle.
    double *l, *q;
    // Energy and its derivatives.
    double *e, *el, *eq, *elq;
    double dl, dq;
  };

  bool gather_bonded = true;
  int *global_num_bonded = NULL;

  int ntables;
  Table *tables;
  int *tabindex;

  void allocate();
  void null_table(Table *);
  void free_table(Table *);
  void read_table(Table *, char *, char *);
  void bcast_table(Table *);

  void param_extract(Table *, char *);
  void ufl_lookup(int, double, double ,double &, double &, double &);
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Unknown table style in angle style table

Self-explanatory.

E: Illegal number of angle table entries

There must be at least 2 table entries.

E: Invalid angle table length

Length must be 2 or greater.

E: Angle table must range from 0 to 180 degrees

Self-explanatory.

E: Cannot open file %s

The specified file cannot be opened.  Check that the path and name are
correct. If the file is a compressed file, also check that the gzip
executable can be found and run.

E: Did not find keyword in table file

Keyword used in pair_coeff command was not found in table file.

E: Invalid keyword in angle table parameters

Self-explanatory.

E: Angle table parameters did not set N

List of angle table parameters must include N setting.

*/
