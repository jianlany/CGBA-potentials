/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Vipin Agrawal (vipin.agrawal@asu.edu)
                        Jay Oswald    (joswald1@asu.edu)
                        Jianlan Ye    (jianlany@asu.edu)
------------------------------------------------------------------------- */
#include "math.h"
#include "stdlib.h"
#include "string.h"
#include "angle_cgangle.h"
#include "atom.h"
#include "neighbor.h"
#include "domain.h"
#include "comm.h"
#include "force.h"
#include "math_const.h"
#include "memory.h"
#include "error.h"
#include <vector>
#include <algorithm>
#include "table_file_reader.h"

using namespace LAMMPS_NS;
using namespace MathConst;

enum { LINEAR, BICUBIC };

#define MAXLINE 1024
#define SMALL 0.001
// uncomment this line to print debug info.
// #define OUTPUT_DEBUG_INFO

/* ----------------------------------------------------------------------
   a simple vector that that contains 4 elements.
------------------------------------------------------------------------- */
class Vec4 {
public:
    Vec4(double a=0.0, double b=0.0, double c=0.0, double d=0.0) 
    : _v{a, b, c, d} {}
    double& operator()(int i) { return _v[i]; }
    double  operator()(int i) const { return _v[i]; }
private:
    double _v[4];
};

/* ----------------------------------------------------------------------
   a simple matrix that that contains 4x4 elements.
------------------------------------------------------------------------- */
class Matrix4 {
public:
    double& operator()(int i, int j) { return _v[4*i+j]; }
    double  operator()(int i, int j) const { return _v[4*i+j]; }
private:
    double _v[16];
};
 
Vec4 operator*(const Matrix4 &A, const Vec4 &x) {
    Vec4 y;
    for (int i=0; i<4; ++i) {
        y(i) = A(i,0)*x(0) + A(i,1)*x(1) + A(i,2)*x(2) + A(i,3)*x(3);
    }
    return y;
}
 
double dot(const Vec4 &x, const Vec4 &y) {
    return x(0)*y(0) + x(1)*y(1) + x(2)*y(2) + x(3)*y(3);
}
 
Matrix4 operator*(const Matrix4 &A, const Matrix4 &B) {
    Matrix4 C;
    for (int i=0; i<4; ++i) {
        for (int j=0; j<4; ++j) {
            C(i,j) = A(i,0)*B(0,j) + A(i,1)*B(1,j)
                   + A(i,2)*B(2,j) + A(i,3)*B(3,j);
        }
    }
    return C;
}
std::ostream& operator<<(std::ostream &o, const Vec4 &A) {
    o << A(0) << " " << A(1) << " " << A(2) << " " << A(3) << "\n";
    return o;
}

std::ostream& operator<<(std::ostream &o, const Matrix4 &A) {
    o << A(0,0) << " " << A(0,1) << " " << A(0,2) << " " << A(0,3) << "\n";
    o << A(1,0) << " " << A(1,1) << " " << A(1,2) << " " << A(1,3) << "\n";
    o << A(2,0) << " " << A(2,1) << " " << A(2,2) << " " << A(2,3) << "\n";
    o << A(3,0) << " " << A(3,1) << " " << A(3,2) << " " << A(3,3) << "\n";
    return o;
}

/* ---------------------------------------------------------------------- */

CgAngle::CgAngle(LAMMPS *lmp) : Angle(lmp)
{
  writedata = 0;
  ntables = 0;
  tables = NULL;
}

/* ---------------------------------------------------------------------- */

CgAngle::~CgAngle()
{
  for (int m = 0; m < ntables; m++) free_table(&tables[m]);
  memory->sfree(tables);
  memory->destroy(global_num_bonded);
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(theta0);
    memory->destroy(tabindex);
  }
}

/* ---------------------------------------------------------------------- */

// Compute hash value for a pair, used for building map with pair as key.
struct pair_hash
{
    template <class T1, class T2>
    std::size_t operator() (const std::pair<T1, T2> &pair) const
    {
        return std::hash<T1>()(pair.first) ^ std::hash<T2>()(pair.second);
    }
};
// Gather the number of atom bonded to an atom from all the processors
void CgAngle::gather_num_bonded(){
  bigint nlocal = atom->nlocal;
  int nbondlist = neighbor->nbondlist;
  int **bondlist = neighbor->bondlist;
  bigint natoms = 0;
  tagint *global_id = atom->tag;
  bigint local_atom_tags[nlocal];
  int newton_bond = force->newton_bond;
  int me, nproc;

  if (!bondlist) error->all(FLERR,"Bond topology list not built, create a bond style, at least use 'bond_style zero'.");

  MPI_Comm_rank(world, &me);
  MPI_Comm_size(world, &nproc);

  MPI_Allreduce(&nlocal,&natoms,1,MPI_LMP_BIGINT,MPI_SUM,world);

  std::vector<bigint> atom_appear_v;
  for (int i = 0; i < nbondlist; i++){
    int i1 = bondlist[i][0];
    int i2 = bondlist[i][1];
    bigint atom1 = global_id[i1];
    bigint atom2 = global_id[i2];
    // If newton_bond is off, then don't keep ghost atoms
    if (newton_bond || i1 < nlocal) atom_appear_v.push_back(atom1);
    if (newton_bond || i2 < nlocal) atom_appear_v.push_back(atom2);
  }

  int num_atom_appear = static_cast<int>(atom_appear_v.size());
  bigint atom_appear[num_atom_appear] = {};
  // Transform std::vector to C array for MPI communication
  std::copy(atom_appear_v.begin(), atom_appear_v.end(), atom_appear);

  // Need all the processors finished casting 
  MPI_Barrier(world);

  // Calculate the total number of entries in atom_appear and broadcast to all processors
  // for future memory allocation
  bigint total_atom_appear = 0;
  MPI_Allreduce(&num_atom_appear, &total_atom_appear, 1, MPI_INT, MPI_SUM, world);

  // Gather the size of atom_appear to root.
  int num_atom_appear_g[nproc] = {};
  MPI_Gather(&num_atom_appear, 1, MPI_INT,
             num_atom_appear_g, 1, MPI_INT, 0, world);

  int recvdispls[nproc] = {};
  for (int i = 0; i < nproc; ++i){
    if (i < nproc-1) recvdispls[i+1] = recvdispls[i] + num_atom_appear_g[i];
  }
  // Gather atom_appear array to root
  bigint atom_appear_global[total_atom_appear] = {};
  MPI_Gatherv(&atom_appear[0], num_atom_appear, MPI_LMP_BIGINT,
              &atom_appear_global[0], num_atom_appear_g, recvdispls, MPI_LMP_BIGINT, 0, world);

  memory->create(global_num_bonded, total_atom_appear, "angle:global_num_bonded");
  // Initializing the allocated memory, just to be safe
  for (bigint i = 0; i < total_atom_appear; ++i){
    global_num_bonded[i] = 0;
  }


  if (me == 0){
    for(bigint i = 0; i<total_atom_appear; ++i){
      global_num_bonded[atom_appear_global[i]]++;
    }
  }
  // Broadcast the global bonded atoms to every processor
  MPI_Bcast(global_num_bonded, total_atom_appear, MPI_INT, 0, world);
  gather_bonded = false;

/* ###############################################################################
Following is a different implementation that only works when newton_bond is off.
Keep it here for future reference.
###############################################################################*/
/*
  // Gather number of atoms of each processor
  int local_natoms[nproc] = {};
  MPI_Gather(&nlocal, 1, MPI_LMP_BIGINT, local_natoms, 1, MPI_LMP_BIGINT, 0, world);
  // Calculate total number of atoms.
  for (int i = 0; i < nproc; ++i){
    natoms += local_natoms[i];
  }
  // Keep track of the global atom id and number of bonded atom
  int num_bonded[nlocal] = {};
  for (bigint i = 0; i < nlocal; ++i){
    local_atom_tags[i] = global_id[i];
    num_bonded[i] = atom->num_bond[i];
    // std::cout << atom->num_bond[i] << '\n';
  }
 //  for(bigint i = 0; i < nlocal; ++i){
 //    std::cout<< "local atom id: " << i << " number of atom bonded: " << num_bonded[i] << '\n';
 // }
  // Calculate memory displacement in receive buffer when mpi_gather
  int recv_displs[nproc] = {};
  for (int i = 0; i < nproc; ++i){
    if (i < nproc-1) recv_displs[i+1] = recv_displs[i] + local_natoms[i];
  }
  // Gather the global atom id from all processor
  bigint atom_tags[natoms] = {};
  int mpisize=0;
  MPI_Type_size(MPI_LMP_BIGINT, &mpisize);
  MPI_Gatherv(&local_atom_tags[0], nlocal, MPI_LMP_BIGINT,
              &atom_tags[0], local_natoms, recv_displs, MPI_LMP_BIGINT,
              0, world);
  // Gather number of bonded atoms from all processor
  int recv_buf[natoms] = {};
  MPI_Gatherv(&num_bonded[0], nlocal, MPI_INT,
              &recv_buf[0], local_natoms, recv_displs, MPI_LMP_BIGINT,
              0, world);
  std::cout << "Allocating memory............";
  memory->create(global_num_bonded, natoms, "angle:global_num_bonded");
  std::cout << "Done\n";
  if (me == 0){
    for(bigint i = 0; i<natoms; ++i){
      global_num_bonded[atom_tags[i]] = recv_buf[i];
    }
  }
*/
}

void CgAngle::compute(int eflag, int vflag)
{
  int i1,i2,i3,n,type;
  double eangle,f1[3],f3[3];
  double delx1,dely1,delz1,delx2,dely2,delz2;
  double fl, ft;
  double rsq1,rsq2,r1,r2,c,s,a,a11,a12,a22;
  double theta,u,qdu; //mdu: minus du, -du/dx=f
  double u1, fl1, ft1, u2, fl2, ft2;
  double vx11,vx12,vy11,vy12,vz11,vz12,vx21,vx22,vy21,vy22,vz21,vz22;
  double b1,b2,aa1,aa2,aa11,aa12,aa21,aa22;

  eangle = 0.0;
  if (eflag || vflag) ev_setup(eflag,vflag);
  else evflag = 0;

  double **x = atom->x;
  double **f = atom->f;
  int **anglelist = neighbor->anglelist;
  int nanglelist = neighbor->nanglelist;
  bigint nlocal = atom->nlocal;
  int newton_bond = force->newton_bond;
  tagint *global_id = atom->tag;
  ///////////////////////
  if (gather_bonded) gather_num_bonded();
  for (n = 0; n < nanglelist; n++) {
    i1 = anglelist[n][0];
    i2 = anglelist[n][1];
    i3 = anglelist[n][2];

    type = anglelist[n][3];

    auto atom1 = global_id[i1];
    auto atom2 = global_id[i2];
    auto atom3 = global_id[i3];

    auto bond1_rep = global_num_bonded[atom1] + global_num_bonded[atom2] - 2;
    auto bond2_rep = global_num_bonded[atom3] + global_num_bonded[atom2] - 2;

    // 1st bond

    delx1 = x[i1][0] - x[i2][0];
    dely1 = x[i1][1] - x[i2][1];
    delz1 = x[i1][2] - x[i2][2];

    rsq1 = delx1*delx1 + dely1*dely1 + delz1*delz1;
    r1 = sqrt(rsq1);

    // 2nd bond
    delx2 = x[i3][0] - x[i2][0];
    dely2 = x[i3][1] - x[i2][1];
    delz2 = x[i3][2] - x[i2][2];

    rsq2 = delx2*delx2 + dely2*dely2 + delz2*delz2;
    r2 = sqrt(rsq2);

    // angle (cos and sin)
    c = delx1*delx2 + dely1*dely2 + delz1*delz2;
    c /= r1*r2;

    if (c > 1.0) c = 1.0;
    if (c < -1.0) c = -1.0;

    s = sqrt(1.0 - c*c);
    if (s < SMALL) s = SMALL;
    s = 1.0/s;

    // tabulated force & energy
    theta = acos(c);

    // force & energy for bond-angle term

    ufl_lookup(type,r1,theta,u1,fl1,ft1);
    ufl_lookup(type,r2,theta,u2,fl2,ft2);
    /*
    ########################################################################
    #                                                                      #
    #           fl AND ft ARE ENERGY DERIVATIVES, NOT FORCES.              #
    #                                                                      #
    ########################################################################
    */

    aa1 = s * ft1;
    aa2 = s * ft2;

    aa11 = aa1 * c / rsq1;
    aa12 = -aa1 / (r1 * r2);
    aa21 = aa2 * c / rsq1;
    aa22 = -aa2 / (r1 * r2);

    vx11 = (aa11 * delx1) + (aa12 * delx2);
    vx12 = (aa21 * delx1) + (aa22 * delx2);
    vy11 = (aa11 * dely1) + (aa12 * dely2);
    vy12 = (aa21 * dely1) + (aa22 * dely2);
    vz11 = (aa11 * delz1) + (aa12 * delz2);
    vz12 = (aa21 * delz1) + (aa22 * delz2);

    aa11 = aa1 * c / rsq2;
    aa21 = aa2 * c / rsq2;

    vx21 = (aa11 * delx2) + (aa12 * delx1);
    vx22 = (aa21 * delx2) + (aa22 * delx1);
    vy21 = (aa11 * dely2) + (aa12 * dely1);
    vy22 = (aa21 * dely2) + (aa22 * dely1);
    vz21 = (aa11 * delz2) + (aa12 * delz1);
    vz22 = (aa21 * delz2) + (aa22 * delz1);
    b1 = fl1 / r1;
    b2 = fl2 / r2;
    // Multiplier 0.5 is because the angle contribution of force is counted twice.
    // Divided by bond_count_global is because a bond is counted multiple times 
    // when it is shared by multiple angle.
    f1[0] = -((vx11 + vx12) + b1*delx1)/bond1_rep;
    f1[1] = -((vy11 + vy12) + b1*dely1)/bond1_rep;
    f1[2] = -((vz11 + vz12) + b1*delz1)/bond1_rep;

    f3[0] = -((vx21 + vx22) + b2*delx2)/bond2_rep;
    f3[1] = -((vy21 + vy22) + b2*dely2)/bond2_rep;
    f3[2] = -((vz21 + vz22) + b2*delz2)/bond2_rep;

    if (eflag) eangle = u1/bond1_rep+u2/bond2_rep;

    // apply force to each of 3 atoms

    if (newton_bond || i1 < nlocal) {
      f[i1][0] += f1[0];
      f[i1][1] += f1[1];
      f[i1][2] += f1[2];
    }

    if (newton_bond || i2 < nlocal) {
      f[i2][0] -= f1[0] + f3[0];
      f[i2][1] -= f1[1] + f3[1];
      f[i2][2] -= f1[2] + f3[2];
    }

    if (newton_bond || i3 < nlocal) {
      f[i3][0] += f3[0];
      f[i3][1] += f3[1];
      f[i3][2] += f3[2];
    }

    if (evflag) ev_tally(i1,i2,i3,nlocal,newton_bond,eangle,f1,f3,
                         delx1,dely1,delz1,delx2,dely2,delz2);
  }
}

/* ---------------------------------------------------------------------- */

void CgAngle::allocate()
{
  allocated = 1;
  int n = atom->nangletypes;

  memory->create(theta0,n+1,"angle:theta0");
  memory->create(tabindex,n+1,"angle:tabindex");

  memory->create(setflag,n+1,"angle:setflag");
  for (int i = 1; i <= n; i++) setflag[i] = 0;
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void CgAngle::settings(int narg, char **arg)
{
  if (narg != 2) error->all(FLERR,"Illegal angle_style command");

  if (strcmp(arg[0],"linear") == 0) {
    tabstyle = LINEAR;
  }
  else if (strcmp(arg[0], "bicubic") == 0) {
    tabstyle = BICUBIC;
  }
  else error->all(FLERR,"Unknown table style in angle style table");
  tablength =utils::inumeric(FLERR,arg[1],false,lmp);
  if (tablength < 4) {
    error->all(FLERR,"Illegal number of bondangle table entries");
  }
  // delete old tables, since cannot just change settings
  for (int m = 0; m < ntables; m++) {
    free_table(&tables[m]);
  }
  memory->sfree(tables);
  if (allocated) {
     memory->destroy(setflag);
     memory->destroy(tabindex);
  }
  allocated = 0;
  ntables = 0;
  tables = NULL;
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void CgAngle::coeff(int narg, char **arg)
{
  if (narg != 3) error->all(FLERR,"Illegal angle_coeff command");
  if (!allocated) allocate();

  tables = (Table *)
    memory->srealloc(tables,(ntables+1)*sizeof(Table),"angle:tables");
  Table *tb = &tables[ntables];
  null_table(tb);

  int me;
  MPI_Comm_rank(world,&me);
  if (me == 0) {
    read_table(tb,arg[1],arg[2]);
  }
  bcast_table(tb);

  // error check on table parameters
  if (tb->ninput <= 1) error->one(FLERR,"Invalid angle table length");

  tb->dq = (tb->q[tb->ninput-1] - tb->q[0]) / (tb->nq-1);
  tb->dl = (tb->l[tb->ninput-1] - tb->l[0]) / (tb->nl-1);

  int ilo,ihi;
  utils::bounds(FLERR, arg[0], 1, atom->nangletypes ,ilo,ihi, error);
  // store ptr to table in tabindex
  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    tabindex[i] = ntables;
    setflag[i] = 1;
    theta0[i] = tb->theta0;
    count++;
  }
  ntables++;

  if (count == 0) error->all(FLERR,"Illegal angle_coeff command");
}

/* ----------------------------------------------------------------------
   return an equilibrium angle
   should not be used, since don't know minimum of tabulated function
------------------------------------------------------------------------- */

double CgAngle::equilibrium_angle(int i)
{
  return theta0[i];
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
 ------------------------------------------------------------------------- */

void CgAngle::write_restart(FILE *fp)
{
  fwrite(&tabstyle,sizeof(int),1,fp);
  fwrite(&tablength,sizeof(int),1,fp);
}

/* ----------------------------------------------------------------------
    proc 0 reads from restart file, bcasts
 ------------------------------------------------------------------------- */

void CgAngle::read_restart(FILE *fp)
{
  if (comm->me == 0) {
    fread(&tabstyle,sizeof(int),1,fp);
    fread(&tablength,sizeof(int),1,fp);
  }
  MPI_Bcast(&tabstyle,1,MPI_INT,0,world);
  MPI_Bcast(&tablength,1,MPI_INT,0,world);

  allocate();
}

/* ---------------------------------------------------------------------- */

double CgAngle::single(int type, int i1, int i2, int i3)
{
  double **x = atom->x;
  double delx1 = x[i1][0] - x[i2][0];
  double dely1 = x[i1][1] - x[i2][1];
  double delz1 = x[i1][2] - x[i2][2];
  domain->minimum_image(delx1,dely1,delz1);
  double r1 = sqrt(delx1*delx1 + dely1*dely1 + delz1*delz1);

  double delx2 = x[i3][0] - x[i2][0];
  double dely2 = x[i3][1] - x[i2][1];
  double delz2 = x[i3][2] - x[i2][2];
  domain->minimum_image(delx2,dely2,delz2);
  double r2 = sqrt(delx2*delx2 + dely2*dely2 + delz2*delz2);
  double c = delx1*delx2 + dely1*dely2 + delz1*delz2;
  c /= r1*r2;
  if (c > 1.0) c = 1.0;
  if (c < -1.0) c = -1.0;

  double theta = acos(c);
  double u1=0.0;
  double u2=0.0;
  double fl1,ft1,fl2,ft2;
  ufl_lookup(type,r1,theta,u1,fl1,ft1);
  ufl_lookup(type,r2,theta,u2,fl2,ft2);

  return u1 + u2;
}

/* ---------------------------------------------------------------------- */

void CgAngle::null_table(Table *tb)
{
  tb->l = tb->q = NULL;
  tb->e = tb->el = tb->eq = tb->elq = NULL;
}

/* ---------------------------------------------------------------------- */

void CgAngle::free_table(Table *tb)
{
  memory->destroy(tb->l);
  memory->destroy(tb->q);
  memory->destroy(tb->e);
  memory->destroy(tb->el);
  memory->destroy(tb->eq);
  memory->destroy(tb->elq);
}

/* ----------------------------------------------------------------------
   read table file, only called by proc 0
------------------------------------------------------------------------- */

void CgAngle::read_table(Table *tb, char *file, char *keyword)
{
  // open file
  TableFileReader reader(lmp, file, "angle");

  char* line = reader.find_section_start(keyword);

  if (!line) {
    error->one(FLERR,"Did not find keyword in table file");
  }

  // read args on 2nd line of section
  // allocate table arrays for file values
  line = reader.next_line();
  param_extract(tb,line);

  memory->create(tb->l, tb->ninput, "angle:l");
  memory->create(tb->q, tb->ninput, "angle:q");
  memory->create(tb->e, tb->ninput, "angle:e");
  memory->create(tb->el, tb->ninput, "angle:el");
  memory->create(tb->eq, tb->ninput, "angle:eq");
  memory->create(tb->elq, tb->ninput, "angle:elq");
  
  // read a,e,f table values from file
  int itmp;
  reader.skip_line();
  for (int i = 0; i < tb->ninput; i++) {
    line = reader.next_line();
    // TODO - bilinear style could exclude elq
    sscanf(line,"%d %lg %lg %lg %lg %lg %lg",
           &itmp, &tb->l[i], &tb->q[i], &tb->e[i], 
           &tb->el[i], &tb->eq[i], &tb->elq[i]);
    tb->q[i] *= MY_PI/180.0;
    tb->eq[i] /= MY_PI/180.0;
    tb->elq[i] /= MY_PI/180.0;
  }
}

/* ----------------------------------------------------------------------
   extract attributes from parameter line in table section
   format of line: NL NQ [nl] [nq]
   where: nl is the number of bond lengths and 
          nq is the number of bond angles
------------------------------------------------------------------------- */

void CgAngle::param_extract(Table *tb, char *line)
{
  tb->theta0 = 180.0;
  tb->ninput = tb->nl = tb->nq = 0;
  char *word = strtok(line," \t\n\r\f");
  while (word) {
    if (strcmp(word,"NL") == 0) {
      word = strtok(NULL," \t\n\r\f");  // SKIP NQ
      word = strtok(NULL," \t\n\r\f");
      tb->nl = atoi(word);
      word = strtok(NULL," \t\n\r\f");
      tb->nq = atoi(word);
      // Calculate number of lines from product.
      tb->ninput = tb->nl * tb->nq;
    }  else {
      error->one(FLERR,"Invalid keyword in angle table parameters");
    }
    word = strtok(NULL," \t\n\r\f");
  }
  if (tb->ninput == 0) error->one(FLERR,"Angle table parameters did not set N");
}

/* ----------------------------------------------------------------------
   broadcast read-in table info from proc 0 to other procs
   this function communicates these values in Table:
     ninput,afile,efile,ffile,theta0
------------------------------------------------------------------------- */

void CgAngle::bcast_table(Table *tb)
{
  MPI_Bcast(&tb->ninput,1,MPI_INT,0,world);
  MPI_Bcast(&tb->nl,1,MPI_INT,0,world);
  MPI_Bcast(&tb->nq,1,MPI_INT,0,world);

  int me;
  MPI_Comm_rank(world,&me);
  if (me > 0) {
    memory->create(tb->l,tb->ninput,"angle:l");
    memory->create(tb->q,tb->ninput,"angle:q");
    memory->create(tb->e,tb->ninput,"angle:e");
    memory->create(tb->el,tb->ninput,"angle:el");
    memory->create(tb->eq,tb->ninput,"angle:eq");
    memory->create(tb->elq,tb->ninput,"angle:elq");
  }

  MPI_Bcast(tb->l,tb->ninput,MPI_DOUBLE,0,world);
  MPI_Bcast(tb->q,tb->ninput,MPI_DOUBLE,0,world);
  MPI_Bcast(tb->e,tb->ninput,MPI_DOUBLE,0,world);
  MPI_Bcast(tb->el,tb->ninput,MPI_DOUBLE,0,world);
  MPI_Bcast(tb->eq,tb->ninput,MPI_DOUBLE,0,world);
  MPI_Bcast(tb->elq,tb->ninput,MPI_DOUBLE,0,world);
  MPI_Bcast(&tb->theta0,1,MPI_DOUBLE,0,world);
}

void CgAngle::ufl_lookup(int type, double l, double q, double &u, double &fl, double &ft)
{
  Table *tb = &tables[tabindex[type]];
#ifdef OUTPUT_DEBUG_INFO
  if ((l > tb->l[tb->ninput-1]) || (l < tb->l[0])){
    std::cout << "l is out of range!!\n" << l << '\n';
  }
  if ((q > tb->q[tb->ninput-1]) || (q < tb->q[0])){
    std::cout << "theta is out of range!!\n" << q/MY_PI*180.0 << '\n';
  }
#endif

  // Map l and q to fit to table grid.
  l = std::min(l, tb->l[tb->ninput-1]);
  l = std::max(l, tb->l[0]);
  q = std::min(q, tb->q[tb->ninput-1]);
  q = std::max(q, tb->q[0]);

  int il = static_cast<int>((l-tb->l[0]) / tb->dl);
  int iq = static_cast<int>((q-tb->q[0]) / tb->dq);
  //   4----3     q
  //   |    |   ^  
  //   |    |   |
  //   1----2   --->  l
  int i1 = iq*tb->nl + il;
  int i2 = i1 + 1;
  int i3 = i1 + tb->nl + 1;
  int i4 = i1 + tb->nl;

  // Fractional coordinates of (l,q) in cell.
  double xi = (l - tb->l[i1]) / tb->dl;
  double eta = (q - tb->q[i1]) / tb->dq;
  if (tabstyle == LINEAR)
  {
    const double N1 = (1.0 - xi) * (1.0 - eta);
    const double N2 = xi * (1.0 - eta);
    const double N3 = xi * eta;
    const double N4 = (1.0 - xi) * eta;

    u  =   N1*tb->e[i1]  + N2*tb->e[i2]  + N3*tb->e[i3]  + N4*tb->e[i4];
    fl =  (N1*tb->el[i1] + N2*tb->el[i2] + N3*tb->el[i3] + N4*tb->el[i4]);
    ft =  (N1*tb->eq[i1] + N2*tb->eq[i2] + N3*tb->eq[i3] + N4*tb->eq[i4]);
  }
  else if (tabstyle == BICUBIC) 
  {
    // Defines interpolation matrices from:
    // https://en.wikipedia.org/wiki/Bicubic_interpolation.
    Matrix4 v1;
    v1(0,0) = 1.0; v1(0,1) = 0.0; v1(0,2) = 0.0; v1(0,3) = 0.0;
    v1(1,0) = 0.0; v1(1,1) = 0.0; v1(1,2) = 1.0; v1(1,3) = 0.0;
    v1(2,0) =-3.0; v1(2,1) = 3.0; v1(2,2) =-2.0; v1(2,3) =-1.0;
    v1(3,0) = 2.0; v1(3,1) =-2.0; v1(3,2) = 1.0; v1(3,3) = 1.0;
    Matrix4 v2;
    v2(0,0) = 1.0; v2(0,1) = 0.0; v2(0,2) =-3.0; v2(0,3) = 2.0;
    v2(1,0) = 0.0; v2(1,1) = 0.0; v2(1,2) = 3.0; v2(1,3) =-2.0;
    v2(2,0) = 0.0; v2(2,1) = 1.0; v2(2,2) =-2.0; v2(2,3) = 1.0;
    v2(3,0) = 0.0; v2(3,1) = 0.0; v2(3,2) =-1.0; v2(3,3) = 1.0;
    // F contains the values from the BA table for the current grid point.
    Matrix4 F;
    // Top left quadrant.
    F(0,0) = tb->e[i1];   F(0,1) = tb->e[i4];
    F(1,0) = tb->e[i2];   F(1,1) = tb->e[i3];
    // Top right quadrant
    F(0,2) = tb->eq[i1]*tb->dq;  F(0,3) = tb->eq[i4]*tb->dq;
    F(1,2) = tb->eq[i2]*tb->dq;  F(1,3) = tb->eq[i3]*tb->dq;
    // Bottom left quadrant.
    F(2,0) = tb->el[i1]*tb->dl;  F(2,1) = tb->el[i4]*tb->dl;
    F(3,0) = tb->el[i2]*tb->dl;  F(3,1) = tb->el[i3]*tb->dl;
    // Bottom right quadrant.
    F(2,2) = tb->elq[i1]*tb->dl*tb->dq; F(2,3) = tb->elq[i4]*tb->dl*tb->dq;
    F(3,2) = tb->elq[i2]*tb->dl*tb->dq; F(3,3) = tb->elq[i3]*tb->dl*tb->dq;

    const Matrix4 A = v1 * F * v2;
    const Vec4 px(1.0, xi,  xi*xi,   xi*xi*xi);
    const Vec4 py(1.0, eta, eta*eta, eta*eta*eta);
    const Vec4 dx(0.0, 1.0, 2.0*xi,  3.0*xi*xi);
    const Vec4 dy(0.0, 1.0, 2.0*eta, 3.0*eta*eta);
    u = dot(px, A*py);
    fl = dot(dx, A*py) / tb->dl;
    ft = dot(px, A*dy) / tb->dq;
  }
#ifdef OUTPUT_DEBUG_INFO
    std::cout << "U = " << u << "\n";
    std::cout << "Fl = " << fl << '\n';
    std::cout << "Fq = " << ft << '\n';
#endif
}

