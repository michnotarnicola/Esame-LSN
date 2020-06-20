/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
#include <string>
//parameters, observables
const int m_props=1000;
int n_props, nbins;
int iv,ik,it,ie,iw, igofr;
double walker[m_props], bin_size;

// averages
double blk_av[m_props],blk_norm,accepted,attempted;
double glob_av[m_props],glob_av2[m_props];
double stima_pot, stima_kin, stima_tot, stima_temp, stima_press,stima_gdir[m_props];
double err_pot, err_kin, err_tot, err_temp, err_press, err_gdir[m_props];

//configuration
const int m_part=108;
double x[m_part],y[m_part],z[m_part],xold[m_part],yold[m_part],zold[m_part];
double vx[m_part],vy[m_part],vz[m_part];

// thermodynamical state
int npart;
double energy,temp,vol,rho,box,rcut;

// simulation
int nstep, iprint, nblk, seed;
double delta;
bool oldconf;

//pigreco
const double pi=3.1415927;

//functions
void Input(void);
void Reset(int);
void Accumulate(void);
void Averages(int);
void Move(void);
void ConfFinal(void);
void ConfXYZ(int);
void Measure(void);
double Force(int, int);
double Pbc(double);
double Error(double, double, int);
/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
