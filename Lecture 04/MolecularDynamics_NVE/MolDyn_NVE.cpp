/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
#include <stdlib.h>     // srand, rand: to generate random number
#include <iostream>     // cin, cout: Standard Input/Output Streams Library
#include <fstream>      // Stream class to both read and write from/to files.
#include <cmath>        // rint, pow
#include <iomanip>
#include "MolDyn_NVE.h"

using namespace std;

/*-------------------------------------------------------------Main-------------------------------------------------------------*/

int main(){ 
  Input();             //Inizialization
  int nconf = 1;
  for(int iblk=1; iblk <= nblk; ++iblk) //Simulation
  {
    cout << "Block number " << iblk << endl;
    Reset(iblk);   //Reset block averages
    for(int istep=1; istep <= nstep; ++istep)
    {
      Move();
      Measure();
      Accumulate(); //Update block averages
      if(istep%iprint == 0) cout << "Number of time-steps: " << istep << endl;
      if(istep%10 == 0){
//        ConfXYZ(nconf);//Write actual configuration in XYZ format //Commented to avoid "filesystem full"! 
        nconf += 1;
      }
    }
    Averages(iblk);   //Print results for current block
  }
  ConfFinal(); //Write final configuration to restart

  return 0;
}

/*-----------------------------------------------------------Functions-----------------------------------------------------------*/

void Input(void){ //Prepare all stuff for the simulation
  ifstream ReadInput,ReadConf;
  double ep, ek, pr, et, vir;

  cout << "Classic Lennard-Jones fluid        " << endl;
  cout << "Molecular dynamics simulation in NVE ensemble  " << endl << endl;
  cout << "Interatomic potential v(r) = 4 * [(1/r)^12 - (1/r)^6]" << endl << endl;
  cout << "The program uses Lennard-Jones units " << endl;

  seed = 1;    //Set seed for random numbers
  srand(seed); //Initialize random number generator
  
  ReadInput.open("input.dat"); //Read input

  ReadInput >> temp;

  ReadInput >> npart;
  cout << "Number of particles = " << npart << endl;

  ReadInput >> rho;
  cout << "Density of particles = " << rho << endl;
  vol = (double)npart/rho;
  cout << "Volume of the simulation box = " << vol << endl;
  box = pow(vol,1.0/3.0);
  cout << "Edge of the simulation box = " << box << endl<< endl;

  ReadInput >> rcut;
  ReadInput >> delta;
  ReadInput >> nstep;
  ReadInput >> iprint;
  ReadInput >> nblk; //Read the number of blocks for blocking averages
  ReadInput >> oldconf; //Read instructions whether or not to read a -dt confiuration

  cout << "The program integrates Newton equations with the Verlet method " << endl;
  cout << "Time step = " << delta << endl;
  cout << "Number of steps = " << nstep << endl << endl;
  cout << "Number of blocks for data averages = " << nblk << endl << endl;
  ReadInput.close();

//Prepare array for measurements
  iv = 0; //Potential energy
  ik = 1; //Kinetic energy
  ie = 2; //Total energy
  it = 3; //Temperature
  iw = 4; //Temperature
  n_props = 5; //Number of observables

//measurement of g(r)
  igofr = 5;
  nbins = 100;
  n_props = n_props + nbins;
  bin_size = (box/2.0)/(double)nbins;


//Read initial configuration
  cout << "Read initial configuration from file config.0 " << endl << endl;
  ReadConf.open("config.0");
  for (int i=0; i<npart; ++i){
    ReadConf >> x[i] >> y[i] >> z[i];
    x[i] = x[i] * box;
    y[i] = y[i] * box;
    z[i] = z[i] * box;
  }
  ReadConf.close();

  if(oldconf==0){    //If an old configuration is not accepted , prepare initial velocities
    cout << "Prepare random velocities with center of mass velocity equal to zero " << endl << endl;
    double sumv[3] = {0.0, 0.0, 0.0};
    for (int i=0; i<npart; ++i){
      vx[i] = rand()/double(RAND_MAX) - 0.5;
      vy[i] = rand()/double(RAND_MAX) - 0.5;
      vz[i] = rand()/double(RAND_MAX) - 0.5;

      sumv[0] += vx[i];
      sumv[1] += vy[i];
      sumv[2] += vz[i];
    }
    for (int idim=0; idim<3; ++idim) sumv[idim] /= (double)npart;

    double sumv2 = 0.0, fs;
    for (int i=0; i<npart; ++i){
      vx[i] = vx[i] - sumv[0];
      vy[i] = vy[i] - sumv[1];
      vz[i] = vz[i] - sumv[2];

      sumv2 += vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i];
    }
    sumv2 /= (double)npart;

    fs = sqrt(3 * temp / sumv2);   // fs = velocity scale factor 
    for (int i=0; i<npart; ++i){
      vx[i] *= fs;
      vy[i] *= fs;
      vz[i] *= fs;

      xold[i] = Pbc(x[i] - vx[i] * delta);
      yold[i] = Pbc(y[i] - vy[i] * delta);
      zold[i] = Pbc(z[i] - vz[i] * delta);
    }
    return;
  }

  else if(oldconf==1) {    //If an old configuration is updloaded:
    //Read pre-initial configuration
    cout << "Read pre-initial configuration from file old.0 " << endl << endl;
    ReadConf.open("old.0");
    for (int i=0; i<npart; ++i){
      ReadConf >> xold[i] >> yold[i] >> zold[i];
      xold[i] = xold[i] * box;
      yold[i] = yold[i] * box;
      zold[i] = zold[i] * box;
    }
    ReadConf.close();

    //Compute r(dt) with the Verlet algorithm : r(dt) = 2*r(0) - r(-dt) + f(0)*dt^2
    Move(); 

    //Compute v(dt/2) and rescale it exploiting v = sqrt(3*T)
    double sumv2 = 0.0, fs;
    for (int i=0; i<npart; ++i){
      vx[i] = Pbc(x[i] - xold[i])/(2.0 * delta);
      vy[i] = Pbc(y[i] - yold[i])/(2.0 * delta);
      vz[i] = Pbc(z[i] - zold[i])/(2.0 * delta);
      sumv2 += vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i];
    }

    sumv2 /= (double)npart;        //<v^2> per particle
    fs = sqrt(3 * temp / sumv2);   // fs = velocity scale factor
    for (int i=0; i<npart; ++i){
      vx[i] *= fs;
      vy[i] *= fs;
      vz[i] *= fs;
      //Keep r(dt) and re-calculate r(0) 
      xold[i] = Pbc(x[i] - vx[i] * delta);
      yold[i] = Pbc(y[i] - vy[i] * delta);
      zold[i] = Pbc(z[i] - vz[i] * delta);
    }
    return;

  }
}


void Move(void){ //Move particles with Verlet algorithm
  double xnew, ynew, znew, fx[m_part], fy[m_part], fz[m_part];

  for(int i=0; i<npart; ++i){ //Force acting on particle i
    fx[i] = Force(i,0);
    fy[i] = Force(i,1);
    fz[i] = Force(i,2);
  }

  for(int i=0; i<npart; ++i){ //Verlet integration scheme

    xnew = Pbc( 2.0 * x[i] - xold[i] + fx[i] * pow(delta,2) );
    ynew = Pbc( 2.0 * y[i] - yold[i] + fy[i] * pow(delta,2) );
    znew = Pbc( 2.0 * z[i] - zold[i] + fz[i] * pow(delta,2) );

    vx[i] = Pbc(xnew - xold[i])/(2.0 * delta);
    vy[i] = Pbc(ynew - yold[i])/(2.0 * delta);
    vz[i] = Pbc(znew - zold[i])/(2.0 * delta);

    xold[i] = x[i];
    yold[i] = y[i];
    zold[i] = z[i];

    x[i] = xnew;
    y[i] = ynew;
    z[i] = znew;
  }
  return;
}


double Force(int ip, int idir){ //Compute forces as -Grad_ip V(r)
  double f=0.0;
  double dvec[3], dr;

  for (int i=0; i<npart; ++i){
    if(i != ip){
      dvec[0] = Pbc( x[ip] - x[i] );  // distance ip-i in pbc
      dvec[1] = Pbc( y[ip] - y[i] );
      dvec[2] = Pbc( z[ip] - z[i] );

      dr = dvec[0]*dvec[0] + dvec[1]*dvec[1] + dvec[2]*dvec[2];
      dr = sqrt(dr);

      if(dr < rcut){
        f += dvec[idir] * (48.0/pow(dr,14) - 24.0/pow(dr,8)); // -Grad_ip V(r)
      }
    }
  }
  
  return f;
}


void Measure(){ //Properties measurement
  int bin;
  double v, t, vij, w, wij;
  double dx, dy, dz, dr;
  v = 0.0; //reset observables
  t = 0.0;
  w = 0.0;

//reset the hystogram of g(r)
  for (int k=igofr; k<igofr+nbins; ++k) walker[k]=0.0;

//cycle over pairs of particles
  for (int i=0; i<npart-1; ++i){
    for (int j=i+1; j<npart; ++j){

     dx = Pbc( xold[i] - xold[j] ); // here I use old configurations [old = r(t)]
     dy = Pbc( yold[i] - yold[j] ); // to be compatible with EKin which uses v(t)
     dz = Pbc( zold[i] - zold[j] ); // => EPot should be computed with r(t)

     dr = dx*dx + dy*dy + dz*dz;
     dr = sqrt(dr);

//update of the histogram of g(r)
     for(bin=0 ; bin<nbins ; bin++){
       if( dr>= bin*bin_size && dr < (bin+1)*bin_size ) walker[igofr+bin]+=2.;
     }

     if(dr < rcut){
       vij = 1.0/pow(dr,12) - 1.0/pow(dr,6);
       wij = 1.0/pow(dr,12) - 0.5/pow(dr,6);

//Potential energy and Virial
       v += vij;
       w += wij;
     }
    }          
  }

//Kinetic energy
  for (int i=0; i<npart; ++i) t += (vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i]);
   
  walker[iv] = 4.0 * v;
  walker[ik] = 0.5 * t;
  walker[it] = 2./3. * walker[ik];
  walker[ie] = walker[iv]+ walker[ik];
  walker[iw] = 48.0 * w / 3.0;


// Instantaneous values
  ofstream Epot, Ekin, Etot, Temp, Press;
  Epot.open("output_epot.dat",ios::app);
  Ekin.open("output_ekin.dat",ios::app);
  Temp.open("output_temp.dat",ios::app);
  Etot.open("output_etot.dat",ios::app);
  Press.open("output_press.dat",ios::app);

  Epot << walker[iv]/(double)npart  << endl;	//Potential energy per particle
  Ekin << walker[ik]/(double)npart  << endl;	//Kinetic energy per particle
  Temp << walker[it]/(double)npart << endl;	//Temperature
  Etot << walker[ie]/(double)npart << endl;	//Total energy per particle
  Press << rho* walker[it]/(double)npart + walker[iw]/vol << endl;	//Pressure

  Epot.close();
  Ekin.close();
  Temp.close();
  Etot.close();
  Press.close();

  return;
}


void Reset(int iblk) //Reset block averages
{
   
  if(iblk == 1)
  {
    for(int i=0; i<n_props; ++i)
    {
       glob_av[i] = 0;
       glob_av2[i] = 0;
    }
  }

  for(int i=0; i<n_props; ++i)
  {
    blk_av[i] = 0;
  }
  blk_norm = 0;
}


void Accumulate(void) //Update block averages
{

  for(int i=0; i<n_props; ++i)
  {
    blk_av[i] = blk_av[i] + walker[i];
  }
  blk_norm = blk_norm + 1.0;
}


void Averages(int iblk) //Print results for current block
{
    
  double r[n_props-igofr], gdir, deltaV;
  ofstream Epot, Ekin, Etot, Temp, Press, Gofr;
  const int wd=12;
    
  Epot.open("ave_epot.dat",ios::app);
  Ekin.open("ave_ekin.dat",ios::app);
  Etot.open("ave_etot.dat",ios::app);
  Temp.open("ave_temp.dat",ios::app);
  Press.open("ave_press.dat",ios::app);
  Gofr.open("ave_gofr.dat",ios::app);
    
  stima_pot = blk_av[iv]/blk_norm/(double)npart; //Potential energy
  glob_av[iv] += stima_pot;
  glob_av2[iv] += stima_pot*stima_pot;
  err_pot=Error(glob_av[iv],glob_av2[iv],iblk);

  stima_kin = blk_av[ik]/blk_norm/(double)npart; //Kinetic energy
  glob_av[ik] += stima_kin;
  glob_av2[ik] += stima_kin*stima_kin;
  err_kin=Error(glob_av[ik],glob_av2[ik],iblk);

  stima_tot = blk_av[ie]/blk_norm/(double)npart; //Total energy
  glob_av[ie] += stima_tot;
  glob_av2[ie] += stima_tot*stima_tot;
  err_tot=Error(glob_av[ie],glob_av2[ie],iblk);

  stima_temp = blk_av[it]/blk_norm/(double)npart; //Temperature
  glob_av[it] += stima_temp;
  glob_av2[it] += stima_temp*stima_temp;
  err_temp=Error(glob_av[it],glob_av2[it],iblk);
    
  stima_press = rho * stima_temp + blk_av[iw]/blk_norm/ vol; //Pressure
  glob_av[iw] += stima_press;
  glob_av2[iw] += stima_press*stima_press;
  err_press=Error(glob_av[iw],glob_av2[iw],iblk);

  for(int ig=igofr; ig <n_props ; ig++){
    int j= ig-igofr;
    r[j]= (ig-igofr+0.5)*bin_size ;
    deltaV= 4.*pi/3. *( pow((r[j]+bin_size),3) - pow(r[j],3) );
    stima_gdir[ig] = blk_av[ig]/blk_norm/(rho*npart*deltaV); //Stima g(r)
    glob_av[ig] += stima_gdir[ig];
    glob_av2[ig] += stima_gdir[ig]*stima_gdir[ig];
    err_gdir[ig]=Error(glob_av[ig],glob_av2[ig],iblk);
  }

//Potential energy per particle
    Epot << setw(wd) << iblk <<  setw(wd) << stima_pot << setw(wd) << glob_av[iv]/(double)iblk << setw(wd) << err_pot << endl;
//Kinetic energy per particle
    Ekin << setw(wd) << iblk <<  setw(wd) << stima_kin << setw(wd) << glob_av[ik]/(double)iblk << setw(wd) << err_kin << endl;
//Total energy per particle
    Etot << setw(wd) << iblk <<  setw(wd) << stima_tot << setw(wd) << glob_av[ie]/(double)iblk << setw(wd) << err_tot << endl;
//Temperature
    Temp << setw(wd) << iblk <<  setw(wd) << stima_temp << setw(wd) << glob_av[it]/(double)iblk << setw(wd) << err_temp << endl;
//Pressure
    Press << setw(wd) << iblk <<  setw(wd) << stima_press << setw(wd) << glob_av[iw]/(double)iblk << setw(wd) << err_press << endl;

//g(r)
    if(iblk==nblk) {
      for(int ig=igofr; ig <n_props ; ig++) Gofr << r[ig-igofr]<< setw(wd) << glob_av[ig]/(double)nblk << setw(wd) << err_gdir[ig] << endl;
    }

    cout << "----------------------------" << endl << endl;

    Epot.close();
    Ekin.close();
    Etot.close();
    Temp.close();
    Press.close();
    Gofr.close();
}


void ConfFinal(void){ //Write final & old final configurations
  ofstream WriteConf;

  cout << "Print final configuration to file config.final " << endl << endl;
  WriteConf.open("config.final");

  for (int i=0; i<npart; ++i){
    WriteConf << x[i]/box << "   " <<  y[i]/box << "   " << z[i]/box << endl;
  }
  WriteConf.close();

  cout << "Print OLD final configuration to file old.final " << endl << endl;
  WriteConf.open("old.final");

  for (int i=0; i<npart; ++i){
    WriteConf << xold[i]/box << "   " <<  yold[i]/box << "   " << zold[i]/box << endl;
  }
  WriteConf.close();

  return;
}


void ConfXYZ(int nconf){ //Write configuration in .xyz format
  ofstream WriteXYZ;

  WriteXYZ.open("frames/config_" + to_string(nconf) + ".xyz");
  WriteXYZ << npart << endl;
  WriteXYZ << "This is only a comment!" << endl;
  for (int i=0; i<npart; ++i){
    WriteXYZ << "LJ  " << Pbc(x[i]) << "   " <<  Pbc(y[i]) << "   " << Pbc(z[i]) << endl;
  }
  WriteXYZ.close();
}


double Pbc(double r){  //Algorithm for periodic boundary conditions with side L=box
    return r - box * rint(r/box);
}


double Error(double sum, double sum2, int iblk)
{
    if( iblk == 1 ) return 0.0;
    else return sqrt((sum2/(double)iblk - pow(sum/(double)iblk,2))/(double)(iblk-1));
}

/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
