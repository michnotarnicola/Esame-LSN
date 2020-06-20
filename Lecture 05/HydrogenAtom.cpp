#include <iostream>
#include <stdlib.h> 
#include <iomanip>
#include <fstream>
#define _USE_MATH_DEFINES
#include <cmath>
#include <array>
#include <algorithm>
#include "random.h"

using namespace std;
Random rnd;
/*--------------------------------------------------------Global Variables-------------------------------------------------------*/
double accept, attempt;
const int M=1000000, Nblk=100, Lblk =M/Nblk, N_step=50;//#simulazioni, #blocchi, lunghezza blocco, numero step Metropolis nel blocco
double x[Lblk], y[Lblk], z[Lblk] , x_new[Lblk], y_new[Lblk], z_new[Lblk], ave_r[Nblk], ave_r2[Nblk];//arrays per campionamento e medie

/*-----------------------------------------------------------Functions-----------------------------------------------------------*/

double errore(double AV, double AV2,int n) {  /*Function per la stima dell'incertezza statistica*/ 
  if(n==0) return 0 ;
  else return sqrt((AV2 - AV*AV)/n) ;
}


double prob1s (double x, double y, double z) {  /*Function per il calcolo della densità di probabiità 1s nel punto x*/ 
  double r= sqrt(x*x+y*y+z*z);
  return 1./M_PI * exp(-2*r);
}


double prob2p (double x, double y, double z) {  /*Function per il calcolo della densità di probabiità 2p nel punto x*/
  double r= sqrt(x*x+y*y+z*z);
  return 1./(32*M_PI)* z*z * exp(-r);
}


void Attempt(int istep, int orbital, int Tyx) { /*Function per generare le mosse del Metropolis*/
  
 if(orbital==0 && Tyx==0) {			//psi 1s T(y|x) uniforme
   x_new[istep]= x[istep] + 2.4*(rnd.Rannyu()-0.5);
   y_new[istep]= y[istep] + 2.4*(rnd.Rannyu()-0.5);
   z_new[istep]= z[istep] + 2.4*(rnd.Rannyu()-0.5);
 }
 if(orbital==0 && Tyx==1) {			//psi 1s T(y|x) gaussiana
   x_new[istep]= rnd.Gauss(x[istep],0.75);
   y_new[istep]= rnd.Gauss(y[istep],0.75);
   z_new[istep]= rnd.Gauss(z[istep],0.75);
 }
 if(orbital==1 && Tyx==0) {			//psi 2p T(y|x) uniforme
   x_new[istep]= x[istep] + 6*(rnd.Rannyu()-0.5);
   y_new[istep]= y[istep] + 6*(rnd.Rannyu()-0.5);
   z_new[istep]= z[istep] + 6*(rnd.Rannyu()-0.5);
 }
 if(orbital==1 && Tyx==1) {			//psi 2p T(y|x) gaussiana
   x_new[istep]= rnd.Gauss(x[istep],1.8);
   y_new[istep]= rnd.Gauss(y[istep],1.8);
   z_new[istep]= rnd.Gauss(z[istep],1.8);
 }
}


void Metropolis (int orbital, int Tyx) {	/*Function per l'algoritmo di Metropolis*/
  double p_old=0,p_new=0,alpha;
  accept=0, attempt=0;

  for(int i=0; i<Lblk; i++) {
    Attempt(i, orbital, Tyx);		//Propongo una mossa  
    if(orbital==0) {
      p_old=prob1s(x[i], y[i], z[i]);
      p_new=prob1s(x_new[i], y_new[i], z_new[i]);
    }
    if(orbital==1) {
      p_old=prob2p(x[i], y[i], z[i]);
      p_new=prob2p(x_new[i], y_new[i], z_new[i]);
    }
    alpha=min(1., p_new/p_old);
    if(rnd.Rannyu() <= alpha) {		//La accetto con probabilità alpha
      x[i]=x_new[i];
      y[i]=y_new[i];
      z[i]=z_new[i];
      accept+=1; 
    }
    attempt+=1;
  }
}


void Equilibration (int orbital, int Tyx) {	/*Function per l'inizializzazione dei dati e la fase di equilibrazione*/

  std::fill_n(ave_r, Nblk, 0);	//Inizializzo gli arrays dei valori medi a 0
  std::fill_n(ave_r2, Nblk, 0);
  
  if(orbital==0) {		//Nel caso 1s inizializzo il sistema in (1,0,0)
    std::fill_n(x, Lblk, 1);
    std::fill_n(y, Lblk, 0);
    std::fill_n(z, Lblk, 0);
  }  
  if(orbital==1) {		//Nel caso 1s inizializzo il sistema in (0,0,4)
    std::fill_n(x, Lblk, 0);
    std::fill_n(y, Lblk, 0);
    std::fill_n(z, Lblk, 4);
  }
 
  int N_eq= 1000; 		//Equilibration steps
  cout << "Fase di equilibrazione ... " << endl; 
  for(int ieq=0; ieq<N_eq; ieq++) Metropolis(orbital, Tyx);
  cout << "Equilibrazione completata" << endl;
}


void Wavefunction (int orbital) {	/*Function per la scrittura su file delle configurazioni campionate (solo nel caso Tyx uniforme)*/
  ofstream psi;
  if(orbital==0) psi.open("psi1s.dat", ios::app);
  if(orbital==1) psi.open("psi2p.dat", ios::app);
  for(int i=0; i<Lblk; i++) psi << x[i] << setw(14)<< y[i] << setw(14) << z[i] << setw(14)<< sqrt(x[i]*x[i] + y[i]*y[i]+ z[i]*z[i])<< endl;
  psi.close();
}


void Measure(int iblk) {	/*Function per il calcolo di <r> e <r^2>*/
  for(int i=0; i<Lblk; i++) ave_r[iblk]+= sqrt(x[i]*x[i] + y[i]*y[i]+ z[i]*z[i]);
  ave_r[iblk]/=Lblk;
  ave_r2[iblk]= ave_r[iblk]* ave_r[iblk];
}


void Averages(int orbital, int Tyx) {	/*Function per le medie cumulate*/
  ofstream write;
  if(orbital==0 && Tyx==0) write.open("rmedio1s_u.dat");
  if(orbital==0 && Tyx==1) write.open("rmedio1s_g.dat");
  if(orbital==1 && Tyx==0) write.open("rmedio2p_u.dat");
  if(orbital==1 && Tyx==1) write.open("rmedio2p_g.dat");

  double stima_r, stima_r2, sigma_r;
  for(int iblk=0; iblk<Nblk; iblk++) {
    stima_r=0, stima_r2=0, sigma_r=0;
    for(int jblk=0; jblk<=iblk; jblk++) {
      stima_r+=ave_r[jblk];
      stima_r2+=ave_r2[jblk];
    }
    stima_r/=(iblk+1);
    stima_r2/=(iblk+1);
    sigma_r= errore(stima_r, stima_r2, iblk);
    write << iblk+1 << setw(14) << ave_r[iblk] << setw(14) << stima_r << setw(14) << sigma_r << endl;
  }   
  write.close();
}

/*-------------------------------------------------------------Main-------------------------------------------------------------*/

int main(){

//Seed per il generatore
   int seed[4];
   int p1, p2;
   ifstream Primes("Primes");
   Primes >> p1 >> p2 ;
   Primes.close();

   ifstream input("seed.in");
   input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
   rnd.SetRandom(seed,p1,p2);
   input.close();
   rnd.SaveSeed();
/*---------------*/

  for(int orbital=0; orbital<=1; orbital++) {		//Orbitale 1s (orbital=0) o orbitale 2p (orbital=1)
   for(int Tyx=0; Tyx<=1; Tyx++) {			//Probabilità di transizione uniforme (Tyx=0) o gaussiana (Tyx=1)
   
     if(orbital==0 && Tyx==0) cout << "Simulazione 1: orbitale 1s. Distribuzione T(y|x) uniforme" << endl << endl;
     if(orbital==0 && Tyx==1) cout << "Simulazione 1: orbitale 1s. Distribuzione T(y|x) gaussiana" << endl << endl;
     if(orbital==1 && Tyx==0) cout << "Simulazione 2: orbitale 2p. Distribuzione T(y|x) uniforme" << endl << endl;
     if(orbital==1 && Tyx==1) cout << "Simulazione 2: orbitale 2p. Distribuzione T(y|x) uniforme" << endl << endl;

     Equilibration(orbital, Tyx);
     for(int iblk=0; iblk<Nblk; iblk++) {
       for(int istep=0; istep<N_step; istep++) Metropolis(orbital, Tyx) ;
       cout << "Accettazione del blocco " << iblk+1 << " =" << accept/attempt << endl;
       if(Tyx==0) Wavefunction(orbital);
       Measure(iblk);
     }
     Averages(orbital, Tyx);
     cout << "--------------------------------" << endl;
   }
  }
  return 0;

}
	
	
