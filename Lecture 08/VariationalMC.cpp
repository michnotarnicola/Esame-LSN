#include <iostream>
#include <fstream>
#include <ostream>
#define _USE_MATH_DEFINES
#include <cmath>
#include <iomanip>
#include <algorithm>
#include <stdlib.h> 
#include "random.h"

using namespace std;
Random rnd;
/*-------------------------------------------------------Global Variables--------------------------------------------------------*/
double mu_tent, sigma_tent;		//Valori tentativi per mu e sigma ottenuti con mosse Metropolis
double mu=1, sigma=1;			//Dato iniziale
double delta=3.5;
double accepted=0, attempted=0;
double hbar=1 , m=1;

/*-----------------------------------------------------------Functions-----------------------------------------------------------*/

double errore(double *AV, double *AV2,int n) { 
   if(n==0) return 0 ;
   else return sqrt(abs( (AV2[n] - AV[n]*AV[n]) )/n) ;
}


double psitrial(double x) {
   return exp(- pow((x-mu_tent),2)/(2*pow(sigma_tent,2) ) ) + exp(- pow((x+mu_tent),2)/(2*pow(sigma_tent,2) ) );
//   return exp(- x*x/2);
}


double second_der_psi (double x) {
  return -1./pow(sigma_tent,2) * psitrial(x) + 1./pow(sigma_tent,4) * ( pow((x-mu_tent),2)* exp(-pow((x-mu_tent),2)/(2*pow(sigma_tent,2))) + pow((x+mu_tent),2)* exp(-pow((x+mu_tent),2)/(2*pow(sigma_tent,2))) );
//  return (x*x-1)*exp(-x*x/2);
}


double E_pot (double x) {
  return pow(x,4) -5./2. * pow(x,2) ;
//  return 1./2*x*x;
}


double* Metropolis (double* x, int Blocklenght) { /*Function che implementa l'algoritmo di Metroplis, proponendo mosse e accettandole*/

  double p_old, p_new, alpha ;
  double xnew[Blocklenght]={0};
  accepted=0, attempted=0;

  for(int i=0; i<Blocklenght ; i++) {
    xnew[i]= x[i] + delta*(rnd.Rannyu()-0.5);	//Genero xnew uniforme in [-0.5*delta +x , 0.5*delta +x]
    p_old = pow(psitrial(x[i]),2);
    p_new = pow(psitrial(xnew[i]),2);
    alpha =min(1., p_new/p_old);
    if(rnd.Rannyu() <= alpha) {
      accepted+=1;
      x[i]=xnew[i];
    }
    attempted+=1;
  }
  return x;
}


double H_average(double* x, int Blocklenght) { /*Function per calcolare <H>*/

  double lambda= hbar*hbar/(2*m);
  double E=0;
  for(int j=0; j<Blocklenght ; j++) E += -lambda*second_der_psi(x[j]) / psitrial(x[j]) + E_pot(x[j]);
  return  E/Blocklenght;
}


void Blockingaverage(double *ave, double *ave2, int Nblk, std::ofstream& myfile) { /* Function per il calcolo delle medie cumulate:
							 			   input -> arrays di r e r^2, # blocchi, nome ofstream */
  const int wd=14;
  double sum_prog[Nblk]={0} ;
  double sum2_prog[Nblk]={0} ;
  double err_prog[Nblk]={0} ;

  for(int i=0; i<Nblk ; i++) {    	
    for(int j=0; j<i+1 ; j++) { 
      sum_prog[i] += ave[j] ;      
      sum2_prog[i] += ave2[j]  ; 
    }
    sum_prog[i]/=(i+1) ;
    sum2_prog[i]/=(i+1) ; 
    err_prog[i] = errore(sum_prog, sum2_prog, i) ; 
    myfile << i+1 << setw(wd) << ave[i] << setw(wd) << sum_prog[i] << setw(wd) << err_prog[i] << endl;
   }
}


void wavefunc (double *x, int Blocklenght, std::ofstream& myfile) { /*Func per la scrittura su file delle coordinate campionate
									input -> arrays di coordinate, L blocco, nome ofstream */

  for (int j=0 ; j<Blocklenght ; j++) myfile << x[j] << endl;
}


double Move () {	/*Function che realizza la simulazione MC (da usare nella procedura di annealing)*/

  int Nblk= 50; 				//Numero blocchi
  int Blocklenght=2000;
  int Eq_step= 1000, Nstep=50;
  double x[Blocklenght]={0};			//Array per campionare la distribuzione di probabilità 
  double H[Nblk]={0};				//Arrays per i valori medi
  double Hmedio=0;

//Equilibrazione    
  std::fill_n(x, Blocklenght, 0); 		//Data la simmetria pari del potenziale inizializzo il sistema nell'origine				
  for(int istep=0 ; istep <Eq_step; istep++) Metropolis(x, Blocklenght);

//Simulazione
  for(int iblk=0 ; iblk < Nblk ; iblk++) {
    for(int istep=0 ; istep < Nstep ; istep++) Metropolis(x, Blocklenght);
    H[iblk]= H_average(x, Blocklenght);
    Hmedio+=H[iblk];
  }
  return Hmedio/Nblk ;
}


void Simulation () {	/*Function che realizza la simulazione FINALE MC (comprende anche la scrittura dati su file)*/

  int Nblk= 50; 				//Numero blocchi
  int Blocklenght=2000;
  int Eq_step= 1000, Nstep=50;
  double x[Blocklenght]={0};			//Array per campionare la distribuzione di probabilità 
  double H[Nblk]={0};				//Arrays per i valori medi
  double H2[Nblk]={0};

//Equilibrazione    
  std::fill_n(x, Blocklenght, 0); 		//Data la simmetria pari del potenziale inizializzo il sistema nell'origine
  cout << "Fase di equilibrazione ... " << endl; 				
  for(int istep=0 ; istep <Eq_step; istep++) Metropolis(x, Blocklenght);
  cout << "Equilibrazione completata " << endl; 
  cout << "Inizio simulazione " << endl; 
  cout << "------------------------------------" << endl;

//Simulazione
  ofstream Energy;
  ofstream Probdensity;
  Energy.open("Hmedio.dat");
  Probdensity.open("Psi_T.dat");

  for(int iblk=0 ; iblk < Nblk ; iblk++) {
    cout << "Blocco numero " << iblk+1 << endl;
    for(int istep=0 ; istep < Nstep ; istep++) Metropolis(x, Blocklenght);
    cout << "Accettazione = " << accepted/attempted << endl;
    cout << "------------------------------------" << endl;

    wavefunc(x, Blocklenght, Probdensity);
    H[iblk]= H_average(x, Blocklenght);
    H2[iblk]= H[iblk]*H[iblk];
  }

  Blockingaverage(H, H2, Nblk, Energy);		
  Probdensity.close(); 
  Energy.close() ;
}

void Annealing() {  		/*Function che implementa il SA*/
  double temp, beta, temp_i=0.1, temp_f=0.01, temp_step=50;
  double delta_mu=0.1, delta_sigma=0.1, Ncycle=50;
  double E_new, E_old , p ;

  ofstream annealing;
  annealing.open("AndamentoH.dat");
  
  mu_tent=mu;							//Con il dato iniziale calcolo E_old
  sigma_tent=sigma;
  E_old= Move(); 

  cout << "-------------Simulated Annealing-------------" << endl;
  for(int itemp=0; itemp <= temp_step; itemp++) {		//Ciclo sulla temperatura fittizia 
    annealing << itemp << setw(14) << E_old << endl;		//Scrivo su file il valore di <H>
  
    temp = temp_i- (temp_i-temp_f)/temp_step * itemp ;
    beta=1/temp;
    
    for(int icycle=0; icycle < Ncycle ; icycle ++) {
      mu_tent= mu + delta_mu*(rnd.Rannyu()-0.5) ; 		//mu_tent in [mu- delta_mu, mu]
      sigma_tent= sigma + delta_sigma*(rnd.Rannyu()-0.5) ;  	//sigma_tent in [sigma- delta_sigma, sigma]
      E_new= Move();
      p=min(1., exp(-beta*(E_new-E_old)) );			//Metropolis!
      if(rnd.Rannyu()<= p) {
        mu=mu_tent ;
        sigma=sigma_tent ;
        E_old=E_new;
      }
    }
    cout << "Step numero " << itemp+1 << " Temperatura: " << temp <<" Mu= " << mu << " Sigma= " << sigma<< endl;
  }
  cout << "-------------------------------------------" << endl;
  annealing.close();
}


/*-------------------------------------------------------------Main-------------------------------------------------------------*/

int main(){

/* Seed per il generatore */
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

  Annealing();
  cout << "Valori minimi trovati :" <<endl;
  cout << "Mu		= " << mu << endl;
  cout << "Sigma	= " << sigma << endl<<endl;
  
  mu_tent=mu;
  sigma_tent=sigma;
  Simulation();
  
  return 0;
  
}
	
	
