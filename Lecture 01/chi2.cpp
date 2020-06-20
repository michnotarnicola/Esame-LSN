#include <iostream>
#include <fstream>
#include <string>
#define _USE_MATH_DEFINES
#include <cmath>
#include <array>
#include <random>
#include "random.h"

using namespace std;
Random rnd;

/*-----------------------------------------------------------Functions-----------------------------------------------------------*/

double* randomgenerator(double *r, int M) { 	/* Function per generare numeri casuali uniformi in [0,1]:
						input -> array per numeri casuali, #lanci */
  for(int i=0; i<M ; i++) r[i]= rnd.Rannyu();
  return r ;
}


int* contatorelanci (double *rand, int nlanci, int m_int, int* counter){/*Function per contare il numero di eventi in ogni sottointervallo:
									   input-> array di numeri casuali nell'esperimento, #lanci ,
										   #sottointervalli */
  for(int i=0; i<m_int ; i++) {    //Fisso un sottointervallo
    for(int j=0; j<nlanci; j++) { //Cerco nell'array tutti gli eventi che cadono in esso
      if(rand[j]<= (i+1.)/m_int && rand[j] > double(i)/m_int) counter[i]+=1 ;
    }
  }
  return counter;
}


double chiquadro(int* counter ,int size, double vatteso) {   /*Function per il calcolo di chi^2:
                     					     input-> array dei conteggi in ogni sottointervallo, #sottointervalli, valore di
								     apettazione di n_i*/
  double chi2=0  ;
  for(int i=0; i<size ; i++) {   
    chi2+= (counter[i] -vatteso)*(counter[i] -vatteso)/vatteso ;
  }
  return chi2 ;
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
/*-----------------*/

  int M=100 ; 			//Numero di intervalli in cui divido [0,1]
  int n=10000 ; 		//Numero di lanci
  int Nstime= 100 ;		//Numero stime di chi^2 da ottenere (i.e. numero esperimenti)
  double vasp=n/M ;		//Valore di aspettazione del numero di eventi in un sottointervallo
  double rand[n*Nstime]={0} ;	//Array per numeri casuali (n per ogni stima da ottnere)
  double rand1[n]={0} ;		//Array ausiliario
  double chi2=0 ;

  randomgenerator(rand, n*Nstime) ;
  ofstream outfile ;
  outfile.open ("chiquadro.dat") ;
   
  for(int i=0; i<Nstime; i++){	  //Ciclo per ottenere una stima di chi^2: fisso l'esperimento i

    int n_i[M]={0} ; 		  //Array per il numero di conteggi in ogni sottointervallo
    for(int j=0; j<n ; j++) {	  //Array che contiene gli n numeri casuali dell'esperimento considerato
      rand1[j]=rand[j+i*n] ;
    }
    contatorelanci(rand1 ,n, M, n_i) ;
    chi2 = chiquadro(n_i ,M, vasp) ;
    outfile << i+1 << " " << chi2 -100 << endl ;
  }

  outfile.close() ;

  return 0;

}
