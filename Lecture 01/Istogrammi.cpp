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

double* randomgenerator(double *r, int M, int distrib) { /* Function per generare numeri casuali uniformi in [0,1]:
							 input -> array per numeri casuali, #lanci, distribuzione da generare */
  const double pi= M_PI ;
  for(int i=0; i<M ; i++) { 
    if(distrib==0) r[i]= rnd.Rannyu() ;
    if(distrib==1) r[i]= tan(pi*(rnd.Rannyu()-0.5)) ;
    if(distrib==2) r[i]= - log(1-rnd.Rannyu()) ;
  }
  return r ; 

}


double* S_Ncalculator (double *rand, int nblock, int L, int Lmax, double* S_N) {   /* Function per il calcolo di S_N in ogni esperimento:
										   input -> array di numeri casuali, numero blocchi, numero di
										            termini da sommare e numero massimo di termini
										   output -> array di S_N in ogni esperimento */
  for(int i=0; i< nblock; i++){    //Fisso il blocco
    for(int j=0; j<L; j++){       //Sommo i primi L numeri nel blocco
      int k=j+i*Lmax ;
      S_N[i]+= rand[k]/L;
    }
  }

  return S_N ;
}


void write (double* S1, double* S2, double* S3, double* S4, int M, basic_string<char> distr) {/*Function per la scrittura dati su file:
                                                                                             input-> arrays di S_N, nome file da scrivere*/
   
  ofstream myfile;
  myfile.open (distr); 
  for(int j=0; j<M; j++) myfile << S1[j] << " " << S2[j] << " " << S3[j] << " " << S4[j] << endl;
  myfile.close();
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

  int M=10000 ;            		 //Numero realizzazioni (blocchi)
  int L[]={1,2,10,100} ;                //Numero lanci per blocco
  double rand[M*100]={0}; 		 //Array per numeri casuali

  for(int k=0; k<3 ; k++) {     	 //Genero le tre distribuzioni una per volta (0=unif , 1=lorentz, 2=esponenziale) 

    double S_N_1[M]={0};
    double S_N_2[M]={0};
    double S_N_3[M]={0};
    double S_N_4[M]={0};

    randomgenerator(rand, M*100, k) ;

    S_Ncalculator (rand, M, L[0], L[3], S_N_1);
    S_Ncalculator (rand, M, L[1], L[3], S_N_2);
    S_Ncalculator (rand, M, L[2], L[3], S_N_3);
    S_Ncalculator (rand, M, L[3], L[3], S_N_4);
    write(S_N_1, S_N_2, S_N_3, S_N_4, M, "Distribuzione" + to_string(k) + ".dat");
  }

  return 0;

}
	
	
