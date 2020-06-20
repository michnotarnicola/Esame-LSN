#include <iostream>
#include <fstream>
#include <string>
#define _USE_MATH_DEFINES
#include <cmath>
#include <array>
#include <random>
#include <algorithm>
#include "random.h"

using namespace std;
Random rnd;

/*-----------------------------------------------------------Functions-----------------------------------------------------------*/

double errore(double *AV, double *AV2,int n) {  /*Function per la stima dell'incertezza statistica: 
						input-> <r> e <r^2> calcolati al blocco n*/
  if(n==0) return 0 ;
  else return sqrt((AV2[n] - AV[n]*AV[n])/n) ;
}


double* boxmuller(double *array, int M, double mu, double sigma) { /* Function per generare numeri casuali N(mu,sigma):
							            input -> array per numeri casuali, #lanci, mu & sigma */
  for(int i=0; i<M ; i++) array[i]=rnd.Gauss(mu, sigma2) ;
  return array;
}


double* meanvalues(double *r ,int N, int L, double *ave, double *ave2) { /* Function per il calcolo di <r> e <r^2> per ogni blocco:
									 input-> array, #blocchi , #lanci per blocco
									 output -> array di r_i e r_i^2*/
  for(int i=0; i<N ; i++) {	   //Fisso un blocco   
    double sum=0  ;		   
    for(int j=0; j<L; j++) {	   //Sommo tutti i valori ottenuti nel blocco i 
      int k = j+i*L  ;
      sum += r[k] ;
    }
    ave[i] = sum/L ;         //r_i 
    ave2[i] = ave[i]*ave[i] ;  //(r_i)^2 
  }
  return ave, ave2 ;
}


void cumulativeaverage(double* avecall, double* ave2call, double* aveput, double* ave2put, int N, const char* nomefile) { 
							 /* Function per il calcolo delle medie cumulate progressive e la scrittura su file:
							 input -> arrays di r_i e r_i^2, #blocchi, nome file su cui scrivere */
  ofstream myfile;
  myfile.open(nomefile);
  double sum_progcall[N]={0}, sum2_progcall[N]={0}, err_progcall[N]={0} ;
  double sum_progput[N]={0} , sum2_progput[N]={0}, err_progput[N]={0} ;

  for(int i=0; i<N ; i++) {    	/*Ciclo per calcolare le medie cumulate progressive: 
				si fissa il blocco i e si va a sommare tutte le stime ottenute sino a quel punto;
				infine si scrivono su file i risultati ottenuti*/
    for(int j=0; j<i+1 ; j++) { 
      sum_progcall[i] += avecall[j] ;      //SUM_{j=0,i} r_j
      sum2_progcall[i] += ave2call[j]  ;   //SUM_{j=0,i} (r_j)^2

      sum_progput[i] += aveput[j] ;      
      sum2_progput[i] += ave2put[j]  ;
    }
    sum_progcall[i]/=(i+1) ; 		//Media cumulativa
    sum2_progcall[i]/=(i+1) ;  		//Media cumulativa dei quadrati
    err_progcall[i] = errore(sum_progcall, sum2_progcall, i) ; 	//Incertezza statistica

    sum_progput[i]/=(i+1) ; 		
    sum2_progput[i]/=(i+1) ;  		
    err_progput[i] = errore(sum_progput, sum2_progput, i) ; 
    
    myfile << i << " " << sum_progcall[i] << " " << err_progcall[i] << " " << sum_progput[i] << " " << err_progput[i] <<endl;
  }
  myfile.close() ;  
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

  int M=10000 ;  		//Numero di simulazioni da generare
  int N= 100; 			//Numero blocchi
  double S0 = 100 ;		//Costo iniziale dell'asset
  double T = 1 ;		//Delivery time
  double K = 100 ;		//Strike price
  double r= 0.1 ;		//Tasso di interesse risk-free
  double sigma = 0.25 ;		//Volatilità
  double Z[M]={0} ;		//Array per variabili casuali gaussiane 
  double S_T[M]={0} ; 		//Array per il valore finale dello spot price

  double C_call[M]={0}, Cc_medio[N]={0}, Cc_medio2[N]={0} ; 	//Arrays per il costo dell'opzione call
  double C_put[M]={0}, Cp_medio[N]={0}, Cp_medio2[N]={0}; 	//Arrays per il costo dell'opzione put

//Primo metodo: valuto direttamente S(T)
  boxmuller(Z, M, 0, 1);     //Genero numeri casuali N(0,1)

  for(int i=0; i<M ; i++) {
    S_T[i]= S0*exp( (r- sigma*sigma/2)*T + sigma* Z[i]* sqrt(T) ) ;
    C_call[i]= exp(-r*T)*max(0.,S_T[i]-K) ;
    C_put[i]= exp(-r*T)*max(0.,K-S_T[i]) ;
  }
   
  meanvalues(C_call, N, M/N, Cc_medio, Cc_medio2);
  meanvalues(C_put, N, M/N, Cp_medio, Cp_medio2);
  cumulativeaverage(Cc_medio, Cc_medio2, Cp_medio, Cp_medio2, N, "Metodo1.dat");

//Secondo metodo: campiono il cammino S(T)
  Z[100]={0} ;

  for(int i=0; i<M ; i++) {
    boxmuller(Z, 100, 0, 1);     //Genero numeri casuali N(0,1)      
    double S_j=S0 ;
    for(int j=0 ; j<100 ; j++) S_j= S_j*exp( (r-sigma*sigma/2)/100. + sigma* Z[j]/10. ) ;
    
    S_T[i]= S_j ;      
    C_call[i]= exp(-r*T)*max(0.,S_T[i]-K) ;
    C_put[i]= exp(-r*T)*max(0.,K-S_T[i]) ;
  }

  meanvalues(C_call, N, M/N, Cc_medio, Cc_medio2);
  meanvalues(C_put, N, M/N, Cp_medio, Cp_medio2);
  cumulativeaverage(Cc_medio, Cc_medio2, Cp_medio, Cp_medio2, N, "Metodo2.dat");

  return 0;

}
	
	
