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

double errore(double AV, double AV2,int n) {  /*Function per la stima dell'incertezza statistica: 
						input-> <r> e <r^2> calcolati al blocco n*/
  if(n==0) return 0 ;
  else return sqrt((AV2 - AV*AV)/n) ;
}


double* randomgenerator(double *r, double *var, int M) { /* Function per generare numeri casuali uniformi in [0,1]:
							 input -> array per numeri casuali, array per varianza, #lanci */
   
  for(int i=0; i<M ; i++) { 
    r[i]=rnd.Rannyu();
    var[i]= (r[i]-0.5)*(r[i]-0.5) ;
  }
  return r , var;
}


double* meanvalues(double *r ,int N, int L, double *ave, double *ave2) { /* Function per il calcolo di <r> e <r^2> per ogni blocco:
									 input-> array di numeri casuali, #blocchi , #lanci per blocco
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


void cumulativeaverage(double* ave,double* ave2, double* var_ave, double* var_ave2,int N) {/* Function per il calcolo delle medie cumulate
											       progressive:
							  				       input -> array di r_i e r_i^2 per <r> e
											        var(r), numero blocchi */

  ofstream myfile;
  myfile.open ("Dati.dat");
  double sum_prog[N]={0} ;
  double sum2_prog[N]={0} ;
  double err_prog[N]={0} ;
  double var_sum_prog[N]={0} ;
  double var_sum2_prog[N]={0} ;
  double var_err_prog[N]={0} ;

  for(int i=0; i<N ; i++) {    	/*Ciclo per calcolare le medie cumulate progressive: 
				si fissa il blocco i e si va a sommare tutte le stime ottenute sino a quel punto;
					infine si scrivono su file i risultati ottenuti*/
    for(int j=0; j<i+1 ; j++) { 
      sum_prog[i] += ave[j] ;      //SUM_{j=0,i} r_j
      sum2_prog[i] += ave2[j]  ;   //SUM_{j=0,i} (r_j)^2

      var_sum_prog[i] += var_ave[j] ;   
      var_sum2_prog[i] += var_ave2[j]  ;  
    }

    sum_prog[i]/=(i+1) ; 		//Media cumulativa
    sum2_prog[i]/=(i+1) ;  		//Media cumulativa dei quadrati
    err_prog[i] = errore(sum_prog[i], sum2_prog[i], i) ; 	//Incertezza statistica

    var_sum_prog[i]/=(i+1) ;
    var_sum2_prog[i]/=(i+1) ; 
    var_err_prog[i] = errore(var_sum_prog[i], var_sum2_prog[i], i) ; 
    
    myfile << i << " " << sum_prog[i]-0.5 << " " << err_prog[i] << " " << var_sum_prog[i]-1/12. << " " << var_err_prog[i] << endl;
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

  int M=100000 ;           	  //Numero totale di lanci
  int N=100    ;                 //Numero di blocchi
  int L= M/N   ;                 //Numero di lanci in ogni blocco

  double r[M]={0}, var[M]={0} ;		  //Arrays per numeri casuali
  double ave[N]={0}, ave2[N]={0}, var_ave[N]={0}, var_ave2[N]={0} ;		  //Arrays per le stime di r , r^2

  randomgenerator(r,var, M);
  meanvalues(r ,N, L, ave, ave2) ;
  meanvalues(var ,N, L, var_ave, var_ave2) ;
  cumulativeaverage(ave, ave2, var_ave, var_ave2, N);

  return 0;

}
	
	
