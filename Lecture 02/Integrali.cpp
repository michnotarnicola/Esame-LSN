#include <iostream>
#include <fstream>
#include <string>
#define _USE_MATH_DEFINES
#include <cmath>
#include <array>
#include <random>
#include <vector>
#include "random.h"

using namespace std;
Random rnd;

/*-----------------------------------------------------------Functions-----------------------------------------------------------*/

double errore(double *AV, double *AV2,int n) {  /*Function per la stima dell'incertezza statistica: 
						input-> <r> e <r^2> calcolati al blocco n*/
  if(n==0) return 0 ;
  else return sqrt((AV2[n] - AV[n]*AV[n])/n) ;
}


double* uniformgenerator(double *r_u, int M) { /*Function per generare M numeri casuali uniformi in [0,1]
                                               input-> array per numeri casuali, numero lanci*/
 
  for(int i=0; i<M ; i++) r_u[i]=rnd.Rannyu();
  return r_u ;
}


double* acceptreject(double *d_x, int M) { /*Function per generare M numeri distributi secondo d_x con il metodo accept-reject
                                           input-> array per numeri casuali, numero lanci*/
 
  int j=0;
  while(d_x[M-1]==0) {  
    double x=rnd.Rannyu();
    double r=rnd.Rannyu();
    if(r<(1.-x*x) && x!=0) {
      d_x[j]=x;
      j+=1 ;
    }
  }  
  return d_x ;
}


double* integralvalue_u(double *r ,int N, int L, double *ave, double *ave2) { /* Function per il calcolo dell'integrale per ogni blocco
									      (campionando una distruzione uniforme):
									      input-> array di numeri casuali, #blocchi , #lanci per blocco
									      output -> array di I_i e I_i^2*/
  const double pi= M_PI ;

  for(int i=0; i<N ; i++) {   		//Fisso un blocco
    double integ=0  ;
    for(int j=0; j<L; j++) {		//Sommo tutti i valori associati al blocco i
      int k = j+i*L  ;
      integ += cos(pi*r[k]/2) ;
      ave[i] = integ*pi/(2*L) ;  	//I_i 
      ave2[i] = ave[i]*ave[i] ;  	//(I_i)^2 
    }	
  }
  return ave, ave2 ;
}


double* integralvalue_is(double *r ,int N, int L, double *ave, double *ave2) { /* Function per il calcolo dell'integrale per ogni blocco 
									       (con l'importance sampling):
									       input-> array di numeri casuali, #blocchi , #lanci per blocco
									       output -> array di I_i e I_i^2*/
  const double pi= M_PI ;
  for(int i=0; i<N ; i++) {   
    double integ=0  ;
    for(int j=0; j<L; j++) {
      int k = j+i*L  ;
      integ += cos(pi*r[k]/2)/(1-r[k]*r[k]) ;
      ave[i] = integ*pi/(3*L) ;    
      ave2[i] = ave[i]*ave[i] ;   
    }
  }
  return ave, ave2 ;
}


void cumulativeaverage (double* integ_u, double* integ2_u, double* integ_d, double* integ2_d, int N) {/* Function per il calcolo delle medie 													      cumulate progressive:
							  				       	      input -> array di r_i e r_i^2 per i due
												      metodi, numero blocchi */
  double sum_prog_u[N]={0} ;
  double sum2_prog_u[N]={0} ;
  double err_prog_u[N]={0} ;
  double sum_prog_d[N]={0} ;
  double sum2_prog_d[N]={0} ;
  double err_prog_d[N]={0} ;

  ofstream myfile;
  myfile.open ("Integrale.dat");

  for(int i=0; i<N ; i++) {    	/*Ciclo per calcolare le medie cumulate progressive: 
				si fissa il blocco i e si va a sommare tutte le stime ottenute sino a quel punto;
				infine si scrivono su file i risultati ottenuti*/
    for(int j=0; j<i+1 ; j++) { 
      sum_prog_u[i] += integ_u[j] ;    //SUM_{j=0,i} I_j
      sum2_prog_u[i] += integ2_u[j]  ;  //SUM_{j=0,i} (I_j)^2
      sum_prog_d[i] += integ_d[j] ;    
      sum2_prog_d[i] += integ2_d[j]  ;  
    }

    sum_prog_u[i]/=(i+1) ; 
    sum2_prog_u[i]/=(i+1) ;  
    err_prog_u[i] = errore(sum_prog_u, sum2_prog_u, i) ; 

    sum_prog_d[i]/=(i+1) ; 
    sum2_prog_d[i]/=(i+1) ; 
    err_prog_d[i] = errore(sum_prog_d, sum2_prog_d, i) ; 
    
    myfile << i << " " << sum_prog_u[i]-1 << " " << err_prog_u[i] << " " << sum_prog_d[i]-1 << " " << err_prog_d[i] << endl;
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

  int M=100000 ;             //Numero totale di lanci
  int N=100    ;             //Numero di blocchi
  int L= M/N   ;             //Numero di lanci in ogni blocco

  double unif[M]={0}  ;      //Arrays per numeri casuali 
  double d_x[M]={0}  ;
  double integ_u[N]={0}, integ2_u[N]={0},integ_d[N]={0} ,integ2_d[N]={0}  ;    //Arrays per le stime di I , I^2 per le 2 tecniche

  uniformgenerator(unif, M) ;
  acceptreject(d_x, M) ;
  integralvalue_u(unif,N,L, integ_u, integ2_u);
  integralvalue_is(d_x,N,L, integ_d, integ2_d);
  cumulativeaverage (integ_u, integ2_u, integ_d, integ2_d, N) ;

  return 0;

}


