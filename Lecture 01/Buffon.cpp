#include <iostream>
#include <fstream>
#include <string>
#define _USE_MATH_DEFINES
#include <cmath>
#include <array>
#include "random.h"

using namespace std;
Random rnd;

/*-----------------------------------------------------------Functions-----------------------------------------------------------*/

double errore(double *AV, double *AV2,int n) {  /*Function per la stima dell'incertezza statistica: 
						input-> <r> e <r^2> calcolati al blocco n*/
  if(n==0) return 0 ;
  else return sqrt((AV2[n] - AV[n]*AV[n])/n) ;
}


double* uniformgenerator(double *r_u, int M) { /*Function per generare M numeri casuali y uniformi in [0,1]
                                               input-> array per numeri casuali, numero lanci*/
 
  for(int i=0; i<M ; i++) r_u[i]= rnd.Rannyu();
  return r_u ;
}


double* acceptreject(double *theta, int M) {/*Function per generare M numeri casuali theta uniformi in [0, pi] con il metodo accept-reject
                                           input-> array per numeri casuali, numero lanci*/
 
  int j=0;
  while(theta[M-1]==0) {  
    double x=2*rnd.Rannyu()-1;
    double y=rnd.Rannyu();
    if(x*x+y*y<1) {
      theta[j]=acos(x/sqrt(x*x+y*y) );
      j+=1 ;
    }
  }  

  return theta ;
}


int* num_hits (double* y, double* theta, int N, int L, double d, double l, int* num_intersez) {/*Function per il conteggio delle intersezioni
                                                                                                per ogni blocco: 
												input-> array di y e theta, #blocchi, #lanci
												        per blocco, d, l, array conteggi*/
   
  for(int i=0; i<N; i++) {  //Fisso il blocco
    for(int j=0 ; j<L ; j++) { //Valuto tutte le misure nel blocco
      int k = j+i*L  ;
      for(int s=0; s<=30; s++){
        if(cos(theta[k]) >=0 && y[k]> s*d && y[k]-l*cos(theta[k]) < s*d ) num_intersez[i]+=1 ;
        else if(cos(theta[k]) <0 && y[k]< s*d && y[k]-l*cos(theta[k]) > s*d )  num_intersez[i]+=1 ;
      }
    }
  }
  return num_intersez ; 
}


double* meanvalues(int *counter,int N, int L, double l, double d, double *ave, double *ave2) { /* Function per il calcolo di <pi> e <pi^2>
                                                                                                per ogni blocco:
							  					input-> array del #intersezioni, #blocchi
									 			output -> array di r_i e r_i^2*/

  for(int i=0; i<N ; i++) {   
    ave[i] = 2.*l*L/(counter[i]*d) ;         //r_i 
    ave2[i] = ave[i]*ave[i] ; 		  //(r_i)^2 
  }
  return ave, ave2 ;
}


void cumulativeaverage(double* ave, double* ave2, int N) {  /* Function per il calcolo delle medie cumulate progressive:
		  					    input -> array di r_i e r_i^2, numero blocchi */
  const double pi= M_PI ;

  ofstream myfile;
  myfile.open ("Pigreco.dat");
  double sum_prog[N]={0}, sum2_prog[N]={0}, err_prog[N]={0} ;

  for(int i=0; i<N ; i++) {    	/*Ciclo per calcolare le medie cumulate progressive: 
				si fissa il blocco i e si va a sommare tutte le stime ottenute sino a quel punto;
				infine si scrivono su file i risultati ottenuti*/
    for(int j=0; j<i+1 ; j++) { 
      sum_prog[i] += ave[j] ;      //SUM_{j=0,i} r_j
      sum2_prog[i] += ave2[j]  ;   //SUM_{j=0,i} (r_j)^2
    }

    sum_prog[i]/=(i+1) ; 		//Media cumulativa
    sum2_prog[i]/=(i+1) ;  		//Media cumulativa dei quadrati
    err_prog[i] = errore(sum_prog, sum2_prog, i) ; 	//Incertezza statistica

    
    myfile << i << " " << sum_prog[i] -pi << " " << err_prog[i] << endl;
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

  double d=1./30 ;
  double l=1./40 ;
  double y[M]={0} ;
  double theta[M]={0} ;
  int num_intersez[N]={0} ;
  double pivalue[N]={0}, pivalue_2[N]={0};

  uniformgenerator(y,M);
  acceptreject(theta, M);
  num_hits (y, theta, N, L, d, l, num_intersez);
  meanvalues(num_intersez ,N, L, l, d, pivalue, pivalue_2) ;
  cumulativeaverage(pivalue, pivalue_2, N);

  return 0;

}
	
