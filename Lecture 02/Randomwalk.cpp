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

int* discreteRW (int *x, int *y, int *z, int M) {  /*Function per il RW discreto, lancio un dado a 6 facce e a seconda del risultato mi muovo
						   avanti/indietro nelle 3 direzioni:
						   input -> arrays per coordinate, numero lanci */

  for(int i=0; i<M; i++) {
    int num=rand()%6 ;
    if(num==0) x[i]+=1 ;
    else if(num==1) x[i]-=1 ;
    else if(num==2) y[i]+=1 ;
    else if(num==3) y[i]-=1 ;
    else if(num==4) z[i]+=1 ;
    else z[i]-=1 ;
  }
  return x, y, z ;
} 


double* contRW (double *x, double *y, double *z, int M) { /*Function per il RW continuo, genero angoli casual theta e phi:
						          input -> arrays per coordinate, numero lanci */

  const double pi= M_PI ;
  for(int i=0; i<M; i++) {
    double phi= 2*pi* rnd.Rannyu() ;
    double theta= acos(1.-2* rnd.Rannyu() ) ;

    x[i] += sin(theta)*cos(phi) ;
    y[i] += sin(theta)*sin(phi) ;
    z[i] += cos(theta) ;
  }
  return x, y, z ;
} 


double* meanvalues(int *x ,int *y, int *z, int N, int L, double *ave, double *ave2) { /* Function per il calcolo di ||r||**2 per ogni
										      blocco (caso discreto):
									 	      input-> arrays coordinate, #blocchi , #lanci per blocco
									 	      output -> array di r_i e r_i^2*/
  for(int i=0; i<N ; i++) {   
    double norm=0  ;
    for(int j=0; j<L; j++) {
      int k = j+i*L  ;
      norm += x[k]*x[k] + y[k]*y[k] + z[k]*z[k] ;
    }
    ave[i] = norm/L ;         //r_i 
    ave2[i] = ave[i]*ave[i] ;  //(r_i)^2 
  }
  return ave, ave2 ;
}

double* meanvalues_cont (double *x ,double *y, double *z, int N, int L, double *ave, double *ave2) { /* Function per il calcolo di ||r||**2 
												     per ogni blocco (caso continuo):
									 			     input-> arrays coordinate, #blocchi ,
													     #lanci per blocco
									 			     output -> array di r_i e r_i^2*/
  for(int i=0; i<N ; i++) {   
    double norma=0.0  ;
    for(int j=0; j<L; j++) {
      int k = j+i*L  ;
      norma += x[k]*x[k] + y[k]*y[k] + z[k]*z[k] ;
    }
    ave[i] = norma/L ;         //r_i 
    ave2[i] = ave[i]*ave[i] ;  //(r_i)^2 
  }
  return ave, ave2 ;
}


void cumulativeaverage(double* ave, double* ave2, int N, std::ofstream& myfile) { 
							 /* Function per il calcolo delle medie dei valor medi e la scrittura su file:
							 input -> arrays di r_i e r_i^2, #blocchi, nome ofstream */
  double sum=0 ;
  double sum2=0 ;
  double err=0 ;

  for(int j=0; j<N ; j++) { 
    sum += ave[j] ;      //SUM_{j=1,N} r_j
    sum2 += ave2[j]  ;   //SUM_{j=1,N} (r_j)^2
  }
  sum/=N ; 		//Media cumulativa
  sum2/=N ;  		//Media cumulativa dei quadrati
  err = sqrt((sum2 - sum*sum)/(N-1)) ; 	//Incertezza statistica
  myfile << sqrt(sum) << " " << 1./(2*sqrt(sum)) * err << endl;
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

  int M=100000 ;            //Numero totale di esperimenti
  int N=100;                //Numero blocchi
  int L=M/N;                //Numero lanci nel blocco
  int Nstep=100    ;        //Numero di passi per esperimento


//Caso 1 : Random Walk discreto

  int x[M]={0}, y[M]={0}, z[M]={0};	     //Posizione del walker in tutti i lanci, fissata nell'origine
  double rmedio[N]={0}, rmedio2[N]={0};

  ofstream Discreto;
  Discreto.open ("RW_discreto.dat");
  Discreto << 0 << " " << 0 << endl; //scrivo su file la distanza iniziale (l'origine) e l'incertezza su essa
  for(int i=0; i<Nstep; i++) { //fisso lo step i=0, ..., 99 (linguaggio C++ per 1, ..., 100)
    discreteRW(x,y,z,M) ;
    meanvalues(x,y,z ,N, L, rmedio, rmedio2) ;
    cumulativeaverage(rmedio, rmedio2, N, Discreto) ;
  }
  Discreto.close() ;


//Caso 2 : Random Walk continuo
  double x_cont[M]={0}, y_cont[M]={0}, z_cont[M]={0};	     
  rmedio[N]={0}, rmedio2[N]={0};
   
  ofstream Continuo;
  Continuo.open ("RW_continuo.dat");
  Continuo << 0 << " " << 0 << endl; //scrivo su file la distanza iniziale (l'origine) e l'incertezza su essa
  for(int i=0; i<Nstep; i++) { //fisso lo step
    contRW(x_cont, y_cont, z_cont, M) ;
    meanvalues_cont(x_cont,y_cont,z_cont ,N, L, rmedio, rmedio2) ;
    cumulativeaverage(rmedio, rmedio2, N, Continuo) ;
  }
  Continuo.close() ;

  return 0;

}
