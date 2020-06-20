#include <iostream>
#include <fstream>
#include <ostream>
#define _USE_MATH_DEFINES
#include <cmath>
#include <iomanip>
#include <algorithm>
#include <array>
#include <vector>
#include <random>
#include "random.h"

using namespace std;

/*-------------------------------------------------------Global Variables-------------------------------------------------------*/
Random rnd;
const int N_cromo=1000, N_cities=32, N_step=N_cromo/5, N_gen=500;
int bestcromo1, bestcromo2, worstcromo1, worstcromo2, finalcromo;

std::vector<std::array<int, N_cities> > Population; 
double xPosition[N_cities], yPosition[N_cities];
int Ordinamento[N_cromo]={0};


/*----------------------------------------------------------Functions-----------------------------------------------------------*/

double Costfunc (std::array<int,N_cities> cromo) { /*Function per il calcolo della funzione costo dato un cromosoma*/
  double L=0;
  int i, iplus1;
  for(int icity=0; icity< N_cities ; icity ++) {
    i= cromo[icity];				//Leggo la sequenza delle città in cromo e calcolo L
    if(icity==N_cities-1) iplus1= cromo[0];	//PBC
    else iplus1= cromo[icity+1];
    L+= pow( (xPosition[i] - xPosition[iplus1]) ,2) + pow( (yPosition[i] - yPosition[iplus1]) ,2);
  }
  return L;
}


bool Check (std::array<int,N_cities> cromo) {  /*Function di controllo che i vincoli sul cromosoma siano rispettati*/
  bool primoelemento=true, ordine=true, citynum=true;
  if(cromo[0]!=0) primoelemento=false;
  for(int icity=0; icity<N_cities; icity++) {
    if(cromo[icity]>=N_cities) citynum=false ;
    for(int jcity=icity+1; jcity<N_cities; jcity++) if(cromo[icity]== cromo[jcity]) ordine=false ;
  }
  return primoelemento && ordine && citynum;
}


std::array<int,N_cities>&  Copy (std::array<int,N_cities> &cromo1, std::array<int,N_cities> &cromo2) {
											/*Function per copiare un cromosoma su un altro*/
  for(int icity=0; icity< N_cities ; icity++) cromo2[icity]= cromo1[icity];
  return cromo1, cromo2;
}


double Averages (double* L) {		/*Function per il calcolo di <L> sulla miglior metà di popolazione*/
  double sum=0;
  for(int icromo=0; icromo< N_cromo/2; icromo++) sum+= L[icromo];
  return sum*2/N_cromo;
}


void Ordering (int igen) {	/*Funzione per riordinare gli individui per fitness crescente e scrivere i valori di L su file*/
  double L[N_cromo];
  double temp;
  int index_min, index_temp;

  for(int icromo=0; icromo<N_cromo; icromo++) {
    L[icromo]=Costfunc(Population[icromo]);	//Calcolo il fitness su tutte le configurazioni
    Ordinamento[icromo]= icromo;		//Inizializzo la sequenza di ordinamento dei cromosomi
  }
 
  for(int i=0; i<N_cromo-1; i++) {		//Selection-sort
    index_min=i;
    for(int j=i+1; j<N_cromo; j++) { 
      if(L[j]<L[i]) {
        index_min=j;
        temp=L[index_min];
        index_temp= Ordinamento[index_min];
        L[index_min]=L[i];
        L[i]=temp;
        Ordinamento[index_min]=Ordinamento[i];
        Ordinamento[i]=index_temp;
      }
    }
  }
 
//Scrittura su file della cost function della miglior stringa e della cost function media
  ofstream bestcost, averagecost;
  bestcost.open("Best_L.dat",ios::app);
  averagecost.open("Ave_L.dat",ios::app);
  bestcost << igen << setw(14) << L[0] << endl;	//L[0] è il valore minimo degli L
  averagecost << igen << setw(14) << Averages(L) << endl;
  bestcost.close();
  averagecost.close();
}


/*------Selection Operator------*/
std::array<int,N_cities>& Selection (std::array<int,N_cities> &Genitore1, std::array<int,N_cities> &Genitore2) {
  double p=5, q=1./p;

  //Selezione dei migliori due individui per la riproduzione
  finalcromo= Ordinamento[0];			//Utile per la funzione BestPath
  bestcromo1= Ordinamento[ int(N_cromo* pow(rnd.Rannyu(),p)) ];
  bestcromo2= Ordinamento[ int(N_cromo* pow(rnd.Rannyu(),p)) ];
  while(bestcromo2==bestcromo1) bestcromo2= Ordinamento[ int(N_cromo* pow(rnd.Rannyu(),p)) ];

  Copy(Population[bestcromo1], Genitore1);
  Copy(Population[bestcromo2], Genitore2);

  //Selezione dei peggiori due individui
  worstcromo1= Ordinamento[ int(N_cromo* pow(rnd.Rannyu(),q)) ];
  worstcromo2= Ordinamento[ int(N_cromo* pow(rnd.Rannyu(),q)) ];
  while(worstcromo2==worstcromo1) worstcromo2= Ordinamento[ int(N_cromo* pow(rnd.Rannyu(),q)) ];

  return Genitore1, Genitore2;
}


/*------Mutation Operators------*/
std::array<int,N_cities> PairPerm (std::array<int,N_cities> cromo) {

  int gene1= rand()%(N_cities-1) +1;		//Scelgo due geni distinti da 2,....,N
  int gene2= rand()%(N_cities-1) +1;
  while(gene2==gene1) gene2= rand()%(N_cities-1) +1;
  int temp=cromo[gene1];
  cromo[gene1]=cromo[gene2];
  cromo[gene2]=temp;
  if(Check(cromo)==false) cout << "Errore nella mutazione: PairPerm!" <<endl;
  return cromo;
}


std::array<int,N_cities>& Shift (std::array<int,N_cities>& cromo) {

  int n_shift= rand()%(N_cities-1);		//Numero di traslazioni in avanti (si trasla a partire dalla seconda città)
  int index=0;
  std::array<int,N_cities> cromo_copy;
  Copy(cromo, cromo_copy);

  for(int icity=n_shift+1; icity<N_cities; icity++) {
    index+=1;
    cromo[icity]= cromo_copy[index];
  }
  for(int icity=1; icity<=n_shift; icity++) {
    index+=1;    
    cromo[icity]= cromo_copy[index];
  }
  if(Check(cromo)==false) cout << "Errore nella mutazione: Shift!" <<endl;
  return cromo;
}


std::array<int,N_cities>& Permutation (std::array<int,N_cities>& cromo) {
  
  int n_perm= rand()%((N_cities-1)/2); /*Lunghezza blocchi da permutare (il primo blocco da dalla seconda città alla n_perm -esima,
					  il secondo dalla N_cities-n_perm alla n_perm)*/
  int index=0;
  std::array<int,N_cities> cromo_copy;
  Copy(cromo, cromo_copy);

  for(int icity=1; icity<=n_perm+1; icity++) {
    cromo[icity]= cromo_copy[N_cities-1-n_perm +index];
    cromo[N_cities-1-n_perm +index]= cromo_copy[icity];
    index+=1;
  }
  for(int icity=n_perm+2; icity<N_cities-n_perm-1; icity++) cromo[icity]= cromo_copy[icity];
  if(Check(cromo)==false) cout << "Errore nella mutazione: Permutation!" <<endl;
  return cromo;
}


std::array<int,N_cities>& Inversion (std::array<int,N_cities>& cromo) {
  
  int n_inv= rand()%(N_cities-1);		//Si inverte un blocco a partire dalla seconda città
  int index=0;
  std::array<int,N_cities> cromo_copy;
  Copy(cromo, cromo_copy);

  for(int icity=1; icity<=n_inv+1; icity++) {
    cromo[icity]= cromo_copy[n_inv+1-icity+1];
    index+=1;
  }
  if(Check(cromo)==false) cout << "Errore nella mutazione: Inversion!" <<endl;
  return cromo;
}


/*------Crossover Operator------*/
std::array<int,N_cities>& Crossover (std::array<int,N_cities> &cromo1, std::array<int,N_cities> &cromo2) {
  
  int cut_1= rand()%(N_cities-1)+1; 	//Posizione dello start e dello stop del cu2
  int cut_2= rand()%(N_cities-1)+1;
  while(cut_1==cut_2) cut_2= rand()%(N_cities-1)+1 ;
  double temp;

  for(int icity=min(cut_1, cut_2); icity<= max(cut_1,cut_2); icity++) {
    temp= cromo1[icity];
    cromo1[icity]= cromo2[icity];
    cromo2[icity]= temp;
  }
  return cromo1, cromo2 ;
}

/*------------------------------*/
void Input () {
  const double pi=M_PI;
  //Scelgo a caso le posizioni delle città e le scrivo su file
  ofstream Citta;
  Citta.open("Citta.dat");
  for(int icity=0; icity<N_cities; icity++) {
//    double theta=2*pi*rnd.Rannyu();				//Caso 1: circonferenza di R=1
//    xPosition[icity]= cos(theta);
//    yPosition[icity]= sin(theta);
    xPosition[icity]= rnd.Rannyu();				//Caso 2: quadrato
    yPosition[icity]= rnd.Rannyu();
    Citta << icity << setw(14) << xPosition[icity] << setw(14) << yPosition[icity] << endl;
  }
  Citta.close();

  //Costruisco la popolazione iniziale
  std::array<int, N_cities> cromo1={0};
  for(int icity=0; icity<N_cities; icity++) cromo1[icity]=icity;
  Population.push_back(cromo1);
  for(int icromo=1; icromo<N_cromo; icromo++) {
    for(int ip=0; ip<100*rnd.Rannyu() ; ip++) PairPerm(cromo1);
    for(int ip=0; ip<10*rnd.Rannyu() ; ip++) Shift(cromo1);
    for(int ip=0; ip<10*rnd.Rannyu() ; ip++) Permutation(cromo1);
    for(int ip=0; ip<10*rnd.Rannyu() ; ip++) Inversion(cromo1);
    Population.push_back(cromo1);
    if(Check(cromo1)==false) cout << "Errore di inizializzazione!" << endl;
  }
  Ordering(0);			//Calcolo la cost function per la generazione igen=0, quella iniziale

}


void Move(int igen) {

  std::array<int,N_cities> Genitore1={0}, Genitore2={0}, Figlio1={0}, Figlio2={0};
  double p_mutaz=0.4, p_cross=0.8;

  for(int istep=1; istep<= N_step; istep++) {
    //Selezione
    Selection(Genitore1, Genitore2);

    //Mutazione  
    if(rnd.Rannyu()<=p_mutaz) PairPerm(Genitore1);
    if(rnd.Rannyu()<=p_mutaz) PairPerm(Genitore2);
    if(rnd.Rannyu()<=p_mutaz) Shift(Genitore1);
    if(rnd.Rannyu()<=p_mutaz) Shift(Genitore2);
    if(rnd.Rannyu()<=p_mutaz) Permutation(Genitore1);
    if(rnd.Rannyu()<=p_mutaz) Permutation(Genitore2);
    if(rnd.Rannyu()<=p_mutaz) Inversion(Genitore1);
    if(rnd.Rannyu()<=p_mutaz) Inversion(Genitore2);

    //Crossover
    Copy(Genitore1, Figlio1);					//Per il momento i figli sono uguali ai genitori
    Copy(Genitore2, Figlio2);
    if(rnd.Rannyu()<=p_cross) {
      Crossover(Figlio1, Figlio2);				//Se faccio il cross modifico i figli
      while(Check(Figlio1)==false && Check(Figlio2)==false) {
        Copy(Genitore1, Figlio1);				//Se il cross non è valido rimetto i figli uguali ai genitori e ci riprovo
        Copy(Genitore2, Figlio2);
        Crossover(Figlio1, Figlio2);
      }
    }

    //Modifico la popolazione finale rimuovendo i due peggiori individui
    Copy(Figlio1, Population[worstcromo1]) ;
    Copy(Figlio2, Population[worstcromo2]) ;
  }
  Ordering(igen);
}


void BestPath () {	/*Function per la scrittura finale del miglior percorso*/
  ofstream bestpath;
  bestpath.open("BestPath.dat");
  for (int icity=0; icity< N_cities ; icity++) bestpath << Population[finalcromo][icity] << endl;
  bestpath.close();
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
  
  Input();
  for(int igen=1; igen<= N_gen; igen++) {
    cout << "Generazione " << igen << endl;
    Move(igen);
    cout << "------------------------------" << endl;
  }
  BestPath();

  return 0;

}
