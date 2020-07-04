#include <cmath>
#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <string>
#include "funzioni.h"

using namespace std;

void setrandom(Random &rnd) {
   int seed[4];
   int p1, p2;
   ifstream Primes("Primes");
   if (Primes.is_open()){
      Primes >> p1 >> p2 ;
   } else cerr << "PROBLEM: Unable to open Primes" << endl;
   Primes.close();

   ifstream input("seed.in");
   string property;
   if (input.is_open()){
      while ( !input.eof() ){
         input >> property;
         if( property == "RANDOMSEED" ){
            input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
            rnd.SetRandom(seed,p1,p2);
         }
      }
      input.close();
   } else cerr << "PROBLEM: Unable to open seed.in" << endl;
}




void calcolacumulate(double *media, double *media_2, double *media_cumulativa, double *media_2_cumulativa, double *errore, int N) {
	for (int i=0; i<N; i++) { //ciclo sui blocchi
		for (int j=0; j<i+1; j++) { //ciclo che mi permette di sommare fino all'i-esimo blocco (sto cumulando) 
			media_cumulativa[i]+=media[j];
			media_2_cumulativa[i]+=media_2[j];
			}
		media_cumulativa[i]=media_cumulativa[i]/(i+1); //media delle medie: devo quindi dividere per il numero di componenti  
		media_2_cumulativa[i]=media_2_cumulativa[i]/(i+1);
		errore[i]=calcolaerrore(media_cumulativa,media_2_cumulativa,i); //associa incertezza a ciascuna media cumulata
		}
	}


double calcolaerrore (double *x, double *y, int i) { //x: gioca il ruolo di media_cumulata; y: media_2_cumulata; i: numero blocchi sommati
	if (i==0) {
		return 0; //se faccio solo una misura (un singolo blocco) non ha senso definire un'incertezza
		}
	return sqrt((y[i]-pow(x[i],2))/i); //formula per l'incertezza
	}



void azzeraarray(double *x, int n) {
	for (int i=0; i<n; i++) {
		x[i]=0.;
		}
	}


void calcolaprezzo(Random &rnd, const double S0,const double t0,const double T,const double K,const double r,const double sigma,double *media, double *media_2,int N, int L, int intervalli,string s) {
	if (s!="call"&&s!="put") { //controlla che la stringa inserita corrisponda ad un'opzione
		cout<<"Errore!"<<endl;
		exit(1);
		}

	int segno;
 
//tra le due opzioni cambia soltanto che nella call serve il massimo tra 0 e S(T)-K; nella put il massimo tra 0 e K-S(T). Computazionalmente differiscono solo per un segno meno 	

	if (s=="call") { 
		segno=1;
		}
	
	if (s=="put") {
		segno=-1;
		}


	for (int i=0; i<N; i++) { //ciclo sul numero dei blocchi
		double somma=0;
		for (int j=0; j<L; j++) { //ciclo sul numero di passi in ogni blocco (nota che il loro prodotto restituisce M)
			double Sprec=S0;
			double ST=S0;
			double dt=(double)(T-t0)/((double) intervalli); //divido in intervalli identici, quindi t_(i+1)-t_i=dt per ogni t
			for (int k=0; k<intervalli; k++) { //formula ricorsiva per trovare S(T)	
				ST=Sprec*exp((r-0.5*sigma*sigma)*dt+sigma*rnd.Gauss(0.,1.)*sqrt(dt));
				Sprec=ST;
				}
			somma=somma+exp(-r*T)*fmax(0.,(segno)*(ST-K)); //formula per la singola stima del costo (l'unica differenza tra call e put è determinata da quel segno) 
			}
		media[i]=somma/L; //def media (per ogni blocco ho una stima del costo dell'opzione)
		media_2[i]=pow(media[i],2); //serve la media di ciascun blocco e il suo valore al quadrato
		}
	}

void stamparisultati(ofstream &out,double *media_cumulativa,double *errore,int L, int N) {
	for (int i=0; i<N; i++) {
		out<<(i+1)*L<<" "<<media_cumulativa[i]<<" "<<errore[i]<<endl;
		} //stampo l'andamento di media cumulativa all'aumentare del numero dei blocchi (quindi anche dei lanci totali) e anche le incertezze
	}
