#include <cmath>
#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
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



void calcolaintegrale1 (Random &rnd, double *media, double *media_2, double *x_random, int M, int N, int L) {
	for (int i=0; i<M; i++) {
		x_random[i]=rnd.Rannyu();
		} //riempio il vettore di numeri casuali (distribuiti uniformemente)
	for (int i=0; i<N; i++) { //ciclo sul numero dei blocchi
		int k; 
		double somma=0;
		for (int j=0; j<L; j++) { //ciclo sul numero di lanci in ogni blocco (nota che il loro prodotto restituisce la dimensione dell'array dei numeri casuali)
			k=j+i*L; //serve per selezionare correttamente la componente dell'array dei numeri casuali
			somma=somma+M_PI/2*cos(M_PI/2.*x_random[k]);
			}
		media[i]=somma/L; //def media (per ogni blocco ho la stima integrale)
		media_2[i]=pow(media[i],2); //serve la media di ciascun blocco e il suo valore al quadrato
		}
	}


void calcolaintegrale2 (Random &rnd, double *media, double *media_2, double *x_random, int M, int N, int L) {
	const double pmax=1920./(pow(M_PI,4)-80*pow(M_PI,2)+1920); //costante utile nei calcoli
	for (int i=0; i<M; i++) { //devo generare numeri casuali secondo la distribuzione data dallo sviluppo al 4° ordine in 0
		double r=rnd.Rannyu();
		x_random[i]=rnd.Rannyu();
		if (r>=1-pow(M_PI*x_random[i],2)/8.+pow(M_PI*x_random[i],4)/384.) { //uso la tecnica del rigetto
			i--; //se il punto non soddisfa il criterio scarto e lo rigenero finché non soddisfa
			}
		}
	for (int i=0; i<N; i++) { //ciclo sul numero dei blocchi
		int k; 
		double somma=0;
		for (int j=0; j<L; j++) { //ciclo sul numero di lanci in ogni blocco (nota che il loro prodotto restituisce la dimensione dell'array dei numeri casuali)
			k=j+i*L; //serve per selezionare correttamente la componente dell'array
			somma=somma+M_PI/2*cos(M_PI/2.*x_random[k])/(pmax*(1-pow(M_PI*x_random[k],2)/8.+pow(M_PI*x_random[k],4)/384.));
			}
		media[i]=somma/L; //def media (ogni blocco --> stima integrale)
		media_2[i]=pow(media[i],2); //serve la media di ciascun blocco e il suo valore al quadrato
		}
	}

void calcolaintegrale3 (Random &rnd, double *media, double *media_2, double *x_random, int M, int N, int L) {
	for (int i=0; i<M; i++) { //devo generare numeri casuali secondo la distribuzione data dallo sviluppo al 1° ordine in 0.5
		double r=rnd.Rannyu();
		x_random[i]=(M_PI+4-sqrt(pow(M_PI+4,2)-16*M_PI*r))/(2*M_PI); //uso il metodo dell'inverso della cumulativa
		}
	for (int i=0; i<N; i++) { //ciclo sul numero dei blocchi
		int k; 
		double somma=0;
		for (int j=0; j<L; j++) { //ciclo sul numero di lanci in ogni blocco (nota che il loro prodotto restituisce la dimensione dell'array dei numeri casuali)
			k=j+i*L; //serve per selezionare correttamente la componente dell'array
			somma=somma+M_PI/2*cos(M_PI/2.*x_random[k])/(1-M_PI/2*(x_random[k]-0.5));
			}
		media[i]=somma/L; //def media
		media_2[i]=pow(media[i],2); //serve la media di ciascun blocco e il suo valore al quadrato
		}
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

double calcolaerrore (double x, double y, int N) { //x: gioca il ruolo di media_cumulata; y: media_2_cumulata; N: numero blocchi sommati
	if (N==0) {
		return 0; //se faccio solo una misura (un singolo blocco) non ha senso definire un'incertezza
		}
	return sqrt((y-pow(x,2))/(N-1)); //formula per l'incertezza
	}


void azzeraarray(double *x, int n) {
	for (int i=0; i<n; i++) {
		x[i]=0.;
		}
	}

void stepdiscreto (Random &rnd, double a, double *x, double *y, double *z, int M) {
	for (int i=0; i<M; i++) { //per ognuno degli M random walk che simulo...
		double r=rnd.Rannyu(0.,6.); //la probabilità è la stessa per tutte le sei direzioni in cui può muoversi 
		int v=(int) r;
		if (v==0) { //scelgo questa convenzione arbitraria
			x[i]=x[i]+a;
			}
		if (v==1) {
			x[i]=x[i]-a;
			}
		if (v==2) {
			y[i]=y[i]+a;
			}
		if (v==3) {
			y[i]=y[i]-a;
			}
		if (v==4) {
			z[i]=z[i]+a;
			}
		if (v==5) {
			z[i]=z[i]-a;
			}
		}
	
//poiché i numeri sono uniformemente distribuiti, in media, considerando un grande numero di random walk, i risultati non dovrebbero cambiare
	}

void stepcontinuo (Random &rnd, double a, double *x, double *y, double *z, int M) {
	for (int i=0; i<M; i++) { //per ciascuno dei M random walk che simulo... 
		double phi=rnd.Rannyu(0.,2*M_PI); //formule per campionare UNIFORMEMENTE l'angolo solido 
		double theta=acos(1-2*rnd.Rannyu());
		x[i]=x[i]+a*sin(theta)*cos(phi);
		y[i]=y[i]+a*sin(theta)*sin(phi);
		z[i]=z[i]+a*cos(theta); //aggiorno le posizioni dopo lo step con le formule delle coordinate sferiche
		}
	}
void studiarandomwalk (double *x, double *y, double *z, double *media, double *media_2, int N, int L) {
	for (int i=0; i<N; i++) { //ciclo sul numero dei blocchi
		int k; 
		double somma=0.;
		for (int j=0; j<L; j++) { //ciclo sul numero di lanci in ogni blocco (nota che il loro prodotto restituisce la dimensione dell'array dei numeri casuali)
			k=j+i*L; //serve per selezionare correttamente il random walk
			somma=somma+lunghezza(x[k],y[k],z[k],0.,0.,0.);
			}
		media[i]=somma/L; //ottengo il valore di |r|^2 per ciascun blocco
		media_2[i]=pow(media[i],2); 
		}
		
	}

void stamparisultati(double *media, double *media_2, int N, ofstream &out, int q) {
	double media_cumulata=0.; //il valore atteso è la media dei valori degli N blocchi 
	double media_2_cumulata=0.; //serve per calcolare l'errore

	for (int i=0; i<N; i++) {
		media_cumulata+=media[i];
		media_2_cumulata+=media_2[i];
		}

	media_cumulata=media_cumulata/N; //calcolo la media degli N blocchi
	media_2_cumulata=media_2_cumulata/N; //serve per l'incertezza

	double errore=calcolaerrore(media_cumulata,media_2_cumulata,N); //calcola l'incertezza associata al valore di media_cumulata (media degli N blocchi) 
	
	out<<q+1<<" "<<sqrt(media_cumulata)<<" "<<0.5*errore/(sqrt(media_cumulata))<<endl; //stampo sqrt(|r|^2) e la sua incertezza (formula data dalla propagazione errori
	}

double lunghezza (double x1, double y1, double z1, double x2, double y2, double z2) { //più propriamente lunghezza al quadrato
	return (pow(x1-x2,2)+pow(y1-y2,2)+pow(z1-z2,2));
	}

