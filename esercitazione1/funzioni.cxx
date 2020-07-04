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



void calcolamedie (double *media, double *media_2, double *x_random, int N, int L) {
	for (int i=0; i<N; i++) { //ciclo sul numero dei blocchi
		int k; 
		double somma=0;
		for (int j=0; j<L; j++) { //ciclo sul numero di lanci in ogni blocco (nota che il loro prodotto restituisce la dimensione dell'array dei numeri casuali)
			k=j+i*L; //serve per selezionare correttamente la componente dell'array dei numeri casuali
			somma=somma+x_random[k];
			}
		media[i]=somma/L; //def media
		media_2[i]=pow(media[i],2); //serve la media di ciascun blocco e il suo valore al quadrato
		}
	}


void calcolavarianze (double *media, double *media_2, double *x_random, int N, int L) {
	for (int i=0; i<N; i++) {
		int k;
		double somma=0;
		for (int j=0; j<L; j++) {
			k=j+i*L;
			somma=somma+pow(x_random[k]-0.5,2); //unica differenza con la funzione precedente 
			}
		media[i]=somma/L; 
		media_2[i]=pow(media[i],2);
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



void azzeraarray(double *x, int n) {
	for (int i=0; i<n; i++) {
		x[i]=0.;
		}
	}




void calcolachiquadro (Random rnd, double *x_chiquadro,double *x_i,double *chiquadro,int n , int m_int, int m_test) {
	for (int i=0; i<m_test; i++) { //ciclo sui 100 valori di chi-quadro desiderati
		for (int j=0; j<n; j++) { //ciclo sugli n numeri casuali utilizzati per ciascun calcolo di chi-quadro
			x_chiquadro[j]=rnd.Rannyu(); //genero il numero casuale
			int k=(int)(x_chiquadro[j]*m_int); //capisco a quale sottointervallo di [0;1) appartiene il numero generato k
			x_i[k]=x_i[k]+1; //aumento l'occorrenza di quel sottointervallo (k={0,1,....,m_int-1} è giusto che sia)
			}
		double somma=0;
		for (int j=0; j<m_int; j++) { //implemto il numeratore della formula del chi-quadro
			somma+=pow(x_i[j]-(double)n/((double)m_int),2); 
			}
		somma=somma/((double)n/((double)m_int));
		chiquadro[i]=somma; //completo la formula chi-quadro
		azzeraarray(x_i,m_int);
		}
	}

void creafiledatiunif(Random rnd, int M, int N, ofstream &out) {
	for (int i=0; i<M; i++) { //ciclo sul numero di dati M che voglio generare
		double S_n=0;
		for (int j=0; j<N; j++) { //ciascun dato è la media di N numeri casuali interi...
			S_n=S_n+(int) rnd.Rannyu(1.,7.);//... distribuiti uniformemente tra 1 e 6. Per implementare ciò basta utilizzare il generatore continuo di numeri casuali (distr. unif.) e ricordarsi che (int) arrotonda per difetto 
			}
		S_n=S_n/N;
		out<<S_n<<endl;
		}
	}

void creafiledatiexp(Random rnd, int M, int N, ofstream &out) { //cambia solo la distribuzione dei dati che voglio generare
	for (int i=0; i<M; i++) {
		double S_n=0;
		for (int j=0; j<N; j++) {
			S_n=S_n+rnd.exp(1.);
			}
		S_n=S_n/N;
		out<<S_n<<endl;
		}
	}

void creafiledatilor(Random rnd, int M, int N, ofstream &out) {
	for (int i=0; i<M; i++) {
		double S_n=0;
		for (int j=0; j<N; j++) {
			S_n=S_n+ rnd.lorentz(0.,1.);
			}
		S_n=S_n/N;
		out<<S_n<<endl;
		}
	}

void creafiledati(Random rnd, int M, int N, ofstream &out) {
	for (int i=0; i<M; i++) {
		double S_n=0;
		for (int j=0; j<N; j++) {
			S_n=S_n+ rnd.exp(1.)+rnd.Rannyu(0.,10.);
			}
		S_n=S_n/(2*N); //2 perché sto facendo la somma di 200 numeri in realtà (100 exp e 100 uniformi)
		out<<S_n<<endl;
		}
	}

void creafiledaticorr(Random rnd, int M, int N, ofstream &out) {
	for (int i=0; i<M; i++) {
		double S_n=0;
		for (int j=0; j<N; j++) {
			double x=rnd.exp(1.);
			S_n=S_n+ x+2*x;
			}
		S_n=S_n/(2*N); //2 per lo stesso motivo precedente
		out<<S_n<<endl;
		}
	}

double lunghezza (double x1, double y1, double x2, double y2) {
	return sqrt(pow(x1-x2,2)+pow(y1-y2,2));
	}

void calcolapi (double *x_ago,double *angolo_ago, double *media, double *media_2, double l,double d,double N,double L) {
	for (int i=0; i<N; i++) { //ciclo sul numero di blocchi
		int k;
		int nhit=0;
		for (int j=0; j<L; j++) { //ciclo sul numero di tiri per blocco
			k=j+i*L; //per selezionare la componente corretta degli array riempiti di posizioni e angoli casuali
			if (x_ago[k]<l/2*(angolo_ago[k])) { //condizione geometrica per avere intersezione riga/ago
				nhit++; 
				}
			}
		media[i]=2*l*L/(nhit*d); //formula per pi di ciascun blocco
		media_2[i]=pow(media[i],2);
		}
	}
