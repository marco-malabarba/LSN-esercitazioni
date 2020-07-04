#include <cmath>
#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include "random.h"
#include "funzioni.h"

using namespace std;

int main () {

Random rnd;
setrandom(rnd);

int M=100000; //numero totale di passi Monte Carlo (un passo=una stima del costo)
int N=100; //numero di blocchi
int L=M/N; 

const double S0=100.; //prezzo iniziale
const double t0=0.; //tempo iniziale
const double T=1.; //tempo termine contratto
const double K=100.; //prezzo pattuito
const double r=0.1; //tasso di interesse
const double sigma=0.25; //volatilità

double *media=new double[N]; //vettore che contiene il costo dell'opzione per ogni blocco
double *media_2=new double[N]; //contiene il costo per ogni blocco al quadrato
double *media_cumulativa=new double[N]; //la componente i-esima contiene la media tra i costi nei blocchi fino all'i-esimo
double *media_2_cumulativa=new double[N]; //come sopra ma al quadrato
double *errore=new double[N]; //vettore che contiene le incertezze di ciascuna delle componenti media_cumulativa  

ofstream out; 
out.open("calldiretta.dat"); 

calcolaprezzo(rnd,S0,t0,T,K,r,sigma,media,media_2,N,L,1,"call"); 
//calcolaprezzo permette di calcolare il costo dell'opzione per ogni blocco; il penultimo input è il numero di intervalli temporali in cui si divide [t0;T] (nella prima parte sarà 1, successivamente 100)  
calcolacumulate(media,media_2,media_cumulativa,media_2_cumulativa,errore,N); //calcolo le medie cumulate e associo gli errori
stamparisultati(out,media_cumulativa,errore,L,N); //stampo sul file il numero di passi, il costo e l'incertezza sul costo 

out.close();

azzeraarray(media_cumulativa,N);
azzeraarray(media_2_cumulativa,N);

out.open("putdiretta.dat");

calcolaprezzo(rnd,S0,t0,T,K,r,sigma,media,media_2,N,L,1,"put");
calcolacumulate(media,media_2,media_cumulativa,media_2_cumulativa,errore,N); //calcolo le medie cumulate e associo gli errori
stamparisultati(out,media_cumulativa,errore,L,N);

out.close();

azzeraarray(media_cumulativa,N);
azzeraarray(media_2_cumulativa,N);

out.open("calldiscreta.dat");

int intervalli=100;

calcolaprezzo(rnd,S0,t0,T,K,r,sigma,media,media_2,N,L,intervalli,"call");
calcolacumulate(media,media_2,media_cumulativa,media_2_cumulativa,errore,N); //calcolo le medie cumulate e associo gli errori
stamparisultati(out,media_cumulativa,errore,L,N);

out.close();

azzeraarray(media_cumulativa,N);
azzeraarray(media_2_cumulativa,N);

out.open("putdiscreta.dat");

calcolaprezzo(rnd,S0,t0,T,K,r,sigma,media,media_2,N,L,intervalli,"put");
calcolacumulate(media,media_2,media_cumulativa,media_2_cumulativa,errore,N); //calcolo le medie cumulate e associo gli errori
stamparisultati(out,media_cumulativa,errore,L,N);

out.close();
return 0;
}

