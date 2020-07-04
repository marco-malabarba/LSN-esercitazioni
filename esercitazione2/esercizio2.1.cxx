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

int M=1000000; //numero totale di lanci
int N=100; //numero di blocchi
int L=M/N;

ofstream out;
out.open("integrale1.dat");
double *x_random=new double[M]; //vettore che conterrà, ogni volta, numeri distribuiti secondo la distribuzione con cui si vuole calcolare l'integrale




double *media=new double[N]; //vettore che contiene stima integrale per ogni blocco
double *media_2=new double[N]; //contiene la stima per ogni blocco al quadrato
double *media_cumulativa=new double[N]; //la componente i-esima contiene la media tra i risultati dei blocchi fino all'i-esimo
double *media_2_cumulativa=new double[N]; //come sopra ma al quadrato
double *errore=new double[N]; //vettore che contiene le incertezze di ciascuna delle componenti media_cumulativa  


calcolaintegrale1(rnd,media,media_2,x_random,M,N,L); //calcola integrale con una distribuzione di probabilità p(x) uniforme

calcolacumulate(media,media_2,media_cumulativa,media_2_cumulativa,errore,N); //calcolo le medie cumulate e associo gli errori

for (int i=0; i<N; i++) {
	out<<(i+1)*L<<" "<<media_cumulativa[i]<<" "<<errore[i]<<endl;
	} //stampo l'andamento di media cumulativa all'aumentare del numero dei blocchi (quindi anche dei lanci totali) e anche le incertezze

out.close();

out.open("integrale2.dat");
azzeraarray(media_cumulativa,N);
azzeraarray(media_2_cumulativa,N);



calcolaintegrale2(rnd,media,media_2,x_random,M,N,L); //come prima ma ora si utilizza una distribuzione derivante dallo sviluppo di Taylor al 4° ordine della funzione
calcolacumulate(media,media_2,media_cumulativa,media_2_cumulativa,errore,N); //come prima


for (int i=0; i<N; i++) {
	out<<(i+1)*L<<" "<<media_cumulativa[i]<<" "<<errore[i]<<endl;
	}
		

out.close();

out.open("integrale3.dat");
azzeraarray(media_cumulativa,N);
azzeraarray(media_2_cumulativa,N);

calcolaintegrale3(rnd,media,media_2,x_random,M,N,L); //uso distribuzione data dallo sviluppo di Taylor al 1°ordine in 0.5
calcolacumulate(media,media_2,media_cumulativa,media_2_cumulativa,errore,N); //come prima

for (int i=0; i<N; i++) {
	out<<(i+1)*L<<" "<<media_cumulativa[i]<<" "<<errore[i]<<endl;
	}

out.close();
	
return 0;
}

