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

int M=10000; //numero totale di lanci
int N=100; //numero di blocchi
int L=M/N;

ofstream out;
out.open("media.dat");
double *x_random=new double[M]; //vettore di numeri casuali

for (int i=0; i<M; i++) {
	x_random[i]=rnd.Rannyu();
	} //riempio il vettore di numeri casuali


double *media=new double[N]; //vettore che contiene <r> (una componente per ogni blocco)
double *media_2=new double[N]; //contiene (<r>)^2
double *media_cumulativa=new double[N]; //la componente i-esima contiene la media tra i risultati dei blocchi fino all'i-esimo
double *media_2_cumulativa=new double[N]; //come sopra ma al quadrato
double *errore=new double[N]; //vettore che contiene le incertezze di ciascuna delle componenti media_cumulativa  


calcolamedie(media,media_2,x_random,N,L); //riempie gli array media e media_2

calcolacumulate(media,media_2,media_cumulativa,media_2_cumulativa,errore,N); //riempio le medie cumulate e 

for (int i=0; i<N; i++) {
	out<<(i+1)*L<<" "<<media_cumulativa[i]-0.5<<" "<<errore[i]<<endl; //sottraggo 0.5 per avere un più facile confronto col valore zero
	} //stampo l'andamento di media cumulativa all'aumentare del numero dei blocchi (quindi anche dei lanci totali) e anche le incertezze

out.close();

out.open("varianza.dat");


azzeraarray(media_cumulativa,N);
azzeraarray(media_2_cumulativa,N);


calcolavarianze(media,media_2,x_random,N,L); //come prima ma calcola il valore medio di ogni blocco di (r-<r>)^2

calcolacumulate(media,media_2,media_cumulativa,media_2_cumulativa,errore,N); //come prima


for (int i=0; i<N; i++) {
	out<<(i+1)*L<<" "<<media_cumulativa[i]-(double)1/12<<" "<<errore[i]<<endl;
	}
		

out.close();	


out.open("chiquadro.dat");

int n=10000; //numeri casuali ad ogni test del chi quadro
int m_int=100; //sottointervalli in cui divido [0;1)
int n_test=100; //calcolo 100 volte il chi quadro

double *x_chiquadro=new double [n]; //sono gli n numeri casuali di ogni test
double *x_i=new double [m_int]; //array della dimensione del numero dei sottointervalli, conterà quante occorrenze ci sono ad ogni sottointervallo
double *chiquadro=new double [n_test]; //registra i valori del chi-quadro

azzeraarray(x_i,m_int);
calcolachiquadro(rnd,x_chiquadro, x_i, chiquadro, n ,m_int,n_test); //calcola i 100 valori del chi quadro 


for (int i=0; i<n_test; i++) {
	out<<i+1<<" "<<chiquadro[i]<<endl;
	} //stampo i 100 valori del chi quadro

out.close();


return 0;
}






	 
