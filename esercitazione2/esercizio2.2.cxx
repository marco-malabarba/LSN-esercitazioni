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
int M=10000; //totale dei random walk
int N=100; //numero di blocchi
int L=M/N;
int step=100; //per ogni random-walk faccio studio fino a 100 step

double a=1.; //lunghezza di ciascun passo
ofstream out;
out.open("RWdiscreto.dat");

double *x=new double [M];
double *y=new double [M];
double *z=new double [M]; 
//per ciascuno degli M random walk salvo la posizione dopo ogni step per poterne calcolare media e incertezza

azzeraarray(x,M);
azzeraarray(y,M);
azzeraarray(z,M);

double *media=new double[N]; //valori di |r|^2 per ciascun blocco
double *media_2=new double[N]; //i suoi componenti sono quelli dell'array sopra al quadrato

for (int q=0; q<step; q++) {
	stepdiscreto(rnd,a,x,y,z,M); //mi sposto di a in una delle direzioni sei direzioni casuali permesse
	studiarandomwalk(x,y,z,media,media_2,N,L); //calcola, per ogni blocco, i valori di |r|^2
	stamparisultati(media,media_2,N,out,q); //fa la media degli |r|^2 di ogni blocco, ne calcola l'incertezza e stampa i risultati
	}

out.close();
out.open("RWcontinuo.dat");


azzeraarray(x,M);
azzeraarray(y,M);
azzeraarray(z,M);

for (int q=0; q<step; q++) {
	stepcontinuo(rnd,a,x,y,z,M); //mi sposto di a(=1) con una distribuzione di angolo solido
	studiarandomwalk(x,y,z,media,media_2,N,L);
	stamparisultati(media,media_2,N,out,q);
	}

out.close();

azzeraarray(x,M);
azzeraarray(y,M);
azzeraarray(z,M);
out.open("RWcontinuo2.dat");


M=100000;
L=M/N;

x=new double [M];
y=new double [M];
z=new double [M]; 

for (int q=0; q<step; q++) {
	stepcontinuo(rnd,a,x,y,z,M); //mi sposto di a con una distribuzione di angolo solido
	studiarandomwalk(x,y,z,media,media_2,N,L);
	stamparisultati(media,media_2,N,out,q);
	}

out.close();

return 0;
}

