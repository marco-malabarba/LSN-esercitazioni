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

double d=1.; //distanza griglia
double l=0.2; //lunghezza ago

//inizio con un metodo più "ingenuo": per stimare pi utilizzo il valore di pi stesso
//parto anche con un piccolo l, poi vedo come cambia la velocità di convergenza per valori più grandi

int M=100000; //numero totale di lanci
int N=100; //numero di blocchi
int L=M/N;

ofstream out;
out.open("pi1.dat");

double *x_ago=new double[M]; //posizione del centro ago dalla riga più vicina (distribuzione uniforme tra 0 e d/2)
double *angolo_ago=new double[M]; //angolo acuto tra ago e linee della griglia (distribuzione uniforme tra 0 e pi/2)

double *media=new double[N]; //definizioni dei successivi array analoghe all'esercizio 1.1 (qui ovviamente ogni componente del vettore conterrà una stima di pi)
double *media_2=new double[N];
double *media_cumulativa=new double[N];
double *media_2_cumulativa=new double[N];
double *errore=new double[N];


for (int i=0; i<M; i++) { //riempio i valori dei due array
	x_ago[i]=rnd.Rannyu(0.,d/2);
	angolo_ago[i]=rnd.Rannyu(0.,M_PI/2); //usato pi!!!!!
	angolo_ago[i]=sin(angolo_ago[i]); //converto subito angolo ---> seno angolo per comodità successiva
	}


calcolapi(x_ago,angolo_ago,media,media_2,l,d,N,L); //per ciascun blocco ho la stima di pi
calcolacumulate(media,media_2,media_cumulativa,media_2_cumulativa,errore,N); //funzione identica a quella esercizio 1.1


for (int i=0; i<N; i++) {
	out<<(i+1)*L<<" "<<media_cumulativa[i]<<" "<<errore[i]<<endl;
	} //stampo l'andamento di pi all'aumentare del numero dei blocchi (quindi anche dei lanci totali) e anche le incertezze
		

out.close();	

//aumento ora l

out.open("pi1_cr.dat"); //cr="convergenza rapida"
l=0.95;
azzeraarray(media,N);
azzeraarray(media_2,N);
azzeraarray(media_cumulativa,N);
azzeraarray(media_2_cumulativa,N);



calcolapi(x_ago,angolo_ago,media,media_2,l,d,N,L);
calcolacumulate(media,media_2,media_cumulativa,media_2_cumulativa,errore,N);
for (int i=0; i<N; i++) {
	out<<(i+1)*L<<" "<<media_cumulativa[i]<<" "<<errore[i]<<endl;
	}
		

out.close();

//provo ora senza utilizzare pi, lascio grande l

M=10000000; //aumento anche M per avere una stima migliore
N=100; //non vario il numero dei blocchi
L=M/N;

double *x=new double[M];
double *y=new double[M];
double *xcm=new double[M];
double *sintheta=new double[M];

out.open("pi2.dat");

azzeraarray(media,N);
azzeraarray(media_2,N);
azzeraarray(media_cumulativa,N);
azzeraarray(media_2_cumulativa,N);


for (int i=0; i<M; i++) {
	x[i]=rnd.Rannyu(-1,1); //genero una coppia di punti x ed y uniformemente distribuiti tra -1 e 1
	y[i]=rnd.Rannyu(-1,1);
	if (pow(x[i],2)+pow(y[i],2)>1) {
		i--; //rigetto il punto (inteso come coppia di coordinate (x;y)) se non appartiene al cerchio inscritto nel quadrato
		}
	else { //se invece il punto è accettabile ho generato una direzione casuale (da intendersi come congiungente origine circonferenza--> punto generato
		xcm[i]=rnd.Rannyu(0,d/2); //posizione del centro ago dalla riga più vicina (distribuzione uniforme tra 0 e d/2)
		sintheta[i]=abs(x[i])/lunghezza(x[i],y[i],0.,0.); //posso calcolare il seno dell'angolo data la direzione casuale generata
		}
	}





calcolapi(xcm,sintheta,media,media_2,l,d,N,L); //basta ora utilizzare la stessa funzione precedente  

calcolacumulate(media,media_2,media_cumulativa,media_2_cumulativa,errore,N);

for (int i=0; i<N; i++) {
	out<<(i+1)*L<<" "<<media_cumulativa[i]<<" "<<errore[i]<<endl;
	}

out.close();
return 0;
}



