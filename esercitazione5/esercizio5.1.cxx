#include <cmath>
#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include "random.h"
#include "funzioni.h"
#include <string>

using namespace std;

int main () {

Random rnd;
setrandom(rnd);

int M=100000000; //numero totale di lanci
int N=100; //numero di blocchi
int L=M/N;


double *x=new double [M];
double *y=new double [M];
double *z=new double [M];
double *r=new double [M];
double T[3];
double *media=new double[N]; //vettore che contiene <r> (una componente per ogni blocco)
double *media_2=new double[N]; //contiene (<r>)^2
double *media_cumulativa=new double[N]; //la componente i-esima contiene la media tra i risultati dei blocchi fino all'i-esimo
double *media_2_cumulativa=new double[N]; //come sopra ma al quadrato
double *errore=new double[N]; //vettore che contiene le incertezze di ciascuna delle componenti media_cumulativa  

double dist=1.2;
double acc=0;
double sigma=0.75;

x[0]=sqrt(3./4.);
y[0]=sqrt(3./4.);
z[0]=sqrt(3./4.);
r[0]=raggio(x[0],y[0],z[0]);



for (int i=0; i<M-1; i++) {
       stepmetropolisunif(1,rnd,T,x,y,z,r,dist,acc,i);
} 

equilibrio((string) "andamento.dat", r);

calcolamedia(r,media,media_2,N,L);

calcolacumulate(media,media_2,media_cumulativa,media_2_cumulativa,errore,N);


stamparisultati((string)"media_100.dat",media_cumulativa,errore,N,acc,M);




azzeraarray(media,N);
azzeraarray(media_2,N);
azzeraarray(media_cumulativa,N);
azzeraarray(media_2_cumulativa,N);

x[0]=20.;
y[0]=20.;
z[0]=20.;
r[0]=raggio(x[0],y[0],z[0]);



for (int i=0; i<M-1; i++) {
        stepmetropolisunif(1,rnd,T,x,y,z,r,dist,acc,i);
} 

equilibrio((string)"andamento_far.dat", r);
calcolamedia(r,media,media_2,N,L);
calcolacumulate(media,media_2,media_cumulativa,media_2_cumulativa,errore,N);

stamparisultati((string)"media_100_far.dat",media_cumulativa,errore,N,acc,M);

azzeraarray(media,N);
azzeraarray(media_2,N);
azzeraarray(media_cumulativa,N);
azzeraarray(media_2_cumulativa,N);


x[0]=0.;
y[0]=0.;
z[0]=5.;
r[0]=raggio(x[0],y[0],z[0]);
dist=3.;


for (int i=0; i<M-1; i++) {
	stepmetropolisunif(2,rnd,T,x,y,z,r,dist,acc,i); //2: #quantico principale
} 




calcolamedia(r,media,media_2,N,L);

calcolacumulate(media,media_2,media_cumulativa,media_2_cumulativa,errore,N);

stamparisultati((string)"media_210.dat",media_cumulativa,errore,N,acc,M);




azzeraarray(media,N);
azzeraarray(media_2,N);
azzeraarray(media_cumulativa,N);
azzeraarray(media_2_cumulativa,N);

x[0]=50.;
y[0]=50.;
z[0]=50.;
r[0]=raggio(x[0],y[0],z[0]);



for (int i=0; i<M-1; i++) {
        stepmetropolisunif(2,rnd,T,x,y,z,r,dist,acc,i);
} 



calcolamedia(r,media,media_2,N,L);

calcolacumulate(media,media_2,media_cumulativa,media_2_cumulativa,errore,N);

stamparisultati((string)"media_210_far.dat",media_cumulativa,errore,N,acc,M);

azzeraarray(media,N);
azzeraarray(media_2,N);
azzeraarray(media_cumulativa,N);
azzeraarray(media_2_cumulativa,N);



x[0]=sqrt(3./4.);
y[0]=sqrt(3./4.);
z[0]=sqrt(3./4.);
r[0]=raggio(x[0],y[0],z[0]);

for (int i=0; i<M-1; i++) {
	stepmetropolisgauss(1,rnd,T,x,y,z,r,sigma,acc,i);
} 



calcolamedia(r,media,media_2,N,L);

calcolacumulate(media,media_2,media_cumulativa,media_2_cumulativa,errore,N);

stamparisultati((string)"media_100_gauss.dat",media_cumulativa,errore,N,acc,M);

azzeraarray(media,N);
azzeraarray(media_2,N);
azzeraarray(media_cumulativa,N);
azzeraarray(media_2_cumulativa,N);

acc=0;
sigma=1.9;
x[0]=0.;
y[0]=0.;
z[0]=5.;
r[0]=raggio(x[0],y[0],z[0]);

for (int i=0; i<M-1; i++) {
        stepmetropolisgauss(2,rnd,T,x,y,z,r,sigma,acc,i);
} 



calcolamedia(r,media,media_2,N,L);

calcolacumulate(media,media_2,media_cumulativa,media_2_cumulativa,errore,N);

stamparisultati((string)"media_210_gauss.dat",media_cumulativa,errore,N,acc,M);


return 0;
}





 

	
