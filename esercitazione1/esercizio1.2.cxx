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
   
int M=10000; //numero di dati per ciascun istogramma

int N=1; //sommo solo un dato distribuito uniformemente
ofstream out;
out.open("uniforme1.dat");

creafiledatiunif(rnd,M,N,out); //funzione che stampa sul file M dati ottenuti dalla media degli N generati casualmente (in questo caso con distribuzione continua
out.close();

N=2; //sommo due dati distribuiti uniformemente
out.open("uniforme2.dat");
creafiledatiunif(rnd,M,N,out);
out.close();

N=10; //sommo dieci dati distribuiti uniformemente
out.open("uniforme10.dat");
creafiledatiunif(rnd,M,N,out);
out.close();

N=100; //sommo cento dati distribuiti uniformemente
out.open("uniforme100.dat");
creafiledatiunif(rnd,M,N,out);
out.close();

//tutto analogo, cambia solo la distribuzione dei numeri generati

N=1;
out.open("exp1.dat");
creafiledatiexp(rnd,M,N,out);
out.close();

N=2;
out.open("exp2.dat");
creafiledatiexp(rnd,M,N,out);
out.close();

N=10;
out.open("exp10.dat");
creafiledatiexp(rnd,M,N,out);
out.close();

N=100;
out.open("exp100.dat");
creafiledatiexp(rnd,M,N,out);
out.close();

N=1;
out.open("lor1.dat");
creafiledatilor(rnd,M,N,out);
out.close();

N=2;
out.open("lor2.dat");
creafiledatilor(rnd,M,N,out);
out.close();

N=10;
out.open("lor10.dat");
creafiledatilor(rnd,M,N,out);
out.close();

N=100;
out.open("lor100.dat");
creafiledatilor(rnd,M,N,out);
out.close();

out.open("somme.dat");
creafiledati(rnd,M,N,out);
out.close();

out.open("corr.dat");
creafiledaticorr(rnd,M,N,out);
out.close();

return 0;

}










