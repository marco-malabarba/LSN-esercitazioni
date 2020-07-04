#include "funzioni.h"


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

double distanza (double x1, double y1, double x2, double y2) {
	return (pow(x2-x1,2)+pow(y2-y1,2));
}

void swap (double &x, double &y) {
	double z=x;
	x=y;
	y=z;
}

int pbc (int x, int tot) {
	if (x<tot&&x>=1) return x;
	if(x>=tot) return x-tot+1; 
	if (x<1) return x+tot-1;
	cerr<<"Errore!"<<endl;
	exit(1);
}

bool cerca(int x, int *y, int c) {
	bool res=false;
	for (int i=0; i<c; i++) {
		if (y[i]==x) res=true;
	}
	return res;
}

void sort (double *x, int n) {
	for (int i=0; i<n; i++) {
		for (int j=i+1; j<n; j++) {
			if (x[j]<x[i]) swap(x[j],x[i]);
		}
	}
}
