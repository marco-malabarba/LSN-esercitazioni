#include <cmath>
#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include "funzioni.h"
#include <string>

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

double p_100 (double x, double y, double z) {
	return (pow(1/sqrt(M_PI)*exp(-raggio(x,y,z)),2));
}

double p_210 (double x, double y, double z) {
	return pow(1/8.*sqrt(2./M_PI)*raggio(x,y,z)*exp(-raggio(x,y,z)/2.)*costheta(x,y,z),2);
}

double raggio (double x, double y, double z) {
	return sqrt(x*x+y*y+z*z);
}

double costheta(double x, double y, double z) {
	return z/raggio(x,y,z);
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



void calcolamedia (double *r,double *media, double *media_2, int N, int L) {
	for (int i=0; i<N; i++) { //ciclo sul numero dei blocchi
		int k; 
		double somma=0;
		for (int j=0; j<L; j++) { 
			k=j+i*L; //serve per selezionare correttamente la componente dell'array dei numeri casuali
			somma=somma+r[k];
		}
		media[i]=somma/L; //def media
		media_2[i]=pow(media[i],2); //serve la media di ciascun blocco e il suo valore al quadrato
	}
}

void stamparisultati (string s, double *media_cumulativa, double *errore, int N, double &acc, int M) {
	int l=(M/N);
	ofstream out(s);
	for (int i=0; i<N; i++) {
		out<<(i+1)*l<<" "<<media_cumulativa[i]<<" "<<errore[i]<<endl;
	} //stampo l'andamento di media cumulativa all'aumentare del numero dei blocchi (quindi anche dei lanci totali) e anche le incertezze

	acc=(double)acc/(double)M;
	cout<<"Probabilità accettazione: "<<acc<<endl;
	acc=0;
	out.close();
}

void stepmetropolisunif (int n, Random &rnd, double T[], double *x, double *y, double *z, double *r, double dist, double &acc, int i) {
	if (n!=1&&n!=2) {
		cout<<"Errore!!!"<<endl;
		exit(1);
	}
	if (n==2) {
		double alfa;
		T[0]=rnd.Rannyu(x[i]-dist,x[i]+dist);
		T[1]=rnd.Rannyu(y[i]-dist,y[i]+dist);   
		T[2]=rnd.Rannyu(z[i]-dist,z[i]+dist);
        	alfa=fmin(1,(p_210(T[0],T[1],T[2])/p_210(x[i],y[i],z[i])));
		if (rnd.Rannyu()<=alfa) {
			x[i+1]=T[0];
			y[i+1]=T[1];
			z[i+1]=T[2];
			acc++;
		}
		else {
			x[i+1]=x[i];
			y[i+1]=y[i];
			z[i+1]=z[i];
		}
		r[i+1]=raggio(x[i+1],y[i+1],z[i+1]);
	}
	if (n==1) {
		double alfa;
		T[0]=rnd.Rannyu(x[i]-dist,x[i]+dist);
		T[1]=rnd.Rannyu(y[i]-dist,y[i]+dist);   
		T[2]=rnd.Rannyu(z[i]-dist,z[i]+dist);
        	alfa=fmin(1,(p_100(T[0],T[1],T[2])/p_100(x[i],y[i],z[i])));
		if (rnd.Rannyu()<=alfa) {
			x[i+1]=T[0];
			y[i+1]=T[1];
			z[i+1]=T[2];
			acc++;
		}
		else {
			x[i+1]=x[i];
			y[i+1]=y[i];
			z[i+1]=z[i];
		}
		r[i+1]=raggio(x[i+1],y[i+1],z[i+1]);
	}
}	

void stepmetropolisgauss (int n, Random &rnd, double T[], double *x, double *y, double *z, double *r, double sigma, double &acc, int i) {
	if (n!=1&&n!=2) {
		cout<<"Errore!!!"<<endl;
		exit(1);
	}
	if (n==2) {
		double alfa;
		T[0]=rnd.Gauss(x[i],sigma);
		T[1]=rnd.Gauss(y[i],sigma);   
		T[2]=rnd.Gauss(z[i],sigma);
        	alfa=fmin(1,(p_210(T[0],T[1],T[2])/p_210(x[i],y[i],z[i])));
		if (rnd.Rannyu()<=alfa) {
			x[i+1]=T[0];
			y[i+1]=T[1];
			z[i+1]=T[2];
			acc++;
		}
		else {
			x[i+1]=x[i];
			y[i+1]=y[i];
			z[i+1]=z[i];
		}
		r[i+1]=raggio(x[i+1],y[i+1],z[i+1]);
	}
	if (n==1) {
		double alfa;
		T[0]=rnd.Gauss(x[i],sigma);
		T[1]=rnd.Gauss(y[i],sigma);   
		T[2]=rnd.Gauss(z[i],sigma);
        	alfa=fmin(1,(p_100(T[0],T[1],T[2])/p_100(x[i],y[i],z[i])));
		if (rnd.Rannyu()<=alfa) {
			x[i+1]=T[0];
			y[i+1]=T[1];
			z[i+1]=T[2];
			acc++;
		}
		else {
			x[i+1]=x[i];
			y[i+1]=y[i];
			z[i+1]=z[i];
		}
		r[i+1]=raggio(x[i+1],y[i+1],z[i+1]);
	}
}	

void equilibrio(string s, double *r) {
  ofstream out(s);
  for (int i=0; i<1000; i++) out<<r[i]<<endl;
  out.close();
}

