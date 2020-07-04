#include <iostream>
#include <cmath>
#include <fstream>
#include <cstdlib>
#include <string>
#include "random.h"
#include "random.cpp"

using namespace std;

double psi (double x, double mu, double sigma);
double potenziale (double x);
double diff_2_psi (double x, double mu, double sigma);
double Hpsi(double x, double mu, double sigma);
void equilibra (double &x, double mu, double sigma, Random &rnd, int step_eq);
void stepmetropolis(double &x, double mu, double sigma, Random &rnd);
void setrandom(Random &rnd);
double calcolavalormediopsi(double &x, double mu, double sigma,Random &rnd, int M);
void cercaparametri (double &mu,double &sigma,double &E_min,double &x,Random &rnd, int M);
void calcolacumulate(double *media, double *media_2, double *media_cumulativa, double *media_2_cumulativa, double *errore, int N);
double calcolaerrore (double *x, double *y, int i);
void stamparisultati (string s, double *media_cumulativa, double *errore, int N);

int main () {

Random rnd;
setrandom(rnd);

double mu=0.;
double sigma=1.;
//double mod_quadro=calcolaintegralemoduloquadro(mu,sigma); //integrale del modulo quadro della PSI_t //non serve per Metropolis!!! 

const double hbar=1.; //scelgo queste unità di misura
const double m=1.; //scelgo queste unità di misura


int M=100000; //numero totale di lanci per ogni coppia (mu;sigma) 
double x=0.;
equilibra (x,mu,sigma,rnd,10000); //suppongo che 10000 step siano sufficienti per l'equilibrazione iniziale
cout<<"Terminata prima fase equilibrazione (effettuata con 10000 step)"<<endl;
cout<<endl;

double E_min=calcolavalormediopsi(x,mu,sigma,rnd,M);


cout<<"Con i parametri:"<<endl;
cout<<"mu="<<mu<<endl;
cout<<"sigma="<<sigma<<endl;
cout<<"Si ottiene un valore medio di H pari a:"<<endl;
cout<<E_min<<endl;	
cout<<endl;


cout<<"Ricerca dei parametri migliori..."<<endl;

cercaparametri(mu,sigma,E_min,x,rnd,M);




cout<<"Con i nuovi parametri:"<<endl;
cout<<"mu="<<mu<<endl;
cout<<"sigma="<<sigma<<endl;
cout<<"Si ottiene un valore medio di H pari a:"<<endl;
cout<<E_min<<endl;	
cout<<endl;

cout<<"Studio del sistema più in dettaglio con questi parametri per valutare anche l'incertezza sull'energia..."<<endl; 
M=100000000;
int N=100; //numero di blocchi (la media a blocchi si utilizzerà solo con i parametri ottimali)
int L=M/N;
double *media=new double[N]; //vettore che contiene <H> (una componente per ogni blocco)
double *media_2=new double[N]; //contiene (<H>)^2
double *media_cumulativa=new double[N]; //la componente i-esima contiene la media tra i risultati dei blocchi fino all'i-esimo
double *media_2_cumulativa=new double[N]; //come sopra ma al quadrato
double *errore=new double[N]; //vettore che contiene le incertezze di ciascuna delle componenti media_cumulativa



for (int i=0; i<N; i++) {
	media[i]=calcolavalormediopsi(x,mu,sigma,rnd,L);
	media_2[i]=media[i]*media[i];
}

cout<<"Medie a blocchi calcolate"<<endl;
calcolacumulate(media,media_2,media_cumulativa,media_2_cumulativa,errore,N);
cout<<"Calcolate le cumulative e gli errori"<<endl;
cout<<endl;

stamparisultati((string) "ground_state.dat",media_cumulativa,errore,N);

double xmin=-3.;
double xmax=3.;
int nbin=(int) ((xmax-xmin)*10.);
double bin_size=(xmax-xmin)/(double) nbin; 
int *counter=new int [nbin];
int valori_hist=10000000;

for (int i=0; i<nbin; i++) counter[i]=0;

for (int i=0; i<valori_hist; i++) {
	stepmetropolis(x,mu,sigma,rnd);
	if (x>=xmin&&x<=xmax) {
		int bin= (int)((x-xmin)/bin_size);
		counter[bin]++;
	}
}

ofstream out;
out.open("histo.dat");
for (int i=0; i<nbin; i++) {
	out<<xmin+i*bin_size+bin_size/2<<" " <<counter[i]/(double) valori_hist <<endl;
}

out.close();
 

return 0;

}

/* NON SERVE PER METROPOLIS (con questa semplice psi_t era calcolabile analiticamente)
double calcolaintegralemoduloquadro (double mu, double sigma) { //integrale del modulo quadro della PSI_t
	return 2*sqrt(M_PI)*sigma*exp(-mu*mu/(sigma*sigma))+2*sqrt(M_PI)*sigma;
}
*/ 

double psi (double x, double mu, double sigma) {
	return exp(-(x-mu)*(x-mu)/(2*sigma*sigma))+exp(-(x+mu)*(x+mu)/(2*sigma*sigma));
}

double potenziale (double x) {
	return pow(x,4)-5./2.*pow(x,2);
}

double diff_2_psi (double x, double mu, double sigma) {
	double primo_fattore=pow(x+mu,2)*exp(-(x+mu)*(x+mu)/(2*sigma*sigma))+pow(x-mu,2)*exp(-(x-mu)*(x-mu)/(2*sigma*sigma));
	primo_fattore=primo_fattore/(sigma*sigma);
	return (primo_fattore-psi(x,mu,sigma))/(sigma*sigma);
}

double Hpsi(double x, double mu, double sigma) {
	return -((1./2.)*diff_2_psi(x,mu,sigma))+(potenziale(x)*psi(x,mu,sigma));
}

void equilibra (double &x, double mu, double sigma, Random &rnd, int step_eq) {
	for (int i=0; i<step_eq; i++) { 
		stepmetropolis(x,mu,sigma,rnd);
	}
	return;
}

void stepmetropolis(double &x, double mu, double sigma, Random &rnd) {
	double x_new=rnd.Rannyu(x-2.,x+2.); //con uno spostamento di 2 ho un'accettazione del 50% circa 
	double alfa=fmin(1.,pow(psi(x_new,mu,sigma)/psi(x,mu,sigma),2)); //psi_t reale--->modulo quadro è equivalente a elevarla al quadrato
	if (rnd.Rannyu()<alfa) { 
		x=x_new;
	}
	return;
}


double calcolavalormediopsi(double &x, double mu, double sigma,Random &rnd, int M) {
	double somma=0.;
	for (int i=0; i<M; i++) {
		stepmetropolis(x,mu,sigma,rnd);
		somma+=Hpsi(x,mu,sigma)/psi(x,mu,sigma);
	}
	return somma/(double)M;
}

void cercaparametri (double &mu,double &sigma,double &E_min,double &x,Random &rnd, int M) {
	for (double mu_new=0.1; mu_new<=1.; mu_new+=0.1) {
		for (double sigma_new=0.1; sigma_new<=1.; sigma_new+=0.1) {
			equilibra (x,mu_new,sigma_new,rnd,1000); //la psi varia di poco, immagino che 500 step per una eventuale riequilibrazione siano sufficienti
			double E_min_new=calcolavalormediopsi(x,mu_new,sigma_new,rnd,M);
			if (E_min_new<E_min) {
				E_min=E_min_new;
				sigma=sigma_new;
				mu=mu_new;
			}
		}
	}

	cout<<"Termine prima fase"<<endl;


	for (double mu_new=mu-0.09; mu_new<=mu+0.09; mu_new+=0.01) {
		for (double sigma_new=sigma-0.09; sigma_new<=sigma+0.09; sigma_new+=0.01) {
		equilibra (x,mu_new,sigma_new,rnd,1000); 
			double E_min_new=calcolavalormediopsi(x,mu_new,sigma_new,rnd,M);
			if (E_min_new<E_min) {
				E_min=E_min_new;
				sigma=sigma_new;
				mu=mu_new;
			}
		}
	}		

	cout<<"Ottenuta precisione sulla prima cifra decimale"<<endl;

	for (double mu_new=mu-0.009; mu_new<=mu+0.009; mu_new+=0.001) {
		for (double sigma_new=sigma-0.009; sigma_new<=sigma+0.009; sigma_new+=0.001) {
			equilibra (x,mu_new,sigma_new,rnd,1000); 
			double E_min_new=calcolavalormediopsi(x,mu_new,sigma_new,rnd,M);
			if (E_min_new<E_min) {
				E_min=E_min_new;
				sigma=sigma_new;
				mu=mu_new;
			}
		}
	}		

	cout<<"Ottenuta precisione sulla seconda cifra decimale"<<endl;

	for (double mu_new=mu-0.0009; mu_new<=mu+0.0009; mu_new+=0.0001) {
		for (double sigma_new=sigma-0.0009; sigma_new<=sigma+0.0009; sigma_new+=0.0001) {
			equilibra (x,mu_new,sigma_new,rnd,1000); 
			double E_min_new=calcolavalormediopsi(x,mu_new,sigma_new,rnd,M);
			if (E_min_new<E_min) {
				E_min=E_min_new;
				sigma=sigma_new;
				mu=mu_new;
			}
		}
	}		

	cout<<"Ottenuta precisione sulla terza cifra decimale"<<endl;


	for (double mu_new=mu-0.00009; mu_new<=mu+0.00009; mu_new+=0.00001) {
		for (double sigma_new=sigma-0.00009; sigma_new<=sigma+0.00009; sigma_new+=0.00001) {
			equilibra (x,mu_new,sigma_new,rnd,1000); 
			double E_min_new=calcolavalormediopsi(x,mu_new,sigma_new,rnd,M);
			if (E_min_new<E_min) {
				E_min=E_min_new;
				sigma=sigma_new;
				mu=mu_new;
			}
		}
	}		


	cout<<"Ottenuta precisione sulla quarta cifra decimale"<<endl;

	for (double mu_new=mu-0.000009; mu_new<=mu+0.000009; mu_new+=0.000001) {
		for (double sigma_new=sigma-0.000009; sigma_new<=sigma+0.000009; sigma_new+=0.000001) {
			equilibra (x,mu_new,sigma_new,rnd,1000); 
			double E_min_new=calcolavalormediopsi(x,mu_new,sigma_new,rnd,M);
			if (E_min_new<E_min) {
				E_min=E_min_new;
				sigma=sigma_new;
				mu=mu_new;
			}
		}
	}		
	
	cout<<"Ottenuta precisione sulla quinta cifra decimale"<<endl;
	cout<<endl;
	return;
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

void stamparisultati (string s, double *media_cumulativa, double *errore, int N) {
	ofstream out(s);
	for (int i=0; i<N; i++) {
		out<<(i+1)<<" "<<media_cumulativa[i]<<" "<<errore[i]<<endl;
	} //stampo l'andamento di media cumulativa all'aumentare del numero dei blocchi e anche le incertezze
	out.close();
}


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

