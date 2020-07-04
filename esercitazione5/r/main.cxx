#include "funzione.cxx"
#include <cmath>
#include <iostream>
#include <random>
#include <fstream>

using namespace std;

double error (double *AV,double *AV2,int n);
double raggio (double x, double y, double z);
double psi(double x, double y, double z);

int main () {



int tiri=10000000;
int blocchi=100;
int L=(int)(tiri/blocchi);

std::default_random_engine generator;
std::uniform_real_distribution<double> distribution(0.0,1.0);


double *ave = new double [blocchi];
double *av2 = new double [blocchi];
double *sum_prog = new double[blocchi];
double *su2_prog = new double[blocchi];
double *err_prog = new double[blocchi];



double X=1;
double Y=1;
double Z=1;
double dist=2.;

for (int i=0; i<1000; i++) {
	double xsucc,ysucc,zsucc,alfa;
	xsucc=X-dist+2*dist*distribution(generator);
	ysucc=Y-dist+2*dist*distribution(generator);
	zsucc=Z-dist+2*dist*distribution(generator);
	alfa=fmin(1.,psi(xsucc,ysucc,zsucc)/psi(X,Y,Z));
	if (distribution(generator)<alfa) {
		X=xsucc;
		Y=ysucc;
		Z=zsucc;
	}
}


for (int i=0; i<blocchi; i++) {
        cout<<"blocco "<<i+1<<endl;
	double summ = 0.;
	for (int j=0; j<L; j++) {
                double xsucc,ysucc,zsucc,alfa;
		xsucc=X-dist+2*dist*distribution(generator);
		ysucc=Y-dist+2*dist*distribution(generator);
		zsucc=Z-dist+2*dist*distribution(generator);
		alfa=fmin(1.,psi(xsucc,ysucc,zsucc)/psi(X,Y,Z));
		if (distribution(generator)<alfa) {
			X=xsucc;
			Y=ysucc;
			Z=zsucc;
		}
		summ += raggio(X,Y,Z);
	}
	ave[i] = summ/L;      
	av2[i] = (ave[i])*(ave[i]);
}

ofstream out("risultati.dat");

for (int i=0; i<blocchi; i++) {
	for (int j=0; j<i+1; j++) { 
		sum_prog[i] += ave[j];
		su2_prog[i] += av2[j]; 
	}
	sum_prog[i]/=(i+1); 
	su2_prog[i]/=(i+1); 
	err_prog[i] = error(sum_prog,su2_prog,i);
	out<<L*(i+1)<<" "<<sum_prog[i]<<" "<<err_prog[i]<<endl;
}

out.close();



return 0;

}

double error (double *AV,double *AV2,int n) {
	if (n==0) return 0.;
	else return sqrt((AV2[n]-AV[n]*AV[n])/n);
} 

double raggio (double x, double y, double z) {
	return sqrt(x*x+y*y+z*z);
}

