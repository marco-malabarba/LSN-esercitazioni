#include <iostream>
#include <cmath>
#include <fstream>
#include "random.h"
#include "popolazione.h"
#include "funzioni.h"
#include <string>

using namespace std;


int main () {

Random rnd;
setrandom(rnd);
int n_cities=32;
int n_individui=100;
double *x = new double[n_cities];
double *y = new double[n_cities];
int ngen=20000;

for (int i=0; i<n_cities; i++) {
	double phi=rnd.Rannyu(0.,2*M_PI);
	x[i]=10*cos(phi);
	y[i]=10*sin(phi);
}

ofstream out1("ave_cerchio.dat");

popolazione pop(x,y,n_cities,n_individui);
pop.creaindividui();

cout<<"Città sulla circonferenza"<<endl;

for (int i=0; i<=ngen; i++) {
	if (i%1000==0) cout<<"Generazione: "<<i<<endl;
	pop.calcolalunghezze();
	pop.calcolaprob();
	double ave=pop.ave_best();
	out1<<i<<" "<<ave<<endl;
	if (i==ngen) {
		int q=pop.cerca_best();
		pop.stampapercorso(q,(string)"cerchio");
	}
	pop.riproduzione();
	pop.mutazioni();
}

cout<<endl;

out1.close();

for (int i=0; i<n_cities; i++) {
	x[i]=rnd.Rannyu(-10.,10.);
	y[i]=rnd.Rannyu(-10.,10.);
}

out1.open("ave_quad.dat");

popolazione pop1(x,y,n_cities,n_individui);
pop1.creaindividui();

cout<<"Città nel quadrato"<<endl;

for (int i=0; i<=ngen; i++) {
	if (i%1000==0) cout<<"Generazione: "<<i<<endl;
	pop1.calcolalunghezze();
	pop1.calcolaprob();
	double ave=pop1.ave_best();
	out1<<i<<" "<<ave<<endl;
	if (i==ngen) {
		int q=pop1.cerca_best();
		pop1.stampapercorso(q,(string)"quad");
	}
	pop1.riproduzione();
	pop1.mutazioni();
}

out1.close();

return 0;
}
 


