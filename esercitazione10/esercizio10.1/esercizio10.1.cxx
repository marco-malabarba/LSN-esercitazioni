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
int n_individui=1;
double *x = new double[n_cities];
double *y = new double[n_cities];


for (int i=0; i<n_cities; i++) {
	double phi=rnd.Rannyu(0.,2*M_PI);
	x[i]=10*cos(phi);
	y[i]=10*sin(phi);
}

ofstream out("best_cerchio.dat");

popolazione pop(x,y,n_cities,n_individui);
pop.creaindividui();
pop.set_p_mut(1.,0.5,0.3,0.3);
cout<<"Città sulla circonferenza"<<endl;


double temp=100.;

for (; temp>4; temp--) {
	cout<<"Temperatura: "<<temp<<endl;
	for (int i=0; i<100; i++) {
		pop.sim_ann(temp);
		out<<temp<<" "<<pop.best()<<endl;
	}
}


for (; temp>=0.1; temp-=0.1) {
	cout<<"Temperatura: "<<temp<<endl;
	for (int i=0; i<100; i++) {
		pop.sim_ann(temp);
		out<<temp<<" "<<pop.best()<<endl;
	}	

}


pop.stampapercorso(pop.cerca_best(),(string)"cerchio");


out.close();

for (int i=0; i<n_cities; i++) {
	x[i]=rnd.Rannyu(-10.,10.);
	y[i]=rnd.Rannyu(-10.,10.);
}

out.open("best_quad.dat");

popolazione pop1(x,y,n_cities,n_individui);
pop1.creaindividui();
pop1.set_p_mut(1.,0.5,0.3,0.3);
cout<<"Città nel quadrato"<<endl;

temp=100.;

for (; temp>4; temp--) {
	cout<<"Temperatura: "<<temp<<endl;
	for (int i=0; i<200; i++) {
		pop1.sim_ann(temp);
		out<<temp<<" "<<pop1.best()<<endl;
	}
}


for (; temp>=0.1; temp-=0.1) {
	cout<<"Temperatura: "<<temp<<endl;
	for (int i=0; i<200; i++) {
		pop1.sim_ann(temp);
		out<<temp<<" "<<pop1.best()<<endl;
	}	

}


pop1.stampapercorso(pop1.cerca_best(),(string)"quad");
out.close();

return 0;
}
	
	
