#include <iostream>
#include <cmath>
#include <fstream>
#include "random.h"
#include "popolazione.h"
#include "funzioni.h"
#include <string>
#include "mpi.h"

using namespace std;
 

 
int main (int argc, char *argv[]) {
Random rnd;
setrandom(rnd,0);
int n_cities=32;
int n_individui=100;
double *x = new double[n_cities];
double *y = new double[n_cities];

int ngen=2000;
int nmigr=10;

for (int i=0; i<n_cities; i++) {
	x[i]=rnd.Rannyu(-10.,10.);
	y[i]=rnd.Rannyu(-10.,10.);
}

int size, rank;
MPI_Init(&argc,&argv);
MPI_Comm_size(MPI_COMM_WORLD, &size);
MPI_Comm_rank(MPI_COMM_WORLD, &rank);
MPI_Status stat1,stat2;
MPI_Request req;
MPI_Status status; 

if (size!=4) {
	cerr<<"Porre size=4!"<<endl;
	exit(1);
}

popolazione pop(x,y,n_cities,n_individui);

ofstream out;

if (rank==0) {
	pop.setrndm(rank);
	out.open("best_0.dat");
}
if (rank==1) {
	pop.setrndm(rank);
	out.open("best_1.dat");
}
if (rank==2) {
	pop.setrndm(rank);
	out.open("best_2.dat");
}

if (rank==3) {
	pop.setrndm(rank);
	out.open("best_3.dat");
}

pop.creaindividui();

for (int i=0; i<=ngen; i++) {
	if (i%1000==0&&rank==0) cout<<"Generazione: "<<i<<endl;
	pop.calcolalunghezze();
	pop.calcolaprob();
	double len_best=pop.best();
	out<<i<<" "<<len_best<<endl;
	if (i==ngen) {
		int q=pop.cerca_best();
		string tit=to_string(rank);
		pop.stampapercorso(q,tit);
	}
	if (i!=0&&i%nmigr==0) {
		int u=pop.cerca_best();
		int n_sw;
		if (rank==0) n_sw= 2*((int) rnd.Rannyu(1.,5.));
		MPI_Bcast(&n_sw,1,MPI_INTEGER,0, MPI_COMM_WORLD);
		int *sw=new int [n_sw];
		
		if (rank==0) {
			for (int j=0; j<n_sw; j+=2) {
				sw[j]=(int) rnd.Rannyu(0.,4.);
				do {
					sw[j+1]=(int) rnd.Rannyu(0.,4.);
				} while (sw[j+1]==sw[j]);
			}

		}
		MPI_Bcast(sw,n_sw,MPI_INTEGER,0, MPI_COMM_WORLD);
		int *send=new int [n_cities];
		int *receive= new int [n_cities];
		for (int j=0; j<n_cities; j++) send[j]=pop.getpercorso(u,j); //u-esimo perché individuo migliore
		for (int j=0; j<n_sw; j+=2) {
			if (rank==sw[j]) {
				int rr=sw[j+1];
				MPI_Isend(send,n_cities,MPI_INTEGER,rr,0,MPI_COMM_WORLD,&req);
				MPI_Recv(receive,n_cities,MPI_INTEGER,rr,0, MPI_COMM_WORLD,&stat2);
				MPI_Wait(&req, &status);
				for (int jj=0;jj<n_cities; jj++) pop.setpercorso(u,jj,receive[jj]); //u-esimo perché il migliore
			}
			
			if (rank==sw[j+1]) {	
				int rr=sw[j]; 
				MPI_Recv(receive,n_cities,MPI_INTEGER,rr,0, MPI_COMM_WORLD,&stat1);
				MPI_Send(send,n_cities,MPI_INTEGER,rr,0,MPI_COMM_WORLD);
				for (int jj=0;jj<n_cities; jj++) pop.setpercorso(u,jj,receive[jj]);
				
			}
			
			
		
			
		} 	 
	}
	pop.riproduzione();
	pop.mutazioni();
}




MPI_Finalize();
  
return 0;
}


