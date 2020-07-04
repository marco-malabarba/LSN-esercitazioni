#include <cmath>
#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <string>
#include "random.h"

using namespace std;

void setrandom(Random &rnd);
double calcolaerrore (double *x, double *y, int i);
void calcolacumulate(double *media, double *media_2, double *media_cumulativa, double *media_2_cumulativa, double *errore, int N);
void azzeraarray(double *x, int n);
void calcolaprezzo (Random &rnd, const double S0,const double t0,const double T,const double K,const double r,const double sigma,double *media, double *media_2,int N, int L, int intervalli,string s);
void stamparisultati(ofstream &out,double *media_cumulativa,double *errore,int L, int N);
