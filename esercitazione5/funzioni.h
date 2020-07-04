#include <cmath>
#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include "random.h"
#include <string>

using namespace std;

void setrandom(Random &rnd);
double p_100 (double x, double y, double z);
double p_210 (double x, double y, double z);
double raggio (double x, double y, double z);
void calcolacumulate(double *media, double *media_2, double *media_cumulativa, double *media_2_cumulativa, double *errore, int N);
double calcolaerrore (double *x, double *y, int i);
double costheta(double x, double y, double z);
void azzeraarray(double *x, int n);
void calcolamedia (double *r,double *media, double *media_2, int N, int L);
void stamparisultati (string s, double *media_cumulativa, double *errore, int N, double &acc, int M);
void stepmetropolisunif (int n, Random &rnd, double T[], double *x, double *y, double *z, double *r, double dist, double &acc, int i);
void stepmetropolisgauss (int n, Random &rnd, double T[], double *x, double *y, double *z, double *r, double sigma, double &acc, int i);
void equilibrio(string s, double *r);
