#include <cmath>
#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include "random.h"

using namespace std;

void setrandom(Random &rnd);
double calcolaerrore (double *x, double *y, int i);
double calcolaerrore (double x, double y, int N);
void calcolacumulate(double *media, double *media_2, double *media_cumulativa, double *media_2_cumulativa, double *errore, int N);
void azzeraarray(double *x, int n);
void calcolaintegrale1 (Random &rnd, double *media, double *media_2, double *x_random, int M, int N, int L);
void calcolaintegrale2 (Random &rnd, double *media, double *media_2, double *x_random, int M, int N, int L);
void calcolaintegrale3 (Random &rnd, double *media, double *media_2, double *x_random, int M, int N, int L);
void stepdiscreto (Random &rnd, double a, double *x, double *y, double *z, int M);
void stepcontinuo (Random &rnd, double a, double *x, double *y, double *z, int M);
void studiarandomwalk (double *x, double *y, double *z, double *media, double *media_2, int N, int L);
void stamparisultati(double *media, double *media_2, int N, ofstream &out, int q);
double lunghezza (double x1, double y1, double z1, double x2, double y2, double z2);
