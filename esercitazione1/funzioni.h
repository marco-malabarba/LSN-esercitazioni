#include <cmath>
#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include "random.h"

using namespace std;

void setrandom(Random &rnd);
double calcolaerrore (double *x, double *y, int i);
void calcolamedie (double *media, double *media_2, double *x_random, int N, int L);
void calcolacumulate(double *media, double *media_2, double *media_cumulativa, double *media_2_cumulativa, double *errore, int N);
void azzeraarray(double *x, int n);
void calcolavarianze (double *media, double *media_2, double *x_random, int N, int L);
void calcolachiquadro (Random rnd,double *x_chiquadro,double *x_i,double *chiquadro,int n , int m_int, int m_test);
void creafiledatiunif(Random rnd, int M, int N, ofstream &out);
void creafiledatiexp(Random rnd, int M, int N, ofstream &out);
void creafiledatilor(Random rnd, int M, int N, ofstream &out);
void creafiledati(Random rnd, int M, int N, ofstream &out);
void creafiledaticorr(Random rnd, int M, int N, ofstream &out);
void calcolapi (double *x_ago,double *angolo_ago, double *media, double *media_2, double l,double d,double N,double L);
double lunghezza (double x1, double y1, double x2, double y2);
