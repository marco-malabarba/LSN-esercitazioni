#ifndef __popolazione_h__
#define __popolazione_h__
#include <iostream>
#include <cstdlib>
#include <string>
#include "random.h"
#include "funzioni.h"

using namespace std;


class popolazione {

public:

popolazione(double *x, double *y, int n_cities, int n_popolazione);
~popolazione();
void setrndm(int k);
void creaindividui();
void stampaindividui();
void calcolalunghezze();
void calcolaprob();
void stampalunghezze();
void stampaprob();
void mutazioni();
void riproduzione();
double ave_best();
int cerca_best();
double best();
void sim_ann(double temp);
void stampapercorso(int q, string k);
void set_p_mut(double p1, double p2, double p3, double p4);
int getcities();
int getindividui();
int getpercorso(int x, int y);
void setpercorso (int x, int y, int z);

private:

int m_n_cities;
int m_n_individui;
double *m_x_city,*m_y_city;
int **m_individui;
double *m_lunghezze;
double *m_prob;
double p_mut1=0.10;
double p_mut2=0.05;
double p_mut3=0.03;
double p_mut4=0.03;
Random rnd;

};

#endif

