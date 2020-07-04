/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
//parameters, observables

using namespace std;

const int m_props=4;
int n_props;
int iv,ik,it,ie;
double stima_pot, stima_kin, stima_etot, stima_temp;

// averages
double acc,att;

//configuration
const int m_part=108;
double x[m_part],y[m_part],z[m_part],xold[m_part],yold[m_part],zold[m_part];
double vx[m_part],vy[m_part],vz[m_part];

// thermodynamical state
int npart;
double energy,temp,vol,rho,box,rcut;

// simulation
int nstep, iprint, seed;
double delta;
bool restart,analizza;
const int blocchi=100;
int L; //numero di misure che compongono un blocco

double media[blocchi][m_props]; //vettore che contiene <Q> (una componente per ogni blocco) Q:variabile termodinamica qualsiasi
double media_2[blocchi][m_props]; //contiene (<Q>)^2
double media_cumulativa[blocchi][m_props]; //la componente i-esima contiene la media tra i risultati dei blocchi fino all'i-esimo
double media_2_cumulativa[blocchi][m_props]; //come sopra ma al quadrato
double errore[blocchi][m_props]; //vettore che contiene le incertezze di ciascuna delle componenti media_cumulativa


int nmisura=1;
int nblocco=0;

//functions
void Input(void);
void Move(void);
void ConfFinal(void);
void ConfXYZ(int);
void Measure(void);
double Force(int, int);
double Pbc(double);
void old(void);
void calcolamedie (bool analizza, int argc, char **argv);
double calcolaerrore (double **x, double **y, int i, int j);
/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
