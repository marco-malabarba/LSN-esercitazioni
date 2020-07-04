/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#include <iostream>
#include <fstream>
#include <ostream>
#include <cmath>
#include <iomanip>
#include "Monte_Carlo_ISING_1D.h"
#include "random.cpp"

using namespace std;

int main()
{ 
  Input(); //Inizialization
  for(int iblk=1; iblk <= nblk; ++iblk) //Simulation
  {
    Reset(iblk);   //Reset block averages
    for(int istep=1; istep <= nstep; ++istep)
    {
      Move(metro);
      Measure(iblk);
      Accumulate(); //Update block averages
    }
    Averages(iblk);   //Print results for current block
  }
  ConfFinal(); //Write final configuration
  Datifinali();
  return 0;
}


void Input(void)
{
  ifstream ReadInput;

  cout << "Classic 1D Ising model             " << endl;
  cout << "Monte Carlo simulation             " << endl << endl;
  cout << "Nearest neighbour interaction      " << endl << endl;
  cout << "Boltzmann weight exp(- beta * H ), beta = 1/T " << endl << endl;
  cout << "The program uses k_B=1 and mu_B=1 units " << endl;

//Read seed for random numbers
   int p1, p2;
   ifstream Primes("Primes");
   Primes >> p1 >> p2 ;
   Primes.close();

   ifstream input("seed.in");
   input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
   rnd.SetRandom(seed,p1,p2);
   input.close();
  
//Read input informations
  ReadInput.open("input.dat");

  ReadInput >> temp;
  beta = 1.0/temp;
  cout << "Temperature = " << temp << endl;

  ReadInput >> nspin;
  cout << "Number of spins = " << nspin << endl;

  ReadInput >> J;
  cout << "Exchange interaction = " << J << endl;

  ReadInput >> h;
  cout << "External field = " << h << endl << endl;
    
  ReadInput >> metro; // if=1 Metropolis else Gibbs

  ReadInput >> nblk;

  ReadInput >> nstep;
  ReadInput >> restart;
  ReadInput >> analizza;
  ReadInput >> verif_eq;

  if(metro==1) cout << "The program perform Metropolis moves" << endl;
  else cout << "The program perform Gibbs moves" << endl;
  cout << "Number of blocks = " << nblk << endl;
  cout << "Number of steps in one block = " << nstep << endl << endl;
  ReadInput.close();


//Prepare arrays for measurements
  iu = 0; //Energy
  ic = 1; //Heat capacity
  im = 2; //Magnetization
  ix = 3; //Magnetic susceptibility
 
  n_props = 4; //Number of observables

//initial configuration
  if (restart==0) {
    for (int i=0; i<nspin; ++i)
    {
      if(rnd.Rannyu() >= 0.5) s[i] = 1;
      else s[i] = -1;
    }
  }

  else {
    ReadInput.open("config.final");
    for (int i=0; i<nspin; i++) ReadInput>>s[i];
  }
//Evaluate energy etc. of the initial configuration
  
  double u=0.;
  for (int i=0; i<nspin; ++i) u += -J * s[i] * s[Pbc(i+1)] - 0.5 * h * (s[i] + s[Pbc(i+1)]);

//Print initial values for the potential energy
  cout << "Initial energy = " << u/(double)nspin << endl;
}


void Move(int metro)
{
  int o;
  double p, energy_old, energy_new, sm;
  double energy_up, energy_down,p_up,p_down;

  for(int i=0; i<nspin; ++i)
  {
  //Select randomly a particle (for C++ syntax, 0 <= o <= nspin-1)
    o = (int)(rnd.Rannyu()*nspin);

    if(metro==1) //Metropolis
    {
    energy_old=Boltzmann(s[o],o);
    energy_new=Boltzmann(-s[o],o);
    double alfa=fmin(1.,exp(-1./temp*(energy_new-energy_old)));
    if (rnd.Rannyu()<=alfa) {
      flip(s[o]);
      accepted++;
      acc_tot++;
      }
    attempted++;
    att_tot++;
// INCLUDE YOUR CODE HERE // fatto
    }
    else //Gibbs sampling
    {
    energy_up=Boltzmann(1,o);
    energy_down=Boltzmann(-1,o);
    p_up=1./(1+exp(-1./temp*(energy_down-energy_up)));
    p_down=1./(1+exp(-1./temp*(energy_up-energy_down)));
    if (rnd.Rannyu()<p_up) s[o]=1;
    else s[o]=-1;
    attempted++;
    accepted++;
// INCLUDE YOUR CODE HERE // fatto
    }
  }
}

double Boltzmann(int sm, int ip)
{
  double ene = -J * sm * ( s[Pbc(ip-1)] + s[Pbc(ip+1)] ) - h * sm;
  return ene;
}

void Measure(int iblk)
{
  int bin;
  double u = 0.0, m = 0.0;
  double csi=0.0, M=0.0 ;

//cycle over spins
  for (int i=0; i<nspin; ++i)
  {
     
     u += -J * s[i] * s[Pbc(i+1)] - 0.5 * h * (s[i] + s[Pbc(i+1)]);
     csi+=s[i];
     M+=s[i];

// INCLUDE YOUR CODE HERE
  }
  csi=pow(csi,2);
  csi=csi*beta;
  walker[iu] = u;
  walker[ic]=u*u;
  walker[im]=M;
  walker[ix]=csi;
// INCLUDE YOUR CODE HERE //fatto!
  
//verifico se siamo all'equilibrio stampando l'energia istantanea per ogni passo Monte Carlo (solo del primo blocco!)

  if (verif_eq==1&&iblk==1) { 
    ofstream out("output.en_ist.0",ios::app);
    out<<(double) u/(double) nspin<<endl;
  }
}


void Reset(int iblk) //Reset block averages
{
   
   if(iblk == 1)
   {
       for(int i=0; i<n_props; ++i)
       {
           glob_av[i] = 0;
           glob_av2[i] = 0;
       }
       acc_tot=0.;
       att_tot=0.;
   }

   for(int i=0; i<n_props; ++i)
   {
     blk_av[i] = 0;
   }
   blk_norm = 0;
   attempted = 0;
   accepted = 0;
}


void Accumulate(void) //Update block averages
{

   for(int i=0; i<n_props; ++i)
   {
     blk_av[i] = blk_av[i] + walker[i];
   }
   blk_norm = blk_norm + 1.0;
}


void Averages(int iblk) //Print results for current block
{
    
   ofstream Ene, Heat, Mag, Chi;
    
    cout << "Block number " << iblk << endl;
    cout << "Acceptance rate " << accepted/attempted << endl << endl;
    
    Ene.open("output.ene.0",ios::app);
    stima_u = blk_av[iu]/(blk_norm-1)/(double)nspin; //Energy
    glob_av[iu]  += stima_u;
    glob_av2[iu] += stima_u*stima_u;
    err_u=Error(glob_av[iu],glob_av2[iu],iblk);
    Ene << iblk <<  " " << stima_u << " " << glob_av[iu]/(double)iblk << " " << err_u << endl;
    Ene.close();

    Heat.open("output.heat.0",ios::app);
    stima_c=blk_av[ic]/(double)(blk_norm-1)-pow(stima_u*nspin,2);
    stima_c=stima_c/nspin;  
    stima_c=beta*beta*stima_c;
    glob_av[ic]  += stima_c;
    glob_av2[ic] += stima_c*stima_c;
    err_c=Error(glob_av[ic],glob_av2[ic],iblk);
    Heat << iblk <<  " " << stima_c << " " << glob_av[ic]/(double)iblk << " " << err_c << endl;
    Heat.close();

    Mag.open("output.mag.0",ios::app);
    stima_m=blk_av[im]/(blk_norm-1)/(double)nspin;
    glob_av[im]  += stima_m;
    glob_av2[im] += stima_m*stima_m;
    err_m=Error(glob_av[im],glob_av2[im],iblk);
    Mag << iblk <<  " " << stima_m << " " << glob_av[im]/(double)iblk << " " << err_m << endl;
    Mag.close();

    Chi.open("output.chi.0",ios::app);
    stima_x=blk_av[ix]/(blk_norm-1)/(double)nspin;
    glob_av[ix]  += stima_x;
    glob_av2[ix] += stima_x*stima_x;
    err_x=Error(glob_av[ix],glob_av2[ix],iblk);
    Chi << iblk <<  " " << stima_x << " " << glob_av[ix]/(double)iblk << " " << err_x << endl;
    Chi.close();

// INCLUDE YOUR CODE HERE // fatto!

    cout << "----------------------------" << endl << endl;
}


void ConfFinal(void)
{
  ofstream WriteConf;

  cout << "Print final configuration to file config.final " << endl << endl;
  WriteConf.open("config.final");
  for (int i=0; i<nspin; ++i)
  {
    WriteConf << s[i] << endl;
  }
  WriteConf.close();

  rnd.SaveSeed();
}

int Pbc(int i)  //Algorithm for periodic boundary conditions
{
    if(i >= nspin) i = i - nspin;
    else if(i < 0) i = i + nspin;
    return i;
}

double Error(double sum, double sum2, int iblk)
{
    if (iblk-1==0) return 0;
    
    return sqrt((sum2/(double)iblk - pow(sum/(double)iblk,2))/(double)(iblk-1));
}

void flip (double &p) {
    p=-p;
    return;
}

void Datifinali (void) {
  if (analizza==0) return;
  if (h==0) {
     string s="ene";
     if (metro==1) s=s+(string)"_metro";
     else s=s+(string)"_gibbs";
     s=s+".dat";
     ofstream out(s,ios::app);
     out<<temp<<" "<<glob_av[iu]/(double)nblk << " " << err_u << endl;
     out.close();
     s="heat";
     if (metro==1) s=s+(string)"_metro";
     else s=s+(string)"_gibbs";
     s=s+".dat";
     out.open(s,ios::app);
     out<<temp<<" "<<glob_av[ic]/(double)nblk << " " << err_c << endl;
     out.close();
     s="chi";
     if (metro==1) s=s+(string)"_metro";
     else s=s+(string)"_gibbs";
     s=s+".dat";
     out.open(s,ios::app);
     out<<temp<<" "<<glob_av[ix]/(double)nblk << " " << err_x << endl;
     out.close();
     if (metro==1) {
        out.open("acc.dat",ios::app);
        out<<temp<<" "<<(double)acc_tot/(double) att_tot<<endl;
     }
     return;
  }
  else {
     string s="mag";
     if (metro==1) s=s+(string)"_metro";
     else s=s+(string)"_gibbs";
     s=s+".dat";
     ofstream out(s,ios::app);
     out<<temp<<" "<<glob_av[im]/(double)nblk << " " << err_m << endl;
     out.close();
     return;
  }
}



/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
