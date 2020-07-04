/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
#include <stdlib.h>     // srand, rand: to generate random number
#include <iostream>     // cin, cout: Standard Input/Output Streams Library
#include <fstream>      // Stream class to both read and write from/to files.
#include <cmath>        // rint, pow
#include "MolDyn_NVE.h"

using namespace std;

int main(int argc, char **argv){ 
 //for (int i=0; i<blocchi; i++) {
    //for (int j=0; j<m_props; j++) {
      //media[i][j]=0;
      //media_2[i][j]=0;
      //media_cumulativa[i][j]=0;
      //media_2_cumulativa[i][j]=0;   
    //}
  //}

  Input();             //Inizialization
  int nconf = 1;
  for(int istep=1; istep <= nstep; ++istep){
     Move();           //Move particles with Verlet algorithm
     if(istep%iprint == 0) cout << "Number of time-steps: " << istep << endl;
     if(istep%1 == 0){
        Measure();     //Properties measurement
        //ConfXYZ(nconf);//Write actual configuration in XYZ format //Commented to avoid "filesystem full"! 
        nconf += 1;
     }
  if (istep==nstep-1) old();
  }
  
  
  calcolamedie(analizza,argc,argv);
  ConfFinal();         //Write final configuration to restart
  return 0;
}


void Input(void){ //Prepare all stuff for the simulation
  ifstream ReadInput,ReadConf;
  double ep, ek, pr, et, vir;

  cout << "Classic Lennard-Jones fluid        " << endl;
  cout << "Molecular dynamics simulation in NVE ensemble  " << endl << endl;
  cout << "Interatomic potential v(r) = 4 * [(1/r)^12 - (1/r)^6]" << endl << endl;
  cout << "The program uses Lennard-Jones units " << endl;

  seed = 1;    //Set seed for random numbers
  srand(seed); //Initialize random number generator
  
  ReadInput.open("input_md.dat"); //Read input

  ReadInput >> temp;

  ReadInput >> npart;
  cout << "Number of particles = " << npart << endl;

  ReadInput >> rho;
  cout << "Density of particles = " << rho << endl;
  vol = (double)npart/rho;
  cout << "Volume of the simulation box = " << vol << endl;
  box = pow(vol,1.0/3.0);
  cout << "Edge of the simulation box = " << box << endl;

  ReadInput >> rcut;
  ReadInput >> delta;
  ReadInput >> nstep;
  ReadInput >> iprint;
  ReadInput >> restart;
  ReadInput >> analizza;
  L=nstep/(blocchi); //faccio una misura solo ogni 10 step
  cout << "The program integrates Newton equations with the Verlet method " << endl;
  cout << "Time step = " << delta << endl;
  cout << "Number of steps = " << nstep << endl << endl;
  ReadInput.close();

//Prepare array for measurements
  iv = 0; //Potential energy
  ik = 1; //Kinetic energy
  ie = 2; //Total energy
  it = 3; //Temperature
  n_props = 4; //Number of observables
  
  bin_size = (box/2.0)/(double)nbin;

if (restart==false) {
//Read initial configuration
  cout << "Read initial configuration from file config.0 " << endl << endl;
  ReadConf.open("config.0");
  for (int i=0; i<npart; ++i){
    ReadConf >> x[i] >> y[i] >> z[i];
    x[i] = x[i] * box;
    y[i] = y[i] * box;
    z[i] = z[i] * box;
  }
  ReadConf.close();

//Prepare initial velocities
   cout << "Prepare random velocities with center of mass velocity equal to zero " << endl << endl;
   double sumv[3] = {0.0, 0.0, 0.0};
   for (int i=0; i<npart; ++i){
     vx[i] = rand()/double(RAND_MAX) - 0.5;
     vy[i] = rand()/double(RAND_MAX) - 0.5;
     vz[i] = rand()/double(RAND_MAX) - 0.5;

     sumv[0] += vx[i];
     sumv[1] += vy[i];
     sumv[2] += vz[i];
   }
   for (int idim=0; idim<3; ++idim) sumv[idim] /= (double)npart;
   double sumv2 = 0.0, fs;
   for (int i=0; i<npart; ++i){
     vx[i] = vx[i] - sumv[0];
     vy[i] = vy[i] - sumv[1];
     vz[i] = vz[i] - sumv[2];

     sumv2 += vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i];
   }
   sumv2 /= (double)npart;

   fs = sqrt(3 * temp / sumv2);   // fs = velocity scale factor 
   for (int i=0; i<npart; ++i){
     vx[i] *= fs;
     vy[i] *= fs;
     vz[i] *= fs;

     xold[i] = Pbc(x[i] - vx[i] * delta);
     yold[i] = Pbc(y[i] - vy[i] * delta);
     zold[i] = Pbc(z[i] - vz[i] * delta);
   }
   return;
}

if (restart==true && temp>0) {
  cout << "Read initial and previous configurations from files old.final and old.0" << endl << endl;
  ReadConf.open("old.final");
  for (int i=0; i<npart; ++i){
    ReadConf >> x[i] >> y[i] >> z[i];
    x[i] = x[i] * box;
    y[i] = y[i] * box;
    z[i] = z[i] * box;
  }
  ReadConf.close();  

  ReadConf.open("old.0");
  for (int i=0; i<npart; ++i){
    ReadConf >> xold[i] >> yold[i] >> zold[i];
    xold[i] = xold[i] * box;
    yold[i] = yold[i] * box;
    zold[i] = zold[i] * box;
  }
  ReadConf.close();  

  double xnew, ynew, znew, fx[m_part], fy[m_part], fz[m_part];

  for(int i=0; i<npart; ++i){ //Force acting on particle i
    fx[i] = Force(i,0);
    fy[i] = Force(i,1);
    fz[i] = Force(i,2);
  }

  for(int i=0; i<npart; ++i){ //Verlet integration scheme

    xnew = Pbc( 2.0 * x[i] - xold[i] + fx[i] * pow(delta,2) );
    ynew = Pbc( 2.0 * y[i] - yold[i] + fy[i] * pow(delta,2) );
    znew = Pbc( 2.0 * z[i] - zold[i] + fz[i] * pow(delta,2) );
    vx[i] = Pbc(xnew - x[i])/(delta);
    vy[i] = Pbc(ynew - y[i])/(delta);
    vz[i] = Pbc(znew - z[i])/(delta);
    xold[i] = x[i];
    yold[i] = y[i];
    zold[i] = z[i];

    x[i] = xnew;
    y[i] = ynew;
    z[i] = znew;
  }
  double t = 0.0;
  for (int i=0; i<npart; ++i) t += 0.5 * (vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i]);
  double fattorescala=temp/(2./3.*t/(double)npart);
  for(int i=0; i<npart; ++i){
    vx[i]*=sqrt(fattorescala);
    vy[i]*=sqrt(fattorescala);
    vz[i]*=sqrt(fattorescala);
    xold[i]=Pbc(x[i]-delta*vx[i]);
    yold[i]=Pbc(y[i]-delta*vy[i]);
    zold[i]=Pbc(z[i]-delta*vz[i]);
  }
  
return;
}

if (restart==true && temp<=0) {
  cout << "REad initial and previous configurations from files old.final and old.0" << endl << endl;
  ReadConf.open("old.final");
  for (int i=0; i<npart; ++i){
    ReadConf >> x[i] >> y[i] >> z[i];
    x[i] = x[i] * box;
    y[i] = y[i] * box;
    z[i] = z[i] * box;
  }
  ReadConf.close();  

  ReadConf.open("old.0");
  for (int i=0; i<npart; ++i){
    ReadConf >> xold[i] >> yold[i] >> zold[i];
    xold[i] = xold[i] * box;
    yold[i] = yold[i] * box;
    zold[i] = zold[i] * box;
  }
  ReadConf.close();
  double xnew, ynew, znew, fx[m_part], fy[m_part], fz[m_part];

  for(int i=0; i<npart; ++i){ //Force acting on particle i
    fx[i] = Force(i,0);
    fy[i] = Force(i,1);
    fz[i] = Force(i,2);
  }

  for(int i=0; i<npart; ++i){ //Verlet integration scheme

    xnew = Pbc( 2.0 * x[i] - xold[i] + fx[i] * pow(delta,2) );
    ynew = Pbc( 2.0 * y[i] - yold[i] + fy[i] * pow(delta,2) );
    znew = Pbc( 2.0 * z[i] - zold[i] + fz[i] * pow(delta,2) );
    vx[i] = Pbc(xnew - x[i])/(delta);
    vy[i] = Pbc(ynew - y[i])/(delta);
    vz[i] = Pbc(znew - z[i])/(delta);
    xold[i] = x[i];
    yold[i] = y[i];
    zold[i] = z[i];

    x[i] = xnew;
    y[i] = ynew;
    z[i] = znew;
  }
return;
}  

cout<<"Errore!!"<<endl;
exit(1);
}


void Move(void){ //Move particles with Verlet algorithm
  double xnew, ynew, znew, fx[m_part], fy[m_part], fz[m_part];

  for(int i=0; i<npart; ++i){ //Force acting on particle i
    fx[i] = Force(i,0);
    fy[i] = Force(i,1);
    fz[i] = Force(i,2);
  }

  for(int i=0; i<npart; ++i){ //Verlet integration scheme

    xnew = Pbc( 2.0 * x[i] - xold[i] + fx[i] * pow(delta,2) );
    ynew = Pbc( 2.0 * y[i] - yold[i] + fy[i] * pow(delta,2) );
    znew = Pbc( 2.0 * z[i] - zold[i] + fz[i] * pow(delta,2) );

    vx[i] = Pbc(xnew - xold[i])/(2.0 * delta);
    vy[i] = Pbc(ynew - yold[i])/(2.0 * delta);
    vz[i] = Pbc(znew - zold[i])/(2.0 * delta);

    xold[i] = x[i];
    yold[i] = y[i];
    zold[i] = z[i];

    x[i] = xnew;
    y[i] = ynew;
    z[i] = znew;
  }
  return;
}

double Force(int ip, int idir){ //Compute forces as -Grad_ip V(r)
  double f=0.0;
  double dvec[3], dr;

  for (int i=0; i<npart; ++i){
    if(i != ip){
      dvec[0] = Pbc( x[ip] - x[i] );  // distance ip-i in pbc
      dvec[1] = Pbc( y[ip] - y[i] );
      dvec[2] = Pbc( z[ip] - z[i] );

      dr = dvec[0]*dvec[0] + dvec[1]*dvec[1] + dvec[2]*dvec[2];
      dr = sqrt(dr);

      if(dr < rcut){
        f += dvec[idir] * (48.0/pow(dr,14) - 24.0/pow(dr,8)); // -Grad_ip V(r)
      }
    }
  }
  
  return f;
}

void Measure(){ //Properties measurement
  int bin;
  double v, t, vij;
  double dx, dy, dz, dr;
  //ofstream Epot, Ekin, Etot, Temp;

  //Epot.open("output_epot.dat",ios::app);
  //Ekin.open("output_ekin.dat",ios::app);
  //Temp.open("output_temp.dat",ios::app);
  //Etot.open("output_etot.dat",ios::app);

  //v = 0.0; //reset observables
  //t = 0.0;

//cycle over pairs of particles
  for (int i=0; i<nbin; i++) walker[i]=0.0;
  for (int i=0; i<npart-1; ++i){
    for (int j=i+1; j<npart; ++j){

     dx = Pbc( xold[i] - xold[j] ); // here I use old configurations [old = r(t)]
     dy = Pbc( yold[i] - yold[j] ); // to be compatible with EKin which uses v(t)
     dz = Pbc( zold[i] - zold[j] ); // => EPot should be computed with r(t)

     dr = dx*dx + dy*dy + dz*dz;
     dr = sqrt(dr);
     bin=(int) (dr/(double)(bin_size));
     if (bin<100) walker[bin]+=2.;

     //if(dr < rcut){
       //vij = 4.0/pow(dr,12) - 4.0/pow(dr,6);

//Potential energy
       //v += vij;
     //}
    }          
  }
//Kinetic energy
  //for (int i=0; i<npart; ++i) t += 0.5 * (vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i]);
   
    //stima_pot = v/(double)npart; //Potential energy per particle
    //stima_kin = t/(double)npart; //Kinetic energy per particle
    //stima_temp = (2.0 / 3.0) * t/(double)npart; //Temperature
    //stima_etot = (t+v)/(double)npart; //Total energy per particle

    //Epot << stima_pot  << endl;
    //media[nblocco][0]+=stima_pot;
    //Ekin << stima_kin  << endl;
    //media[nblocco][1]+=stima_kin;
    //Temp << stima_temp << endl;
    //media[nblocco][2]+=stima_temp;
    //Etot << stima_etot << endl;
    //media[nblocco][3]+=stima_etot;
    //if (nmisura%L==0) {
      //media[nblocco][0]=media[nblocco][0]/L;
      //media[nblocco][1]=media[nblocco][1]/L;
      //media[nblocco][2]=media[nblocco][2]/L;
      //media[nblocco][3]=media[nblocco][3]/L;
      //nblocco++;
    //}
    //Epot.close();
    //Ekin.close();
    //Temp.close();
    //Etot.close();
    //nmisura++;
  for (int i=0; i<nbin; i++) {
    double valore=walker[i]/(rho*npart*4./3.*M_PI*(pow((bin_size)*(i+1.),3.)-pow(bin_size*(i),3.)));
    glob_av[i]+=valore;
    glob_av2[i]+=valore*valore;
  }
   
  return;
}


void ConfFinal(void){ //Write final configuration
  ofstream WriteConf,writeold;

  cout << "Print final configuration to files config.final and old.final" << endl << endl;
  WriteConf.open("config_md.final");
  writeold.open("old.final");

  for (int i=0; i<npart; ++i){
    WriteConf << x[i]/box << "   " <<  y[i]/box << "   " << z[i]/box << endl;
    writeold << x[i]/box << "   " <<  y[i]/box << "   " << z[i]/box << endl;
  }
  WriteConf.close();
  writeold.close();
  return;
}

void old (void) {
  ofstream out;
  out.open("old.0");
  for (int i=0; i<npart; ++i){
    out << x[i]/box << "   " <<  y[i]/box << "   " << z[i]/box << endl;
  }
  out.close();
  return;
}

void ConfXYZ(int nconf){ //Write configuration in .xyz format
  ofstream WriteXYZ;

  WriteXYZ.open("frames/config_" + to_string(nconf) + ".xyz");
  WriteXYZ << npart << endl;
  WriteXYZ << "This is only a comment!" << endl;
  for (int i=0; i<npart; ++i){
    WriteXYZ << "LJ  " << Pbc(x[i]) << "   " <<  Pbc(y[i]) << "   " << Pbc(z[i]) << endl;
  }
  WriteXYZ.close();
}

double Pbc(double r){  //Algorithm for periodic boundary conditions with side L=box
    return r - box * rint(r/box);
}


void calcolamedie (bool analizza, int argc, char **argv) {
  if (analizza==false) return;
  //for (int i=0; i<blocchi; i++) {
    //for (int j=0; j<m_props; j++) {
       //media_2[i][j]=pow(media[i][j],2);
    //}
  //}
  //for (int k=0; k<m_props; k++) {
    //for (int i=0; i<blocchi; i++) { //ciclo sui blocchi
      //for (int j=0; j<i+1; j++) { //ciclo che mi permette di sommare fino all'i-esimo blocco (sto cumulando) 
        //media_cumulativa[i][k]+=media[j][k];
        //media_2_cumulativa[i][k]+=media_2[j][k];
      //}
      //media_cumulativa[i][k]=media_cumulativa[i][k]/(i+1); //media delle medie: devo quindi dividere per il numero di componenti  
      //media_2_cumulativa[i][k]=media_2_cumulativa[i][k]/(i+1);
      //if (i==0) {
	//errore[i][k]=0; //se faccio solo una misura (un singolo blocco) non ha senso definire un'incertezza
      //}
      //else errore[i][k]=sqrt((media_2_cumulativa[i][k]-pow(media_cumulativa[i][k],2))/i); //formula per l'incertezza
    //}
  //}
  /*string title;
  ofstream out;
  if (argc==1) {
    title=string("average_epot") + string(".dat");
    out.open(title);
    for (int i=0; i<blocchi; i++) {
      out<<L*(i+1)<<" "<<media_cumulativa[i][0]<<" "<<errore[i][0]<<endl;
    }
    out.close();
    title=string("average_ekin") + string(".dat");
    out.open(title);
    for (int i=0; i<blocchi; i++) {
      out<<L*(i+1)<<" "<<media_cumulativa[i][1]<<" "<<errore[i][1]<<endl;
    }
    out.close(); 
    title=string("average_temp") + string(".dat");
    out.open(title);
    for (int i=0; i<blocchi; i++) {
      out<<L*(i+1)<<" "<<media_cumulativa[i][2]<<" "<<errore[i][2]<<endl;
    }
    out.close();
    title=string("average_etot") + string(".dat");
    out.open(title);
    for (int i=0; i<blocchi; i++) {
      out<<L*(i+1)<<" "<<media_cumulativa[i][3]<<" "<<errore[i][3]<<endl;
    }
    out.close();
  }
  if (argc>1) {
    title=string("average_epot") + string("_") + string(argv[1])+ string(".dat");
    out.open(title);
    for (int i=0; i<blocchi; i++) {
      out<<L*(i+1)<<" "<<media_cumulativa[i][0]<<" "<<errore[i][0]<<endl;
    }
    out.close();
    title=string("average_ekin") + string("_") + string(argv[1])+ string(".dat");
    out.open(title);
    for (int i=0; i<blocchi; i++) {
      out<<L*(i+1)<<" "<<media_cumulativa[i][1]<<" "<<errore[i][1]<<endl;
    }
    out.close(); 
    title=string("average_temp") + string("_") + string(argv[1])+ string(".dat");
    out.open(title);
    for (int i=0; i<blocchi; i++) {
      out<<L*(i+1)<<" "<<media_cumulativa[i][2]<<" "<<errore[i][2]<<endl;
    }
    out.close();
    title=string("average_etot") + string("_") + string(argv[1])+ string(".dat");
    out.open(title);
    for (int i=0; i<blocchi; i++) {
      out<<L*(i+1)<<" "<<media_cumulativa[i][3]<<" "<<errore[i][3]<<endl;
    }
    out.close();
  }*/
  
  ofstream out ("output.gave_md.0");
  for (int i=0; i<nbin; i++) {
    glob_av[i]=glob_av[i]/L
    glob_av2[i]=glob_av2[i]/L*L
    errorgdir=Error(glob_av[i],glob_av2[i],blocchi);
    out<<i*bin_size<<" "<<glob_av[i]/(double)(blocchi)<<" "<<errorgdir<<endl;
  }
  out.close();
return;
}
 

double Error(double sum, double sum2, int iblk)
{
    if( iblk == 1 ) return 0.0;
    else if (((sum2/(double)iblk - pow(sum/(double)iblk,2))/(double)(iblk-1))<0) return -1;
    return sqrt((sum2/(double)iblk - pow(sum/(double)iblk,2))/(double)(iblk-1));
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
