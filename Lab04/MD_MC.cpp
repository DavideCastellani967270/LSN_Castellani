/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/


/*
Author: Davide Castellani
Email: davide.castellani1@studenti.unimi.it
Matr.: 967270


EXERCISE 04

AIM OF THE CODE:
This code's purpuse is to simulate the dynamics of a classical N-body system, modelling the two-particle interaction with the adiabatic Lennard-Jones potential.
Given a starting crystalline structure to avoid superposition and normally distributed random starting velocities, the system is evolved through time in
two possible methods: by computing the forces acting on each particle and individually updating their position (Molecular Dynamics, MD) or by using a
Metropolis agorithm to propose a new step (Monte Carlo, MC). While doing so, multiple thermodynamical quantities are measured (potential energy, kinetic
energy, total energy, temperature and pressure) in order to compute their progressive average values and uncertainties as a functions of the increasing 
number of blocks.
In particular it will be simulated the three states of matter (solid, liquid and gas) of nobel gas argon (Ar) at a fixed density, temperature and cutoff radius.
Note that is expected a steep drop in temperature during the initial few iterations due to the low-entropy crystalline starting configuration, so equilibration
time to get the system's temperature to the desired value is needed.
  

KEY WORDS:
  - Statistical mechanics
  - Molecular dynamics
  - Integration algorithms: the Verlet's algorithm
  - Periodic boundary conditions
  - Lennard-Jones fluid
  - Statistical uncertainties in Monte Carlo simulations: the blocking method
*/

#include <iostream>
#include <fstream>
#include <ostream>
#include <cmath>
#include <iomanip>
#include "MD_MC.h"

using namespace std;

int main(){

  Input();                                        //inizialization
  int nconf = 1;
  for(int iblk=1; iblk <= nblk+1; iblk++){        //cycles over the nblk+1 blocks: the +1 is for equilibration as here a full block is skipped before taking measurements
    Reset(iblk);                                  //resets block averages
    
    for(int istep=1; istep <= nstep; istep++){    //cycles over the steps in each block
      Move();                                     //updates the positions
      Measure();                                  //measures the observables
      Accumulate();                               //updates the block averages
      if(istep%10 == 0){
        //ConfXYZ(nconf);                         //writes actual configuration in XYZ format: commented to avoid "filesystem full"! 
        nconf += 1;
      }
    }

    if(iblk != 1)             //skips the first block in order to let the system equilibrate
      Averages(iblk-1);       //prints the results for current block
  }
  ConfFinal();                //writes the final configuration

  return 0;
}


//this method initialises the system with the starting configurations and velocities and then prints the initial observables values
void Input(void){
  ifstream ReadInput, ReadConf, ReadVelocity, Primes, Seed;

  //prints some informations about the code
  cout << "Classic Lennard-Jones fluid        " << endl;
  cout << "MD(NVE) / MC(NVT) simulation       " << endl << endl;
  cout << "Interatomic potential v(r) = 4 * [(1/r)^12 - (1/r)^6]" << endl << endl;
  cout << "Boltzmann weight exp(- beta * sum_{i<j} v(r_ij) ), beta = 1/T " << endl << endl;
  cout << "The program uses Lennard-Jones units " << endl;

  //reads two primes to set up the random numbers generator (RNG)
  int p1, p2;
  Primes.open("Primes");
  Primes >> p1 >> p2 ;
  Primes.close();

  //read input informations
  ReadInput.open("input.in");

  ReadInput >> iNVET;       // 0 = molecular dynamics (MD), 1 = Monte Carlo (MC)  ->  reads if the user wants to do a MD or MC simulation
  ReadInput >> restart;     // 0 = no restart, 1 = restart                        ->  reads if the user wants to restart the simulation from a previous run

  if(restart)
    Seed.open("seed.out");                                //opens the seed.out file in order to set the RNG to the exact configuration of the last run
  else 
    Seed.open("seed.in");                                 //opens the seed.in file in order to set the RNG for a new run
  Seed >> seed[0] >> seed[1] >> seed[2] >> seed[3];
  rnd.SetRandom(seed,p1,p2);                              //sets up the RNG
  Seed.close();

  ReadInput >> temp;                                      //reads the INITIAL temperature (in reduced units), prints its value and calculates the beta factor
  beta = 1.0/temp;
  cout << "Initial temperature = " << temp << endl;

  ReadInput >> npart;                                     //reads the number of particles and prints it
  cout << "Number of particles = " << npart << endl;

  ReadInput >> rho;                                       //reads the number density of particles and prints it
  cout << "Density of particles = " << rho << endl;
  vol = (double)npart/rho;                                //computes the volume of the simulation box and the lenght of its edge
  box = pow(vol,1.0/3.0);
  cout << "Volume of the simulation box = " << vol << endl;
  cout << "Edge of the simulation box = " << box << endl;

  ReadInput >> rcut;                                      //reads the cutoff radius of the interatomic interaction and prints it
  cout << "Cutoff of the interatomic potential = " << rcut << endl << endl;
    
  ReadInput >> delta;                                     //reads and prints the moving time parameter (in reduced units, dt^* = dt sqrt( epsilon/ m*sigma^2 ))

  ReadInput >> nblk;                                      //reads and prints the number of blocks

  ReadInput >> nstep;                                     //reads and prints the number of steps in each block

  cout << "The program perform Metropolis moves with uniform translations" << endl;
  cout << "Moves parameter = " << delta << endl;
  cout << "Number of blocks = " << nblk << endl;
  cout << "Number of steps in one block = " << nstep << endl << endl;
  ReadInput.close();

  //prepare arrays for measurements
  iv = 0;      //array position referring to potential energy
  it = 1;      //array position referring to temperature
  ik = 2;      //array position referring to kinetic energy
  ie = 3;      //array position referring to total energy
  ip = 4;      //array position referring to pressure
  n_props = 5; //total number of observables

  //read initial configuration and velocities
  cout << "Read initial configuration" << endl << endl;
  if(restart) {
    ReadConf.open("config.out");                  //opens the file containing the final configuration of the previous run
    ReadVelocity.open("velocity.out");            //opens the file containing the final velocities of the previous run

    for (int i=0; i<npart; ++i)
      ReadVelocity >> vx[i] >> vy[i] >> vz[i];    //reads the old velocities
  
  } else {
    ReadConf.open("config.in");                   //opens the file containing the initial configuration (crystal structure, to avoid superpositions)

    cout << "Prepare velocities with center of mass velocity equal to zero" << endl;
    double sumv[3] = {0.0, 0.0, 0.0};             //array that will contain the drift velocities (mean velocities in each direction)

    for (int i=0; i<npart; ++i) {
      vx[i] = rnd.Gauss(0., sqrt(temp));          //gives a normally distributed velocities accordingly to Maxwell/Boltzmann statistics (sigma = sqrt(T))
      vy[i] = rnd.Gauss(0., sqrt(temp));
      vz[i] = rnd.Gauss(0., sqrt(temp));

      sumv[0] += vx[i];                   //adds the ith particle velocities components to previous ones
      sumv[1] += vy[i];
      sumv[2] += vz[i];
    }

    for (int idim=0; idim<3; ++idim)
      sumv[idim] /= (double)npart;        //computes the average velocity in each component (drift velocities)
    
    double sumv2 = 0.0, fs;
    for (int i=0; i<npart; ++i){
      vx[i] = vx[i] - sumv[0];            //subtracts the drift velocities for each particle to get the velocities related to the system's center of mass
      vy[i] = vy[i] - sumv[1];
      vz[i] = vz[i] - sumv[2];

      sumv2 += vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i];       //computes the squared modulus of the ith particle velocity and adds it to the ones of the previous particles
    }

    sumv2 /= (double)npart;               //computes < |v|^2 >
    fs = sqrt(3 * temp / sumv2);          //and the velocity scale factor fs
    cout << "velocity scale factor: " << fs << endl << endl;

    for (int i=0; i<npart; ++i){
      vx[i] *= fs;                        //scales the velocities of each particle
      vy[i] *= fs;
      vz[i] *= fs;
    }
  }

  for (int i=0; i<npart; ++i) {
    ReadConf >> x[i] >> y[i] >> z[i];     //reads the configuration (in box side units) from the previously selected file
    x[i] = Pbc( x[i] * box );             //scales the configuration (in reduced units) while checking the periodic boundary conditions
    y[i] = Pbc( y[i] * box );
    z[i] = Pbc( z[i] * box );
  }

  ReadConf.close();

  //find the old configuration from the initial one: this is an essential information when making a new move, evolving the system
  for (int i=0; i<npart; ++i) {
    if(iNVET) {         //in MC, to draw the next configuration only the present one is needed, so here it simply copies the starting positions
      xold[i] = x[i];
      yold[i] = y[i];
      zold[i] = z[i];
    }
    else                //in MD, the integration scheme used (the Verlet's algorithm) requires the notion of both the present and the previous positions
    {
      xold[i] = Pbc(x[i] - vx[i] * delta);    //compute old configuration from a linear uniform motion assumption (and checks PBCs)
      yold[i] = Pbc(y[i] - vy[i] * delta);
      zold[i] = Pbc(z[i] - vz[i] * delta);
    }
  }
  
  //evaluates the properties of the initial configuration
  Measure();

  //prints the initial values for measured properties
  cout << "Initial potential energy = " << walker[iv]/(double)npart << endl;
  cout << "Initial temperature      = " << walker[it] << endl;
  cout << "Initial kinetic energy   = " << walker[ik]/(double)npart << endl;
  cout << "Initial total energy     = " << walker[ie]/(double)npart << endl;
  cout << "Initial pressure         = " << walker[ip] << endl << endl;

  return;
}


//this method evolves the present configuration into a new one (not necessarily different in MC)
void Move(){
  int o;
  double p, energy_old, energy_new;
  double xnew, ynew, znew;

  if(iNVET) {   //Monte Carlo (NVT) move

    for(int i=0; i<npart; ++i){
      
      o = (int)(rnd.Rannyu() * npart);            //selects randomly a particle (for C++ syntax, 0 <= o <= npart-1)

      energy_old = Boltzmann(x[o],y[o],z[o],o);   //computes the energy of the oth particle in the old configuration

      x[o] = Pbc( x[o] + delta*(rnd.Rannyu() - 0.5) );    //computes the energy that the oth particle would have its new proposed position (after checking PBCs)
      y[o] = Pbc( y[o] + delta*(rnd.Rannyu() - 0.5) );
      z[o] = Pbc( z[o] + delta*(rnd.Rannyu() - 0.5) );

      energy_new = Boltzmann(x[o],y[o],z[o],o);

      //Metropolis algorithm's test
      p = exp(beta * (energy_old-energy_new));      //Boltzmann's weight is used to get the acceptance probability
      if(p >= rnd.Rannyu()) {

        xold[o] = x[o];                             //updates the old position with the proposed one
        yold[o] = y[o];
        zold[o] = z[o];

        accepted = accepted + 1.0;                  //keeps track of the acceptance rate of the algorithm
      } else {

        x[o] = xold[o];                             //the proposed move hasn't been accepted, so the old position is overwritten
        y[o] = yold[o];
        z[o] = zold[o];
      }
      attempted = attempted + 1.0;
    }
  } else    //Molecular Dynamics (NVE) move
  {
    double fx[m_part], fy[m_part], fz[m_part];      //vectors containing the forces in xyz directions for every particle (m_part = 100 is the maximum number of particles)

    for(int i=0; i<npart; ++i){     //computes the components of the force acting on particle i
      fx[i] = Force(i,0);
      fy[i] = Force(i,1);
      fz[i] = Force(i,2);
    }

    //Verlet integration scheme
    for(int i=0; i<npart; ++i){ 

      xnew = Pbc( 2.0 * x[i] - xold[i] + fx[i] * pow(delta,2) );
      ynew = Pbc( 2.0 * y[i] - yold[i] + fy[i] * pow(delta,2) );
      znew = Pbc( 2.0 * z[i] - zold[i] + fz[i] * pow(delta,2) );

      vx[i] = Pbc(xnew - xold[i])/(2.0 * delta);
      vy[i] = Pbc(ynew - yold[i])/(2.0 * delta);
      vz[i] = Pbc(znew - zold[i])/(2.0 * delta);

      xold[i] = x[i];       //replaces the old positions with the present ones
      yold[i] = y[i];
      zold[i] = z[i];

      x[i] = xnew;          //replaces the present positions with the new ones
      y[i] = ynew;
      z[i] = znew;

      accepted = accepted + 1.0;          //this method doesn't propose moves, so the acceptance rate must be 1
      attempted = attempted + 1.0;
    }
  }
  return;
}


//this method computes the energy of the particle ip in a LJ potential within a cutoff radius
double Boltzmann(double xx, double yy, double zz, int ip) {

  double ene=0.0;
  double dx, dy, dz, dr;

  for (int i=0; i<npart; ++i){            //cycles over every particle but the one whose energy is being calculated
    if(i != ip){

      dx = Pbc(xx - x[i]);                //computes the distance ip-i in PBC
      dy = Pbc(yy - y[i]);
      dz = Pbc(zz - z[i]);

      dr = dx*dx + dy*dy + dz*dz;
      dr = sqrt(dr);

      if(dr < rcut) {                             //checks if the particle i is within the cutoff radius
        ene += 1.0/pow(dr,12) - 1.0/pow(dr,6);    //adds the energy to the previously calculated
      }
    }
  }

  return 4.0*ene;
}


//this method computes the force in direction idir (x=0, y=1, z=2) acting on the ip particle as -Grad_ip V(r)
double Force(int ip, int idir){

  double f=0.0;
  double dvec[3], dr;                     //dvec = vector of the distances

  for (int i=0; i<npart; ++i){            //cycles over every particle but the one the force that is being calculated acts on
    if(i != ip){
      dvec[0] = Pbc( x[ip] - x[i] );      //computes the distance ip-i in PBCs
      dvec[1] = Pbc( y[ip] - y[i] );
      dvec[2] = Pbc( z[ip] - z[i] );

      dr = dvec[0]*dvec[0] + dvec[1]*dvec[1] + dvec[2]*dvec[2];
      dr = sqrt(dr);

      if(dr < rcut){                                                //checks if the particle i is within the cutoff radius
        f += dvec[idir] * (48.0/pow(dr,14) - 24.0/pow(dr,8));       //adds -Grad_ip V(r) to the previously calculated 
      }
    }
  }
 
  return f;
}


//This method measures the thermodynamical properties in a given configuration
void Measure() {

  double v = 0.0, kin=0.0, p=0.0;         //support variables
  double dx, dy, dz, dr;

  for (int i=0; i<npart-1; ++i){          //cycles over every possible pair of particles
    for (int j=i+1; j<npart; ++j){ 

      dx = Pbc(x[i] - x[j]);              //computes the distance between the particles i and j in PBCs
      dy = Pbc(y[i] - y[j]);
      dz = Pbc(z[i] - z[j]);

      dr = dx*dx + dy*dy + dz*dz;
      dr = sqrt(dr);

      if(dr < rcut){                            //checks if the jth particle is within the cutoff radius
        v += 1.0/pow(dr,12) - 1.0/pow(dr,6);    //computes 1/4 * the potential energy between i-j and adds it to the previously calculated
        p += 1.0/pow(dr,12) - 0.5/pow(dr,6);    //computes the term to be averaged in the pressure formula
      }
    }          
  }

  for (int i=0; i<npart; ++i)
    kin += 0.5 * (vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i]);                   //measures kinetic energy

  walker[iv] = 4.0 * v;                                                       //current value of potential energy
  walker[ik] = kin;                                                           //current value of kinetic energy
  walker[it] = (2.0 / 3.0) * kin/(double)npart;                               //current value of temperature
  walker[ie] = 4.0 * v + kin;                                                 //current value of total energy
  walker[ip] = rho * walker[it] + 48./(3.*vol) * (p/( (double)npart ));       //current value of pressure

  return;
}


//Given a block, this method resets its averages, normalization and acceptance rate
void Reset(int iblk) {
   
  if(iblk == 1) {                       //if iblk is the first block, sets to 0 the global averages of every observable
    for(int i=0; i<n_props; ++i) {

      glob_av[i] = 0;
      glob_av2[i] = 0;
    }
  }

  for(int i=0; i<n_props; ++i) {
    blk_av[i] = 0;                    //resets the block average of every observable
  }
   
  blk_norm = 0;
  attempted = 0;
  accepted = 0;
}


//This method updates the block averages
void Accumulate(void) {

  for(int i=0; i<n_props; ++i){
    blk_av[i] = blk_av[i] + walker[i];     //adds the currently measured values to the ones already computed
  }

  blk_norm = blk_norm + 1.0;               //the block normalization is the number of times the measures were made in that block
}


//This method prints the resulting averages and error for the current block
void Averages(int iblk) {
    
  ofstream Epot, Ekin, Etot, Temp, Press;
  const int wd=12;
  
  cout << "Block number " << iblk << endl;                            //prints block number and acceptance rate
  cout << "Acceptance rate " << accepted/attempted << endl << endl;
    
  Epot.open("output_epot.dat",ios::app);            //opens the output files
  Ekin.open("output_ekin.dat",ios::app);
  Temp.open("output_temp.dat",ios::app);
  Etot.open("output_etot.dat",ios::app);
  Press.open("output_press.dat",ios::app);
    
  stima_pot = blk_av[iv]/blk_norm/(double)npart;    //average potential energy of the block per particle
  glob_av[iv] += stima_pot;                         //sum of the average potential energies per particle
  glob_av2[iv] += stima_pot*stima_pot;
  err_pot=Error(glob_av[iv],glob_av2[iv],iblk);     //computes the statistical error
    
  stima_kin = blk_av[ik]/blk_norm/(double)npart;    //average kinetic energy of the block per particle
  glob_av[ik] += stima_kin;                         //sum of the average kinetic energies per particle
  glob_av2[ik] += stima_kin*stima_kin;
  err_kin=Error(glob_av[ik],glob_av2[ik],iblk);     //computes the statistical error

  stima_etot = blk_av[ie]/blk_norm/(double)npart;   //average total energy of the block per particle
  glob_av[ie] += stima_etot;                        //sum of the average total energies per particle
  glob_av2[ie] += stima_etot*stima_etot;
  err_etot=Error(glob_av[ie],glob_av2[ie],iblk);    //computes the statistical error

  stima_temp = blk_av[it]/blk_norm;                 //average temperature of the block
  glob_av[it] += stima_temp;                        //sum of the average temperatures
  glob_av2[it] += stima_temp*stima_temp;
  err_temp=Error(glob_av[it],glob_av2[it],iblk);    //computes the statistical error

  stima_press = blk_av[ip]/blk_norm;                //average pressure of the block
  glob_av[ip] += stima_press;                       //sum of the average pressures
  glob_av2[ip] += stima_press*stima_press;
  err_press=Error(glob_av[ip],glob_av2[ip],iblk);   //computes the statistical error

  //prints the block number, last estimate of potential energy per particle, its progressive global average and statistical error
  Epot << setw(wd) << iblk <<  setw(wd) << stima_pot << setw(wd) << glob_av[iv]/(double)iblk << setw(wd) << err_pot << endl;
  //prints the block number, last estimate of kinetic energy per particle, its progressive global average and statistical error
  Ekin << setw(wd) << iblk <<  setw(wd) << stima_kin << setw(wd) << glob_av[ik]/(double)iblk << setw(wd) << err_kin << endl;
  //prints the block number, last estimate of total energy per particle, its progressive global average and statistical error
  Etot << setw(wd) << iblk <<  setw(wd) << stima_etot << setw(wd) << glob_av[ie]/(double)iblk << setw(wd) << err_etot << endl;
  //prints the block number, last estimate of temperature, its progressive global average and statistical error
  Temp << setw(wd) << iblk <<  setw(wd) << stima_temp << setw(wd) << glob_av[it]/(double)iblk << setw(wd) << err_temp << endl;
  //prints the block number, last estimate of pressure, its progressive global average and statistical error
  Press << setw(wd) << iblk <<  setw(wd) << stima_press << setw(wd) << glob_av[ip]/(double)iblk << setw(wd) << err_press << endl;


  cout << "----------------------------" << endl << endl;

  Epot.close();
  Ekin.close();
  Etot.close();
  Temp.close();
  Press.close();
}


//this method prints the final configuration
void ConfFinal(void){
  ofstream WriteConf, WriteVelocity, WriteSeed;

  cout << "Print final configuration to file config.out" << endl << endl;
  WriteConf.open("config.out");
  WriteVelocity.open("velocity.out");

  for (int i=0; i<npart; ++i){
    WriteConf << x[i]/box << "   " <<  y[i]/box << "   " << z[i]/box << endl;   //prints the final xyz cordinates of every particle normalized by the box length 
    WriteVelocity << vx[i] << "   " <<  vy[i] << "   " << vz[i] << endl;        //prints the final velocity components of every particle
  }

  WriteConf.close();
  WriteVelocity.close();

  rnd.SaveSeed();           //prints the final value of the seed so that it will be possible to restart the simulation in the exact configuration it ended
}


//writes the nconf configuration in .xyz format (the one needed to visualize the simulation in Ovito)
void ConfXYZ(int nconf){ 
  ofstream WriteXYZ;

  WriteXYZ.open("frames/config_" + to_string(nconf) + ".xyz");
  WriteXYZ << npart << endl;
  WriteXYZ << "This is only a comment!" << endl;
  for (int i=0; i<npart; ++i){
    WriteXYZ << "LJ  " << Pbc(x[i]) << "   " <<  Pbc(y[i]) << "   " << Pbc(z[i]) << endl;
  }
  WriteXYZ.close();
}

//algorithm for checking the periodic boundary conditions with side L=box
double Pbc(double r){
  return r - box * rint(r/box);   //given an imput distance r, it calculates how many box lengths r is, rounds it to the nearest intiger and subtracts it from r
}

//function that computes the standard deviation of the mean of iblk blocks
double Error(double sum, double sum2, int iblk){
  return sqrt(fabs(sum2/(double)iblk - pow(sum/(double)iblk,2))/(double)iblk);
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