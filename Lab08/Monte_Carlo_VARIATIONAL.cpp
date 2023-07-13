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


EXERCISE 08

AIM OF THE CODE:
This code's purpuse is to variationally optimize the ground state of a single quantum particle in a one dimensional (1D) space confined by the following 
external potential: V(x)= x^4 - 5/2 x^2
Due to the fact that this potential is not analytically solvable, in order to obtain an approximate wave function for the ground state, the implemented
variational Monte Carlo method uses a trial wave function, parametrized by a set of two variational parameters: sigma (the standard deviation) and mu (the
mean value). Therefore, once the kinetic term has been analitically calculated, the expectation value for the Hamiltonian is sampled using a Metropolis 
algorithm and than the optimal parameters are extracted
  

KEY WORDS:
  - Optimization problems
  - Variational Monte Carlo
  - Simulated annealing
  - Statistical uncertainties in Monte Carlo simulations: the blocking method
*/

#include <iostream>
#include <fstream>
#include <ostream>
#include <cmath>
#include <iomanip>
#include "Monte_Carlo_VARIATIONAL.h"

using namespace std;

int main(){ 
  Input();                               //Inizialization

  Hamiltonian H(mu, sigma, delta, nblk, nstep);

  ofstream Ham, Param;

  Ham.open("ham.out",ios::app);
  Param.open("param.out",ios::app);

  while(temp > temp_min){
    SAstep ++;
  
    old_mu = mu;
    old_sigma = sigma;
    Old_eigen = H.CalcEigenvalue();

    ChangeParam();                //Change parameters
    H.SetParam(new_mu, new_sigma);
    New_eigen =  H.CalcEigenvalue();

    p = exp( -beta * ( New_eigen[0] - Old_eigen[0] ) );

    if (p >= rnd.Rannyu()){
      mu = new_mu;
      sigma = new_sigma;
      Eigen = New_eigen;

    } else{
      mu = old_mu;
      sigma = old_sigma;
      Eigen = Old_eigen;
    }

    H.SetParam(mu, sigma);

    //prints the temperature, the current SA step, the energy final global average and statistical error
    Ham <<  setw(wd) << temp <<  setw(wd) << SAstep  <<  setw(wd) << Eigen[0] <<  setw(wd) << Eigen[1] << endl;
    //prints the temperature, the current SA step and the parameters
    Param <<  setw(wd) << temp <<  setw(wd) << SAstep <<  setw(wd) << mu <<  setw(wd) << sigma << endl;

    temp *= cool_rate;
    beta = 1./temp;
  }

  Ham.close();
  Param.close();

  last++;
  Eigen= H.CalcEigenvalue();

  H.Histo_Psi2();

  H.ConfFinal();                           //Write final configuration

  return 0;
}


//this method initialises the system with the starting configurations and velocities and then prints the initial observables values
void Input(void){
  ifstream ReadInput, ReadConf, Primes, Seed;

  //prints some informations about the code
  cout << "Single quantum particle in 1D potential        " << endl;
  cout << "Variational MC(NVT) simulation       " << endl << endl;
  cout << "Potential v(r) = x^4 - 5/2 * x^2" << endl << endl;
  cout << "Weight: module squared of the wavefunction " << endl << endl;

  //reads two primes to set up the random numbers generator (RNG)
  int p1, p2;
  Primes.open("Primes");
  Primes >> p1 >> p2 ;
  Primes.close();

  //read input informations
  ReadInput.open("input.in");

  ReadInput >> restart;     // 0 = no restart, 1 = restart

  if(restart)
    Seed.open("seed.out");                                //opens the seed.out file in order to set the RNG to the exact configuration of the last run
  else 
    Seed.open("seed.in");                                 //opens the seed.in file in order to set the RNG for a new run
  Seed >> seed[0] >> seed[1] >> seed[2] >> seed[3];
  rnd.SetRandom(seed,p1,p2);                              //sets up the RNG
  Seed.close();

  ReadInput >> mu;                                        //reads and prints the centers of the starting gaussians                                    
  cout << "Starting gaussian centers = +-" << mu << endl;

  ReadInput >> sigma;                                     //reads and prints the standard deviations of the starting gaussians
  cout << "Starting gaussian standard deviation = " << sigma << endl;
    
  ReadInput >> delta;                                     //reads the moving time parameter (in reduced units, dt^* = dt sqrt( epsilon/ m*sigma^2 ))
  ReadInput >> delta_par;                                 //reads the moving gaussian centers and standard deviation parameter 

  ReadInput >> cool_rate;                                 //reads and prints the cooling rate
  cout << "Cooling rate: T' = " << cool_rate << " T" << endl;

  ReadInput >> temp;                                      //reads and prints the starting temperature
  cout << "Starting temperature (LJ units) = " << temp << endl;
  beta = 1./temp;

  ReadInput >> temp_min;                                  //reads and prints the temperature that stops the simulated annealing algorithm
  cout << "Stopping temperature (LJ units) = " << temp_min << endl;

  ReadInput >> nblk;                                      //reads and prints the number of blocks

  ReadInput >> nstep;                                     //reads and prints the number of steps in each block

  cout << "The program perform Metropolis moves with uniform translations" << endl;
  cout << "Moves parameter = " << delta << endl;
  cout << "Number of blocks = " << nblk << endl;
  cout << "Number of steps in one block = " << nstep << endl << endl;
  ReadInput.close();

  //Read initial configuration
  cout << "Read initial configuration" << endl << endl;
  if(restart){
    ReadConf.open("config.out");
  } else {
    ReadConf.open("config.in");
  }

  ReadConf >> x;        //reads the configuration
  ReadConf.close();


  return;
}

//this method moves with a Metropolis algorithm the parameters
void ChangeParam(){

  new_mu = mu + delta_par*(rnd.Rannyu() - 0.5);
  new_sigma = sigma + delta_par*(rnd.Rannyu() - 0.5);

  return;
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