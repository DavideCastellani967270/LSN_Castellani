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


EXERCISE 06

AIM OF THE CODE:
This code's purpuse is to simulate a 1D Ising model where the only interactions considered are the ones with the adjacent spins (with periodic boundary
conditions). More precisely, given in input the maximum temperature and the model parameters J (the exchange interaction) and h (the external magnetic field), 
it performs 10000 Monte Carlo steps per block using the M(RT)^2 or the Gibb's sampling algorithm and repeats the process for 20 smaller temperature values.
While doing so, multiple observable quantities are measured (internal energy (h=0), heat capacity (h=0), magnetic susceptibility (h=0)) and magnetization (h!=0)
and computes their progressive average values and uncertainties as a functions of the increasing number of blocks and temperature.

KEY WORDS:
  - Statistical mechanics
  - Ising model
  - Metropolis algorithm
  - Gibb's sampling method
  - Periodic boundary conditions
  - Statistical uncertainties in Monte Carlo simulations: the blocking method
*/


#include <iostream>
#include <fstream>
#include <ostream>
#include <cmath>
#include <iomanip>
#include "Monte_Carlo_ISING_1D.h"

using namespace std;

int main(){ 

  Input();                                      //inizialization

  do{
    cout << "Temperature: " << temp << endl << endl;

    for(int iblk=1; iblk <= nblk; ++iblk) {       //cycles over the nblk blocks
      Reset(iblk);                                //reset block averages

      for(int istep=1; istep <= nstep; ++istep){  //cycles over the steps in each block
        Move(metro);                              //updates the positions
        Measure();                                //measures the observables
        Accumulate();                             //updates block averages
      }

      Averages(iblk);     //prints results for current block
    }
    ConfFinal();          //writes the final configuration

    temp -= dt;           //updates the temperature
    beta = 1.0/temp;      //updates the beta parameter

  } while (temp >= min_temp);

  return 0;
}


//this method initialises the system with the starting configurations and then prints the initial energy value
void Input(void){
  ifstream ReadInput;

  //prints some informations about the code
  cout << "Classic 1D Ising model             " << endl;
  cout << "Monte Carlo simulation             " << endl << endl;
  cout << "Nearest neighbour interaction      " << endl << endl;
  cout << "Boltzmann weight exp(- beta * H ), beta = 1/T " << endl << endl;
  cout << "The program uses k_B=1 and mu_B=1 units " << endl;

  //reads two primes to set up the random numbers generator (RNG)
  int p1, p2;
  ifstream Primes("Primes");
  Primes >> p1 >> p2 ;
  Primes.close();

  ifstream input("seed.in");                            //opens the seed.in file in order to set the RNG for a new run
  input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
  rnd.SetRandom(seed,p1,p2);                            //sets up the RNG
  input.close();
  
  //read input informations
  ReadInput.open("input.dat");

  ReadInput >> temp;                                    //reads the INITIAL temperature (in reduced units), prints its value and calculates the beta factor
  beta = 1.0/temp;
  dt = (temp - min_temp)/ntemp;                         //computes the updating temperature step lenght and prints its value
  cout << "Temperature = " << temp << endl;
  cout << "Temperature step dT = " << dt << endl;

  ReadInput >> nspin;                                   //reads the number of particles and prints it
  cout << "Number of spins = " << nspin << endl;

  ReadInput >> J;                                       //reads the the exchange interaction modulus and prints it
  cout << "Exchange interaction = " << J << endl;

  ReadInput >> h;                                       //reads the the external magnetic field modulus and prints it
  cout << "External field = " << h << endl << endl;
    
  ReadInput >> metro;                                   // 1 = Metropolis, anything else = Gibbs

  ReadInput >> nblk;                                    //reads and prints the number of blocks

  ReadInput >> nstep;                                   //reads and prints the number of steps in each block

  if(metro==1) cout << "The program perform Metropolis moves" << endl;
  else cout << "The program perform Gibbs moves" << endl;
  cout << "Number of blocks = " << nblk << endl;
  cout << "Number of steps in one block = " << nstep << endl << endl;
  ReadInput.close();


  //prepare arrays for measurements
  iu = 0;       //array position referring to internal energy
  ic = 1;       //array position referring to heat capacity
  im = 2;       //array position referring to magnetization
  ix = 3;       //array position referring to magnetic susceptibility
  n_props = 4;  //Number of observables

  //initial configuration: all spin are put up
  for (int i=0; i<nspin; ++i){
    s[i] = 1;
  }
  
  //evaluates the properties of the initial configuration
  Measure();

  //prints initial values of the internal energy per spin
  cout << "Initial energy per spin = " << walker[iu]/(double)nspin << endl;
}


//this method evolves the present configuration into a new one (not necessarily different)
void Move(int metro){
  int o;
  double p, energy_old, energy_new, r;

  for(int i=0; i<nspin; ++i){
  
    o = (int)(rnd.Rannyu()*nspin);              //randomly selects a particle to be flipped

    if(metro==1) //Metropolis
    {
      energy_old = Boltzmann( s[o], o );              //stores old energy
      energy_new = Boltzmann( -s[o], o );             //evaluates new energy after the spin flip
      p = exp( -beta * (energy_new - energy_old)  );  //in the Metropolis algorithm, the acceptance probability is given by the Boltzmann weight evaluated with the energy difference

      if ( p > 1 )               //if the new configuration has less energy, the proposed move is always accepted
        p = 1;
      
      if (p > rnd.Rannyu()){     //if the extracted number is smaller than p, the proposed move is accepted else it is rejected    
        s[o] = -s[o];            //flips the spin
        accepted ++;
      }
      attempted ++;
    }
    else //Gibbs sampling
    {
      p = 1. / (1 + exp( - 2 * beta * (J  *( s[Pbc(o-1)] + s[Pbc(o+1)] ) + h) ) );   //calculates the acceptance probability of imposing the spin +1

      r = rnd.Rannyu();         //generates random number
      if (p > r)                //assign the spin
        s[o] = +1;            
      else
        s[o] = -1;              //if the spin +1 move is not accepted, Gibbs sampling automatically flips it to be -1

      accepted ++;
      attempted ++;
    }
  }
}

//this method computes the contribuition of the spin ip with value sm to the energy in the Ising model: 
//since the energy difference in the exponent of Boltzmann weight cancels out the common terms, this is all is needed to calculate it
double Boltzmann(int sm, int ip){
  double ene = -J * sm * ( s[Pbc(ip-1)] + s[Pbc(ip+1)] ) - h * sm;        //the only interactions considered are the ones with the adjacent spins
  return ene;
}

//this method measures the observables in a given configuration
void Measure(){
  double u = 0.0, m = 0.0;

  for (int i=0; i<nspin; ++i){                                            //cycles over spins
    u += -J * s[i] * s[Pbc(i+1)] - 0.5 * h * (s[i] + s[Pbc(i+1)]);        //computes the internal energy with nearest neighbors interaction
    m += s[i];                                                            //magnetization is proportional to the average spin value
  }

  walker[iu] = u;             //current value of total internal energy
  walker[ic] = u*u;           //current value of (total internal energy)^2   ->   to compute heat capacity
  walker[im] = m;             //current value of the sum of all spins        ->   to compute magnetization
  walker[ix] = m*m;           //current value of the (sum of all spins)^2    ->   to compute magnetic susceptibility
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


//this method prints the resulting averages and error for the current block
void Averages(int iblk){
    
  ofstream Ene, Heat, Mag, Chi;
  const int wd=12;
    
  cout << "Block number " << iblk << endl;                            //prints block number and acceptance rate
  cout << "Acceptance rate " << accepted/attempted << endl << endl;
    
  if(metro==1){                        //opens the output files checking the type of sampling method and the presence of an external magnetic field
    if(h == 0){
      Ene.open("M.output.ene.0",ios::app);
      Heat.open("M.output.heat.0",ios::app);
      Chi.open("M.output.chi.0",ios::app);
    } else{
      Mag.open("M.output.mag.0",ios::app);
    }
  } else{
    if(h == 0){
    Ene.open("G.output.ene.0",ios::app);
    Heat.open("G.output.heat.0",ios::app);
    Chi.open("G.output.chi.0",ios::app);
    } else{
      Mag.open("G.output.mag.0",ios::app);
    }
  }
    
  if(h==0){
    //Energy
    stima_u = blk_av[iu]/blk_norm/(double)nspin;      //average internal energy of the block per spin
    glob_av[iu]  += stima_u;                          //sum of the average internal energies per spin
    glob_av2[iu] += stima_u*stima_u;
    err_u=Error(glob_av[iu],glob_av2[iu],iblk);       //computes the statistical error
    //prints the block number, the temperature, the last estimate of internal energy per spin, its progressive global average and statistical error
    Ene << setw(wd) << iblk <<  setw(wd) << temp << setw(wd) << stima_u << setw(wd) << glob_av[iu]/(double)iblk << setw(wd) << err_u << endl;
    Ene.close();

    //Heat capacity
    stima_c = pow(beta, 2) * (blk_av[ic]/blk_norm - pow(blk_av[iu]/blk_norm, 2)) /(double)nspin;      //average heat capacity of the block per spin
    glob_av[ic]  += stima_c;                        //sum of the average heat capacities per spin
    glob_av2[ic] += stima_c*stima_c;
    err_c=Error(glob_av[ic],glob_av2[ic],iblk);     //computes the statistical error
    //prints the block number, the temperature, the last estimate of heat capacity per spin, its progressive global average and statistical error
    Heat << setw(wd) << iblk <<  setw(wd) << temp << setw(wd) << stima_c << setw(wd) << glob_av[ic]/(double)iblk << setw(wd) << err_c << endl;
    Heat.close();

    //Magnetic susceptibility
    stima_x = beta * blk_av[ix]/blk_norm/(double)nspin;      //average magnetic susceptibility of the block per spin
    glob_av[ix]  += stima_x;                                 //sum of the average magnetic susceptibilities per spin
    glob_av2[ix] += stima_x*stima_x;
    err_x=Error(glob_av[ix],glob_av2[ix],iblk);              //computes the statistical error
    //prints the block number, the temperature, the last estimate of magnetic susceptibility per spin, its progressive global average and statistical error
    Chi << setw(wd) << iblk <<  setw(wd) << temp << setw(wd) << stima_x << setw(wd) << glob_av[ix]/(double)iblk << setw(wd) << err_x << endl;
    Chi.close();

  } else {
    //Magnetization
    stima_m = blk_av[im]/blk_norm/(double)nspin;      //average magnetization of the block per spin
    glob_av[im]  += stima_m;                          //sum of the average magnetizations per spin
    glob_av2[im] += stima_m*stima_m;
    err_m=Error(glob_av[im],glob_av2[im],iblk);       //computes the statistical error
    //prints the block number, the temperature, the last estimate of magnetization per spin, its progressive global average and statistical error
    Mag << setw(wd) << iblk <<  setw(wd) << temp << setw(wd) << stima_m << setw(wd) << glob_av[im]/(double)iblk << setw(wd) << err_m << endl;
    Mag.close();
  }

  cout << "----------------------------" << endl << endl;
}


//this method prints the final configuration
void ConfFinal(void){
  ofstream WriteConf;

  cout << "Print final configuration to file config.final " << endl << endl;
  WriteConf.open("config.final");

  for (int i=0; i<nspin; ++i){
    WriteConf << s[i] << endl;        //prints the final list of spins
  }

  WriteConf.close();

  rnd.SaveSeed();             //prints the final value of the seed so that it will be possible to restart the simulation in the exact configuration it ended
}


//algorithm for checking the periodic boundary conditions
int Pbc(int i){
  if(i >= nspin) i = i - nspin;     //if the spin index i = nspin + k with 0 <= k <= nspin, Pbc makes i = k
  else if(i < 0) i = i + nspin;     //if the spin index i = -k with 0 <= k <= nspin, Pbc makes i = nspin + k

  return i;
}

//function that computes the standard deviation of the mean of iblk blocks
double Error(double sum, double sum2, int iblk){
  if(iblk==1) return 0.0;
  else return sqrt((sum2/(double)iblk - pow(sum/(double)iblk,2))/(double)(iblk-1));
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