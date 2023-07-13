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


EXERCISE 02.2

AIM OF THE CODE:
   This code simulates 10^4 random walks of 100 steps each in two different types of lattice: the first is cubic, which means that every step is taken
   in one of the xyz direction (positive or negative), while the second one makes steps in a random direction in space, selecting an azimut and a polar angle.
   These RW are divided in 100 blocks and for each one it has been computed the sqrt( <|r_j|^2> ) > (the sqrt of the average distance^2 from the origin at
   step j); then their total average value (considering all blocks) and its statistical uncertainty in function of the step number are saved in a file so
   that they can be plotted.

KEY WORDS:
   - Random walks
   - Statistical uncertainties in Monte Carlo simulations: the blocking method
*/


#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <vector>
#include <cstdlib>
#include "random.h"

using namespace std;

//vector from the origin in 3D
struct xyz_vector {
   double x, y, z;
};


//function that computes the standard deviation of the mean of n+1 blocks
double error ( double ave , double ave2, int n){
    if (n==0)
        return 0;
    else
        return sqrt((ave2 - ave*ave)/n);
}

//function that computes the distance from the origin squared
double distance2 ( xyz_vector r ){
   return r.x*r.x + r.y*r.y + r.z*r.z;
}

//function that makes a step with fixed lenght in a randomly selected direction (given the lattice type)
xyz_vector new_step ( Random& rnd, xyz_vector r, double length, string dir_type ){

   double dir[3] = {0., 0., 0.};                   //vector containing lenghts of the step in each xyz direction

   if ( dir_type == "dis" ){                       //checks the lattice type: "dis" = discrete, everything else will be interpreted to be continuum

      int dis_dir = int ( rnd.Rannyu(0, 3) );      //randommly selects the direction of the step: x = [0,1), y = [1,2) and z = [2,3)
      int dis_vers = rnd.CoinToss();               //randommly selects if it will be a forward (+1) or backward (-1) step

      dir[dis_dir] = length * dis_vers;
   
   } else{
      
      double phi = rnd.Rannyu( 0, 2* M_PI );       //random azimuth
      double theta = rnd.Polar_angle();            //random polar angle

      dir[0] = length * sin(theta) * cos(phi);
      dir[1] = length * sin(theta) * sin(phi);
      dir[2] = length * cos(theta);
   }

   r.x += dir[0];
   r.y += dir[1];
   r.z += dir[2];

   return r;
}
 


int main (int argc, char *argv[]){

   //calling and initializing the (pseudo-)random number generator
   Random rnd;
   int seed[4];
   int p1, p2;
   ifstream Primes("Primes");
   if (Primes.is_open()){
      Primes >> p1 >> p2 ;
   } else cerr << "PROBLEM: Unable to open Primes" << endl;
   Primes.close();

   ifstream input("seed.in");
   string property;
   if (input.is_open()){
      while ( !input.eof() ){
         input >> property;
         if( property == "RANDOMSEED" ){
            input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
            rnd.SetRandom(seed,p1,p2);
         }
      }
      input.close();
   } else cerr << "PROBLEM: Unable to open seed.in" << endl;



   int M_rw = 10000;       //total number of Random Walks
   int step = 100;         //number of steps for each RW
   int b = 100;            //total number of blocks in which the RW are divided
   int L = int (M_rw/b);   //number of RW for each block


   //*** Discete lattice ***
   vector <double> dis_ave_r(step);             //vector which in position i is stored the sum of |r_i|^2 of every RW in a fixed block

   vector <double> dis_ave_blocks(step);        //vector which in position i is stored the sums of sqrt( <|r_i|^2> ) of every block in a fixed step i
   vector <double> dis_ave2_blocks(step);       //vector which in position i is stored the sums of <|r_i|^2> of every block in a fixed step i
   vector <double> dis_err_blocks(step);

   xyz_vector dis_r = {0., 0., 0.};             //discrete case position vector

   ofstream dis_out;
   dis_out.open("../2.discrete_rw.out");


   //*** Continuous lattice ***
   vector <double> con_ave_r(step);

   vector <double> con_ave_blocks(step);
   vector <double> con_ave2_blocks(step);
   vector <double> con_err_blocks(step);

   xyz_vector con_r = {0., 0., 0.};             //continuum case position vector

   ofstream con_out;
   con_out.open("../2.continuous_rw.out");


   for (int k = 0; k < b; k++){                 //cycles over the blocks
      for (int i = 0; i < L; i++){              //fixed the block, cycles over the RW in it

         dis_r = {0., 0., 0.};                  //resets to the origin
         con_r = {0., 0., 0.};

         for (int j = 0; j < step; j++){        //fixed the RW, cycles over the steps in it

            dis_r = new_step( rnd, dis_r, 1., "dis" );      //updates the position with a new step
            con_r = new_step( rnd, con_r, 1., "con" );

            dis_ave_r[j] += distance2( dis_r );             //computes |r_j|^2 (distance^2 from the origin at step j of Rw i in block k)
            con_ave_r[j] += distance2( con_r );             //and adds it to the equivalents of the others RWs (0,1,...,i-1) in the same block k
         }
      }
      
      for (int j = 0; j < step; j++){

         dis_ave_blocks[j] += sqrt( dis_ave_r[j]/L );    //evaluates sqrt( <|r_j|^2> ) (average distance from the origin at step j in block k)
         con_ave_blocks[j] += sqrt( con_ave_r[j]/L );    //and adds it to the equivalents of the others blocks (0,1,...,k-1)

         dis_ave2_blocks[j] += dis_ave_r[j]/L;     //evaluates <|r_j|^2> (average distance^2 from the origin at step j in block k)
         con_ave2_blocks[j] += con_ave_r[j]/L;     //and adds it to the equivalents of the others blocks (0,1,...,k-1)
         
         dis_ave_r[j] = 0.;                  //resets <|r_j|^2> before changing block
         con_ave_r[j] = 0.;
      }
   }

   for (int j=0; j<step; j++){

      dis_ave_blocks[j] /= b;                //evaluates < sqrt( <|r_j|^2> ) > (average of the average distance from the origin at step j)
      dis_ave2_blocks[j] /= b;               //evaluates < <|r_j|^2> > (average of the average distance^2 from the origin at step j)
               
      dis_err_blocks[j] = error( dis_ave_blocks[j], dis_ave2_blocks[j], b-1 );

      dis_out << dis_ave_blocks[j] << "\t" << dis_err_blocks[j] << endl;


      con_ave_blocks[j] /= b;
      con_ave2_blocks[j] /= b;
               
      con_err_blocks[j] = error( con_ave_blocks[j], con_ave2_blocks[j], b-1 );

      con_out << con_ave_blocks[j] << "\t" << con_err_blocks[j] << endl;
   }

   dis_out.close();
   con_out.close();

   rnd.SaveSeed();
   return 0;
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