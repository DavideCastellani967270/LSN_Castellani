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


EXERCISE 01.3

AIM OF THE CODE:
   This code's purpose is to calculate the value of pi simulating the Buffon's experiment. Infact, it is possible to get an estimate of it from random 
   throws of a needle (with lenght L) onto a plane with equally spaced lines (d > L is the distance): the probability that the needle intersects one of
   the lines is P = 2L / (pi * d), from which pi = 2L / (P * d).
   Let the lines be the y=n*d family, with n intiger: their translational simmetry (continuous along the x axis and descrete in y) allows us to simplify the 
   system in the one dimensional finite segment of the y axis [-d/2, d/2). Here the lines are replaced with only the y=0 point and the needle with a random
   y_center in [-d/2, d/2) and the uppermost (or lowermost) y value from the center that the needle touches (y_extr).
   Having (y_center >= 0) and (y_center - y_extr <= 0) or (y_center < 0) and (y_center + y_extr >= 0) immediately means that the needle intersects the y=0 
   line somewhere along the x axis. A blocking average has been used to evaluate pi and its statistical error.

KEY WORDS:
   - Inversion of the cumulative distribution
   - Statistical uncertainties in Monte Carlo simulations: the blocking method
*/


#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cstdlib>
#include <cmath>
#include "random.h"

using namespace std;

struct needle {
   double y_center;              //y value of the needle center
   double y_extr;                //uppermost (lowermost) y value from the center that the needle touches 
};

//function that computes the standard deviation of the mean of n+1 blocks
double error ( double ave , double ave2, int n){
    if (n==0)
        return 0;
    else
        return sqrt((ave2 - ave*ave)/n);
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

   //Defining parameters
   double d = 5.;             //spacing between the lines on the plane
   double L = 1.;              //lenght of the needle
   int Nhit = 0;              //number of intersections

   int Nblocks = 100;              //total number of blocks
   int Npi = 100;                  //number of estimates of pi in each bloch
   int Nthrows = 10000;             //number of throws in each estimate of pi

   vector<double> ave;       //Vector of the <Pi> of each block
   vector<double> ave2;      //Vector of the <Pi>^2 of each block

   vector<double> sum_prog(Nblocks);         //Vector of the progressive sums
   vector<double> sum2_prog(Nblocks);
   vector<double> err_prog(Nblocks);

   double x, y, r2, Pi, piSum=0, piSum2 = 0;

   //Generating throws
   needle ndl;

   for (int j=0; j < Nblocks; j++){             //cycles over the blocks
      piSum =0;
      piSum2 =0;

      for (int i =0; i < Npi; i++){               //cycles over the estimates of pi in each bloch
         Nhit =0;
            
         for (int k=0; k < Nthrows; k++){           //cycles over the throws in each estimate of pi
                        
            ndl.y_center = rnd.Rannyu( -0.5*d, 0.5*d );        //the fundamental block of the plane is the spacing d with 0.5d on each side

            if ( abs(ndl.y_center) <= L/2. ){         //speeds up running time by not extracting the direction of the needle if it's impossible for it to intersect the line
               do{                     
                  x = rnd.Rannyu();                   //draws a random point in the first quarter of the unit circle (origin aside) and computes 
                  y = rnd.Rannyu();                   //sin(theta) = y/sqrt(x^2 + y^2)

                  r2 = x*x + y*y;
               } while (r2 > 1 || r2 == 0);

               ndl.y_extr = L/2. * y/sqrt(r2);           //0 = horizontal, 0.5*L = vertical
   
               if (ndl.y_center <= 0 && ndl.y_center + ndl.y_extr >= 0)    //the center of the needle is under the line but the top is above it
                  Nhit += 1;
               if (ndl.y_center > 0 && ndl.y_center - ndl.y_extr <= 0)     //the center of the needle is above the line but the bottom is under it
                  Nhit += 1;
            }
         }

         if (Nhit > 0){
            Pi = 2*L*Nthrows/(Nhit*d);         //evaluates pi
            
            piSum += Pi;
            piSum2 += Pi*Pi;

         } else
            cout << "No hit has been registered in block number " << j+1 << endl;
      }
      
      ave.push_back( piSum/Npi );
      ave2.push_back( ave[j]*ave[j] );
   }

   ofstream out;
   out.open( "../3.Pi.out" );

   //printing the final values and errors using the blocking method
   for (int i=0; i<Nblocks; i++){

      if (i == 0){
         sum_prog[0] = ave[0];
         sum2_prog[0] = ave2[0];
      } else {
         sum_prog[i] = sum_prog[i-1]*i + ave[i];         //computes the sum of the first i+1 terms summing the i+1th term to i*(the mean value of the first i terms)
         sum2_prog[i] = sum2_prog[i-1]*i + ave2[i];
      }

      sum_prog[i] = sum_prog[i]/(i+1);            //normalizes the sum to get the mean value
      sum2_prog[i] = sum2_prog[i]/(i+1);
      err_prog[i] = error(sum_prog[i], sum2_prog[i], i);

      out << sum_prog[i] << " " << err_prog[i] << endl;
   }

   out.close();

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
