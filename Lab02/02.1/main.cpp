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


EXERCISE 02.1

AIM OF THE CODE:
   This code's purpose is to compute nuumerically the value of the integral of f(x) = g(x)p(x) = pi/2 * cos(pi*x/2) between 0 and 1 using two different
   sampling approaches. The first one considers p(x)=1 and draws uniform numbers from [0,1) and substitutes them into f(x) to get the average value;
   the second one approximates f(x) with a similar normalized distribution function d(x) = 3/2 * (1 - x^2) whose cumulative is invertible, and
   calculates the integral as the average value of g(x)p(x) / d(x) where x are drawn from d(x). This process, called importance sampling, guarantees
   the convergence of the integral to be faster.

KEY WORDS:
   - Monte Carlo integration
   - Importance sampling
   - Inversion of the cumulative distribution
   - Statistical uncertainties in Monte Carlo simulations: the blocking method
*/


#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <cstdlib>
#include "random.h"

using namespace std;

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


   // *** integral by uniform sampling ***

   int M=1000000;             //total number of throws
   int N=100;                 //number of blocks
   int L=int(M/N);            //number of throws in each block

   vector<double> ave_I(N);            //vector of the <I> of each block
   vector<double> ave2_I(N);           //vector of the <I>^2 of each block

   vector<double> sum_prog(N);         //vector of the progressive sums
   vector<double> sum2_prog(N);
   vector<double> err_prog(N);

   double x =0.;

   ofstream out;
   out.open( "../1.uniform.out" );

   for (int i=0; i<N; i++){               //cycles over the blocks

      for (int j=0; j<L; j++){                                    //cycles over the throws in each block
         ave_I[i] += M_PI/2 * cos(M_PI * rnd.Rannyu()/2 );        //sum of f(x) where x is a uniformly generated random number
      }
      ave_I[i] /= L;                      //normalizes to get the average value
      ave2_I[i] = ave_I[i] * ave_I[i];
   }

   //printing the final values and errors using the blocking method
   for (int i=0; i<N; i++){

      if (i == 0){
         sum_prog[0] = ave_I[0];
         sum2_prog[0] = ave2_I[0];
      } else {
         sum_prog[i] = sum_prog[i-1]*i + ave_I[i];         //computes the sum of the first i+1 terms summing the i+1th term to i*(the mean value of the first i terms)
         sum2_prog[i] = sum2_prog[i-1]*i + ave2_I[i];
      }

      sum_prog[i] = sum_prog[i]/(i+1);
      sum2_prog[i] = sum2_prog[i]/(i+1);
      err_prog[i] = error(sum_prog[i], sum2_prog[i], i);

      out << sum_prog[i] << " " << err_prog[i] << endl;
   }

   out.close();


   // *** integral by importance sampling ***

   out.open( "../1.sampling.out" );

   for (int i=0; i<N; i++){

      ave_I[i]=0.;            //resets the average values
      ave2_I[i]=0.;

      for (int j=0; j<L; j++){
         x = rnd.Function02_1();

         ave_I[i] += M_PI/3. * ( cos(M_PI * x/2.) ) / (1.-x*x) ;     //sum of g(x)p(x)/d(x) where x is random number drawn from d(x)
      }
      ave_I[i] /= L;
      ave2_I[i] = ave_I[i] * ave_I[i];
   }



   for (int i=0; i<N; i++){
      sum_prog[i]=0.;
      sum2_prog[i]=0.;

      if (i == 0){
         sum_prog[0] = ave_I[0];
         sum2_prog[0] = ave2_I[0];
      } else {
         sum_prog[i] = sum_prog[i-1]*i + ave_I[i];         //computes the sum of the first i+1 terms summing the i+1th term to i*(the mean value of the first i terms)
         sum2_prog[i] = sum2_prog[i-1]*i + ave2_I[i];
      }

      sum_prog[i] = sum_prog[i]/(i+1);
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
