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


EXERCISE 03.1

AIM OF THE CODE:
   Following Black and Scholes theory of option pricing, this code computes the European Call-option and Put-option prices via 10^6 Monte Carlo simulations.
   The final asset prices are calculated both by sampling it directly (from time 0 -> T) and discreetly (from time 0 -> Dt -> 2Dt -> ... -> T).
   It has been used data blocking to compute the options averages and their statistical error.

KEY WORDS:
   - Brownian motion
   - Stochastic differential equations
   - Geometric Brownian motion
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


//function that computes the standard deviation of the mean of n+1 blocks
double error (double ave, double ave2, int n){
    if (n==0)
        return 0;
    else
        return sqrt((ave2 - ave*ave)/n);
}

//function that returns the maximum between two imput numbers
double max (double a, double b){
   if (a>=b)
      return a;
   else
      return b;
}

//function that directly returns the asset price at an instant "time" given its initial value S0 and the GBM parameters mu and sigma
double final_asset_price ( Random& rnd, double S0, double mu, double sigma, double time ){
   double W = rnd.Gauss(0, time);
   return S0 * exp( (mu - (sigma*sigma /2))*time + sigma*W );
}

//function that returns the asset price after an increment in time Dt given its previous value Sti and the GBM parameters mu and sigma
double discrete_asset_price ( Random& rnd, double Sti, double mu, double sigma, double Dt ){
   double Z = rnd.Gauss(0, 1);
   return Sti * exp( (mu - (sigma*sigma /2))*Dt + sigma*Z*sqrt(Dt) );
}

//function that returns the call option price given the risk-free interest rate r, the delivery time T and the asset price at that time SiT
double call_option_price ( double r, double T, double SiT , double K ){
   return exp(-r*T)*max(0., SiT- K);
}

//function that returns the put option price given the risk-free interest rate r, the delivery time T and the asset price at that time SiT
double put_option_price ( double r, double T, double SiT , double K ){
   return exp(-r*T)*max(0., K-SiT);
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


   const double S0 = 100.;                //asset price at t=0
   const double T = 1.;                   //delivery time
   const double K = 100.;                 //strike price
   const double r = 0.1;                  //risk-free interest rate
   const double vol = 0.25;               //volatility

   const int B = 100;                     //number of blocks
   const int N_asset_prices = 1E6;        //number of estimated asset prices
   const int L = int (N_asset_prices/B);  //number of estimated asset prices per block

   double ST;

   vector<double> Call_ave(B);               //vectors containing the averages of each block
   vector<double> Put_ave(B);


   //*** DIRECT SAMPLING OF THE FINAL ASSET PRICE ***

   for (int j=0; j < B; j++){                    //cycles over the blocks
      for (int i = 0; i < L; i++){                    //cycles over the number of estimates of asset prices per block

         ST = final_asset_price( rnd, S0, r, vol, T );

         Call_ave[j] += call_option_price( r, T, ST, K );         //evaluates the ith call option price and adds it to the previous estimates of the jth block
         Put_ave[j] += put_option_price( r, T, ST, K );
      }

      Call_ave[j] /= L;                //computes the average call option price for each block
      Put_ave[j] /= L;
   }

   //support variables for blocking average computation
   double sum_C=0., ave_C=0., sum_C2=0., ave_C2=0.;
   double sum_P=0., ave_P=0., sum_P2=0., ave_P2=0.;

   ofstream out_C;
   out_C.open("../1.unic_call.out");

   ofstream out_P;
   out_P.open("../1.unic_put.out");

   for (int j=0; j<B; j++){
      sum_C += Call_ave[j];               //adds the jth block average call option price to the previous (0,1,...,j-1)
      sum_P += Put_ave[j];

      sum_C2 += Call_ave[j]*Call_ave[j];              //adds the (jth block average call option price)^2 to the previous (0,1,...,j-1)
      sum_P2 += Put_ave[j]*Put_ave[j];
      
      ave_C = sum_C / (j+1);              //finds the progressive average of the first j+1 block average call option prices (0,1,...,j-1,j)
      ave_P = sum_P / (j+1);

      ave_C2 = sum_C2 / (j+1);            //finds the progressive average of the first j+1 (block average call option prices)^2 (0,1,...,j-1,j)
      ave_P2 = sum_P2 / (j+1);

      out_C << ave_C << " " <<  error(ave_C, ave_C2, j) << endl;        //prints the progressive average and its statistical error
      out_P << ave_P << " " <<  error(ave_P, ave_P2, j) << endl;
   }

   out_C.close();
   out_P.close();


   //*** DISCRETE SAMPLING OF THE ASSET PRICE PATH ***

   const int N_intervals = 100;           //number of intervals into which delivery time is divided
   double Dt = T / N_intervals;           //lenght of those intervals

   for (int j=0; j < B; j++){          //cycles over the blocks

      Call_ave[j] = 0;                 //resets the vectors used in the previous part of the code
      Put_ave[j] = 0;

      for (int i = 0; i < L; i++){           //cycles over the number of estimates of asset prices per block

         ST = S0;                      //resets the starting asset price (t=0)
         
         for (int k = 0; k < N_intervals; k++){
            ST = discrete_asset_price( rnd, ST, r, vol, Dt );        //updates ST from T -> T+Dt 
         }

         Call_ave[j] += call_option_price( r, T, ST, K );         //evaluates the ith call option price and adds it to the previous estimates of the jth block
         Put_ave[j] += put_option_price( r, T, ST, K );
      }

      Call_ave[j] /= L;                  //computes the average call option price for each block
      Put_ave[j] /= L;
   }

   sum_C =0, sum_C2 =0;                  //resets the support variables
   sum_P =0, sum_P2 =0;

   out_C.open("../1.disc_call.out");
   out_P.open("../1.disc_put.out");

   for (int j=0; j<B; j++){
      sum_C += Call_ave[j];                     //adds the jth block average call option price to the previous (0,1,...,j-1)
      sum_P += Put_ave[j];

      sum_C2 += Call_ave[j]*Call_ave[j];        //adds the (jth block average call option price)^2 to the previous (0,1,...,j-1)
      sum_P2 += Put_ave[j]*Put_ave[j];
      
      ave_C = sum_C / (j+1);                    //finds the progressive average of the first j+1 block average call option prices (0,1,...,j-1,j)
      ave_P = sum_P / (j+1);

      ave_C2 = sum_C2 / (j+1);                  //finds the progressive average of the first j+1 (block average call option prices)^2 (0,1,...,j-1,j)
      ave_P2 = sum_P2 / (j+1);

      out_C << ave_C << " " <<  error(ave_C, ave_C2, j) << endl;           //prints the progressive average and its statistical error
      out_P << ave_P << " " <<  error(ave_P, ave_P2, j) << endl;
   }

   out_C.close();
   out_P.close();

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