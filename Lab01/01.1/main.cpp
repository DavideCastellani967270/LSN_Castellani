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


EXERCISE 01.1

AIM OF THE CODE:
    This code's purpose is to evaluate the integral of r between [0,1] as the mean value of random numbers extracted uniformly in that interval.
    It also computes its variance as the mean of the squared differences between the extracted numbers and the expected value 0.5.
    While doing so, the statistical error (standard deviation of the mean) for both observables has been determined using the blocking method.
    In part 3, I test the random number generator with a chi squared test. If the drawing distribution is truly random, one should expect its
    value to be approximately the same of the number of sub-divisions of the interval (in this case, 100).

KEY WORDS:
    - Central limit theorem
    - Mean value, variance and standard deviation of the mean
    - Statistical uncertainties in Monte Carlo simulations: the blocking method
    - Pseudo and true randomness
*/



#include <iostream>
#include <fstream>
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


    int M=100000;              //Total number of throws
    int N=100;                 //Number of blocks
    int L=int(M/N);            //Number of throws in each block, please use for M a multiple of N

    vector<double> ave;      //Vector of the averages <r> in each block
    vector<double> ave2;     //Vector of the <r>^2 in each block
    
    vector<double> ave_sigma;       //Vector of the <variance> of each block
    vector<double> ave2_sigma;      //Vector of the <variance>^2 of each block

    vector<double> sum_prog(N);         //Vector of the progressive sums
    vector<double> sum2_prog(N);
    vector<double> err_prog(N);

    double r, sum, sum_var;


    for (int i=0; i<N; i++){            //cycles through the blocks
        sum=0.;
        sum_var=0.;

        for (int j=0; j<L; j++){        //cycles through the number of throws in each block

            r = rnd.Rannyu();
            sum += r;
            sum_var += (r-0.5)*(r-0.5);       
        }

        ave.push_back(sum /L);              //saves the estimate value of <r> of the current block
        ave2.push_back(ave[i]*ave[i]);      //saves the estimate value of <r>^2 of the current block

        ave_sigma.push_back(sum_var /L);                    //saves the estimate value of <(r-0.5)^2> of the current block
        ave2_sigma.push_back(ave_sigma[i]*ave_sigma[i]);    //saves the estimate value of <(r-0.5)^2>^2 of the current block
    }


    ofstream out;
	out.open("../1.sum_r_prog.out");

    //printing the final values and errors of the integral using the blocking method
    for (int i=0; i<N; i++){

        if (i == 0){
            sum_prog[0] = ave[0];
            sum2_prog[0] = ave2[0];
        } else{
            sum_prog[i] = sum_prog[i-1]*i + ave[i];         //computes the sum of the first i+1 terms summing the i+1th term to i*(the mean value of the first i terms)
            sum2_prog[i] = sum2_prog[i-1]*i + ave2[i];
        }

        sum_prog[i] = sum_prog[i]/(i+1);            //normalizes the sum to get the mean value
        sum2_prog[i] = sum2_prog[i]/(i+1);
        err_prog[i] = error(sum_prog[i], sum2_prog[i], i);

        out << sum_prog[i] << " " << err_prog[i] << endl;
    }

    out.close();

    //printing the final values and errors of the variance using the blocking method
    out.open("../1.sum_sigma_prog.out");

    for (int i=0; i<N; i++){

        if (i == 0){
            sum_prog[0] = ave_sigma[0];
            sum2_prog[0] = ave2_sigma[0];
        } else{
            sum_prog[i] = sum_prog[i-1]*i + ave_sigma[i];
            sum2_prog[i] = sum2_prog[i-1]*i + ave2_sigma[i];
        }

        sum_prog[i] = sum_prog[i]/(i+1);
        sum2_prog[i] = sum2_prog[i]/(i+1);
        err_prog[i] = error(sum_prog[i], sum2_prog[i], i);

        out << sum_prog[i] << " " << err_prog[i] << endl;
    }

    out.close();


    //********* PART 3 *********

    M=100;                          //number of bins
    N=10000;                        //number of throws

    double chi2=0.;

    vector<int> counter (M, 0);     //vector that keeps track of the counts n_i in each bin

    out.open("../1.chi2j.out");

    for (int j=0; j<100; j++){

        chi2=0.;                                //resets the variable and the vector
        counter.assign(counter.size(), 0);

        for (int i=0; i<N; i++){                //extracts the numbers and updates the counter
            r = int( 100*rnd.Rannyu() );
            counter[r] += 1;
        }

        for (int i=0; i<M; i++){
            chi2 += pow( counter[i] - double(N/M), 2 ) / double(N/M);       //computes the chi2 values
        }

        out << chi2 << endl;
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