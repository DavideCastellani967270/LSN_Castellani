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


EXERCISE 01.2

AIM OF THE CODE:
   This code's purpose is to simulate three different dice by extracting intigers between 1 and 6.
   The first one is a normal die (extracts with uniform probability); the second one is an exponential die and the third one a lorentzian die.
   Then, for each die, it has been computed and printed 10000 mean values of N = {1, 2, 10, 100} sigle throws in order to check the central limit theorem
   (at least for the standard and exponential dice).

KEY WORDS:
   - Inversion of the cumulative distribution
   - Central limit theorem
*/


#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cstdlib>
#include "random.h"

using namespace std;
 
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


   // STANDARD DIE
   ofstream out;
   out.open( "../2.standard_die.out" );

   vector<int> N {1, 2, 10, 100};               //vector containing the number of single die throws to be averaged
   double SN = 0.;

   for (int k=0; k<4; k++){               //cycles through the vector N
   
      for (int i=0; i<10000; i++){

         SN = 0.;

         for (int j=0; j<N[k]; j++)
            SN += int( rnd.Rannyu(1, 7) );         //sums N[k] unifomly extracted intigers between 1 and 6
      
         SN = SN/N[k];
         out << SN << endl;
      }
   }

   out.close();


   // EXPONENTIAL DIE
   out.open( "../2.exponential_die.out" );
   int result = 0;

   for (int k=0; k<4; k++){
   
      for (int i=0; i<10000; i++){

         SN = 0.;

         for (int j=0; j<N[k]; j++){
            do{
               result = int( rnd.Exp(1)+1 );    //the exponential distribution draws numbers from [0, +infty]: for this reason they have to be shifted (+1)
            } while ( result >= 7 );            //and casted. This isn't a huge loss in efficiency thanks to the fast decay of the exponential function.
               
            SN += result;
         }

         SN = SN/N[k];
         out << SN << endl;
      }
   }

   out.close();


   // LORENTZIAN DIE
   out.open( "../2.lorentzian_die.out" );

   for (int k=0; k<4; k++){
   
      for (int i=0; i<10000; i++){

         SN = 0.;

         for (int j=0; j<N[k]; j++){
            do{                                                //the lorentzian distribution draws numbers from [-infity, +infty]: so I take their absolute
               result = int( abs(rnd.Lorentz(0, 1)) +1 );      //value, shift it by +1 and delete the number if it is >=7. Just like the exponential, the 
            } while ( result >= 7 );                           //Cauchy-Lorentz is fast decaying, so the efficiency isn't much affected.
               
            SN += result;
         }

         SN = SN/N[k];
         out << SN << endl;
      }
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
