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


EXERCISE 09

AIM OF THE CODE:
This code's purpuse is to numerically find the solution to the two-dimensional travelling salesman problem, which consists in finding the shortest path possible
that connects N points on a plane with a closed line. Firt, it generates 34 randomly distibuted cities on a circumference (to which we know the exat solution)
or in a square; than it creates a population of possible solutions to the problem and measures their fitness, which means "how well" they solve it.
In order to get to the unique solution, the code performs a number of changes in the population (generations) following the recombination algorithms found in
genetics and sorts the best fitting paths.
  

KEY WORDS:
  - Optimization problems
  - Genetic alghorithms
*/


#include <iostream>
#include <fstream>
#include <ostream>
#include <cmath>
#include <iomanip>
#include "GA_Traveling_Salesman.h"

using namespace std;

//functions
void Input( void );
void Print( const territory&, const population& );


int main(){

  //Creation of starting population and territory
  Input();
  territory terr (terr_type);
  population pop;

  //for every generation  
  for (int i = 0; i < Ngen; i++){
    pop.NewGen();                   //updates the less fit half of the population
    pop.CalcLenght( terr );         //computes the lenght of every chromosome in the population
    pop.CalcFitness();              //sorts the chromosomes from the shortest to the longest
  }

  Print ( terr, pop );              //prints the cities coordinates in the order they appear in the best path

  return 0;
}


//initializes the random number generator and reads the input parameters
void Input( void ){
  //reads two primes to set up the random numbers generator (RNG)
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
            input >> seed[0] >> seed[1] >> seed[2] >> seed[3];         //sets up the RNG
            rnd.SetRandom(seed,p1,p2);
         }
      }
      input.close();
   } else cerr << "PROBLEM: Unable to open seed.in" << endl;

  ifstream ReadInput;
  ReadInput.open("input.in");

  ReadInput >> terr_type;       //reads the type territory: 0 -> cities randomly arranged on a circumference,  !=0 -> cities randomly arranged in a square
  ReadInput >> Ncities;         //reads the number of cities
  ReadInput >> pop_size;        //reads the number of chromosomes that makes a population
  ReadInput >> Ngen;            //reads the number of generations the program simulates
  ReadInput >> pselect;         //reads the exponent in the chromosome selection algorithm
  ReadInput >> prob_cross;      //reads the probability of doing a crossover
  ReadInput >> prob_mute;       //reads the probability of doing a mutation

  cout << "This program solves the travelling salesman problem for" << endl;
  cout << Ncities << " cities randomly distributed ";
  if (terr_type)
    cout << "on a circumference." << endl << endl;
  else
    cout << "in a square." << endl << endl;
  cout << "Population size: " << pop_size << endl;
  cout << "Generations: " << Ngen << endl << endl;

  ReadInput.close();
}

//this function prints the cities coordinates in the order they appear in the best path
void Print( const territory& terr, const population& pop ){

  ofstream PrintOutput;
  if (terr_type == 0)                           //checks the type of the territory so that two separate files are created
    PrintOutput.open("C.output.out");
  else
    PrintOutput.open("S.output.out");

  cout << "Best path:" << endl;  

  chromosome bestchromo = pop.GetChromo(0);     //copies the best path
  city c;

  for (int i = 0; i < Ncities; i++){            //prints the cities in the order that says the best chromosome
    cout << bestchromo.GetOrder(i) << ", " ;    //prints on screen the best path by writing the order in which the cities appear

    c = terr.GetCity( bestchromo.GetOrder(i) ); //gets the coordinates of the city in the position reported in the position i of the chromosome
    PrintOutput << c.x << "\t" << c.y << endl;
  }

  PrintOutput.close();
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