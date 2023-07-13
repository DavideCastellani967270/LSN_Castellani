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


EXERCISE 10

AIM OF THE CODE:
This code's purpuse is to numerically find the solution to the two-dimensional travelling salesman problem, which consists in finding the shortest path possible
that connects N points on a plane with a closed line. Firt, it generates 34 randomly distibuted cities on a circumference (to which we know the exat solution)
or in a square; than it creates a population of possible solutions to the problem and measures their fitness, which means "how well" they solve it.
In order to get to the unique solution, the code performs a number of changes in the population (generations) following the recombination algorithms found in
genetics and sorts the best fitting paths.
  

KEY WORDS:
  - Optimization problems
  - Genetic alghorithms
  - Parallel Computing
  - Message Passing Interface
*/


#include <iostream>
#include <fstream>
#include <ostream>
#include <cmath>
#include <iomanip>
#include "Parallel_TSP.h"
#include <mpi.h>

using namespace std;

//functions
void Input( int );
void Print( const territory&, const population& );
vector<int>& newSchedule ( int );


int main(int argc, char* argv[]){

  //parallelisation
  int size, rank;
  MPI_Init(&argc, &argv);                   //initializes the MPI enviroment

  MPI_Comm_size(MPI_COMM_WORLD, &size);     //comunicates to all processes how many they are
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);     //comunicates to all processes which from [0,...,size-1] they are

  Input( rank );                      //initializes the simulation
  
  //Creation of the territory (worldwide) and the starting population (restricted to the continent)
  territory terr ( "American_capitals.dat" );
  population pop;

  //support variables and vectors
  chromosome bestchromo;
  vector<int> send_bestpath(Ncities), receive_bestpath(Ncities);
  vector<int> all_bestpaths(size * Ncities);
  int receiver;

  //parallel simulation with migrations
  for (int imigr = 0; imigr < Nmigr; imigr ++){         //cycles over the number of migrations

    for (int igen = 0; igen < Ngen_settled; igen ++){   //cycles over the settled generations
      pop.NewGen();                                     //updates the population
      pop.CalcLenght( terr );                           //computes the lenght of every chromosome in the population
      pop.CalcFitness();                                //sorts the chromosomes from the shortest to the longest
    }
    
    bestchromo = pop.GetChromo(0);              //every continent gets its best chromosome
    send_bestpath = bestchromo.GetPath();       //its best path is extracted so that it can be sent
    
    schedule = newSchedule( size );             //new schedule (worldwide) is created
    receiver = schedule[rank];                  //every continent reads who is the receiver of its message

    //every continent sends its path send_bestpath, containg Ncities intigers, to the receiver and receives receive_bestpath, containg Ncities intigers, from an unspecified source
    MPI_Sendrecv(send_bestpath.data(), Ncities, MPI_INT, receiver, 0,
                 receive_bestpath.data(), Ncities, MPI_INT, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    bestchromo.SetPath(receive_bestpath);       //every continent sets the received path
    pop.SetChromo( 0, bestchromo );             //every continent sets the new imported chromosome as the best one
    
  }

  MPI_Barrier(MPI_COMM_WORLD);           //waits until every continent has sent&received the path from last migration

  //gathering of the best paths
  pop.CalcLenght( terr );                //computes the lenght of every chromosome in the population
  pop.CalcFitness();                     //sorts from best to worst

  bestchromo = pop.GetChromo(0);         //gets the best chromosome
  send_bestpath = bestchromo.GetPath();  //gets the best path so that it can be sent

  //every continent sends its path send_bestpath, containg Ncities intigers, to the continent ranked 0 that collects them in all_bestpaths
  MPI_Gather(send_bestpath.data(), Ncities, MPI_INT,
             all_bestpaths.data(), Ncities, MPI_INT,
             0, MPI_COMM_WORLD);

  MPI_Barrier(MPI_COMM_WORLD);          //waits until every continent has sent its path

  
  if (rank == 0){
    pop_size = size;
    population bestpop;                           //continent 0 creates a population formed by the best individuals

    for(int i = 0; i < size; i++){
      bestchromo = bestpop.GetChromo(i);

      for (int j = 0; j < Ncities; j++){
        bestchromo.SetOrder(j, all_bestpaths[i * Ncities + j]);     //copies the best paths
      }

      bestpop.SetChromo(i, bestchromo);
      bestpop.GetChromo(i).Check();
    }


    bestpop.CalcLenght( terr );                //computes the lenght of every chromosome in the population
    bestpop.CalcFitness();                     //sorts from best to worst

    Print ( terr, bestpop );                   //prints the cities coordinates in the order they appear in the best path
  }

  MPI_Finalize();

  return 0;
}


//initializes the random number generator and reads the input parameters
void Input( int rank ){

  //reads two primes to set up the random numbers generator (RNG)
  ifstream Primes("Primes");
  if (Primes.is_open()){
    Primes >> p1 >> p2 ;
    for (int i =0; i < rank; i++)
      Primes >> pskip >> pskip;
    Primes >> p3 >> p4;

  } else cerr << "PROBLEM: Unable to open Primes" << endl;
  Primes.close();

  ifstream input("seed.in");
  string property;
  if (input.is_open()){
    while ( !input.eof() ){
      input >> property;
      if( property == "RANDOMSEED" ){
        input >> world_seed[0] >> world_seed[1] >> world_seed[2] >> world_seed[3];         //sets up the RNG
        world_rnd.SetRandom(world_seed,p1,p2);

        for (int j =0; j < 4; j++){
          seed[j] = world_seed[j] + int(world_rnd.Rannyu(0, 100)) *rank;        //randommy changes the seed in function of the rank
        }
          
        rnd.SetRandom(seed,p3,p4);
      }
    }
    input.close();
  } else cerr << "PROBLEM: Unable to open seed.in" << endl;
  

  ifstream ReadInput;
  ReadInput.open("input.in");

  ReadInput >> pop_size;        //reads the number of chromosomes that makes a population
  ReadInput >> Ngen;            //reads the number of generations the program simulates
  ReadInput >> Nmigr;           //reads the number of migrations

  Ngen_settled = (int) Ngen/Nmigr;  //reads the number of generations before migrating

  ReadInput >> pselect;         //reads the exponent in the chromosome selection algorithm
  ReadInput >> prob_cross;      //reads the probability of doing a crossover
  ReadInput >> prob_mute;       //reads the probability of doing a mutation

  ReadInput.close();
}

//this method creates a sending schedule
vector<int>& newSchedule ( int Ncontinents ){
  schedule.clear();                       //deletes previous schedule

  for (int i = 0; i < Ncontinents; i++){     
    schedule.push_back(i);                 //creates a ordered schedule [0,1,2,3,...,size-1] -> every continent migrates to itself
  }

  int pos, i_value;
  for (int i=1; i < Ncontinents; i++){                  //scrumbles the schedule
    
    i_value = schedule[i];                              //saves the old continent label in the ith position

    do{ pos = int( world_rnd.Rannyu(0, Ncontinents) );     //extracts a random continent for the ith continent to be swapped with
    } while ( i == schedule[pos] );                         //checks that every message will be sent to a new continent

    schedule[i] = schedule[pos];          //swaps the two continents
    schedule[pos] = i_value;
  }

  return schedule;
}


//this function prints the cities coordinates in the order they appear in the best path
void Print( const territory& terr, const population& pop ){

  ofstream PrintOutput; 
  PrintOutput.open("../AC.output.out");

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