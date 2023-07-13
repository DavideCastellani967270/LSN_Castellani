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
*/


#include <iostream>
#include <fstream>
#include <ostream>
#include <cmath>
#include <iomanip>
#include <vector>

using namespace std;

//Random numbers
#include "random.h"
int world_seed[4];
int seed[4];
int p1, p2, p3, p4, pskip;
Random rnd, world_rnd;

//parameters, observables and support variables
int Ncities, terr_type, pop_size, Ngen, Nmigr, Ngen_settled;               //terr_type = 0 -> cities randomly on a circumference, terr_type !=0 -> cities randomly in a square
double pselect, prob_cross, prob_mute;
double angle, sum_order, sum2_order, prog_length;
vector<int> Msequence;
vector<int> Fsequence;
vector<int> schedule;

//cities are defined by their cartesian coordinates
struct city {
  double x;
  double y;
};



//classes
#ifndef __salesman_classes__
#define __salesman_classes__


//this class contains a vector of cities randomly generated: the order in which they appear defines their labels
class territory {

  public:

  //constructor
  territory( int terr_type ){

    if(terr_type == 0){                       //circumference
      for (int i = 0; i < Ncities; i++){

        angle = world_rnd.Rannyu(0, 2*M_PI);                          //generates a random angle in radiants
        m_Territory.push_back( {cos(angle), sin(angle)} );      //computes the coordinates of the city
      }
    } else {                                    //square
      for (int i = 0; i < Ncities; i++){
          
        m_Territory.push_back( {world_rnd.Rannyu(), world_rnd.Rannyu()} );  //generates the xy coordinates by drawing two uniform numbers between [0,1)
      }
    }
  }



  //methods
  //this method returns the city labeled with the intiger pos (after checking pos is within the territory vector)
  city GetCity( int pos ) const {
    if (pos < Ncities )
      return m_Territory[pos];
    else{
      cerr << "Error: Nonexistent city called." << endl;
      exit(-3);
    }
  }


  protected:
  vector<city> m_Territory;
};



//this class contains a permutation of the city labels (exept the first), therefore designing a path for the travelling salesman to follow. For every chromosome,
//the class also saves its lenght in order to calculate its fitness.
class chromosome {
  public:

    //constructor
    chromosome() {

      for (int i = 0; i < Ncities; i++){     
        m_Order.push_back(i);                 //creates a ordered chromosome [0,1,2,3,...,Ncities-1]
      }

      int pos, i_value;
      for (int i=1; i < Ncities; i++){        //scrumbles the chromosome (but the first digit)
        pos = int( rnd.Rannyu(1, Ncities) );  //extracts a random position for the ith city to be swapped with
        i_value = m_Order[i];                 //saves the old city label in the ith position

        m_Order[i] = m_Order[pos];            //swaps the two cities
        m_Order[pos] = i_value;
      }

      m_Length = 0.;                          //initializes the path lenght to 0
    }


    //methods
    //check function is mandatory to verify the chromosome is still a good one
    void Check() {
      if (m_Order[0] != 0) {                                  //checks that the starting city has't changed
        cerr << "Starting city changed! EXITING" << endl;
        exit(-1);
      }

      sum_order = 0;                                          //checks that the sum of the labels hasn't changed
      sum2_order = 0;                                         //checks that the sum of the labels^2 hasn't changed

      for (int i = 1; i < Ncities; i++) {
        sum_order += m_Order[i];
        sum2_order += m_Order[i] * m_Order[i];
      }

      if (sum_order != 0.5 * Ncities * (Ncities - 1)) {
        cerr << "Cities order numbers don't add up! EXITING" << endl;
        cerr << "Expected sum: " << 0.5 * Ncities * (Ncities - 1) << endl;
        cerr << "Calculated sum: " << sum_order << endl;
        exit(-2);
      }
      /*
      if (sum2_order != (1. / 6.) * Ncities * (Ncities - 1) * (2 * Ncities - 1)) {
        cerr << "The sum of squared cities order numbers is incorrect! EXITING" << endl;
        cerr << "Expected sum of squares: " << (1. / 6.) * Ncities * (Ncities - 1) * (2 * Ncities - 1) << endl;
        cerr << "Calculated sum of squares: " << sum2_order << endl;
        exit(-2);
      }
      */
    }

    //stores the chromosome length
    void SetLength(double L){
      m_Length = L;
    }

    //returns the chromosome length
    double GetLength() const{
      return m_Length;
    }

    //returns the city vector
    vector<int> GetPath() const{
      return m_Order;
    }

    //this method sets sets the path of a chromosome
    void SetPath(vector<int> newpath) {
      m_Order = newpath;
      Check();
    }

    //returns the label of the city in position pos of the chromosome
    int GetOrder(int pos) const{
      return m_Order[pos];
    }

    //changes the label of the city in position pos of the chromosome with new_ord
    void SetOrder(int pos, int new_ord){
      m_Order[pos] = new_ord;
    }

    //overriding of the "=" operator so that it can copy chromosomes
    chromosome& operator=(const chromosome& other) {
      if (this != &other) {
        m_Length = other.GetLength();         //copies in this chromosome the other's chormosome lenght
      }

      for(int i = 1; i<Ncities; i++){
        m_Order[i] = other.GetOrder(i);       //copies in this chromosome the other's chormosome city labels order
      }
      return *this;
    }

    //MUTATIONS
    //this method swaps two cities within the same chromosome
    void PairShuffle() {
      int i = int( rnd.Rannyu(1, Ncities) );    //randomly selects the first city to swap (excluding the starting one)
      int j = int( rnd.Rannyu(1, Ncities) );    //randomly selects the second city to swap (excluding the starting one)

      int i_value = m_Order[i];                 //saves the old city label in the ith position

      m_Order[i] = m_Order[j];                  //swaps the cities
      m_Order[j] = i_value;
    }

    //this method shifts the first m contiguous cities (excluding the starting one) of n positions within the same chromosome
    void Shift() {
      vector<int> OldChromo;                    //support cromosome

      for (int i=0; i<Ncities; i++){
        OldChromo.push_back(m_Order[i]);        //saves the old chromosome
      }

      int m = int (rnd.Rannyu(1, Ncities-1));   //extracts the number of cities to be shifted, from 1 (only the city in position 1) to Ncities-2 (every city but the first and last)
      int n = int (rnd.Rannyu(1, Ncities-m));   //extracts the shift length, from 1 (only the next city goes before) to Ncities-m (all the next cities go before)

      for (int i=1; i<=m+n; i++){               //only the cities in positions [1,m+n] are affected from a shift
        if (i <= m)                             
          m_Order[i+n] = OldChromo[i];          //the city labels in positions [1,m] are copied in the forwards shifted positions [1+n,m+n]
        else
          m_Order[i-m] = OldChromo[i];          //the city labels in positions [m+1, Ncities-1] are copied in the backwards shifted positions [1, Ncities-1 -m]
      }
    }

    //this method swaps two m-sized blocks of cities within the same chromosome
    void BlockShuffle() {
      int m = int (rnd.Rannyu(2, Ncities/2));   //extracts the size of the blocks to be swapped, from 2 to Ncities/2 - 1
      int j = int (rnd.Rannyu(m+1, Ncities-m+1));//extracts the position of the second block, from m+1 (the second block starts right after the first) to Ncities-m (the second block is the last m digits)

      int i_value;
      for (int i = 1; i <= m; i++){
        i_value = m_Order[i];                   //saves the old chromosome's ith value of the first block

        m_Order[i] = m_Order[j];                //swaps the ith value with the jth
        m_Order[j] = i_value;
        j++;                                    //updates the j
      }      
    }

    //this method flips a block that starts
    void Flip() {   
      int i = int (rnd.Rannyu(1, Ncities-1));  //extracts the starting position of the block, from 1 to Ncities-2 (to make the block lenght at least =2)
      int f = int (rnd.Rannyu(i+1, Ncities));    //extracts the ending position of the block, from the next after i to the last possible Ncities -1

      reverse(m_Order.begin() + i, m_Order.begin() + f);
    }


  protected:
    vector<int> m_Order;
    double m_Length;
};



//this class is a vector of chromosomes: it represents all the genetic makeup of a generation
class population : public chromosome {
  public:

  //constructors
  population(){

    for(int i = 0; i < pop_size; i++){
      chromosome chromo;
      m_Pop.push_back( chromo );
    }
  }


  //methods
  //this method returns the chromosome in the position pos of the population
  chromosome GetChromo(int pos) const{
    return m_Pop[pos];
  }

  void SetChromo(int pos, chromosome& newchromo ){
    m_Pop[pos] = newchromo;
  }


  //for each chromosome of the population, this method, given a territory, reads the travelling order stored in it and computes the total lenght of the resulting path
  void CalcLenght( const territory& terr ){

    for(int j=0; j < pop_size; j++){          //cycles over the chromosomes in the population
      prog_length = 0.;                   

      for (int i = 0; i < Ncities-1; i++){    //cycles over the cities in the chromosome
        prog_length += Distance( terr.GetCity( m_Pop[j].GetOrder(i) ), terr.GetCity( m_Pop[j].GetOrder(i+1) ) );        //distance between the city i and i+1 of the chromosome j
      }

      prog_length += Distance( terr.GetCity( m_Pop[j].GetOrder(Ncities-1) ), terr.GetCity( m_Pop[j].GetOrder(0) ) );    //distance between the last and the first city of the chromosome j
      m_Pop[j].SetLength(prog_length);        //sets the path leght of chromosome j
    }
  }

  //this method computes the distance between two cities (two points in the cartesian plane)
  double Distance( city A, city B ){
    return sqrt( pow((A.x - B.x), 2) + pow((A.y - B.y), 2) );
  }


  //this is the control operation used in the sorting algorithm: in this way the chromosomes will be arranged in increasing order of path lenght
  static bool CompareLength(const chromosome& chrom1, const chromosome& chrom2) {
    return ( chrom1.GetLength() < chrom2.GetLength() );
  }

  //this method reorganises the population sorting the chromosomes by increasing path lenght
  void CalcFitness(){
    sort( m_Pop.begin(), m_Pop.end(), CompareLength );
  }



  //selection operator: this method defines the probability law with which a father and a mother will be chosen to do a crossover
  int Select(){
    return int( pop_size * pow(rnd.Rannyu(), pselect) );    //the bigger pselect the more probable is to select the first chromosomes of the population
  }


  //this method completes a crossover operation: given two poorly fit chromosomes, selects a father, a mother and draws a random cut in the chromosome. Then,
  //it copies into the bad_chromosomes the parents rearranged accordingly to the order in wich the cut off cities appear in the parter.
  void Crossover(int bad_chromo1, int bad_chromo2) {
    int mother = Select();                                  //selects a mother
    int father = Select();                                  //selects a father

    chromosome MoChro = m_Pop[mother];                      //copies the mother so that she will not be changed
    chromosome FaChro = m_Pop[father];                      //copies the father so that he will not be changed

    int cut = int(rnd.Rannyu(1, Ncities - 1));              //cuts go from 1 (all the cities except the starting one) to Ncities -2 (only the last two)

    Msequence.clear();                                      //clears the Msequence vector before use
    Fsequence.clear();                                      //clears the Fsequence vector before use

    //stores in the Msequence/Fsequence vector the positions in the father/mother's chromosome of the cities that have been cut off from the mother/father
    for (int i = cut; i < Ncities; i++) {                   //cycles over the cut off cities

      for (int j = 1; j < Ncities; j++) {                   //cycles over the parter chromosome
        if (Msequence.size() < i-cut +1){                   //this contol makes the search stop once the ith city of the mother has been found
          if ( MoChro.GetOrder(i) == FaChro.GetOrder(j) ) { //checks if the city in position j of the father is the same city as the one in position i of the mother
            Msequence.push_back(j);                         //saves the position j
          }
        }
            
        if (Fsequence.size() < i-cut +1){                   //this contol makes the search stop once the ith city of the father has been found
          if ( FaChro.GetOrder(i) == MoChro.GetOrder(j) ) { //checks if the city in position j of the mother is the same city as the one in position i of the father
            Fsequence.push_back(j);                         //saves the position j
          }
        }
      }
    }

    if (Msequence.size() != Fsequence.size()){              //checks if the sizes of the vectors are the same
      cerr << "Error, different sizes for crossover" << endl;
      exit(-10);
    }
    if (Msequence.size() != Ncities - cut){                 //checks if the size of the vectors is the correct one
      cerr << "error, incorrect sizes for crossover" << endl;
      exit(-11);
    }

    //sorts the positions vectors so as to get the order with which the cut off cities appear in the parter's chromosome
    sort(Msequence.begin(), Msequence.end());   
    sort(Fsequence.begin(), Fsequence.end());

    //substitutes the positions in Msequence/Fsequence vectors with the cities in the father/mother's chromosome indexed by them
    for (int i = 0; i<Msequence.size(); i++){
      Msequence[i] = FaChro.GetOrder(Msequence[i]);
      Fsequence[i] = MoChro.GetOrder(Fsequence[i]);
    }

    int k = 0;
    //copies into bad_chromo1/bad_chromo2 the mother/father rearranged accordingly to the order in wich the cut off cities appear in the father/mother
    for (int i = 1; i < Ncities; i++) {

      if (i < cut) {                                        //before the cut, the parents city labels are copied
        m_Pop[bad_chromo1].SetOrder(i, MoChro.GetOrder(i));
        m_Pop[bad_chromo2].SetOrder(i, FaChro.GetOrder(i));
      } else {                                              //from the cut on, the reordered city labels are copied
        m_Pop[bad_chromo1].SetOrder(i, (Msequence[k]));
        m_Pop[bad_chromo2].SetOrder(i, (Fsequence[k]));
        k++;
      }
    }
  }

  //this method takes mutant chromosome and the random value prob in [0,prob_mute) and decides which mutations to apply
  void Mutate( int mutant, double prob ){
    if (prob < prob_mute/2){
        if (prob < prob_mute/4){              //if prob in [0, prob_mute* 1/4) -> does a pair shuffle
            m_Pop[mutant].PairShuffle();
        } else {                              //if prob in [prob_mute* 1/4, prob_mute* 1/2) -> does a shift
            m_Pop[mutant].Shift();
        }
    } else{
        if (prob < 3*prob_mute/4){            //if prob in [prob_mute* 1/2, prob_mute* 3/4) -> does a block shuffle
            m_Pop[mutant].BlockShuffle();
        } else {                              //if prob in [prob_mute* 3/4, prob_mute) -> does a flip
            m_Pop[mutant].Flip();
        }
    }
  }



  //this method updates the generation by proposing a crossover and/or a mutation to the 3/4 of the current population
  void NewGen(){
    double r_mute1, r_mute2, r_cross;

    for (int i = int(pop_size/10); i<pop_size-1; i=i+2){     //cycles over pairs of the least fitting chromosomes
      r_cross = rnd.Rannyu();

      if (r_cross < prob_cross){
        Crossover(i, i+1);          //crossover accepted
      }

      r_mute1 = rnd.Rannyu();
      r_mute2 = rnd.Rannyu();

      if (r_mute1 < prob_mute){
        Mutate( i, r_mute1 );       //mutation accepted for chromosome i
      }
      if (r_mute2 < prob_mute){
        Mutate( i+1, r_mute2 );     //mutation accepted for chromosome i+1
      }

      m_Pop[i].Check();             //checks if the chromosome are still correctly defined
      m_Pop[i+1].Check();
    }
  }


  protected:
    vector<chromosome> m_Pop;
};


#endif // __salesman_classes__




/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
