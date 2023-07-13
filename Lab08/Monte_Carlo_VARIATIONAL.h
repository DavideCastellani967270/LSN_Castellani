/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/


#include <iostream>
#include <fstream>
#include <ostream>
#include <cmath>
#include <iomanip>
#include <vector>

using namespace std;


//Random numbers
#include "random.h"
int seed[4];
Random rnd;

//configuration
double x;
const double xmax=3.;
const int nbins=100;
double bin_size = 2*xmax/nbins;

//parameters, observables
double Psi2[nbins];   
double  p;
vector<double> Old_eigen(2), New_eigen(2), Eigen(2);

// thermodynamical state
double sigma, mu, temp, temp_min, beta;
double old_sigma, old_mu, new_sigma, new_mu;

// simulation
int nstep, nblk, restart;
int last = 0;
double delta, delta_par, cool_rate;
int SAstep = 0;
const int wd=12; 


//functions
void Input(void);
void ChangeParam(void);


#ifndef __Hamiltonian__
#define __Hamiltonian__

class Hamiltonian{

    public:

    //constructor
    Hamiltonian(double mu, double sigma, double delta, int nblk, int nstep ){

        m_mu = mu;
        m_sigma = sigma;
        m_delta = delta;
        m_nblk = nblk;
        m_nstep = nstep;
    }



    //methods
    vector<double> CalcEigenvalue(void){
        vector<double> EneandErr(2);

        for(int iblk=1; iblk <= nblk; iblk++){  //Simulation
            Reset(iblk);                         //Reset block averages
        
            for(int istep=1; istep <= nstep; istep++){
                Move();
                Measure();
                Accumulate();                      //Update block averages
            }

            Averages(iblk);       //Print results for current block
            EneandErr[0] = eigen;
            EneandErr[1] = err_ham;
        }

        return EneandErr;
    }
     
    //Given a block, this method resets its averages, normalization and acceptance rate
    void Reset(int iblk) {
   
        if(iblk == 1){         //if iblk is the first block, sets to 0 the global averages of every observable
            glob_av = 0;
            glob_av2 = 0;
        }

        blk_av = 0;
        blk_norm = 0;

        attempted = 0;
        accepted = 0;
    }

    //this method evolves the present configuration into a new one (not necessarily different in MC)
    void Move(void){
        double p, psi2_old, psi2_new;
        double xnew;

        //Old
        psi2_old = Prob_density(x);             //computes the value of the probability distribution we're sampling in the old x

        //New
        xnew = x + m_delta*(rnd.Rannyu() - 0.5);  //proposes a randomly selected new position
        psi2_new = Prob_density(xnew);          //computes the value of the probability distribution we're sampling in the new x

        //Metropolis test
        p = psi2_new/psi2_old;                  //the acceptance probability is the ratio between the two probability densities

        if(p >= rnd.Rannyu())  {
            x = xnew;                             //updates the old position with the proposed one
            accepted ++;            //keeps track of the acceptance rate of the algorithm
        } 
  
        attempted ++;
        return;
    }

    //this method evaluates in x the square modulus of the wave function = probability density function
    double Prob_density(double x){
        return pow( exp( - pow( x - m_mu, 2) / (2.*m_sigma*m_sigma) ) + exp( - pow( (x + m_mu), 2) / (2.*m_sigma*m_sigma) ), 2 );     //the normalization is not needed thanks to the fact that we'll use the Metropolis algorithm 
    }
    
    //this method measures the properties in a given configuration, here H psi / psi
    void Measure(void) {
        double Epot = pow(x, 4) - 2.5 * pow(x, 2);
        double Ekin = exp( -(x-m_mu)*(x-m_mu)/(2.0* m_sigma*m_sigma) ) * ( (x-m_mu)*(x-m_mu)/(m_sigma * m_sigma) - 1 )/(m_sigma*m_sigma)
                  + exp( -(x+m_mu)*(x+m_mu)/(2.0* m_sigma*m_sigma) ) * ( (x+m_mu)*(x+m_mu)/(m_sigma * m_sigma) - 1 )/(m_sigma*m_sigma);
        Ekin = - 0.5 * Ekin /( exp( - (x-m_mu)*(x-m_mu)/(2.0 * m_sigma*m_sigma) ) + exp( - (x+m_mu)*(x+m_mu)/(2.0 * m_sigma*m_sigma) ) );    //hbar, m = 1

        walker = (Epot + Ekin);
        return;
    }

    //This method updates the block averages
    void Accumulate(void) {
        blk_av = blk_av + walker;         //adds the currently measured values to the ones already computed
        blk_norm = blk_norm + 1.0;        //the block normalization is the number of times the measures were made in that block
    }

    //This method prints the resulting averages and error for the current block
    void Averages(int iblk) {

        stima_ham = blk_av/blk_norm;                    //average total energy of the block
        glob_av += stima_ham;                           //sum of the average total energies
        glob_av2 += stima_ham*stima_ham;

        if (iblk == nblk ){                             //for all temperatures, only the final block results are printed
            cout << "Temperature " << temp << endl;
            cout << "Acceptance rate " << (double) accepted/attempted << endl << endl;
            cout << "----------------------------" << endl << endl;
        }

        eigen = glob_av/(double)iblk;
        err_ham = Error(glob_av,glob_av2,iblk);           //computes the statistical error   
    
        if(last){
            ofstream BestConf;
            BestConf.open("FinalHam.out", ios::app);
            BestConf << iblk << "  " << stima_ham << "  " << eigen << "  " << err_ham << endl;
            BestConf.close();
        }
    }

    //error
    double Error(double sum, double sum2, int iblk){
        return sqrt(fabs(sum2/(double)iblk - pow(sum/(double)iblk,2))/(double)iblk);
    }



    //histogram of psi^2
    void Histo_Psi2(void){
        int norm=0;

        for(int istep=1; istep <= nstep; istep++){
            Move();
            if(x > - xmax && x < xmax){
                Psi2[ int((x + xmax) / bin_size) ] += 1;   //filling the histo
                norm ++;
            } 
        }
        norm = norm/nbins;

        //Print 
        ofstream Psi2_hist;
        Psi2_hist.open("Psi2_Histogram.out");
        for(int i=0; i<=nbins; i++){
            Psi2_hist << -xmax + bin_size * i << " \t" << Psi2[i]/(2*xmax * norm) << endl;
        }
        Psi2_hist.close();

        return;
    }


    //metod for changing the parameters
    void SetParam(double newmu, double newsigma){
        m_mu = newmu;
        m_sigma = newsigma;
    }


    //this method prints the final configuration
    void ConfFinal(void){
        ofstream WriteConf, WriteSeed;

        cout << "Print final position to file config.out" << endl << endl;

        WriteConf.open("config.out");
        WriteConf << x << endl;
        WriteConf.close();

        rnd.SaveSeed();                    //prints the final value of the seed so that it will be possible to restart the simulation in the exact configuration it ended
    }

    //writes the nconf configuration in .xyz format (the one needed to visualize the simulation in Ovito)
    void ConfXYZ(int nconf){ 
        ofstream WriteXYZ;

        WriteXYZ.open("frames/config_" + to_string(nconf) + ".xyz");
        WriteXYZ << "This is only a comment!" << endl;
  
        WriteXYZ << "LJ  " << x << "   " << 0 << "   " << 0 << endl;

        WriteXYZ.close();
    }


    private: 

    // averages
    double blk_av, blk_norm, accepted, attempted;
    double stima_ham, glob_av, glob_av2;
    double eigen;
    double err_ham;

    // simulation
    int m_nstep, m_nblk;
    double m_delta;

    //parameters, observables
    double m_mu, m_sigma;
    double walker;
};




#endif

/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
