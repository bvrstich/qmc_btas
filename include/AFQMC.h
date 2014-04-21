#ifndef AFQMC_H
#define AFQMC_H

#include <iostream>
#include <iomanip>
#include <complex>
#include <vector>

using std::complex;
using std::vector;

#include "Propagator.h"

class AFQMC {

   public:
   
      //constructor with input trialwavefunction
      AFQMC(const MPS< complex<double> > &,double,int);
      
      //Destructor
      virtual ~AFQMC();
      
      //Let the walkers propagate for steps steps
      void walk(int);

      //Propagate my population of walkers for 1 timestep. Return the sum of the coeff of my walkers.
      double Propagate();
      
      //Control the population of walkers based on scaling * weight
      void PopulationControl(double);

      //Calculate the single walker projected energies, update the energy history, calculate the fluctuation metric, and the total projected energy
      complex<double> gEP();

      //Write the projected energy, target energy
      void write(int,double,double);

      //Setup the walkers
      void SetupWalkers();

   private:
      
      //!The Trotter decomposition of the Hamiltonian
      Trotter *trotter;

      //!number of trotter terms
      int n_trot;
      
      //!The total desired number of walkers
      int Nw;
      
      //!The imaginary time step size (positive)
      double dtau;

      //propagator
      Propagator *P;
      
      //!Trial wfn 
      MPS< complex<double> > Psi0;
      
      //!The walkers
      std::vector<Walker*> walker;
      
};

#endif
