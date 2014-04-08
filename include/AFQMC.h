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
      AFQMC(const MPO<complex<double>,Quantum> &ham,const Trotter &Trotter,const MPS<complex<double>,Quantum> &Psi0,int DW, int Nwalkers,double dtau);
      
      //Destructor
      ~AFQMC();
      
      //Let the walkers propagate for steps steps
      void Walk(const int steps);

      //Propagate my population of walkers for 1 timestep. Return the sum of the coeff of my walkers.
      double Propagate();
      
      //Control the population of walkers based on scaling * weight
      void PopulationControl(double scaling);

      //Calculate the single walker projected energies, update the energy history, calculate the fluctuation metric, and the total projected energy
      complex<double> gEP();

      //Write the projected energy, target energy
      void write(const int step,const int nwalkers, const double projectedEnergy, const double targetEnergy);

      //Setup the walkers
      void SetupWalkers();

   private:
      
      //!The Trotter decomposition of the Hamiltonian
      Trotter *trotter;

      //!The MPO form of the hamiltonian
      MPO<complex<double>,Quantum> ham;

      //!number of trotter terms
      int n_trot;
      
      //!The MPS truncation dimension of the walkers
      int DW;
      
      //!The total desired number of walkers
      int Nwalkers;
      
      //!The imaginary time step size (positive)
      double dtau;

      //propagator
      Propagator *P;
      
      //!Trial wfn 
      MPS<complex<double>,Quantum> Psi0;
      
      //!The walkers
      std::vector<Walker*> theWalkers;

      //!length of the chain
      int L;

      //!physical dimension
      int d;
      
};

#endif
