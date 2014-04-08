#ifndef WALKER_H
#define WALKER_H

#include <iostream>
#include <iomanip>
#include <complex>
#include <vector>

using std::complex;
using std::vector;

#include "SpinQuantum.h"
namespace btas { typedef SpinQuantum Quantum; }; 

#include "MPSblas.h"
#include "Propagator.h"

using namespace btas;
using namespace mpsxx;

class Trotter;

class Walker {

   public:
   
      //Constructor copying an MPSstate
      Walker(const MPS<complex<double>,Quantum> &, double weight,int n_trot);
      
      //Constructor copying an entire Walker
      Walker(const Walker &walker);
      
      //Destructor
      virtual ~Walker();
      
      //Return the walker weight
      double gWeight() const;
      
      //Return the walker overlap
      complex<double> gOverlap() const;

      const MPS<complex<double>,Quantum> &gPsiW() const;

      MPS<complex<double>,Quantum> &gPsiW();

      //Set the overlap with the trial wfn
      void sOverlap(const MPS<complex<double>,Quantum> &Psi0);

      void sEL(const MPO<complex<double>,Quantum> &O,const MPS<complex<double>,Quantum> &Psi0);

      void sEL(complex<double> );

      void sVL(const Trotter &trotter,const MPS<complex<double>,Quantum> &Psi0);

      int gn_trot() const;

      complex<double> gEL() const;

      const std::vector< complex<double> > &gVL() const;

      complex<double> gVL(int,int) const;

      void multWeight(double);

      void sWeight(double);

      void propagate(const Propagator &P);
 
   private:
   
      //The MPS state
      MPS<complex<double>,Quantum> PsiW;

      //nr of trotter terms
      int n_trot;
      
      //The walker overlap with the trial wfn
      complex<double> overlap;

      //!local energy
      complex<double> EL;

      //!local auxiliary operators: <PsiT|v|phi>/<PsiT|phi>
      std::vector< complex<double> > VL;
      
      //The walker weight
      double weight;
      
};

#endif
