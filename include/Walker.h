#ifndef WALKER_H
#define WALKER_H

#include <iostream>
#include <iomanip>
#include <complex>
#include <vector>

using std::complex;
using std::vector;

using namespace btas;

class Propagator;
class Trotter;

/**
 * class definition of Walker, made to describe the product state walkers. An array of L size-2 vector. Each site represents a rotation of the spin.
 */
class Walker : public vector< ZArray<1> > {

   public:
   
      //Constructor copying an MPSstate
      Walker(int,int,double weight,int n_trot);
      
      //Constructor copying an entire Walker
      Walker(const Walker &walker);
      
      //Destructor
      virtual ~Walker();

      //Return the walker weight
      double gWeight() const;
      
      //Return the walker overlap
      complex<double> gOverlap() const;

      int gn_trot() const;

      complex<double> gEL() const;

      const std::vector< complex<double> > &gVL() const;

      complex<double> gVL(int,int) const;

      const std::vector< std::vector< complex<double> > > &gauxvec() const;

      void fill_Random();

      complex<double> calc_overlap(const MPS< complex<double> > &mps) const;

      void propagate(const Propagator &P);

      void normalize();

      void fill_xyz();

      const ZArray<1> &gVxyz(int,int) const;

      //Set the overlap with the trial wfn
      void sOverlap(const MPS< complex<double> > &Psi0);

      void sVL(const Trotter &trotter,const MPS< complex<double> > &Psi0);

      void multWeight(double);

      void sWeight(double);
 
      void sEL(complex<double> );

     //void sEL(const MPO<complex<double>,Quantum> &O,const MPS<complex<double>,Quantum> &Psi0);

  private:
   
      //nr of trotter terms
      int n_trot;
      
      //The walker overlap with the trial wfn
      complex<double> overlap;

      //!local energy
      complex<double> EL;

      //!local auxiliary operators: <PsiT|v|phi>/<PsiT|phi>
      std::vector< complex<double> > VL;

      //!intermediate storage for calculation of auxiliary operator expectation values
      std::vector< std::vector< complex<double> > > auxvec;

      //!Sx,Sy and Sz operator application to walker
      std::vector< ZArray<1> > Vxyz;
      
      //The walker weight
      double weight;
      
};

#endif
