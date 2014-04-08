#ifndef TROTTER_H
#define TROTTER_H

#include <iostream>
#include <iomanip>
#include <complex>
#include <vector>

#include <btas/DENSE/DArray.h>
#include <btas/DENSE/ZArray.h>

class Trotter {

   public:
   
      //Constructor
      Trotter(int,const DArray<2> &, double);

      //copy constructor
      Trotter(const Trotter &trot_copy);
      
      //Destructor
      ~Trotter();
      
      double gdtau() const;
      
      const ZArray<2> &gV() const;
      
      int gn_trot() const;

      int gL() const;

      int gd() const;

      const MPO<complex<double>,Quantum> &gV_op(int,int) const;

      void fill_prop(MPO<complex<double>,Quantum> &P,int k,int r,complex<double> x);

      const QSZArray<2,Quantum> &gMx(int) const;

      const QSZArray<2,Quantum> &gMy(int) const;

      const QSZArray<2,Quantum> &gMz(int) const;


   private:
   
      //The chain length
      int L;
      
      //!the transformation matrix
      ZArray<2> V;

      //!physical dimension
      int d;

      //!mpo's corresponding to the auxiliary field operators
      std::vector< MPO<complex<double>,Quantum > > V_op;

      //!the number of non-zero eigenvalues of the coupling matrix: equals the number of trotter product terms
      int n_trot;

      //!The time step
      double dtau;

      //!eigenvector of Sx |Sx><Sx|
      std::vector< QSZArray<2,Quantum> > Mx;
      
      //!eigenvector of Sy |Sy><Sy|
      std::vector< QSZArray<2,Quantum> > My;

      //!eigenvector of Sz |Sz><Sz|
      std::vector< QSZArray<2,Quantum> > Mz;

};

#endif
