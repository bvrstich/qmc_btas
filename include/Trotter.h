#ifndef TROTTER_H
#define TROTTER_H

#include <iostream>
#include <iomanip>
#include <complex>
#include <vector>

class Trotter {

   public:
   
      //Constructor
      Trotter(double);

      //copy constructor
      Trotter(const Trotter &trot_copy);
      
      //Destructor
      ~Trotter();
      
      double gdtau() const;
      
      const TArray<complex<double>,2> &gV() const;
      
      int gn_trot() const;

      const TArray<complex<double>,2> &gMx(int) const;

      const TArray<complex<double>,2> &gMy(int) const;

      const TArray<complex<double>,2> &gMz(int) const;

   private:
   
      //!the transformation matrix
      TArray<complex<double>,2> V;

      //!the number of non-zero eigenvalues of the coupling matrix: equals the number of trotter product terms
      int n_trot;

      //!The time step
      double dtau;

      //!eigenvector of Sx |Sx><Sx|
      std::vector< TArray<complex<double>,2> > Mx;
      
      //!eigenvector of Sy |Sy><Sy|
      std::vector< TArray<complex<double>,2> > My;

      //!eigenvector of Sz |Sz><Sz|
      std::vector< TArray<complex<double>,2> > Mz;

};

#endif
