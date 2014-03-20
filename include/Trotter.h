#ifndef TROTTER_H
#define TROTTER_H

#include <iostream>
#include <iomanip>

#include <btas/DENSE/DArray.h>
#include <btas/DENSE/ZArray.h>

class Trotter {

   public:
   
      //Constructor
      Trotter(const DArray<2> &J, double dtau);
      
      //Destructor
      ~Trotter();
      
      double gdtau() const;
      
      const ZArray<2> &gV() const;
      
      int gn_trot() const;

      int gL() const;

   private:
   
      //The chain length
      int L;
      
      //!the transformation matrix
      ZArray<2> V;

      //!the number of non-zero eigenvalues of the coupling matrix: equals the number of trotter product terms
      int n_trot;

      //The time step
      double dtau;
      
};

#endif
