#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <complex>
#include <vector>

#include "include.h"

/**
 * constructor, diagonalizes the coupling matrix J to make the auxiliary field operators V(r,k)
 * @param J coupling matrix
 * @param dtau time step
 */
Trotter::Trotter(const DArray<2> &J,double dtau){

   L = J.shape(0);

   this->dtau = dtau;

   //diagonalize the coupling matrix
   DArray<1> Jeig;
   DArray<2> Uv;

   Syev('V','U',J,Jeig,Uv);

   n_trot = 0;

   for(int i = 0;i < L;++i){

      if(fabs(Jeig(i)) > 1.0e-14)
         n_trot++;

   }

   //now transform the elements with dtau and the eigenvalues for the propagator to form the transformation V
   V.resize(n_trot,L);

   int k = 0;

   for(int i = 0;i < L;++i){

      if(fabs(Jeig(i)) > 1.0e-14){

         complex<double> tmp =  std::sqrt( (complex<double>)-Jeig(i) * dtau);

         for(int j = 0;j < L;++j)
            V(k,j) = tmp * Uv(i,j);

         ++k;

      }

   }

}

/**
 * destructor
 */
Trotter::~Trotter(){

}

/**
 * @return the timestep
 */
double Trotter::gdtau() const {

   return dtau;

}
      
/**
 * @return the auxiliary field matrix
 */
const ZArray<2> &Trotter::gV() const {

   return V;

}

/**
 * @return the number of trotter terms
 */
int Trotter::gn_trot() const {

   return n_trot;

}

/**
 * @return the length of the chain
 */
int Trotter::gL() const {

   return L;

}
