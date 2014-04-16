#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <complex>
#include <vector>

using std::cout;
using std::endl;

#include "include.h"

/**
 * constructor, diagonalizes the coupling matrix J to make the auxiliary field operators V(r,k)
 * @param d physical dimension
 * @param J coupling matrix
 * @param dtau time step
 */
Trotter::Trotter(int d_in,const DArray<2> &J,double dtau){

   L = J.shape(0);

   this->dtau = dtau;
   this->d = d_in;

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
            V(k,j) = tmp * Uv(j,i);

         ++k;

      }

   }

   //Sx
   DArray<1> eig;
   ZArray<2> U;

   Heev('V', 'U', Heisenberg::gSx(), eig, U);

   Mx.resize(d);

   for(int i = 0;i < d;++i){

      Mx[i].resize(d,d);

      for(int j = 0;j < d;++j)
         for(int k = 0;k < d;++k)
            Mx[i](j,k) = U(i,j) * std::conj(U(i,k));

   }

   //Sy
   Heev('V', 'U', Heisenberg::gSy(), eig, U);

   My.resize(d);

   for(int i = 0;i < d;++i){

      My[i].resize(d,d);

      for(int j = 0;j < d;++j)
         for(int k = 0;k < d;++k)
            My[i](j,k) = U(i,j) * std::conj(U(i,k));

   }

   //Sz
   Mz.resize(d);

   for(int i = 0;i < d;++i){

      Mz[i].resize(d,d);
      Mz[i] = 0.0;

      Mz[i](i,i) = 1.0;

   }

}

/** 
 * copy constructor
 */
Trotter::Trotter(const Trotter &trot_copy){

   dtau = trot_copy.gdtau();
   V = trot_copy.gV();
   n_trot = trot_copy.gn_trot();
   L = trot_copy.gL();
   d = trot_copy.gd();

   Mx.resize(d);
   My.resize(d);
   Mz.resize(d);

   for(int i = 0;i < d;++i){

      Mx[i] = trot_copy.gMx(i);
      My[i] = trot_copy.gMy(i);
      Mz[i] = trot_copy.gMz(i);

   }

}

/**
 * destructor
 */
Trotter::~Trotter(){ }

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

/**
 * @return the physical dimension
 */
int Trotter::gd() const {

   return d;

}

/**
 * @return the eigenvector 'matrix' |Sx><Sx| corresponding to the i'th eigenvalue ranked from small to large 
 */
const ZArray<2> &Trotter::gMx(int i) const {

   return Mx[i];

}

/**
 * @return the eigenvector 'matrix' |Sy><Sy| corresponding to the i'th eigenvalue ranked from small to large 
 */
const ZArray<2> &Trotter::gMy(int i) const {

   return My[i];

}

/**
 * @return the eigenvector 'matrix' |Sz><Sz| corresponding to the i'th eigenvalue ranked from small to large 
 */
const ZArray<2> &Trotter::gMz(int i) const {

   return Mz[i];

}
