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
 * @param dtau time step
 */
Trotter::Trotter(double dtau){

   int L = Global::gL();
   int d = Global::gd();

   this->dtau = dtau;

   //diagonalize the coupling matrix
   DArray<1> Jeig;
   DArray<2> Uv;

   Syev('V','U',Global::gJ(),Jeig,Uv);

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
            Mx[i](j,k) = U(j,i) * std::conj(U(k,i));

   }

   //Sy
   Heev('V', 'U', Heisenberg::gSy(), eig, U);

   My.resize(d);

   for(int i = 0;i < d;++i){

      My[i].resize(d,d);

      for(int j = 0;j < d;++j)
         for(int k = 0;k < d;++k)
            My[i](j,k) = U(j,i) * std::conj(U(k,i));

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

   int d = Global::gd();

   dtau = trot_copy.gdtau();
   V = trot_copy.gV();
   n_trot = trot_copy.gn_trot();

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
