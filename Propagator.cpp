#include <iostream>
#include <cmath>
#include <vector>
#include <complex>
#include <fstream>

using std::ostream;
using std::ofstream;
using std::ifstream;
using std::cout;
using std::endl;
using std::ios;
using std::vector;
using std::complex;

#include "include.h"

/**
 * standard constructor: sets the size of vector, and initializes the QSZArray objects to the correct quantumnumbers and dimensions
 * @param L size
 */
Propagator::Propagator(int L_in,int d_in) : vector< QSZArray<2,Quantum> > (L_in) {

   L = L_in;
   d = d_in;

   x = 0.0;
   k = 0;
   r = 0;
  
   //physical indices
   Qshapes<Quantum> qp;

   for(int i = 0;i < d;++i)
      qp.push_back(Quantum::zero());

   Dshapes dp;

   for(int i = 0;i < d;++i)
      dp.push_back(1);

   for(int i = 0;i < L;++i)
      (*this)[i].resize(Quantum::zero(),make_array(qp,-qp),make_array(dp,dp));
 
}

/**
 * copy constructor
 * @param prop_copy input Propagator object
 */
Propagator::Propagator(const Propagator &prop_copy) : vector< QSZArray<2,Quantum> > (prop_copy) {

   L = prop_copy.gL();
   d = prop_copy.gd();

   x = prop_copy.gx();
   k = prop_copy.gk();
   r = prop_copy.gr();
   
}

/**
 * destructor
 */
Propagator::~Propagator() { }

/**
 * @return the auxiliary field variable x
 */
complex<double> Propagator::gx() const {

   return x;

}

/**
 * @return the trotter index k
 */
int Propagator::gk() const {

   return k;

}

/**
 * @return the type of operator r (0=x,1=y,2=z)
 */
int Propagator::gr() const {

   return r;

}

/**
 * @return the length of the chain
 */
int Propagator::gL() const {

   return L;

}

/**
 * @return the physical dimension
 */
int Propagator::gd() const {

   return d;

}
/**
 * set the variables to fill the propagator on
 */
void Propagator::set(complex<double> x_in,int k_in,int r_in){

   x = x_in;
   k = k_in;
   r = r_in;

}

/**
 * fill the propagator using a Trotter input object
 */
void Propagator::fill(const Trotter &trotter) {

   if(r == 0){//x

      for(int site = 0;site < L;++site){

         (*this)[site] = trotter.gMx(0);

         double m = 0.5 * ( d - 1.0 );

         Scal( exp(x * trotter.gV()(k,site) * m ) , (*this)[site] );

         for(int i = 1;i < d;++i){

            m++;

            Axpy( exp(  x * trotter.gV()(k,site) * m ) , trotter.gMx(i), (*this)[site] );

         }

      }

   }
   else if(r == 1){//y

      for(int site = 0;site < L;++site){

         (*this)[site] = trotter.gMy(0);

         double m = 0.5 * ( d - 1.0 );

         Scal( exp(x * trotter.gV()(k,site) * m ) , (*this)[site] );

         for(int i = 1;i < d;++i){

            m++;

            Axpy( exp(  x * trotter.gV()(k,site) * m ) , trotter.gMy(i), (*this)[site] );

         }

      }

   }
   else{//z

      for(int site = 0;site < L;++site){

         (*this)[site] = trotter.gMz(0);

         double m = 0.5 * ( d - 1.0 );

         Scal( exp(x * trotter.gV()(k,site) * m ) , (*this)[site] );

         for(int i = 1;i < d;++i){

            m++;

            Axpy( exp(  x * trotter.gV()(k,site) * m ) , trotter.gMz(i), (*this)[site] );

         }

      }

   }

}
