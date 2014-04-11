#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <complex>

using std::cout;
using std::endl;
using std::vector;
using std::complex;
using std::ofstream;
using std::ifstream;

#include "include.h"

/** 
 * constructor: just sets the length of the vector, nothing is allocates or initialized
 * @param L_in length of the chain
 */
template<typename T>
MPS<T>::MPS(int L_in) : vector< TArray<T,3> >(L_in) { }

/** 
 * standard constructor: just takes in
 * @param L_in length of the chain
 * @param d_in physical dimension of the sites
 * @param D_in virtual max bond dimension
 * allocates the tensors and fills them randomly
 */
template<typename T>
MPS<T>::MPS(int L_in,int d_in,int D_in) : vector< TArray<T,3> >(L_in) {

   D = D_in;
   d = d_in;

   vector<int> vdim(L_in + 1);

   vdim[0] = 1;

   for(int i = 1;i < L_in;++i){

      int tmp = vdim[i - 1] * d;

      if(tmp < D)
         vdim[i] = tmp;
      else 
         vdim[i] = D;

   }

   vdim[L_in] = 1;

   for(int i = L_in - 1;i > 0;--i){

      int tmp = vdim[i + 1] * d;

      if(tmp < vdim[i])
         vdim[i] = tmp;

   }

   for(int i = 0;i < this->size();++i){

      (*this)[i].resize(d,vdim[i],vdim[i+1]);
      (*this)[i].generate(Global::rgen<T>);

   }

}

/** 
 * @param filename inputfile
 */
template<typename T>
MPS<T>::MPS(const char *filename) : vector< TArray<T,3> >() {

   ifstream in(filename);

   int L;

   in >> L >> D >> d;

   this->resize(L);

   vector<int> vdim(L + 1);

   for(int i = 0;i <= L;++i)
      in >> i >> vdim[i];

   int dim;

   for(int i = 0;i < L;++i){

      (*this)[i].resize(d,vdim[i],vdim[i+1]);

      in >> i >> dim;

      int teller;

      for(int s = 0;s < d;++s)
         for(int k = 0;k < vdim[i + 1];++k)
            for(int j = 0;j < vdim[i];++j)
               in >> i >> teller >> (*this)[i](s,j,k);
            

   }

}


/**
 * copy constructor
 */
template<typename T>
MPS<T>::MPS(const MPS<T> &mps_copy) : vector< TArray<T,3> >(mps_copy) {

   D = mps_copy.gD();

}

/**
 * empty destructor
 */
template<typename T>
MPS<T>::~MPS(){ }

/**
 * @return virtual dimension of the MPS
 */
template<typename T>
int MPS<T>::gD() const {

   return D;

}

/**
 * scale the MPS with a constant factor
 * @param alpha scalingfactor
 */
template<>
void MPS<double>::scal(double alpha){

   int sign;

   if(alpha > 0)
      sign = 1;
   else
      sign = -1;

   alpha = pow(fabs(alpha),1.0/(double)this->size());

   Scal(sign * alpha,(*this)[0]);

   for(int i = 1;i < this->size();++i)
      Scal(alpha,(*this)[i]);

}

/**
 * scale the MPS with a constant factor
 * @param alpha scalingfactor
 */
template<>
void MPS< complex<double> >::scal(complex<double> alpha){

   alpha = pow(fabs(alpha),1.0/(complex<double>)this->size());

   Scal(alpha,(*this)[0]);

   for(int i = 1;i < this->size();++i)
      Scal(alpha,(*this)[i]);

}

template MPS<double>::MPS(int,int,int);
template MPS< complex<double> >::MPS(int,int,int);

template MPS<double>::MPS(const char *);
template MPS< complex<double> >::MPS(const char *);

template MPS<double>::MPS(int);
template MPS< complex<double> >::MPS(int);

template MPS<double>::MPS(const MPS<double> &);
template MPS< complex<double> >::MPS(const MPS< complex<double> > &);

template MPS<double>::~MPS();
template MPS< complex<double> >::~MPS();

template int MPS<double>::gD() const;
template int MPS< complex<double> >::gD() const;
