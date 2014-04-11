#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <complex>
#include <vector>

#include "include.h"

using std::cout;
using std::endl;

/**
 * construct a Walker object
 * @param L length of the chain
 * @param d physical dimension
 * @param weight input weight
 * @param n_trot number of trotter terms
 */
Walker::Walker(int L,int d,double weight,int n_trot) : std::vector< ZArray<1> >(L){

   this->weight = weight;
   this->n_trot = n_trot;

   for(int i = 0;i < L;++i)
      (*this)[i].resize(d);

   VL.resize(3*n_trot);

   auxvec.resize(n_trot);

   for(int i = 0;i < n_trot;++i)//for x,y and z components
      auxvec[i].resize(3);

}

/**
 * copy constructor
 */
Walker::Walker(const Walker &walker) : std::vector< ZArray<1> >(walker) {

   weight = walker.gWeight();
   n_trot = walker.gn_trot();
   overlap = walker.gOverlap();
   EL = walker.gEL();
   VL = walker.gVL();
   auxvec = walker.gauxvec();

}

/**
 * destructor
 */
Walker::~Walker(){ }

/** 
 * @return the number of trotter terms
 */
int Walker::gn_trot() const {

   return n_trot;

}

/** 
 * @return the weight corresponding to the walker
 */
double Walker::gWeight() const{

   return weight; 

}

/** 
 * @return the overlap of the walker with the Trial
 */
complex<double> Walker::gOverlap() const{

   return overlap; 

}

/** 
 * @return the local energy
 */
complex<double> Walker::gEL() const{

   return EL; 

}

/** 
 * @return the shifts, i.e. the vector containing the projected expectation values of the auxiliary field operators
 */
const std::vector< complex<double> > &Walker::gVL() const{

   return VL; 

}

/** 
 * @return intermediate object for calcuation of expectation values of auxiliary fields
 */
const std::vector< std::vector< complex<double> > > &Walker::gauxvec() const{

   return auxvec; 

}


/**
 * calculate the overlap between the D=1 walker and a large D MPS
 */
complex<double> Walker::calc_overlap(const MPS< complex<double> > &mps) const {

   complex<double> one(1.0,0.0);
   complex<double> zero(0.0,0.0);

   ZArray<2> E;

   Gemv(CblasTrans,one,mps[0],(*this)[0],zero,E);

   ZArray<2> E_tmp;
   ZArray<2> tmp;

   for(int i = 1;i < mps.size();++i){

      tmp.clear();
      Gemv(CblasTrans,one,mps[i],(*this)[i],zero,tmp);

      E_tmp.clear();
      Gemm(CblasNoTrans,CblasNoTrans,one,E,tmp,zero,E_tmp);

      E = std::move(E_tmp);

   }

   return E(0,0);

}

/**
 * calculate the overlap with the trial, Psi0
 * @param Psi0 the trial
 */
void Walker::sOverlap(const MPS<complex<double> > &Psi0){

   complex<double> prev_overlap = overlap;

   //just for stability, doesn't do anything physical
   this->normalize();

   overlap = this->calc_overlap(Psi0);

   //phase free projection
   weight *= std::max(0.0,cos(std::arg(overlap/prev_overlap)));

}

/** 
 * set the Local Energy with a number
 */
void Walker::sEL(complex<double> EL_in){

   EL = EL_in;

}

/** 
 * set the shifts, the auxiliary field projected expectation values: <PsiT|v|PsiW>/<PsiT|PsiW> --> overlap has to be set
 * @param trotter object containing the information about the auxiliary field operators
 * @param Psi0 trial wavefunction 
 */
void Walker::sVL(const Trotter &trotter,const MPS< complex<double> > &Psi0){

   complex<double> one(1.0,0.0);
   complex<double> zero(0.0,0.0);

   int L = this->size();
   int d = trotter.gd();

   //start by constructing renormalized operator for overlap of Psi0 and walker state
   std::vector< ZArray<2> > ro(L - 1);

   //last site
   Gemv(CblasTrans,one,Psi0[L - 1],(*this)[L - 1],zero,ro[L - 2]);

   ZArray<2> tmp;

   for(int i = L - 2;i > 0;--i){

      tmp.clear();
      Gemv(CblasTrans,one,Psi0[i],(*this)[i],zero,tmp);

      Gemm(CblasNoTrans,CblasNoTrans,one,tmp,ro[i],zero,ro[i - 1]);

   }

   //now fill the auxvec
   ZArray<1> vec(d);

   //first site = 0

   //x
   Gemv(CblasNoTrans,one,trotter.gSx(),(*this)[0],zero,vec);

   tmp.clear();
   Gemv(CblasTrans,one,Psi0[0],vec,zero,tmp);

   auxvec[0][0] = Dot(tmp,ro[0]);

   //y
   Gemv(CblasNoTrans,one,trotter.gSy(),(*this)[0],zero,vec);

   tmp.clear();
   Gemv(CblasTrans,one,Psi0[0],vec,zero,tmp);

   auxvec[0][1] = Dot(tmp,ro[0]);

   //z
   Gemv(CblasNoTrans,one,trotter.gSz(),(*this)[0],zero,vec);

   tmp.clear();
   Gemv(CblasTrans,one,Psi0[0],vec,zero,tmp);

   auxvec[0][2] = Dot(tmp,ro[0]);

   //make the left operator
   ro[0].clear();
   Gemv(CblasTrans,one,Psi0[0],(*this)[0],zero,ro[0]);

   ZArray<2> tmp2;

   //calculate auxvec for middle sites
   for(int i = 1;i < L - 1;++i){

      //x
      Gemv(CblasNoTrans,one,trotter.gSx(),(*this)[i],zero,vec);

      tmp.clear();
      Gemv(CblasTrans,one,Psi0[i],vec,zero,tmp);

      tmp2.clear();
      Gemm(CblasNoTrans,CblasNoTrans,one,ro[i - 1],tmp,zero,tmp2);

      auxvec[i][0] = Dot(tmp2,ro[i]);

      //y
      Gemv(CblasNoTrans,one,trotter.gSy(),(*this)[i],zero,vec);

      tmp.clear();
      Gemv(CblasTrans,one,Psi0[i],vec,zero,tmp);

      tmp2.clear();
      Gemm(CblasNoTrans,CblasNoTrans,one,ro[i - 1],tmp,zero,tmp2);

      auxvec[i][1] = Dot(tmp2,ro[i]);

      //z
      Gemv(CblasNoTrans,one,trotter.gSz(),(*this)[i],zero,vec);

      tmp.clear();
      Gemv(CblasTrans,one,Psi0[i],vec,zero,tmp);

      tmp2.clear();
      Gemm(CblasNoTrans,CblasNoTrans,one,ro[i - 1],tmp,zero,tmp2);

      auxvec[i][2] = Dot(tmp2,ro[i]);

      //make the left operator
      tmp.clear();
      Gemv(CblasTrans,one,Psi0[i],(*this)[i],zero,tmp);

      ro[i].clear();
      Gemm(CblasNoTrans,CblasNoTrans,one,ro[i - 1],tmp,zero,ro[i]);

   }

   //last site = 0

   //x
   Gemv(CblasNoTrans,one,trotter.gSx(),(*this)[L - 1],zero,vec);

   tmp.clear();
   Gemv(CblasTrans,one,Psi0[L - 1],vec,zero,tmp);

   auxvec[L - 1][0] = Dot(tmp,ro[L - 2]);

   //y
   Gemv(CblasNoTrans,one,trotter.gSy(),(*this)[L - 1],zero,vec);

   tmp.clear();
   Gemv(CblasTrans,one,Psi0[L - 1],vec,zero,tmp);

   auxvec[L - 1][1] = Dot(tmp,ro[L - 2]);

   //z
   Gemv(CblasNoTrans,one,trotter.gSz(),(*this)[L - 1],zero,vec);

   tmp.clear();
   Gemv(CblasTrans,one,Psi0[L - 1],vec,zero,tmp);

   auxvec[L - 1][2] = Dot(tmp,ro[L - 2]);

   //Done with the Sx, Sy and Sz on every site: ready to construct auxiliary expectation values quickly!

}

/**
 * muliply the weight by a factor
 */
void Walker::multWeight(double factor){

   weight *= factor; 

}

/**
 * muliply the weight by a factor
 */
void Walker::sWeight(double new_weight){

   weight = new_weight;

}

/**
 * Apply the propagator (a D=1 MPO) to the current state of the walker. The walker is changed when calling this function
 * @param P the propagator
 */
void Walker::propagate(const Propagator &P){

   ZArray<1> tmp;

   complex<double> one(1.0,0.0);
   complex<double> zero(0.0,0.0);

   enum {j,k,l,m};

   for(int i = 0;i < this->size();++i){

      tmp.clear();

      Gemv(CblasNoTrans,one,P[i],(*this)[i],zero,tmp);

      (*this)[i] = std::move(tmp);

   }

}

/**
 * @return the projected expectation value of the shift corresponding to trotter indices k and r
 */
complex<double> Walker::gVL(int k,int r) const {

   return VL[r*n_trot + k];

}

/** 
 * mostly for debugging purposes, fill the vector with random complex numbers
 */
void Walker::fill_Random(){

   for(int i = 0;i < this->size();++i)
      (*this)[i].generate(Global::rgen< complex<double> >);

}

/** 
 * normalize the walker state, will not physically effect the afqmc walk, just keeps the numbers from blowing up.
 */
void Walker::normalize(){

   for(int i = 0;i < this->size();++i)
      Normalize((*this)[i]);

}

/** 
 * set the Local Energy: overlap has to be set first!
 * @param O mpo containing Hamiltonian
 * @param Psi0 trial wavefunction 
 */
/*
   void Walker::sEL(const MPO<complex<double>,Quantum> &O,const MPS<complex<double>,Quantum> &Psi0){

   EL = mpsxx::inprod(mpsxx::Left,PsiW,O,Psi0) / overlap;

   }
 */
