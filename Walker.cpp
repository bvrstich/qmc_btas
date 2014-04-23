#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <complex>
#include <vector>

#include "include.h"

using std::cout;
using std::ifstream;
using std::ofstream;
using std::endl;

/**
 * empty constructor
 */
Walker::Walker() : std::vector< TArray<complex<double>,1> >(Global::gL()){ }

/**
 * construct a Walker object
 * @param n_trot number of trotter terms
 */
Walker::Walker(int n_trot_in) : std::vector< TArray<complex<double>,1> >(Global::gL()){

   int L = Global::gL(); 
   int d = Global::gd(); 

   this->n_trot = n_trot_in;

   for(int i = 0;i < L;++i)
      (*this)[i].resize(d);

   Vxyz.resize(3*L);

   for(int i = 0;i < L;++i)
      for(int r = 0;r < 3;++r)//regular vector
         Vxyz[3*i + r].resize(d);

   VL.resize(3*n_trot);

}

/**
 * copy constructor
 */
Walker::Walker(const Walker &walker) : std::vector< TArray<complex<double>,1> >(walker) {

   this->weight = walker.gWeight();

   n_trot = walker.gn_trot();
   overlap = walker.gOverlap();
   EL = walker.gEL();
   VL = walker.gVL();
   Vxyz = walker.Vxyz;

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
 * calculate the overlap between the D=1 walker and a large D MPS
 */
complex<double> Walker::calc_overlap(const MPS< complex<double> > &mps) const {

   int L = Global::gL();
   int d = Global::gd();

   complex<double> one(1.0,0.0);
   complex<double> zero(0.0,0.0);

   blas::gemv(CblasRowMajor, CblasTrans, d, d, one, mps[0].data(), d, (*this)[0].data(), 1, zero, Global::LO[0].data(), 1);

   for(int i = 1;i < L;++i){

      int Ldim = mps[i].shape(1);
      int Rdim = mps[i].shape(2);

      blas::gemv(CblasRowMajor, CblasTrans, d, Global::gemv_list[i], one, mps[i].data(), Global::gemv_list[i], (*this)[i].data(), 1, zero, Global::loc[i].data(), 1);

      blas::gemv(CblasRowMajor, CblasTrans, Ldim, Rdim, one, Global::loc[i].data(), Rdim, Global::LO[i - 1].data(), 1, zero, Global::LO[i].data(), 1);

   }

   return Global::LO[L - 1](0);

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

   int L = Global::gL();
   int d = Global::gd();

   //start by constructing renormalized operator for overlap of Psi0 and walker state

   //last site
   blas::gemv(CblasRowMajor, CblasTrans, d, d, one, Psi0[L - 1].data(), d, (*this)[L - 1].data(), 1, zero, Global::RO[L - 1].data(), 1);

   for(int i = L - 2;i > 0;--i){

      int Ldim = Psi0[i].shape(1);
      int Rdim = Psi0[i].shape(2);

      blas::gemv(CblasRowMajor, CblasTrans, d, Global::gemv_list[i], one, Psi0[i].data(), Global::gemv_list[i], (*this)[i].data(), 1, zero, Global::loc[i].data(), 1);

      blas::gemv(CblasRowMajor, CblasNoTrans, Ldim, Rdim, one, Global::loc[i].data(), Rdim, Global::RO[i + 1].data(), 1, zero, Global::RO[i].data(), 1);

   }

   //first site = 0
   for(int r = 0;r < 3;++r){

      blas::gemv(CblasRowMajor, CblasTrans, d, d, one, Psi0[0].data(), d, Vxyz[r].data(), 1, zero, Global::loc[0].data(), 1);

      Global::auxvec[0][r] = blas::dot(d,Global::loc[0].data(),1,Global::RO[1].data(),1);

   }

   //make the left operator
   blas::gemv(CblasRowMajor, CblasTrans, d, d, one, Psi0[0].data(), d, (*this)[0].data(), 1, zero, Global::LO[0].data(), 1);

   //calculate Global::auxvec for middle sites
   for(int i = 1;i < L - 1;++i){

      int Ldim = Psi0[i].shape(1);
      int Rdim = Psi0[i].shape(2);

      for(int r = 0;r < 3;++r){

         blas::gemv(CblasRowMajor, CblasTrans, d, Global::gemv_list[i], one, Psi0[i].data(), Global::gemv_list[i], Vxyz[i*3 + r].data(), 1, zero, Global::loc[i].data(), 1);

         blas::gemv(CblasRowMajor, CblasTrans, Ldim, Rdim, one, Global::loc[i].data(), Rdim, Global::LO[i - 1].data(), 1, zero, Global::LO[i].data(), 1);

         Global::auxvec[i][r] = blas::dot(Rdim,Global::LO[i].data(),1,Global::RO[i + 1].data(),1);

      }

      //make the left operator
      blas::gemv(CblasRowMajor, CblasTrans, d, Global::gemv_list[i], one, Psi0[i].data(), Global::gemv_list[i], (*this)[i].data(), 1, zero, Global::loc[i].data(), 1);

      blas::gemv(CblasRowMajor, CblasTrans, Ldim, Rdim, one, Global::loc[i].data(), Rdim, Global::LO[i - 1].data(), 1, zero, Global::LO[i].data(), 1);

   }

   //last site = L-1
   for(int r = 0;r < 3;++r){

      blas::gemv(CblasRowMajor, CblasTrans, d, d, one, Psi0[L - 1].data(), d, Vxyz[(L - 1)*3 + r].data(), 1, zero, Global::loc[L - 1].data(), 1);

      Global::auxvec[L - 1][r] = blas::dot(d,Global::LO[L - 2].data(),1,Global::loc[L - 1].data(),1);

   }

   //Done with the Sx, Sy and Sz on every site: ready to construct auxiliary expectation values quickly!
   for(int k = 0;k < n_trot;++k)
      for(int r = 0;r < 3;++r){

         VL[r*n_trot + k] = Global::auxvec[0][r] * trotter.gV()(k,0);

         for(int i = 1;i < L;++i)
            VL[r*n_trot + k] += Global::auxvec[i][r] * trotter.gV()(k,i);

         VL[r*n_trot + k] /= overlap;

      }

}

/**
 * muliply the weight by a factor
 */
void Walker::multWeight(double factor){

   weight *= factor; 

}

/**
 * set new weight
 */
void Walker::sWeight(double new_weight){

   weight = new_weight;

}

/**
 * Apply the propagator (a D=1 MPO) to the current state of the walker. The walker is changed when calling this function
 * @param P the propagator
 */
void Walker::propagate(const Propagator &P){

   int d = Global::gd();

   complex<double> one(1.0,0.0);
   complex<double> zero(0.0,0.0);

   enum {j,k,l,m};

   TArray<complex<double>,1> tmp(d);

   for(int i = 0;i < this->size();++i){

      blas::gemv(CblasRowMajor,CblasNoTrans, d, d, one, P[i].data(), d, (*this)[i].data(), 1, zero, tmp.data(), 1);

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

   int d = Global::gd();

   for(int i = 0;i < this->size();++i){

      double nrm = blas::nrm2(d, (*this)[i].data(), 1);

      blas::scal(d, static_cast< complex<double> >(1.0/nrm), (*this)[i].data(), 1);

   }

}

/** 
 * set the Local Energy: overlap has to be set first!
 * @param Psi0 trial wavefunction 
 */
void Walker::sEL(const MPS< complex<double> > &Psi0){

   EL = Heisenberg::energy(Psi0,*this) / overlap;

}

/**
 * fill the Sx,Sy and Sz rotated walker vectors
 */
void Walker::fill_xyz() {

   complex<double> one(1.0,0.0);
   complex<double> zero(0.0,0.0);

   int L = this->size();
   int d = Global::gd();

   for(int i = 0;i < L;++i){

      blas::gemv(CblasRowMajor,CblasNoTrans, d, d, one, Heisenberg::gSx().data(), d, (*this)[i].data(), 1, zero, Vxyz[i*3].data(), 1);
      blas::gemv(CblasRowMajor,CblasNoTrans, d, d, one, Heisenberg::gSy().data(), d, (*this)[i].data(), 1, zero, Vxyz[i*3 + 1].data(), 1);
      blas::gemv(CblasRowMajor,CblasNoTrans, d, d, one, Heisenberg::gSz().data(), d, (*this)[i].data(), 1, zero, Vxyz[i*3 + 2].data(), 1);

   }

}

/**
 * @return the xyz vector on site [i] of type r = 0,1,2 (x,y,z)
 */
const TArray<complex<double>,1> &Walker::gVxyz(int i,int r) const {

   return Vxyz[3*i + r];

}

/** 
 * fill the walker from file
 * @param filename location of inputfile
 */
void Walker::read(const char *filename){

   ifstream in(filename);

   for(int i = 0;i < this->size();++i)
      for(int s = 0;s < (*this)[i].shape(0);++s)
         in >> i >> s >> (*this)[i](s);

}

/**
 * copy the necessary information of a Walker object into this
 * @param walker_copy input Walker
 */
void Walker::copy_essential(const Walker &walker_copy){

   EL = walker_copy.gEL();
   overlap = walker_copy.gOverlap();
   weight = walker_copy.gWeight();

   int L = Global::gL();
   int d = Global::gd();

   for(int i = 0;i < L;++i)
      blas::copy(d, walker_copy[i].data(), 1, (*this)[i].data(), 1);

}
