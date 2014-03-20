#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <complex>
#include <vector>

#include "include.h"

/**
 * construct a Walker object
 * @param PsiW input walker MPS
 * @param weight input weight
 * @param n_trot number of trotter terms
 */
Walker::Walker(const MPS<complex<double>,Quantum> &PsiW, double weight,int n_trot){

   this->PsiW = PsiW;
   this->weight = weight;
   this->n_trot = n_trot;

   VL.resize(3*n_trot);

}

/**
 * copy constructor
 */
Walker::Walker(const Walker &walker){

   PsiW = walker.gPsiW();
   weight = walker.gWeight();
   n_trot = walker.gn_trot();
   overlap = walker.gOverlap();
   EL = walker.gEL();
   VL = walker.gVL();

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
 * @return the walker MPS : const version
 */
const MPS<complex<double>,Quantum> &Walker::gPsiW() const {
   
   return PsiW; 
   
}

/**
 * @return the walker MPS 
 */
MPS<complex<double>,Quantum> &Walker::gPsiW(){
   
   return PsiW; 
   
}

/**
 * calculate the overlap with the trial, Psi0
 * @param Psi0 the trial
 */
void Walker::sOverlap(const MPS<complex<double>,Quantum> &Psi0){

   complex<double> prev_overlap = overlap;

   normalize(PsiW);

   overlap = mpsxx::dot(mpsxx::Left,PsiW,Psi0);

   //phase free projection
   weight *= std::max(0.0,cos(std::arg(overlap/prev_overlap)));

}

/** 
 * set the Local Energy: overlap has to be set first!
 * @param O mpo containing Hamiltonian
 * @param Psi0 trial wavefunction 
 */
void Walker::sEL(const MPO<complex<double>,Quantum> &O,const MPS<complex<double>,Quantum> &Psi0){

   EL = mpsxx::inprod(mpsxx::Left,PsiW,O,Psi0) / overlap;

}

/** 
 * set the Local Energy with a number
 */
void Walker::sEL(complex<double> EL_in){

   EL = EL_in;

}
