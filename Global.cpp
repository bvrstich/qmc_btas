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

#include "include.h"

Random Global::RN;

int Global::Lx;
int Global::Ly;
int Global::L;
double Global::J2;
int Global::j2;
DArray<2> Global::J;
int Global::d;
int Global::D;
bool Global::pbc;


//!function which generates random complex numbers uniformly on a square of side 2
template<>
complex<double> Global::rgen(){ 

   return complex<double>(2.0*RN() - 1.0,2.0*RN() - 1.0); 

}

//!function which generates uniform random numbers between [-1:1]
template<>
double Global::rgen(){ 

   return 2.0*RN() - 1.0;

}

/**
 * initialize the global variables
 * @param Lx_in input x dim
 * @param Ly_in input y dim
 * @param j2_in strength of the J2 coupling parameter
 * @param d_in physical dimension 
 * @param D_in virtual dimension of the trial
 * @param pbc true if periodic boundary conditions are assumed
 */
void Global::init(int Lx_in,int Ly_in,int j2_in,int d_in,int D_in,bool pbc_in){

   Lx = Lx_in;
   Ly = Ly_in;
   L = Lx*Ly;

   j2 = j2_in;
   J2 = 0.1*(double)j2;

   d = d_in;
   D = D_in;

   pbc = pbc_in;

   J.resize(L,L);

   //coupling matrix:
   coupling::J1J2_2D(pbc,J2,J);

}

/**
 * @return physical dimension
 */
int Global::gd(){

   return d;

}

/**
 * @return max virtual bond dimension of trial
 */
int Global::gD(){

   return D;

}

/**
 * @return max virtual bond dimension of trial
 */
const DArray<2> &Global::gJ(){

   return J;

}

/**
 * @return the x dim
 */
int Global::gLx() {

   return Lx;

}

/**
 * @return the y dim
 */
int Global::gLy() {

   return Ly;

}

/**
 * @return the length of the chain
 */
int Global::gL() {

   return L;

}

/**
 * @return integer value of the coupling parameter
 */
int Global::gj2() {

   return j2;

}


/**
 * @return J2 double value of the coupling paramter
 */
double Global::gJ2() {

   return J2;

}
