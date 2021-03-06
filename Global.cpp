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
TArray<double,2> Global::J;
int Global::d;
int Global::D;
bool Global::pbc;

std::vector< TArray<complex<double>,1> > Global::LO;
std::vector< TArray<complex<double>,1> > Global::RO;

vector<int> Global::gemv_list;

vector< TArray<complex<double>,2> > Global::loc;

vector< Walker > Global::backup_walker;

std::vector< std::vector< complex<double> > > Global::auxvec;

int Global::omp_num_threads;

/**
 * initialize the storage on dimensions of input MPS
 */
void Global::init_storage(const MPS< complex<double> > &mps){

   LO.resize(omp_num_threads * L);
   RO.resize(omp_num_threads * L);

   for(int proc = 0;proc < omp_num_threads;++proc)
      for(int i = 0;i < L;++i){

      LO[proc*L + i].resize(mps[i].shape(2));
      RO[proc*L + i].resize(mps[i].shape(1));

   }

   //make the dimension list for the gemv routine
   gemv_list.resize(L);

   for(int i = 0;i < L;++i)
      gemv_list[i] = mps[i].shape(1)*mps[i].shape(2);

   loc.resize(omp_num_threads * L);

   for(int proc = 0;proc < omp_num_threads;++proc)
      for(int i = 0;i < L;++i)
         loc[proc*L + i].resize(mps[i].shape(1),mps[i].shape(2));

   backup_walker.resize(omp_num_threads);

   for(int proc = 0;proc < omp_num_threads;++proc){

      backup_walker[proc].resize(L);

      for(int i = 0;i < L;++i)
         backup_walker[proc][i].resize(d);

   }

   auxvec.resize(omp_num_threads * L);

   for(int i = 0;i < auxvec.size();++i)//for x,y and z components
      auxvec[i].resize(3);

}

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

#ifdef _OPENMP
   omp_num_threads = omp_get_max_threads();
#else
   omp_num_threads = 1;
#endif

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
const TArray<double,2> &Global::gJ(){

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

/**
 * @return the number of omp threads
 */
int Global::gomp_num_threads(){

   return omp_num_threads;

}
