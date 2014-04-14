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

DArray<2> Heisenberg::J;

int Heisenberg::L;

int Heisenberg::d;

ZArray<2> Heisenberg::Sx;

ZArray<2> Heisenberg::Sy;

ZArray<2> Heisenberg::Sz;

std::vector< std::vector<int> > Heisenberg::in;

std::vector< std::vector<int> > Heisenberg::out;

/**
 * @param d_in physical dimension
 * @param J_in coupling matrix
 */
void Heisenberg::init(int d_in,const DArray<2> &J_in){

   L = J_in.shape(0);
   d = d_in;

   J = J_in;

   //Sx
   Sx.resize(d,d);

   Sx(0,0) = complex<double>(0.0,0.0);
   Sx(0,1) = complex<double>(0.5,0.0);
   Sx(1,0) = complex<double>(0.5,0.0);
   Sx(1,1) = complex<double>(0.0,0.0);

   //Sy
   Sy.resize(d,d);

   Sy(0,0) = complex<double>(0.0,0.0);
   Sy(0,1) = complex<double>(0.0,0.5);
   Sy(1,0) = complex<double>(0.0,-0.5);
   Sy(1,1) = complex<double>(0.0,0.0);

   //Sz
   Sz.resize(d,d);

   Sz(0,0) = complex<double>(-0.5,0.0);
   Sz(0,1) = complex<double>(0.0,0.0);
   Sz(1,0) = complex<double>(0.0,0.0);
   Sz(1,1) = complex<double>(0.5,0.0);

   //out going!
   out.resize(L);

   for(int i = 0;i < L;++i){

      for(int j = i + 1;j < L;++j)
         if( fabs(J(i,j)) > 1.0e-15 )
            out[i].push_back(j);

   }

    //incoming!
   in.resize(L);

   for(int i = 1;i < L;++i){

      for(int j = 0;j < i;++j)
         if( fabs(J(i,j)) > 1.0e-15 )
            in[i].push_back(j);

   }

   for(int i = 0;i < L;++i){

      cout << i << "\t|\t" << in[i] << "\t" << out[i] << endl;

   }

}

/**
 * @return the Sx operator
 */
const ZArray<2> &Heisenberg::gSx() {

   return Sx;

}

/**
 * @return the Sy operator
 */
const ZArray<2> &Heisenberg::gSy() {

   return Sy;

}

/**
 * @return the Sz operator
 */
const ZArray<2> &Heisenberg::gSz() {

   return Sz;

}
