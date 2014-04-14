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

vector<int> Heisenberg::no;

vector< vector<int> > Heisenberg::to(L);

vector< vector<int> > Heisenberg::nc(L);

int Heisenberg::max_ro;

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

   no.resize(L);

   to.resize(L);
   nc.resize(L);

   //site = 0
   no[0] = 1;

   to[0].resize(no[0]);
   nc[0].resize(no[0]);

   to[0][0] = 0;
   nc[0][0] = out[0].size();

   for(int i = 1;i < L;++i){

      no[i] = no[i - 1] + 1;

      to[i] = to[i - 1];
      nc[i] = nc[i - 1];

      to[i].push_back(i);
      nc[i].push_back(out[i].size());

      //now check if anything has reduced
      for(int j = 0;j < in[i].size();++j){

         int type = in[i][j];

         for(int k = 0;k < to[i].size();++k)
            if(type == to[i][k])
               nc[i][k]--;

      }

      int flag = 0;
      int indx = -1;

      while(flag == 0){

         for(int j = 0;j < nc[i].size();++j)
            if(nc[i][j] == 0){

               flag = 1;
               indx = j;
               break;

            }

         if(flag == 1){

            no[i]--;
            nc[i].erase(nc[i].begin() + indx);
            to[i].erase(to[i].begin() + indx);

            flag = 0;

         }
         else
            flag = 2;

      }


   }

   max_ro = 0;

   for(int i = 0;i < L;++i)
      if(no[i] > max_ro)
         max_ro = no[i];

   cout << max_ro << endl;

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

/**
 * @return the expectation value of the energy
 */
complex<double> Heisenberg::energy(const MPS< complex<double> > &mps,const Walker &walker){

   return complex<double> (0.0,0.0);

}
