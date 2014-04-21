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

ZArray<2> Heisenberg::Sx;

ZArray<2> Heisenberg::Sy;

ZArray<2> Heisenberg::Sz;

std::vector< std::vector<int> > Heisenberg::in;

std::vector< std::vector<int> > Heisenberg::out;

vector<int> Heisenberg::no;

vector< vector<int> > Heisenberg::to;

vector< vector<int> > Heisenberg::nc;

int Heisenberg::max_ro;

vector< ZArray<2> > Heisenberg::I;

vector< ZArray<2> > Heisenberg::C;

vector< vector< ZArray<2> > > Heisenberg::ro;

vector< ZArray<2> > Heisenberg::loc;

vector< vector<int*> > Heisenberg::job_cont;

vector< vector<int> > Heisenberg::job_close;

/**
 * @param d_in physical dimension
 * @param J_in coupling matrix
 */
void Heisenberg::init(){

   int d = Global::gd();
   int L = Global::gL();

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
         if( fabs(Global::gJ()(i,j)) > 1.0e-15 )
            out[i].push_back(j);

   }

   //incoming!
   in.resize(L);

   for(int i = 1;i < L;++i){

      for(int j = 0;j < i;++j)
         if( fabs(Global::gJ()(i,j)) > 1.0e-15 )
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

   //what is the max amount of renormalized operators?
   max_ro = 0;

   for(int i = 0;i < L;++i)
      if(no[i] > max_ro)
         max_ro = no[i];

   //get the instructions
   job_close.resize(L - 2);

   for(int i = 0;i < L - 2;++i){

      for(int j = 0;j < in[i + 1].size();++j){

         for(int k = 0;k < to[i].size();++k)
            if(in[i + 1][j] == to[i][k])
               job_close[i].push_back(k);

      }

   }

   job_cont.resize(L - 2);

   for(int i = 0;i < L - 2;++i){

      for(int j = 0;j < to[i].size();++j)
         for(int k = 0;k < to[i + 1].size();++k){

            if(to[i][j] == to[i + 1][k]){

               int *pair = new int [2];

               pair[0] = j;
               pair[1] = k;

               job_cont[i].push_back(pair);

            }

         }

   }

}

/**
 * remove some memory I allocated dynamically just for the fun of it
 */
void Heisenberg::clear(){

   for(int i = 0;i < job_cont.size();++i)
      for(int j = 0;j < job_cont[i].size();++j)
         delete [] job_cont[i][j];

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

   int L = Global::gL();
   int d = Global::gd();

   complex<double> one(1.0,0.0);
   complex<double> zero(0.0,0.0);

   //first site
   //Gemv(CblasTrans,one,mps[0],walker[0],zero,I[0]);
   blas::gemv(CblasRowMajor, CblasTrans, d, d, one, mps[0].data(), d, walker[0].data(), 1, zero, I[0].data(), 1);

   for(int r = 0;r < 3;++r)
      Gemv(CblasTrans,one,mps[0],walker.gVxyz(0,r),zero,ro[0][r]);

   C[0] = 0.0;

   //middle sites
   for(int i = 1;i < L - 1;++i){

      //continue from previous site:
      Gemv(CblasTrans,one,mps[i],walker[i],zero,loc[i]);

      //I
      Gemm(CblasNoTrans,CblasNoTrans,one,I[i - 1],loc[i],zero,I[i]);

      //C
      Gemm(CblasNoTrans,CblasNoTrans,one,C[i - 1],loc[i],zero,C[i]);

      //transfer previous operators
      for(int j = 0;j < job_cont[i-1].size();++j){

         int input = job_cont[i - 1][j][0];
         int output = job_cont[i - 1][j][1];

         for(int r = 0;r < 3;++r)
            Gemm(CblasNoTrans,CblasNoTrans,one,ro[i - 1][input*3 + r],loc[i],zero,ro[i][output*3 + r]);

      }

      //start up new operators and close down couples
      for(int r = 0;r < 3;++r){

         Gemv(CblasTrans,one,mps[i],walker.gVxyz(i,r),zero,loc[i]);

         //start up
         Gemm(CblasNoTrans,CblasNoTrans,one,I[i - 1],loc[i],zero,ro[i][3*(no[i] - 1) + r]);

         //close down
         for(int j = 0;j < job_close[i - 1].size();++j){

            int input = job_close[i - 1][j];

            Gemm(CblasNoTrans,CblasNoTrans,one,ro[i - 1][input*3 + r],loc[i],one,C[i]);

         }

      }

   }

   //now get the result
   Gemv(CblasTrans,one,mps[L - 1],walker[L - 1],zero,loc[L - 1]);

   complex<double> val = Dot(C[L - 2],loc[L - 1]);

   for(int r = 0;r < 3;++r){

      Gemv(CblasTrans,one,mps[L - 1],walker.gVxyz(L - 1,r),zero,loc[L - 1]);

      for(int j = 0;j < no[L - 2];++j)
         val += Dot(ro[L - 2][j*3 + r],loc[L - 1]);

   }

   return val;

}

/**
 * initialize the storages needed for the evaluation of the energy, so that this doesn't have to be done every time during the program
 */
void Heisenberg::init_storage(const MPS< complex<double> > &mps){

   int L = Global::gL();

   I.resize(L - 1);

   C.resize(L - 1);

   for(int i = 0;i < L - 1;++i){

      I[i].resize(1,mps[i].shape(2));
      C[i].resize(1,mps[i].shape(2));

   }

   ro.resize(L - 1);

   for(int i = 0;i < L - 1;++i){

      ro[i].resize(3*no[i]);//xyz 

      for(int j = 0;j < 3*no[i];++j)
         ro[i][j].resize(1,mps[i].shape(2));

   }

   loc.resize(L);

   for(int i = 0;i < L - 1;++i)
      loc[i].resize(mps[i].shape(1),mps[i].shape(2));

}
