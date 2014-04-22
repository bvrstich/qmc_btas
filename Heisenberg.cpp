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

vector<int> Heisenberg::gemv_list;

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
   blas::gemv(CblasRowMajor, CblasTrans, d, d, one, mps[0].data(), d, walker[0].data(), 1, zero, I[0].data(), 1);

   for(int r = 0;r < 3;++r)
      blas::gemv(CblasRowMajor, CblasTrans, d, d, one, mps[0].data(), d, walker.gVxyz(0,r).data(), 1, zero, ro[0][r].data(), 1);

   C[0] = 0.0;

   //middle sites
   for(int i = 1;i < L - 1;++i){

      int Ldim = mps[i].shape(1);
      int Rdim = mps[i].shape(2);

      //continue from previous site:
      blas::gemv(CblasRowMajor, CblasTrans, d, gemv_list[i], one, mps[i].data(), gemv_list[i], walker[i].data(), 1, zero, loc[i].data(), 1);

      //I
      blas::gemv(CblasRowMajor, CblasTrans, Ldim, Rdim, one, loc[i].data(), Rdim, I[i - 1].data(), 1, zero, I[i].data(), 1);

      //C
      blas::gemv(CblasRowMajor, CblasTrans, Ldim, Rdim, one, loc[i].data(), Rdim, C[i - 1].data(), 1, zero, C[i].data(), 1);

      //transfer previous operators
      for(int j = 0;j < job_cont[i-1].size();++j){

         int input = job_cont[i - 1][j][0];
         int output = job_cont[i - 1][j][1];

         for(int r = 0;r < 3;++r)
            blas::gemv(CblasRowMajor, CblasTrans, Ldim, Rdim, one, loc[i].data(), Rdim, ro[i - 1][input*3 + r].data(), 1, zero, ro[i][output*3 + r].data(), 1);

      }

      //start up new operators and close down couples
      for(int r = 0;r < 3;++r){

         blas::gemv(CblasRowMajor, CblasTrans, d, gemv_list[i], one, mps[i].data(), gemv_list[i], walker.gVxyz(i,r).data(), 1, zero, loc[i].data(), 1);

         //start up
         blas::gemv(CblasRowMajor, CblasTrans, Ldim, Rdim, one, loc[i].data(), Rdim, I[i - 1].data(), 1, zero, ro[i][3*(no[i] - 1) + r].data(), 1);

         //close down
         for(int j = 0;j < job_close[i - 1].size();++j){

            int input = job_close[i - 1][j];

            blas::gemv(CblasRowMajor, CblasTrans, Ldim, Rdim, one, loc[i].data(), Rdim, ro[i - 1][input*3 + r].data(), 1, one, C[i].data(), 1);

         }

      }

   }

   //now get the result
   blas::gemv(CblasRowMajor, CblasTrans, d, d, one, mps[L - 1].data(), d, walker[L - 1].data(), 1, zero, loc[L - 1].data(), 1);

   complex<double> val = blas::dot(d,C[L - 2].data(),1,loc[L - 1].data(),1);

   for(int r = 0;r < 3;++r){

      blas::gemv(CblasRowMajor, CblasTrans, d, d, one, mps[L - 1].data(), d, walker.gVxyz(L - 1,r).data(), 1, zero, loc[L - 1].data(), 1);

      for(int j = 0;j < no[L - 2];++j)
         val += blas::dot(d,ro[L - 2][j*3 + r].data(),1,loc[L - 1].data(),1);

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

   for(int i = 0;i < L;++i)
      loc[i].resize(mps[i].shape(1),mps[i].shape(2));

   //make the dimension list for the gemv routine
   gemv_list.resize(L);

   for(int i = 0;i < L;++i)
      gemv_list[i] = mps[i].shape(1)*mps[i].shape(2);

}
