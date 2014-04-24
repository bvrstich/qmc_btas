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

TArray<complex<double>,2> Heisenberg::Sx;

TArray<complex<double>,2> Heisenberg::Sy;

TArray<complex<double>,2> Heisenberg::Sz;

std::vector< std::vector<int> > Heisenberg::in;

std::vector< std::vector<int> > Heisenberg::out;

vector<int> Heisenberg::no;

vector< vector<int> > Heisenberg::to;

vector< vector<int> > Heisenberg::nc;

int Heisenberg::max_ro;

vector< TArray<complex<double>,2> > Heisenberg::I;

vector< TArray<complex<double>,2> > Heisenberg::C;

vector< vector< TArray<complex<double>,2> > > Heisenberg::ro;

vector< vector<int*> > Heisenberg::job_cont;

vector< vector<int> > Heisenberg::job_close;

vector< vector< complex<double> > > Heisenberg::close_coupling;

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

   close_coupling.resize(L - 1);

   for(int i = 0;i < L - 2;++i)
      for(int j = 0;j < job_close[i].size();++j)
         close_coupling[i].push_back(Global::gJ()(to[i][job_close[i][j]],i+1));

   for(int i = 0;i < in[L - 1].size();++i)
      close_coupling[L - 2].push_back(Global::gJ()(in[L - 1][i],L - 1));

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
const TArray<complex<double>,2> &Heisenberg::gSx() {

   return Sx;

}

/**
 * @return the Sy operator
 */
const TArray<complex<double>,2> &Heisenberg::gSy() {

   return Sy;

}

/**
 * @return the Sz operator
 */
const TArray<complex<double>,2> &Heisenberg::gSz() {

   return Sz;

}

/**
 * @return the expectation value of the energy
 */
complex<double> Heisenberg::energy(const MPS< complex<double> > &mps,const Walker &walker){

   int L = Global::gL();
   int d = Global::gd();

#ifdef _OPENMP
   int myID = omp_get_thread_num();
#else
   int myID = 0;
#endif

   complex<double> one(1.0,0.0);
   complex<double> zero(0.0,0.0);

   //first site
   blas::gemv(CblasRowMajor, CblasTrans, d, d, one, mps[0].data(), d, walker[0].data(), 1, zero, I[myID*(L - 1) + 0].data(), 1);

   for(int r = 0;r < 3;++r)
      blas::gemv(CblasRowMajor, CblasTrans, d, d, one, mps[0].data(), d, walker.gVxyz(0,r).data(), 1, zero, ro[myID*(L - 1) + 0][r].data(), 1);

   C[myID*(L - 1) + 0] = 0.0;

   //middle sites
   for(int i = 1;i < L - 1;++i){

      int Ldim = mps[i].shape(1);
      int Rdim = mps[i].shape(2);

      //continue from previous site:
      blas::gemv(CblasRowMajor, CblasTrans, d, Global::gemv_list[i], one, mps[i].data(), Global::gemv_list[i], walker[i].data(), 1, zero, Global::loc[myID*L + i].data(), 1);

      //I
      blas::gemv(CblasRowMajor, CblasTrans, Ldim, Rdim, one, Global::loc[myID*L + i].data(), Rdim, I[myID*(L - 1) + i - 1].data(), 1, zero, I[myID*(L - 1) + i].data(), 1);

      //C
      blas::gemv(CblasRowMajor, CblasTrans, Ldim, Rdim, one, Global::loc[myID*L + i].data(), Rdim, C[myID*(L - 1) + i - 1].data(), 1, zero, C[myID*(L - 1) + i].data(), 1);

      //transfer previous operators
      for(int j = 0;j < job_cont[i-1].size();++j){

         int input = job_cont[i - 1][j][0];
         int output = job_cont[i - 1][j][1];

         for(int r = 0;r < 3;++r)
            blas::gemv(CblasRowMajor, CblasTrans, Ldim, Rdim, one, Global::loc[myID*L + i].data(), Rdim, ro[myID*(L - 1) + i - 1][input*3 + r].data(), 1, zero, ro[myID*(L - 1) + i][output*3 + r].data(), 1);

      }

      //start up new operators and close down couples
      for(int r = 0;r < 3;++r){

         blas::gemv(CblasRowMajor, CblasTrans, d, Global::gemv_list[i], one, mps[i].data(), Global::gemv_list[i], walker.gVxyz(i,r).data(), 1, zero, Global::loc[myID*L + i].data(), 1);

         //start up
         blas::gemv(CblasRowMajor, CblasTrans, Ldim, Rdim, one, Global::loc[myID*L + i].data(), Rdim, I[myID*(L - 1) + i - 1].data(), 1, zero, ro[myID*(L - 1) + i][3*(no[i] - 1) + r].data(), 1);

         //close down
         for(int j = 0;j < job_close[i - 1].size();++j){

            int input = job_close[i - 1][j];

            blas::gemv(CblasRowMajor, CblasTrans, Ldim, Rdim, close_coupling[i - 1][j], Global::loc[myID*L + i].data(), Rdim, ro[myID*(L - 1) + i - 1][input*3 + r].data(), 1, one, C[myID*(L - 1) + i].data(), 1);

         }

      }

   }

   //now get the result
   blas::gemv(CblasRowMajor, CblasTrans, d, d, one, mps[L - 1].data(), d, walker[L - 1].data(), 1, zero, Global::loc[myID*L + L - 1].data(), 1);

   complex<double> val = blas::dot(d,C[myID*(L - 1) + L - 2].data(),1,Global::loc[myID*L + L - 1].data(),1);

   for(int r = 0;r < 3;++r){

      blas::gemv(CblasRowMajor, CblasTrans, d, d, one, mps[L - 1].data(), d, walker.gVxyz(L - 1,r).data(), 1, zero, Global::loc[myID*L + L - 1].data(), 1);

      for(int j = 0;j < no[L - 2];++j)
         val += close_coupling[L - 2][j] * blas::dot(d,ro[myID*(L - 1) + L - 2][j*3 + r].data(),1,Global::loc[myID*L + L - 1].data(),1);

   }

   return val;

}

/**
 * initialize the storages needed for the evaluation of the energy, so that this doesn't have to be done every time during the program
 */
void Heisenberg::init_storage(const MPS< complex<double> > &mps){

   int L = Global::gL();
   int nomp = Global::gomp_num_threads();

   I.resize( nomp * (L - 1) );

   C.resize( nomp * (L - 1) );

   for(int proc = 0;proc < nomp;++proc)
      for(int i = 0;i < L - 1;++i){

         I[proc * (L - 1) + i].resize(1,mps[i].shape(2));
         C[proc * (L - 1) + i].resize(1,mps[i].shape(2));

      }

   ro.resize( nomp * (L - 1) );

   for(int proc = 0;proc < nomp;++proc)
      for(int i = 0;i < L - 1;++i){

         ro[proc * (L - 1) + i].resize(3*no[i]);//xyz 

         for(int j = 0;j < 3*no[i];++j)
            ro[proc * (L - 1) + i][j].resize(1,mps[i].shape(2));

      }

}
