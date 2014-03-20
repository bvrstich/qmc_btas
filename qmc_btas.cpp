#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <complex>
#include <vector>

using std::cout;
using std::endl;
using std::ofstream;
using std::ifstream;
using std::complex;
using std::vector;

#include "SpinQuantum.h"
namespace btas { typedef SpinQuantum Quantum; }; // Defined as default quantum number class

#include "include.h"

using namespace btas;
using namespace mpsxx;

Random RN;

complex<double> rgen_complex() { return complex<double>(2.0*RN() - 1.0,2.0*RN() - 1.0); }
double rgen_real() { return 2.0*RN() - 1.0; }
void read(MPS<complex<double>,Quantum> &,const char *);

int main(int argc,char *argv[]){

   cout.precision(15);
   srand(time(NULL));

   int L = 16;
   int d = 2;
   
   DArray<2> J(L,L);
   coupling::J1J2_2D(true,0.0,J);

   double dtau = 0.01;

   Trotter trotter(J,dtau);

   return 0;

}

void read(MPS<complex<double>,Quantum> &A,const char *filename){

   ifstream in(filename);

   int L;
   int D,d;

   //get the length
   in >> L >> D >> d;

   A.resize(L);

   //get the virtual dimensions
   vector<int> vdim(L + 1);

   for(int i = 0;i < L + 1;++i)
      in >> i >> vdim[i];

   Qshapes<Quantum> qp;
   qp.push_back(Quantum::zero());
   qp.push_back(Quantum::zero());

   Dshapes dp;
   dp.push_back(1);
   dp.push_back(1);

   Qshapes<Quantum> qz;
   qz.push_back(Quantum::zero());

   int tmp,tmp2;

   //allocate the matrices
   for(int i = 0;i < L;++i){

      in >> tmp >> tmp2;

      Dshapes id;
      id.push_back(vdim[i]);

      Dshapes od;
      od.push_back(vdim[i + 1]);

      A[i].resize(Quantum::zero(),make_array(qz,qp,qz),make_array(id,dp,od));

      for(int s = 0;s < 2;++s){

         STArray<complex<double>,3>::iterator it = A[i].find(make_array(0,s,0));

         int I;
         int index;
         double value;

            for(int dk = 0;dk < vdim[i + 1];++dk)
               for(int dj = 0;dj < vdim[i];++dj){

                  in >> i >> index >> value;

                  (*it->second)(dj,0,dk) = complex<double>(value,0.0);

               }

      }

   }

}
