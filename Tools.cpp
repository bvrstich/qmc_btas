#include <iostream>
#include <iomanip>
#include <fstream>

using std::cout;
using std::endl;
using std::ostream;
using std::ifstream;

#include "include.h"

namespace Tools {

   //!A random object: to be used globally
   Random RN;

   //!function which generates random complex numbers uniformly on a square of side 2
   complex<double> rgen_complex(){ 

      return complex<double>(2.0*RN() - 1.0,2.0*RN() - 1.0); 

   }

   //!function which generates uniform random numbers between [-1:1]
   double rgen_real(){ 

     return 2.0*RN() - 1.0;

   }

   /**
    * read in an mps from file "filename" and put it in MPS A
    * @param A MPS will be constructed when calling this function
    * @param filename the inputfile
    */
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
}
