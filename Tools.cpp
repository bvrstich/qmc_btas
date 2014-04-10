#include <iostream>
#include <iomanip>
#include <fstream>

using std::cout;
using std::endl;
using std::ostream;
using std::ifstream;

#include "include.h"

namespace Tools {

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

      MPS<double,Quantum> B(A.size());

      load_mpx(B,filename); 

      cout << B << endl;

   }

}
