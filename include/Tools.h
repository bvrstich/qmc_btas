#ifndef TOOLS_H
#define TOOLS_H

#include <iostream>
#include <iomanip>
#include <complex>

using namespace btas;
using namespace mpsxx;
using std::complex;

namespace Tools {

   //function which generates random complex numbers uniformly on a square of side 2
   complex<double> rgen_complex(); 

   //function which generates uniform random numbers between [-1:1]
   double rgen_real(); 
      
   //reads in a MPS from file
   void read(MPS<complex<double>,Quantum> &,const char *);

};

#endif
