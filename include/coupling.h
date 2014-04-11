#ifndef COUPLING_H
#define COUPLING_H

#include <iostream>
#include <iomanip>
#include <complex>

using namespace btas;
using std::complex;

namespace coupling {

   void J1J2_1D(bool,double,DArray<2> &);

   void J1J2_2D(bool,double,DArray<2> &);

   void rand(DArray<2> &);

};

#endif
