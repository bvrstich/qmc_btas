#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <complex>

using std::cout;
using std::endl;
using std::ofstream;
using std::ifstream;
using std::complex;

#include "include.h"

using namespace btas;

Random RN;

complex<double> rgen_complex() { return complex<double>(2.0*RN() - 1.0,2.0*RN() - 1.0); }
double rgen_real() { return 2.0*RN() - 1.0; }

int main(int argc,char *argv[]){

   cout.precision(15);
   srand(time(NULL));

   int D = 10;
   int d = 2;

   ZArray<4> M_a(D,D,D,D);
   M_a.generate(rgen_complex);

   ZArray<4> M_b(D,D,D,D);
   M_b.generate(rgen_complex);

   Znormalize(M_a);
   Zorthogonalize(M_a,M_b);

   cout << Zdotc(M_a,M_b) << endl;

   return 0;

}
