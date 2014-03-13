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

#include "SpinQuantum.h"
namespace btas { typedef SpinQuantum Quantum; }; // Defined as default quantum number class

#include "include.h"

using namespace btas;
using namespace mpsxx;

Random RN;

complex<double> rgen_complex() { return complex<double>(2.0*RN() - 1.0,2.0*RN() - 1.0); }
double rgen_real() { return 2.0*RN() - 1.0; }

int main(int argc,char *argv[]){

   cout.precision(15);
   srand(time(NULL));

   int L = 4;

   int D = 10;
   int d = 2;
   
   Qshapes<Quantum> qp;
   Dshapes dp;

   physical(d,qp,dp);

   MPS< complex<double> , Quantum > A = create< complex<double> , Quantum >(L,Quantum::zero(),qp,dp,D,rgen_complex);

   DArray<2> J(L,L);
   J.generate(rgen_real);

   MPO< complex<double> , Quantum > O = heisenberg<Quantum>(d,J,0.0);

   cout << O << endl;

   return 0;

}
