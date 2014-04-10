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

int main(int argc,char *argv[]){

   cout.precision(15);
   srand(time(NULL));

   int L = 16;

   int D = 16;
   int d = 2;

   //coupling matrix:
   DArray<2> J(L,L);
   coupling::J1J2_2D(true,0.0,J);

   //hamiltonian
   MPO<complex<double>,Quantum> ham = heisenberg<Quantum>(false,J,0.0);

   //state
   MPS<complex<double>,Quantum> A(L);
   Tools::read(A,"input/J1J2/4x4/J2=0.0/Psi0");

   return 0;

}
