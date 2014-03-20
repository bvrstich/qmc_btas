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
