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

int main(int argc,char *argv[]){

   cout.precision(15);
   srand(time(NULL));

   int L = atoi(argv[1]);
   int d = atoi(argv[2]);
   int D = atoi(argv[3]);

   //coupling matrix:
   DArray<2> J(L,L);
   coupling::J1J2_2D(true,0.0,J);

   char filename[200];
   sprintf(filename,"/home/bright/bestanden/programmas/dmrg/J1J2/4x4/J2=0.0/Psi0/DT=%d.mps",D);

   MPS< complex<double> > mps(filename);

   Walker walker(L,d,1.0,L);
   walker.fill_Random();

   //overlap
   cout << walker.calc_overlap(mps) << endl;

   return 0;

}
