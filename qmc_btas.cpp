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

   int L = 20;

   int D = 10;
   int d = 2;

   Qshapes<Quantum> qp;

   int m = -d + 1;

   while(m < d){

      qp.push_back(Quantum(m));

      m += 2;

   }

   MPS<Quantum> A = create(L,Quantum::zero(),qp,D,rgen_real);

   save_mpx(A,"/tmp/input/A");

   normalize(A);

   cout << compress(A,mpsxx::Left,0) << endl;
   cout << compress(A,mpsxx::Right,0) << endl;

   MPS<Quantum> B = create(L,Quantum::zero(),qp,D,rgen_real);

   save_mpx(B,"/tmp/input/B");

   normalize(B);

   cout << compress(B,mpsxx::Left,0) << endl;
   cout << compress(B,mpsxx::Right,0) << endl;

   DArray<2> J(L,L);

   J.generate(rgen_real);

   MPO<Quantum> O = heisenberg<Quantum>(J,0.0);

   cout << inprod(mpsxx::Left,A,O,B) << endl;

   save_mpx(O,"/tmp/input/MPO");

   cout << compress(O,mpsxx::Left,0) << endl;
   cout << compress(O,mpsxx::Right,0) << endl;

   cout << inprod(mpsxx::Left,A,O,B) << endl;

   return 0;

}
