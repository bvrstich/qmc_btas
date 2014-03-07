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

   int L = 10;

   int D = 10;
   int d = 2;

   Qshapes<Quantum> qp;

   int m = -d + 1;

   while(m < d){

      qp.push_back(Quantum(m));

      m += 2;

   }

   MPS<double,Quantum> A = create< double >(L,Quantum::zero(),qp,D,rgen_real);
   normalize(A);

   save_mpx(A,"/tmp/input/A");

   MPS<double,Quantum> B = create< double >(L,Quantum::zero(),qp,D,rgen_real);
   normalize(B);

   save_mpx(B,"/tmp/input/B");

   cout << dot(mpsxx::Left,A,B) << endl;

   cout << compress(A,mpsxx::Left,0) << endl;
   cout << compress(A,mpsxx::Right,0) << endl;

   cout << compress(A,mpsxx::Left,0) << endl;
   cout << compress(A,mpsxx::Right,0) << endl;

   cout << compress(B,mpsxx::Left,0) << endl;
   cout << compress(B,mpsxx::Right,0) << endl;

   cout << compress(B,mpsxx::Left,0) << endl;
   cout << compress(B,mpsxx::Right,0) << endl;

   cout << dot(mpsxx::Left,A,B) << endl;

   return 0;

}
