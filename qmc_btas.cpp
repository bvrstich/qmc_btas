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
using namespace Tools;

int main(int argc,char *argv[]){

   cout.precision(15);
   srand(time(NULL));

   int L = 16;
   int d = 2;

   int DW = 4;

   //coupling matrix:
   DArray<2> J(L,L);
   coupling::J1J2_2D(true,0.0,J);

   //hamiltonian
   MPO<complex<double>,Quantum> ham = heisenberg<Quantum>(d,J,0.0);

   //make the trotter
   double dtau = 0.01;
   Trotter trotter(d,J,dtau);

   MPS<complex<double>,Quantum> Psi0;
   Tools::read(Psi0,"/home/bright/bestanden/programmas/dmrg/J1J2/4x4/J2=0.0/Psi0/DT=256.mps");

   cout << inprod(mpsxx::Left,Psi0,ham,Psi0) << endl;
/*
   //read in trial wavefunction
   int nwalkers = 1000;

   AFQMC afqmc(ham,trotter,Psi0,DW,nwalkers,dtau);

   afqmc.Walk(1);
*/
   return 0;

}
