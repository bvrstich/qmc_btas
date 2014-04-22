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

#include "include.h"

using namespace btas;

int main(int argc,char *argv[]){

   cout.precision(15);
   srand(time(NULL));

   int L = atoi(argv[1]);
   int d = atoi(argv[2]);
   int D = atoi(argv[3]);
   int j2 = atoi(argv[4]);

   Global::init(L,L,j2,d,D,true);

   Heisenberg::init();

   //read in the trial state
   char filename[200];
   sprintf(filename,"input/J1J2/%dx%d/J2=0.%d/Psi0/DT=%d.mps",L,L,j2,D);

   MPS< complex<double> > mps(filename);

   //intialize some storage
   Heisenberg::init_storage(mps);
   Global::init_storage(mps);

   double dtau = 0.01;
   int Nw = 1000;

   AFQMC afqmc(mps,dtau,Nw);
   afqmc.walk(10);

   Heisenberg::clear();

   return 0;

}
