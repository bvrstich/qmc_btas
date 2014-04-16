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

   //coupling matrix:
   DArray<2> J(L,L);
   coupling::J1J2_2D(true,0.0,J);

   Heisenberg::init(d,J);

   char filename[200];
   sprintf(filename,"input/J1J2/4x4/J2=0.0/Psi0/DT=%d.mps",D);

   MPS< complex<double> > mps(filename);

   double dtau = 0.01;

   Trotter trotter(d,J,dtau);

   Walker walker(L,d,1.0,trotter.gn_trot());
   sprintf(filename,"input/J1J2/4x4/J2=0.0/PsiW/DT=%d.mps",D);

   walker.read(filename);

   cout << walker.calc_overlap(mps) << endl;
/*
   walker.fill_xyz();

   walker.sVL(trotter,mps);

   for(int i = 0;i < trotter.gn_trot();++i)
      for(int r = 0;r < 3;++r)
         cout << i << "\t" << r << "\t|\t" << walker.gVL(i,r) << endl;
   cout << endl;

   Heisenberg::init_storage(mps);
   cout << Heisenberg::energy(mps,walker) << endl;
*/
   Heisenberg::clear();

   return 0;

}
