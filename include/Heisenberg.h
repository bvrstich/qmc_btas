#ifndef HEISENBERG_H
#define HEISENBERG_H

#include <iostream>
#include <iomanip>
#include <complex>
#include <vector>

#include <btas/DENSE/DArray.h>
#include <btas/DENSE/ZArray.h>

/**
 * class written for the efficient evaluation of the energy overlap of an MPS with a D=1 walker.
 */
class Heisenberg {

   public:
   
      //Constructor
      static void init(int,const DArray<2> &J);

      static const ZArray<2> &gSx();

      static const ZArray<2> &gSy();

      static const ZArray<2> &gSz();

   private:

      //!coupling matrix
      static DArray<2> J;

      //!bookkeeping, site interactions coming in
      static std::vector< std::vector<int> > in;

      //!bookkeeping, site interactions going out
      static std::vector< std::vector<int> > out;
   
      //The chain length
      static int L;
      
      //!physical dimension
      static int d;

      //!Sx matrix
      static ZArray<2> Sx;

      //!Sy matrix
      static ZArray<2> Sy;

      //!Sz matrix
      static ZArray<2> Sz;

};

#endif
