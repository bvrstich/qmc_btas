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
      static void init();

      static void clear();

      static const ZArray<2> &gSx();

      static const ZArray<2> &gSy();

      static const ZArray<2> &gSz();

      static complex<double> energy(const MPS< complex<double> > &,const Walker &);

      static void init_storage(const MPS< complex<double> > &);

   private:

      //!bookkeeping, site interactions coming in
      static std::vector< std::vector<int> > in;

      //!bookkeeping, site interactions going out
      static std::vector< std::vector<int> > out;

      //!Sx matrix
      static ZArray<2> Sx;

      //!Sy matrix
      static ZArray<2> Sy;

      //!Sz matrix
      static ZArray<2> Sz;

      //!number of operators required on site
      static vector<int> no;

      //!type of the different operators required on site
      static vector< vector<int> > to;

      //!number of connecting sites to the right of this site with this operator
      static vector< vector<int> > nc;

      //!the max number of operators
      static int max_ro;

      //!some storage for the energy evaluations: identity
      static vector< ZArray<2> > I;

      //!some storage for the energy evaluations: closed terms
      static vector< ZArray<2> > C;

      //!some storage for the energy evaluations: connecting terms
      static vector< vector< ZArray<2> > > ro;

      //!some storage for the energy evaluations: local overlap
      static vector< ZArray<2> > loc;

      //!operations to be performed at every site: input-output
      static vector< vector<int*> > job_cont;

      //!operations to be performed at every site: input
      static vector< vector<int> > job_close;

      static vector<int> gemv_list;

};

#endif
