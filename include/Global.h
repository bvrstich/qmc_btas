#ifndef GLOBAL_H
#define GLOBAL_H

#include <iostream>
#include <fstream>
#include <vector>
#include <complex>

using std::ostream;
using std::vector;
using std::complex;

template<typename T>
class MPS;

class Walker;

using namespace btas;

#include "Random.h"

/**
 * @author Brecht Verstichel
 * @date 26-03-2014\n\n
 * Class which stores some globar variables, storage, and functions on them for use in program
 */
class Global {

   public:

      //!a Random class object
      static Random RN;

      template<typename T>
         static T rgen();

      static void init(int,int,int,int,int,bool);

      static int gLx();

      static int gLy();

      static int gL();

      static double gJ2();

      static int gj2();

      static int gD();

      static int gd();

      static const TArray<double,2> &gJ();

      static void init_storage(const MPS< complex<double> > &);

      static void clear();

      //!public storage for contractions: for the overlap
      static vector< TArray<complex<double>,1> > LO;
      static vector< TArray<complex<double>,1> > RO;

      //!some storage: local overlap
      static vector< TArray<complex<double>,2> > loc;

      static vector<int> gemv_list;

      static vector< Walker > backup_walker;

      //!intermediate storage for calculation of auxiliary operator expectation values
      static vector< vector< complex<double> > > auxvec;

      static int gomp_num_threads();

   private:

      //!flag for periodic or non-periodic boundary conditions
      static bool pbc;

      //!physical dimension
      static int d;

      //!max bond dimension of the trial state
      static int D;

      //!coupling matrix of the Hamiltonian
      static TArray<double,2> J;

      //!x dimension
      static int Lx;

      //!y dimension
      static int Ly;

      //!length of the chain
      static int L;

      //!J2 coupling parameter
      static double J2;

      //!for output to files
      static int j2;

      //!omp stuff
      static int omp_num_threads;

};

#endif

/* vim: set ts=3 sw=3 expandtab :*/
