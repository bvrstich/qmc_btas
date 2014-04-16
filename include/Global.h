#ifndef GLOBAL_H
#define GLOBAL_H

#include <iostream>
#include <fstream>
#include <vector>

using std::ostream;
using std::vector;

using namespace btas;

#include "Random.h"

/**
 * @author Brecht Verstichel
 * @date 26-03-2014\n\n
 * Class which stores some globar variables and functions on them for use in program
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

      static const DArray<2> &gJ();

   private:

      //!flag for periodic or non-periodic boundary conditions
      static bool pbc;

      //!physical dimension
      static int d;

      //!max bond dimension of the trial state
      static int D;

      //!coupling matrix of the Hamiltonian
      static DArray<2> J;

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

};

#endif

/* vim: set ts=3 sw=3 expandtab :*/
