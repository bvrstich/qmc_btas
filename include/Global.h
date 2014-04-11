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

};

#endif

/* vim: set ts=3 sw=3 expandtab :*/
