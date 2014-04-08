#ifndef PROPAGATOR_H
#define PROPAGATOR_H

#include <iostream>
#include <fstream>
#include <vector>
#include <complex>

using std::ostream;
using std::vector;
using std::complex;

/**
 * @author Brecht Verstichel
 * @date 21-03-2014\n\n
 * This class Propagator is a class written for the construction and the action of the imaginary time propagator 
 * on a MPS walker.
 */
class Propagator : public vector< QSZArray<2,Quantum> > {

   public:
      
      //constructor
      Propagator(int,int);

      //copy constructor
      Propagator(const Propagator &);

      //destructor
      virtual ~Propagator();

      using vector::operator=;

      using vector::operator[];

      complex<double> gx() const;

      int gk() const;

      int gr() const;

      int gL() const;

      int gd() const;

      void set(complex<double>,int,int);

      void fill(const Trotter &);

   private:

      //!random auxiliary field variable - shift
      complex<double> x;

      //!trotter index
      int k;

      //!type of operator (0=x,1=y or 2=z)
      int r;

      //length of the chain
      int L;

      //!physical dimension
      int d;

};

#endif

/* vim: set ts=3 sw=3 expandtab :*/
