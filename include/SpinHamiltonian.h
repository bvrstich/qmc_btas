#ifndef SPINHAMILTONIAN_H
#define SPINHAMILTONIAN_H

#include <iostream>
#include <iomanip>
#include <complex>

using namespace btas;
using namespace mpsxx;
using std::complex;

//some functions which initialize an MPO to a certian Hamiltonian

template<class Q>
MPO<complex<double>,Q> heisenberg(int d,const DArray<2> &J,double B);

template<class Q>
void physical(int d,Qshapes<Q> &qp,Dshapes &dp);

template<class Q>
void insert_id(QSZArray<4,Q> &,int row,int col,complex<double> val);

template<class Q>
void insert_Sz(QSZArray<4,Q> &,int row,int col,complex<double> val);

template<class Q>
void insert_Sy(QSZArray<4,Q> &,int row,int col,complex<double> val);

template<class Q>
void insert_Sx(QSZArray<4,Q> &,int row,int col,complex<double> val);

#endif
