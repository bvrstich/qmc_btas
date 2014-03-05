#ifndef SPINHAMILTONIAN_H
#define SPINHAMILTONIAN_H

#include <iostream>
#include <iomanip>

using namespace btas;
using namespace mpsxx;

//some functions which initialize an MPO to a certian Hamiltonian
template<class Q>
MPO<Q> ising(int,int,double,double);

template<class Q>
MPO<Q> XY(int,int,double,double);

template<class Q>
MPO<Q> heisenberg(int,int,double,double);

template<class Q>
MPO<Q> heisenberg(const DArray<2> &J,double B);

template<class Q>
MPO<Q> raise(int,int);

template<class Q>
MPO<Q> lower(int,int);

template<class Q>
MPO<Q> Sz(int,int);

template<class Q>
void physical(int d,Qshapes<Q> &qp);

template<class Q>
void insert_id(QSDArray<4,Q> &,int row,int col,double val);

template<class Q>
void insert_Sz(QSDArray<4,Q> &,int row,int col,double val);

template<class Q>
void insert_Sp(QSDArray<4,Q> &,int row,int col,double val);

template<class Q>
void insert_Sm(QSDArray<4,Q> &,int row,int col,double val);

#endif
