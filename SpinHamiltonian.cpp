#include <iostream>
#include <iomanip>

using std::cout;
using std::endl;
using std::ostream;

#include "include.h"

/**
 * @param d local dimension
 * @param qp Qshapes object containing the local quantumnumbers on output, input is destroyed
 */
template<class Q>
void physical(int d,Qshapes<Q> &qp){

   qp.clear();

   int m = -d + 1;

   while(m < d){

      qp.push_back(Q(m));

      m += 2;

   }

}

/**
 * initialize the MPO<Q> to represent a nearest-neighbour ising Hamiltonian on a lattice of size L and with coupling constant J
 * @param L length of the chain
 * @param d local dimension: i.e. defines the size of the local spins
 * @param J coupling constant
 * @param B magnetic fieldstrength in z direction
 */
template<class Q>
MPO<Q> ising(int L,int d,double J,double B){

   double sz = 0.5 * (d - 1.0);//size of local spin

   MPO<Q> mpo(L);

   //physical indices
   Qshapes<Quantum> qp;

   physical(d,qp);

   Qshapes<Quantum> qz;
   qz.push_back(Quantum::zero());

   //incoming
   Qshapes<Quantum> qi;
   qi.push_back(Quantum::zero());//I has spin 0
   qi.push_back(Quantum::zero());//Sz has spin 0
   qi.push_back(Quantum::zero());//B has spin 0

   //outgoing
   Qshapes<Quantum> qo;
   qo.push_back(Quantum::zero());//I has spin 0
   qo.push_back(Quantum::zero());//Sz has spin 0
   qo.push_back(Quantum::zero());//B has spin 0

   TVector<Qshapes<Quantum>,4> qshape = make_array(qz,qp,-qp,qo);

   //initialize the quantumnumbers of the MPO<Q>
   mpo[0].resize(Quantum::zero(),qshape);

   qshape = make_array(qi,qp,-qp,qo);

   for(int i = 1;i < L-1;++i)
      mpo[i].resize(Quantum::zero(),qshape);

   qshape = make_array(qi,qp,-qp,-qz);

   mpo[L-1].resize(Quantum::zero(),qshape);

   //first the identity
   double mz = -sz;

   for(int m = 0;m < d;++m){

      // set block elements
      DArray<4> I_op(1, 1, 1, 1);//identity
      I_op = 1.0;

      DArray<4> Sz_op(1, 1, 1, 1);//identity
      Sz_op = mz;

      DArray<4> B_op(1, 1, 1, 1);//identity
      B_op = -B*mz;

      DArray<4> J_op(1, 1, 1, 1);//identity
      J_op = J*mz;

      mpo[0].insert(shape(0,m,m,0),I_op);
      mpo[0].insert(shape(0,m,m,1),Sz_op);
      mpo[0].insert(shape(0,m,m,2),B_op);

      for(int i = 1;i < L - 1;++i){

         mpo[i].insert(shape(0,m,m,0),I_op);
         mpo[i].insert(shape(0,m,m,1),Sz_op);
         mpo[i].insert(shape(0,m,m,2),B_op);
         mpo[i].insert(shape(1,m,m,2),J_op);
         mpo[i].insert(shape(2,m,m,2),I_op);

      }

      mpo[L-1].insert(shape(0,m,m,0),B_op);
      mpo[L-1].insert(shape(1,m,m,0),J_op);
      mpo[L-1].insert(shape(2,m,m,0),I_op);

      mz += 1.0;

   }
/*
   //merge everything together
   TVector<Qshapes<Quantum>,1> qmerge;
   TVector<Dshapes,1> dmerge;

   qmerge[0] = mpo[0].qshape(3);
   dmerge[0] = mpo[0].dshape(3);

   QSTmergeInfo<1> info(qmerge,dmerge);

   QSDArray<4> tmp;
   QSTmerge(mpo[0],info,tmp);

   mpo[0] = tmp;

   for(int i = 1;i < L - 1;++i){

      //first merge the row
      qmerge[0] = mpo[i].qshape(0);
      dmerge[0] = mpo[i].dshape(0);

      info.reset(qmerge,dmerge);

      tmp.clear();

      QSTmerge(info,mpo[i],tmp);

      //then merge the column
      qmerge[0] = tmp.qshape(3);
      dmerge[0] = tmp.dshape(3);

      info.reset(qmerge,dmerge);

      mpo[i].clear();

      QSTmerge(tmp,info,mpo[i]);

   }

   //only merge row for i = L - 1
   qmerge[0] = mpo[L - 1].qshape(0);
   dmerge[0] = mpo[L - 1].dshape(0);

   info.reset(qmerge,dmerge);

   tmp.clear();

   QSTmerge(info,mpo[L - 1],tmp);

   mpo[L - 1] = tmp;
*/

   return mpo;

}

/**
 * initialize the MPO<Q> to represent a Sz operator
 * @param L length of the chain
 * @param d local dimension: i.e. defines the size of the local spins
 */
template<class Q>
MPO<Q> Sz(int L,int d){

   double sz = 0.5 * (d - 1.0);//size of local spin

   MPO<Q> mpo(L);

   //physical indices
   Qshapes<Quantum> qp;
   physical(d,qp);

   Qshapes<Quantum> qz;
   qz.push_back(Quantum::zero());//Sz has spin 0

   //incoming
   Qshapes<Quantum> qi;
   qi.push_back(Quantum::zero());//Sz has spin 0
   qi.push_back(Quantum::zero());//I has spin 0

   //outgoing
   Qshapes<Quantum> qo;
   qo.push_back(Quantum::zero());//Sz has spin 0
   qo.push_back(Quantum::zero());//I has spin 0

   TVector<Qshapes<Quantum>,4> qshape = make_array(qz,qp,-qp,qo);

   //initialize the quantumnumbers of the MPO<Q>
   mpo[0].resize(Quantum::zero(),qshape);

   qshape = make_array(qi,qp,-qp,qo);

   for(int i = 1;i < L-1;++i)
      mpo[i].resize(Quantum::zero(),qshape);

   qshape = make_array(qi,qp,-qp,qz);

   mpo[L-1].resize(Quantum::zero(),qshape);

   double mz = -sz;

   for(int m = 0;m < d;++m){

      // set block elements
      DArray<4> I_op(1, 1, 1, 1);//identity
      I_op = 1.0;

      DArray<4> Sz_op(1, 1, 1, 1);//Sz
      Sz_op = mz;

      mpo[0].insert(shape(0,m,m,0),I_op);
      mpo[0].insert(shape(0,m,m,1),Sz_op);

      for(int i = 1;i < L - 1;++i){

         mpo[i].insert(shape(0,m,m,0),I_op);
         mpo[i].insert(shape(0,m,m,1),Sz_op);
         mpo[i].insert(shape(1,m,m,1),I_op);

      }

      mpo[L-1].insert(shape(0,m,m,0),Sz_op);
      mpo[L-1].insert(shape(1,m,m,0),I_op);

      mz += 1.0;

   }

   //merge everything together
   TVector<Qshapes<Quantum>,1> qmerge;
   TVector<Dshapes,1> dmerge;

   qmerge[0] = mpo[0].qshape(3);
   dmerge[0] = mpo[0].dshape(3);

   QSTmergeInfo<1> info(qmerge,dmerge);

   QSDArray<4,Q> tmp;
   QSTmerge(mpo[0],info,tmp);

   mpo[0] = tmp;

   for(int i = 1;i < L - 1;++i){

      //first merge the row
      qmerge[0] = mpo[i].qshape(0);
      dmerge[0] = mpo[i].dshape(0);

      info.reset(qmerge,dmerge);

      tmp.clear();

      QSTmerge(info,mpo[i],tmp);

      //then merge the column
      qmerge[0] = tmp.qshape(3);
      dmerge[0] = tmp.dshape(3);

      info.reset(qmerge,dmerge);

      mpo[i].clear();

      QSTmerge(tmp,info,mpo[i]);

   }

   //only merge row for i = L - 1
   qmerge[0] = mpo[L - 1].qshape(0);
   dmerge[0] = mpo[L - 1].dshape(0);

   info.reset(qmerge,dmerge);

   tmp.clear();

   QSTmerge(info,mpo[L - 1],tmp);

   mpo[L - 1] = tmp;

   return mpo;

}

/**
 * initialize the MPO<Q> to represent a spin raising operator
 * @param L length of the chain
 * @param d local dimension: i.e. defines the size of the local spins
 */
template<class Q>
MPO<Q> raise(int L,int d){

   double sz = 0.5 * (d - 1.0);//size of local spin

   MPO<Q> mpo(L);

   //physical indices
   Qshapes<Quantum> qp;

   physical(d,qp);

   Qshapes<Quantum> qz;
   qz.push_back(Quantum::zero());

   Qshapes<Quantum> qt;
   qt.push_back(Quantum(2));

   //incoming
   Qshapes<Quantum> qi;
   qi.push_back(Quantum::zero());//I has spin 0
   qi.push_back(Quantum(2));//S+ has spin 2

   //outgoing
   Qshapes<Quantum> qo;
   qo.push_back(Quantum::zero());//I has spin 0
   qo.push_back(Quantum(-2));//S+ has spin 2

   TVector<Qshapes<Quantum>,4> qshape = make_array(qz,qp,-qp,qo);

   //initialize the quantumnumbers of the MPO<Q>
   mpo[0].resize(Quantum::zero(),qshape);

   qshape = make_array(qi,qp,-qp,qo);

   for(int i = 1;i < L-1;++i)
      mpo[i].resize(Quantum::zero(),qshape);

   qshape = make_array(qi,qp,-qp,-qt);

   mpo[L-1].resize(Quantum::zero(),qshape);

   //first the identity
   double mz = -sz;

   for(int m = 0;m < d;++m){

      // set block elements
      DArray<4> I_op(1, 1, 1, 1);//identity
      I_op = 1.0;

      mpo[0].insert(shape(0,m,m,0),I_op);

      for(int i = 1;i < L - 1;++i){

         mpo[i].insert(shape(0,m,m,0),I_op);
         mpo[i].insert(shape(1,m,m,1),I_op);

      }

      mpo[L-1].insert(shape(1,m,m,0),I_op);

      mz += 1.0;

   }

   //then the raising operator
   mz = -sz;

   for(int m = 0; m < d - 1; ++m) {

      // set block elements
      DArray<4> Sp(1, 1, 1, 1);
      Sp = std::sqrt( (sz - mz) * (sz + mz + 1.0) );

      // insert blocks
      mpo[0].insert(shape(0,m+1,m,1),Sp);

      for(int i = 1; i < L-1; ++i) 
         mpo[i].insert(shape(0,m+1,m,1),Sp);

      mpo[L-1].insert(shape(0,m+1,m,0),Sp);

      mz  += 1.0;

   }

   return mpo;

}

/**
 * initialize the MPO<Q> to represent a spin lowering operator
 * @param L length of the chain
 * @param d local dimension: i.e. defines the size of the local spins
 */
template<class Q>
MPO<Q> lower(int L,int d){

   double sz = 0.5 * (d - 1.0);//size of local spin

   MPO<Q> mpo(L);

   //physical indices
   Qshapes<Quantum> qp;

   physical(d,qp);

   Qshapes<Quantum> qz;
   qz.push_back(Quantum::zero());

   Qshapes<Quantum> qt;
   qt.push_back(Quantum(-2));

   //incoming
   Qshapes<Quantum> qi;
   qi.push_back(Quantum::zero());//I has spin 0
   qi.push_back(Quantum(-2));//S+ has spin 2

   //outgoing
   Qshapes<Quantum> qo;
   qo.push_back(Quantum::zero());//I has spin 0
   qo.push_back(Quantum(2));//S+ has spin 2

   TVector<Qshapes<Quantum>,4> qshape = make_array(qz,qp,-qp,qo);

   //initialize the quantumnumbers of the MPO<Q>
   mpo[0].resize(Quantum::zero(),qshape);

   qshape = make_array(qi,qp,-qp,qo);

   for(int i = 1;i < L-1;++i)
      mpo[i].resize(Quantum::zero(),qshape);

   qshape = make_array(qi,qp,-qp,-qt);

   mpo[L-1].resize(Quantum::zero(),qshape);

   //first the identity
   double mz = -sz;

   for(int m = 0;m < d;++m){

      // set block elements
      DArray<4> I_op(1, 1, 1, 1);//identity
      I_op = 1.0;

      mpo[0].insert(shape(0,m,m,0),I_op);

      for(int i = 1;i < L - 1;++i){

         mpo[i].insert(shape(0,m,m,0),I_op);
         mpo[i].insert(shape(1,m,m,1),I_op);

      }

      mpo[L-1].insert(shape(1,m,m,0),I_op);

      mz += 1.0;

   }

   //then the lowering operator
   mz = sz;

   for(int m = d - 1;m > 0;--m) {

      // set block elements
      DArray<4> Sm(1, 1, 1, 1);
      Sm = std::sqrt( (sz + mz) * (sz - mz + 1.0) );

      // insert blocks
      mpo[0].insert(shape(0,m-1,m,1),Sm);

      for(int i = 1; i < L-1; ++i) 
         mpo[i].insert(shape(0,m-1,m,1),Sm);

      mpo[L-1].insert(shape(0,m-1,m,0),Sm);

      mz  -= 1.0;

   }

   return mpo;

}

/**
 * initialize the MPO<Q> to represent a nearest-neighbour ising Hamiltonian on a lattice of size L and with coupling constant J
 * @param L length of the chain
 * @param d local dimension: i.e. defines the size of the local spins
 * @param J coupling constant between S+S-
 * @param B magnetic fieldstrength in z direction
 */
template<class Q>
MPO<Q> XY(int L,int d,double J,double B){

   double sz = 0.5 * (d - 1.0);//size of local spin

   MPO<Q> mpo(L);

   //physical indices
   Qshapes<Quantum> qp;

   physical(d,qp);

   Qshapes<Quantum> qz;
   qz.push_back(Quantum::zero());

   //incoming
   Qshapes<Quantum> qi;
   qi.push_back(Quantum::zero());//I has spin 0
   qi.push_back(Quantum(2));//S- has spin -2
   qi.push_back(Quantum(-2));//S+ has spin +2
   qi.push_back(Quantum::zero());//B has spin 0

   //outgoing
   Qshapes<Quantum> qo;
   qo.push_back(Quantum::zero());//I has spin 0
   qo.push_back(Quantum(-2));//S+ has spin 2
   qo.push_back(Quantum(2));//S- has spin 2
   qo.push_back(Quantum::zero());//B has spin 0

   TVector<Qshapes<Quantum>,4> qshape = make_array(qz,qp,-qp,qo);

   //initialize the quantumnumbers of the MPO<Q>
   mpo[0].resize(Quantum::zero(),qshape);

   qshape = make_array(qi,qp,-qp,qo);

   for(int i = 1;i < L-1;++i)
      mpo[i].resize(Quantum::zero(),qshape);

   qshape = make_array(qi,qp,-qp,qz);

   mpo[L-1].resize(Quantum::zero(),qshape);

   //first the identity and local term
   double mz = -sz;

   for(int m = 0;m < d;++m){

      // set block elements
      DArray<4> I_op(1, 1, 1, 1);//identity
      I_op = 1.0;

      DArray<4> B_op(1, 1, 1, 1);//identity
      B_op = -B*mz;

      mpo[0].insert(shape(0,m,m,0),I_op);
      mpo[0].insert(shape(0,m,m,3),B_op);

      for(int i = 1;i < L - 1;++i){

         mpo[i].insert(shape(0,m,m,0),I_op);
         mpo[i].insert(shape(0,m,m,3),B_op);
         mpo[i].insert(shape(3,m,m,3),I_op);

      }

      mpo[L-1].insert(shape(0,m,m,0),B_op);
      mpo[L-1].insert(shape(3,m,m,0),I_op);

      mz += 1.0;

   }

   //then the raising operator
   mz = -sz;

   for(int m = 0; m < d - 1; ++m) {

      // set block elements
      DArray<4> Sp(1, 1, 1, 1);
      Sp = std::sqrt( J*(sz - mz) * (sz + mz + 1.0) );

      // insert blocks
      mpo[0].insert(shape(0,m+1,m,1),Sp);

      for(int i = 1;i < L - 1;++i){

         mpo[i].insert(shape(0,m+1,m,1),Sp);
         mpo[i].insert(shape(2,m+1,m,3),Sp);

      }

      mpo[L-1].insert(shape(2,m+1,m,0),Sp);

      mz  += 1.0;

   }

   //then the lowering operator
   mz = sz;

   for(int m = d-1; m > 0;--m) {

      // set block elements
      DArray<4> Sm(1, 1, 1, 1);
      Sm = std::sqrt( J*(sz + mz) * (sz - mz + 1.0) );

      // insert blocks
      mpo[0].insert(shape(0,m-1,m,2),Sm);

      for(int i = 1;i < L - 1;++i){

         mpo[i].insert(shape(0,m-1,m,2),Sm);
         mpo[i].insert(shape(1,m-1,m,3),Sm);

      }

      mpo[L-1].insert(shape(1,m-1,m,0),Sm);

      mz -= 1.0;

   }


   //merge everything together
   TVector<Qshapes<Quantum>,1> qmerge;
   TVector<Dshapes,1> dmerge;

   qmerge[0] = mpo[0].qshape(3);
   dmerge[0] = mpo[0].dshape(3);

   QSTmergeInfo<1> info(qmerge,dmerge);

   QSDArray<4,Q> tmp;
   QSTmerge(mpo[0],info,tmp);

   mpo[0] = tmp;

   for(int i = 1;i < L - 1;++i){

      //first merge the row
      qmerge[0] = mpo[i].qshape(0);
      dmerge[0] = mpo[i].dshape(0);

      info.reset(qmerge,dmerge);

      tmp.clear();

      QSTmerge(info,mpo[i],tmp);

      //then merge the column
      qmerge[0] = tmp.qshape(3);
      dmerge[0] = tmp.dshape(3);

      info.reset(qmerge,dmerge);

      mpo[i].clear();

      QSTmerge(tmp,info,mpo[i]);

   }

   //only merge row for i = L - 1
   qmerge[0] = mpo[L - 1].qshape(0);
   dmerge[0] = mpo[L - 1].dshape(0);

   info.reset(qmerge,dmerge);

   tmp.clear();

   QSTmerge(info,mpo[L - 1],tmp);

   mpo[L - 1] = tmp;

   return mpo;

}

template<class Q>
void insert_id(QSDArray<4,Q> &O,int row,int col,double val){

   DArray<4> Ip(1,1,1,1);
   Ip = val;

   O.insert(shape(row,0,0,col),Ip);
   O.insert(shape(row,1,1,col),Ip);

}

template<class Q>
void insert_Sz(QSDArray<4,Q> &O,int row,int col,double val){

   DArray<4> Ip(1,1,1,1);

   Ip = -0.5 * val;
   O.insert(shape(row,0,0,col),Ip);

   Ip = 0.5 * val;
   O.insert(shape(row,1,1,col),Ip);

}

template<class Q>
void insert_Sp(QSDArray<4,Q> &O,int row,int col,double val){

   DArray<4> Ip(1,1,1,1);

   Ip = val;
   O.insert(shape(row,1,0,col),Ip);

}

template<class Q>
void insert_Sm(QSDArray<4,Q> &O,int row,int col,double val){

   DArray<4> Ip(1,1,1,1);

   Ip = val;
   O.insert(shape(row,0,1,col),Ip);

}

/**
 * initialize the MPO<Q> to represent a nearest-neighbour anisotropic Heisenberg Hamiltonian on a lattice of size L and with coupling constant Jz for the z spins and Jxy for XY
 * inside a magnetic field B which defines the z direction. nearest neighbour interaction is repulsive, i.e. Jz and Jxy > 0.0!
 * @param L length of the chain
 * @param d local dimension: i.e. defines the size of the local spins
 * @param Jz coupling constant > 0 
 * @param Jxy coupling constant > 0
 * @param B magnetic fieldstrength
 */
template<class Q>
MPO<Q> heisenberg(int L,int d,double J,double B){

   double sz = 0.5 * (d - 1.0);//size of local spin

   MPO<Q> mpo(L);

   //physical indices
   Qshapes<Quantum> qp;

   physical(d,qp);

   Qshapes<Quantum> qz;
   qz.push_back(Quantum::zero());

   //incoming
   Qshapes<Quantum> qi;
   qi.push_back(Quantum::zero());//I has spin 0
   qi.push_back(Quantum(-2));//S- has spin -2
   qi.push_back(Quantum(2));//S+ has spin +2
   qi.push_back(Quantum::zero());//Sz has spin 0
   qi.push_back(Quantum::zero());//B has spin 0

   //outgoing
   Qshapes<Quantum> qo;
   qo.push_back(Quantum::zero());//I has spin 0
   qo.push_back(Quantum(2));//S+ has spin 2
   qo.push_back(Quantum(-2));//S- has spin 2
   qo.push_back(Quantum::zero());//Sz has spin 0
   qo.push_back(Quantum::zero());//B has spin 0

   //initialize the quantumnumbers of the MPO<Q>
   mpo[0].resize(Quantum::zero(),make_array(qz,qp,-qp,qo));

   for(int i = 1;i < L-1;++i)
      mpo[i].resize(Quantum::zero(),make_array(qi,qp,-qp,qo));

   mpo[L-1].resize(Quantum::zero(),make_array(qi,qp,-qp,-qz));

   //set the input: first site
   insert_id(mpo[0],0,0,1.0);
   insert_Sm(mpo[0],0,1,sqrt(0.5 * J));
   insert_Sp(mpo[0],0,2,sqrt(0.5 * J));
   insert_Sz(mpo[0],0,3,sqrt(J));
   insert_Sz(mpo[0],0,4,B);

   //middle sites
   for(int i = 1;i < L - 1;++i){

      //first row
      insert_id(mpo[i],0,0,1.0);
      insert_Sm(mpo[i],0,1,sqrt(0.5 * J));
      insert_Sp(mpo[i],0,2,sqrt(0.5 * J));
      insert_Sz(mpo[i],0,3,sqrt(J));
      insert_Sz(mpo[i],0,4,B);

      //last col
      insert_Sp(mpo[i],1,4,sqrt(0.5 * J));
      insert_Sm(mpo[i],2,4,sqrt(0.5 * J));
      insert_Sz(mpo[i],3,4,sqrt(J));
      insert_id(mpo[i],4,4,1.0);

   }

   //last site
   insert_Sz(mpo[L-1],0,0,B);
   insert_Sp(mpo[L-1],1,0,sqrt(0.5 * J));
   insert_Sm(mpo[L-1],2,0,sqrt(0.5 * J));
   insert_Sz(mpo[L-1],3,0,sqrt(J));
   insert_id(mpo[L-1],4,0,1.0);

   //merge everything together
   TVector<Qshapes<Quantum>,1> qmerge;
   TVector<Dshapes,1> dmerge;

   qmerge[0] = mpo[0].qshape(3);
   dmerge[0] = mpo[0].dshape(3);

   QSTmergeInfo<1> info(qmerge,dmerge);

   QSDArray<4,Q> tmp;
   QSTmerge(mpo[0],info,tmp);

   mpo[0] = tmp;

   for(int i = 1;i < L - 1;++i){

      //first merge the row
      qmerge[0] = mpo[i].qshape(0);
      dmerge[0] = mpo[i].dshape(0);

      info.reset(qmerge,dmerge);

      tmp.clear();

      QSTmerge(info,mpo[i],tmp);

      //then merge the column
      qmerge[0] = tmp.qshape(3);
      dmerge[0] = tmp.dshape(3);

      info.reset(qmerge,dmerge);

      mpo[i].clear();

      QSTmerge(tmp,info,mpo[i]);

   }

   //only merge row for i = L - 1
   qmerge[0] = mpo[L - 1].qshape(0);
   dmerge[0] = mpo[L - 1].dshape(0);

   info.reset(qmerge,dmerge);

   tmp.clear();

   QSTmerge(info,mpo[L - 1],tmp);

   mpo[L - 1] = tmp;

   return mpo;

}

/**
 * initialize the MPO<Q> to represent a Heisenberg model with general interaction matrix J_[ij]
 * @param J DArray<2> object containing the couplings
 * @param B magnetic fieldstrength
 */
template<class Q>
MPO<Q> heisenberg(const DArray<2> &J,double B){

   int L = J.shape(0);

   MPO<Q> mpo(L);

   //physical indices
   Qshapes<Quantum> qp;
   physical(2,qp);

   Qshapes<Quantum> qz;
   qz.push_back(Quantum::zero());

   //incoming
   Qshapes<Quantum> qi;
   Qshapes<Quantum> qo;

   qo.push_back(Quantum::zero());//I has spin 0
   qo.push_back(Quantum(2));//S+ has spin 2
   qo.push_back(Quantum(-2));//S- has spin 2
   qo.push_back(Quantum::zero());//Sz has spin 0
   qo.push_back(Quantum::zero());//B has spin 0

   //initialize the quantumnumbers of the MPO<Q>
   mpo[0].resize(Quantum::zero(),make_array(qz,qp,-qp,qo));

   //set the input: first site
   insert_id(mpo[0],0,0,1.0);
   insert_Sm(mpo[0],0,1,1.0);
   insert_Sp(mpo[0],0,2,1.0);
   insert_Sz(mpo[0],0,3,1.0);
   insert_Sz(mpo[0],0,4,B);

   //before middle site
   for(int i = 1;i < L/2;++i){

      qi = -qo;
      qo.clear();

      qo.push_back(Quantum::zero());//I has spin 0

      for(int row = 0;row < i + 1;++row)
         qo.push_back(Quantum(2));//S+ has spin 2

      for(int row = 0;row < i + 1;++row)
         qo.push_back(Quantum(-2));//S- has spin 2

      for(int row = 0;row < i + 1;++row)
         qo.push_back(Quantum::zero());//Sz has spin 0

      qo.push_back(Quantum::zero());//B has spin 0

      mpo[i].resize(Quantum::zero(),make_array(qi,qp,-qp,qo));

      //first row
      insert_id(mpo[i],0,0,1.0);
      insert_Sm(mpo[i],0,i + 1,1.0);
      insert_Sp(mpo[i],0,2*i + 2,1.0);
      insert_Sz(mpo[i],0,3*i + 3,1.0);

      //pass on stuff from previous sites: S-
      for(int col = 1;col < i + 1;++col)
         insert_id(mpo[i],col,col,1.0);

      //pass on stuff from previous sites: S+
      for(int col = i + 2;col < 2*i + 2;++col)
         insert_id(mpo[i],col - 1,col,1.0);

      //pass on stuff from previous sites: Sz
      for(int col = 2*i + 3;col < 3*i + 3;++col)
         insert_id(mpo[i],col - 2,col,1.0);

      //last col: close down everyting
      insert_Sz(mpo[i],0,3*i + 4,B);

      int row = 1;

      //S-
      for(int min_ri = 0;min_ri < i;++min_ri){

         insert_Sp(mpo[i],row,3*i + 4,0.5 * J(min_ri,i));
         ++row;

      }

      //S+
      for(int pl_ri = 0;pl_ri < i;++pl_ri){

         insert_Sm(mpo[i],row,3*i + 4,0.5 * J(pl_ri,i));
         ++row;

      }

      //Sz
      for(int z_ri = 0;z_ri < i;++z_ri){

         insert_Sz(mpo[i],row,3*i + 4,J(z_ri,i));
         ++row;

      }

      insert_id(mpo[i],row,3*i + 4,1.0);

   }

   //on the flip site: L/2
   {

      qi = -qo;
      qo.clear();

      qo.push_back(Quantum::zero());//C has spin 0

      for(int row = L/2 + 1;row < L;++row)
         qo.push_back(Quantum(2));//S+ has spin 2

      for(int row = L/2 + 1;row < L;++row)
         qo.push_back(Quantum(-2));//S- has spin 2

      for(int row = L/2 + 1;row < L;++row)
         qo.push_back(Quantum::zero());//Sz has spin 0

      qo.push_back(Quantum::zero());//B has spin 0

      mpo[L/2].resize(Quantum::zero(),make_array(qi,qp,-qp,qo));

      //first row
      insert_Sz(mpo[L/2],0,0,B);//B

      int col = 1;

      for(int i = L/2 + 1;i < L;++i){

         insert_Sm(mpo[L/2],0,col,0.5 * J(L/2,i));
         ++col;

      }

      int col_p = col;

      for(int i = L/2 + 1;i < L;++i){

         insert_Sp(mpo[L/2],0,col,0.5 * J(L/2,i));
         ++col;

      }

      int col_z = col;

      for(int i = L/2 + 1;i < L;++i){

         insert_Sz(mpo[L/2],0,col,J(L/2,i));
         ++col;

      }

      insert_id(mpo[L/2],0,col,1.0);

      int row = 1;

      //now fill the rest: pass on S- with id's
      for(int i = 0;i < L/2;++i){

         insert_Sp(mpo[L/2],row,0,0.5 * J(i,L/2));

         col = 1;

         for(int j = L/2 + 1;j < L;++j){

            insert_id(mpo[L/2],row,col,0.5 * J(i,j));
            ++col;

         }

         ++row;

      }

      //now fill the rest: pass on S+ with id's
      for(int i = 0;i < L/2;++i){

         insert_Sm(mpo[L/2],row,0,0.5 * J(i,L/2));

         col = col_p;

         for(int j = L/2 + 1;j < L;++j){

            insert_id(mpo[L/2],row,col,0.5 * J(i,j));
            ++col;

         }

         ++row;

      }

      //now fill the rest: pass on Sz with id's
      for(int i = 0;i < L/2;++i){

         insert_Sz(mpo[L/2],row,0,J(i,L/2));

         col = col_z;

         for(int j = L/2 + 1;j < L;++j){

            insert_id(mpo[L/2],row,col,J(i,j));
            ++col;

         }

         ++row;

      }

      insert_id(mpo[L/2],row,0,1.0);//1

   }

   //after the flip site
   for(int i = L/2 + 1;i < L - 1;++i){

      qi = -qo;
      qo.clear();

      qo.push_back(Quantum::zero());//C has spin 0

      for(int row = i + 1;row < L;++row)
         qo.push_back(Quantum(2));//S+ has spin 2

      for(int row = i + 1;row < L;++row)
         qo.push_back(Quantum(-2));//S- has spin 2

      for(int row = i + 1;row < L;++row)
         qo.push_back(Quantum::zero());//Sz has spin 0

      qo.push_back(Quantum::zero());//B has spin 0

      mpo[i].resize(Quantum::zero(),make_array(qi,qp,-qp,qo));

      //first col
      insert_id(mpo[i],0,0,1.0);

      insert_Sp(mpo[i],1,0,1.0);

      int row = 1;

      for(int j = i;j < L;++j)
         ++row;

      int row_p = row;

      insert_Sm(mpo[i],row,0,1.0);

      for(int j = i;j < L;++j)
         ++row;

      int row_z = row;

      insert_Sz(mpo[i],row,0,1.0);

      for(int j = i;j < L;++j)
         ++row;

      insert_Sz(mpo[i],row,0,B);

      //id's!
      int col = 1;
      row = 2;

      for(int j = i + 1;j < L;++j){

         insert_id(mpo[i],row,col,1.0);

         ++row;
         ++col;

      }

      //Sp
      row = row_p + 1;

      for(int j = i + 1;j < L;++j){

         insert_id(mpo[i],row,col,1.0);

         ++row;
         ++col;

      }

      //Sz
      row = row_z + 1;

      for(int j = i + 1;j < L;++j){

         insert_id(mpo[i],row,col,1.0);

         ++row;
         ++col;

      }

      //last row
      col = 1;

      for(int j = i + 1;j < L;++j){

         insert_Sm(mpo[i],row,col,0.5*J(i,j));
         ++col;

      }

      for(int j = i + 1;j < L;++j){

         insert_Sp(mpo[i],row,col,0.5*J(i,j));
         ++col;

      }

      for(int j = i + 1;j < L;++j){

         insert_Sz(mpo[i],row,col,J(i,j));
         ++col;

      }

      insert_id(mpo[i],row,col,1.0);

   }

   //last site
   qi = -qo;
   mpo[L - 1].resize(Quantum::zero(),make_array(qi,qp,-qp,qz));

   insert_id(mpo[L - 1],0,0,1.0);
   insert_Sp(mpo[L - 1],1,0,1.0);
   insert_Sm(mpo[L - 1],2,0,1.0);
   insert_Sz(mpo[L - 1],3,0,1.0);
   insert_Sz(mpo[L - 1],4,0,B);

   //merge everything together
   TVector<Qshapes<Quantum>,1> qmerge;
   TVector<Dshapes,1> dmerge;

   qmerge[0] = mpo[0].qshape(3);
   dmerge[0] = mpo[0].dshape(3);

   QSTmergeInfo<1> info(qmerge,dmerge);

   QSDArray<4,Q> tmp;
   QSTmerge(mpo[0],info,tmp);

   mpo[0] = tmp;

   for(int i = 1;i < L - 1;++i){

      //first merge the row
      qmerge[0] = mpo[i].qshape(0);
      dmerge[0] = mpo[i].dshape(0);

      info.reset(qmerge,dmerge);

      tmp.clear();

      QSTmerge(info,mpo[i],tmp);

      //then merge the column
      qmerge[0] = tmp.qshape(3);
      dmerge[0] = tmp.dshape(3);

      info.reset(qmerge,dmerge);

      mpo[i].clear();

      QSTmerge(tmp,info,mpo[i]);

   }

   //only merge row for i = L - 1
   qmerge[0] = mpo[L - 1].qshape(0);
   dmerge[0] = mpo[L - 1].dshape(0);

   info.reset(qmerge,dmerge);

   tmp.clear();

   QSTmerge(info,mpo[L - 1],tmp);

   mpo[L - 1] = tmp;

   return mpo;

}

template void physical<Quantum>(int d,Qshapes<Quantum> &);
template MPO<Quantum> ising<Quantum>(int L,int d,double J,double B);
template MPO<Quantum> Sz<Quantum>(int L,int d);
template MPO<Quantum> raise<Quantum>(int L,int d);
template MPO<Quantum> lower<Quantum>(int L,int d);
template MPO<Quantum> XY<Quantum>(int L,int d,double J,double B);
template MPO<Quantum> heisenberg<Quantum>(int L,int d,double J,double B);
template MPO<Quantum> heisenberg<Quantum>(const DArray<2> &,double B);
