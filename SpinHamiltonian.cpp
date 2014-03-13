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
void physical(int d,Qshapes<Q> &qp,Dshapes &dp){

   qp.clear();
   dp.clear();

   qp.push_back(Q::zero());
   dp.push_back(d);

}

template<class Q>
void insert_id(QSZArray<4,Q> &O,int row,int col,complex<double> val){

   ZArray<4> Ip(1,1,1,1);
   Ip = val;

   O.insert(shape(row,0,0,col),Ip);
   O.insert(shape(row,1,1,col),Ip);

}

template<class Q>
void insert_Sz(QSZArray<4,Q> &O,int row,int col,complex<double> val){

   ZArray<4> Ip(1,1,1,1);

   Ip = -0.5 * val;
   O.insert(shape(row,0,0,col),Ip);

   Ip = 0.5 * val;
   O.insert(shape(row,1,1,col),Ip);

}

template<class Q>
void insert_Sy(QSZArray<4,Q> &O,int row,int col,complex<double> val){

   ZArray<4> Ip(1,1,1,1);

   Ip = val * complex<double>(0.0,1.0);
   O.insert(shape(row,0,1,col),Ip);

   Ip = val * complex<double>(0.0,-1.0);
   O.insert(shape(row,1,0,col),Ip);

}

template<class Q>
void insert_Sx(QSZArray<4,Q> &O,int row,int col,complex<double> val){

   ZArray<4> Ip(1,1,1,1);

   Ip = val;
   O.insert(shape(row,0,1,col),Ip);

   Ip = val;
   O.insert(shape(row,1,0,col),Ip);

}

/**
 * initialize the MPO<Q> to represent a Heisenberg model with general interaction matrix J_[ij]
 * @param J DArray<2> object containing the couplings
 * @param B magnetic fieldstrength
 */
template<class Q>
MPO<complex<double>,Q> heisenberg(int d,const DArray<2> &J,double B){

   int L = J.shape(0);

   MPO<complex<double>,Q> mpo(L);

   //physical indices
   Qshapes<Q> qp;
   Dshapes dp;

   physical(d,qp,dp);

   Qshapes<Q> qz;
   qz.push_back(Q::zero());

   //initialize the quantumnumbers of the MPO<Q>
   mpo[0].resize(Q::zero(),make_array(qz,qp,-qp,qz));

   //set the input: first site
   /*
   insert_id(mpo[0],0,0,1.0);
   insert_Sx(mpo[0],0,1,1.0);
   insert_Sy(mpo[0],0,2,1.0);
   insert_Sz(mpo[0],0,3,1.0);
   insert_Sz(mpo[0],0,4,B);
   */
/*
   //before middle site
   for(int i = 1;i < L/2;++i){

      qi = -qo;
      qo.clear();

      qo.push_back(Q::zero());//I has spin 0

      for(int row = 0;row < i + 1;++row)
         qo.push_back(Q(2));//S+ has spin 2

      for(int row = 0;row < i + 1;++row)
         qo.push_back(Q(-2));//S- has spin 2

      for(int row = 0;row < i + 1;++row)
         qo.push_back(Q::zero());//Sz has spin 0

      qo.push_back(Q::zero());//B has spin 0

      mpo[i].resize(Q::zero(),make_array(qi,qp,-qp,qo));

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

      qo.push_back(Q::zero());//C has spin 0

      for(int row = L/2 + 1;row < L;++row)
         qo.push_back(Q(2));//S+ has spin 2

      for(int row = L/2 + 1;row < L;++row)
         qo.push_back(Q(-2));//S- has spin 2

      for(int row = L/2 + 1;row < L;++row)
         qo.push_back(Q::zero());//Sz has spin 0

      qo.push_back(Q::zero());//B has spin 0

      mpo[L/2].resize(Q::zero(),make_array(qi,qp,-qp,qo));

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

      qo.push_back(Q::zero());//C has spin 0

      for(int row = i + 1;row < L;++row)
         qo.push_back(Q(2));//S+ has spin 2

      for(int row = i + 1;row < L;++row)
         qo.push_back(Q(-2));//S- has spin 2

      for(int row = i + 1;row < L;++row)
         qo.push_back(Q::zero());//Sz has spin 0

      qo.push_back(Q::zero());//B has spin 0

      mpo[i].resize(Q::zero(),make_array(qi,qp,-qp,qo));

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
   mpo[L - 1].resize(Q::zero(),make_array(qi,qp,-qp,qz));

   insert_id(mpo[L - 1],0,0,1.0);
   insert_Sp(mpo[L - 1],1,0,1.0);
   insert_Sm(mpo[L - 1],2,0,1.0);
   insert_Sz(mpo[L - 1],3,0,1.0);
   insert_Sz(mpo[L - 1],4,0,B);

   //merge everything together
   TVector<Qshapes<Q>,1> qmerge;
   TVector<Dshapes,1> dmerge;

   qmerge[0] = mpo[0].qshape(3);
   dmerge[0] = mpo[0].dshape(3);

   QSTmergeInfo<1> info(qmerge,dmerge);

   QSTArray< complex<double> ,4,Q> tmp;
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

template void physical<Quantum>(int d,Qshapes<Quantum> &,Dshapes &dp);
template MPO<complex<double>,Quantum> heisenberg<Quantum>(int d,const DArray<2> &,double B);
