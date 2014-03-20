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

   if(abs(val) > 1.0e-15){

      ZArray<4> Ip(1,1,1,1);
      Ip = val;

      O.insert(shape(row,0,0,col),Ip);
      O.insert(shape(row,1,1,col),Ip);

   }

}

template<class Q>
void insert_Sz(QSZArray<4,Q> &O,int row,int col,complex<double> val){

   if(abs(val) > 1.0e-15){

      ZArray<4> Ip(1,1,1,1);

      Ip = -0.5 * val;
      O.insert(shape(row,0,0,col),Ip);

      Ip = 0.5 * val;
      O.insert(shape(row,1,1,col),Ip);

   }

}

template<class Q>
void insert_Sy(QSZArray<4,Q> &O,int row,int col,complex<double> val){

   if(abs(val) > 1.0e-15){

      ZArray<4> Ip(1,1,1,1);

      Ip = val * complex<double>(0.0,0.5);
      O.insert(shape(row,0,1,col),Ip);

      Ip = val * complex<double>(0.0,-0.5);
      O.insert(shape(row,1,0,col),Ip);

   }

}

template<class Q>
void insert_Sx(QSZArray<4,Q> &O,int row,int col,complex<double> val){

   if(abs(val) > 1.0e-15){

      ZArray<4> Ip(1,1,1,1);

      Ip = 0.5 * val;
      O.insert(shape(row,0,1,col),Ip);

      Ip = 0.5 * val;
      O.insert(shape(row,1,0,col),Ip);

   }

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

   for(int i = 0;i < d;++i)
      qp.push_back(Q::zero());

   Qshapes<Q> qz;
   qz.push_back(Q::zero());

   Qshapes<Q> qi;
   Qshapes<Q> qo;

   for(int i = 0;i < 5;++i)
      qo.push_back(Q::zero());//I has spin 0

   //initialize the quantumnumbers of the MPO<Q>
   mpo[0].resize(Q::zero(),make_array(qz,qp,-qp,qo));

   //set the input: first site

   insert_id(mpo[0],0,0,1.0);
   insert_Sx(mpo[0],0,1,1.0);
   insert_Sy(mpo[0],0,2,1.0);
   insert_Sz(mpo[0],0,3,1.0);
   insert_Sz(mpo[0],0,4,B);

   //before middle site
   for(int i = 1;i < L/2;++i){

      qi = -qo;
      qo.clear();

      int od = (i + 1)*3 + 2;

      for(int j = 0;j < od;++j)
         qo.push_back(Q::zero());

      mpo[i].resize(Q::zero(),make_array(qi,qp,-qp,qo));

      //first row
      insert_id(mpo[i],0,0,1.0);
      insert_Sx(mpo[i],0,i + 1,1.0);
      insert_Sy(mpo[i],0,1 + (i + 1) + i,1.0);
      insert_Sz(mpo[i],0,1 + 2*(i + 1) + i,1.0);
      insert_Sz(mpo[i],0,1 + 3*(i + 1),B);

      //pass on stuff from previous sites: Sx
      for(int col = 1;col <= i;++col)
         insert_id(mpo[i],col,col,1.0);

      //pass on stuff from previous sites: Sy
      for(int col = 1 + (i + 1);col < 2*(i + 1);++col)
         insert_id(mpo[i],col - 1,col,1.0);

      //pass on stuff from previous sites: Sz
      for(int col = 1 + 2*(i + 1);col < 3*(i + 1);++col)
         insert_id(mpo[i],col - 2,col,1.0);

      //last col: close down everyting
      int row = 1;

      //Sx
      for(int xi = 0;xi < i;++xi){

         insert_Sx(mpo[i],row,od - 1,J(xi,i));
         ++row;

      }

      //Sy
      for(int yi = 0;yi < i;++yi){

         insert_Sy(mpo[i],row,od - 1,J(yi,i));
         ++row;

      }

      //Sz
      for(int zi = 0;zi < i;++zi){

         insert_Sz(mpo[i],row,od - 1,J(zi,i));
         ++row;

      }

      insert_id(mpo[i],row,od - 1,1.0);

   }

   //on the flip site: L/2
   {

      qi = -qo;
      qo.clear();

      int od = 2 + 3*(L/2 - 1);//number of outgoing states

      for(int i = 0;i < od;++i)
         qo.push_back(Q::zero());

      mpo[L/2].resize(Q::zero(),make_array(qi,qp,-qp,qo));

      //first row
      insert_Sz(mpo[L/2],0,0,B);//B

      int col = 1;

      for(int xi = L/2 + 1;xi < L;++xi){

         insert_Sx(mpo[L/2],0,col,J(L/2,xi));
         ++col;

      }

      for(int yi = L/2 + 1;yi < L;++yi){

         insert_Sy(mpo[L/2],0,col, J(L/2,yi));
         ++col;

      }

      for(int zi = L/2 + 1;zi < L;++zi){

         insert_Sz(mpo[L/2],0,col,J(L/2,zi));
         ++col;

      }

      insert_id(mpo[L/2],0,col,1.0);

      int row = 1;

      //now fill the rest: pass on Sx with id's
      for(int xi = 0;xi < L/2;++xi){

         insert_Sx(mpo[L/2],row,0,J(xi,L/2));

         col = 1;

         for(int xj = L/2 + 1;xj < L;++xj){

            insert_id(mpo[L/2],row,col,J(xi,xj));
            ++col;

         }

         ++row;

      }

      //now fill the rest: pass on Sy with id's
      for(int yi = 0;yi < L/2;++yi){

         insert_Sy(mpo[L/2],row,0,J(yi,L/2));

         col = L/2;

         for(int yj = L/2 + 1;yj < L;++yj){

            insert_id(mpo[L/2],row,col,J(yi,yj));
            ++col;

         }

         ++row;

      }

      //now fill the rest: pass on Sz with id's
      for(int zi = 0;zi < L/2;++zi){

         insert_Sz(mpo[L/2],row,0,J(zi,L/2));

         col = L - 1;

         for(int zj = L/2 + 1;zj < L;++zj){

            insert_id(mpo[L/2],row,col,J(zi,zj));
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

      int od = 2 + 3*(L - i - 1);

      for(int j = 0;j < od;++j)
         qo.push_back(Q::zero());//C has spin 0

      mpo[i].resize(Q::zero(),make_array(qi,qp,-qp,qo));

      //first col
      insert_id(mpo[i],0,0,1.0);
      insert_Sx(mpo[i],1,0,1.0);
      insert_Sy(mpo[i],1 + L - i,0,1.0);
      insert_Sz(mpo[i],1 + 2*(L - i),0,1.0);
      insert_Sz(mpo[i],1 + 3*(L - i),0,B);

      //id's!
      int col = 1;
      int row = 2;

      for(int j = i + 1;j < L;++j){

         insert_id(mpo[i],row,col,1.0);

         ++row;
         ++col;

      }

      //Sp
      row = 2 + L - i;

      for(int j = i + 1;j < L;++j){

         insert_id(mpo[i],row,col,1.0);

         ++row;
         ++col;

      }

      //Sz
      row = 2 + 2*(L - i);

      for(int j = i + 1;j < L;++j){

         insert_id(mpo[i],row,col,1.0);

         ++row;
         ++col;

      }

      //last row
      col = 1;

      for(int j = i + 1;j < L;++j){

         insert_Sx(mpo[i],row,col,J(i,j));
         ++col;

      }

      for(int j = i + 1;j < L;++j){

         insert_Sy(mpo[i],row,col,J(i,j));
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
   insert_Sx(mpo[L - 1],1,0,1.0);
   insert_Sy(mpo[L - 1],2,0,1.0);
   insert_Sz(mpo[L - 1],3,0,1.0);
   insert_Sz(mpo[L - 1],4,0,B);
/*
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
