#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <complex>
#include <vector>

#include "include.h"

using std::cout;
using std::endl;
using std::complex;
using std::vector;
using std::ofstream;
using std::ifstream;
using std::ios;

/**
 * constructor of the AFQMC object, takes input parameters that define the QMC walk.
 * @param trotter trotter terms for the interaction
 * @param Psi0_in input trialwavefunction
 * @param Nwalkers number of Walker states
 * @param dtau time step of each evolution
 */
AFQMC::AFQMC(const MPS< complex<double> > &Psi0_in,double dtau_in,int Nw_in){

   this->dtau = dtau_in;
   this->Psi0 = Psi0_in;

   this->trotter = new Trotter(dtau);

   n_trot = trotter->gn_trot();

   this->Nw = Nw_in;

   SetupWalkers();

   P = new Propagator();

}

AFQMC::~AFQMC(){

   for(int i = 0;i < walker.size();++i)
      delete walker[i];

   delete trotter;
   delete P;

}

/**
 * initialize the walkers
 */
void AFQMC::SetupWalkers(){

   walker.resize(Nw);

   walker[0] = new Walker(n_trot);

   //read in the initial walker state
   char filename[200];
   sprintf(filename,"input/J1J2/%dx%d/J2=0.%d/PsiW/DT=%d.mps",Global::gLx(),Global::gLy(),Global::gj2(),Global::gD());

   walker[0]->read(filename);

   walker[0]->sOverlap(Psi0);

   walker[0]->fill_xyz();

   walker[0]->sEL(Psi0);

   walker[0]->sVL(*trotter,Psi0);
   walker[0]->sWeight(1.0);

   //now copy this to all the other walkers
   for(int cnt = 1;cnt < walker.size();cnt++)
      walker[cnt] = new Walker(*walker[0]);

}

void AFQMC::walk(const int steps){

   complex<double> EP = gEP();

   char filename[200];
   sprintf(filename,"output/J1J2/%dx%d/J2=0.%d/DT=%d.txt",Global::gLx(),Global::gLy(),Global::gj2(),Global::gD());

   ofstream output(filename,ios::trunc);

   output << "#Step\t\tE_P\t\tE_T\t" << endl;
   output.close();

#ifdef _DEBUG
   cout << "Energy at start = " << EP << endl;
   cout << "---------------------------------------------------------" << endl;
#endif

   for(int step = 1;step <= steps;step++){

      //Propagate the walkers of each rank separately --> no MPI in that function
      double wsum = Propagate();

      //Form the total sum of the walker weights and calculate the scaling for population control
      double avgw = wsum / (double)walker.size();

      double scaling = Nw / wsum;

      double ET = log(scaling)/dtau;

      EP = gEP();

#ifdef _DEBUG
      cout << "        Step = " << step << endl;
      cout << "   # walkers = " << walker.size() << endl;
      cout << " avg(weight) = " << avgw << endl;
      cout << "         E_P = " << EP << endl;
      cout << "         E_T = " << ET << endl;
      cout << "---------------------------------------------------------" << endl;
#endif

      write(step,std::real(EP), ET);

      //Based on scaling, first control the population on each rank separately, and then balance the population over the ranks (uses MPI)
      PopulationControl(scaling);

      double min_en = 0.0;
      double min_ov = 1.0;

      for(int i = 0;i < walker.size();++i){

         if(min_en > std::real(walker[i]->gEL()))
            min_en = std::real(walker[i]->gEL());

         if(min_ov > std::abs(walker[i]->gOverlap()))
            min_ov = std::abs(walker[i]->gOverlap());

      }

#ifdef _DEBUG
      cout << "Minimal Energy:\t" << min_en << endl;
      cout << "Minimal Overlap:\t" << min_ov << endl;
#endif

   }

}

/**
 * Here the trotter terms, propagator terms are applied to every walker individually.
 */
double AFQMC::Propagate(){

   double sum = 0.0;
   double width = sqrt(2.0/dtau);

   int num_rej = 0;

   for(int i = 0;i < walker.size();i++){

      Global::backup_walker.copy_essential(*walker[i]);

      //now loop over the auxiliary fields:
      for(int k = 0;k < n_trot;++k)
         for(int r = 0;r < 3;++r){

            double x = Global::RN.normal();

            complex<double> shift = walker[i]->gVL(k,r);

            //set the values
            P->set(x + shift,k,r);

            //and fill the propagator
            P->fill(*trotter);

            //and apply it to the walker:
            walker[i]->propagate(*P);

         }

      walker[i]->sOverlap(Psi0);

      complex<double> prev_EL = walker[i]->gEL();

      walker[i]->fill_xyz();

      walker[i]->sEL(Psi0);

      complex<double> EL = walker[i]->gEL();

      if( (std::real(EL) < std::real(prev_EL) - width) || (std::real(EL) > std::real(prev_EL) + width) ){//very rare event, will cause numerical unstability

         num_rej++;

         //copy the state back!
         walker[i]->copy_essential(Global::backup_walker);

      }
      else{//go on

         double scale = exp(-0.5 * dtau * std::real(EL + prev_EL));

         walker[i]->multWeight(scale);

         walker[i]->sVL(*trotter,Psi0);

      }

      sum += walker[i]->gWeight();

   }

   return sum;

}

/**
 * redistribute the weights to stabilize the walk, keep the population in check
 */
void AFQMC::PopulationControl(double scaling){

   double minw = 1.0;
   double maxw = 1.0;

   double sum = 0.0;

   for(int i = 0;i < walker.size();i++){

      walker[i]->multWeight(scaling);

      double weight = walker[i]->gWeight();

      if(weight < minw)
         minw = weight;

      if(weight > maxw)
         maxw = weight;

      if (weight < 0.25){ //Energy doesn't change statistically

         int nCopies = (int) ( weight + Global::RN());

         if(nCopies == 0){

#ifdef _DEBUG
            cout << "Walker with weight " << weight << " will be deleted." << endl;
#endif

            delete walker[i];

            walker.erase(walker.begin() + i);

         }
         else
            walker[i]->sWeight(1.0);

      }

      if(weight > 1.5){ //statically energy doesn't change

         int nCopies =(int) ( weight + Global::RN());
         double new_weight = weight / (double) nCopies;

         walker[i]->sWeight(new_weight);

#ifdef _DEBUG
         cout << "Walker with weight " << weight << " will be " << nCopies << " times copied with weight " << new_weight << "." << endl;
#endif

         for(int n = 1;n < nCopies;++n){

            Walker *nw = new Walker(*walker[i]);

            walker.push_back(nw);

         }

      }

      sum += weight;

   }

#ifdef _DEBUG
   cout << endl;
   cout << "total weight:\t" << sum << endl;
   cout << endl;

   cout << "The min. encountered weight is " << minw << " ." << endl;
   cout << "The max. encountered weight is " << maxw << " ." << endl;
#endif

}

/**
 * @return the total projected energy of the walkers at a certain timestep
 */
complex<double> AFQMC::gEP(){

   complex<double> projE_num = 0.0;
   complex<double> projE_den = 0.0;

   for(int wi = 0;wi < walker.size();wi++){

      complex<double> w_loc_en = walker[wi]->gEL(); // <Psi_T | H | walk > / <Psi_T | walk >

      //For the projected energy
      projE_num   += walker[wi]->gWeight() * w_loc_en;
      projE_den += walker[wi]->gWeight();

   }

   complex<double> EP = 0.0;

   EP = projE_num / projE_den;

   return EP;

}

void AFQMC::write(const int step,const double EP, const double ET){

   char filename[200];
   sprintf(filename,"output/J1J2/%dx%d/J2=0.%d/DT=%d.txt",Global::gLx(),Global::gLy(),Global::gj2(),Global::gD());

   ofstream output(filename,ios::app);
   output.precision(10);
   output << step << "\t\t" << walker.size() << "\t" << EP << "\t\t" << ET << endl;
   output.close();

}
