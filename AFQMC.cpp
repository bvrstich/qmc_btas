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

   walker[0] = new Walker(Global::gL(),Global::gd(),1.0,n_trot);

   //read in the initial walker state
   char filename[200];
   sprintf(filename,"input/J1J2/%dx%d/J2=0.%d/PsiW/DT=%d.mps",Global::gLx(),Global::gLy(),Global::gj2(),Global::gD());

   walker[0]->read(filename);
/*
   walker[0]->sOverlap(Psi0);
   walker[0]->sEL(ham,Psi0);
   walker[0]->sVL(*trotter,Psi0);

   //now copy this to all the other walkers
   for(int cnt = 1;cnt < theWalkers.size();cnt++)
      theWalkers[cnt] = new Walker(*theWalkers[0]);
*/
}
/*
void AFQMC::Walk(const int steps){

   complex<double> EP = gEP();

   char filename[200];
   sprintf(filename,"output/J1J2/%dx%d/J2=0.%d/DT=%d.txt",Global::gLx(),Global::gLy(),Global::gj2(),mps.gD());

   ofstream output(filename,ios::trunc);

   output << "#Step\t\tE_P\t\tE_T\t" << endl;
   output.close();

   cout << "Energy at start = " << EP << endl;
   cout << "---------------------------------------------------------" << endl;

   for (int step=1;step <= steps;step++){

      //Propagate the walkers of each rank separately --> no MPI in that function
      double wsum = Propagate();

      //Form the total sum of the walker weights and calculate the scaling for population control
      double avgw = wsum / (double)theWalkers.size();

      double scaling = Nwalkers / wsum;

      double ET = log(scaling)/dtau;

      EP = gEP();

      cout << "        Step = " << step << endl;
      cout << "   # walkers = " << theWalkers.size() << endl;
      cout << " avg(weight) = " << avgw << endl;
      cout << "         E_P = " << EP << endl;
      cout << "         E_T = " << ET << endl;
      cout << "---------------------------------------------------------" << endl;

      write(step,theWalkers.size(),std::real(EP), ET);

      //Based on scaling, first control the population on each rank separately, and then balance the population over the ranks (uses MPI)
      PopulationControl(scaling);

      double min_en = 0.0;
      double min_ov = 1.0;

      for(int i = 0;i < theWalkers.size();++i){

         if(min_en > std::real(theWalkers[i]->gEL()))
            min_en = std::real(theWalkers[i]->gEL());

         if(min_ov > std::abs(theWalkers[i]->gOverlap()))
            min_ov = std::abs(theWalkers[i]->gOverlap());

      }

      cout << "Minimal Energy:\t" << min_en << endl;
      cout << "Minimal Overlap:\t" << min_ov << endl;

   }

}
*/
/**
 * Here the trotter terms, propagator terms are applied to every walker individually.
 */
 /*
double AFQMC::Propagate(){

   double sum = 0.0;

   for(int walker=0; walker < theWalkers.size(); walker++){

      cout << walker << endl;

      //now loop over the auxiliary fields:
      for(int k = 0;k < n_trot;++k)
         for(int r = 0;r < 3;++r){

            double x = Tools::RN.normal();

            complex<double> shift = theWalkers[walker]->gVL(k,r);

            //set the values
            P->set(x + shift,k,r);

            //and fill the propagator
            //      P->fill(*trotter);

            //and apply it to the walker
            //       theWalkers[walker]->propagate(*P);

         }

      theWalkers[walker]->sOverlap(Psi0);

      complex<double> prev_EL = theWalkers[walker]->gEL();

      //  theWalkers[walker]->sEL(ham,Psi0);

      complex<double> EL = theWalkers[walker]->gEL();

      double scale = exp(-0.5 * dtau * std::real(EL + prev_EL));

      theWalkers[walker]->multWeight(scale);

      //   theWalkers[walker]->sVL(*trotter,Psi0);

      sum += theWalkers[walker]->gWeight();

   }

   return sum;

}
*/
/*
void AFQMC::PopulationControl(const double scaling){

   double minw = 1.0;
   double maxw = 1.0;

   double sum = 0.0;

   for(int walker = 0;walker < theWalkers.size();walker++){

      theWalkers[walker]->multWeight(scaling);

      double weight = theWalkers[walker]->gWeight();

      if(weight < minw)
         minw = weight;

      if(weight > maxw)
         maxw = weight;

      if (weight < 0.25){ //Energy doesn't change statistically

         int nCopies = (int) ( weight + Tools::RN());

         if(nCopies == 0){

            cout << "Walker with weight " << weight << " will be deleted." << endl;

            delete theWalkers[walker];

            theWalkers.erase(theWalkers.begin() + walker);

         }
         else
            theWalkers[walker]->sWeight(1.0);

      }

      if(weight > 1.5){ //statically energy doesn't change

         int nCopies =(int) ( weight + Tools::RN());
         double new_weight = weight / (double) nCopies;

         theWalkers[walker]->sWeight(new_weight);

         cout << "Walker with weight " << weight << " will be " << nCopies << " times copied with weight " << new_weight << "." << endl;

         for(int i = 1;i < nCopies;++i){

            Walker *nw = new Walker(*theWalkers[walker]);

            theWalkers.push_back(nw);

         }

      }

      sum += weight;

   }

   cout << endl;
   cout << "total weight:\t" << sum << endl;
   cout << endl;

   cout << "The min. encountered weight is " << minw << " ." << endl;
   cout << "The max. encountered weight is " << maxw << " ." << endl;

}

complex<double> AFQMC::gEP(){

   complex<double> projE_num = 0.0;
   complex<double> projE_den = 0.0;

   for(int walker = 0;walker < theWalkers.size();walker++){

      const complex<double> walkerEnergy = theWalkers[walker]->gEL(); // <Psi_T | H | walk > / <Psi_T | walk >

      //For the projected energy
      projE_num   += theWalkers[walker]->gWeight() * walkerEnergy;
      projE_den += theWalkers[walker]->gWeight();

   }

   complex<double> EP = 0.0;

   EP = projE_num / projE_den;

   return EP;

}

void AFQMC::write(const int step,const int nwalkers,const double EP, const double ET){

   char filename[100];
   sprintf(filename,"output/J1J2/4x4/J2=0.0/DW%d.txt",DW);

   ofstream output(filename,ios::app);
   output.precision(10);
   output << step << "\t\t" << nwalkers << "\t" << EP << "\t\t" << ET << endl;
   output.close();

}
*/
