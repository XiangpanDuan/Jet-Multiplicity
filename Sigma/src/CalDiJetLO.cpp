#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>
#include <iomanip>
#include <cstdlib>
#include <chrono>
#include "DiJetLO.h"


namespace CalDiJetLO{

  double SigQuark(double range[], size_t dim, void *params){
    DiJetLO *jetp = (DiJetLO*) params;
    // std::string type="quark";
    jetp->setpT(range[0]);
    jetp->setLOParameters(range[1],range[2]);

    return jetp->SigmaLOQuark();
  }

  double SigGluon(double range[], size_t dim, void *params){
    DiJetLO *jetp = (DiJetLO*) params;
    // std::string type="gluon";
    jetp->setpT(range[0]);
    jetp->setLOParameters(range[1],range[2]);
    
    return jetp->SigmaLOGluon();
  }


  //------------------------------------------------------------
  //Cross section calculation: d^2sigma/dpTdy
  void CalSigma(){

    //Record start time
    auto starttime = std::chrono::high_resolution_clock::now();

    //Output files
    std::stringstream ss[3];
    std::string OutputString[3];
    std::ofstream OutputFile[3];
    for(int i=0; i<3; i++){
      ss[i] << "../Output/dsigma_dpTdy_" << Input::Ecm << "GeV_type" << i << ".dat";
      OutputString[i]=ss[i].str();
      OutputFile[i].open(OutputString[i], std::ofstream::out);
    }


    //Jet initialization
    Particle *Parton = new Particle(Input::name);
    DiJetLO  *Jet    = new DiJetLO(Input::Ecm,Parton);

    //Set initial conditions
    Jet->setpTScale(std::pow(2.,Input::scale));
    Jet->setNf(Input::nf);
    Jet->setLambdaQCD(Input::nloop);

    //Range of jet rapidity
    double rapmin3=-Input::rap3;
    double rapmax3= Input::rap3;
    double rapmin4=-Input::rap4;
    double rapmax4= Input::rap4;
    double drap3=rapmax3-rapmin3;  //observed jet rapidity

    //Final observed jet momentum
    int    pTnum=Input::pTnum;
    double pTmin=Input::pTmin;
    double pTmax=Input::pTmax;
    double pTbin=(pTmax-pTmin)/pTnum;
    for(int k=0; k<pTnum; k++){
      double pT=pTmin+pTbin/2.+pTbin*k;
      double pTbinmin=pT-pTbin/2.;
      double pTbinmax=pT+pTbin/2.;

      //Range of the integration: from rangemin to rangemax
      //pp collisions
      size_t dim=3;
      double rangemin[]={pTbinmin, rapmin3, rapmin4};
      double rangemax[]={pTbinmax, rapmax3, rapmax4};
      std::cout << pT << "  " << dim << std::endl;
      
      size_t calls = Input::calls;
      Jet->setCalls(calls);

      double resq,errq;
      double resg,errg;

      //Quark cross section
      Jet->setupMC(&SigQuark, dim, rangemin, rangemax, Jet);
      Jet->calculateMC(resq,errq);
      Jet->cleanupMC();
      // std::cout << "quark d^2sigma/dpTdy = " <<  resq/(pTbin*drap3)*0.001 << " nb with err = " << errq/(pTbin*drap3)*0.001 << " nb." << std::endl;
      OutputFile[1] << pT << "   " << resq/(pTbin*drap3)*0.001 << "   " << errq/(pTbin*drap3)*0.001 << std::endl;

      //Gluon cross section
      Jet->setupMC(&SigGluon, dim, rangemin, rangemax, Jet);
      Jet->calculateMC(resg,errg);
      Jet->cleanupMC();
      // std::cout << "gluon d^2sigma/dpTdy = " <<  resg/(pTbin*drap3)*0.001 << " nb with err = " << errg/(pTbin*drap3)*0.001 << " nb." << std::endl;
      OutputFile[2] << pT << "   " << resg/(pTbin*drap3)*0.001 << "   " << errg/(pTbin*drap3)*0.001 << std::endl;

      //Total cross section
      std::cout << "total d^2sigma/dpTdy = " <<  (resq+resg)/(pTbin*drap3)*0.001 << " nb with err = " << (errq+errg)/(pTbin*drap3)*0.001 << " nb." << std::endl;
      OutputFile[0] << pT << "   " << (resq+resg)/(pTbin*drap3)*0.001 << "   " << (errq+errg)/(pTbin*drap3)*0.001 << std::endl;


      //Record end time
      auto endtime=std::chrono::high_resolution_clock::now();
      std::chrono::duration<double> duration=endtime-starttime;
      std::cout << "Function execution time: " << duration.count() << " seconds" << std::endl;
    }
    delete Jet;

    for(int i=0; i<3; i++){
      OutputFile[i].close();
    }
    
  }

}



//####################################################################################################
//Main function
int main(){

  std::cout << "######################################################################" << std::endl;
  std::cout << "# Leading-Order Cross Section Calculation in pp collisions           #" << std::endl;
  std::cout << "######################################################################" << std::endl;

  if(MODE!=0){
    std::cerr << "Error in pp collisions!!!" << std::endl;
    std::exit(EXIT_FAILURE);
  }


  CalDiJetLO::CalSigma();

  return 0;
}