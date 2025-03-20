#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>
#include <iomanip>
#include <cstdlib>
#include <chrono>
#include "InclusiveJetLO.h"


namespace CalculateLO{

  double SigmaQuark(double range[], size_t dim, void *params){
    DiJetLO *jetq = (DiJetLO*) params;
    // std::string type="quark";
    jetq->setpT(range[0]);
    jetq->setLOParameters(range[1],range[2]);

    return jetq->SigmaLOQuark();
  }

  double SigmaGluon(double range[], size_t dim, void *params){
    DiJetLO *jetg = (DiJetLO*) params;
    // std::string type="gluon";
    jetg->setpT(range[0]);
    jetg->setLOParameters(range[1],range[2]);

    return jetg->SigmaLOGluon();
  }

  double dSigmaQuark(double range[], size_t dim, void *params){
    DiJetLO *jetq = (DiJetLO*) params;
    // std::string type="quark";
    jetq->setLOParameters(range[0],range[1]);

    return jetq->SigmaLOQuark();
  }

  double dSigmaGluon(double range[], size_t dim, void *params){
    DiJetLO *jetg = (DiJetLO*) params;
    // std::string type="gluon";
    jetg->setLOParameters(range[0],range[1]);

    return jetg->SigmaLOGluon();
  }


  //------------------------------------------------------------
  //Cross section calculation: sigma
  void CalSigmaLO(){

    //Record start time
    auto starttime = std::chrono::high_resolution_clock::now();

    //Output files
    std::stringstream ss[3];
    std::string OutputString[3];
    std::ofstream OutputFile[3];
    for(int i=0; i<3; i++){
      if(Input::rap4==10.) ss[i] << "../Output/sigma_" << Input::Ecm << "GeV_inclusivejet_type" << i << ".dat";
      if(Input::rap4!=10.) ss[i] << "../Output/sigma_" << Input::Ecm << "GeV_dijet_type" << i << ".dat";
      OutputString[i]=ss[i].str();
      OutputFile[i].open(OutputString[i], std::ofstream::out);
      OutputFile[i] << "# pT(GeV) sigma(pb) error(pb)" << std::endl;
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
    // double drap3=rapmax3-rapmin3;  //observed jet rapidity

    //Final observed jet momentum
    int    pTnum=Input::pTnum;
    double pTmin=Input::pTmin;
    double pTmax=Input::pTmax;
    double pTbin=(pTmax-pTmin)/pTnum;
    //pT loop
    for(int k=0; k<pTnum; k++){
      double pT=pTmin+pTbin*k;
      double pTbinmin=pT-pTbin/2.;
      double pTbinmax=pT+pTbin/2.;

      //Range of the integration: from rangemin to rangemax
      //pp collisions
      size_t dim=3;
      double rangemin[]={pTbinmin, rapmin3, rapmin4};
      double rangemax[]={pTbinmax, rapmax3, rapmax4};
      std::cout << "pT=" << pT << "GeV,   dim=" << dim << std::endl;

      size_t calls = Input::calls;
      Jet->setCalls(calls);

      double resq,errq;
      double resg,errg;

      //Quark cross section
      Jet->setupMC(&SigmaQuark, dim, rangemin, rangemax, Jet);
      Jet->calculateMC(resq,errq);
      Jet->cleanupMC();
      // std::cout << "quark sigma = " << resq << " pb with err = " << errq << " pb." << std::endl;
      OutputFile[1] << pT << "   " << resq << "   " << errq << std::endl;

      //Gluon cross section
      Jet->setupMC(&SigmaGluon, dim, rangemin, rangemax, Jet);
      Jet->calculateMC(resg,errg);
      Jet->cleanupMC();
      // std::cout << "gluon sigma = " << resg << " pb with err = " << errg << " pb." << std::endl;
      OutputFile[2] << pT << "   " << resg << "   " << errg << std::endl;

      //Total cross section
      std::cout << "total sigma = " << resq+resg << " pb with err = " << errq+errg << " pb." << std::endl;
      OutputFile[0] << pT << "   "  << resq+resg << "   " << errq+errg << std::endl;


      //Record end time
      auto endtime=std::chrono::high_resolution_clock::now();
      std::chrono::duration<double> duration=endtime-starttime;
      std::cout << "Execution time: " << duration.count() << " seconds" << std::endl;

    }//End of pT loop

    delete Jet;

    //Close output files
    for(int i=0; i<3; i++){
      OutputFile[i].close();
    }

  }


  //------------------------------------------------------------
  //Differential cross section calculation: dsigma/dpT
  void CalDSigmaLO(){

    //Record start time
    auto starttime = std::chrono::high_resolution_clock::now();

    //Output files
    std::stringstream ss[3];
    std::string OutputString[3];
    std::ofstream OutputFile[3];
    for(int i=0; i<3; i++){
      if(Input::rap4==10.) ss[i] << "../Output/dsigma_dpT_" << Input::Ecm << "GeV_inclusivejet_type" << i << ".dat";
      if(Input::rap4!=10.) ss[i] << "../Output/dsigma_dpT_" << Input::Ecm << "GeV_dijet_type" << i << ".dat";
      OutputString[i]=ss[i].str();
      OutputFile[i].open(OutputString[i], std::ofstream::out);
      OutputFile[i] << "# pT(GeV) dsigma/dpT(pb/GeV) error(pb/GeV)" << std::endl;
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
    // double drap3=rapmax3-rapmin3;  //observed jet rapidity

    //Final observed jet momentum
    int    pTnum=Input::pTnum;
    double pTmin=Input::pTmin;
    double pTmax=Input::pTmax;
    double pTbin=(pTmax-pTmin)/pTnum;
    //pT loop
    for(int k=0; k<pTnum; k++){
      double pT=pTmin+pTbin*k;

      //Range of the integration: from rangemin to rangemax
      //pp collisions
      size_t dim=2;
      double rangemin[]={rapmin3, rapmin4};
      double rangemax[]={rapmax3, rapmax4};
      std::cout << "pT=" << pT << "GeV,   dim=" << dim << std::endl;

      size_t calls = Input::calls;
      Jet->setCalls(calls);

      double resq,errq;
      double resg,errg;

      //Quark cross section
      Jet->setpT(pT);
      Jet->setupMC(&dSigmaQuark, dim, rangemin, rangemax, Jet);
      Jet->calculateMC(resq,errq);
      Jet->cleanupMC();
      // std::cout << "quark dsigma/dpT = " <<  resq << " pb with err = " << errq << " pb." << std::endl;
      OutputFile[1] << pT << "   " << resq << "   " << errq << std::endl;

      //Gluon cross section
      Jet->setpT(pT);
      Jet->setupMC(&dSigmaGluon, dim, rangemin, rangemax, Jet);
      Jet->calculateMC(resg,errg);
      Jet->cleanupMC();
      // std::cout << "gluon dsigma/dpT = " <<  resg << " pb with err = " << errg << " pb." << std::endl;
      OutputFile[2] << pT << "   " << resg << "   " << errg << std::endl;

      //Total cross section
      std::cout << "total dsigma/dpT = " << resq+resg << " pb with err = " << errq+errg << " pb." << std::endl;
      OutputFile[0] << pT << "   " << resq+resg << "   " << errq+errg << std::endl;


      //Record end time
      auto endtime=std::chrono::high_resolution_clock::now();
      std::chrono::duration<double> duration=endtime-starttime;
      std::cout << "Execution time: " << duration.count() << " seconds" << std::endl;

    }//End of pT loop

    delete Jet;

    //Close output files
    for(int i=0; i<3; i++){
      OutputFile[i].close();
    }

  }


  // //------------------------------------------------------------
  // //Test
  // void Test(){
  //   QCD *Test = new QCD();
  //   for(int k=0; k<25; k++){
  //     double pT=150.+100.*k;
  //     double RR=0.4;
  //     double alphas=Test->alphas(pT*RR);
  //     std::cout << pT << " " << pT*RR << " " << alphas << std::endl;
  //   }
  //   delete Test;
  // }


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

  // CalculateLO::CalSigmaLO();
  CalculateLO::CalDSigmaLO();

  // CalculateLO::Test();

  return 0;
}