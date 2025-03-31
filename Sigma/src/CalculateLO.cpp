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
    InclusiveJetLO *jetq = (InclusiveJetLO*) params;
    // std::string type="quark";
    jetq->setpT(range[0]);
    jetq->setParametersLO(range[1],range[2]);

    return jetq->getSigmaLOQuark();
  }

  double SigmaGluon(double range[], size_t dim, void *params){
    InclusiveJetLO *jetg = (InclusiveJetLO*) params;
    // std::string type="gluon";
    jetg->setpT(range[0]);
    jetg->setParametersLO(range[1],range[2]);

    return jetg->getSigmaLOGluon();
  }

  double dSigmaQuark(double range[], size_t dim, void *params){
    InclusiveJetLO *jetq = (InclusiveJetLO*) params;
    // std::string type="quark";
    jetq->setParametersLO(range[0],range[1]);

    return jetq->getSigmaLOQuark();
  }

  double dSigmaGluon(double range[], size_t dim, void *params){
    InclusiveJetLO *jetg = (InclusiveJetLO*) params;
    // std::string type="gluon";
    jetg->setParametersLO(range[0],range[1]);

    return jetg->getSigmaLOGluon();
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
      ss[i] << "../Output/sigma_" << Input::Ecm << "GeV_inclusivejet_type" << i << ".dat";
      OutputString[i]=ss[i].str();
      OutputFile[i].open(OutputString[i], std::ofstream::out);
      OutputFile[i] << "# pT(GeV) sigma(pb) error(pb)" << std::endl;
    }


    //Jet initialization
    Particle *Parton = new Particle(Input::name);
    InclusiveJetLO *Jet = new InclusiveJetLO(Input::Ecm,Parton);

    //Set initial condition
    Jet->setpTScale(std::pow(2.0,Input::scale));
    Jet->setNf(Input::nf);
    Jet->setLambdaQCD(Input::nloop);

    //Range of jet rapidity
    double rap3min=-Input::rap3;
    double rap3max= Input::rap3;
    double rap4min=-Input::rap4;
    double rap4max= Input::rap4;
    // double drap3=rap3max-rap3min;  //observed jet rapidity

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
      double rangemin[]={pTbinmin, rap3min, rap4min};
      double rangemax[]={pTbinmax, rap3max, rap4max};
      std::cout << "pT = " << pT << " GeV,   dim = " << dim << std::endl;

      size_t calls=Input::calls;
      Jet->setCalls(calls);

      double resq,errq;
      double resg,errg;

      //Quark cross section
      Jet->setupMC(&SigmaQuark, dim, rangemin, rangemax, Jet);
      Jet->calculateMC(resq,errq);
      Jet->freeMC();
      OutputFile[1] << pT << "   " << resq << "   " << errq << std::endl;
      // std::cout << "Quark sigma = " << resq << " pb with error = " << errq << " pb." << std::endl;

      //Gluon cross section
      Jet->setupMC(&SigmaGluon, dim, rangemin, rangemax, Jet);
      Jet->calculateMC(resg,errg);
      Jet->freeMC();
      OutputFile[2] << pT << "   " << resg << "   " << errg << std::endl;
      // std::cout << "Gluon sigma = " << resg << " pb with error = " << errg << " pb." << std::endl;

      //Total cross section
      OutputFile[0] << pT << "   "  << resq+resg << "   " << std::sqrt(errq*errq+errg*errg) << std::endl;
      std::cout << "Total sigma = " << resq+resg << " pb with error = " << std::sqrt(errq*errq+errg*errg) << " pb." << std::endl;


      //Record end time
      auto endtime=std::chrono::high_resolution_clock::now();
      std::chrono::duration<double> duration=endtime-starttime;
      std::cout << "Execution time: " << duration.count() << " seconds." << std::endl;

    }//End of pT loop

    delete Jet;

    //Close output files
    for(int i=0; i<3; i++){
      OutputFile[i].close();
    }

  }


  //------------------------------------------------------------
  //Differential cross section calculation: dsigma/dpT
  void CaldSigmaLO(){

    //Record start time
    auto starttime = std::chrono::high_resolution_clock::now();

    //Output files
    std::stringstream ss[3];
    std::string OutputString[3];
    std::ofstream OutputFile[3];
    for(int i=0; i<3; i++){
      ss[i] << "../Output/dsigma_dpT_" << Input::Ecm << "GeV_inclusivejet_type" << i << ".dat";
      OutputString[i]=ss[i].str();
      OutputFile[i].open(OutputString[i], std::ofstream::out);
      OutputFile[i] << "# pT(GeV) dsigma/dpT(pb/GeV) error(pb/GeV)" << std::endl;
    }


    //Jet initialization
    Particle *Parton = new Particle(Input::name);
    InclusiveJetLO *Jet = new InclusiveJetLO(Input::Ecm,Parton);

    //Set initial condition
    Jet->setpTScale(std::pow(2.0,Input::scale));
    Jet->setNf(Input::nf);
    Jet->setLambdaQCD(Input::nloop);

    //Range of jet rapidity
    double rap3min=-Input::rap3;
    double rap3max= Input::rap3;
    double rap4min=-Input::rap4;
    double rap4max= Input::rap4;
    // double drap3=rap3max-rap3min;  //observed jet rapidity

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
      double rangemin[]={rap3min, rap4min};
      double rangemax[]={rap3max, rap4max};
      std::cout << "pT = " << pT << " GeV,   dim = " << dim << std::endl;

      size_t calls = Input::calls;
      Jet->setCalls(calls);

      double resq,errq;
      double resg,errg;

      //Quark cross section
      Jet->setpT(pT);
      Jet->setupMC(&dSigmaQuark, dim, rangemin, rangemax, Jet);
      Jet->calculateMC(resq,errq);
      Jet->freeMC();
      OutputFile[1] << pT << "   " << resq << "   " << errq << std::endl;
      // std::cout << "Quark dsigma/dpT = " << resq << " pb with error = " << errq << " pb." << std::endl;

      //Gluon cross section
      Jet->setpT(pT);
      Jet->setupMC(&dSigmaGluon, dim, rangemin, rangemax, Jet);
      Jet->calculateMC(resg,errg);
      Jet->freeMC();
      OutputFile[2] << pT << "   " << resg << "   " << errg << std::endl;
      // std::cout << "Gluon dsigma/dpT = " << resg << " pb with error = " << errg << " pb." << std::endl;

      //Total cross section
      OutputFile[0] << pT << "   " << resq+resg << "   " << std::sqrt(errq*errq+errg*errg) << std::endl;
      std::cout << "Total dsigma/dpT = " << resq+resg << " pb with error = " << std::sqrt(errq*errq+errg*errg) << " pb." << std::endl;


      //Record end time
      auto endtime=std::chrono::high_resolution_clock::now();
      std::chrono::duration<double> duration=endtime-starttime;
      std::cout << "Execution time: " << duration.count() << " seconds." << std::endl;

    }//End of pT loop

    delete Jet;

    //Close output files
    for(int i=0; i<3; i++){
      OutputFile[i].close();
    }

  }


}



//####################################################################################################
//Main function
int main(){

  std::cout << "######################################################################" << std::endl;
  std::cout << "# Leading Order Cross Section Calculation in pp collisions           #" << std::endl;
  std::cout << "######################################################################" << std::endl;

  if(MODE!=0){
    std::cerr << "Error in pp collisions!!!" << std::endl;
    std::exit(EXIT_FAILURE);
  }

  // CalculateLO::CalSigmaLO();
  CalculateLO::CaldSigmaLO();


  return 0;
}