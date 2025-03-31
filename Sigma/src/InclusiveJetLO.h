#ifndef INCLUSIVEJETLO_H
#define INCLUSIVEJETLO_H

#include <iostream>
#include <iomanip>
#include <cmath>
#include <gsl/gsl_math.h>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_vegas.h>
// #define MC_MISER
// #include <gsl/gsl_monte_miser.h>
#include "QCD.h"
#include "DefineInput.cpp"


class InclusiveJetLO{

 private:

  Particle *_parton;                      //particle class
  int    _nf;                             //quarks flavors
  double _unitGeV2pb;                     //GeV^2 pb from PDG 2024
  double _preXsection;                    //prefactor for cross secton calculation
  double _Ecm;                            //_Ecm=sqrt{_s} in CoM frame
  double _s;                              //collisional energy square
  double _mP2;                            //parton mass square
  double _y3c,_y4d;                       //rapidity of final jet
  double _x1,_x2,_shat,_that,_uhat;       //kinematics for hard scattering
  double _pdf1aplus[30],_pdf1aminus[30];  //PDFs with momentum fraction x1
  double _pdf2bplus[30],_pdf2bminus[30];  //PDFs with momentum fraction x2
  double _Dz3c[30],_Dz4d[30];             //JFFs for different parton numbering scheme

  //Variables for Monte Carlo
  const gsl_rng_type *_Type;
  gsl_rng *_rng;
#ifdef MC_MISER
  gsl_monte_miser_state *_mcStatus;
#else
  gsl_monte_vegas_state *_mcStatus;
#endif
  gsl_monte_function _mcFun;
  size_t _dim, _calls;
  double *_rangemin, *_rangemax;


 protected:

  QCD *_qcd;            //QCD class
  double _lambdaQCD;    //Λ_QCD with quark flavors dependence
  double _alphasPDF;    //AlphaS from PDF
  double _pT,_pTscale;  //initial jet momentum and scale to control the error range


 public:

  InclusiveJetLO();
  InclusiveJetLO(const double Ecm, Particle *Parton);
  virtual ~InclusiveJetLO();

  //pT and pT scale
  inline void   setpT(const double pT) {_pT=pT;}
  inline double getpT() const {return _pT;}
  inline void   setpTScale(const double pTscale) {_pTscale=pTscale;}
  //QCD parameters
  void setNf(const int nf);
  void setLambdaQCD(const int nloop);

  //Kinematic calculations
  void MomentumFractions(const double y3c, const double y4d);
  void Mandelstam(const double y3c, const double y4d);
  void setKinematics(const double y3c, const double y4d);
  void getKinematics();
  void getFourMomenta();

  //Kinematic cuts
  inline bool KinematicsQcut() const {return _x1>0.0 && _x1<1.0 && _x2>0.0 && _x2<1.0;}  //return whether it is allowed kinematics

  //PDFs: parton distribution functions
  void   setPDF1();
  void   setPDF2();
  double getPDF1(const int i);
  double getPDF2(const int i);
  //JFFs: jet fragmentation functions
  void   setJFF3();
  void   setJFF4();
  double getJFF3(const int i);
  double getJFF4(const int i);

  //Prefactor for cross section calculation at leading order
  void   setXsection();
  inline double getXsection() const {return _preXsection;}
  //Amplitudes squared
  inline double M2qqp2qqp(const double &s, const double &t, const double &u);    //qq'->qq'
  inline double M2qq2qq(const double &s, const double &t, const double &u);      //qq->qq
  inline double M2qqb2qpqpb(const double &s, const double &t, const double &u);  //qqbar->q'q'bar
  inline double M2qqb2qqb(const double &s, const double &t, const double &u);    //qqbar->qqbar
  inline double M2qqb2gg(const double &s, const double &t, const double &u);     //qqbar->gg
  inline double M2gg2qqb(const double &s, const double &t, const double &u);     //gg->qqbar
  inline double M2gq2gq(const double &s, const double &t, const double &u);      //gq->gq
  inline double M2gg2gg(const double &s, const double &t, const double &u);      //gg->gg

  //Leading order cross section for different channels
  double SigmaLOqqp2qqp();
  double SigmaLOqq2qq();
  double SigmaLOqqb2qpqpb();
  double SigmaLOqqb2qqb();
  double SigmaLOqqb2gg();
  double SigmaLOgg2qqb();
  double SigmaLOgq2gq();
  double SigmaLOgq2gq_Q();
  double SigmaLOgq2gq_G();
  double SigmaLOgg2gg();
  //Leading order kinematic parameters
  void   setParametersLO(const double y3c, const double y4d);
  //Leading order cross section
  double getSigmaLO();
  double getSigmaLOQuark();
  double getSigmaLOGluon();


  //Monte Carlo
  void setCalls(const size_t calls);
  void setupMC(double (*func)(double *, size_t, void *), size_t dim, double *rangemin, double *rangemax, void *para);
  void calculateMC(double &res, double &err);
  void freeMC();

};

#endif  //INCLUSIVEJETLO_H
