#ifndef DIJETLO_H
#define DIJETLO_H

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


class DiJetLO{

 private:

  Particle *_parton;
  unsigned int _nf;                       //flavors of quarks, 3 by default
  double _preXsection;                    //prefactor for calculating cross secton
  double _Ecm;                            //_Ecm=sqrt{_s} in CoM frame
  double _s;                              //collisional energy square
  double _mP2;                            //parton mass square
  double _unitGeV2pb;                     //GeV^2 pb from PDG 2024
  double _y3c,_y4d;                       //rapidity of final jet
  double _x1,_x2,_shat,_that,_uhat;       //kinematics for hard scattering
  double _pdf1aplus[30],_pdf1aminus[30];  //PDF with momentum fraction x1
  double _pdf2bplus[30],_pdf2bminus[30];  //PDF with momentum fraction x2
  double _Dz3c[30],_Dz4d[30];             //JFF for different parton numbering scheme

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

  QCD *_qcd;
  double _lambdaQCD;    //Λ_QCD depends on quark flavour _nf
  double _pT,_pTscale;  //initial jet momentum and scale to control the errer bar


 public:

  DiJetLO();
  DiJetLO(const double Ecm, Particle *Parton);
  virtual ~DiJetLO();

  //pT and pT scale
  inline void setpT(const double pT) {_pT=pT;}
  inline void setpTScale(const double pTscale) {_pTscale=pTscale;}
  //QCD
  void setNf(const unsigned int nf);
  void setLambdaQCD(const unsigned int nloop);

  //Kinematic calculations
  void MomentumFractions(const double y3c, const double y4d);
  void Mandelstam(const double y3c, const double y4d);
  void setKinematics(const double y3c, const double y4d);
  void getKinematics();
  void getFourMomenta();

  //Kinematic cuts
  inline bool KinematicsQcut() const {return _x1>0.0 && _x1<1.0 && _x2>0.0 && _x2<1.0;}  //return whether it is allowed kinematics

  //PDF: Parton distribution function
  void   setPDF1();
  void   setPDF2();
  double getPDF1(const int i);
  double getPDF2(const int i);
  //JFF: Jet fragmentation function
  void   setJetFF3();
  void   setJetFF4();
  double getJetFF3(const int i);
  double getJetFF4(const int i);

  //The prefactor*M^2 gives DiJetLO cross section
  double alphasPDF(const double pT);
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

  //Differential cross section at leading order
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
  void   setLOParameters(const double y3c, const double y4d);
  //Leading order cross section
  double SigmaLO();
  double SigmaLOQuark();
  double SigmaLOGluon();


  //Monte Carlo
  void setCalls(const size_t calls);
  void setupMC(double (*func)(double *, size_t, void *), size_t dim, double *rangemin, double *rangemax, void *para);
  void calculateMC(double &res, double &err);
  void cleanupMC();

};

#endif  //DIJETLO_H
