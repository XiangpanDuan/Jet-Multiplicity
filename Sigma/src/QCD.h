#ifndef QCD_H
#define QCD_H

#include <string>
#include <cmath>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_zeta.h>
#include "LHAPDF/LHAPDF.h"
#include "Particle.h"


class QCD{

 private:

  LHAPDF::PDF *_pdf;               //download LHAPDF from https://lhapdf.hepforge.org/install.html
  Particle *_Z0;
  int    _nf;                      //flavors of quarks, 3 by default
  double _lambdaQCD;
  double _TF;
  double _Nc,_CA,_CF;              //number of colors
  double _beta0,_beta1,_beta2;     //beta function coefficient
  double _gcusp0,_gcusp1,_gcusp2;  //_gcusp=gamma^cusp
  double _gHq0,_gHg0,_gHq1,_gHg1;  //non-cusp anomalous dims for the hard function


 public:

  QCD();
  QCD(const int nf);
  ~QCD();
  void setQCDInitialCondition();
  void setNf(const int nf);
  void setQCDParameters();

  //Color factor
  inline int    Nf() const {return _nf;}
  inline double TF() const {return _TF;}
  inline double Nc() const {return _Nc;}
  inline double CA() const {return _CA;}
  inline double CF() const {return _CF;}

  //PDF, alphas, and Λ_QCD
  inline double pdf(const int f, const double x, const double Q) const {return _pdf->xfxQ(f,x,Q)/x;}  //PDF as a function of flavor with momentum fraction x and factorization scale Q
  inline double alphas(const double Q) const {return _pdf->alphasQ(Q);}                               //strong coupling constant alphas as a function of Q
  double LambdaQCD(const int nloop);

  //Splitting function
  double Pqg(const double z);
  double Pqg(const double z, const double x, const double Q);
  double Pgq(const double z);
  double Pgq(const double z, const double x, const double Q);
  double Pgqb(const double z, const double x, const double Q);
  double Pgg(const double z, const double x, const double Q);
  double Pqq(const int flavor, const double z, const double x, const double Q);

};

#endif  //QCD_H
