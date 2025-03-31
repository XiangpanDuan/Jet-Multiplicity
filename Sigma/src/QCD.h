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

  LHAPDF::PDF *_pdf;               //LHAPDF from https://lhapdf.hepforge.org/install.html
  Particle *_Z0;                   //parton Z0 class
  int    _nf;                      //quark flavors
  double _lambdaQCD;               //Λ_QCD
  double _TF;                      //color factor
  double _Nc,_CA,_CF;              //color factor
  double _beta0,_beta1,_beta2;     //beta function coefficient
  double _gcusp0,_gcusp1,_gcusp2;  //cusp anomalous dimensions (gcusp=gamma^cusp)
  double _gHq0,_gHg0,_gHq1,_gHg1;  //non-cusp anomalous dimensions for the hard function


 public:

  QCD();
  QCD(const int nf);
  ~QCD();
  void setInitialConditionQCD();
  void setNfQCD(const int nf);
  void setParametersQCD();

  //Quark flavors and color factor
  inline int    Nf() const {return _nf;}
  inline double TF() const {return _TF;}
  inline double Nc() const {return _Nc;}
  inline double CA() const {return _CA;}
  inline double CF() const {return _CF;}

  //PDFs, alphas, and Λ_QCD
  inline double PDFs(const int f, const double x, const double Q) const {return _pdf->xfxQ(f,x,Q)/x;}  //PDFs with momentum fraction x and factorization scale Q
  inline double AlphaS(const double Q) const {return _pdf->alphasQ(Q);}                                //strong coupling alphas as a function of Q
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
