#include "QCD.h"


QCD::QCD(){
  setQCDInitialCondition();
  setNf(3);
}

QCD::QCD(const int nf){
  setQCDInitialCondition();
  setNf(nf);
}

QCD::~QCD(){
  delete _pdf;
  delete _Z0;
}

void QCD::setQCDInitialCondition(){
  _TF=0.5;
  _Nc=3.0; _CA=3.0; _CF=4.0/3.0;

  _Z0 = new Particle("Z0");

  //Download the pdfset from https://lhapdf.hepforge.org/pdfsets.html and put it into the folder which can be found out by lhapdf-config --datadir
  const int imem=0;  //load file member
  _pdf = LHAPDF::mkPDF("CT18NLO",imem);  //or CT14nlo
  // _pdf = LHAPDF::mkPDF("NNPDF23_lo_as_0130_qed",imem);  //same PDF with the default MadGraph model
}

void QCD::setNf(const int nf){
  _nf=nf;
  setQCDParameters();
  // std::cout << "_nf=" << _nf << std::endl;
}

void QCD::setQCDParameters(){
  //beta_i
  _beta0=11.0*_Nc/3.0-2.0*_nf/3.0;
  _beta1=34.0*_Nc*_Nc/3.0-10.0*_Nc*_nf/3.0-2.0*_CF*_nf;
  _beta2=2857.0*_Nc*_Nc*_Nc/54.0\
        +_TF*_nf*(2.0*_CF*_CF-205.0*_CF*_Nc/9.0-1415.0*_Nc*_Nc/27.0)\
        +_TF*_TF*_nf*_nf*(44.0*_CF/9.0+158.0*_Nc/27.0);

  //gamma^cusp_i
  _gcusp0=4.0;
  _gcusp1=(268.0/9.0-4.0*M_PI*M_PI/3.0)*_Nc-40.0*_nf/9.0;
  _gcusp2=_Nc*_Nc*(490.0/3.0-536.0*M_PI*M_PI/27.0+44.0*pow(M_PI,4.0)/45.0+88.0*gsl_sf_zeta(3.0)/3.0)\
         +_Nc*_TF*_nf*(-1672.0/27.0+160.0*M_PI*M_PI/27.0-224.0*gsl_sf_zeta(3.0)/3.0)\
         +_CF*_TF*_nf*(-220.0/3.0+64.0*gsl_sf_zeta(3.0)) - 64.0*_TF*_TF*_nf*_nf/27.0;

  //non-cusp anomalous dimension for the hard function
  _gHq0=-3.0*_CF;
  _gHq1=_CF*_CF*(-1.5+2.0*M_PI*M_PI-24.0*gsl_sf_zeta(3.0))\
       +_CF*_Nc*(-961.0/54.0-11.0*M_PI*M_PI/6.0+26.0*gsl_sf_zeta(3.0))\
       +_CF*_TF*_nf*(130.0/27.0+2.0*M_PI*M_PI/3.0);
  _gHg0=-_beta0;
  _gHg1=_Nc*_Nc*(-692.0/27.0+11.0*M_PI*M_PI/18.0+2.0*gsl_sf_zeta(3.0))\
       +_Nc*_TF*_nf*(256.0/27.0-2.0*M_PI*M_PI/9.0)+4.0*_CF*_TF*_nf;
}


//####################################################################################################
//In 153-154 pages of [S. Navas et al. (Particle Data Group), Phys. Rev. D 110, 030001 (2024)]
//Λ_QCD
double QCD::LambdaQCD(const int nloop){
  //LO with 1-loop β-function coefficient
  // double alphaS_mZ2=0.118;
  // double b0=_beta0/(4.*M_PI);
  // _lambdaQCD=_Z0->Mass()/std::pow(std::exp(1./(b0*alphaS_mZ2)),0.5);
  if(nloop==1){
    if(_nf==3){_lambdaQCD=0.2457484;}
    if(_nf==4){_lambdaQCD=0.1530858;}
    if(_nf==5){_lambdaQCD=0.0878275;}
  }
  //NLO with 2-loop β-function coefficient
  if(nloop==2){
    if(_nf==3){_lambdaQCD=0.7376440;}
    if(_nf==4){_lambdaQCD=0.4362115;}
    if(_nf==5){_lambdaQCD=0.2262347;}
  }

  return _lambdaQCD;
}


//####################################################################################################
//Splitting function
double QCD::Pqg(const double z){
  double res=0.0;
  if(z>0. && z<=1.) res=z*z-z+0.5;
  return res;
}

double QCD::Pqg(const double z, const double x, const double Q){
  double res=0.0;
  if(z>x && z<1.) res=Pqg(z)*pdf(21,x/z,Q)/z;
  return res;
}

double QCD::Pgq(const double z){
  double res=0.0;
  if(z>0. && z<1.) res=_CF*(z-2.+2./z);
  return res;
}

double QCD::Pgq(const double z, const double x, const double Q){
  double res=0.0;
  if(z>x && z<1.){
    for(int i=1; i<=_nf; i++){
      res+=pdf(i,x/z,Q);
    }
    res*=Pgq(z)/z;
  }
  return res;
}

double QCD::Pgqb(const double z, const double x, const double Q){
  double res=0.0;
  if(z>x && z<1.){
    for(int i=1; i<=_nf; i++){
      res+=pdf(-i,x/z,Q);
    }
    res*=Pgq(z)/z;
  }
  return res;
}

double QCD::Pqq(const int flavor, const double z, const double x, const double Q){
  double res=0.0;
  if(z>0. && z<1.){
    double qsumx=pdf(flavor,x,Q);
    res=1.5*qsumx;  //make sure the range for z integration is from 0 to 1
    double qsum=0.0;
    if(z>x) qsum=pdf(flavor,x/z,Q);
    res+=((1.+z*z)*qsum/z-2.*qsumx)/((1.-z));
    res*=_CF;
  }
  return res;
}

double QCD::Pgg(const double z, const double x, const double Q){
  double res=0.0;
  if(z>0. && z<1.){
    res=(11.*_Nc/6.-_nf/3.)*pdf(21,x,Q);  //make sure the range for z integration is from 0 to 1
    if(z>x){
      res+=2.*_Nc*((z*z-z+1.)*(z*z-z+1.)*pdf(21,x/z,Q)/(z*z)-pdf(21,x,Q))/(1.-z);
    }
    else{
      res+=-2.*_Nc*pdf(21,x,Q)/(1.-z);
    }
  }
  return res;
}
