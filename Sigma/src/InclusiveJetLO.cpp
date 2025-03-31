#include "InclusiveJetLO.h"


InclusiveJetLO::InclusiveJetLO(){
  _nf=3;
  _lambdaQCD=0.2457484;
  _qcd = new QCD();
  _pTscale=1.0;
  _unitGeV2pb=3.893793721e8;  //1 (ℏ𝑐)^2 = 3.893793721e8 GeV^2 pb from PDG 2024
  _calls=50000;
  _Ecm=5020.;  //GeV
  _s=_Ecm*_Ecm;
  _parton = new Particle();
  _mP2=0.0;
}

InclusiveJetLO::InclusiveJetLO(const double Ecm, Particle *Parton){
  _nf=3;
  _lambdaQCD=0.2457484;
  _qcd = new QCD();
  _pTscale=1.0;
  _unitGeV2pb=3.893793721e8;
  _calls=50000;
  _Ecm=Ecm;
  _s=_Ecm*_Ecm;
  _parton=Parton;
  _mP2=_parton->Mass()*_parton->Mass();
}

InclusiveJetLO::~InclusiveJetLO(){
  delete _qcd;
  delete _parton;
}


//####################################################################################################
//QCD parameters setting
//Quark flavors
void InclusiveJetLO::setNf(const int nf){
  _nf=nf;
  _qcd->setNfQCD(nf);
}

//LambdaQCD for n-loop β-function coefficient
void InclusiveJetLO::setLambdaQCD(const int nloop){
  _lambdaQCD=_qcd->LambdaQCD(nloop);
  // std::cout << "lambdaQCD=" << _lambdaQCD << std::endl;
}


//####################################################################################################
//Kinematics calculations in 248-250 pages of [QCD and Collider Physics] and J. F. Owens paper for [Large Momentum Transfer Production of Direct Photons, Jets, and Particles]
//Calculate Bjorken x's from Ecm, pT, and the rapidities of the two produced particles
void InclusiveJetLO::MomentumFractions(const double y3c, const double y4d){
  //Massless partons
  double xT=2.*_pT/std::sqrt(_s);
  _x1=(1./2.)*xT*(std::exp( y3c)+std::exp( y4d));
  _x2=(1./2.)*xT*(std::exp(-y3c)+std::exp(-y4d));
  // //Massive partons
  // double mT3=std::sqrt(_mP2+_pT*_pT);   //triggered parton or gamma
  // double mJ=0.;                         //If jet is a q/g/gamma jet, jet mass is equal to 0 (mJ=0.), except for the heavy quark or Z/W parton.
  // double mT4=std::sqrt(mJ*mJ+_pT*_pT);  //observed jet
  // _x1=(1./std::sqrt(_s))*(mT3*std::exp( y3c)+mT4*std::exp( y4d));
  // _x2=(1./std::sqrt(_s))*(mT3*std::exp(-y3c)+mT4*std::exp(-y4d));
}

//Mandelstam variables
void InclusiveJetLO::Mandelstam(const double y3c, const double y4d){
  //Satisfy s+t+u=m1^2+m2^2+m3^2+m4^2=0 for massless partons
  _shat= _x1*_x2*_s;
  _that=-_x1*_pT*std::sqrt(_s)*std::exp(-y3c);  //equal to _that=-_x2*_pT*std::sqrt(_s)*std::exp(-y4d);
  _uhat=-_x2*_pT*std::sqrt(_s)*std::exp( y3c);  //equal to _uhat=-_x1*_pT*std::sqrt(_s)*std::exp(-y4d);
  // std::cout << "pT=" << _pT << std::endl;
  // getKinematics();
  // _that=(-1./2.)*_shat*(1.-std::tanh((y3c-y4d)/2.));
  // _uhat=(-1./2.)*_shat*(1.+std::tanh((y3c-y4d)/2.));
  // _that=-_shat*1./(std::exp(y3c-y4d)+1.);
  // _uhat=-_shat*1./(std::exp(-(y3c-y4d))+1.);
  // getKinematics();
}

void InclusiveJetLO::setKinematics(const double y3c, const double y4d){
  _y3c=y3c;
  _y4d=y4d;
  MomentumFractions(y3c,y4d);
  Mandelstam(y3c,y4d);
}

void InclusiveJetLO::getKinematics(){
  std::cout << std::setprecision(16) << "x1=" << _x1 << ", x2=" << _x2 << ", s=" << _shat << ", t=" << _that << ", u=" << _uhat << ", s+t+u=" << _shat+_that+_uhat << std::endl;
}

//Four momenta of the incoming partons and outgoing particles
void InclusiveJetLO::getFourMomenta(){
  //mT=sqrt(_parton->Mass()*_parton->Mass()+_pT*_pT),
  //p1=_x1*(_Ecm/2., 0, 0,  _Ecm/2.)
  //p2=_x2*(_Ecm/2., 0, 0, -_Ecm/2.)
  //p3=(mT*cosh(y3c), pT*cos(phi3), pT*sin(phi3), mT*sinh(y3c))
  //p4=(mT*cosh(y4d), pT*cos(phi4), pT*sin(phi4), mT*sinh(y4d))
  std::cout << 0.5*_x1*_Ecm << " " << 0.0 << " " << 0.0 << " " << -0.5*_x1*_Ecm << std::endl;
  std::cout << 0.5*_x2*_Ecm << " " << 0.0 << " " << 0.0 << " " << -0.5*_x2*_Ecm << std::endl;
  std::cout << std::sqrt(_pT*_pT+_mP2)*std::cosh(_y3c) << " " << 0.0 << " " << _pT << " "  << std::sqrt(_pT*_pT+_mP2)*std::sinh(_y3c) << std::endl;
  std::cout << _pT*std::cosh(_y4d) << " " << 0.0 << " " << -_pT <<  " " << _pT*std::sinh(_y4d) << std::endl;
 }


//####################################################################################################
//PDFs: parton distribution functions
void InclusiveJetLO::setPDF1(){
  double Q=_pT*_pTscale;
  if(KinematicsQcut()){
    for(int i=1; i<=_nf; i++){
      _pdf1aplus[i] =_qcd->PDFs( i,_x1,Q);
      _pdf1aminus[i]=_qcd->PDFs(-i,_x1,Q);
    }
    _pdf1aplus[21]=_qcd->PDFs(21,_x1,Q);
  }
}
void InclusiveJetLO::setPDF2(){
  double Q=_pT*_pTscale;
  if(KinematicsQcut()){
    for(int i=1; i<=_nf; i++){
      _pdf2bplus[i] =_qcd->PDFs( i,_x2,Q);
      _pdf2bminus[i]=_qcd->PDFs(-i,_x2,Q);
    }
    _pdf2bplus[21]=_qcd->PDFs(21,_x2,Q);
  }
}

double InclusiveJetLO::getPDF1(const int i){
  double  pdfval=0.0;
  if(i>0) pdfval=_pdf1aplus[i];
  if(i<0) pdfval=_pdf1aminus[-i];
  return  pdfval;
}
double InclusiveJetLO::getPDF2(const int i){
  double  pdfval=0.0;
  if(i>0) pdfval=_pdf2bplus[i];
  if(i<0) pdfval=_pdf2bminus[-i];
  return  pdfval;
}


//####################################################################################################
//JFFs: jet fragmentation functions
void InclusiveJetLO::setJFF3(){
  _Dz3c[21]=1.0;
  for(int i=1; i<=_nf; i++){
    _Dz3c[i]=1.0;
  }
}
void InclusiveJetLO::setJFF4(){
  _Dz4d[21]=1.0;
  for(int i=1; i<=_nf; i++){
    _Dz4d[i]=1.0;
  }
}

double InclusiveJetLO::getJFF3(const int i){
  return 1.0;
  // return _Dz3c[std::abs(i)];
}
double InclusiveJetLO::getJFF4(const int i){
  return 1.0;
  // return _Dz4d[std::abs(i)];
}


//####################################################################################################
//In 248 and 250 pages of [QCD and Collider Physics] and J. F. Owens paper for [Large Momentum Transfer Production of Direct Photons, Jets, and Particles]
//Function prefactor for "dsigma/dpT=2*pi*pT*alphas^2/(x1*x2*s^2)*M2*PDF1(x2)*PDF2(x2)" without PDFs and M2 at LO
void InclusiveJetLO::setXsection(){
  double Q=_pT*_pTscale;
  _preXsection=2.*M_PI*_pT*_qcd->AlphaS(Q)*_qcd->AlphaS(Q)/(_x1*_x2*_s*_s);
}


//####################################################################################################
//See Table 7.1 in 249 pages of [QCD and Collider Physics] and Table 1 in 470 pages of J. F. Owens paper for [Large Momentum Transfer Production of Direct Photons, Jets, and Particles]
//Amplitudes squared: M^2/(gs^4) averaged and summed over the spin and color indices repectively in the intial and final states
inline double InclusiveJetLO::M2qqp2qqp(const double &s, const double &t, const double &u){
  //qq'->qq' and qqbar'->qqbar'
  return (4./9.)*((s*s+u*u)/(t*t));
}

inline double InclusiveJetLO::M2qq2qq(const double &s, const double &t, const double &u){
  //qq->qq
  return (4./9.)*((s*s+u*u)/(t*t)+(s*s+t*t)/(u*u))-(8./27.)*((s*s)/(u*t));
}

inline double InclusiveJetLO::M2qqb2qpqpb(const double &s, const double &t, const double &u){
  //qqbar->q'q'bar
  return (4./9.)*((t*t+u*u)/(s*s));
}

inline double InclusiveJetLO::M2qqb2qqb(const double &s, const double &t, const double &u){
  //qqbar->qqbar
  return (4./9.)*((s*s+u*u)/(t*t)+(t*t+u*u)/(s*s))-(8./27.)*((u*u)/(s*t));
}

inline double InclusiveJetLO::M2qqb2gg(const double &s, const double &t, const double &u){
  //qqbar->gg
  return (32./27.)*((t*t+u*u)/(t*u))-(8./3.)*((t*t+u*u)/(s*s));
}

inline double InclusiveJetLO::M2gg2qqb(const double &s, const double &t, const double &u){
  //gg->qqbar
  return (1./6.)*((t*t+u*u)/(t*u))-(3./8.)*((t*t+u*u)/(s*s));
}

inline double InclusiveJetLO::M2gq2gq(const double &s, const double &t, const double &u){
  //gq->gq
  return -(4./9.)*((s*s+u*u)/(s*u))+(u*u+s*s)/(t*t);
}

inline double InclusiveJetLO::M2gg2gg(const double &s, const double &t, const double &u){
  //gg->gg
  return (9./2.)*(3.-(t*u)/(s*s)-(s*u)/(t*t)-(s*t)/(u*u));
}


//####################################################################################################
//Differential cross section of inclusive jet at LO (2->2)
//dsigma/dpT=2*π*pT*⍺s^2/(x1*x2*s^2)*M2*pdf1(x1)*pdf2(x2) in pb/GeV
//qq'->qq' and qq'->q'q (t,u channel)
double InclusiveJetLO::SigmaLOqqp2qqp(){
  double res=0.0;
  if(KinematicsQcut()){
    double dsi3j4=0.0;
    double dsj4i3=0.0;
    for(int i=(-_nf); i<=_nf; i++){
      if(i==0) continue;
      for(int j=(-_nf); j<=_nf; j++){
        if(j==0) continue;
        if(std::abs(i)==std::abs(j)) continue;
        //qq'->q3q'4
        dsi3j4+=getPDF1(i)*getPDF2(j)*getJFF3(i)*getJFF4(j);
        //qq'->q'4q3, exchange final partons (t->u,u->t)
        dsj4i3+=getPDF1(i)*getPDF2(j)*getJFF3(j)*getJFF4(i);
      }
    }
    res=(dsi3j4*M2qqp2qqp(_shat,_that,_uhat)+dsj4i3*M2qqp2qqp(_shat,_uhat,_that))*getXsection();
  }
  return res*_unitGeV2pb;  //return with unit pb
}

//qq->qq (t,u channel)
double InclusiveJetLO::SigmaLOqq2qq(){
  double res=0.0;
  if(KinematicsQcut()){
    double dsi3i4=0.0;
    for(int i=(-_nf); i<=_nf; i++){
      if(i==0) continue;
      //qq->q3q4
      dsi3i4+=getPDF1(i)*getPDF2(i)*getJFF3(i)*getJFF4(i)*(1./2.);  //1./2. from identical particles in the final state
    }
    dsi3i4*=2.0;  //2.0 for q3 and q4 inclusive jet in the final state
    res=dsi3i4*M2qq2qq(_shat,_that,_uhat)*getXsection();
  }
  return res*_unitGeV2pb;  //return with unit pb
}

//qqbar->q'q'bar (s channel)
double InclusiveJetLO::SigmaLOqqb2qpqpb(){
  double res=0.0;
  if(KinematicsQcut()){
    double dsj3j4=0.0;
    for(int i=(-_nf); i<=_nf; i++){
      if(i==0) continue;
      for(int j=(-_nf); j<=_nf; j++){
        if(j==0) continue;
        if(std::abs(i)==std::abs(j)) continue;
        //qqbar->q'3q'bar4
        dsj3j4+=getPDF1(i)*getPDF2(-i)*getJFF3(j)*getJFF4(-j);
      }
    }
    res=dsj3j4*M2qqb2qpqpb(_shat,_that,_uhat)*getXsection();
  }
  return res*_unitGeV2pb;  //return with unit pb
}

//qqbar->qqbar and qqbar->qbarq (s,t,u channel)
double InclusiveJetLO::SigmaLOqqb2qqb(){
  double res=0.0;
  if(KinematicsQcut()){
    double dsi3i4=0.0;
    double dsi4i3=0.0;
    for(int i=(-_nf); i<=_nf; i++){
      if(i==0) continue;
      //qqbar->q3qbar4
      dsi3i4+=getPDF1(i)*getPDF2(-i)*getJFF3( i)*getJFF4(-i);
      //qqbar->qbar3q4, exchange final partons (t->u,u->t)
      dsi4i3+=getPDF1(i)*getPDF2(-i)*getJFF3(-i)*getJFF4( i);
    }
    res=(dsi3i4*M2qqb2qqb(_shat,_that,_uhat)+dsi4i3*M2qqb2qqb(_shat,_uhat,_that))*getXsection();
  }
  return res*_unitGeV2pb;  //return with unit pb
}

//qqbar->gg (s,t,u channel)
double InclusiveJetLO::SigmaLOqqb2gg(){
  double res=0.0;
  if(KinematicsQcut()){
    double dsg3g4=0.0;
    for(int i=(-_nf); i<=_nf; i++){
      if(i==0) continue;
      //qqbar->g3g4
      dsg3g4+=getPDF1(i)*getPDF2(-i)*getJFF3(21)*getJFF4(21)*(1./2.);  //1./2. from identical particles in the final state
    }
    dsg3g4*=2.0;  //2.0 for g3 and g4 inclusive jet in the final state
    res=dsg3g4*M2qqb2gg(_shat,_that,_uhat)*getXsection();
  }
  return res*_unitGeV2pb;  //return with unit pb
}

//gg->qqbar (s,t,u channel)
double InclusiveJetLO::SigmaLOgg2qqb(){
  double res=0.0;
  if(KinematicsQcut()){
    double dsj3j4=0.0;
    for(int j=(-_nf); j<=_nf; j++){
      if(j==0) continue;
      //gg->qqbar
      dsj3j4+=getPDF1(21)*getPDF2(21)*getJFF3(j)*getJFF4(-j);
    }
    res=dsj3j4*M2gg2qqb(_shat,_that,_uhat)*getXsection();
  }
  return res*_unitGeV2pb;  //return with unit pb
}

//gq->gq and qg->gq, qg->qg and gq->qg (s,t,u channel)
double InclusiveJetLO::SigmaLOgq2gq(){
  double res=0.0;
  if(KinematicsQcut()){
    double ds =0.0;
    double ds_=0.0;
    for(int i=(-_nf); i<=_nf; i++){
      if(i==0) continue;
      //gq->g3q4 and qg->q3g4
      ds +=getPDF1(21)*getPDF2( i)*getJFF3(21)*getJFF4( i);
      ds +=getPDF1( i)*getPDF2(21)*getJFF3( i)*getJFF4(21);
      //qg->g3q4 and gq->q3g4, exchange final partons (t->u,u->t)
      ds_+=getPDF1( i)*getPDF2(21)*getJFF3(21)*getJFF4( i);
      ds_+=getPDF1(21)*getPDF2( i)*getJFF3( i)*getJFF4(21);
    }
    res=(ds*M2gq2gq(_shat,_that,_uhat)+ds_*M2gq2gq(_shat,_uhat,_that))*getXsection();
  }
  return res*_unitGeV2pb;  //return with unit pb
}
//Gluon jet (3): observed jet
//gq->gq and qg->gq (s,t,u channel)
double InclusiveJetLO::SigmaLOgq2gq_G(){
  double res=0.0;
  if(KinematicsQcut()){
    double dsg =0.0;
    double dsg_=0.0;
    for(int i=(-_nf); i<=_nf; i++){
      if(i==0) continue;
      //gq->g3q4
      dsg +=getPDF1(21)*getPDF2( i)*getJFF3(21)*getJFF4( i);
      //qg->g3q4, exchange final partons (t->u,u->t)
      dsg_+=getPDF1( i)*getPDF2(21)*getJFF3(21)*getJFF4( i);
    }
    res=(dsg*M2gq2gq(_shat,_that,_uhat)+dsg_*M2gq2gq(_shat,_uhat,_that))*getXsection();
  }
  return res*_unitGeV2pb;  //return with unit pb
}
//Quark jet (3): observed jet
//qg->qg and gq->qg (s,t,u channel)
double InclusiveJetLO::SigmaLOgq2gq_Q(){
  double res=0.0;
  if(KinematicsQcut()){
    double dsq =0.0;
    double dsq_=0.0;
    for(int i=(-_nf); i<=_nf; i++){
      if(i==0) continue;
      //qg->q3g4
      dsq +=getPDF1( i)*getPDF2(21)*getJFF3( i)*getJFF4(21);
      //gq->q3g4, exchange final partons (t->u,u->t)
      dsq_+=getPDF1(21)*getPDF2( i)*getJFF3( i)*getJFF4(21);
    }
    res=(dsq*M2gq2gq(_shat,_that,_uhat)+dsq_*M2gq2gq(_shat,_uhat,_that))*getXsection();
  }
  return res*_unitGeV2pb;  //return with unit pb
}

//gg->gg (s,t,u channel)
double InclusiveJetLO::SigmaLOgg2gg(){
  double res=0.0;
  if(KinematicsQcut()){
    double dsg3g4=0.0;
    //gg->g3g4
    dsg3g4+=getPDF1(21)*getPDF2(21)*getJFF3(21)*getJFF4(21)*(1./2.);  //1./2. from identical particles in the final state
    dsg3g4*=2.0;  //2.0 for g3 and g4 inclusive jet in the final state

    res=dsg3g4*M2gg2gg(_shat,_that,_uhat)*getXsection();
  }
  return res*_unitGeV2pb;  //return with unit pb
}


//####################################################################################################
//Leading order kinematic parameters
void InclusiveJetLO::setParametersLO(const double y3c, const double y4d){
  setKinematics(y3c,y4d);
  setPDF1(); setPDF2();
  setJFF3(); setJFF4();
  setXsection();
}

//Leading order cross section
//All jets
double InclusiveJetLO::getSigmaLO(){
  double res=0.0;
  if(KinematicsQcut()){
    res=SigmaLOqqp2qqp()+SigmaLOqq2qq()+SigmaLOqqb2qpqpb()+SigmaLOqqb2qqb()+SigmaLOqqb2gg()+SigmaLOgq2gq()+SigmaLOgg2qqb()+SigmaLOgg2gg();
  }
  return res;
}

//Quark jets
double InclusiveJetLO::getSigmaLOQuark(){
  double res=0.0;
  if(KinematicsQcut()){
    res=SigmaLOqqp2qqp()+SigmaLOqq2qq()+SigmaLOqqb2qpqpb()+SigmaLOqqb2qqb()+SigmaLOgq2gq_Q()+SigmaLOgg2qqb();
  }
  return res;
}

//Gluon jets
double InclusiveJetLO::getSigmaLOGluon(){
  double res=0.0;
  if(KinematicsQcut()){
    res=SigmaLOqqb2gg()+SigmaLOgq2gq_G()+SigmaLOgg2gg();
  }
  return res;
}


//####################################################################################################
//MC integration
void InclusiveJetLO::setCalls(const size_t calls){
  _calls=calls;
}

void InclusiveJetLO::setupMC(double (*func)(double *, size_t, void *), size_t dim, double *rangemin, double *rangemax, void *para){
  gsl_rng_env_setup();
  _Type=gsl_rng_default;
  _rng=gsl_rng_alloc(_Type);

  _dim=dim;
  _rangemin=rangemin; _rangemax=rangemax;
  _mcFun.f=func; _mcFun.dim=dim; _mcFun.params=para;

#ifdef MC_MISER
  _mcStatus=gsl_monte_miser_alloc(dim);
  std::cout << "Miser is used ..." << std::endl;
#else
  _mcStatus=gsl_monte_vegas_alloc(dim);
  std::cout << "Vegas is used ..." << std::endl;
#endif
}

void InclusiveJetLO::calculateMC(double &res, double &err){
#ifdef MC_MISER
  gsl_monte_miser_integrate(&_mcFun, _rangemin, _rangemax, _dim, _calls, _rng, _mcStatus, &res, &err);
#else
  gsl_monte_vegas_integrate(&_mcFun, _rangemin, _rangemax, _dim, _calls, _rng, _mcStatus, &res, &err);
#endif
}

void InclusiveJetLO::freeMC(){
  gsl_rng_free(_rng);
#ifdef MC_MISER
  gsl_monte_miser_free(_mcStatus);
#else
  gsl_monte_vegas_free(_mcStatus);
#endif
}
