#include "DiJetLO.h"


DiJetLO::DiJetLO(){
  _qcd = new QCD();
  _nf=_qcd->Nf();
  _lambdaQCD=_qcd->LambdaQCD(1);
  _pTscale=1.;
  _unitGeV2pb=3.893793721e8;  //1 (ℏ𝑐)^2 = 3.893793721e8 GeV^2 pb from PDG 2024
  _calls=50000;
  _Ecm=5.02e3;  //GeV
  _s=_Ecm*_Ecm;
  _parton = new Particle();
  _mP2=0.0;
}

DiJetLO::DiJetLO(const double Ecm, Particle *Parton){
  _qcd = new QCD();
  _nf=_qcd->Nf();
  _lambdaQCD=_qcd->LambdaQCD(1);
  _pTscale=1.;
  _unitGeV2pb=3.893793721e8;
  _calls=50000;
  _Ecm=Ecm;
  _s=_Ecm*_Ecm;
  _parton=Parton;
  _mP2=_parton->Mass()*_parton->Mass();
}

DiJetLO::~DiJetLO(){
  delete _qcd;
  delete _parton;
}


//####################################################################################################
//QCD
//Quark flavor
void DiJetLO::setNf(const unsigned int nf){
  _nf=nf;
  _qcd->setNf(nf);
}

//LambdaQCD for loop
void DiJetLO::setLambdaQCD(const unsigned int nloop){
  _lambdaQCD=_qcd->LambdaQCD(nloop);  //nloop=1 for LO
  // std::cout << "lambdaQCD=" << _lambdaQCD << std::endl;
}


//####################################################################################################
//Kinematics calculation in 248-250 pages of [QCD and Collider Physics]
//Calculate Bjorken x's from Ecm, pT and the rapidities of the two produced particles
void DiJetLO::MomentumFractions(const double y3c, const double y4d){
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
void DiJetLO::Mandelstam(const double y3c, const double y4d){
  //satisfy s+t+u=m1^2+m2^2+m3^2+m4^2=0 for massless partons
  _shat= _x1*_x2*_s;
  _that=-_x1*_Ecm*_pT*std::exp(-y3c);
  _uhat=-_x1*_Ecm*_pT*std::exp(-y4d);  //equal to _uhat=-_x2*_Ecm*_pT*std::exp( y3c);
  // std::cout << "pT=" << _pT << std::endl;
  // getKinematics();
  // _that=(-1./2.)*_shat*(1.-std::tanh((y3c-y4d)/2.));
  // _uhat=(-1./2.)*_shat*(1.+std::tanh((y3c-y4d)/2.));
  // _that=-_shat*1./(std::exp(y3c-y4d)+1.);
  // _uhat=-_shat*1./(std::exp(-(y3c-y4d))+1.);
  // getKinematics();
}

void DiJetLO::setKinematics(const double y3c, const double y4d){
  _y3c=y3c;
  _y4d=y4d;
  MomentumFractions(y3c,y4d);
  Mandelstam(y3c,y4d);
}

void DiJetLO::getKinematics(){
  std::cout << std::setprecision(16) << "x1=" << _x1 << ", x2=" << _x2 << ", s=" << _shat << ", t=" << _that << ", u=" << _uhat << ", s+t+u=" << _shat+_that+_uhat << std::endl;
}

//Four momenta of the incoming partons and outgoing particles
void DiJetLO::getFourMomenta(){
  //Four-momentum of incoming and outgoing particles
  //mT=sqrt(_parton->Mass()*_parton->Mass()+_pT*_pT),
  //p1=_x1*(_Ecm/2., 0, 0,  _Ecm/2.),
  //p2=_x2*(_Ecm/2., 0, 0, -_Ecm/2.),
  //p3=(mT*cosh(y3c), pT*cos(phi3), pT*sin(phi3), mT*sinh(y3c)),
  //p4=(mT*cosh(y4d), pT*cos(phi4), pT*sin(phi4), mT*sinh(y4d)).
  std::cout << 0.5*_x1*_Ecm << " " << 0.0 << " " << 0.0 << " " << -0.5*_x1*_Ecm << std::endl;
  std::cout << 0.5*_x2*_Ecm << " " << 0.0 << " " << 0.0 << " " << -0.5*_x2*_Ecm << std::endl;
  std::cout << std::sqrt(_pT*_pT+_mP2)*std::cosh(_y3c) << " " << 0.0 << " " << _pT << " "  << std::sqrt(_pT*_pT+_mP2)*std::sinh(_y3c) << std::endl;
  std::cout << _pT*std::cosh(_y4d) << " " << 0.0 << " " << -_pT <<  " " << _pT*std::sinh(_y4d) << std::endl;
 }


//####################################################################################################
//PDF: parton distribution function
void DiJetLO::setPDF1(){
  if(KinematicsQcut()){
    for(unsigned int i=1; i<=_nf; i++){
      _pdf1aplus[i] =_qcd->pdf( i,_x1,_pT*_pTscale);
      _pdf1aminus[i]=_qcd->pdf(-i,_x1,_pT*_pTscale);
    }
    _pdf1aplus[21]=_qcd->pdf(21,_x1,_pT*_pTscale);
  }
}
void DiJetLO::setPDF2(){
  if(KinematicsQcut()){
    for(unsigned int i=1; i<=_nf; i++){
      _pdf2bplus[i] =_qcd->pdf( i,_x2,_pT*_pTscale);
      _pdf2bminus[i]=_qcd->pdf(-i,_x2,_pT*_pTscale);
    }
    _pdf2bplus[21]=_qcd->pdf(21,_x2,_pT*_pTscale);
  }
}

double DiJetLO::getPDF1(const int i){
  double  pdfval=0.0;
  if(i>0) pdfval=_pdf1aplus[i];
  if(i<0) pdfval=_pdf1aminus[-i];
  return  pdfval;
}
double DiJetLO::getPDF2(const int i){
  double  pdfval=0.0;
  if(i>0) pdfval=_pdf2bplus[i];
  if(i<0) pdfval=_pdf2bminus[-i];
  return  pdfval;
}


//####################################################################################################
//JFF: Jet fragmentation function
void DiJetLO::setJetFF3(){
  _Dz3c[21]=1.0;
  for(unsigned int i=1; i<=_nf; i++){
    _Dz3c[i]=1.0;
  }
}
void DiJetLO::setJetFF4(){
  _Dz4d[21]=1.0;
  for(unsigned int i=1; i<=_nf; i++){
    _Dz4d[i]=1.0;
  }
}

double DiJetLO::getJetFF3(const int i){
  return 1.0;
  // return _Dz3c[std::abs(i)];
}
double DiJetLO::getJetFF4(const int i){
  return 1.0;
  // return _Dz4d[std::abs(i)];
}


//####################################################################################################
//In 248 and 250 pages of [QCD and Collider Physics]
//Function prefactor for "dsigma=2*pi*pT*alphas^2/(x1*x2*s^2)*M2*pdf1(x2)*pdf2(x2)" without pdfs and M2
double DiJetLO::alphasPDF(const double pT){
    return _qcd->alphas(pT*_pTscale);
}

void DiJetLO::setXsection(){
  _preXsection=2.*M_PI*_pT*alphasPDF(_pT)*alphasPDF(_pT)/(_x1*_x2*_s*_s);  //DiJet LO
}


//####################################################################################################
//See Table 7.1 in 249 pages of [QCD and Collider Physics]
//Amplitudes squared: M^2/(gs^4) averaged and summed over the spin and color indices repectively in the intial and final states
inline double DiJetLO::M2qqp2qqp(const double &s, const double &t, const double &u){
  //qq'->qq' and qqbar'->qqbar'
  return (4./9.)*((s*s+u*u)/(t*t));
}

inline double DiJetLO::M2qq2qq(const double &s, const double &t, const double &u){
  //qq->qq
  return (4./9.)*((s*s+u*u)/(t*t)+(s*s+t*t)/(u*u))-(8./27.)*(s*s/(u*t));
}

inline double DiJetLO::M2qqb2qpqpb(const double &s, const double &t, const double &u){
  //qqbar->q'q'bar
  return (4./9.)*((t*t+u*u)/(s*s));
}

inline double DiJetLO::M2qqb2qqb(const double &s, const double &t, const double &u){
  //qqbar->qqbar
  return (4./9.)*((s*s+u*u)/(t*t)+(t*t+u*u)/(s*s))-(8./27.)*(u*u/(s*t));
}

inline double DiJetLO::M2qqb2gg(const double &s, const double &t, const double &u){
  //qqbar->gg
  return (32./27.)*((t*t+u*u)/(t*u))-(8./3.)*((t*t+u*u)/(s*s));
}

inline double DiJetLO::M2gg2qqb(const double &s, const double &t, const double &u){
  //gg->qqbar
  return (1./6.)*((t*t+u*u)/(t*u))-(3./8.)*((t*t+u*u)/(s*s));
}

inline double DiJetLO::M2gq2gq(const double &s, const double &t, const double &u){
  //gq->gq
  return -(4./9.)*((s*s+u*u)/(s*u))+(u*u+s*s)/(t*t);
}

inline double DiJetLO::M2gg2gg(const double &s, const double &t, const double &u){
  //gg->gg
  return (9./2.)*(3.-(t*u)/(s*s)-(s*u)/(t*t)-(s*t)/(u*u));
}


//####################################################################################################
//Differential cross sectoin of Dijet at LO (2->2)
//dsigma=2*pi*pT*alphas^2/(x1*x2*s^2)*M2*pdf1(x2)*pdf2(x2) in pb/GeV
//qq'->qq' and qq'->q'q
double DiJetLO::SigmaLOqqp2qqp(){
  double res=0.0;
  if(KinematicsQcut()){
    double dsqqp =0.0;
    double dsqqp_=0.0;
    for(unsigned int i=1; i<=_nf; i++){
      for(unsigned int j=1; j<=_nf; j++){
        if(i==j) continue;
        //qq'->qq'
        dsqqp +=getPDF1( i)*getPDF2( j)*getJetFF3( i)*getJetFF4( j);
        dsqqp +=getPDF1( i)*getPDF2(-j)*getJetFF3( i)*getJetFF4(-j);
        dsqqp +=getPDF1(-i)*getPDF2( j)*getJetFF3(-i)*getJetFF4( j);
        dsqqp +=getPDF1(-i)*getPDF2(-j)*getJetFF3(-i)*getJetFF4(-j);
        //qq'->q'q, exchange parton (t->u,u->t)
        dsqqp_+=getPDF1( i)*getPDF2( j)*getJetFF3( j)*getJetFF4( i);
        dsqqp_+=getPDF1( i)*getPDF2(-j)*getJetFF3(-j)*getJetFF4( i);
        dsqqp_+=getPDF1(-i)*getPDF2( j)*getJetFF3( j)*getJetFF4(-i);
        dsqqp_+=getPDF1(-i)*getPDF2(-j)*getJetFF3(-j)*getJetFF4(-i);
      }
    }
    res=(dsqqp*M2qqp2qqp(_shat,_that,_uhat)+dsqqp_*M2qqp2qqp(_shat,_uhat,_that))*getXsection();
    // std::cout << "res1=" << res << std::endl;
  }
  return res*_unitGeV2pb;  //return with unit pb
}

//qq->qq
double DiJetLO::SigmaLOqq2qq(){
  double res=0.0;
  if(KinematicsQcut()){
    double dsqq=0.0;
    for(unsigned int i=1; i<=_nf; i++){
      //qq->qq
      dsqq+=getPDF1( i)*getPDF2( i)*getJetFF3( i)*getJetFF4( i);
      dsqq+=getPDF1(-i)*getPDF2(-i)*getJetFF3(-i)*getJetFF4(-i);
    }
    res=dsqq*M2qq2qq(_shat,_that,_uhat)*getXsection();
    // std::cout << "res2=" << res << std::endl;
  }
  return res*_unitGeV2pb;  //return with unit pb
}

//qqbar->q'q'bar
double DiJetLO::SigmaLOqqb2qpqpb(){
  double res=0.0;
  if(KinematicsQcut()){
    double dsqpqpb=0.0;
    for(unsigned int i=1; i<=_nf; i++){
      for(unsigned int j=1; j<=_nf; j++){
        if(i==j) continue;
        //qqbar->q'q'bar
        dsqpqpb+=getPDF1( i)*getPDF2(-i)*getJetFF3( j)*getJetFF4(-j);
        dsqpqpb+=getPDF1( i)*getPDF2(-i)*getJetFF3(-j)*getJetFF4( j);
        dsqpqpb+=getPDF1(-i)*getPDF2( i)*getJetFF3( j)*getJetFF4(-j);
        dsqpqpb+=getPDF1(-i)*getPDF2( i)*getJetFF3(-j)*getJetFF4( j);
      }
    }
    res=dsqpqpb*M2qqb2qpqpb(_shat,_that,_uhat)*getXsection();
    // std::cout << "res3=" << res << std::endl;
  }
  return res*_unitGeV2pb;  //return with unit pb
}

//qqbar->qqbar and qqbar->qbarq
double DiJetLO::SigmaLOqqb2qqb(){
  double res=0.0;
  if(KinematicsQcut()){
    double dsqqb =0.0;
    double dsqqb_=0.0;
    for(unsigned int i=1; i<=_nf; i++){
      //qqbar->qqbar
      dsqqb +=getPDF1( i)*getPDF2(-i)*getJetFF3( i)*getJetFF4(-i);
      dsqqb +=getPDF1(-i)*getPDF2( i)*getJetFF3(-i)*getJetFF4( i);
      //qqbar->qbarq, exchange parton (t->u,u->t)
      dsqqb_+=getPDF1( i)*getPDF2(-i)*getJetFF3(-i)*getJetFF4( i);
      dsqqb_+=getPDF1(-i)*getPDF2( i)*getJetFF3( i)*getJetFF4(-i);
    }
    res=(dsqqb*M2qqb2qqb(_shat,_that,_uhat)+dsqqb_*M2qqb2qqb(_shat,_uhat,_that))*getXsection();
    // std::cout << "res4=" << res << std::endl;
  }
  return res*_unitGeV2pb;  //return with unit pb
}

//qqbar->gg
double DiJetLO::SigmaLOqqb2gg(){
  double res=0.0;
  if(KinematicsQcut()){
    double dsqqbg=0.0;
    for(unsigned int i=1; i<=_nf; i++){
      //qqbar->gg
      dsqqbg+=getPDF1( i)*getPDF2(-i)*getJetFF3(21)*getJetFF4(21);
      dsqqbg+=getPDF1(-i)*getPDF2( i)*getJetFF3(21)*getJetFF4(21);
    }
    res=dsqqbg*M2qqb2gg(_shat,_that,_uhat)*getXsection();
    // std::cout << "res5=" << res << std::endl;
  }
  return res*_unitGeV2pb;  //return with unit pb
}

//gg->qqbar
double DiJetLO::SigmaLOgg2qqb(){
  double res=0.0;
  if(KinematicsQcut()){
    double dsggqqb=0.0;
    for(unsigned int j=1; j<=_nf; j++){
      //gg->qqbar
      dsggqqb+=getPDF1(21)*getPDF2(21)*getJetFF3( j)*getJetFF4(-j);
      dsggqqb+=getPDF1(21)*getPDF2(21)*getJetFF3(-j)*getJetFF4( j);
    }
    res=dsggqqb*M2gg2qqb(_shat,_that,_uhat)*getXsection();
    // std::cout << "res6=" << res << std::endl;
  }
  return res*_unitGeV2pb;  //return with unit pb
}

//gq(qg)->gq(qg) and gq(qg)->qg(gq)
double DiJetLO::SigmaLOgq2gq(){
  double res=0.0;
  if(KinematicsQcut()){
    double dsgq =0.0;
    double dsgq_=0.0;
    for(unsigned int i=1; i<=_nf; i++){
      //gq(qg)->gq(qg)
      dsgq +=getPDF1(21)*getPDF2( i)*getJetFF3(21)*getJetFF4( i);
      dsgq +=getPDF1(21)*getPDF2(-i)*getJetFF3(21)*getJetFF4(-i);
      dsgq +=getPDF1( i)*getPDF2(21)*getJetFF3( i)*getJetFF4(21);
      dsgq +=getPDF1(-i)*getPDF2(21)*getJetFF3(-i)*getJetFF4(21);
      //gq(qg)->qg(gq), exchange parton (t->u,u->t)
      dsgq_+=getPDF1(21)*getPDF2( i)*getJetFF3( i)*getJetFF4(21);
      dsgq_+=getPDF1(21)*getPDF2(-i)*getJetFF3(-i)*getJetFF4(21);
      dsgq_+=getPDF1( i)*getPDF2(21)*getJetFF3(21)*getJetFF4( i);
      dsgq_+=getPDF1(-i)*getPDF2(21)*getJetFF3(21)*getJetFF4(-i);
    }
    res=(dsgq*M2gq2gq(_shat,_that,_uhat)+dsgq_*M2gq2gq(_shat,_uhat,_that))*getXsection();
    // std::cout << "res7=" << res << std::endl;
  }
  return res*_unitGeV2pb;  //return with unit pb
}

//Quark jet (3): observed jet
//qg->qg and gq->qg
double DiJetLO::SigmaLOgq2gq_Q(){
  double res=0.0;
  if(KinematicsQcut()){
    double dsgq1 =0.0;
    double dsgq1_=0.0;
    for(unsigned int i=1; i<=_nf; i++){
      //qg->qg
      dsgq1 +=getPDF1( i)*getPDF2(21)*getJetFF3( i)*getJetFF4(21);
      dsgq1 +=getPDF1(-i)*getPDF2(21)*getJetFF3(-i)*getJetFF4(21);
      //gq->qg, exchange parton (t->u,u->t)
      dsgq1_+=getPDF1(21)*getPDF2( i)*getJetFF3( i)*getJetFF4(21);
      dsgq1_+=getPDF1(21)*getPDF2(-i)*getJetFF3(-i)*getJetFF4(21);
    }
    res=(dsgq1*M2gq2gq(_shat,_that,_uhat)+dsgq1_*M2gq2gq(_shat,_uhat,_that))*getXsection();
    // std::cout << "Quark: res7=" << res << std::endl;
  }
  return res*_unitGeV2pb;  //return with unit pb
}
//Gluon jet (3): observed jet
//gq->gq and qg->gq
double DiJetLO::SigmaLOgq2gq_G(){
  double res=0.0;
  if(KinematicsQcut()){
    double dsgq2 =0.0;
    double dsgq2_=0.0;
    for(unsigned int i=1; i<=_nf; i++){
      //gq->gq
      dsgq2 +=getPDF1(21)*getPDF2( i)*getJetFF3(21)*getJetFF4( i);
      dsgq2 +=getPDF1(21)*getPDF2(-i)*getJetFF3(21)*getJetFF4(-i);
      //qg->gq, exchange parton (t->u,u->t)
      dsgq2_+=getPDF1( i)*getPDF2(21)*getJetFF3(21)*getJetFF4( i);
      dsgq2_+=getPDF1(-i)*getPDF2(21)*getJetFF3(21)*getJetFF4(-i);
    }
    res=(dsgq2*M2gq2gq(_shat,_that,_uhat)+dsgq2_*M2gq2gq(_shat,_uhat,_that))*getXsection();
    // std::cout << "Gluon: res7=" << res << std::endl;
  }
  return res*_unitGeV2pb;  //return with unit pb
}

//gg->gg
double DiJetLO::SigmaLOgg2gg(){
  double res=0.0;
  if(KinematicsQcut()){
    double dsgggg=0.0;
    //gg->gg
    dsgggg+=getPDF1(21)*getPDF2(21)*getJetFF3(21)*getJetFF4(21);
    
    res=dsgggg*M2gg2gg(_shat,_that,_uhat)*getXsection();
    // std::cout << "res8=" << res << std::endl;
  }
  return res*_unitGeV2pb;  //return with unit pb
}


//####################################################################################################
//Leading order kinematic parameters
void DiJetLO::setLOParameters(const double y3c, const double y4d){
  setKinematics(y3c,y4d);
  setPDF1(); setPDF2();
  setJetFF3(); setJetFF4();
  setXsection();
}

//Leading order cross section
double DiJetLO::SigmaLO(){
  double res=0.0;
  if(KinematicsQcut()){
    res=SigmaLOqqp2qqp()+SigmaLOqq2qq()+SigmaLOqqb2qpqpb()+SigmaLOqqb2qqb()+SigmaLOqqb2gg()+SigmaLOgq2gq()+SigmaLOgg2qqb()+SigmaLOgg2gg();
  }
  return res;
}

double DiJetLO::SigmaLOQuark(){
  double res=0.0;
  if(KinematicsQcut()){
    res=SigmaLOqqp2qqp()+SigmaLOqq2qq()+SigmaLOqqb2qpqpb()+SigmaLOqqb2qqb()+SigmaLOgq2gq_Q()+SigmaLOgg2qqb();
  }
  return res;
}

double DiJetLO::SigmaLOGluon(){
  double res=0.0;
  if(KinematicsQcut()){
    res=SigmaLOqqb2gg()+SigmaLOgq2gq_G()+SigmaLOgg2gg();
  }
  return res;
}


//####################################################################################################
//MC integration
void DiJetLO::setCalls(const size_t calls){
  _calls=calls;
}

void DiJetLO::setupMC(double (*func)(double *, size_t, void *), size_t dim, double *rangemin, double *rangemax, void *para){
  gsl_rng_env_setup ();
  _Type = gsl_rng_default;
  _rng  = gsl_rng_alloc (_Type);

  _dim=dim;
  _rangemin=rangemin; _rangemax=rangemax;
  _mcFun.f=func; _mcFun.dim=dim; _mcFun.params=para;

#ifdef MC_MISER
  _mcStatus = gsl_monte_miser_alloc (dim);
  std::cout << "miser is used..." << std::endl;
#else
  _mcStatus = gsl_monte_vegas_alloc (dim);
  std::cout << "vegas is used..." << std::endl;
#endif
}

void DiJetLO::calculateMC(double &res, double &err){  //calculate the total cross section in pb unit
#ifdef MC_MISER
  gsl_monte_miser_integrate (&_mcFun, _rangemin, _rangemax, _dim, _calls, _rng, _mcStatus, &res, &err);
#else
  gsl_monte_vegas_integrate (&_mcFun, _rangemin, _rangemax, _dim, _calls, _rng, _mcStatus, &res, &err);
#endif
}

void DiJetLO::cleanupMC(){
  gsl_rng_free (_rng);
#ifdef MC_MISER
  gsl_monte_miser_free (_mcStatus);
#else
  gsl_monte_vegas_free (_mcStatus);
#endif
}

