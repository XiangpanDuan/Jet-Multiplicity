#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_bessel.h>


//Mean multiplicity
//Formula 5.51(P124) for LLA and 7.34(P174) for MLLA in [Basics of Perturbative QCD]
int main()
{
  int    mode=0;    //0:LL; 1:MLL
  double Q0=1.43;   //perturbative sacle
  double Rsize=0.4;

  std::stringstream ss;
  if(mode==0) ss << "../Output/LL_Multi_Qmed"  << Q0 << "_R" << Rsize << ".dat";
  if(mode==1) ss << "../Output/MLL_Multi_Qmed" << Q0 << "_R" << Rsize << ".dat";
  std::string OutputString=ss.str();
  std::ofstream OutputFile;
  OutputFile.open(OutputString);
  OutputFile << "# pT <multiQ> <multiG>" << std::endl;

  double LambdaQCD=0.245748;
  double CF=4./3.;
  // double CA=3.;
  double Nc=3.;
  unsigned int nf=3;
  double par=1.0;
  // double varLL=1.0;   //LL
  // double varMLL=1.0;  //MLL
  double a=11./3.*Nc+2.*nf/(3.*Nc*Nc);
  double b=11./3.*Nc-2./3.*nf;
  double A=std::sqrt(16.*Nc/b);
  double B=a/b;
  double ypT,lambdapT,YpT;
  double x1,x2;
  double multiQ,multiG;
  double pT,Q;
  int    pTnum=250;
  double pTmin=100.;
  double pTbin=10.;
  // int    pTnum=100;
  // double pTmin=4.;
  // double pTbin=1.;
  for(int i=0; i<pTnum; i++){
    pT=pTmin+i*pTbin;
    Q=pT*Rsize;  //R dependence
    ypT=std::log(Q/Q0);
    lambdapT=std::log(Q0/LambdaQCD);
    YpT=ypT+lambdapT;
    x1=A*std::sqrt(YpT);
    x2=A*std::sqrt(lambdapT);
    if(mode==0){
      // multiQ=varQ*x1*(gsl_sf_bessel_I1(x1)*gsl_sf_bessel_K0(x2)+gsl_sf_bessel_K1(x1)*gsl_sf_bessel_I0(x2));  //LL
      multiG=par*x1*(gsl_sf_bessel_I1(x1)*gsl_sf_bessel_K0(x2)+gsl_sf_bessel_K1(x1)*gsl_sf_bessel_I0(x2));  //LL
      multiQ=(CF/Nc)*(multiG-1.)+1.;
    }
    if(mode==1){
      // multiQ=par*x1*std::pow(x2/x1,B)*(gsl_sf_bessel_Inu(B+1,x1)*gsl_sf_bessel_Knu(B,x2)+gsl_sf_bessel_Knu(B+1,x1)*gsl_sf_bessel_Inu(B,x2));  //MLL
      multiG=par*x1*std::pow(x2/x1,B)*(gsl_sf_bessel_Inu(B+1,x1)*gsl_sf_bessel_Knu(B,x2)+gsl_sf_bessel_Knu(B+1,x1)*gsl_sf_bessel_Inu(B,x2));  //MLL
      multiQ=(CF/Nc)*(multiG-1.)+1.;
    }
    OutputFile << pT << " " << multiQ << " " << multiG << std::endl;
  }
  OutputFile.close();


  // //
  // double TR=0.5;
  // double b0=(11.*Nc-4.*nf*TR)/(12.*M_PI);
  // double a0=1./4.+5.*nf/(54.*M_PI*b0);
  // double varLL=0.016965;    //LLA  from pythia8 simulation
  // double varMLL=0.0403164;  //MLLA from pythia8 simulation
  // for(int i=0; i<pTnum; i++){
  //   pT=pTmin+i*pTbin;
  //   Q=pT*Rsize;
  //   double alphaS=1./(b0*std::log(Q*Q/(LambdaQCD*LambdaQCD)));
  //   if(mode==0){
  //     multiQ=varQ*varLL*std::exp(std::pow(2.*Nc/(M_PI*b0)*std::log(Q*Q/(LambdaQCD*LambdaQCD)),0.5));  //LLA
  //     multiG=varG*varLL*std::exp(std::pow(2.*Nc/(M_PI*b0)*std::log(Q*Q/(LambdaQCD*LambdaQCD)),0.5));  //LLA
  //   }
  //   if(mode==1){
  //     multiQ=varQ*varMLL*std::exp(std::pow(2.*Nc/(M_PI*b0)*std::log(Q*Q/(LambdaQCD*LambdaQCD)),0.5)+a0*std::log(alphaS));  //MLLA
  //     multiG=varG*varMLL*std::exp(std::pow(2.*Nc/(M_PI*b0)*std::log(Q*Q/(LambdaQCD*LambdaQCD)),0.5)+a0*std::log(alphaS));  //MLLA
  //   }
  // }


  return 0;
}