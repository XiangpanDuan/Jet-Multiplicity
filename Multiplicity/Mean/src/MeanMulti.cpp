#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_bessel.h>


//Mean multiplicity
//Formula 5.51(P124) for DLA and 7.34(P174) for MDLA in [Basics of Perturbative QCD]
int main()
{
  int    mode=0;  //0:DLA; 1:MDLA
  double Q0=0.5;  //infrared cutoff
  double Rsize=0.4;

  std::stringstream ss;
  if(mode==0) ss << "../Output/DLA_Multi_Q"  << Q0 << "_R" << Rsize << ".dat";
  if(mode==1) ss << "../Output/MDLA_Multi_Q" << Q0 << "_R" << Rsize << ".dat";
  std::string OutputString=ss.str();
  std::ofstream OutputFile;
  OutputFile.open(OutputString);
  OutputFile << "# pT <multiQ> <multiG>" << std::endl;

  double par=1.0;
  double LambdaQCD=0.2457484;
  double CF=4./3.;
  // double CA=3.0;
  double Nc=3.0;
  int    nf=3;
  double a=11./3.*Nc+2.*nf/(3.*Nc*Nc);
  double b=11./3.*Nc-2./3.*nf;
  double A=std::sqrt(16.*Nc/b);
  double B=a/b;
  double ypT,lambdapT,YpT;
  double x1,x2;
  double multiQ,multiG;
  double pT,Q;
  const int    pTnum=25;
  const double pTmin=100.;
  const double pTmax=2600.;
  const double pTbin=(pTmax-pTmin)/pTnum;
  for(int i=0; i<pTnum; i++){
    pT=pTmin+i*pTbin;
    Q=pT*Rsize;  //R dependence
    ypT=std::log(Q/Q0);
    lambdapT=std::log(Q0/LambdaQCD);
    YpT=ypT+lambdapT;
    x1=A*std::sqrt(YpT);
    x2=A*std::sqrt(lambdapT);
    if(mode==0){
      //DLA
      multiG=par*x1*(gsl_sf_bessel_I1(x1)*gsl_sf_bessel_K0(x2)+gsl_sf_bessel_K1(x1)*gsl_sf_bessel_I0(x2));
      multiQ=(CF/Nc)*(multiG-1.)+1.;
    }
    if(mode==1){
      //MDLA
      multiG=par*x1*std::pow(x2/x1,B)*(gsl_sf_bessel_Inu(B+1,x1)*gsl_sf_bessel_Knu(B,x2)+gsl_sf_bessel_Knu(B+1,x1)*gsl_sf_bessel_Inu(B,x2));
      multiQ=(CF/Nc)*(multiG-1.)+1.;
    }

    OutputFile << pT << " " << multiQ << " " << multiG << std::endl;
  }

  OutputFile.close();


  return 0;
}