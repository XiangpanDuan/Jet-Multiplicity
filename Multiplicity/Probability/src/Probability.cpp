#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <cmath>
#include <gsl/gsl_integration.h>
#include "LHAPDF/LHAPDF.h"


double pT,Q;
const int    pTnum=25;
const double pTmin=100.;
const double pTmax=2600.;
const double pTbin=(pTmax-pTmin)/pTnum;
const int    Pnum=1000;
const int    num=5000;
std::vector<std::vector<double>> Pq(pTnum, std::vector<double>(Pnum, 0.0));
std::vector<std::vector<double>> Pg(pTnum, std::vector<double>(Pnum, 0.0));
std::vector<std::vector<double>> S (pTnum, std::vector<double>(Pnum, 0.0));
std::vector<std::vector<std::vector<double>>> PTable(pTnum, std::vector<std::vector<double>>(Pnum, std::vector<double>(num+1, 0.0)));
std::vector<std::vector<std::vector<double>>> STable(pTnum, std::vector<std::vector<double>>(Pnum, std::vector<double>(num+1, 0.0)));
double Rsize=0.4;  //jet cone size
double Q0=1.0;     //infrared cutoff
// double alphas=0.1;
double Nc=3.0;
double CA=3.0;
double CF=4./3.;
double ymin,ymax,ybin;
double Pqsum,Pgsum;
double nMulq,nMulg;
int    ival,jval;


//------------------------------------------------------------
//Alphas(kT^2)
LHAPDF::PDF *pdf = LHAPDF::mkPDF("CT18NLO",0);
double AlphaS(const double yp){
  double kT=Q0*std::exp(yp);
  return pdf->alphasQ(kT);
}

//------------------------------------------------------------
//Linear interpolation
double LinearInterpolation(const double yp){
  int ylow=(int)(std::floor(yp/ybin));
  int yhigh=ylow+1;
  if(ylow<0 || yhigh>num) return 0.0;
  return PTable[ival][jval][ylow]+(PTable[ival][jval][yhigh]-PTable[ival][jval][ylow])*(yp-ylow*ybin)/ybin;
}

//------------------------------------------------------------
//Integral table (0,y')
double Function(const double yp, void *params){
  (void) (params);
  double gamma0square=2.*Nc/M_PI*AlphaS(yp);
  // double gamma0square=2.*Nc/M_PI*alphas;
  double Pval=LinearInterpolation(yp);
  return (ymax-yp)*gamma0square*Pval;
}

//------------------------------------------------------------
//Integral (0,y)
double FunctionBase(const double yp, void *params){
  (void) (params);
  double gamma0square=2.*Nc/M_PI*AlphaS(yp);
  // double gamma0square=2.*Nc/M_PI*alphas;
  return (ymax-yp)*gamma0square;
}


//------------------------------------------------------------
//Calculate probability
void CalProbability(){

  //------------------------------
  //Output files
  std::stringstream ss[3];
  ss[0] << "../Output/DLA_Probability_Q" << Q0 << "_R" << Rsize << "_AlphaS" << "_Multi.dat";
  ss[1] << "../Output/DLA_Probability_Q" << Q0 << "_R" << Rsize << "_AlphaS" << "_Quark.dat";
  ss[2] << "../Output/DLA_Probability_Q" << Q0 << "_R" << Rsize << "_AlphaS" << "_Gluon.dat";
  std::string OutputString[3];
  std::ofstream OutputFile[3];
  for(int i=0; i<3; i++){
    OutputString[i]=ss[i].str();
    OutputFile[i].open(OutputString[i]);
    OutputFile[i] << "# pT";
    if(i==0){OutputFile[i] << " <MultiQ>" << " <MultiG>" << std::endl;}
    if(i!=0){
      for(int j=0; j<Pnum; j++){
        OutputFile[i] << " P" << j+1;
      }
      OutputFile[i] << " <Multi>" << std::endl;
    }
  }

  //------------------------------
  //pT loop
  for(int i=0; i<pTnum; i++){
    ival=i;
    pT=pTmin+pTbin*i;
    Q=pT*Rsize;  //jet cone size dependence
    ybin=std::log(Q/Q0)/num;

    OutputFile[0] << pT << " ";
    OutputFile[1] << pT << " ";
    OutputFile[2] << pT << " ";
    Pqsum=0.0;
    Pgsum=0.0;
    nMulq=0.0;
    nMulg=0.0;

    double res,err;

    //------------------------------
    //Probability loop
    for(int j=0; j<Pnum; j++){
      jval=j;

      //------------------------------
      //Probability integral table
      if(j==0) PTable[i][j][0]=1.0;  //the minimum interpolation for P1(0) in the table
      if(j!=0) PTable[i][j][0]=0.0;  //the minimum interpolation for Pn(0) in the table
      for(int k=1; k<(num+1); k++){
        ymin=0.0;
        ymax=ybin*k;  //update current y value
        if(j==0){
          gsl_function FunBase;
          FunBase.function = &FunctionBase;
          FunBase.params = nullptr;
          gsl_integration_workspace *workspace = gsl_integration_workspace_alloc(10000);
          res=0.0; err=0.0;
          gsl_integration_qag(&FunBase, ymin, ymax, 1e-6, 1e-6, 10000, 6, workspace, &res, &err);
          // gsl_integration_qags(&FunBase, ymin, ymax, 1e-6, 1e-6, 10000, workspace, &res, &err);
          gsl_integration_workspace_free(workspace);
          PTable[i][j][k]=std::exp(-(CA/Nc)*res);
        }
        if(j>0){
          for(int jj=0; jj<j; jj++){
            PTable[i][j][k]+=(double)(jj+1)/j*PTable[i][j-1-jj][k]*STable[i][jj][k]*(CA/Nc);
          }
        }
        // if(k==1 && (j==0 || j==1)) std::cout << pT << " " << PTable[i][j][k] << std::endl;
        gsl_function Fun;
        Fun.function = &Function;
        Fun.params = nullptr;
        gsl_integration_workspace *wspace = gsl_integration_workspace_alloc(10000);
        res=0.0; err=0.0;
        gsl_integration_qag(&Fun, ymin, ymax, 1e-6, 1e-6, 10000, 6, wspace, &res, &err);
        // gsl_integration_qags(&Fun, ymin, ymax, 1e-6, 1e-6, 10000, wspace, &res, &err);
        gsl_integration_workspace_free(wspace);
        STable[i][j][k]=res;  //not consider k=0

      }//End of probability integral table


      //------------------------------
      //Probability value
      ymin=0.0;
      ymax=std::log(Q/Q0);
      if(j==0){
        gsl_function FBase;
        FBase.function = &FunctionBase;
        FBase.params = nullptr;
        gsl_integration_workspace *wsp = gsl_integration_workspace_alloc(10000);
        res=0.0; err=0.0;
        gsl_integration_qag(&FBase, ymin, ymax, 1e-6, 1e-6, 10000, 6, wsp, &res, &err);
        // gsl_integration_qags(&FBase, ymin, ymax, 1e-6, 1e-6, 10000, wsp, &res, &err);
        gsl_integration_workspace_free(wsp);
        Pq[i][j]=std::exp(-(CF/Nc)*res);
        Pg[i][j]=std::exp(-(CA/Nc)*res);
      }
      if(j>0){
        for(int jj=0; jj<j; jj++){
          Pq[i][j]+=(double)(jj+1)/j*Pq[i][j-1-jj]*S[i][jj]*(CF/Nc);
          Pg[i][j]+=(double)(jj+1)/j*Pg[i][j-1-jj]*S[i][jj]*(CA/Nc);
        }
      }
      gsl_function F;
      F.function = &Function;
      F.params = nullptr;
      gsl_integration_workspace *ws = gsl_integration_workspace_alloc(10000);
      res=0.0; err=0.0;
      gsl_integration_qag(&F, ymin, ymax, 1e-6, 1e-6, 10000, 6, ws, &res, &err);
      // gsl_integration_qags(&F, ymin, ymax, 1e-6, 1e-6, 10000, ws, &res, &err);
      gsl_integration_workspace_free(ws);
      S[i][j]=res;


      //------------------------------
      //Output files
      nMulq+=(j+1)*Pq[i][j];
      nMulg+=(j+1)*Pg[i][j];
      OutputFile[1] << Pq[i][j] << " ";
      OutputFile[2] << Pg[i][j] << " ";
      Pqsum+=Pq[i][j];
      Pgsum+=Pg[i][j];

    }//End of probability loop

    OutputFile[0] << nMulq << " " << nMulg << std::endl;
    OutputFile[1] << nMulq << std::endl;
    OutputFile[2] << nMulg << std::endl;

    //Test P1+P2+P3+...+Pn=1
    std::cout << "pT=" << pT << ",  Pqsum=" << Pqsum << ",  Pgsum=" << Pgsum << std::endl;

  }//End of pT loop

  for(int i=0; i<3; i++){
    OutputFile[i].close();
  }

}



//####################################################################################################
//Main function
int main(){

  std::cout << "################################################################################" << std::endl;
  std::cout << "# Jet Multiplicity Distribution P(n) in Double Logarithmic Approximation (DLA) #" << std::endl;
  std::cout << "################################################################################" << std::endl;


  CalProbability();


  return 0;
}

