#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <cmath>
#include <gsl/gsl_integration.h>


double pT,Q;
double pTmin=100.;
double pTmax=2600.;
const  int pTnum=250;
double pTbin=(pTmax-pTmin)/pTnum;
const  int Pnum=200;
const  int num=1000;
std::vector<std::vector<double>> Pq(pTnum, std::vector<double>(Pnum, 0.0));
std::vector<std::vector<double>> Pg(pTnum, std::vector<double>(Pnum, 0.0));
std::vector<std::vector<double>> S (pTnum, std::vector<double>(Pnum, 0.0));
std::vector<std::vector<std::vector<double>>> PTable(pTnum, std::vector<std::vector<double>>(Pnum, std::vector<double>(num, 0.0)));
std::vector<std::vector<std::vector<double>>> STable(pTnum, std::vector<std::vector<double>>(Pnum, std::vector<double>(num, 0.0)));
double Rsize=0.4;
double Q0=1.43;  //perturbative sacle
double LambdaQCD=0.245748;
double lambda=std::log(Q0/LambdaQCD);
double Nc=3.;
double CA=3.;
double CF=4./3.;
unsigned int nf=3;
double b=11./3.*Nc-2./3.*nf;
double ymin,ymax,ybin;
double Pqsum,Pgsum;
double nMulq,nMulg;
int    ival,jval;

//------------------------------------------------------------
int Factorial(const int n){
  if(n<=1) return 1;        //0!=1, 1!=1
  return n*Factorial(n-1);  //n!
}

//------------------------------------------------------------
double LinearInterpolation(const double yp){
  int ylow=(int)(std::floor(yp/ybin));
  int yhigh=ylow+1;
  if(ylow<0 || yhigh>=num) return 0.0;
  return PTable[ival][jval][ylow]+(PTable[ival][jval][yhigh]-PTable[ival][jval][ylow])*(yp-ylow*ybin)/ybin;
}

//------------------------------------------------------------
double Function(const double yp, void *params){
  (void) (params);
  double gamma0square=4.*Nc/(b*(yp+lambda));
  double Pval=LinearInterpolation(yp);
  return (ymax-yp)*gamma0square*Pval;
}

//------------------------------------------------------------
double FunctionBase(const double yp, void *params){
  (void) (params);
  double gamma0square=4.*Nc/(b*(yp+lambda));
  return (ymax-yp)*gamma0square;
}


//------------------------------------------------------------
void CalProbability(){

  //Output files
  std::stringstream ss[3];
  ss[0] << "../Output/LL_Probability_Qmed" << Q0 << "_R" << Rsize << "_Multi.dat";
  ss[1] << "../Output/LL_Probability_Qmed" << Q0 << "_R" << Rsize << "_Quark.dat";
  ss[2] << "../Output/LL_Probability_Qmed" << Q0 << "_R" << Rsize << "_Gluon.dat";
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

    //Probability loop
    for(int j=0; j<Pnum; j++){
      jval=j;

      double res,err;

      //Probability integral table
      for(int k=0; k<num; k++){
        ymin=0.0;
        ymax=ybin*(k+1);  //update current y value
        if(j==0){
          gsl_function FunBase;
          FunBase.function = &FunctionBase;
          FunBase.params = nullptr;
          gsl_integration_workspace *workspace = gsl_integration_workspace_alloc(10000);
          res=0.0; err=0.0;
          gsl_integration_qag(&FunBase, ymin, ymax, 0, 1e-3, 10000, 6, workspace, &res, &err);
          // gsl_integration_qags(&FunBase, ymin, ymax, 0, 1e-3, 10000, workspace, &res, &err);
          gsl_integration_workspace_free(workspace);
          PTable[i][j][k]=std::exp(-(CA/Nc)*res);
        }
        if(j>0){
          for(int jj=0; jj<j; jj++){
            PTable[i][j][k]+=(double)(jj+1)/j*PTable[i][j-1-jj][k]*STable[i][jj][k]*(CA/Nc);
          }
        }
        // if(j>0){
        //   for(int jj=0; jj<j; jj++){
        //     PTable[i][j][k]+=(double)(jj+1)/j*PTable[i][j-1-jj][k]*STable[i][jj][k]/Factorial(jj+1)*(CA/Nc);
        //   }
        // }
        gsl_function Fun;
        Fun.function = &Function;
        Fun.params = nullptr;
        gsl_integration_workspace *workspace = gsl_integration_workspace_alloc(10000);
        res=0.0; err=0.0;
        gsl_integration_qag(&Fun, ymin, ymax, 0, 1e-3, 10000, 6, workspace, &res, &err);
        // gsl_integration_qags(&Fun, ymin, ymax, 0, 1e-3, 10000, workspace, &res, &err);
        gsl_integration_workspace_free(workspace);
        STable[i][j][k]=res;
        if(res<0) std::cout << "i=" << i << ", j=" << j << ", k=" << k << ", S=" << STable[i][j][k] << std::endl;
      }


      //Probability value
      ymin=0.0;
      ymax=std::log(Q/Q0);
      if(j==0){
        gsl_function FBase;
        FBase.function = &FunctionBase;
        FBase.params = nullptr;
        gsl_integration_workspace *ws = gsl_integration_workspace_alloc(10000);
        res=0.0; err=0.0;
        gsl_integration_qag(&FBase, ymin, ymax, 0, 1e-3, 10000, 6, ws, &res, &err);
        // gsl_integration_qags(&FBase, ymin, ymax, 0, 1e-3, 10000, ws, &res, &err);
        gsl_integration_workspace_free(ws);
        Pq[i][j]=std::exp(-(CF/Nc)*res);
        Pg[i][j]=std::exp(-(CA/Nc)*res);
      }
      if(j>0){
        for(int jj=0; jj<j; jj++){
          Pq[i][j]+=(double)(jj+1)/j*Pq[i][j-1-jj]*S[i][jj]*(CF/Nc);
          Pg[i][j]+=(double)(jj+1)/j*Pg[i][j-1-jj]*S[i][jj]*(CA/Nc);
        }
      }
      // if(j>0){
      //   for(int jj=0; jj<j; jj++){
      //     Pq[i][j]+=(double)(jj+1)/j*Pq[i][j-1-jj]*S[i][jj]/Factorial(jj+1)*(CF/Nc);
      //     Pg[i][j]+=(double)(jj+1)/j*Pg[i][j-1-jj]*S[i][jj]/Factorial(jj+1)*(CA/Nc);
      //   }
      // }
      gsl_function F;
      F.function = &Function;
      F.params = nullptr;
      gsl_integration_workspace *wsp = gsl_integration_workspace_alloc(10000);
      res=0.0; err=0.0;
      gsl_integration_qag(&F, ymin, ymax, 0, 1e-3, 10000, 6, wsp, &res, &err);
      // gsl_integration_qags(&F, ymin, ymax, 0, 1e-3, 10000, wsp, &res, &err);
      gsl_integration_workspace_free(wsp);
      S[i][j]=res;


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
int main(){

  CalProbability();

  return 0;
}
