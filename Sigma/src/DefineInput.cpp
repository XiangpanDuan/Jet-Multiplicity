#include <string>

#define  MODE   0                    //0:pp collisions


namespace Input{

//Initial condition
const double Ecm=13000.;             //GeV, collisional energy in CoM frame
const std::string name="parton";     //particle name in hard scattering: parton,d,u,s,c,b,t,g,photon,Z0,Wp,Higgs
const int nf=3;                      //quark flavors
const int nloop=1;                   //Λ_QCD calculation with n-loop β-function coefficient
const int scale=0;                   //0,+1,-1: pT scale (2^scale*pT) to control the pT error range in PDFs and AlphaS from LHAPDF

//Monte Carlo calls
const size_t calls=50000;

//Range of integration
const int    pTnum=25;               //number of transverse momentum of final observed jet
const double pTmin=100.;             //GeV, minimum transverse momentum of final observed jet
const double pTmax=2600.;            //GeV, maximal transverse momentum of final observed jet
const double rap3 =2.1;              //rapidity of final observed jet
const double rap4 =2.1;              //rapidity of final another  jet, rap4=∞(infinite) means to observe the inclusive jet

}


////////////////////////////////////////////////////////////
//Instruction
//MODE=0     "make" or "make pp.exe"