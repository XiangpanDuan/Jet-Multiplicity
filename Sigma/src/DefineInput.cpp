#include <string>

#define  MODE   0                    //0:pp collisions


namespace Input{

//Initial conditions
const double Ecm=13000.;             //GeV, collisional energy in CoM frame
const std::string name="parton";     //particle name in hard scattering: parton,d,u,s,c,b,t,g,photon,Z0,Wp,Higgs
const int scale=0;                   //0,+1,-1: pT scale to control the pT errer bar in pdf and alphas from LHAPDF
const unsigned int nf=3;             //quark flavors
const unsigned int nloop=1;          //Λ_QCD calculation (1:LO; 2:NLO)

//Monte Carlo calls
const size_t calls=50000;

//Range of integration
const int    pTnum=250;              //number of transverse momentum of final observed jet
const double pTmin=95.;              //GeV, minimum transverse momentum of final observed jet
const double pTmax=2595.;            //GeV, maximal transverse momentum of final observed jet
const double rap3 =2.1;              //rapidity of final observed jet
const double rap4 =2.1;              //rapidity of final another  jet

}


////////////////////////////////////////////////////////////
//Instruction
//MODE=0     "make" or "make pp.exe"