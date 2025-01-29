#include <iostream>
#include <fstream>
#include "Pythia8/Pythia.h"

using namespace Pythia8;


int main(){

    std::ofstream OutputFile;
    // //7 TeV {3.5,9.5,22,35}
    // OutputFile.open("../Output/ATLAS_DiJet_Hadron_7000GeV_35GeV_Tune21.dat");
    // OutputFile.open("../Output/ATLAS_DiJet_Parton_7000GeV_35GeV_Tune21.dat");
    // //13 TeV {90.,290.,590.,980.,1580.}
    OutputFile.open("../Output/ATLAS_DiJet_Hadron_13000GeV_1580GeV_Tune21.dat");
    // OutputFile.open("../Output/ATLAS_DiJet_Parton_13000GeV_980GeV_Tune21.dat");


    //------------------------------------------------------------
    //Pythia initialization
    //Generator
    Pythia pythia;
    // //Shorthand for the event record in pythia
    // Event &event = pythia.event;

    //Tuning parameters
    // pythia.readString("Tune:pp = 5" );  //Tune 4C
    // pythia.readString("Tune:pp = 14");  //Monash 2013 tune to both e^+e^- and pp/pbarp data
    // pythia.readString("Tune:pp = 18");  //CMS Tune MonashStar
    pythia.readString("Tune:pp = 21");     //ATLAS A14 central tune with NNPDF2.3LO

    // //PDF selection
    // pythia.readString("PDF:pSet = 13");

    //Initial condition
    // pythia.readString("Beams:eCM = 7000.");   //GeV, collision CM energy
    // pythia.readString("Beams:eCM = 8000.");   //GeV, collision CM energy
    pythia.readString("Beams:eCM = 13000.");  //GeV, collision CM energy
    pythia.readString("Beams:idA = 2212");    //2212 is proton
    pythia.readString("Beams:idB = 2212");

    //Inclusive jet
    pythia.readString("HardQCD:all = off");
    pythia.readString("PhaseSpace:pTHatMin = 1580.");  //minimal pT scale
    // pythia.readString("SoftQCD:all = on");  //switch for the group of all soft QCD processes
    //Dijet (Light quarks and gluons)
    pythia.readString("HardQCD:gg2gg = on");
    pythia.readString("HardQCD:gg2qqbar = on");
    pythia.readString("HardQCD:qg2qg = on");
    pythia.readString("HardQCD:qq2qq = on");
    pythia.readString("HardQCD:qqbar2gg = on");
    pythia.readString("HardQCD:qqbar2qqbarNew = on");
    pythia.readString("HardQCD:nQuarkNew = 3");
    // //Heavy-Flavor
    // // pythia.readString("HardQCD:gg2ccbar = on");
    // // pythia.readString("HardQCD:qqbar2ccbar = on");
    // pythia.readString("HardQCD:hardccbar = on");
    // // pythia.readString("HardQCD:gg2bbbar = on");
    // // pythia.readString("HardQCD:qqbar2bbbar = on");
    // pythia.readString("HardQCD:hardbbbar = on");
    // // Gamma-Jet
    // pythia.readString("PromptPhoton:qg2qgamma = on");
    // pythia.readString("PromptPhoton:qqbar2ggamma = on");
    // pythia.readString("PromptPhoton:gg2ggamma = on");

    // //Process level
    // pythia.readString("ProcessLevel:all = off");  //the trick

    //Parton level
    pythia.readString("PartonLevel:all = on");
    // pythia.readString("PartonLevel:MPI = off");
    // pythia.readString("PartonLevel:ISR = on");
    // pythia.readString("PartonLevel:FSR = on");
    // pythia.readString("PartonLevel:FSRinProcess = on");
    // pythia.readString("PartonLevel:FSRinResonances = on");
    // pythia.readString("PartonLevel:earlyResDec = off");
    // pythia.readString("TimeShower:pTmin = 2.0");  //parton shower cut-off pT for QCD emissions

    //Hadron level
    pythia.readString("HadronLevel:all = on");  //off to only get partons from hard process
    // pythia.readString("HadronLevel:Hadronize = on");
    // pythia.readString("HadronLevel:Decay = off");

    //Prevent unstable particles from decaying
    pythia.readString("ParticleDecays:limitTau0 = on");
    pythia.readString("ParticleDecays:tau0Max = 10");
    //Random number
    pythia.readString("Random:setSeed = on");
    pythia.readString("Random:seed = 0");  //a value 0 gives a random seed based on the time
    //Initialization
    pythia.init();


    //------------------------------------------------------------
    //Event loop
    int nEvents=1000000;
    for(int iEvent=0; iEvent<nEvents; iEvent++){
        //Generate next event and skip it if error
        if(!pythia.next()) continue;

        //Output the number of event particles
        int num=0;
        for(int i=0; i<pythia.event.size(); i++){
            if(pythia.event[i].isFinal()) num+=1;
        }
        OutputFile << "#    " << iEvent+1 << "   " << num << "   " << pythia.info.sigmaGen(0) << "   " << pythia.info.sigmaErr(0) << std::endl;  //cross section, in units of mb
        //Hard process 3+4->5+6
        OutputFile << "##   " << pythia.event[3].status() << "   " << pythia.event[3].id() << "   " << pythia.event[3].px() << "   " << pythia.event[3].py() << "   " << pythia.event[3].pz() << "   " << pythia.event[3].e() << "   " << pythia.event[3].m() << "   " << pythia.event[3].eta() << "   " << pythia.event[3].phi() << std::endl;
        OutputFile << "##   " << pythia.event[4].status() << "   " << pythia.event[4].id() << "   " << pythia.event[4].px() << "   " << pythia.event[4].py() << "   " << pythia.event[4].pz() << "   " << pythia.event[4].e() << "   " << pythia.event[4].m() << "   " << pythia.event[4].eta() << "   " << pythia.event[4].phi() << std::endl;
        OutputFile << "##   " << pythia.event[5].status() << "   " << pythia.event[5].id() << "   " << pythia.event[5].px() << "   " << pythia.event[5].py() << "   " << pythia.event[5].pz() << "   " << pythia.event[5].e() << "   " << pythia.event[5].m() << "   " << pythia.event[5].eta() << "   " << pythia.event[5].phi() << std::endl;
        OutputFile << "##   " << pythia.event[6].status() << "   " << pythia.event[6].id() << "   " << pythia.event[6].px() << "   " << pythia.event[6].py() << "   " << pythia.event[6].pz() << "   " << pythia.event[6].e() << "   " << pythia.event[6].m() << "   " << pythia.event[6].eta() << "   " << pythia.event[6].phi() << std::endl;

        //Particle loop
        int    index=0,id;
        double px,py,pz,e,m;
        for(int iParti=0; iParti<pythia.event.size(); iParti++){
            if(pythia.event[iParti].isFinal()){
                index+=1;
                id=pythia.event[iParti].id();
                px=pythia.event[iParti].px();
                py=pythia.event[iParti].py();
                pz=pythia.event[iParti].pz();
                e =pythia.event[iParti].e();
                m =pythia.event[iParti].m();
                OutputFile << index << "   " << id << "   " << px << "   " << py << "   " << pz << "   " << e << "   " << m << std::endl;
            }
        }
    }
    OutputFile.close();

    //Statistics on event generation
    pythia.stat();
    std::cout << "#    " << nEvents << "   " << pythia.info.sigmaGen() << "   " << pythia.info.sigmaErr() << std::endl;  //cross section, in units of mb

    return 0;
}
