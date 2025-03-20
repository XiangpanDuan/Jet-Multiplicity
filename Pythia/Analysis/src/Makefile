CXX = g++
CXXFLAGS = -std=c++11 -Wall -O3

IROOT = ${shell root-config --cflags}
LROOT = ${shell root-config --libs --glibs}
IFASTJET = ${shell fastjet-config --cxxflags --plugins}
LFASTJET = ${shell fastjet-config --libs --plugins}
# IFASTJET = ${shell /Users/duanxiangpan/Desktop/Work/Model/Fastjet/fastjet-install/bin/fastjet-config --cxxflags --plugins}
# LFASTJET = ${shell /Users/duanxiangpan/Desktop/Work/Model/Fastjet/fastjet-install/bin/fastjet-config --libs --plugins}
# FASTJET = `/Users/duanxiangpan/Desktop/Work/Model/Fastjet/fastjet-install/bin/fastjet-config --cxxflags --libs --plugins`
# IFASTJET = -I/Users/duanxiangpan/Desktop/Work/Model/Fastjet/fastjet-install/include
# LFASTJET = -L/Users/duanxiangpan/Desktop/Work/Model/Fastjet/fastjet-install/lib -lfastjettools -lfastjet
IFJCONTRIB = -I/Users/duanpanpan/Desktop/Work/Model/FastJet/fjcontrib-install/include
LFJCONTRIB = -L/Users/duanpanpan/Desktop/Work/Model/FastJet/fjcontrib-install/lib -lEnergyCorrelator

# OBJ = Jet_Multiplicity_Dijet_ATLAS_7000GeV
# OBJ = Jet_Multiplicity_Dijet_ATLAS_13000GeV
OBJ = Jet_ECF_Multiplicity_Dijet_ATLAS_13000GeV
# OBJ = Jet_Sigma_InclusiveJet_ATLAS_5020GeV

TARGET = Analysis.exe
all: $(TARGET)

$(OBJ).o: $(OBJ).cpp
	$(CXX) $(CXXFLAGS) $(IROOT) $(IFASTJET) $(IFJCONTRIB) -c $(OBJ).cpp

$(TARGET): $(OBJ).o
	$(CXX) $(CXXFLAGS) $^ $(LROOT) $(LFASTJET) $(LFJCONTRIB) -o $@


clean:
	rm -f *.o *.*~ *.exe
