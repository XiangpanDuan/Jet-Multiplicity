CXX = g++
CXXFLAGS = -std=c++11 -Wall -O3

IPDF = ${shell lhapdf-config --cflags}
LPDF = ${shell lhapdf-config --libs}
IGSL = ${shell gsl-config --cflags}
LGSL = ${shell gsl-config --libs}
# IPDF = -I/opt/homebrew/Cellar/lhapdf/6.5.3/include
# LPDF = -L/opt/homebrew/Cellar/lhapdf/6.5.3/lib -lLHAPDF
# IGSL = -I/Users/duanpanpan/Desktop/Work/Model/GSL/gsl-install/include
# LGSL = -L/Users/duanpanpan/Desktop/Work/Model/GSL/gsl-install/lib -lgsl -lgslcblas -lm


TARGET = pp.exe
all: $(TARGET)


Particle.o: Particle.cpp Particle.h
	$(CXX) $(CXXFLAGS) -c Particle.cpp

QCD.o: QCD.cpp QCD.h
	$(CXX) $(CXXFLAGS) $(IPDF) $(IGSL) -c QCD.cpp

DiJetLO.o: DiJetLO.cpp DiJetLO.h
	$(CXX) $(CXXFLAGS) $(IPDF) $(IGSL) -c DiJetLO.cpp

#Calculate Leading-Order cross section
CalDiJetLO.o: CalDiJetLO.cpp
	$(CXX) $(CXXFLAGS) $(IPDF) $(IGSL) -c CalDiJetLO.cpp

pp.exe: CalDiJetLO.o DiJetLO.o QCD.o Particle.o
	$(CXX) $(CXXFLAGS) $^ $(LPDF) $(LGSL) -o $@



clean:
	rm -f *.o *.*~ *.exe
