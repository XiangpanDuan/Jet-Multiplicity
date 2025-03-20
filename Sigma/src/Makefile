CXX = g++
CXXFLAGS = -std=c++11 -Wall -O3

IPDF = ${shell lhapdf-config --cflags}
LPDF = ${shell lhapdf-config --libs}
IGSL = ${shell gsl-config --cflags}
LGSL = ${shell gsl-config --libs}
# IPDF = -I/opt/homebrew/Cellar/lhapdf/6.5.3/include
# LPDF = -L/opt/homebrew/Cellar/lhapdf/6.5.3/lib -lLHAPDF
# IGSL = -I/Users/duanxiangpan/Desktop/Work/Model/GSL/gsl-install/include
# LGSL = -L/Users/duanxiangpan/Desktop/Work/Model/GSL/gsl-install/lib -lgsl -lgslcblas -lm


TARGET = pp.exe
all: $(TARGET)


Particle.o: Particle.cpp Particle.h
	$(CXX) $(CXXFLAGS) -c Particle.cpp

QCD.o: QCD.cpp QCD.h
	$(CXX) $(CXXFLAGS) $(IPDF) $(IGSL) -c QCD.cpp

InclusiveJetLO.o: InclusiveJetLO.cpp InclusiveJetLO.h
	$(CXX) $(CXXFLAGS) $(IPDF) $(IGSL) -c InclusiveJetLO.cpp

#Calculate Leading-Order cross section
CalculateLO.o: CalculateLO.cpp
	$(CXX) $(CXXFLAGS) $(IPDF) $(IGSL) -c CalculateLO.cpp

pp.exe: CalculateLO.o InclusiveJetLO.o QCD.o Particle.o
	$(CXX) $(CXXFLAGS) $^ $(LPDF) $(LGSL) -o $@



clean:
	rm -f *.o *.*~ *.exe
