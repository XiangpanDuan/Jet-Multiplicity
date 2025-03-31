CXX = g++
CXXFLAGS = -std=c++11 -Wall -O3

IGSL = ${shell gsl-config --cflags}
LGSL = ${shell gsl-config --libs}
IPDF = ${shell lhapdf-config --cflags}
LPDF = ${shell lhapdf-config --libs}
# IGSL = -I/Users/duanxiangpan/Desktop/Work/Model/GSL/gsl-install/include
# LGSL = -L/Users/duanxiangpan/Desktop/Work/Model/GSL/gsl-install/lib -lgsl -lgslcblas -lm
# IPDF = -I/opt/homebrew/Cellar/lhapdf/6.5.3/include
# LPDF = -L/opt/homebrew/Cellar/lhapdf/6.5.3/lib -lLHAPDF


TARGET = Probability.exe
all: $(TARGET)

Probability.o: Probability.cpp
	$(CXX) $(CXXFLAGS) $(IGSL) $(IPDF) -c Probability.cpp

$(TARGET): Probability.o
	$(CXX) $(CXXFLAGS) $^ $(LGSL) $(LPDF) -o $@


clean:
	rm -f *.o *.*~ *.exe
