CXX = g++
CXXFLAGS = -std=c++11 -Wall -O3

IGSL = ${shell gsl-config --cflags}
LGSL = ${shell gsl-config --libs}
# IGSL = -I/Users/duanxiangpan/Desktop/Work/Model/GSL/gsl-install/include
# LGSL = -L/Users/duanxiangpan/Desktop/Work/Model/GSL/gsl-install/lib -lgsl -lgslcblas -lm


TARGET = Probability.exe
all: $(TARGET)

Probability.o: Probability.cpp
	$(CXX) $(CXXFLAGS) $(IGSL) -c Probability.cpp

$(TARGET): Probability.o
	$(CXX) $(CXXFLAGS) $^ $(LGSL) -o $@


clean:
	rm -f *.o *.*~ *.exe
