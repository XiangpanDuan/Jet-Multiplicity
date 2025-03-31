CXX = g++
CXXFLAGS = -std=c++11 -Wall -O3

IGSL = ${shell gsl-config --cflags}
LGSL = ${shell gsl-config --libs}
# IGSL = -I/Users/duanxiangpan/Desktop/Work/Model/GSL/gsl-install/include
# LGSL = -L/Users/duanxiangpan/Desktop/Work/Model/GSL/gsl-install/lib -lgsl -lgslcblas -lm


TARGET = MeanMulti_nf.exe
all: $(TARGET)

MeanMulti_nf.o: MeanMulti_nf.cpp
	$(CXX) $(CXXFLAGS) $(IGSL) -c MeanMulti_nf.cpp

$(TARGET): MeanMulti_nf.o
	$(CXX) $(CXXFLAGS) MeanMulti_nf.o $(LGSL) -o $@
	

clean:
	rm -f *.o *.*~ *.exe
