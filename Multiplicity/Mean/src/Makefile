CXX = g++
CXXFLAGS = -std=c++11 -Wall -O3

IGSL = ${shell gsl-config --cflags}
LGSL = ${shell gsl-config --libs}
# IGSL = -I/Users/duanxiangpan/Desktop/Work/Model/GSL/gsl-install/include
# LGSL = -L/Users/duanxiangpan/Desktop/Work/Model/GSL/gsl-install/lib -lgsl -lgslcblas -lm


TARGET = MeanMulti.exe
all: $(TARGET)

MeanMulti.o: MeanMulti.cpp
	$(CXX) $(CXXFLAGS) $(IGSL) -c MeanMulti.cpp

$(TARGET): MeanMulti.o
	$(CXX) $(CXXFLAGS) MeanMulti.o $(LGSL) -o $@
	

clean:
	rm -f *.o *.*~ *.exe
