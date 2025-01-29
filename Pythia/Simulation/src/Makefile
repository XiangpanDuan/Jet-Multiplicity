CXX = g++
CXXFLAGS = -std=c++11 -Wall -O3

IPYTHIA = ${shell pythia8-config --cxxflags}
LPYTHIA = ${shell pythia8-config --libs}
# IPYTHIA = -I/Users/duanxiangpan/Desktop/Work/Model/PYTHIA8/pythia8310/include
# LPYTHIA = -L/Users/duanxiangpan/Desktop/Work/Model/PYTHIA8/pythia8310/lib -lpythia8 -ldl


OBJ = Pythia_Multiplicity_ATLAS
# OBJ = Pythia_Sigma_ATLAS

TARGET =  Calculate.exe
all: $(TARGET)

${OBJ}.o: ${OBJ}.cpp
	$(CXX) $(CXXFLAGS) $(IPYTHIA) -c ${OBJ}.cpp

$(TARGET): ${OBJ}.o
	$(CXX) $(CXXFLAGS) $^ $(LPYTHIA) -o $@


clean:
	rm -f *.o *.*~ *.exe
