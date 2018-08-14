
CXXFLAGS = -std=c++14

all: src/xinvert.cpp src/main.cpp
	cd bin; $(CXX) $(CXXFLAGS) -c ../src/xinvert.cpp
	cd bin; $(CXX) -$(CXXFLAGS) -c ../src/main.cpp
	$(CXX) $(CXXFLAGS) -oinvert bin/xinvert.o bin/main.o

clean:
	-rm invert bin/*.o
