CXX=g++
CXXFLAGS=-O2 -std=c++17
LIBS=-lm

all : makerobin

%: %.cpp
	$(CXX) $(CXXFLAGS) -o $@ $< $(LIBS)

clean: 
	rm -f *.o makerobin

fileclean:
	rm -f *.obj
