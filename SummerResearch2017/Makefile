# Configure compiler settings
CXX=clang++
CXXFLAGS= -g -std=c++11 -Werror -O3
# The object files for the program.
OFILES = Balanced_Network.o Neuron.o
# The header files for the program
HFILES = Balanced_Network.h Neuron.h Statistics.h
# UnitTest++ keeps its object files in this directory.
# UNITTEST_LIB = -lUnitTest++

all: main

# This target builds your main program.
main: $(HFILES) $(OFILES) main.o
	$(CXX) $(CXXFLAGS) -o $@ main.o $(OFILES)

# This target describes how to compile a .o file from a .cpp file.
%.o: %.cpp $(HFILES)
	$(CXX) $(CXXFLAGS) -c -o $@ $<

# This target deletes the temporary files we have built.
.PHONY: clean all
clean:
	rm -f *.o
	rm -f main
