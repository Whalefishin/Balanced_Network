# Configure compiler settings
CXX=clang++
CXXFLAGS=-g -std=c++11 -Werror
# The object files for the program.
OFILES = \
	goal.o \
	gui/cpipe.o \
	gui/subprocess.o \
	gui/dataReader.o \
	gui/dataWriter.o \
	gui/railwayDataReader.o \
	gui/railwayDataWriter.o \
	railwayGUI.o \
	ioUtils.o \
	railwayGame.o
# The header files for the program
HFILES = $(shell ls -1 *.h gui/*.h)

# UnitTest++ keeps its object files in this library.
UNITTEST_LIB = -lUnitTest++

all: railway tests manualTests 

# This target builds your main program.
railway: $(HFILES) $(OFILES) main.o
	$(CXX) $(CXXFLAGS) -o $@ main.o $(OFILES)


# This target builds your tests.
tests: tests.o $(HFILES) $(OFILES)
	$(CXX) $(CXXFLAGS) -o $@ tests.o $(OFILES)  $(UNITTEST_LIB)

# This target builds your manual testing program.
manualTests: manualTests.o $(HFILES) $(OFILES)
	$(CXX) $(CXXFLAGS) -o $@ manualTests.o $(OFILES)

# This target describes how to compile a .o file from a .cpp file.
%.o: %.cpp $(HFILES)
	$(CXX) $(CXXFLAGS) -c -o $@ $<

# This target deletes the temporary files we have built.
.PHONY: clean
clean:
	rm -f *.o
	rm -f gui/*.o
	rm -f railway
	rm -f tests
	rm -f manualTests
