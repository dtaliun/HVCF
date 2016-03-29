GTESTLIB = /Users/dtaliun/Documents/gtest/googletest/cmake/test/
GTESTINCS = /Users/dtaliun/Documents/gtest/googletest/include/
HDF5LIB=/Users/dtaliun/Documents/hdf5-1.10.0-pre2/hdf5/lib
HDF5INCS=/Users/dtaliun/Documents/hdf5-1.10.0-pre2/hdf5/include

AUXLIBS = ../../auxc/FileReader/src/*.o \
		../../auxc/MiniVCF/src/*.o
AUXDIRS = ../../auxc/FileReader/src \
		../../auxc/MiniVCF/src

APPLIBS = ../src/*.o
APPDIRS = ../src

CXX = g++
CXXFLAGS = -std=c++11 -Wall -L/opt/local/lib/ -lz -L$(GTESTLIB) -L$(HDF5LIB) -lhdf5 -lgtest -larmadillo
INCS = -I$(GTESTINCS) -I$(HDF5INCS)

OBJECTS = HVCFTest.o Main_TestAll.o

.PHONY: all auxlibs applibs

all: auxlibs applibs testAll

auxlibs:
	@for auxdir in $(AUXDIRS); do \
		(cd $${auxdir} && make -j 4) || exit 1; \
	done
	
applibs:
	@for appdir in $(APPDIRS); do \
		(cd $${appdir} && make -j 4) || exit 1; \
	done

testAll: $(OBJECTS)
	$(CXX) $(CXXFLAGS) $(INCS) -o $@ $(OBJECTS) $(AUXLIBS) $(APPLIBS)
	
.cpp.o:
	$(CXX) $(CXXFLAGS) $(INCS) -c -o $@ $<
	
clean:
	rm -f testAll *.o  ../src/*.o