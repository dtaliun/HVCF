GTESTLIB=/Users/dtaliun/Documents/gtest/googletest/googletest/build/
GTESTINCS=/Users/dtaliun/Documents/gtest/googletest/googletest/include/
HDF5LIB=/Users/dtaliun/Documents/hdf5-1.10.0/hdf5/lib
HDF5INCS=/Users/dtaliun/Documents/hdf5-1.10.0/hdf5/include
BLOSCLIB=/Users/dtaliun/Documents/c-blosc-1.7.1/blosc
BLOSCINCS=/Users/dtaliun/Documents/c-blosc-1.7.1/blosc

AUXLIBS = ../../auxc/FileReader/src/*.o \
		../../auxc/MiniVCF/src/*.o
AUXDIRS = ../../auxc/FileReader/src \
		../../auxc/MiniVCF/src

APPLIBS = ../src/*.o
APPDIRS = ../src

BLOSCLIBS = ../src/blosc/*.o
BLOSCDIRS = ../src/blosc

CXX = g++
CXXFLAGS = -std=c++11 -Wall -L/opt/local/lib/ -L$(GTESTLIB) -L$(HDF5LIB) -L$(BLOSCLIB)
LIBS = -lz -lhdf5 -lblosc -larmadillo -lgtest
INCS = -I$(GTESTINCS) -I$(HDF5INCS) -I$(BLOSCINCS)

OBJECTS = HVCFTestReadWrite.o HVCFTestLD.o Main_TestAll.o

.PHONY: all blosclibs auxlibs applibs

all: blosclibs auxlibs applibs testAll

blosclibs:
	@for bloscdir in $(BLOSCDIRS); do \
		(cd $${bloscdir} && make -j 4) || exit 1; \
	done

auxlibs:
	@for auxdir in $(AUXDIRS); do \
		(cd $${auxdir} && make -j 4) || exit 1; \
	done
	
applibs:
	@for appdir in $(APPDIRS); do \
		(cd $${appdir} && make -j 4) || exit 1; \
	done

testAll: $(OBJECTS)
	$(CXX) $(CXXFLAGS) $(INCS) -o $@ $(OBJECTS) $(BLOSCLIBS) $(AUXLIBS) $(APPLIBS) $(LIBS)
	
.cpp.o:
	$(CXX) $(CXXFLAGS) $(INCS) -c -o $@ $<
	
clean:
	rm -f testAll *.o  ../src/*.o ../src/blosc/*.o
