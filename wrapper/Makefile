HDF5LIB=/Users/dtaliun/Documents/hdf5-1.10.0/hdf5/lib
HDF5INCS=/Users/dtaliun/Documents/hdf5-1.10.0/hdf5/include
BLOSCLIB=/Users/dtaliun/Documents/c-blosc-1.7.1/blosc
BLOSCINCS=/Users/dtaliun/Documents/c-blosc-1.7.1/blosc
BOOSTPYTHONINCS=/Users/dtaliun/Documents/boost_1_60_0/local_build/include/
BOOSTPYTHONLIB=/Users/dtaliun/Documents/boost_1_60_0/local_build/lib/
PYTHONINCS=/System/Library/Frameworks/Python.framework/Versions/2.7/include/python2.7/

AUXLIBS = ../../auxc/FileReader/src/*.o \
		../../auxc/MiniVCF/src/*.o
AUXDIRS = ../../auxc/FileReader/src \
		../../auxc/MiniVCF/src
	
BLOSCLIBS = ../src/blosc/*.o
BLOSCDIRS = ../src/blosc	
	
APPLIBS = ../src/*.o
APPDIRS = ../src

CXX = g++
CXXFLAGS = -std=c++11 -O3 -Wall -L$(BOOSTPYTHONLIB) -L$(HDF5LIB) -L$(BLOSCLIB)
INCS = -I$(PYTHONINCS) -I$(BOOSTPYTHONINCS) -I$(HDF5INCS) -I$(BLOSCINCS)
LIBS = -lz -lhdf5 -lblosc -larmadillo -lboost_python -lpython2.7

OBJECTS = PyHVCF.o

.PHONY: all blosclibs auxlibs applibs

all: blosclibs auxlibs applibs PyHVCF.so

PyHVCF.so: $(OBJECTS)
	$(CXX) $(CXXFLAGS) -shared -fPIC -o $@ $(OBJECTS) $(BLOSCLIBS) $(AUXLIBS) $(APPLIBS) $(LIBS)
	
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

.cpp.o:
	$(CXX) $(CXXFLAGS) $(INCS) -c -o $@ $<
		
clean:
	rm -f *.o *.so
