HDF5LIB=/Users/dtaliun/Documents/hdf5-1.10.0-alpha1/hdf5/lib
HDF5INCS=/Users/dtaliun/Documents/hdf5-1.10.0-alpha1/hdf5/include

CXX = g++
CXXFLAGS = -std=c++11 -Wall -O3 -L$(HDF5LIB)
INCS = -I$(HDF5INCS)

all: HVCF.o HVCFException.o HVCFOpenException.o HVCFCloseException.o HVCFWriteException.o

.cpp.o:
	$(CXX) $(CXXFLAGS) $(INCS) -c -o $@ $<
	
clean: