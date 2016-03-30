HDF5LIB=/Users/dtaliun/Documents/hdf5-1.10.0-pre2/hdf5/lib
HDF5INCS=/Users/dtaliun/Documents/hdf5-1.10.0-pre2/hdf5/include
BLOSCLIB=/Users/dtaliun/Documents/c-blosc-1.7.1/blosc
BLOSCINCS=/Users/dtaliun/Documents/c-blosc-1.7.1/blosc

CXX = g++
CXXFLAGS = -std=c++11 -Wall -O3 -L$(HDF5LIB) -L$(BLOSCLIB)
INCS = -I$(HDF5INCS) -I$(BLOSCINCS)

all: HVCFException.o \
	HVCFOpenException.o \
	HVCFCloseException.o \
	HVCFWriteException.o \
	HVCFReadException.o \
	HDF5Identifier.o \
	HDF5FileIdentifier.o \
	HDF5GroupIdentifier.o \
	HDF5DatasetIdentifier.o \
	HDF5DataspaceIdentifier.o \
	HDF5DatatypeIdentifier.o \
	HDF5PropertyIdentifier.o \
	WriteBuffer.o \
	HVCF.o

.cpp.o:
	$(CXX) $(CXXFLAGS) $(INCS) -c -o $@ $<
	
clean: