HDF5LIB=/Users/dtaliun/Documents/hdf5-1.10.0-alpha1/hdf5/lib
HDF5INCS=/Users/dtaliun/Documents/hdf5-1.10.0-alpha1/hdf5/include

CXX = g++
CXXFLAGS = -std=c++11 -Wall -O3 -L$(HDF5LIB)
INCS = -I$(HDF5INCS)

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
	IOBuffer.o \
	HVCF.o

.cpp.o:
	$(CXX) $(CXXFLAGS) $(INCS) -c -o $@ $<
	
clean: