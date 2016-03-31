HDF5LIB=/Users/dtaliun/Documents/hdf5-1.10.0/hdf5/lib
HDF5INCS=/Users/dtaliun/Documents/hdf5-1.10.0/hdf5/include
BLOSCLIB=/Users/dtaliun/Documents/c-blosc-1.7.1/blosc
BLOSCINCS=/Users/dtaliun/Documents/c-blosc-1.7.1/blosc

INCS = -I$(HDF5INCS) -I$(BLOSCINCS)

all: $(patsubst %.c,%.o,$(wildcard *.c))

.c.o:
	gcc -Wall -O3 $(INCS) -c -o $@ $<
	
clean:
	rm -f *.o
