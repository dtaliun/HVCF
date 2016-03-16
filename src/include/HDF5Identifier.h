#ifndef SRC_HDF5IDENTIFIER_H_
#define SRC_HDF5IDENTIFIER_H_

#include <iostream>
#include <limits>
#include "hdf5.h"
#include "HVCFException.h"

using namespace std;

namespace sph_umich_edu {

class HDF5Identifier {
protected:
	hid_t identifier;

public:
	HDF5Identifier();
	virtual ~HDF5Identifier() throw (HVCFException);

	void set(hid_t identifier);
	hid_t get();
	virtual void close() throw (HVCFException) = 0;
	hid_t release();

};

}

#endif
