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

	HDF5Identifier(const HDF5Identifier& hdf5identifier) = delete;
	HDF5Identifier& operator=(const HDF5Identifier& hdf5identifier) = delete;
	HDF5Identifier(HDF5Identifier&& hdf5identifier) = delete;
	HDF5Identifier& operator=(HDF5Identifier&& hdf5identifier) = delete;

	hid_t operator=(hid_t identifier); //for assignment: avoids using set() every time.
	operator hid_t() const; //for implicit/explicit conversion: avoids using get() every time.

	void set(hid_t identifier);
	hid_t get() const;
	virtual void close() throw (HVCFException) = 0;
	hid_t release();

};

}

#endif
