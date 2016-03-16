#ifndef SRC_HDF5PROPERTYIDENTIFIER_H_
#define SRC_HDF5PROPERTYIDENTIFIER_H_

#include "HDF5Identifier.h"

using namespace std;

namespace sph_umich_edu {

class HDF5PropertyIdentifier: public HDF5Identifier {
public:
	HDF5PropertyIdentifier();
	virtual ~HDF5PropertyIdentifier() throw (HVCFException);

	using HDF5Identifier::operator=;

	void close() throw (HVCFException);
};

}

#endif
