#ifndef SRC_HDF5DATASETIDENTIFIER_H_
#define SRC_HDF5DATASETIDENTIFIER_H_

#include "HDF5Identifier.h"

using namespace std;

namespace sph_umich_edu {

class HDF5DatasetIdentifier: public HDF5Identifier {
public:
	HDF5DatasetIdentifier();
	virtual ~HDF5DatasetIdentifier() throw (HVCFException);

	using HDF5Identifier::operator=;

	void close() throw (HVCFException);
};

}

#endif
