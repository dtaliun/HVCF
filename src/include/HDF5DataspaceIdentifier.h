#ifndef SRC_HDF5DATASPACEIDENTIFIER_H_
#define SRC_HDF5DATASPACEIDENTIFIER_H_

#include "HDF5Identifier.h"

using namespace std;

namespace sph_umich_edu {

class HDF5DataspaceIdentifier: public HDF5Identifier {
public:
	HDF5DataspaceIdentifier();
	virtual ~HDF5DataspaceIdentifier() noexcept;

	using HDF5Identifier::operator=;

	void close() throw (HVCFException);
};

}

#endif
