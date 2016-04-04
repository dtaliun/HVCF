#ifndef SRC_HDF5GROUPIDENTIFIER_H_
#define SRC_HDF5GROUPIDENTIFIER_H_

#include "HDF5Identifier.h"

using namespace std;

namespace sph_umich_edu {

class HDF5GroupIdentifier: public HDF5Identifier {
public:
	HDF5GroupIdentifier();
	virtual ~HDF5GroupIdentifier() noexcept;

	using HDF5Identifier::operator=;

	void close() throw (HVCFException);
};

}

#endif
