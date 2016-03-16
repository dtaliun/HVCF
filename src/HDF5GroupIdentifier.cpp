#include "include/HDF5GroupIdentifier.h"

namespace sph_umich_edu {

HDF5GroupIdentifier::HDF5GroupIdentifier() {

}

HDF5GroupIdentifier::~HDF5GroupIdentifier() {
	close();
}

void HDF5GroupIdentifier::close() throw (HVCFException) {
	if (identifier >= 0) {
		if (H5Gclose(identifier) < 0) {
			throw HVCFException(__FILE__, __FUNCTION__, __LINE__, "Error while closing HDF5 group identifier.");
		}
		this->identifier = numeric_limits<hid_t>::min();
	}
}

}
