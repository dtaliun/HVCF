#include "include/HDF5PropertyIdentifier.h"

namespace sph_umich_edu {

HDF5PropertyIdentifier::HDF5PropertyIdentifier() {

}

HDF5PropertyIdentifier::~HDF5PropertyIdentifier() noexcept {
	try {
		close();
	} catch (std::exception &e) {
		// do not propagate any exceptions
	}
}

void HDF5PropertyIdentifier::close() throw (HVCFException) {
	if (identifier >= 0) {
		if (H5Pclose(identifier) < 0) {
			throw HVCFException(__FILE__, __FUNCTION__, __LINE__, "Error while closing HDF5 property identifier.");
		}
		this->identifier = numeric_limits<hid_t>::min();
	}
}

}
