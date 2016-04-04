#include "include/HDF5DataspaceIdentifier.h"

namespace sph_umich_edu {

HDF5DataspaceIdentifier::HDF5DataspaceIdentifier() {

}

HDF5DataspaceIdentifier::~HDF5DataspaceIdentifier() noexcept {
	try {
		close();
	} catch (std::exception &e) {
		// do not propagate any exceptions
	}
}

void HDF5DataspaceIdentifier::close() throw (HVCFException) {
	if (identifier >= 0) {
		if (H5Sclose(identifier) < 0) {
			throw HVCFException(__FILE__, __FUNCTION__, __LINE__, "Error while closing HDF5 dataspace identifier.");
		}
		this->identifier = numeric_limits<hid_t>::min();
	}
}

}
