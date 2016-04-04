#include "include/HDF5DatasetIdentifier.h"

namespace sph_umich_edu {

HDF5DatasetIdentifier::HDF5DatasetIdentifier() {

}

HDF5DatasetIdentifier::~HDF5DatasetIdentifier() noexcept {
	try {
		close();
	} catch (std::exception &e) {
		// do not propagate any exceptions
	}
}

void HDF5DatasetIdentifier::close() throw (HVCFException) {
	if (identifier >= 0) {
		if (H5Dclose(identifier) < 0) {
			throw HVCFException(__FILE__, __FUNCTION__, __LINE__, "Error while closing HDF5 dataset identifier.");
		}
		this->identifier = numeric_limits<hid_t>::min();
	}
}

}
