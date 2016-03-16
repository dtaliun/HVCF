#include "include/HDF5FileIdentifier.h"

namespace sph_umich_edu {

HDF5FileIdentifier::HDF5FileIdentifier() {

}

HDF5FileIdentifier::~HDF5FileIdentifier() {
	cout << "Closing file identifier" << endl;
	close();
}

void HDF5FileIdentifier::close() throw (HVCFException) {
	if (identifier >= 0) {
		if (H5Fclose(identifier) < 0) {
			throw HVCFException(__FILE__, __FUNCTION__, __LINE__, "Error while closing HDF5 file identifier.");
		}
		this->identifier = numeric_limits<hid_t>::min();
		cout << "closed" << endl;
	}
}

}
