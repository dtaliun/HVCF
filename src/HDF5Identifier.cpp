#include "include/HDF5Identifier.h"

namespace sph_umich_edu {

HDF5Identifier::HDF5Identifier() : identifier(numeric_limits<hid_t>::min()) {

}

HDF5Identifier::~HDF5Identifier() throw (HVCFException) {

}

void HDF5Identifier::set(hid_t identifier) {
	this->identifier = identifier;
}

hid_t HDF5Identifier::get() {
	return identifier;
}

hid_t HDF5Identifier::release() {
	hid_t identifier = this->identifier;
	this->identifier = numeric_limits<hid_t>::min();
	return identifier;
}

}
