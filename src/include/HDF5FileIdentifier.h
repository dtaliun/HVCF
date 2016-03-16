#ifndef SRC_HDF5FILEIDENTIFIER_H_
#define SRC_HDF5FILEIDENTIFIER_H_

#include "HDF5Identifier.h"

using namespace std;

namespace sph_umich_edu {

class HDF5FileIdentifier: public HDF5Identifier {
public:
	HDF5FileIdentifier();
	virtual ~HDF5FileIdentifier() throw (HVCFException);

	void close() throw (HVCFException);
};

}

#endif
