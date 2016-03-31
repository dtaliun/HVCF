#ifndef HVCFCREATEEXCEPTION_H_
#define HVCFCREATEEXCEPTION_H_

#include "HVCFException.h"

using namespace std;

namespace sph_umich_edu {

class HVCFCreateException : public HVCFException {
public:
	HVCFCreateException(const char* source_file, const char* function, unsigned int line, const char* message);
	virtual ~HVCFCreateException();
};

}

#endif
