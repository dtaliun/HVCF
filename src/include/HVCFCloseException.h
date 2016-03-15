#ifndef HVCFCLOSEEXCEPTION_H_
#define HVCFCLOSEEXCEPTION_H_

#include "HVCFException.h"

using namespace std;

namespace sph_umich_edu {

class HVCFCloseException : public HVCFException {
public:
	HVCFCloseException(const char* source_file, const char* function, unsigned int line, const char* message);
	virtual ~HVCFCloseException();
};

}

#endif
