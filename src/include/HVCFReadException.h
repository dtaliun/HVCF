#ifndef HVCFREADEXCEPTION_H_
#define HVCFREADEXCEPTION_H_

#include "HVCFException.h"

using namespace std;

namespace sph_umich_edu {

class HVCFReadException : public HVCFException {
public:
	HVCFReadException(const char* source_file, const char* function, unsigned int line, const char* message);
	virtual ~HVCFReadException();
};

}

#endif
