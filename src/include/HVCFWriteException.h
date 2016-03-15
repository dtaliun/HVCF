#ifndef HVCFWRITEEXCEPTION_H_
#define HVCFWRITEEXCEPTION_H_

#include "HVCFException.h"

using namespace std;

namespace sph_umich_edu {

class HVCFWriteException : public HVCFException {
public:
	HVCFWriteException(const char* source_file, const char* function, unsigned int line, const char* message);
	virtual ~HVCFWriteException();
};

}

#endif
