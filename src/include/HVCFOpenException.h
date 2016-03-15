#ifndef HVCFOPENEXCEPTION_H_
#define HVCFOPENEXCEPTION_H_

#include "HVCFException.h"

using namespace std;

namespace sph_umich_edu {

class HVCFOpenException : public HVCFException {
public:
	HVCFOpenException(const char* source_file, const char* function, unsigned int line, const char* message);
	virtual ~HVCFOpenException();
};

}

#endif
