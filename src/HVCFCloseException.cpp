#include "include/HVCFCloseException.h"

namespace sph_umich_edu {

HVCFCloseException::HVCFCloseException(
		const char* source_file, const char* function, unsigned int line, const char* message) :
				HVCFException(source_file, function, line, message) {

}

HVCFCloseException::~HVCFCloseException() {

}

}
