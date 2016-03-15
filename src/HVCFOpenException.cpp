#include "include/HVCFOpenException.h"

namespace sph_umich_edu {

HVCFOpenException::HVCFOpenException(
		const char* source_file, const char* function, unsigned int line, const char* message) :
				HVCFException(source_file, function, line, message) {

}

HVCFOpenException::~HVCFOpenException() {

}

}
