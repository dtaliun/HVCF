#include "include/HVCFCreateException.h"

namespace sph_umich_edu {

HVCFCreateException::HVCFCreateException(
		const char* source_file, const char* function, unsigned int line, const char* message) :
				HVCFException(source_file, function, line, message) {

}

HVCFCreateException::~HVCFCreateException() {

}

}
