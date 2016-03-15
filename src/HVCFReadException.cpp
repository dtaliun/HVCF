#include "include/HVCFReadException.h"

namespace sph_umich_edu {

HVCFReadException::HVCFReadException(
		const char* source_file, const char* function, unsigned int line, const char* message) :
				HVCFException(source_file, function, line, message) {

}

HVCFReadException::~HVCFReadException() {

}

}
