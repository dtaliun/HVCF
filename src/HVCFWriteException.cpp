#include "include/HVCFWriteException.h"

namespace sph_umich_edu {

HVCFWriteException::HVCFWriteException(
		const char* source_file, const char* function, unsigned int line, const char* message) :
				HVCFException(source_file, function, line, message) {

}

HVCFWriteException::~HVCFWriteException() {

}

}
