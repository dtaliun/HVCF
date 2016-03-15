#include "include/HVCFException.h"

namespace sph_umich_edu {

HVCFException::HVCFException(
		const char* source_file, const char* function, unsigned int line, const char* message):
				source_file(source_file), function(function), line(line), message(message) {

}

HVCFException::~HVCFException() {

}

const char* HVCFException::what() const noexcept {
	return message.c_str();
}

const string& HVCFException::get_source_file_name() {
	return source_file;
}

const string& HVCFException::get_function_name() {
	return function;
}

unsigned int HVCFException::get_line() {
	return line;
}

const string& HVCFException::get_message() {
	return message;
}

}
