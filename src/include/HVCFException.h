#ifndef SRC_HVCFEXCEPTION_H_
#define SRC_HVCFEXCEPTION_H_

#include <iostream>
#include <sstream>
#include <exception>
#include <stdexcept>

using namespace std;

namespace sph_umich_edu {

class HVCFException : public std::exception {
protected:
	string source_file;
	string function;
	unsigned int line;
	string message;

public:
	HVCFException(const char* source_file, const char* function, unsigned int line, const char* message);
	virtual ~HVCFException();

	virtual const char* what() const noexcept;

	const string& get_source_file_name();
	const string& get_function_name();
	unsigned int get_line();
	const string& get_message();
};

}

#endif
