#ifndef SRC_INCLUDE_IOBUFFER_H_
#define SRC_INCLUDE_IOBUFFER_H_

#include <iostream>
#include <memory>

#include "../../../auxc/MiniVCF/src/include/VCFReader.h"
#include "HVCFWriteException.h"

using namespace std;

namespace sph_umich_edu {

class IOBuffer {
private:
	unsigned int max_variants;
	unsigned int n_samples;
	unsigned int n_haplotypes;

	unique_ptr<unsigned char[]> haplotypes;
	unique_ptr<char*[]> names;

	unsigned int n_variants;

public:
	IOBuffer(unsigned int max_variants, unsigned int n_samples);
	virtual ~IOBuffer();

	void add_variant(const Variant& variant) throw (HVCFWriteException);
	void reset();

	unsigned int get_max_variants() const;
	unsigned int get_n_variants() const;
	unsigned int get_n_haplotypes() const;
	const unsigned char* get_haplotypes_buffer() const;
	char* const* get_names_buffer() const;
	bool is_full() const;
	bool is_empty() const;
};

}

#endif
