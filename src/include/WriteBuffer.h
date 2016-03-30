#ifndef SRC_INCLUDE_WRITEBUFFER_H_
#define SRC_INCLUDE_WRITEBUFFER_H_

#include <iostream>
#include <memory>

#include "../../../auxc/MiniVCF/src/include/VCFReader.h"
#include "HVCFWriteException.h"
#include "Types.h"

using namespace std;

namespace sph_umich_edu {

class WriteBuffer {
private:
	unsigned int max_variants;
	unsigned int n_samples;
	unsigned int n_haplotypes;

	unique_ptr<unsigned char[]> haplotypes;
	unique_ptr<variants_entry_type[]> variants;

	unsigned int n_variants;

public:
	WriteBuffer(unsigned int max_variants, unsigned int n_samples);
	virtual ~WriteBuffer();

	void add_variant(const Variant& variant) throw (HVCFWriteException);
	void reset();

	unsigned int get_max_variants() const;
	unsigned int get_n_variants() const;
	unsigned int get_n_haplotypes() const;
	const unsigned char* get_haplotypes_buffer() const;
	const variants_entry_type* get_variants_buffer() const;
	bool is_full() const;
	bool is_empty() const;
};

}

#endif
