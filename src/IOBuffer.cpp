#include "include/IOBuffer.h"

namespace sph_umich_edu {

IOBuffer::IOBuffer(unsigned int max_variants, unsigned int n_samples):
		max_variants(max_variants),
		n_samples(n_samples),
		n_haplotypes(n_samples + n_samples),
		haplotypes(nullptr),
		n_variants(0u) {

	haplotypes = unique_ptr<unsigned char[]>(new unsigned char[n_haplotypes * max_variants]);

}

IOBuffer::~IOBuffer() {
	cout << "Destroy IO BUffer!" << endl;
}

void IOBuffer::add_variant(const Variant& variant) throw (HVCFWriteException) {
	if (n_variants >= max_variants) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Memory buffer overflow while writing.");
	}

	for (unsigned int s = 0u; s < variant.get_n_samples(); ++s) {
		if (variant.get_genotype(s).get_alleles().size() != 2) {
			throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while writing variant to memory buffer.");
		}
		haplotypes[n_variants * n_haplotypes + 2u * s] = static_cast<unsigned char>(variant.get_genotype(s).get_alleles().at(0));
		haplotypes[n_variants * n_haplotypes + 2u * s + 1u] = static_cast<unsigned char>(variant.get_genotype(s).get_alleles().at(1));
	}

	++n_variants;
}

void IOBuffer::reset() {
	n_variants = 0u;
}

unsigned int IOBuffer::get_max_variants() const {
	return max_variants;
}

unsigned int IOBuffer::get_n_variants() const {
	return n_variants;
}

unsigned int IOBuffer::get_n_haplotypes() const {
	return n_haplotypes;
}

const unsigned char* IOBuffer::get_haplotypes_buffer() const {
	return haplotypes.get();
}

bool IOBuffer::is_full() const {
	return (n_variants >= max_variants);
}

}
