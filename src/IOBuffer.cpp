#include "include/IOBuffer.h"

namespace sph_umich_edu {

IOBuffer::IOBuffer(unsigned int max_variants, unsigned int n_samples):
		max_variants(max_variants),
		n_samples(n_samples),
		n_haplotypes(n_samples + n_samples),
		haplotypes(nullptr),
		names(nullptr),
		n_variants(0u) {

	haplotypes = unique_ptr<unsigned char[]>(new unsigned char[n_haplotypes * max_variants]{});
	names = unique_ptr<char*[]>(new char*[max_variants]{nullptr});
}

IOBuffer::~IOBuffer() {
	cout << "Destroy IO BUffer!" << endl;
	for (unsigned int i = 0u; i < max_variants; ++i) {
		if (names[i] != nullptr) {
			free(names[i]);
			names[i] = nullptr;
		}
	}
}

void IOBuffer::add_variant(const Variant& variant) throw (HVCFWriteException) {
	if (n_variants >= max_variants) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Memory buffer overflow while writing.");
	}

	if (variant.get_alt().get_values().size() != 1) { // Support only bi-allelic (for computing LD it is file, but must be extended).
		cout << "Warning!" << endl;
		return;
	}

	for (unsigned int s = 0u; s < variant.get_n_samples(); ++s) {
		if (variant.get_genotype(s).get_alleles().size() != 2) { // Support only HUMAN chromosomes 1-22 (should be extened for special case of chr Y).
			throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while writing variant to memory buffer.");
		}
		haplotypes[n_variants * n_haplotypes + 2u * s] = static_cast<unsigned char>(variant.get_genotype(s).get_alleles().at(0));
		haplotypes[n_variants * n_haplotypes + 2u * s + 1u] = static_cast<unsigned char>(variant.get_genotype(s).get_alleles().at(1));
	}

	unique_ptr<char[]> name = unique_ptr<char[]>(
			new char[variant.get_chrom().get_text().length() +
					 variant.get_pos().get_text().length() +
					 variant.get_ref().get_text().length() +
					 variant.get_alt().get_text().length() + 4u]{});

	name[0] = '\0';
	strcat(name.get(), variant.get_chrom().get_text().c_str());
	strcat(name.get(), ":");
	strcat(name.get(), variant.get_pos().get_text().c_str());
	strcat(name.get(), "_");
	strcat(name.get(), variant.get_ref().get_text().c_str());
	strcat(name.get(), "/");
	strcat(name.get(), variant.get_alt().get_text().c_str());

	names[n_variants] = name.release();

	++n_variants;
}

void IOBuffer::reset() {
	for (unsigned int i = 0u; i < n_variants; ++i) {
		if (names[i] != nullptr) {
			free(names[i]);
			names[i] = nullptr;
		}
	}

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

char* const* IOBuffer::get_names_buffer() const {
	return names.get();
}

bool IOBuffer::is_full() const {
	return (n_variants >= max_variants);
}

bool IOBuffer::is_empty() const {
	return (n_variants == 0u);
}

}
