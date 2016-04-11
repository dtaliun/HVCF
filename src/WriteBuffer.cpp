#include "include/WriteBuffer.h"

namespace sph_umich_edu {

WriteBuffer::WriteBuffer(unsigned int max_variants, unsigned int n_samples):
		max_variants(max_variants),
		n_samples(n_samples),
		n_haplotypes(n_samples + n_samples),
		haplotypes(nullptr),
		variants(nullptr),
		n_variants(0u),
		n_flushed_variants(0u) {

	haplotypes = unique_ptr<unsigned char[]>(new unsigned char[n_haplotypes * max_variants]{});
	variants = unique_ptr<variants_entry_type[]>(new variants_entry_type[max_variants]{});
	for (unsigned int i = 0u; i < max_variants; ++i) {
		variants[i].name = nullptr;
		variants[i].ref = nullptr;
		variants[i].alt = nullptr;
	}

	flushed_haplotypes = unique_ptr<unsigned char[]>(new unsigned char[n_haplotypes * max_variants]{});
	flushed_variants = unique_ptr<variants_entry_type[]>(new variants_entry_type[max_variants]{});
	for (unsigned int i = 0u; i < max_variants; ++i) {
		flushed_variants[i].name = nullptr;
		flushed_variants[i].ref = nullptr;
		flushed_variants[i].alt = nullptr;
	}
}

WriteBuffer::~WriteBuffer() {
	for (unsigned int i = 0u; i < max_variants; ++i) {
		if (variants[i].name != nullptr) {
			delete variants[i].name;
			variants[i].name = nullptr;
		}
		if (variants[i].ref != nullptr) {
			delete variants[i].ref;
			variants[i].ref = nullptr;
		}
		if (variants[i].alt != nullptr) {
			delete variants[i].alt;
			variants[i].alt = nullptr;
		}
		if (flushed_variants[i].name != nullptr) {
			delete flushed_variants[i].name;
			flushed_variants[i].name = nullptr;
		}
		if (flushed_variants[i].ref != nullptr) {
			delete flushed_variants[i].ref;
			flushed_variants[i].ref = nullptr;
		}
		if (flushed_variants[i].alt != nullptr) {
			delete flushed_variants[i].alt;
			flushed_variants[i].alt = nullptr;
		}
	}
}

void WriteBuffer::add_variant(const Variant& variant) throw (HVCFWriteException) {
	if (n_variants >= max_variants) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Memory buffer overflow while writing.");
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

	unique_ptr<char[]> ref = unique_ptr<char[]>(new char[variant.get_ref().get_text().length() + 1u]{});
	unique_ptr<char[]> alt = unique_ptr<char[]>(new char[variant.get_alt().get_text().length() + 1u]{});

	name[0] = '\0';
	strcat(name.get(), variant.get_chrom().get_text().c_str());
	strcat(name.get(), ":");
	strcat(name.get(), variant.get_pos().get_text().c_str());
	strcat(name.get(), "_");
	strcat(name.get(), variant.get_ref().get_text().c_str());
	strcat(name.get(), "/");
	strcat(name.get(), variant.get_alt().get_text().c_str());

	ref[0] = '\0';
	strcat(ref.get(), variant.get_ref().get_text().c_str());

	alt[0] = '\0';
	strcat(alt.get(), variant.get_alt().get_text().c_str());

	variants[n_variants].name = name.release();
	variants[n_variants].ref = ref.release();
	variants[n_variants].alt = alt.release();
	variants[n_variants].position = variant.get_pos().get_value();

	++n_variants;
}

tuple<const unsigned char*, const variants_entry_type*, unsigned int, unsigned int> WriteBuffer::flush() {
	for (unsigned int i = 0u; i < n_flushed_variants; ++i) {
		if (flushed_variants[i].name != nullptr) {
			delete flushed_variants[i].name;
			flushed_variants[i].name = nullptr;
		}
		if (flushed_variants[i].ref != nullptr) {
			delete flushed_variants[i].ref;
			flushed_variants[i].ref = nullptr;
		}
		if (flushed_variants[i].alt != nullptr) {
			delete flushed_variants[i].alt;
			flushed_variants[i].alt = nullptr;
		}
	}

	n_flushed_variants = n_variants;
	flushed_haplotypes.swap(haplotypes);
	flushed_variants.swap(variants);

	n_variants = 0u;

	return std::make_tuple(flushed_haplotypes.get(), flushed_variants.get(), n_flushed_variants, n_haplotypes);
}

unsigned int WriteBuffer::get_max_variants() const {
	return max_variants;
}

unsigned int WriteBuffer::get_n_samples() const {
	return n_samples;
}

bool WriteBuffer::is_full() const {
	return (n_variants >= max_variants);
}

bool WriteBuffer::is_empty() const {
	return (n_variants == 0u);
}

}
