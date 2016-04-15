#ifndef SRC_INCLUDE_TYPES_H_
#define SRC_INCLUDE_TYPES_H_

#include <string>
#include <map>
#include "hdf5.h"

#include "HDF5DatasetIdentifier.h"

using namespace std;

namespace sph_umich_edu {

typedef struct {
	unsigned long long int ull_value_1;
	unsigned long long int ull_value_2;
	hsize_t offset_1;
	hsize_t offset_2;
	hsize_t bucket_offset;
	hsize_t bucket_size;
} interval_index_entry_type;

typedef struct {
	hsize_t bucket_offset;
	hsize_t bucket_size;
} hash_index_entry_type;

typedef struct {
	char* string_value;
	hsize_t offset;
} string_index_entry_type;

typedef struct {
	unsigned long long int ull_value;
	hsize_t offset;
} ull_index_entry_type;

typedef struct {
	char* name;
	char* ref;
	char* alt;
	unsigned long long int position;
} variants_entry_type;

typedef struct {
	char* name;
	hsize_t offset_1;
	hsize_t offset_2;
} subsets_entry_type;

typedef struct {
	vector<tuple<hsize_t, hsize_t, hsize_t>> chunks; // offset_1 (start), offset_2 (end), size (offset_2 - offset_1 + 1)
	hsize_t n_samples;
} subsets_cache_entry;

typedef struct {
	HDF5DatasetIdentifier names_index_id;
	HDF5DatasetIdentifier names_index_buckets_id;
	unordered_map<string, subsets_cache_entry> subsets;
} samples_cache_entry;

typedef struct {
	HDF5DatasetIdentifier names_index_id;
	HDF5DatasetIdentifier names_index_buckets_id;
	HDF5DatasetIdentifier intervals_index_id;
	HDF5DatasetIdentifier intervals_index_buckets_id;
	HDF5DatasetIdentifier variants_id;
	HDF5DatasetIdentifier haplotypes_id;
} chromosomes_cache_entry;

typedef struct VariantQueryResult {
	string name;
	string ref;
	string alt;
	unsigned long long int position;

	VariantQueryResult(const char* name, const char* ref, const char* alt, unsigned long long int position) :
		name(name), ref(ref), alt(alt), position(position) {
	}

	bool operator==(VariantQueryResult const& result) const { // needed for boost.python
		return (position == result.position && position == result.position && name.compare(result.name) == 0);
	}
} variant_query_result;

typedef struct FrequencyQueryResult {
	string name;
	string ref;
	string alt;
	unsigned long long int position;
	double ref_af;
	double alt_af;

	FrequencyQueryResult(const char* name, const char* ref, const char* alt, unsigned long long int position, double ref_af, double alt_af) :
		name(name), ref(ref), alt(alt), position(position), ref_af(ref_af), alt_af(alt_af) {
	}

	bool operator==(FrequencyQueryResult const& result) const { // needed for boost.python
		return (position == result.position && position == result.position && name.compare(result.name) == 0);
	}
} frequency_query_result;

typedef struct LDQueryResult{
	string name1;
	unsigned long long int position1;
	string name2;
	unsigned long long int position2;
	double r;
	double rsquare;

	LDQueryResult() : name1(""), position1(0ul), name2(""), position2(0ul), r(0.0), rsquare(0.0) {

	}

	LDQueryResult(const char* name1, unsigned long int position1,
			const char* name2, unsigned long int position2, double r, double rsquare) :
				name1(name1), position1(position1), name2(name2), position2(position2), r(r), rsquare(rsquare) {

	}

	bool operator==(LDQueryResult const& result) const { // needed for boost.python
		return (position1 == result.position1 && position2 == result.position2 &&
				name1.compare(result.name1) == 0 && name2.compare(result.name2) == 0);
	}

} ld_query_result;

typedef struct VariantHaplotypesQueryResult {
	string sample;
	unsigned char allele1;
	unsigned char allele2;

	VariantHaplotypesQueryResult() : sample(""), allele1(0), allele2(0) {

	}

	VariantHaplotypesQueryResult(const char* sample, unsigned char allele1, unsigned char allele2) :
		sample(sample), allele1(allele1), allele2(allele2) {
	}

	bool operator==(VariantHaplotypesQueryResult const& result) const { // needed for boost.python
		return (sample.compare(result.sample) == 0);
	}
} variant_haplotypes_query_result;

typedef struct SampleHaplotypesQueryResult {
	string name;
	unsigned long long position;
	unsigned char allele1;
	unsigned char allele2;

	SampleHaplotypesQueryResult() : name(""), position(0ul), allele1(0), allele2(0) {

	}

	SampleHaplotypesQueryResult(const char* name, unsigned long long int position, unsigned char allele1, unsigned char allele2) :
		name(name), position(position), allele1(allele1), allele2(allele2) {

	}

	bool operator==(SampleHaplotypesQueryResult const& result) const { // needed for boost.python
		return (position == result.position && name.compare(result.name) == 0);
	}

} sample_haplotypes_query_result;

}

#endif
