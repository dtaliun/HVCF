#ifndef SRC_INCLUDE_TYPES_H_
#define SRC_INCLUDE_TYPES_H_

#include <string>
#include "hdf5.h"

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
	unsigned long long int position;
} variants_entry_type;

typedef struct {
	char* name;
	hsize_t offset_1;
	hsize_t offset_2;
} subsets_entry_type;

typedef struct VariantsPair{
	string name1;
	unsigned long long int position1;
	string name2;
	unsigned long long int position2;
	double r;
	double rsquare;

	VariantsPair(const char* name1, unsigned long int position1,
			const char* name2, unsigned long int position2, double r, double rsquare) :
				name1(name1), position1(position1), name2(name2), position2(position2), r(r), rsquare(rsquare) {

	}

	bool operator==(VariantsPair const& pair) const { // needed for boost.python
		return (position1 == pair.position1 && position2 == pair.position2);
	}

} variants_pair;

}

#endif
