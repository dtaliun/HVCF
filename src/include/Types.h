#ifndef SRC_INCLUDE_TYPES_H_
#define SRC_INCLUDE_TYPES_H_

#include "hdf5.h"

namespace sph_umich_edu {

typedef struct {
	unsigned long long int ull_value_1;
	unsigned long long int ull_value_2;
	hsize_t offset_1;
	hsize_t offset_2;
	hobj_ref_t obj_reference;
} interval_index_key_type;

typedef struct {
	char* string_value;
	hsize_t offset;
} string_index_key_type;

typedef struct {
	unsigned long long int ull_value;
	hsize_t offset;
} ull_index_key_type;

typedef struct {
	char* name;
	unsigned long long int position;
} variants_entry_type;

}

#endif
