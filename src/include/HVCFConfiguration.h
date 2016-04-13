#ifndef SRC_INCLUDE_HVCFCONFIGURATION_H_
#define SRC_INCLUDE_HVCFCONFIGURATION_H_

#include "hdf5.h"

namespace sph_umich_edu {

class HVCFConfiguration {
public:
	static constexpr char GZIP_COMPRESSION[] = "GZIP";
	static constexpr char BLOSC_LZ4HC_COMPRESSION[] = "BLOSC_LZ4HC";

	unsigned int n_variants_hash_buckets;
	unsigned int n_samples_hash_buckets;
	unsigned int max_variants_in_interval_bucket;
	hsize_t variants_chunk_size;
	hsize_t samples_chunk_size;
	const char* compression;
	unsigned int compression_level;
	size_t metadata_cache_initial_size;
	size_t metadata_cache_min_size;
	size_t metadata_cache_max_size;
	size_t sieve_buffer_max_size;
	size_t chunk_cache_n_slots;
	size_t chunk_cache_size;

	HVCFConfiguration();
	virtual ~HVCFConfiguration();
};

}

#endif
