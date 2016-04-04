#ifndef SRC_INCLUDE_HVCFCONFIGURATION_H_
#define SRC_INCLUDE_HVCFCONFIGURATION_H_

#include "hdf5.h"

namespace sph_umich_edu {

class HVCFConfiguration {
public:
	static constexpr char GZIP_COMPRESSION[] = "GZIP";
	static constexpr char BLOSC_LZ4HC_COMPRESSION[] = "BLOSC_LZ4HC";

	unsigned int n_hash_buckets;
	unsigned int max_variants_in_interval_bucket;
	hsize_t variants_chunk_size;
	hsize_t samples_chunk_size;
	const char* compression;
	unsigned int compression_level;

	HVCFConfiguration();
	virtual ~HVCFConfiguration();
};

}

#endif
