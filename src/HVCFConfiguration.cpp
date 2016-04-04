#include "include/HVCFConfiguration.h"

namespace sph_umich_edu {

constexpr char HVCFConfiguration::GZIP_COMPRESSION[];
constexpr char HVCFConfiguration::BLOSC_LZ4HC_COMPRESSION[];

HVCFConfiguration::HVCFConfiguration() {
	n_hash_buckets = 10000;
	max_variants_in_interval_bucket = 1000;
	variants_chunk_size = 500;
	samples_chunk_size = 1000;
	compression = HVCFConfiguration::GZIP_COMPRESSION;
//	compression = HVCFConfiguration::BLOSC_LZ4HC_COMPRESSION;
	compression_level = 9;
}

HVCFConfiguration::~HVCFConfiguration() {

}

}
