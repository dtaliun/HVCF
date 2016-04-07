#include "include/HVCFConfiguration.h"

namespace sph_umich_edu {

constexpr char HVCFConfiguration::GZIP_COMPRESSION[];
constexpr char HVCFConfiguration::BLOSC_LZ4HC_COMPRESSION[];

HVCFConfiguration::HVCFConfiguration() {
	n_hash_buckets = 10000;
	max_variants_in_interval_bucket = 1000;
	variants_chunk_size = 1000;
	samples_chunk_size = 100;
	compression = HVCFConfiguration::GZIP_COMPRESSION;
//	compression = HVCFConfiguration::BLOSC_LZ4HC_COMPRESSION;
	compression_level = 9;
	metadata_cache_initial_size = 64 * 1024 * 1024;
	metadata_cache_min_size = 8 * 1024 * 1024;
	metadata_cache_max_size = 128 * 1024 * 1024;
	sieve_buffer_max_size = 64 * 1024 * 1024; // assuming max hyperslab = 10000 (variants) * 5000 (variants) * sizeof(char)
	chunk_cache_n_slots = 10000; // following HDF5 documentation, 100x more than max number of chunks in cache (i.e. 100 * 100)
	chunk_cache_size = 32 * 1024 * 1024; // if we want to hold up to 100 chunks in cache for haplotypes (i.e. not less than 100 * variants_chunk_size * 2 * samples_chunk_size * sizeof(char))
}

HVCFConfiguration::~HVCFConfiguration() {

}

}
