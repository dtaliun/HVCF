#include "include/HVCF.h"

namespace sph_umich_edu {

constexpr char HVCF::CHROMOSOMES_GROUP[];
constexpr char HVCF::SAMPLES_GROUP[];
constexpr char HVCF::VARIANTS_DATASET[];
constexpr char HVCF::HAPLOTYPES_DATASET[];
constexpr char HVCF::SAMPLE_NAMES_DATASET[];
constexpr char HVCF::SAMPLE_SUBSETS_DATASET[];
constexpr char HVCF::VARIABLE_LENGTH_STRING_TYPE[];
constexpr char HVCF::VARIANTS_ENTRY_TYPE[];
constexpr char HVCF::SUBSETS_ENTRY_TYPE[];
constexpr char HVCF::STRING_INDEX_ENTRY_TYPE[];
constexpr char HVCF::INTERVAL_INDEX_ENTRY_TYPE[];
constexpr char HVCF::HASH_INDEX_ENTRY_TYPE[];
constexpr char HVCF::ULL_INDEX_ENTRY_TYPE[];
constexpr char HVCF::NAMES_INDEX_GROUP[];
constexpr char HVCF::INTERVALS_INDEX_GROUP[];
constexpr char HVCF::INTERVALS_INDEX[];
constexpr char HVCF::HASH_INDEX[];
constexpr char HVCF::INDEX_BUCKETS[];


HVCF::HVCF() : HVCF(HVCFConfiguration()) {

}

HVCF::HVCF(const HVCFConfiguration& configuration) {
//	Disables automatic HDF5 error stack printing to stderr when function call returns negative value.
//	H5Eset_auto(H5E_DEFAULT, nullptr, nullptr);

	N_VARIANTS_HASH_BUCKETS = configuration.n_variants_hash_buckets;
	N_SAMPLES_HASH_BUCKETS = configuration.n_samples_hash_buckets;
	MAX_VARIANTS_IN_INTERVAL_BUCKET = configuration.max_variants_in_interval_bucket;
	VARIANTS_CHUNK_SIZE = configuration.variants_chunk_size;
	SAMPLES_CHUNK_SIZE = configuration.samples_chunk_size;
	COMPRESSION = configuration.compression;
	COMPRESSION_LEVEL = configuration.compression_level;
	METADATA_CACHE_INITIAL_SIZE = configuration.metadata_cache_initial_size;
	METADATA_CACHE_MIN_SIZE = configuration.metadata_cache_min_size;
	METADATA_CACHE_MAX_SIZE = configuration.metadata_cache_max_size;
	SIEVE_BUFFER_MAX_SIZE = configuration.sieve_buffer_max_size;
	CHUNK_CACHE_N_SLOTS = configuration.chunk_cache_n_slots;
	CHUNK_CACHE_SIZE = configuration.chunk_cache_size;

//  Register Blosc
	if (strcmp(COMPRESSION, HVCFConfiguration::BLOSC_LZ4HC_COMPRESSION) == 0) {
		char* blosc_version = nullptr;
		char* blosc_date = nullptr;
		if (register_blosc(&blosc_version, &blosc_date) < 0) {
			// throw some expression
		}
	}
}

HVCF::~HVCF() {

}

hid_t HVCF::create_variants_entry_memory_datatype() throw (HVCFCreateException) {
	HDF5DatatypeIdentifier native_string_datatype_id;
	HDF5DatatypeIdentifier memory_datatype_id;

	if ((native_string_datatype_id = H5Topen(file_id, VARIABLE_LENGTH_STRING_TYPE, H5P_DEFAULT)) < 0) {
		throw HVCFCreateException(__FILE__, __FUNCTION__, __LINE__, "Error while opening committed datatype.");
	}

	if ((memory_datatype_id = H5Tcreate(H5T_COMPOUND, sizeof(variants_entry_type))) < 0) {
		throw HVCFCreateException(__FILE__, __FUNCTION__, __LINE__, "Error while creating compound datatype.");
	}

	if (H5Tinsert(memory_datatype_id, "name", HOFFSET(variants_entry_type, name), native_string_datatype_id) < 0) {
		throw HVCFCreateException(__FILE__, __FUNCTION__, __LINE__, "Error while adding new member to compound datatype.");
	}

	if (H5Tinsert(memory_datatype_id, "ref", HOFFSET(variants_entry_type, ref), native_string_datatype_id) < 0) {
		throw HVCFCreateException(__FILE__, __FUNCTION__, __LINE__, "Error while adding new member to compound datatype.");
	}

	if (H5Tinsert(memory_datatype_id, "alt", HOFFSET(variants_entry_type, alt), native_string_datatype_id) < 0) {
		throw HVCFCreateException(__FILE__, __FUNCTION__, __LINE__, "Error while adding new member to compound datatype.");
	}

	if (H5Tinsert(memory_datatype_id, "position", HOFFSET(variants_entry_type, position), H5T_NATIVE_ULLONG) < 0) {
		throw HVCFCreateException(__FILE__, __FUNCTION__, __LINE__, "Error while adding new member to compound datatype.");
	}

	return memory_datatype_id.release();
}

hid_t HVCF::create_subsets_entry_memory_datatype() throw (HVCFCreateException) {
	HDF5DatatypeIdentifier native_string_datatype_id;
	HDF5DatatypeIdentifier memory_datatype_id;

	if ((native_string_datatype_id = H5Topen(file_id, VARIABLE_LENGTH_STRING_TYPE, H5P_DEFAULT)) < 0) {
		throw HVCFCreateException(__FILE__, __FUNCTION__, __LINE__, "Error while opening committed datatype.");
	}

	if ((memory_datatype_id = H5Tcreate(H5T_COMPOUND, sizeof(subsets_entry_type))) < 0) {
		throw HVCFCreateException(__FILE__, __FUNCTION__, __LINE__, "Error while creating compound datatype.");
	}

	if (H5Tinsert(memory_datatype_id, "name", HOFFSET(subsets_entry_type, name), native_string_datatype_id) < 0) {
		throw HVCFCreateException(__FILE__, __FUNCTION__, __LINE__, "Error while adding new member to compound datatype.");
	}

	if (H5Tinsert(memory_datatype_id, "offset_1", HOFFSET(subsets_entry_type, offset_1), H5T_NATIVE_HSIZE) < 0) {
		throw HVCFCreateException(__FILE__, __FUNCTION__, __LINE__, "Error while adding new member to compound datatype.");
	}

	if (H5Tinsert(memory_datatype_id, "offset_2", HOFFSET(subsets_entry_type, offset_2), H5T_NATIVE_HSIZE) < 0) {
		throw HVCFCreateException(__FILE__, __FUNCTION__, __LINE__, "Error while adding new member to compound datatype.");
	}

	return memory_datatype_id.release();
}

hid_t HVCF::create_ull_index_entry_memory_datatype() throw (HVCFCreateException) {
	HDF5DatatypeIdentifier memory_datatype_id;

	if ((memory_datatype_id = H5Tcreate(H5T_COMPOUND, sizeof(ull_index_entry_type))) < 0) {
		throw HVCFCreateException(__FILE__, __FUNCTION__, __LINE__, "Error while creating compound datatype.");
	}

	if (H5Tinsert(memory_datatype_id, "ull_value", HOFFSET(ull_index_entry_type, ull_value), H5T_NATIVE_ULLONG) < 0) {
		throw HVCFCreateException(__FILE__, __FUNCTION__, __LINE__, "Error while adding new member to compound datatype.");
	}

	if (H5Tinsert(memory_datatype_id, "offset", HOFFSET(ull_index_entry_type, offset), H5T_NATIVE_HSIZE) < 0) {
		throw HVCFCreateException(__FILE__, __FUNCTION__, __LINE__, "Error while adding new member to compound datatype.");
	}

	return memory_datatype_id.release();
}

hid_t HVCF::create_string_index_entry_memory_datatype() throw (HVCFCreateException) {
	HDF5DatatypeIdentifier native_string_datatype_id;
	HDF5DatatypeIdentifier memory_datatype_id;

	if ((native_string_datatype_id = H5Topen(file_id, VARIABLE_LENGTH_STRING_TYPE, H5P_DEFAULT)) < 0) {
		throw HVCFCreateException(__FILE__, __FUNCTION__, __LINE__, "Error while opening committed datatype.");
	}

	if ((memory_datatype_id = H5Tcreate(H5T_COMPOUND, sizeof(string_index_entry_type))) < 0) {
		throw HVCFCreateException(__FILE__, __FUNCTION__, __LINE__, "Error while creating compound datatype.");
	}

	if (H5Tinsert(memory_datatype_id, "string_value", HOFFSET(string_index_entry_type, string_value), native_string_datatype_id) < 0) {
		throw HVCFCreateException(__FILE__, __FUNCTION__, __LINE__, "Error while adding new member to compound datatype.");
	}

	if (H5Tinsert(memory_datatype_id, "offset", HOFFSET(string_index_entry_type, offset), H5T_NATIVE_HSIZE) < 0) {
		throw HVCFCreateException(__FILE__, __FUNCTION__, __LINE__, "Error while adding new member to compound datatype.");
	}

	return memory_datatype_id.release();
}

hid_t HVCF::create_interval_index_entry_memory_datatype() throw (HVCFCreateException) {
	HDF5DatatypeIdentifier memory_datatype_id;

	if ((memory_datatype_id = H5Tcreate(H5T_COMPOUND, sizeof(interval_index_entry_type))) < 0) {
		throw HVCFCreateException(__FILE__, __FUNCTION__, __LINE__, "Error while creating compound datatype.");
	}

	if (H5Tinsert(memory_datatype_id, "ull_value_1", HOFFSET(interval_index_entry_type, ull_value_1), H5T_NATIVE_ULLONG) < 0) {
		throw HVCFCreateException(__FILE__, __FUNCTION__, __LINE__, "Error while adding new member to compound datatype.");
	}

	if (H5Tinsert(memory_datatype_id, "ull_value_2", HOFFSET(interval_index_entry_type, ull_value_2), H5T_NATIVE_ULLONG) < 0) {
		throw HVCFCreateException(__FILE__, __FUNCTION__, __LINE__, "Error while adding new member to compound datatype.");
	}

	if (H5Tinsert(memory_datatype_id, "offset_1", HOFFSET(interval_index_entry_type, offset_1), H5T_NATIVE_HSIZE) < 0) {
		throw HVCFCreateException(__FILE__, __FUNCTION__, __LINE__, "Error while adding new member to compound datatype.");
	}

	if (H5Tinsert(memory_datatype_id, "offset_2", HOFFSET(interval_index_entry_type, offset_2), H5T_NATIVE_HSIZE) < 0) {
		throw HVCFCreateException(__FILE__, __FUNCTION__, __LINE__, "Error while adding new member to compound datatype.");
	}

	if (H5Tinsert(memory_datatype_id, "bucket_offset", HOFFSET(interval_index_entry_type, bucket_offset), H5T_NATIVE_HSIZE) < 0) {
		throw HVCFCreateException(__FILE__, __FUNCTION__, __LINE__, "Error while adding new member to compound datatype.");
	}

	if (H5Tinsert(memory_datatype_id, "bucket_size", HOFFSET(interval_index_entry_type, bucket_size), H5T_NATIVE_HSIZE) < 0) {
		throw HVCFCreateException(__FILE__, __FUNCTION__, __LINE__, "Error while adding new member to compound datatype.");
	}

	return memory_datatype_id.release();
}

hid_t HVCF::create_hash_index_entry_memory_datatype() throw (HVCFCreateException) {
	HDF5DatatypeIdentifier memory_datatype_id;

	if ((memory_datatype_id = H5Tcreate(H5T_COMPOUND, sizeof(hash_index_entry_type))) < 0) {
		throw HVCFCreateException(__FILE__, __FUNCTION__, __LINE__, "Error while creating compound datatype.");
	}

	if (H5Tinsert(memory_datatype_id, "bucket_offset", HOFFSET(hash_index_entry_type, bucket_offset), H5T_NATIVE_HSIZE) < 0) {
		throw HVCFCreateException(__FILE__, __FUNCTION__, __LINE__, "Error while adding new member to compound datatype.");
	}

	if (H5Tinsert(memory_datatype_id, "bucket_size", HOFFSET(hash_index_entry_type, bucket_size), H5T_NATIVE_HSIZE) < 0) {
		throw HVCFCreateException(__FILE__, __FUNCTION__, __LINE__, "Error while adding new member to compound datatype.");
	}

	return memory_datatype_id.release();
}

void HVCF::initialize_ull_index_buckets(hid_t group_id, const char* index_group_name) throw (HVCFWriteException) {
	HDF5GroupIdentifier index_group_id;
	HDF5DataspaceIdentifier dataspace_id;
	HDF5DatasetIdentifier dataset_id;
	HDF5PropertyIdentifier dataset_property_id;
	HDF5DatatypeIdentifier index_entry_type_id;

	hsize_t initial_dims[1]{0};
	hsize_t maximum_dims[1]{H5S_UNLIMITED};
	hsize_t chunk_dims[1]{1000};

	if ((index_entry_type_id = H5Topen(file_id, ULL_INDEX_ENTRY_TYPE, H5P_DEFAULT)) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while opening datatype.");
	}

	if ((index_group_id = H5Gopen(group_id, index_group_name, H5P_DEFAULT)) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while opening group.");
	}

	if ((dataspace_id = H5Screate_simple(1, initial_dims, maximum_dims)) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while creating dataspace.");
	}

	if ((dataset_property_id = H5Pcreate(H5P_DATASET_CREATE)) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while creating dataset property.");
	}

	if (strcmp(COMPRESSION, HVCFConfiguration::GZIP_COMPRESSION) == 0) {
		if ((H5Pset_chunk(dataset_property_id, 1, chunk_dims) < 0) || (H5Pset_deflate(dataset_property_id, COMPRESSION_LEVEL) < 0)) {
			throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while setting dataset properties.");
		}
	} else if (strcmp(COMPRESSION, HVCFConfiguration::BLOSC_LZ4HC_COMPRESSION) == 0) {
		unsigned int cd_values[7];
		cd_values[4] = COMPRESSION_LEVEL;
		// 0 -- shuffle not active, 1 -- shuffle active
		cd_values[5] = 1;
		// Compressor to use
		cd_values[6] = BLOSC_LZ4HC; // does better but slower compression. decompression is still very fast.
		if ((H5Pset_chunk(dataset_property_id, 1, chunk_dims) < 0) || (H5Pset_filter(dataset_property_id, FILTER_BLOSC, H5Z_FLAG_OPTIONAL, 7, cd_values) < 0)) {
			throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while setting dataset properties.");
		}
	} else {
		if (H5Pset_chunk(dataset_property_id, 1, chunk_dims) < 0) {
			throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while setting dataset properties.");
		}
	}

	if ((dataset_id = H5Dcreate(index_group_id, INDEX_BUCKETS, index_entry_type_id, dataspace_id, H5P_DEFAULT, dataset_property_id, H5P_DEFAULT)) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while creating dataset.");
	}
}

void HVCF::initialize_string_index_buckets(hid_t group_id, const char* index_group_name) throw (HVCFWriteException) {
	HDF5GroupIdentifier index_group_id;
	HDF5DataspaceIdentifier dataspace_id;
	HDF5DatasetIdentifier dataset_id;
	HDF5PropertyIdentifier dataset_property_id;
	HDF5DatatypeIdentifier index_entry_type_id;

	hsize_t initial_dims[1]{0};
	hsize_t maximum_dims[1]{H5S_UNLIMITED};
	hsize_t chunk_dims[1]{1000};

	if ((index_entry_type_id = H5Topen(file_id, STRING_INDEX_ENTRY_TYPE, H5P_DEFAULT)) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while opening datatype.");
	}

	if ((index_group_id = H5Gopen(group_id, index_group_name, H5P_DEFAULT)) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while opening group.");
	}

	if ((dataspace_id = H5Screate_simple(1, initial_dims, maximum_dims)) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while creating dataspace.");
	}

	if ((dataset_property_id = H5Pcreate(H5P_DATASET_CREATE)) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while creating dataset property.");
	}

	if (strcmp(COMPRESSION, HVCFConfiguration::GZIP_COMPRESSION) == 0) {
		if ((H5Pset_chunk(dataset_property_id, 1, chunk_dims) < 0) || (H5Pset_deflate(dataset_property_id, COMPRESSION_LEVEL) < 0)) {
			throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while setting dataset properties.");
		}
	} else if (strcmp(COMPRESSION, HVCFConfiguration::BLOSC_LZ4HC_COMPRESSION) == 0) {
		unsigned int cd_values[7];
		cd_values[4] = COMPRESSION_LEVEL;
		// 0 -- shuffle not active, 1 -- shuffle active
		cd_values[5] = 1;
		// Compressor to use
		cd_values[6] = BLOSC_LZ4HC; // does better but slower compression. decompression is still very fast.
		if ((H5Pset_chunk(dataset_property_id, 1, chunk_dims) < 0) || (H5Pset_filter(dataset_property_id, FILTER_BLOSC, H5Z_FLAG_OPTIONAL, 7, cd_values) < 0)) {
			throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while setting dataset properties.");
		}
	} else {
		if (H5Pset_chunk(dataset_property_id, 1, chunk_dims) < 0) {
			throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while setting dataset properties.");
		}
	}

	if ((dataset_id = H5Dcreate(index_group_id, INDEX_BUCKETS, index_entry_type_id, dataspace_id, H5P_DEFAULT, dataset_property_id, H5P_DEFAULT)) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while creating dataset.");
	}
}

void HVCF::cache_names_index_bucket(hid_t group_id, hash_index_entry_type& hash_index_entry, vector<string_index_entry_type>& bucket, vector<string_index_entry_type>& buckets_cache) throw (HVCFWriteException) {
	HDF5GroupIdentifier index_group_id;
	HDF5DatasetIdentifier dataset_id;
	HDF5DataspaceIdentifier file_dataspace_id;

	hsize_t file_dims[1]{0};
	hsize_t mem_dims[1]{bucket.size()};

	// BEGIN: sort bucket.
	sort(bucket.begin(), bucket.end(),
			[] (const string_index_entry_type& f, const string_index_entry_type& s) -> bool {
				return (strcmp(f.string_value, s.string_value) < 0);
	});
	// END: sort bucket.

	// BEGIN: get file offset.
	if ((index_group_id = H5Gopen(group_id, NAMES_INDEX_GROUP, H5P_DEFAULT)) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while opening group.");
	}

	if ((dataset_id = H5Dopen(index_group_id, INDEX_BUCKETS, H5P_DEFAULT)) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while opening dataset.");
	}

	if ((file_dataspace_id = H5Dget_space(dataset_id)) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while getting dataspace.");
	}

	if (H5Sget_simple_extent_dims(file_dataspace_id, file_dims, nullptr) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while getting dataspace dimensions.");
	}
	// END: get file offset.

	// BEGIN: set bucket location and write to buffer
	hash_index_entry.bucket_offset = file_dims[0] + buckets_cache.size();
	hash_index_entry.bucket_size = mem_dims[0];

	buckets_cache.insert(buckets_cache.end(), bucket.begin(), bucket.end());
	// END: set bucket location and write to buffer
}

void HVCF::write_names_index_buckets(hid_t group_id, vector<string_index_entry_type>& buckets_cache) throw (HVCFWriteException) {
	if (buckets_cache.size() == 0u) {
		return;
	}

	HDF5GroupIdentifier index_group_id;
	HDF5DatasetIdentifier dataset_id;
	HDF5DataspaceIdentifier file_dataspace_id;
	HDF5DataspaceIdentifier memory_dataspace_id;

	hsize_t file_dims[1]{0};
	hsize_t mem_dims[1]{buckets_cache.size()};
	hsize_t offset[1]{0};

	if ((memory_dataspace_id = H5Screate_simple(1, mem_dims, nullptr)) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while creating memory dataspace.");
	}

	// BEGIN: get file offset.
	if ((index_group_id = H5Gopen(group_id, NAMES_INDEX_GROUP, H5P_DEFAULT)) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while opening group.");
	}

	if ((dataset_id = H5Dopen(index_group_id, INDEX_BUCKETS, H5P_DEFAULT)) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while opening dataset.");
	}

	if ((file_dataspace_id = H5Dget_space(dataset_id)) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while getting dataspace.");
	}

	if (H5Sget_simple_extent_dims(file_dataspace_id, file_dims, nullptr) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while getting dataspace dimensions.");
	}

	file_dataspace_id.close();
	// END: get file offset.

	offset[0] = file_dims[0];
	file_dims[0] += mem_dims[0];

	// BEGIN: set new file extent and write
	if (H5Dset_extent(dataset_id, file_dims) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while setting dataset dimensions.");
	}

	if ((file_dataspace_id = H5Dget_space(dataset_id)) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while getting dataspace.");
	}

	if (H5Sselect_hyperslab(file_dataspace_id, H5S_SELECT_SET, offset, nullptr, mem_dims, nullptr) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while making selection in dataspace.");
	}

	if (H5Dwrite(dataset_id, string_index_entry_memory_datatype_id, memory_dataspace_id, file_dataspace_id, H5P_DEFAULT, buckets_cache.data()) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while writing to dataset.");
	}
	// END: set new file extent and write

	// BEGIN: deallocate char* allocated by HDF5
	if (H5Dvlen_reclaim(string_index_entry_memory_datatype_id, memory_dataspace_id, H5P_DEFAULT, buckets_cache.data()) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while reclaiming HDF5 memory.");
	}
	// END: deallocate char* allocated by HDF5

	buckets_cache.clear();
}

void HVCF::write_names_index(hid_t group_id, const hash_index_entry_type* hash_index_entries, unsigned int n_hash_index_entries) throw (HVCFWriteException) {
	if (n_hash_index_entries == 0u) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Empty list of offsets");
	}

	HDF5DatatypeIdentifier hash_index_entry_type_id;

	HDF5GroupIdentifier index_group_id;
	HDF5DatasetIdentifier dataset_id;
	HDF5DataspaceIdentifier file_dataspace_id;
	HDF5DataspaceIdentifier memory_dataspace_id;
	HDF5PropertyIdentifier dataset_property_id;

	hsize_t dims[1]{n_hash_index_entries};
	hsize_t chunk_dims[1]{100};

	if (dims[0] < chunk_dims[0]) {
		chunk_dims[0] = dims[0];
	}

	// BEGIN: open comitted datatype.
	if ((hash_index_entry_type_id = H5Topen(file_id, HASH_INDEX_ENTRY_TYPE, H5P_DEFAULT)) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while opening datatype.");
	}
	// END: open comitted datatype.

	if ((index_group_id = H5Gopen(group_id, NAMES_INDEX_GROUP, H5P_DEFAULT)) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while opening group.");
	}

	if ((file_dataspace_id = H5Screate_simple(1, dims, NULL)) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while creating dataspace.");
	}

	if ((dataset_property_id = H5Pcreate(H5P_DATASET_CREATE)) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while creating dataset property.");
	}

	if (strcmp(COMPRESSION, HVCFConfiguration::GZIP_COMPRESSION) == 0) {
		if ((H5Pset_chunk(dataset_property_id, 1, chunk_dims) < 0) || (H5Pset_deflate(dataset_property_id, COMPRESSION_LEVEL) < 0)) {
			throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while setting dataset properties.");
		}
	} else if (strcmp(COMPRESSION, HVCFConfiguration::BLOSC_LZ4HC_COMPRESSION) == 0) {
		unsigned int cd_values[7];
		cd_values[4] = COMPRESSION_LEVEL;
		// 0 -- shuffle not active, 1 -- shuffle active
		cd_values[5] = 1;
		// Compressor to use
		cd_values[6] = BLOSC_LZ4HC; // does better but slower compression. decompression is still very fast.
		if ((H5Pset_chunk(dataset_property_id, 1, chunk_dims) < 0) || (H5Pset_filter(dataset_property_id, FILTER_BLOSC, H5Z_FLAG_OPTIONAL, 7, cd_values) < 0)) {
			throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while setting dataset properties.");
		}
	} else {
		if (H5Pset_chunk(dataset_property_id, 1, chunk_dims) < 0) {
			throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while setting dataset properties.");
		}
	}

	if ((dataset_id = H5Dcreate(index_group_id, HASH_INDEX, hash_index_entry_type_id, file_dataspace_id, H5P_DEFAULT, dataset_property_id, H5P_DEFAULT)) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while creating dataset.");
	}

	if (H5Dwrite(dataset_id, hash_index_entry_memory_datatype_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, hash_index_entries) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while writing to dataset.");
	}
}

void HVCF::cache_intervals_index_bucket(hid_t group_id, interval_index_entry_type& interval_index_entry, vector<ull_index_entry_type>& bucket, vector<ull_index_entry_type>& buckets_cache) throw (HVCFWriteException) {
	HDF5GroupIdentifier index_group_id;
	HDF5DatasetIdentifier dataset_id;
	HDF5DataspaceIdentifier file_dataspace_id;

	hsize_t file_dims[1]{0};
	hsize_t mem_dims[1]{bucket.size()};

	// BEGIN: sort bucket's entries.
	sort(bucket.begin(), bucket.end(),
			[] (const ull_index_entry_type& f, const ull_index_entry_type& s) -> bool {
				return (f.ull_value < s.ull_value);
	});
	// END: sort bucket's entries.

	// BEGIN: get file offset.
	if ((index_group_id = H5Gopen(group_id, INTERVALS_INDEX_GROUP, H5P_DEFAULT)) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while opening group.");
	}

	if ((dataset_id = H5Dopen(index_group_id, INDEX_BUCKETS, H5P_DEFAULT)) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while opening dataset.");
	}

	if ((file_dataspace_id = H5Dget_space(dataset_id)) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while getting dataspace.");
	}

	if (H5Sget_simple_extent_dims(file_dataspace_id, file_dims, nullptr) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while getting dataspace dimensions.");
	}
	// END: get file offset.

	// BEGIN: set bucket location and write to buffer.
	interval_index_entry.ull_value_1 = bucket[0].ull_value;
	interval_index_entry.ull_value_2 = bucket[bucket.size() - 1].ull_value;
	interval_index_entry.offset_1 = bucket[0].offset;
	interval_index_entry.offset_2 = bucket[bucket.size() - 1].offset;
	interval_index_entry.bucket_offset = file_dims[0] + buckets_cache.size();
	interval_index_entry.bucket_size = mem_dims[0];

	buckets_cache.insert(buckets_cache.end(), bucket.begin(), bucket.end());
	// END: set bucket location and write to buffer.
}

void HVCF::write_intervals_index_buckets(hid_t group_id, vector<ull_index_entry_type>& buckets_cache) throw (HVCFWriteException) {
	if (buckets_cache.size() == 0u) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Empty list of offsets");
	}

	HDF5GroupIdentifier index_group_id;
	HDF5DatasetIdentifier dataset_id;
	HDF5DataspaceIdentifier file_dataspace_id;
	HDF5DataspaceIdentifier memory_dataspace_id;

	hsize_t file_dims[1]{0};
	hsize_t mem_dims[1]{buckets_cache.size()};
	hsize_t offset[1]{0};

	if ((memory_dataspace_id = H5Screate_simple(1, mem_dims, nullptr)) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while creating memory dataspace.");
	}

	// BEGIN: get file offset.
	if ((index_group_id = H5Gopen(group_id, INTERVALS_INDEX_GROUP, H5P_DEFAULT)) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while opening group.");
	}

	if ((dataset_id = H5Dopen(index_group_id, INDEX_BUCKETS, H5P_DEFAULT)) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while opening dataset.");
	}

	if ((file_dataspace_id = H5Dget_space(dataset_id)) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while getting dataspace.");
	}

	if (H5Sget_simple_extent_dims(file_dataspace_id, file_dims, nullptr) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while getting dataspace dimensions.");
	}
	file_dataspace_id.close();
	// END: get file offset.

	offset[0] = file_dims[0];
	file_dims[0] += mem_dims[0];

	// BEGIN: set new file extent and write
	if (H5Dset_extent(dataset_id, file_dims) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while setting dataset dimensions.");
	}

	if ((file_dataspace_id = H5Dget_space(dataset_id)) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while getting dataspace.");
	}

	if (H5Sselect_hyperslab(file_dataspace_id, H5S_SELECT_SET, offset, nullptr, mem_dims, nullptr) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while making selection in dataspace.");
	}

	if (H5Dwrite(dataset_id, ull_index_entry_memory_datatype_id, memory_dataspace_id, file_dataspace_id, H5P_DEFAULT, buckets_cache.data()) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while writing to dataset.");
	}
	// END: set new file extent and write

	buckets_cache.clear();
}

void HVCF::write_intervals_index(hid_t chromosome_group_id, const interval_index_entry_type* interval_index_entries, unsigned int n_interval_index_entries) throw (HVCFWriteException) {
	if (n_interval_index_entries == 0u) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Empty list of offsets");
	}

	HDF5DatatypeIdentifier interval_index_entry_type_id;
	HDF5DatatypeIdentifier interval_index_entry_memory_datatype_id;

	HDF5GroupIdentifier index_group_id;
	HDF5DatasetIdentifier dataset_id;
	HDF5DataspaceIdentifier file_dataspace_id;
	HDF5DataspaceIdentifier memory_dataspace_id;
	HDF5PropertyIdentifier dataset_property_id;

	hsize_t dims[1]{n_interval_index_entries};
	hsize_t chunk_dims[1]{100};

	if (dims[0] < chunk_dims[0]) {
		chunk_dims[0] = dims[0];
	}

	// BEGIN: open comitted datatype.
	if ((interval_index_entry_type_id = H5Topen(file_id, INTERVAL_INDEX_ENTRY_TYPE, H5P_DEFAULT)) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while opening datatype.");
	}
	// END: open comitted datatype.

	// BEGIN: create memory datatype for index.
	try {
		interval_index_entry_memory_datatype_id = create_interval_index_entry_memory_datatype();
	} catch (HVCFCreateException &e) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while creating memory datatypes.");
	}
	// END: create memory datatype index.

	if ((index_group_id = H5Gopen(chromosome_group_id, INTERVALS_INDEX_GROUP, H5P_DEFAULT)) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while opening group.");
	}

	if ((file_dataspace_id = H5Screate_simple(1, dims, NULL)) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while creating dataspace.");
	}

	if ((dataset_property_id = H5Pcreate(H5P_DATASET_CREATE)) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while creating dataset property.");
	}

	if (strcmp(COMPRESSION, HVCFConfiguration::GZIP_COMPRESSION) == 0) {
		if ((H5Pset_chunk(dataset_property_id, 1, chunk_dims) < 0) || (H5Pset_deflate(dataset_property_id, COMPRESSION_LEVEL) < 0)) {
			throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while setting dataset properties.");
		}
	} else if (strcmp(COMPRESSION, HVCFConfiguration::BLOSC_LZ4HC_COMPRESSION) == 0) {
		unsigned int cd_values[7];
		cd_values[4] = COMPRESSION_LEVEL;
		// 0 -- shuffle not active, 1 -- shuffle active
		cd_values[5] = 1;
		// Compressor to use
		cd_values[6] = BLOSC_LZ4HC; // does better but slower compression. decompression is still very fast.
		if ((H5Pset_chunk(dataset_property_id, 1, chunk_dims) < 0) || (H5Pset_filter(dataset_property_id, FILTER_BLOSC, H5Z_FLAG_OPTIONAL, 7, cd_values) < 0)) {
			throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while setting dataset properties.");
		}
	} else {
		if (H5Pset_chunk(dataset_property_id, 1, chunk_dims) < 0) {
			throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while setting dataset properties.");
		}
	}

	if ((dataset_id = H5Dcreate(index_group_id, INTERVALS_INDEX, interval_index_entry_type_id, file_dataspace_id, H5P_DEFAULT, dataset_property_id, H5P_DEFAULT)) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while creating dataset.");
	}

	if (H5Dwrite(dataset_id, interval_index_entry_memory_datatype_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, interval_index_entries) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while writing to dataset.");
	}
}

void HVCF::read_variant_names_into_bucket(hid_t chromosome_group_id, const vector<hsize_t>& offsets, vector<string_index_entry_type>& bucket) throw (HVCFWriteException) {
	if (offsets.size() == 0u) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Empty list of offsets");
	}

	HDF5DatasetIdentifier dataset_id;
	HDF5DataspaceIdentifier file_dataspace_id;
	HDF5DataspaceIdentifier memory_dataspace_id;

	hsize_t mem_dims[1]{offsets.size()};

	variants_entry_type variants_buffer[offsets.size()];

	// BEGIN: read variants.
	if ((dataset_id = H5Dopen(chromosome_group_id, VARIANTS_DATASET, H5P_DEFAULT)) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while opening dataset.");
	}

	if ((file_dataspace_id = H5Dget_space(dataset_id)) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while getting dataspace.");
	}

	if ((memory_dataspace_id = H5Screate_simple(1, mem_dims, nullptr)) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while creating memory dataspace.");
	}

	if (H5Sselect_elements(file_dataspace_id, H5S_SELECT_SET, offsets.size(), offsets.data()) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while making selection in dataspace.");
	}

	if (H5Dread(dataset_id, variants_entry_memory_datatype_id, memory_dataspace_id, file_dataspace_id, H5P_DEFAULT, variants_buffer) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while reading from dataset");
	}
	// END: read variants.

	bucket.resize(offsets.size());
	for (unsigned int i = 0u; i < offsets.size(); ++i) {
		bucket[i].string_value = variants_buffer[i].name;
		bucket[i].offset = offsets[i];
	}
}

void HVCF::read_sample_names_into_bucket(hid_t samples_group_id, const vector<hsize_t>& offsets, vector<string_index_entry_type>& bucket) throw (HVCFWriteException) {
	if (offsets.size() == 0u) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Empty list of offsets");
	}

	HDF5DatasetIdentifier dataset_id;
	HDF5DataspaceIdentifier file_dataspace_id;
	HDF5DataspaceIdentifier memory_dataspace_id;

	hsize_t mem_dims[1]{offsets.size()};

	char* sample_names_buffer[offsets.size()];

	// BEGIN: read samples.
	if ((dataset_id = H5Dopen(samples_group_id, SAMPLE_NAMES_DATASET, H5P_DEFAULT)) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while opening dataset.");
	}

	if ((file_dataspace_id = H5Dget_space(dataset_id)) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while getting dataspace.");
	}

	if ((memory_dataspace_id = H5Screate_simple(1, mem_dims, nullptr)) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while creating memory dataspace.");
	}

	if (H5Sselect_elements(file_dataspace_id, H5S_SELECT_SET, offsets.size(), offsets.data()) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while making selection in dataspace.");
	}

	if (H5Dread(dataset_id, native_string_datatype_id, memory_dataspace_id, file_dataspace_id, H5P_DEFAULT, sample_names_buffer) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while reading from dataset");
	}
	// END: read samples.

	bucket.resize(offsets.size());
	for (unsigned int i = 0u; i < offsets.size(); ++i) {
		bucket[i].string_value = sample_names_buffer[i];
		bucket[i].offset = offsets[i];
	}
}

void HVCF::read_positions_into_bucket(hid_t chromosome_group_id, const vector<hsize_t>& offsets, vector<ull_index_entry_type>& bucket) throw (HVCFWriteException) {
	if (offsets.size() == 0u) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Empty list of offsets");
	}

	HDF5DatasetIdentifier dataset_id;
	HDF5DataspaceIdentifier file_dataspace_id;
	HDF5DataspaceIdentifier memory_dataspace_id;

	hsize_t mem_dims[1]{offsets.size()};

	variants_entry_type variants_buffer[offsets.size()];

	// BEGIN: read variants.
	if ((dataset_id = H5Dopen(chromosome_group_id, VARIANTS_DATASET, H5P_DEFAULT)) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while opening dataset.");
	}

	if ((file_dataspace_id = H5Dget_space(dataset_id)) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while getting dataspace.");
	}

	if ((memory_dataspace_id = H5Screate_simple(1, mem_dims, nullptr)) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while creating memory dataspace.");
	}

	if (H5Sselect_elements(file_dataspace_id, H5S_SELECT_SET, offsets.size(), offsets.data()) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while making selection in dataspace.");
	}

	if (H5Dread(dataset_id, variants_entry_memory_datatype_id, memory_dataspace_id, file_dataspace_id, H5P_DEFAULT, variants_buffer) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while reading from dataset");
	}
	// END: read variants.

	bucket.resize(offsets.size());
	for (unsigned int i = 0u; i < offsets.size(); ++i) {
		bucket[i].ull_value = variants_buffer[i].position;
		bucket[i].offset = offsets[i];
	}
}

void HVCF::write_haplotypes(hid_t group_id, const unsigned char* buffer, unsigned int n_variants, unsigned int n_haplotypes) throw (HVCFWriteException) {
	HDF5DatasetIdentifier dataset_id;
	HDF5DataspaceIdentifier file_dataspace_id;
	HDF5DataspaceIdentifier memory_dataspace_id;

	hsize_t mem_dims[2]{n_variants, n_haplotypes};
	hsize_t file_dims[2]{0, 0};
	hsize_t file_offset[2]{0, 0};

	if ((dataset_id = H5Dopen(group_id, HAPLOTYPES_DATASET, H5P_DEFAULT)) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while opening dataset.");
	}

	if ((file_dataspace_id = H5Dget_space(dataset_id)) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while getting dataspace.");
	}

	if (H5Sget_simple_extent_dims(file_dataspace_id, file_dims, nullptr) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while getting dataspace dimensions.");
	}

	file_dataspace_id.close();

	file_offset[0] = file_dims[0];
	file_offset[1] = 0;
	file_dims[0] += mem_dims[0];

	if (H5Dset_extent(dataset_id, file_dims) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while setting dataset dimensions.");
	}

	if ((file_dataspace_id = H5Dget_space(dataset_id)) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while getting dataspace.");
	}

	if ((memory_dataspace_id = H5Screate_simple(2, mem_dims, nullptr)) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while creating memory dataspace.");
	}

	if (H5Sselect_hyperslab(file_dataspace_id, H5S_SELECT_SET, file_offset, nullptr, mem_dims, nullptr) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while making selection in dataspace.");
	}

	if (H5Dwrite(dataset_id, H5T_NATIVE_UCHAR, memory_dataspace_id, file_dataspace_id, H5P_DEFAULT, buffer) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while writing to dataset.");
	}
}

void HVCF::write_variants(hid_t group_id, const variants_entry_type* buffer, unsigned int n_variants) throw (HVCFWriteException) {
	HDF5DatasetIdentifier dataset_id;
	HDF5DataspaceIdentifier file_dataspace_id;
	HDF5DatatypeIdentifier variants_entry_memory_datatype_id;
	HDF5DataspaceIdentifier memory_dataspace_id;

	hsize_t mem_dims[1]{n_variants};
	hsize_t file_dims[1]{0};
	hsize_t file_offset[1]{0};

	// BEGIN: get current dimensions.
	if ((dataset_id = H5Dopen(group_id, VARIANTS_DATASET, H5P_DEFAULT)) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while opening dataset.");
	}

	if ((file_dataspace_id = H5Dget_space(dataset_id)) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while getting dataspace.");
	}

	if (H5Sget_simple_extent_dims(file_dataspace_id, file_dims, nullptr) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while getting dataspace dimensions.");
	}

	file_dataspace_id.close();
	// END: get current dimensions.

	// BEGIN: set new dimensions.
	file_offset[0] = file_dims[0];
	file_dims[0] += mem_dims[0];

	if (H5Dset_extent(dataset_id, file_dims) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while setting dataset dimensions.");
	}
	// END: set new dimensions.

	// BEGIN: create memory datatype.
	try {
		variants_entry_memory_datatype_id = create_variants_entry_memory_datatype();
	} catch (HVCFCreateException &e) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while creating memory datatypes.");
	}
	// END: create memory datatype.

	// BEGIN: write.
	if ((file_dataspace_id = H5Dget_space(dataset_id)) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while getting dataspace.");
	}

	if ((memory_dataspace_id = H5Screate_simple(1, mem_dims, nullptr)) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while creating memory dataspace.");
	}

	if (H5Sselect_hyperslab(file_dataspace_id, H5S_SELECT_SET, file_offset, nullptr, mem_dims, nullptr) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while making selection in dataspace.");
	}

	if (H5Dwrite(dataset_id, variants_entry_memory_datatype_id, memory_dataspace_id, file_dataspace_id, H5P_DEFAULT, buffer) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while writing to dataset.");
	}
	// END: write.
}

hid_t HVCF::create_sample_names_dataset(hid_t group_id, hsize_t chunk_size) throw (HVCFWriteException) {
	HDF5DataspaceIdentifier dataspace_id;
	HDF5DatasetIdentifier dataset_id;
	HDF5PropertyIdentifier dataset_property_id;
	HDF5DatatypeIdentifier native_string_datatype_id;

	hsize_t initial_dims[1]{0};
	hsize_t maximum_dims[1]{H5S_UNLIMITED};
	hsize_t chunk_dims[1]{chunk_size};

	if ((native_string_datatype_id = H5Topen(file_id, VARIABLE_LENGTH_STRING_TYPE, H5P_DEFAULT)) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while opening datatype.");
	}

	if ((dataspace_id = H5Screate_simple(1, initial_dims, maximum_dims)) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while creating dataspace.");
	}

	if ((dataset_property_id = H5Pcreate(H5P_DATASET_CREATE)) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while creating dataset property.");
	}

	if (strcmp(COMPRESSION, HVCFConfiguration::GZIP_COMPRESSION) == 0) {
		if ((H5Pset_chunk(dataset_property_id, 1, chunk_dims) < 0) || (H5Pset_deflate(dataset_property_id, COMPRESSION_LEVEL) < 0)) {
			throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while setting dataset properties.");
		}
	} else if (strcmp(COMPRESSION, HVCFConfiguration::BLOSC_LZ4HC_COMPRESSION) == 0) {
		unsigned int cd_values[7];
		cd_values[4] = COMPRESSION_LEVEL;
		// 0 -- shuffle not active, 1 -- shuffle active
		cd_values[5] = 1;
		// Compressor to use
		cd_values[6] = BLOSC_LZ4HC; // does better but slower compression. decompression is still very fast.
		if ((H5Pset_chunk(dataset_property_id, 1, chunk_dims) < 0) || (H5Pset_filter(dataset_property_id, FILTER_BLOSC, H5Z_FLAG_OPTIONAL, 7, cd_values) < 0)) {
			throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while setting dataset properties.");
		}
	} else {
		if (H5Pset_chunk(dataset_property_id, 1, chunk_dims) < 0) {
			throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while setting dataset properties.");
		}
	}

	if ((dataset_id = H5Dcreate(group_id, SAMPLE_NAMES_DATASET, native_string_datatype_id, dataspace_id, H5P_DEFAULT, dataset_property_id, H5P_DEFAULT)) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while creating dataset.");
	}

	return dataset_id.release();
}

hid_t HVCF::create_sample_subsets_dataset(hid_t group_id, hsize_t chunk_size) throw (HVCFWriteException) {
	HDF5DataspaceIdentifier dataspace_id;
	HDF5DatasetIdentifier dataset_id;
	HDF5PropertyIdentifier dataset_property_id;
	HDF5DatatypeIdentifier subsets_entry_type_id;

	hsize_t initial_dims[1]{0};
	hsize_t maximum_dims[1]{H5S_UNLIMITED};
	hsize_t chunk_dims[1]{chunk_size};

	if ((subsets_entry_type_id = H5Topen(file_id, SUBSETS_ENTRY_TYPE, H5P_DEFAULT)) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while opening datatype.");
	}

	if ((dataspace_id = H5Screate_simple(1, initial_dims, maximum_dims)) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while creating dataspace.");
	}

	if ((dataset_property_id = H5Pcreate(H5P_DATASET_CREATE)) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while creating dataset property.");
	}

	if (strcmp(COMPRESSION, HVCFConfiguration::GZIP_COMPRESSION) == 0) {
		if ((H5Pset_chunk(dataset_property_id, 1, chunk_dims) < 0) || (H5Pset_deflate(dataset_property_id, COMPRESSION_LEVEL) < 0)) {
			throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while setting dataset properties.");
		}
	} else if (strcmp(COMPRESSION, HVCFConfiguration::BLOSC_LZ4HC_COMPRESSION) == 0) {
		unsigned int cd_values[7];
		cd_values[4] = COMPRESSION_LEVEL;
		// 0 -- shuffle not active, 1 -- shuffle active
		cd_values[5] = 1;
		// Compressor to use
		cd_values[6] = BLOSC_LZ4HC; // does better but slower compression. decompression is still very fast.
		if ((H5Pset_chunk(dataset_property_id, 1, chunk_dims) < 0) || (H5Pset_filter(dataset_property_id, FILTER_BLOSC, H5Z_FLAG_OPTIONAL, 7, cd_values) < 0)) {
			throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while setting dataset properties.");
		}
	} else {
		if (H5Pset_chunk(dataset_property_id, 1, chunk_dims) < 0) {
			throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while setting dataset properties.");
		}
	}

	if ((dataset_id = H5Dcreate(group_id, SAMPLE_SUBSETS_DATASET, subsets_entry_type_id, dataspace_id, H5P_DEFAULT, dataset_property_id, H5P_DEFAULT)) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while creating dataset.");
	}

	return dataset_id.release();
}

hid_t HVCF::create_haplotypes_dataset(hid_t group_id, hsize_t variants_chunk_size, hsize_t samples_chunk_size) throw (HVCFWriteException) {
	HDF5DataspaceIdentifier dataspace_id;
	HDF5DatasetIdentifier dataset_id;
	HDF5PropertyIdentifier dataset_property_id;

	hsize_t n_samples = get_n_samples();
	hsize_t n_haplotypes = 2 * n_samples;
	hsize_t haplotypes_chunk_size = 2 * samples_chunk_size;

	if (haplotypes_chunk_size > n_haplotypes) {
		haplotypes_chunk_size = n_haplotypes;
	}

	hsize_t initial_dims[2]{0, n_haplotypes};
	hsize_t maximum_dims[2]{H5S_UNLIMITED, n_haplotypes};
	hsize_t chunk_dims[2]{variants_chunk_size, haplotypes_chunk_size};

	if ((dataspace_id = H5Screate_simple(2, initial_dims, maximum_dims)) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while creating dataspace.");
	}

	if ((dataset_property_id = H5Pcreate(H5P_DATASET_CREATE)) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while creating dataset property.");
	}

	if (strcmp(COMPRESSION, HVCFConfiguration::GZIP_COMPRESSION) == 0) {
		if ((H5Pset_chunk(dataset_property_id, 2, chunk_dims) < 0) || (H5Pset_deflate(dataset_property_id, COMPRESSION_LEVEL) < 0)) {
			throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while setting dataset properties.");
		}
	} else if (strcmp(COMPRESSION, HVCFConfiguration::BLOSC_LZ4HC_COMPRESSION) == 0) {
		unsigned int cd_values[7];
		cd_values[4] = COMPRESSION_LEVEL;
		// 0 -- shuffle not active, 1 -- shuffle active
		cd_values[5] = 1;
		// Compressor to use
		cd_values[6] = BLOSC_LZ4HC; // does better but slower compression. decompression is still very fast.
		if ((H5Pset_chunk(dataset_property_id, 2, chunk_dims) < 0) || (H5Pset_filter(dataset_property_id, FILTER_BLOSC, H5Z_FLAG_OPTIONAL, 7, cd_values) < 0)) {
			throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while setting dataset properties.");
		}
	} else {
		if (H5Pset_chunk(dataset_property_id, 2, chunk_dims) < 0) {
			throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while setting dataset properties.");
		}
	}

	if ((dataset_id = H5Dcreate(group_id, HAPLOTYPES_DATASET, H5T_NATIVE_UCHAR, dataspace_id, H5P_DEFAULT, dataset_property_id, H5P_DEFAULT)) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while creating dataset.");
	}

	return dataset_id.release();
}

hid_t HVCF::create_variants_dataset(hid_t group_id, hsize_t chunk_size) throw (HVCFWriteException) {
	HDF5DatatypeIdentifier datatype_id;
	HDF5DataspaceIdentifier dataspace_id;
	HDF5DatasetIdentifier dataset_id;
	HDF5PropertyIdentifier dataset_property_id;

	hsize_t initial_dims[1]{0};
	hsize_t maximum_dims[1]{H5S_UNLIMITED};
	hsize_t chunk_dims[1]{chunk_size};

	if ((datatype_id = H5Topen(file_id, VARIANTS_ENTRY_TYPE, H5P_DEFAULT)) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while opening datatype.");
	}

	if ((dataspace_id = H5Screate_simple(1, initial_dims, maximum_dims)) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while creating dataspace.");
	}

	if ((dataset_property_id = H5Pcreate(H5P_DATASET_CREATE)) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while creating dataset property.");
	}

	if (strcmp(COMPRESSION, HVCFConfiguration::GZIP_COMPRESSION) == 0) {
		if ((H5Pset_chunk(dataset_property_id, 1, chunk_dims) < 0) || (H5Pset_deflate(dataset_property_id, COMPRESSION_LEVEL) < 0)) {
			throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while setting dataset properties.");
		}
	} else if (strcmp(COMPRESSION, HVCFConfiguration::BLOSC_LZ4HC_COMPRESSION) == 0) {
		unsigned int cd_values[7];
		cd_values[4] = COMPRESSION_LEVEL;
		// 0 -- shuffle not active, 1 -- shuffle active
		cd_values[5] = 1;
		// Compressor to use
		cd_values[6] = BLOSC_LZ4HC; // does better but slower compression. decompression is still very fast.
		if ((H5Pset_chunk(dataset_property_id, 1, chunk_dims) < 0) || (H5Pset_filter(dataset_property_id, FILTER_BLOSC, H5Z_FLAG_OPTIONAL, 7, cd_values) < 0)) {
			throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while setting dataset properties.");
		}
	} else {
		if (H5Pset_chunk(dataset_property_id, 1, chunk_dims) < 0) {
			throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while setting dataset properties.");
		}
	}

	if ((dataset_id = H5Dcreate(group_id, VARIANTS_DATASET, datatype_id, dataspace_id, H5P_DEFAULT, dataset_property_id, H5P_DEFAULT)) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while creating dataset.");
	}

	return dataset_id.release();
}

void HVCF::create_chromosome_indices(hid_t chromosome_group_id) throw (HVCFWriteException) {
	if ((H5Lexists(chromosome_group_id, NAMES_INDEX_GROUP, H5P_DEFAULT) > 0) ||
			(H5Lexists(chromosome_group_id, INTERVALS_INDEX_GROUP, H5P_DEFAULT) > 0)) {
		return;
	}

	HDF5GroupIdentifier names_index_group_id;
	HDF5GroupIdentifier intervals_index_group_id;
	HDF5DatasetIdentifier dataset_id;
	HDF5DataspaceIdentifier file_dataspace_id;
	HDF5DatatypeIdentifier memory_datatype_id;
	HDF5DataspaceIdentifier memory_dataspace_id;

	const unsigned int read_chunk_size = 100000;

	hsize_t file_dims[1]{0};
	hsize_t file_offset[1]{0};
	hsize_t mem_dims[1]{read_chunk_size};

	variants_entry_type buffer[read_chunk_size];

	unordered_map<unsigned int, vector<hsize_t>> names_index_buckets;
	auto names_index_buckets_it = names_index_buckets.end();

	// index construction for positions is based on assumption that variants dataset are already ordered by position.
	// we group positions into intervals. every interval will have not more than MAX_VARIANTS_PER_INTERVAL positions.
	vector<vector<hsize_t>> intervals;

	std::hash<string> hash_function;
	size_t hash_value;

	// BEGIN: create groups for indices.
	if ((names_index_group_id = H5Gcreate(chromosome_group_id, NAMES_INDEX_GROUP, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT)) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while creating group.");
	}

	if ((intervals_index_group_id = H5Gcreate(chromosome_group_id, INTERVALS_INDEX_GROUP, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT)) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while creating group.");
	}
	// END: create groups for indices.

	// BEGIN: get dataset dimensions.
	if ((dataset_id = H5Dopen(chromosome_group_id, VARIANTS_DATASET, H5P_DEFAULT)) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while opening dataset.");
	}

	if ((file_dataspace_id = H5Dget_space(dataset_id)) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while getting dataspace.");
	}

	if (H5Sget_simple_extent_dims(file_dataspace_id, file_dims, nullptr) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while getting dataspace dimensions.");
	}
	// END: get dataset dimensions.

	// BEGIN: create memory datatype.
	try {
		memory_datatype_id = create_variants_entry_memory_datatype();
	} catch (HVCFCreateException &e) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while creating memory datatypes.");
	}
	// END: create memory datatype.

	// BEGIN: read data in chunks.
	while (file_offset[0] < file_dims[0]) {
		if (file_offset[0] + read_chunk_size > file_dims[0]) {
			mem_dims[0] = file_dims[0] - file_offset[0];
		} else {
			mem_dims[0] = read_chunk_size;
		}

		if ((memory_dataspace_id = H5Screate_simple(1, mem_dims, nullptr)) < 0) {
			throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while creating memory dataspace.");
		}

		if (H5Sselect_hyperslab(file_dataspace_id, H5S_SELECT_SET, file_offset, NULL, mem_dims, NULL) < 0) {
			throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while making selection in dataspace.");
		}

		if (H5Dread(dataset_id, memory_datatype_id, memory_dataspace_id, file_dataspace_id, H5P_DEFAULT, buffer) < 0) {
			throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while reading from dataset");
		}

		for (unsigned int i = 0; i < mem_dims[0]; ++i) {
			hash_value = hash_function(buffer[i].name) % N_VARIANTS_HASH_BUCKETS;
			names_index_buckets_it = names_index_buckets.find(hash_value);
			if (names_index_buckets_it == names_index_buckets.end()) {
				names_index_buckets_it = names_index_buckets.emplace(hash_value, vector<hsize_t>()).first;
			}
			names_index_buckets_it->second.push_back(file_offset[0]);

			if ((intervals.size() == 0) || (intervals.back().size() >= MAX_VARIANTS_IN_INTERVAL_BUCKET)) {
				intervals.emplace_back();
			}
			intervals.back().push_back(file_offset[0]);

			file_offset[0] += 1;
		}

		if (H5Dvlen_reclaim(memory_datatype_id, memory_dataspace_id, H5P_DEFAULT, buffer) < 0) {
			throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while reclaiming HDF5 memory.");
		}

		memory_dataspace_id.close();
	}
	// END: read data in chunks.

	// BEGIN: write indices on disk.
	initialize_string_index_buckets(chromosome_group_id, NAMES_INDEX_GROUP);

	hash_index_entry_type hash_index_keys[N_VARIANTS_HASH_BUCKETS];
	vector<string_index_entry_type> hash_bucket;
	vector<string_index_entry_type> hash_buckets_cache;

	for (unsigned int i = 0u; i < N_VARIANTS_HASH_BUCKETS; ++i) {
		hash_index_keys[i].bucket_offset = 0;
		hash_index_keys[i].bucket_size = 0;
	}
	hash_bucket.reserve(1000);
	hash_buckets_cache.reserve(100001);

	for (auto&& bucket_offsets : names_index_buckets) {
		read_variant_names_into_bucket(chromosome_group_id, bucket_offsets.second, hash_bucket);
		cache_names_index_bucket(chromosome_group_id, hash_index_keys[bucket_offsets.first], hash_bucket, hash_buckets_cache);
		if (hash_buckets_cache.size() > 100000) {
			write_names_index_buckets(chromosome_group_id, hash_buckets_cache);
		}
	}
	write_names_index_buckets(chromosome_group_id, hash_buckets_cache);
	write_names_index(chromosome_group_id, hash_index_keys, N_VARIANTS_HASH_BUCKETS);

	initialize_ull_index_buckets(chromosome_group_id, INTERVALS_INDEX_GROUP);

	interval_index_entry_type interval_index_keys[intervals.size()];
	vector<ull_index_entry_type> intervals_bucket;
	vector<ull_index_entry_type> intervals_buckets_cache;
	intervals_bucket.reserve(1000);
	intervals_buckets_cache.reserve(100001);

	for (unsigned int i = 0; i < intervals.size(); ++i) {
		read_positions_into_bucket(chromosome_group_id, intervals[i], intervals_bucket);
		cache_intervals_index_bucket(chromosome_group_id, interval_index_keys[i], intervals_bucket, intervals_buckets_cache);
		if (intervals_buckets_cache.size() > 100000) {
			write_intervals_index_buckets(chromosome_group_id, intervals_buckets_cache);
		}
	}
	write_intervals_index_buckets(chromosome_group_id, intervals_buckets_cache);
	write_intervals_index(chromosome_group_id, interval_index_keys, intervals.size());
	// END: write indices on disk.
}

void HVCF::create_samples_indices() throw (HVCFWriteException) {
	if (H5Lexists(samples_group_id, NAMES_INDEX_GROUP, H5P_DEFAULT) > 0) {
		return;
	}

	HDF5GroupIdentifier names_index_group_id;
	HDF5DatasetIdentifier samples_all_dataset_id;
	HDF5DataspaceIdentifier file_dataspace_id;
	HDF5DatatypeIdentifier native_string_datatype_id;
	hsize_t file_dims[1]{0};

	if ((names_index_group_id = H5Gcreate(samples_group_id, NAMES_INDEX_GROUP, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT)) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while creating group.");
	}

	if ((native_string_datatype_id = H5Topen(file_id, VARIABLE_LENGTH_STRING_TYPE, H5P_DEFAULT)) < 0) {
		throw HVCFReadException(__FILE__, __FUNCTION__, __LINE__, "Error while opening datatype.");
	}

	if ((samples_all_dataset_id = H5Dopen(samples_group_id, SAMPLE_NAMES_DATASET, H5P_DEFAULT)) < 0) {
		throw HVCFReadException(__FILE__, __FUNCTION__, __LINE__, "Error while opening dataset.");
	}

	if ((file_dataspace_id = H5Dget_space(samples_all_dataset_id)) < 0) {
		throw HVCFReadException(__FILE__, __FUNCTION__, __LINE__, "Error while getting dataspace.");
	}

	if (H5Sget_simple_extent_dims(file_dataspace_id, file_dims, nullptr) < 0) {
		throw HVCFReadException(__FILE__, __FUNCTION__, __LINE__, "Error while getting dataspace dimensions.");
	}

	char* buffer[file_dims[0]];

	unordered_map<unsigned int, vector<hsize_t>> names_index_buckets;
	auto names_index_buckets_it = names_index_buckets.end();

	std::hash<string> hash_function;
	size_t hash_value;

	if (H5Dread(samples_all_dataset_id, native_string_datatype_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, buffer) < 0) {
		throw HVCFReadException(__FILE__, __FUNCTION__, __LINE__, "Error while reading from dataset");
	}

	for (unsigned int i = 0; i < file_dims[0]; ++i) {
		hash_value = hash_function(buffer[i]) % N_SAMPLES_HASH_BUCKETS;
		names_index_buckets_it = names_index_buckets.find(hash_value);
		if (names_index_buckets_it == names_index_buckets.end()) {
			names_index_buckets_it = names_index_buckets.emplace(hash_value, vector<hsize_t>()).first;
		}
		names_index_buckets_it->second.push_back(i);
	}

	if (H5Dvlen_reclaim(native_string_datatype_id, file_dataspace_id, H5P_DEFAULT, buffer) < 0) {
		throw HVCFReadException(__FILE__, __FUNCTION__, __LINE__, "Error while reclaiming HDF5 memory.");
	}

	// BEGIN: write indices on disk.
	initialize_string_index_buckets(samples_group_id, NAMES_INDEX_GROUP);

	hash_index_entry_type hash_index_keys[N_SAMPLES_HASH_BUCKETS];

	vector<string_index_entry_type> bucket;
	vector<string_index_entry_type> buckets_buffer;

	for (unsigned int i = 0u; i < N_SAMPLES_HASH_BUCKETS; ++i) {
		hash_index_keys[i].bucket_offset = 0;
		hash_index_keys[i].bucket_size = 0;
	}
	bucket.reserve(100);
	buckets_buffer.reserve(100);

	for (auto&& bucket_offsets : names_index_buckets) {
		read_sample_names_into_bucket(samples_group_id, bucket_offsets.second, bucket);
		cache_names_index_bucket(samples_group_id, hash_index_keys[bucket_offsets.first], bucket, buckets_buffer);
		if (buckets_buffer.size() > 100) {
			write_names_index_buckets(samples_group_id, buckets_buffer);
		}
	}
	write_names_index_buckets(samples_group_id, buckets_buffer);
	write_names_index(samples_group_id, hash_index_keys, N_SAMPLES_HASH_BUCKETS);
}

void HVCF::create_indices() throw (HVCFWriteException) {
	create_samples_indices();
	for (auto&& chromosome : chromosomes) {
		create_chromosome_indices(chromosome.second->get());
	}
}

void HVCF::write_samples(const vector<string>& samples) throw (HVCFWriteException) {
	if (H5Lexists(samples_group_id, SAMPLE_NAMES_DATASET, H5P_DEFAULT) > 0) {
		return;
	}

	HDF5DatasetIdentifier dataset_id;
	HDF5DataspaceIdentifier file_dataspace_id;
	HDF5DataspaceIdentifier memory_dataspace_id;
	HDF5DatatypeIdentifier native_string_datatype_id;
	HDF5DatatypeIdentifier subsets_entry_memory_datatype_id;

	hsize_t mem_dims[1]{samples.size()};
	hsize_t file_offset[1]{0};
	hsize_t file_dims[1]{samples.size()};
	const char* buffer[samples.size()];

	for (unsigned int i = 0u; i < samples.size(); ++i) {
		buffer[i] = samples.at(i).c_str();
	}

	if ((native_string_datatype_id = H5Topen(file_id, VARIABLE_LENGTH_STRING_TYPE, H5P_DEFAULT)) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while opening datatype.");
	}

	if ((dataset_id = create_sample_names_dataset(samples_group_id, 1000)) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while creating dataset.");
	}

	if (H5Dset_extent(dataset_id, file_dims) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while setting dataset dimensions.");
	}

	if ((file_dataspace_id = H5Dget_space(dataset_id)) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while getting dataspace.");
	}

	if ((memory_dataspace_id = H5Screate_simple(1, mem_dims, nullptr)) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while creating memory dataspace.");
	}

	if (H5Sselect_hyperslab(file_dataspace_id, H5S_SELECT_SET, file_offset, nullptr, mem_dims, nullptr) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while making selection in dataspace.");
	}

	if (H5Dwrite(dataset_id, native_string_datatype_id, memory_dataspace_id, file_dataspace_id, H5P_DEFAULT, buffer) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while writing to dataset.");
	}

	dataset_id.close();
	memory_dataspace_id.close();
	file_dataspace_id.close();

	mem_dims[0] = 1;
	file_offset[0] = 0;
	file_dims[0] = 1;
	subsets_entry_type subsets_entry;

	subsets_entry.name = const_cast<char*>("ALL");
	subsets_entry.offset_1 = 0;
	subsets_entry.offset_2 = samples.size() - 1;

	try {
		subsets_entry_memory_datatype_id  = create_subsets_entry_memory_datatype();
	} catch (HVCFCreateException &e) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while creating memory datatype.");
	}

	if ((dataset_id = create_sample_subsets_dataset(samples_group_id, 5000)) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while creating dataset.");
	}

	if (H5Dset_extent(dataset_id, file_dims) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while setting dataset dimensions.");
	}

	if ((file_dataspace_id = H5Dget_space(dataset_id)) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while getting dataspace.");
	}

	if ((memory_dataspace_id = H5Screate_simple(1, mem_dims, nullptr)) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while creating memory dataspace.");
	}

	if (H5Sselect_hyperslab(file_dataspace_id, H5S_SELECT_SET, file_offset, nullptr, mem_dims, nullptr) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while making selection in dataspace.");
	}

	if (H5Dwrite(dataset_id, subsets_entry_memory_datatype_id, memory_dataspace_id, file_dataspace_id, H5P_DEFAULT, &subsets_entry) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while writing to dataset.");
	}
}

void HVCF::flush_write_buffer(future<void>& async_write) throw (HVCFWriteException) {
	auto buffers_it = write_buffers.end();
	for (auto&& entry : chromosomes) {
		buffers_it = write_buffers.find(entry.first);
		if (!buffers_it->second->is_empty()) {
			if (async_write.valid()) {
				async_write.wait();
			}
			auto flushed = buffers_it->second->flush();
			write_haplotypes(entry.second->get(), std::get<0>(flushed), std::get<2>(flushed), std::get<3>(flushed));
			write_variants(entry.second->get(), std::get<1>(flushed), std::get<2>(flushed));
		}
	}
}

void HVCF::write_variant(const Variant& variant, future<void>& async_write) throw (HVCFWriteException) {
	const string& chromosome = variant.get_chrom().get_value();
	auto chromosomes_it = chromosomes.end();
	auto buffers_it = write_buffers.end();

	if (variant.get_alt().get_values().size() != 1) { // Support only bi-allelic (for computing LD it is fine, but must be extended).
		return;
	}

	if (chromosomes.count(chromosome) == 0) {
		chromosomes_it = chromosomes.emplace(chromosome, std::move(unique_ptr<HDF5GroupIdentifier>(new HDF5GroupIdentifier()))).first;
		buffers_it = write_buffers.emplace(chromosome, std::move(unique_ptr<WriteBuffer>(new WriteBuffer(100000, get_n_samples())))).first;
		chromosomes_it->second->set(create_chromosome_group(chromosome));
	} else {
		chromosomes_it = chromosomes.find(chromosome);
		buffers_it = write_buffers.find(chromosome);
	}

	if (buffers_it->second->is_full()) {
		if (async_write.valid()) {
			async_write.wait();
		}

		auto write = [&](
				hid_t group_id,
				const unsigned char* haplotypes,
				const variants_entry_type* variants,
				unsigned int n_variants,
				unsigned int n_haplotypes) -> void {
			write_haplotypes(group_id, haplotypes, n_variants, n_haplotypes);
			write_variants(group_id, variants, n_variants);
		};

		auto flushed = buffers_it->second->flush();

		async_write = async(std::launch::async,
				write, chromosomes_it->second->get(), std::get<0>(flushed), std::get<1>(flushed), std::get<2>(flushed), std::get<3>(flushed));
	}

	buffers_it->second->add_variant(variant);
}

hid_t HVCF::create_chromosome_group(const string& name) throw (HVCFWriteException) {
	HDF5GroupIdentifier group_id;
	HDF5DatasetIdentifier dataset_id;

	if (name.length() == 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Empty dataspace name.");
	}

	if ((group_id = H5Gcreate(chromosomes_group_id, name.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT)) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while creating group.");
	}

	dataset_id = create_haplotypes_dataset(group_id, VARIANTS_CHUNK_SIZE, SAMPLES_CHUNK_SIZE);
	dataset_id.close();

	dataset_id = create_variants_dataset(group_id, VARIANTS_CHUNK_SIZE);
	dataset_id.close();

	return group_id.release();
}

void HVCF::load_chromosomes_cache() throw (HVCFReadException) {
	HDF5GroupIdentifier index_group_id;
	auto chromosomes_cache_it = chromosomes_cache.end();

	chromosomes_cache.clear();

	for (auto&& chromosome : chromosomes) {
		chromosomes_cache_it = chromosomes_cache.emplace(chromosome.first, std::move(unique_ptr<chromosomes_cache_entry>(new chromosomes_cache_entry()))).first;

		if ((index_group_id = H5Gopen(chromosome.second->get(), NAMES_INDEX_GROUP, H5P_DEFAULT)) < 0) {
			throw HVCFReadException(__FILE__, __FUNCTION__, __LINE__, "Error while opening group.");
		}
		chromosomes_cache_it->second->names_index_id.set(H5Dopen(index_group_id, HASH_INDEX, H5P_DEFAULT));
		if (chromosomes_cache_it->second->names_index_id.get() < 0) {
			throw HVCFReadException(__FILE__, __FUNCTION__, __LINE__, "Error while opening dataset.");
		}
		chromosomes_cache_it->second->names_index_buckets_id.set(H5Dopen(index_group_id, INDEX_BUCKETS, H5P_DEFAULT));
		if (chromosomes_cache_it->second->names_index_buckets_id.get() < 0) {
			throw HVCFReadException(__FILE__, __FUNCTION__, __LINE__, "Error while opening dataset.");
		}
		index_group_id.close();

		if ((index_group_id = H5Gopen(chromosome.second->get(), INTERVALS_INDEX_GROUP, H5P_DEFAULT)) < 0) {
			throw HVCFReadException(__FILE__, __FUNCTION__, __LINE__, "Error while opening group.");
		}
		chromosomes_cache_it->second->intervals_index_id.set(H5Dopen(index_group_id, INTERVALS_INDEX, H5P_DEFAULT));
		if (chromosomes_cache_it->second->intervals_index_id.get() < 0) {
			throw HVCFReadException(__FILE__, __FUNCTION__, __LINE__, "Error while opening dataset.");
		}
		chromosomes_cache_it->second->intervals_index_buckets_id.set(H5Dopen(index_group_id, INDEX_BUCKETS, H5P_DEFAULT));
		if (chromosomes_cache_it->second->intervals_index_buckets_id.get() < 0) {
			throw HVCFReadException(__FILE__, __FUNCTION__, __LINE__, "Error while opening dataset.");
		}
		index_group_id.close();

		chromosomes_cache_it->second->variants_id.set(H5Dopen(chromosome.second->get(), VARIANTS_DATASET, H5P_DEFAULT));
		if (chromosomes_cache_it->second->variants_id.get() < 0) {
			throw HVCFReadException(__FILE__, __FUNCTION__, __LINE__, "Error while opening dataset.");
		}

		chromosomes_cache_it->second->haplotypes_id.set(H5Dopen(chromosome.second->get(), HAPLOTYPES_DATASET, H5P_DEFAULT));
		if (chromosomes_cache_it->second->haplotypes_id.get() < 0) {
			throw HVCFReadException(__FILE__, __FUNCTION__, __LINE__, "Error while opening dataset.");
		}
	}
}

void HVCF::load_samples_cache() throw (HVCFReadException) {
	HDF5GroupIdentifier index_group_id;
	HDF5DatasetIdentifier dataset_id;
	HDF5DataspaceIdentifier file_dataspace_id;
	HDF5DataspaceIdentifier memory_dataspace_id;
	HDF5DatatypeIdentifier subsets_entry_memory_datatype_id;

	hsize_t file_dims[1]{0};

	samples_cache.names_index_id.close();
	samples_cache.names_index_buckets_id.close();
	samples_cache.subsets.clear();

	if ((index_group_id = H5Gopen(samples_group_id, NAMES_INDEX_GROUP, H5P_DEFAULT)) < 0) {
		throw HVCFReadException(__FILE__, __FUNCTION__, __LINE__, "Error while opening group.");
	}

	if ((samples_cache.names_index_id = H5Dopen(index_group_id, HASH_INDEX, H5P_DEFAULT)) < 0) {
		throw HVCFReadException(__FILE__, __FUNCTION__, __LINE__, "Error while opening dataset.");
	}

	if ((samples_cache.names_index_buckets_id = H5Dopen(index_group_id, INDEX_BUCKETS, H5P_DEFAULT)) < 0) {
		throw HVCFReadException(__FILE__, __FUNCTION__, __LINE__, "Error while opening dataset.");
	}

	try {
		subsets_entry_memory_datatype_id  = create_subsets_entry_memory_datatype();
	} catch (HVCFCreateException &e) {
		throw HVCFReadException(__FILE__, __FUNCTION__, __LINE__, "Error while creating memory datatype.");
	}

	if ((dataset_id = H5Dopen(samples_group_id, SAMPLE_SUBSETS_DATASET, H5P_DEFAULT)) < 0) {
		throw HVCFReadException(__FILE__, __FUNCTION__, __LINE__, "Error while creating dataset.");
	}

	if ((file_dataspace_id = H5Dget_space(dataset_id)) < 0) {
		throw HVCFReadException(__FILE__, __FUNCTION__, __LINE__, "Error while getting dataspace.");
	}

	if (H5Sget_simple_extent_dims(file_dataspace_id, file_dims, nullptr) < 0) {
		throw HVCFReadException(__FILE__, __FUNCTION__, __LINE__, "Error while getting dataspace dimensions.");
	}

	subsets_entry_type subsets_entry_buffer[file_dims[0]];

	if (H5Dread(dataset_id, subsets_entry_memory_datatype_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, subsets_entry_buffer) < 0) {
		throw HVCFReadException(__FILE__, __FUNCTION__, __LINE__, "Error while reading from dataset");
	}

	auto subsets_cache_it = samples_cache.subsets.end();
	hsize_t offset_1 = 0;
	hsize_t offset_2 = 0;
	hsize_t size = 0;
	for (hsize_t i = 0; i < file_dims[0]; ++i) {
		subsets_cache_it = samples_cache.subsets.emplace(subsets_entry_buffer[i].name, subsets_cache_entry()).first;
		offset_1 = subsets_entry_buffer[i].offset_1;
		offset_2 = subsets_entry_buffer[i].offset_2;
		size = offset_2 - offset_1 + 1u;
		subsets_cache_it->second.chunks.emplace_back(offset_1, offset_2, size);
		subsets_cache_it->second.n_samples += size;
	}

	if (H5Dvlen_reclaim(subsets_entry_memory_datatype_id, file_dataspace_id, H5P_DEFAULT, subsets_entry_buffer) < 0) {
		throw HVCFReadException(__FILE__, __FUNCTION__, __LINE__, "Error while reclaiming HDF5 memory.");
	}
}

void HVCF::load_cache() throw (HVCFReadException) {
	load_samples_cache();
	load_chromosomes_cache();
}

void HVCF::create(const string& name) throw (HVCFWriteException) {
	HDF5DatatypeIdentifier datatype_id;
	HDF5PropertyIdentifier file_access_property_id;
	H5AC_cache_config_t config;
	size_t native_string_datatype_size = 0;

	this->name = name;

	if ((file_access_property_id = H5Pcreate(H5P_FILE_ACCESS)) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while creating file access property.");
	}

	config.version = H5AC__CURR_CACHE_CONFIG_VERSION;
	if (H5Pget_mdc_config(file_access_property_id, &config) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while getting metadata cache configuration.");
	}

	config.set_initial_size = true;
	config.initial_size = METADATA_CACHE_INITIAL_SIZE;
	config.min_size = METADATA_CACHE_MIN_SIZE;
	config.max_size = METADATA_CACHE_MAX_SIZE;

	if (H5Pset_mdc_config(file_access_property_id, &config) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while setting metadata cache configuration.");
	}

	if (H5Pset_sieve_buf_size(file_access_property_id, SIEVE_BUFFER_MAX_SIZE) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while setting sieve buffer size.");
	}

	if (H5Pset_cache(file_access_property_id, 0, CHUNK_CACHE_N_SLOTS, CHUNK_CACHE_SIZE, 0.75) < 0) {
		throw HVCFOpenException(__FILE__, __FUNCTION__, __LINE__, "Error while setting cache parameters.");
	}

	if ((file_id = H5Fcreate(this->name.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, file_access_property_id)) < 0) {
		throw HVCFOpenException(__FILE__, __FUNCTION__, __LINE__, "Error while creating file.");
	}

//	VERBOSE
//	HDF5PropertyIdentifier file_access_property_id_test;
//	H5AC_cache_config_t cache_config;
//	size_t sieve_buffer_size, cache_nslots, cache_nbytes;
//	double w0;
//	cache_config.version = H5AC__CURR_CACHE_CONFIG_VERSION;
//	file_access_property_id_test = H5Fget_access_plist(file_id);
//	H5Pget_mdc_config(file_access_property_id_test, &cache_config);
//	H5Pget_sieve_buf_size(file_access_property_id_test, &sieve_buffer_size);
//	H5Pget_cache(file_access_property_id_test, NULL, &cache_nslots, &cache_nbytes, &w0);
//	cout << "Created file configuration:" << endl;
//	cout << "\tInitital size: " << cache_config.initial_size << endl;
//	cout << "\tSieve buffer size: " << sieve_buffer_size << endl;
//	cout << "\tCache slots: " << cache_nslots << endl;
//	cout << "\tCache size: " << cache_nbytes << endl;

	// BEGIN: create and commit variable length string datatype.
	if (((datatype_id = H5Tcopy(H5T_C_S1)) < 0) || (H5Tset_size(datatype_id, H5T_VARIABLE) < 0)) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while creating variable length string datatype.");
	}

	if (H5Tcommit(file_id, VARIABLE_LENGTH_STRING_TYPE, datatype_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while committing variable length string datatype.");
	}
	datatype_id.close();
	// END: create and commit variable length string datatype.

	if ((native_string_datatype_id = H5Topen(file_id, VARIABLE_LENGTH_STRING_TYPE, H5P_DEFAULT)) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while opening datatype.");
	}

	if ((native_string_datatype_size = H5Tget_size(native_string_datatype_id)) == 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while obtaining datatype size.");
	}

	// BEGIN: create and commit compound datatype.
	if ((datatype_id = H5Tcreate(H5T_COMPOUND, sizeof(unsigned long long int) + sizeof(hsize_t))) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while creating compound datatype.");
	}

	if (H5Tinsert(datatype_id, "ull_value", 0, H5T_NATIVE_ULLONG) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while adding new member to compound datatype.");
	}

	if (H5Tinsert(datatype_id, "offset", sizeof(unsigned long long int), H5T_NATIVE_HSIZE) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while adding new member to compound datatype.");
	}

	if (H5Tcommit(file_id, ULL_INDEX_ENTRY_TYPE, datatype_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while committing compound datatype.");
	}
	datatype_id.close();
	// END: create and commit compound datatype.

	// BEGIN: create and commit compound datatype.
	if ((datatype_id = H5Tcreate(H5T_COMPOUND, native_string_datatype_size + sizeof(hsize_t))) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while creating compound datatype.");
	}

	if (H5Tinsert(datatype_id, "string_value", 0, native_string_datatype_id) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while adding new member to compound datatype.");
	}

	if (H5Tinsert(datatype_id, "offset", native_string_datatype_size, H5T_NATIVE_HSIZE) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while adding new member to compound datatype.");
	}

	if	(H5Tcommit(file_id, STRING_INDEX_ENTRY_TYPE, datatype_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while committing compound datatype.");
	}
	datatype_id.close();
	// END: create and commit compound datatype.

	// BEGIN: create and commit compound datatype.
	if ((datatype_id = H5Tcreate(H5T_COMPOUND, 2 * sizeof(unsigned long long int) + 4 * sizeof(hsize_t))) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while creating compound datatype.");
	}

	if (H5Tinsert(datatype_id, "ull_value_1", 0, H5T_NATIVE_ULLONG) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while adding new member to compound datatype.");
	}

	if (H5Tinsert(datatype_id, "ull_value_2", sizeof(unsigned long long int), H5T_NATIVE_ULLONG) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while adding new member to compound datatype.");
	}

	if (H5Tinsert(datatype_id, "offset_1", 2 * sizeof(unsigned long long int), H5T_NATIVE_HSIZE) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while adding new member to compound datatype.");
	}

	if (H5Tinsert(datatype_id, "offset_2", 2 * sizeof(unsigned long long int) + sizeof(hsize_t), H5T_NATIVE_HSIZE) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while adding new member to compound datatype.");
	}

	if (H5Tinsert(datatype_id, "bucket_offset", 2 * sizeof(unsigned long long int) + 2 * sizeof(hsize_t), H5T_NATIVE_HSIZE) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while adding new member to compound datatype.");
	}

	if (H5Tinsert(datatype_id, "bucket_size", 2 * sizeof(unsigned long long int) + 3 * sizeof(hsize_t), H5T_NATIVE_HSIZE) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while adding new member to compound datatype.");
	}

	if	(H5Tcommit(file_id, INTERVAL_INDEX_ENTRY_TYPE, datatype_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while committing compound datatype.");
	}
	datatype_id.close();
	// END: create and commit compound datatype.

	// BEGIN: create and commit compound datatype.
	if ((datatype_id = H5Tcreate(H5T_COMPOUND, 2 * sizeof(hsize_t))) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while creating compound datatype.");
	}

	if (H5Tinsert(datatype_id, "bucket_offset", 0, H5T_NATIVE_HSIZE) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while adding new member to compound datatype.");
	}

	if (H5Tinsert(datatype_id, "bucket_size", sizeof(hsize_t), H5T_NATIVE_HSIZE) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while adding new member to compound datatype.");
	}

	if	(H5Tcommit(file_id, HASH_INDEX_ENTRY_TYPE, datatype_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while committing compound datatype.");
	}
	datatype_id.close();
	// END: create and commit compound datatype.

	// BEGIN: create and commit compound datatype.
	if ((datatype_id = H5Tcreate(H5T_COMPOUND, 3u * native_string_datatype_size + sizeof(unsigned long long int))) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while creating compound datatype.");
	}

	if (H5Tinsert(datatype_id, "name", 0, native_string_datatype_id) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while adding new member to compound datatype.");
	}

	if (H5Tinsert(datatype_id, "ref", native_string_datatype_size, native_string_datatype_id) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while adding new member to compound datatype.");
	}

	if (H5Tinsert(datatype_id, "alt", 2u * native_string_datatype_size, native_string_datatype_id) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while adding new member to compound datatype.");
	}

	if (H5Tinsert(datatype_id, "position", 3u * native_string_datatype_size, H5T_NATIVE_ULLONG) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while adding new member to compound datatype.");
	}

	if	(H5Tcommit(file_id, VARIANTS_ENTRY_TYPE, datatype_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while committing compound datatype.");
	}
	datatype_id.close();
	// END: create and commit compound datatype.

	// BEGIN: create and commit compound datatype.
	if ((datatype_id = H5Tcreate(H5T_COMPOUND, native_string_datatype_size + 2 * sizeof(hsize_t))) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while creating compound datatype.");
	}

	if (H5Tinsert(datatype_id, "name", 0, native_string_datatype_id) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while adding new member to compound datatype.");
	}

	if (H5Tinsert(datatype_id, "offset_1", native_string_datatype_size, H5T_NATIVE_HSIZE) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while adding new member to compound datatype.");
	}

	if (H5Tinsert(datatype_id, "offset_2", native_string_datatype_size + sizeof(hsize_t), H5T_NATIVE_HSIZE) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while adding new member to compound datatype.");
	}

	if	(H5Tcommit(file_id, SUBSETS_ENTRY_TYPE, datatype_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while committing compound datatype.");
	}
	datatype_id.close();
	// END: create and commit compound datatype.

	if ((samples_group_id = H5Gcreate(file_id, SAMPLES_GROUP, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT)) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while creating group.");
	}

	if ((chromosomes_group_id = H5Gcreate(file_id, CHROMOSOMES_GROUP, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT)) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while creating group.");
	}

	try {
		hash_index_entry_memory_datatype_id = create_hash_index_entry_memory_datatype();
		string_index_entry_memory_datatype_id = create_string_index_entry_memory_datatype();
		ull_index_entry_memory_datatype_id = create_ull_index_entry_memory_datatype();
		interval_index_entry_memory_datatype_id = create_interval_index_entry_memory_datatype();
		variants_entry_memory_datatype_id = create_variants_entry_memory_datatype();
	} catch (HVCFCreateException &e) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while creating memory datatypes.");
	}
}

void HVCF::open(const string& name) throw (HVCFOpenException) {
	HDF5PropertyIdentifier file_access_property_id;
	H5AC_cache_config_t config;

	this->name = name;

	if ((file_access_property_id = H5Pcreate(H5P_FILE_ACCESS)) < 0) {
		throw HVCFOpenException(__FILE__, __FUNCTION__, __LINE__, "Error while creating file access property.");
	}

	config.version = H5AC__CURR_CACHE_CONFIG_VERSION;
	if (H5Pget_mdc_config(file_access_property_id, &config) < 0) {
		throw HVCFOpenException(__FILE__, __FUNCTION__, __LINE__, "Error while getting metadata cache configuration.");
	}

	config.set_initial_size = true;
	config.initial_size = METADATA_CACHE_INITIAL_SIZE;
	config.min_size = METADATA_CACHE_MIN_SIZE;
	config.max_size = METADATA_CACHE_MAX_SIZE;

	if (H5Pset_mdc_config(file_access_property_id, &config) < 0) {
		throw HVCFOpenException(__FILE__, __FUNCTION__, __LINE__, "Error while setting metadata cache configuration.");
	}

	if (H5Pset_sieve_buf_size(file_access_property_id, SIEVE_BUFFER_MAX_SIZE) < 0) {
		throw HVCFOpenException(__FILE__, __FUNCTION__, __LINE__, "Error while setting sieve buffer size.");
	}

	if (H5Pset_cache(file_access_property_id, 0, CHUNK_CACHE_N_SLOTS, CHUNK_CACHE_SIZE, 0.75) < 0) {
		throw HVCFOpenException(__FILE__, __FUNCTION__, __LINE__, "Error while setting cache parameters.");
	}

	if ((file_id = H5Fopen(this->name.c_str(), H5F_ACC_RDONLY, file_access_property_id)) < 0) {
		throw HVCFOpenException(__FILE__, __FUNCTION__, __LINE__, "Error while opening file.");
	}

//	VERBOSE
//	HDF5PropertyIdentifier file_access_property_id_test;
//	H5AC_cache_config_t cache_config;
//	size_t sieve_buffer_size, cache_nslots, cache_nbytes;
//	double w0;
//	cache_config.version = H5AC__CURR_CACHE_CONFIG_VERSION;
//	file_access_property_id_test = H5Fget_access_plist(file_id);
//	H5Pget_mdc_config(file_access_property_id_test, &cache_config);
//	H5Pget_sieve_buf_size(file_access_property_id_test, &sieve_buffer_size);
//	H5Pget_cache(file_access_property_id_test, NULL, &cache_nslots, &cache_nbytes, &w0);
//	cout << "Created file configuration:" << endl;
//	cout << "\tInitital size: " << cache_config.initial_size << endl;
//	cout << "\tSieve buffer size: " << sieve_buffer_size << endl;
//	cout << "\tCache slots: " << cache_nslots << endl;
//	cout << "\tCache size: " << cache_nbytes << endl;

	if ((samples_group_id = H5Gopen(file_id, SAMPLES_GROUP, H5P_DEFAULT)) < 0) {
		throw HVCFOpenException(__FILE__, __FUNCTION__, __LINE__, "Error while opening group.");
	}

	if ((chromosomes_group_id = H5Gopen(file_id, CHROMOSOMES_GROUP, H5P_DEFAULT)) < 0) {
		throw HVCFOpenException(__FILE__, __FUNCTION__, __LINE__, "Error while opening group.");
	}

	H5G_info_t chromosomes_group_info;
	H5O_info_t object_info;
	hsize_t object_name_length = 0;
	unique_ptr<char[]> object_name = nullptr;
	auto chromosomes_it = chromosomes.end();

	if (H5Gget_info(chromosomes_group_id, &chromosomes_group_info) < 0) {
		throw HVCFOpenException(__FILE__, __FUNCTION__, __LINE__, "Error while getting group information.");
	}

	for (hsize_t i = 0; i < chromosomes_group_info.nlinks; ++i) {
		if (H5Oget_info_by_idx(chromosomes_group_id, ".", H5_INDEX_NAME, H5_ITER_INC, i, &object_info, 0) < 0) {
			throw HVCFOpenException(__FILE__, __FUNCTION__, __LINE__, "Error while getting object information.");
		}

		if (object_info.type != H5O_type_t::H5O_TYPE_GROUP) {
			continue;
		}

		object_name_length = H5Lget_name_by_idx(chromosomes_group_id, ".", H5_INDEX_NAME, H5_ITER_INC, i, NULL, 0, H5P_DEFAULT);
		if (object_name_length < 0) {
			throw HVCFOpenException(__FILE__, __FUNCTION__, __LINE__, "Error while getting group name size.");
		}

		object_name = unique_ptr<char[]>(new char[object_name_length + 1]{});

		if (H5Lget_name_by_idx(chromosomes_group_id, ".", H5_INDEX_NAME, H5_ITER_INC, i, object_name.get(), object_name_length + 1, H5P_DEFAULT) < 0) {
			throw HVCFOpenException(__FILE__, __FUNCTION__, __LINE__, "Error while getting group name.");
		}

		chromosomes_it = chromosomes.emplace(object_name.get(), std::move(unique_ptr<HDF5GroupIdentifier>(new HDF5GroupIdentifier()))).first;
		chromosomes_it->second->set(H5Gopen(chromosomes_group_id, chromosomes_it->first.c_str(), H5P_DEFAULT));
		if (chromosomes_it->second->get() < 0) {
			throw HVCFOpenException(__FILE__, __FUNCTION__, __LINE__, "Error while opening group.");
		}
	}

	if ((native_string_datatype_id = H5Topen(file_id, VARIABLE_LENGTH_STRING_TYPE, H5P_DEFAULT)) < 0) {
		throw HVCFOpenException(__FILE__, __FUNCTION__, __LINE__, "Error while opening datatype.");
	}

	try {
		hash_index_entry_memory_datatype_id = create_hash_index_entry_memory_datatype();
		string_index_entry_memory_datatype_id = create_string_index_entry_memory_datatype();
		ull_index_entry_memory_datatype_id = create_ull_index_entry_memory_datatype();
		interval_index_entry_memory_datatype_id = create_interval_index_entry_memory_datatype();
		variants_entry_memory_datatype_id = create_variants_entry_memory_datatype();
	} catch (HVCFCreateException &e) {
		throw HVCFOpenException(__FILE__, __FUNCTION__, __LINE__, "Error while creating memory datatypes.");
	}

	try {
		load_cache();
	} catch (HVCFReadException &e) {
		throw HVCFOpenException(__FILE__, __FUNCTION__, __LINE__, "Error while creating cache.");
	}
}

void HVCF::close() throw (HVCFCloseException) {
	samples_cache.subsets.clear();
	samples_cache.names_index_id.close();
	samples_cache.names_index_buckets_id.close();
	chromosomes_cache.clear();
	chromosomes.clear();
	chromosomes_group_id.close();
	samples_group_id.close();
	native_string_datatype_id.close();
	hash_index_entry_memory_datatype_id.close();
	string_index_entry_memory_datatype_id.close();
	interval_index_entry_memory_datatype_id.close();
	ull_index_entry_memory_datatype_id.close();
	variants_entry_memory_datatype_id.close();
	file_id.close();
	name.clear();
}

void HVCF::import_vcf(const string& name) throw (HVCFWriteException) {
	VCFReader vcf;

	try {
		future<void> async_write;

		vcf.open(name);
		write_samples(std::move(vcf.get_variant().get_samples()));
		while (vcf.read_next_variant()) {
			write_variant(vcf.get_variant(), async_write);
		}
		flush_write_buffer(async_write);

		vcf.close();
	} catch (ReaderException &e) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while reading VCF file.");
	} catch (VCFException &e) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while reading VCF file.");
	}

	create_indices();
	load_cache();
}

void HVCF::create_sample_subset(const string& name, const std::vector<string>& samples) throw (HVCFWriteException) {
	HDF5DatasetIdentifier dataset_id;
	HDF5DataspaceIdentifier file_dataspace_id;
	HDF5DataspaceIdentifier memory_dataspace_id;
	HDF5DatatypeIdentifier subsets_entry_memory_datatype_id;

	map<string, unsigned int> all_samples;
	auto all_samples_it = all_samples.end();

	hsize_t offests_buffer[samples.size()];
	vector<tuple<hsize_t, hsize_t>> squeezed_offsets;

	for (auto&& sample : get_samples()) {
		all_samples.emplace(std::move(sample), all_samples.size());
	}

	for (unsigned int i = 0u; i < samples.size(); ++i) {
		all_samples_it = all_samples.find(samples[i]);
		if (all_samples_it == all_samples.end()) {
			throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Sample not found.");
		}
		offests_buffer[i] = static_cast<hsize_t>(all_samples_it->second);
	}

	std::sort(offests_buffer, offests_buffer + samples.size(), std::less_equal<hsize_t>());

	hsize_t start = offests_buffer[0];
	hsize_t end = offests_buffer[0];
	for (unsigned int i = 1; i < samples.size(); ++i) {
		if (offests_buffer[i] - offests_buffer[i - 1] > 1) {
			squeezed_offsets.emplace_back(start, end);
			start = offests_buffer[i];
			end = offests_buffer[i];
		} else {
			end = offests_buffer[i];
		}
	}
	squeezed_offsets.emplace_back(start, end);

	hsize_t mem_dims[1]{squeezed_offsets.size()};
	hsize_t file_dims[1]{0};
	hsize_t file_offset[1]{0};


	subsets_entry_type populations_entry_buffer[squeezed_offsets.size()];
	for (unsigned int i = 0; i < squeezed_offsets.size(); ++i) {
		populations_entry_buffer[i].name = const_cast<char*>(name.c_str());
		populations_entry_buffer[i].offset_1 = get<0>(squeezed_offsets[i]);
		populations_entry_buffer[i].offset_2 = get<1>(squeezed_offsets[i]);
	}

	try {
		subsets_entry_memory_datatype_id  = create_subsets_entry_memory_datatype();
	} catch (HVCFCreateException &e) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while creating memory datatype.");
	}

	if ((dataset_id = H5Dopen(samples_group_id, SAMPLE_SUBSETS_DATASET, H5P_DEFAULT)) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while creating dataset.");
	}

	if ((file_dataspace_id = H5Dget_space(dataset_id)) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while getting dataspace.");
	}

	if (H5Sget_simple_extent_dims(file_dataspace_id, file_dims, nullptr) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while getting dataspace dimensions.");
	}

	file_dataspace_id.close();

	file_offset[0] = file_dims[0];
	file_dims[0] += mem_dims[0];

	if (H5Dset_extent(dataset_id, file_dims) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while setting dataset dimensions.");
	}

	if ((file_dataspace_id = H5Dget_space(dataset_id)) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while getting dataspace.");
	}

	if ((memory_dataspace_id = H5Screate_simple(1, mem_dims, nullptr)) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while creating memory dataspace.");
	}

	if (H5Sselect_hyperslab(file_dataspace_id, H5S_SELECT_SET, file_offset, nullptr, mem_dims, nullptr) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while making selection in dataspace.");
	}

	if (H5Dwrite(dataset_id, subsets_entry_memory_datatype_id, memory_dataspace_id, file_dataspace_id, H5P_DEFAULT, populations_entry_buffer) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while writing to dataset.");
	}

	try {
		load_samples_cache();
	} catch (HVCFReadException &e) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while creating subset cache.");
	}
}

hsize_t HVCF::get_n_samples() throw (HVCFReadException) {
	HDF5DatasetIdentifier samples_all_dataset_id;
	HDF5DataspaceIdentifier file_dataspace_id;
	hsize_t file_dims[1]{0};

	if ((samples_all_dataset_id = H5Dopen(samples_group_id, SAMPLE_NAMES_DATASET, H5P_DEFAULT)) < 0) {
		throw HVCFReadException(__FILE__, __FUNCTION__, __LINE__, "Error while opening dataset.");
	}

	if ((file_dataspace_id = H5Dget_space(samples_all_dataset_id)) < 0) {
		throw HVCFReadException(__FILE__, __FUNCTION__, __LINE__, "Error while getting dataspace.");
	}

	if (H5Sget_simple_extent_dims(file_dataspace_id, file_dims, nullptr) < 0) {
		throw HVCFReadException(__FILE__, __FUNCTION__, __LINE__, "Error while getting dataspace dimensions.");
	}

	return file_dims[0];
}

vector<string> HVCF::get_samples() throw (HVCFReadException) {
	HDF5DatasetIdentifier samples_all_dataset_id;
	HDF5DataspaceIdentifier file_dataspace_id;
	HDF5DatatypeIdentifier native_string_datatype_id;
	hsize_t file_dims[1]{0};

	if ((native_string_datatype_id = H5Topen(file_id, VARIABLE_LENGTH_STRING_TYPE, H5P_DEFAULT)) < 0) {
		throw HVCFReadException(__FILE__, __FUNCTION__, __LINE__, "Error while opening datatype.");
	}

	if ((samples_all_dataset_id = H5Dopen(samples_group_id, SAMPLE_NAMES_DATASET, H5P_DEFAULT)) < 0) {
		throw HVCFReadException(__FILE__, __FUNCTION__, __LINE__, "Error while opening dataset.");
	}

	if ((file_dataspace_id = H5Dget_space(samples_all_dataset_id)) < 0) {
		throw HVCFReadException(__FILE__, __FUNCTION__, __LINE__, "Error while getting dataspace.");
	}

	if (H5Sget_simple_extent_dims(file_dataspace_id, file_dims, nullptr) < 0) {
		throw HVCFReadException(__FILE__, __FUNCTION__, __LINE__, "Error while getting dataspace dimensions.");
	}

	char* buffer[file_dims[0]];
	vector<string> samples;

	if (H5Dread(samples_all_dataset_id, native_string_datatype_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, buffer) < 0) {
		throw HVCFReadException(__FILE__, __FUNCTION__, __LINE__, "Error while reading from dataset");
	}

	for (unsigned int i = 0; i < file_dims[0]; ++i) {
		samples.emplace_back(buffer[i]);
	}

	if (H5Dvlen_reclaim(native_string_datatype_id, file_dataspace_id, H5P_DEFAULT, buffer) < 0) {
		throw HVCFReadException(__FILE__, __FUNCTION__, __LINE__, "Error while reclaiming HDF5 memory.");
	}

	return samples;
}

unsigned int HVCF::get_n_sample_subsets() throw (HVCFReadException) {
	HDF5DatasetIdentifier dataset_id;
	HDF5DataspaceIdentifier file_dataspace_id;
	HDF5DataspaceIdentifier memory_dataspace_id;
	HDF5DatatypeIdentifier populations_entry_memory_datatype_id;

	hsize_t file_dims[1]{0};

	try {
		populations_entry_memory_datatype_id  = create_subsets_entry_memory_datatype();
	} catch (HVCFCreateException &e) {
		throw HVCFReadException(__FILE__, __FUNCTION__, __LINE__, "Error while creating memory datatype.");
	}

	if ((dataset_id = H5Dopen(samples_group_id, SAMPLE_SUBSETS_DATASET, H5P_DEFAULT)) < 0) {
		throw HVCFReadException(__FILE__, __FUNCTION__, __LINE__, "Error while creating dataset.");
	}

	if ((file_dataspace_id = H5Dget_space(dataset_id)) < 0) {
		throw HVCFReadException(__FILE__, __FUNCTION__, __LINE__, "Error while getting dataspace.");
	}

	if (H5Sget_simple_extent_dims(file_dataspace_id, file_dims, nullptr) < 0) {
		throw HVCFReadException(__FILE__, __FUNCTION__, __LINE__, "Error while getting dataspace dimensions.");
	}

	subsets_entry_type populations_entry_buffer[file_dims[0]];

	if (H5Dread(dataset_id, populations_entry_memory_datatype_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, populations_entry_buffer) < 0) {
		throw HVCFReadException(__FILE__, __FUNCTION__, __LINE__, "Error while reading from dataset");
	}

	set<string> populations;

	for (hsize_t i = 0; i < file_dims[0]; ++i) {
		populations.emplace(populations_entry_buffer[i].name);
	}

	if (H5Dvlen_reclaim(populations_entry_memory_datatype_id, file_dataspace_id, H5P_DEFAULT, populations_entry_buffer) < 0) {
		throw HVCFReadException(__FILE__, __FUNCTION__, __LINE__, "Error while reclaiming HDF5 memory.");
	}

	return populations.size();
}

vector<string> HVCF::get_sample_subsets() throw (HVCFReadException) {
	HDF5DatasetIdentifier dataset_id;
	HDF5DataspaceIdentifier file_dataspace_id;
	HDF5DataspaceIdentifier memory_dataspace_id;
	HDF5DatatypeIdentifier subsets_entry_memory_datatype_id;

	hsize_t file_dims[1]{0};

	try {
		subsets_entry_memory_datatype_id  = create_subsets_entry_memory_datatype();
	} catch (HVCFCreateException &e) {
		throw HVCFReadException(__FILE__, __FUNCTION__, __LINE__, "Error while creating memory datatype.");
	}

	if ((dataset_id = H5Dopen(samples_group_id, SAMPLE_SUBSETS_DATASET, H5P_DEFAULT)) < 0) {
		throw HVCFReadException(__FILE__, __FUNCTION__, __LINE__, "Error while creating dataset.");
	}

	if ((file_dataspace_id = H5Dget_space(dataset_id)) < 0) {
		throw HVCFReadException(__FILE__, __FUNCTION__, __LINE__, "Error while getting dataspace.");
	}

	if (H5Sget_simple_extent_dims(file_dataspace_id, file_dims, nullptr) < 0) {
		throw HVCFReadException(__FILE__, __FUNCTION__, __LINE__, "Error while getting dataspace dimensions.");
	}

	subsets_entry_type subsets_entry_buffer[file_dims[0]];

	if (H5Dread(dataset_id, subsets_entry_memory_datatype_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, subsets_entry_buffer) < 0) {
		throw HVCFReadException(__FILE__, __FUNCTION__, __LINE__, "Error while reading from dataset");
	}

	set<string> subsets;

	for (hsize_t i = 0; i < file_dims[0]; ++i) {
		subsets.emplace(subsets_entry_buffer[i].name);
	}

	if (H5Dvlen_reclaim(subsets_entry_memory_datatype_id, file_dataspace_id, H5P_DEFAULT, subsets_entry_buffer) < 0) {
		throw HVCFReadException(__FILE__, __FUNCTION__, __LINE__, "Error while reclaiming HDF5 memory.");
	}

	vector<string> result(subsets.begin(), subsets.end());

	return result;
}

unsigned int HVCF::get_n_samples_in_subset(const string& name) throw (HVCFReadException) {
	HDF5DatasetIdentifier dataset_id;
	HDF5DataspaceIdentifier file_dataspace_id;
	HDF5DataspaceIdentifier memory_dataspace_id;
	HDF5DatatypeIdentifier subsets_entry_memory_datatype_id;

	hsize_t file_dims[1]{0};

	try {
		subsets_entry_memory_datatype_id  = create_subsets_entry_memory_datatype();
	} catch (HVCFCreateException &e) {
		throw HVCFReadException(__FILE__, __FUNCTION__, __LINE__, "Error while creating memory datatype.");
	}

	if ((dataset_id = H5Dopen(samples_group_id, SAMPLE_SUBSETS_DATASET, H5P_DEFAULT)) < 0) {
		throw HVCFReadException(__FILE__, __FUNCTION__, __LINE__, "Error while creating dataset.");
	}

	if ((file_dataspace_id = H5Dget_space(dataset_id)) < 0) {
		throw HVCFReadException(__FILE__, __FUNCTION__, __LINE__, "Error while getting dataspace.");
	}

	if (H5Sget_simple_extent_dims(file_dataspace_id, file_dims, nullptr) < 0) {
		throw HVCFReadException(__FILE__, __FUNCTION__, __LINE__, "Error while getting dataspace dimensions.");
	}

	subsets_entry_type subsets_entry_buffer[file_dims[0]];

	if (H5Dread(dataset_id, subsets_entry_memory_datatype_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, subsets_entry_buffer) < 0) {
		throw HVCFReadException(__FILE__, __FUNCTION__, __LINE__, "Error while reading from dataset");
	}

	unsigned int n_samples = 0u;

	for (hsize_t i = 0; i < file_dims[0]; ++i) {
		if (strcmp(subsets_entry_buffer[i].name, name.c_str()) == 0) {
			n_samples += subsets_entry_buffer[i].offset_2 - subsets_entry_buffer[i].offset_1 + 1u;
		}
	}

	if (H5Dvlen_reclaim(subsets_entry_memory_datatype_id, file_dataspace_id, H5P_DEFAULT, subsets_entry_buffer) < 0) {
		throw HVCFReadException(__FILE__, __FUNCTION__, __LINE__, "Error while reclaiming HDF5 memory.");
	}

	return n_samples;
}

vector<string> HVCF::get_samples_in_subset(const string& name) throw (HVCFReadException) {
	HDF5DatasetIdentifier dataset_id;
	HDF5DataspaceIdentifier file_dataspace_id;
	HDF5DataspaceIdentifier memory_dataspace_id;
	HDF5DatatypeIdentifier subsets_entry_memory_datatype_id;
	HDF5DatatypeIdentifier native_string_datatype_id;

	hsize_t file_dims[1]{0};

	try {
		subsets_entry_memory_datatype_id  = create_subsets_entry_memory_datatype();
	} catch (HVCFCreateException &e) {
		throw HVCFReadException(__FILE__, __FUNCTION__, __LINE__, "Error while creating memory datatype.");
	}

	if ((native_string_datatype_id = H5Topen(file_id, VARIABLE_LENGTH_STRING_TYPE, H5P_DEFAULT)) < 0) {
		throw HVCFReadException(__FILE__, __FUNCTION__, __LINE__, "Error while opening datatype.");
	}

	if ((dataset_id = H5Dopen(samples_group_id, SAMPLE_SUBSETS_DATASET, H5P_DEFAULT)) < 0) {
		throw HVCFReadException(__FILE__, __FUNCTION__, __LINE__, "Error while creating dataset.");
	}

	if ((file_dataspace_id = H5Dget_space(dataset_id)) < 0) {
		throw HVCFReadException(__FILE__, __FUNCTION__, __LINE__, "Error while getting dataspace.");
	}

	if (H5Sget_simple_extent_dims(file_dataspace_id, file_dims, nullptr) < 0) {
		throw HVCFReadException(__FILE__, __FUNCTION__, __LINE__, "Error while getting dataspace dimensions.");
	}

	subsets_entry_type subsets_entry_buffer[file_dims[0]];

	if (H5Dread(dataset_id, subsets_entry_memory_datatype_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, subsets_entry_buffer) < 0) {
		throw HVCFReadException(__FILE__, __FUNCTION__, __LINE__, "Error while reading from dataset");
	}

	vector<tuple<hsize_t, hsize_t>> population_offsets;
	unsigned int n_samples = 0u;

	for (hsize_t i = 0; i < file_dims[0]; ++i) {
		if (strcmp(subsets_entry_buffer[i].name, name.c_str()) == 0) {
			population_offsets.emplace_back(subsets_entry_buffer[i].offset_1, subsets_entry_buffer[i].offset_2);
			n_samples += subsets_entry_buffer[i].offset_2 - subsets_entry_buffer[i].offset_1 + 1u;
		}
	}

	if (H5Dvlen_reclaim(subsets_entry_memory_datatype_id, file_dataspace_id, H5P_DEFAULT, subsets_entry_buffer) < 0) {
		throw HVCFReadException(__FILE__, __FUNCTION__, __LINE__, "Error while reclaiming HDF5 memory.");
	}

	file_dataspace_id.close();
	dataset_id.close();

	hsize_t mem_dims[1]{n_samples};
	hsize_t file_offset[1]{0};
	hsize_t counts[1]{0};
	char* samples_buffer[n_samples];

	if ((dataset_id = H5Dopen(samples_group_id, SAMPLE_NAMES_DATASET, H5P_DEFAULT)) < 0) {
		throw HVCFReadException(__FILE__, __FUNCTION__, __LINE__, "Error while opening dataset.");
	}

	if ((file_dataspace_id = H5Dget_space(dataset_id)) < 0) {
		throw HVCFReadException(__FILE__, __FUNCTION__, __LINE__, "Error while getting dataspace.");
	}

	if ((memory_dataspace_id = H5Screate_simple(1, mem_dims, nullptr)) < 0) {
		throw HVCFReadException(__FILE__, __FUNCTION__, __LINE__, "Error while creating memory dataspace.");
	}

	if (H5Sselect_none(file_dataspace_id) < 0) {
		throw HVCFReadException(__FILE__, __FUNCTION__, __LINE__, "Error while making selection in dataspace.");
	}

	for (unsigned int i = 0; i < population_offsets.size(); ++i) {
		file_offset[0] = get<0>(population_offsets[i]);
		counts[0] = get<1>(population_offsets[i]) - get<0>(population_offsets[i]) + 1u;

		if (H5Sselect_hyperslab(file_dataspace_id, H5S_SELECT_OR, file_offset, NULL, counts, NULL) < 0) {
			throw HVCFReadException(__FILE__, __FUNCTION__, __LINE__, "Error while making selection in dataspace.");
		}
	}

	if (H5Dread(dataset_id, native_string_datatype_id, memory_dataspace_id, file_dataspace_id, H5P_DEFAULT, samples_buffer) < 0) {
		throw HVCFReadException(__FILE__, __FUNCTION__, __LINE__, "Error while reading from dataset.");
	}

	vector<string> samples;

	for (unsigned int i = 0; i < n_samples; ++i) {
		samples.emplace_back(samples_buffer[i]);
	}

	if (H5Dvlen_reclaim(native_string_datatype_id, memory_dataspace_id, H5P_DEFAULT, samples_buffer) < 0) {
		throw HVCFReadException(__FILE__, __FUNCTION__, __LINE__, "Error while reclaiming HDF5 memory.");
	}

	return samples;
}

unsigned int HVCF::get_n_chromosomes() const {
	return chromosomes.size();
}

vector<string> HVCF::get_chromosomes() const {
	vector<string> names;
	for (auto&& chromosome : chromosomes) {
		names.emplace_back(chromosome.first);
	}
	return names;
}

bool HVCF::has_chromosome(const string& chromosome) const {
	return chromosomes.count(chromosome) > 0;
}

unsigned long long int HVCF::get_chromosome_start(const string& chromosome) const throw (HVCFReadException) {
	auto chromosomes_cache_it = chromosomes_cache.find(chromosome);
	if (chromosomes_cache_it == chromosomes_cache.end()) {
		return 0;
	}

	HDF5DataspaceIdentifier dataspace_id;
	HDF5DataspaceIdentifier memory_dataspace_id;

	hsize_t file_dims[1]{0};

	if ((dataspace_id = H5Dget_space(chromosomes_cache_it->second->intervals_index_id)) < 0) {
		throw HVCFReadException(__FILE__, __FUNCTION__, __LINE__, "Error while getting dataspace.");
	}

	if (H5Sget_simple_extent_dims(dataspace_id, file_dims, NULL) < 0) {
		throw HVCFReadException(__FILE__, __FUNCTION__, __LINE__, "Error while getting dataspace dimensions.");
	}

	interval_index_entry_type interval_index_entries_buffer[file_dims[0]];

	if (H5Dread(chromosomes_cache_it->second->intervals_index_id, interval_index_entry_memory_datatype_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, interval_index_entries_buffer) < 0) {
		throw HVCFReadException(__FILE__, __FUNCTION__, __LINE__, "Error while reading from dataset.");
	}

	return interval_index_entries_buffer[0].ull_value_1;
}

unsigned long long int HVCF::get_chromosome_end(const string& chromosome) const throw (HVCFReadException) {
	auto chromosomes_cache_it = chromosomes_cache.find(chromosome);
	if (chromosomes_cache_it == chromosomes_cache.end()) {
		return 0;
	}

	HDF5DataspaceIdentifier dataspace_id;
	HDF5DataspaceIdentifier memory_dataspace_id;

	hsize_t file_dims[1]{0};

	if ((dataspace_id = H5Dget_space(chromosomes_cache_it->second->intervals_index_id)) < 0) {
		throw HVCFReadException(__FILE__, __FUNCTION__, __LINE__, "Error while getting dataspace.");
	}

	if (H5Sget_simple_extent_dims(dataspace_id, file_dims, NULL) < 0) {
		throw HVCFReadException(__FILE__, __FUNCTION__, __LINE__, "Error while getting dataspace dimensions.");
	}

	interval_index_entry_type interval_index_entries_buffer[file_dims[0]];

	if (H5Dread(chromosomes_cache_it->second->intervals_index_id, interval_index_entry_memory_datatype_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, interval_index_entries_buffer) < 0) {
		throw HVCFReadException(__FILE__, __FUNCTION__, __LINE__, "Error while reading from dataset.");
	}

	return interval_index_entries_buffer[file_dims[0] - 1].ull_value_2;
}

hsize_t HVCF::get_n_variants() const throw (HVCFReadException) {
	hsize_t total = 0;

	for (auto&& entry : chromosomes) {
		total += get_n_variants_in_chromosome(entry.first);
	}

	return total;
}

hsize_t HVCF::get_n_variants_in_chromosome(const string& chromosome) const throw (HVCFReadException) {
	auto chromosomes_it = chromosomes.find(chromosome);

	if (chromosomes_it == chromosomes.end()) {
		return 0;
	}

	HDF5DatasetIdentifier dataset_id;
	HDF5DataspaceIdentifier file_dataspace_id;

	hsize_t file_dims[1]{0};

	if ((dataset_id = H5Dopen(chromosomes_it->second->get(), VARIANTS_DATASET, H5P_DEFAULT)) < 0) {
		throw HVCFReadException(__FILE__, __FUNCTION__, __LINE__, "Error while opening dataset.");
	}

	if ((file_dataspace_id = H5Dget_space(dataset_id)) < 0) {
		throw HVCFReadException(__FILE__, __FUNCTION__, __LINE__, "Error while getting dataspace.");
	}

	if (H5Sget_simple_extent_dims(file_dataspace_id, file_dims, nullptr) < 0) {
		throw HVCFReadException(__FILE__, __FUNCTION__, __LINE__, "Error while getting dataspace dimensions.");
	}

	return file_dims[0];
}

long long int HVCF::get_variant_offset_by_position_eq(const string& chromosome, unsigned long long int position) throw (HVCFReadException) {
	auto chromosomes_cache_it = chromosomes_cache.find(chromosome);
	if (chromosomes_cache_it == chromosomes_cache.end()) {
		return -1;
	}

	HDF5DataspaceIdentifier dataspace_id;
	HDF5DataspaceIdentifier memory_dataspace_id;

	hsize_t file_dims[1]{0};

	if ((dataspace_id = H5Dget_space(chromosomes_cache_it->second->intervals_index_id)) < 0) {
		throw HVCFReadException(__FILE__, __FUNCTION__, __LINE__, "Error while getting dataspace.");
	}

	if (H5Sget_simple_extent_dims(dataspace_id, file_dims, NULL) < 0) {
		throw HVCFReadException(__FILE__, __FUNCTION__, __LINE__, "Error while getting dataspace dimensions.");
	}

	interval_index_entry_type interval_index_entries_buffer[file_dims[0]];

	if (H5Dread(chromosomes_cache_it->second->intervals_index_id, interval_index_entry_memory_datatype_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, interval_index_entries_buffer) < 0) {
		throw HVCFReadException(__FILE__, __FUNCTION__, __LINE__, "Error while reading from dataset.");
	}

	dataspace_id.close();

	interval_index_entry_type search_interval_index_entry;
	search_interval_index_entry.ull_value_2 = position;

	auto interval_index_entry = lower_bound(interval_index_entries_buffer, interval_index_entries_buffer + file_dims[0], search_interval_index_entry,
			[] (const interval_index_entry_type& f, const interval_index_entry_type& s) -> bool {
				return (f.ull_value_2 < s.ull_value_2);
			});

	if (interval_index_entry == interval_index_entries_buffer + file_dims[0]) {
		return -1;
	}

	if (position < interval_index_entry->ull_value_1) {
		return -1;
	}

	hsize_t mem_dims[1] {interval_index_entry->bucket_size};
	hsize_t offset[1] {interval_index_entry->bucket_offset};

	ull_index_entry_type ull_index_entries_buffer[interval_index_entry->bucket_size];

	if ((memory_dataspace_id = H5Screate_simple(1, mem_dims, nullptr)) < 0) {
		throw HVCFReadException(__FILE__, __FUNCTION__, __LINE__, "Error while creating memory dataspace.");
	}

	if ((dataspace_id = H5Dget_space(chromosomes_cache_it->second->intervals_index_buckets_id)) < 0) {
		throw HVCFReadException(__FILE__, __FUNCTION__, __LINE__, "Error while getting dataspace.");
	}

	if (H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, offset, NULL, mem_dims, NULL) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while making selection in dataspace.");
	}

	if (H5Dread(chromosomes_cache_it->second->intervals_index_buckets_id, ull_index_entry_memory_datatype_id, memory_dataspace_id, dataspace_id, H5P_DEFAULT, ull_index_entries_buffer) < 0) {
		throw HVCFReadException(__FILE__, __FUNCTION__, __LINE__, "Error while reading from dataset.");
	}

	ull_index_entry_type search_ull_index_entry;
	search_ull_index_entry.ull_value = position;

	auto ull_index_entry = lower_bound(ull_index_entries_buffer, ull_index_entries_buffer + mem_dims[0], search_ull_index_entry,
			[] (const ull_index_entry_type& f, const ull_index_entry_type& s) -> bool {
				return (f.ull_value < s.ull_value);
			});

	if ((ull_index_entry == ull_index_entries_buffer + mem_dims[0]) || (ull_index_entry->ull_value != position)) {
		return -1;
	}

	return ull_index_entry->offset;
}

long long int HVCF::get_variant_offset_by_position_ge(const string& chromosome, unsigned long long int position) throw (HVCFReadException) {
	auto chromosomes_cache_it = chromosomes_cache.find(chromosome);
	if (chromosomes_cache_it == chromosomes_cache.end()) {
		return -1;
	}

	HDF5DataspaceIdentifier dataspace_id;
	HDF5DataspaceIdentifier memory_dataspace_id;

	hsize_t file_dims[1]{0};

	if ((dataspace_id = H5Dget_space(chromosomes_cache_it->second->intervals_index_id)) < 0) {
		throw HVCFReadException(__FILE__, __FUNCTION__, __LINE__, "Error while getting dataspace.");
	}

	if (H5Sget_simple_extent_dims(dataspace_id, file_dims, NULL) < 0) {
		throw HVCFReadException(__FILE__, __FUNCTION__, __LINE__, "Error while getting dataspace dimensions.");
	}

	interval_index_entry_type interval_index_entries_buffer[file_dims[0]];

	if (H5Dread(chromosomes_cache_it->second->intervals_index_id, interval_index_entry_memory_datatype_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, interval_index_entries_buffer) < 0) {
		throw HVCFReadException(__FILE__, __FUNCTION__, __LINE__, "Error while reading from dataset.");
	}

	dataspace_id.close();

	interval_index_entry_type search_interval_index_entry;
	search_interval_index_entry.ull_value_2 = position;

	auto interval_index_entry = lower_bound(interval_index_entries_buffer, interval_index_entries_buffer + file_dims[0], search_interval_index_entry,
			[] (const interval_index_entry_type& f, const interval_index_entry_type& s) -> bool {
				return (f.ull_value_2 < s.ull_value_2);
			});

	if (interval_index_entry == interval_index_entries_buffer + file_dims[0]) {
		return -1;
	}

	if (position <= interval_index_entry->ull_value_1) {
		return interval_index_entry->offset_1;
	}

	hsize_t offset[1] {interval_index_entry->bucket_offset};
	hsize_t mem_dims[1] {interval_index_entry->bucket_size};

	ull_index_entry_type ull_index_entries_buffer[interval_index_entry->bucket_size];

	if ((memory_dataspace_id = H5Screate_simple(1, mem_dims, nullptr)) < 0) {
		throw HVCFReadException(__FILE__, __FUNCTION__, __LINE__, "Error while creating memory dataspace.");
	}

	if ((dataspace_id = H5Dget_space(chromosomes_cache_it->second->intervals_index_buckets_id)) < 0) {
		throw HVCFReadException(__FILE__, __FUNCTION__, __LINE__, "Error while getting dataspace.");
	}

	if (H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, offset, NULL, mem_dims, NULL) < 0) {
		throw HVCFReadException(__FILE__, __FUNCTION__, __LINE__, "Error while making selection in dataspace.");
	}

	if (H5Dread(chromosomes_cache_it->second->intervals_index_buckets_id, ull_index_entry_memory_datatype_id, memory_dataspace_id, dataspace_id, H5P_DEFAULT, ull_index_entries_buffer) < 0) {
		throw HVCFReadException(__FILE__, __FUNCTION__, __LINE__, "Error while reading from dataset.");
	}

	ull_index_entry_type search_ull_index_entry;
	search_ull_index_entry.ull_value = position;

	auto ull_index_entry = lower_bound(ull_index_entries_buffer, ull_index_entries_buffer + mem_dims[0], search_ull_index_entry,
			[] (const ull_index_entry_type& f, const ull_index_entry_type& s) -> bool {
				return (f.ull_value < s.ull_value);
			});

	if (ull_index_entry == ull_index_entries_buffer + mem_dims[0]) {
		return -1;
	}

	return ull_index_entry->offset;
}

long long int HVCF::get_variant_offset_by_position_le(const string& chromosome, unsigned long long int position) throw (HVCFReadException) {
	auto chromosomes_cache_it = chromosomes_cache.find(chromosome);
	if (chromosomes_cache_it == chromosomes_cache.end()) {
		return -1;
	}

	HDF5DataspaceIdentifier dataspace_id;
	HDF5DataspaceIdentifier memory_dataspace_id;

	hsize_t file_dims[1]{0};

	if ((dataspace_id = H5Dget_space(chromosomes_cache_it->second->intervals_index_id)) < 0) {
		throw HVCFReadException(__FILE__, __FUNCTION__, __LINE__, "Error while getting dataspace.");
	}

	if (H5Sget_simple_extent_dims(dataspace_id, file_dims, NULL) < 0) {
		throw HVCFReadException(__FILE__, __FUNCTION__, __LINE__, "Error while getting dataspace dimensions.");
	}

	interval_index_entry_type interval_index_entries_buffer[file_dims[0]];

	if (H5Dread(chromosomes_cache_it->second->intervals_index_id, interval_index_entry_memory_datatype_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, interval_index_entries_buffer) < 0) {
		throw HVCFReadException(__FILE__, __FUNCTION__, __LINE__, "Error while reading from dataset.");
	}

	dataspace_id.close();

	interval_index_entry_type search_interval_index_entry;
	search_interval_index_entry.ull_value_2 = position;

	auto interval_index_entry = lower_bound(interval_index_entries_buffer, interval_index_entries_buffer + file_dims[0], search_interval_index_entry,
			[] (const interval_index_entry_type& f, const interval_index_entry_type& s) -> bool {
				return (f.ull_value_2 < s.ull_value_2);
			});

	if ((interval_index_entry == interval_index_entries_buffer) && (position < interval_index_entry->ull_value_1)) {
		return -1;
	}

	if (interval_index_entry == interval_index_entries_buffer + file_dims[0]) {
		return interval_index_entries_buffer[file_dims[0] - 1].offset_2;
	}

	if (position < interval_index_entry->ull_value_1) {
		return (--interval_index_entry)->offset_2;
	}

	if (position == interval_index_entry->ull_value_1) {
		return interval_index_entry->offset_1;
	}

	hsize_t offset[1] {interval_index_entry->bucket_offset};
	hsize_t mem_dims[1] {interval_index_entry->bucket_size};

	ull_index_entry_type ull_index_entries_buffer[interval_index_entry->bucket_size];

	if ((memory_dataspace_id = H5Screate_simple(1, mem_dims, nullptr)) < 0) {
		throw HVCFReadException(__FILE__, __FUNCTION__, __LINE__, "Error while creating memory dataspace.");
	}

	if ((dataspace_id = H5Dget_space(chromosomes_cache_it->second->intervals_index_buckets_id.get())) < 0) {
		throw HVCFReadException(__FILE__, __FUNCTION__, __LINE__, "Error while getting dataspace.");
	}

	if (H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, offset, NULL, mem_dims, NULL) < 0) {
		throw HVCFReadException(__FILE__, __FUNCTION__, __LINE__, "Error while making selection in dataspace.");
	}

	if (H5Dread(chromosomes_cache_it->second->intervals_index_buckets_id.get(), ull_index_entry_memory_datatype_id, memory_dataspace_id, dataspace_id, H5P_DEFAULT, ull_index_entries_buffer) < 0) {
		throw HVCFReadException(__FILE__, __FUNCTION__, __LINE__, "Error while reading from dataset.");
	}

	ull_index_entry_type search_ull_index_entry;
	search_ull_index_entry.ull_value = position;

	auto ull_index_entry = lower_bound(ull_index_entries_buffer, ull_index_entries_buffer + mem_dims[0], search_ull_index_entry,
			[] (const ull_index_entry_type& f, const ull_index_entry_type& s) -> bool {
				return (f.ull_value < s.ull_value);
			});

	if ((ull_index_entry == ull_index_entries_buffer) && (position < ull_index_entry->ull_value)) {
		return -1;
	}

	if (ull_index_entry == ull_index_entries_buffer + mem_dims[0]) {
		return ull_index_entries_buffer[mem_dims[0] - 1].offset;
	}

	if (position < ull_index_entry->ull_value) {
		return (--ull_index_entry)->offset;
	}

	return ull_index_entry->offset;
}

long long int HVCF::get_variant_offset_by_name(const string& chromosome, const string& name) throw (HVCFReadException) {
	auto chromosomes_cache_it = chromosomes_cache.find(chromosome);
	if (chromosomes_cache_it == chromosomes_cache.end()) {
		return -1;
	}

	HDF5DataspaceIdentifier dataspace_id;
	HDF5DataspaceIdentifier memory_dataspace_id;

	std::hash<string> hash_function;

	hsize_t hash_value = hash_function(name) % N_VARIANTS_HASH_BUCKETS;

	hsize_t offset[1]{hash_value};
	hsize_t mem_dims[1]{1};

	hash_index_entry_type hash_index_entry;

	if ((dataspace_id = H5Dget_space(chromosomes_cache_it->second->names_index_id)) < 0) {
		throw HVCFReadException(__FILE__, __FUNCTION__, __LINE__, "Error while getting dataspace.");
	}

	if ((memory_dataspace_id = H5Screate_simple(1, mem_dims, nullptr)) < 0) {
		throw HVCFReadException(__FILE__, __FUNCTION__, __LINE__, "Error while creating memory dataspace.");
	}

	if (H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, offset, NULL, mem_dims, NULL) < 0) {
		throw HVCFReadException(__FILE__, __FUNCTION__, __LINE__, "Error while making selection in dataspace.");
	}

	if (H5Dread(chromosomes_cache_it->second->names_index_id, hash_index_entry_memory_datatype_id, memory_dataspace_id, dataspace_id, H5P_DEFAULT, &hash_index_entry) < 0) {
		throw HVCFReadException(__FILE__, __FUNCTION__, __LINE__, "Error while reading from dataset");
	}

	dataspace_id.close();
	memory_dataspace_id.close();

	if ((hash_index_entry.bucket_offset == 0) && (hash_index_entry.bucket_size == 0)) {
		return -1;
	}

	offset[0] = hash_index_entry.bucket_offset;
	mem_dims[0] = hash_index_entry.bucket_size;

	string_index_entry_type string_index_entries_buffer[hash_index_entry.bucket_size];

	if ((dataspace_id = H5Dget_space(chromosomes_cache_it->second->names_index_buckets_id)) < 0) {
		throw HVCFReadException(__FILE__, __FUNCTION__, __LINE__, "Error while getting dataspace.");
	}

	if ((memory_dataspace_id = H5Screate_simple(1, mem_dims, nullptr)) < 0) {
		throw HVCFReadException(__FILE__, __FUNCTION__, __LINE__, "Error while creating memory dataspace.");
	}

	if (H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, offset, NULL, mem_dims, NULL) < 0) {
		throw HVCFReadException(__FILE__, __FUNCTION__, __LINE__, "Error while making selection in dataspace.");
	}

	if (H5Dread(chromosomes_cache_it->second->names_index_buckets_id, string_index_entry_memory_datatype_id, memory_dataspace_id, dataspace_id, H5P_DEFAULT, string_index_entries_buffer) < 0) {
		throw HVCFReadException(__FILE__, __FUNCTION__, __LINE__, "Error while reading from dataset.");
	}

	string_index_entry_type search_string_index_entry;
	search_string_index_entry.string_value = (char*)name.c_str();

	auto string_index_entry = lower_bound(string_index_entries_buffer, string_index_entries_buffer + mem_dims[0], search_string_index_entry,
			[] (const string_index_entry_type& f, const string_index_entry_type& s) -> bool {
				return (strcmp(f.string_value, s.string_value) < 0);
			});

	long long int variant_offset = 0;

	if ((string_index_entry == string_index_entries_buffer + mem_dims[0]) || (strcmp(string_index_entry->string_value, name.c_str())) != 0) {
		variant_offset = -1;
	} else {
		variant_offset = string_index_entry->offset;
	}

	if (H5Dvlen_reclaim(string_index_entry_memory_datatype_id, memory_dataspace_id, H5P_DEFAULT, string_index_entries_buffer) < 0) {
		throw HVCFReadException(__FILE__, __FUNCTION__, __LINE__, "Error while reclaiming HDF5 memory.");
	}

	return variant_offset;
}

long long int HVCF::get_sample_offset(const string& name) throw (HVCFReadException) {
	HDF5DataspaceIdentifier dataspace_id;
	HDF5DataspaceIdentifier memory_dataspace_id;

	std::hash<string> hash_function;

	hsize_t hash_value = hash_function(name) % N_SAMPLES_HASH_BUCKETS;

	hsize_t offset[1]{hash_value};
	hsize_t mem_dims[1]{1};

	hash_index_entry_type hash_index_entry;

	if ((dataspace_id = H5Dget_space(samples_cache.names_index_id)) < 0) {
		throw HVCFReadException(__FILE__, __FUNCTION__, __LINE__, "Error while getting dataspace.");
	}

	if ((memory_dataspace_id = H5Screate_simple(1, mem_dims, nullptr)) < 0) {
		throw HVCFReadException(__FILE__, __FUNCTION__, __LINE__, "Error while creating memory dataspace.");
	}

	if (H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, offset, NULL, mem_dims, NULL) < 0) {
		throw HVCFReadException(__FILE__, __FUNCTION__, __LINE__, "Error while making selection in dataspace.");
	}

	if (H5Dread(samples_cache.names_index_id, hash_index_entry_memory_datatype_id, memory_dataspace_id, dataspace_id, H5P_DEFAULT, &hash_index_entry) < 0) {
		throw HVCFReadException(__FILE__, __FUNCTION__, __LINE__, "Error while reading from dataset");
	}

	dataspace_id.close();
	memory_dataspace_id.close();

	if ((hash_index_entry.bucket_offset == 0) && (hash_index_entry.bucket_size == 0)) {
		return -1;
	}

	offset[0] = hash_index_entry.bucket_offset;
	mem_dims[0] = hash_index_entry.bucket_size;

	string_index_entry_type string_index_entries_buffer[hash_index_entry.bucket_size];

	if ((dataspace_id = H5Dget_space(samples_cache.names_index_buckets_id)) < 0) {
		throw HVCFReadException(__FILE__, __FUNCTION__, __LINE__, "Error while getting dataspace.");
	}

	if ((memory_dataspace_id = H5Screate_simple(1, mem_dims, nullptr)) < 0) {
		throw HVCFReadException(__FILE__, __FUNCTION__, __LINE__, "Error while creating memory dataspace.");
	}

	if (H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, offset, NULL, mem_dims, NULL) < 0) {
		throw HVCFReadException(__FILE__, __FUNCTION__, __LINE__, "Error while making selection in dataspace.");
	}

	if (H5Dread(samples_cache.names_index_buckets_id, string_index_entry_memory_datatype_id, memory_dataspace_id, dataspace_id, H5P_DEFAULT, string_index_entries_buffer) < 0) {
		throw HVCFReadException(__FILE__, __FUNCTION__, __LINE__, "Error while reading from dataset.");
	}

	string_index_entry_type search_string_index_entry;
	search_string_index_entry.string_value = (char*)name.c_str();

	auto string_index_entry = lower_bound(string_index_entries_buffer, string_index_entries_buffer + mem_dims[0], search_string_index_entry,
			[] (const string_index_entry_type& f, const string_index_entry_type& s) -> bool {
				return (strcmp(f.string_value, s.string_value) < 0);
			});

	long long int variant_offset = 0;

	if ((string_index_entry == string_index_entries_buffer + mem_dims[0]) || (strcmp(string_index_entry->string_value, name.c_str())) != 0) {
		variant_offset = -1;
	} else {
		variant_offset = string_index_entry->offset;
	}

	if (H5Dvlen_reclaim(string_index_entry_memory_datatype_id, memory_dataspace_id, H5P_DEFAULT, string_index_entries_buffer) < 0) {
		throw HVCFReadException(__FILE__, __FUNCTION__, __LINE__, "Error while reclaiming HDF5 memory.");
	}

	return variant_offset;
}

void HVCF::compute_ld(const string& chromosome, const string& subset, unsigned long long int start_position, unsigned long long end_position, vector<variants_pair>& result) throw (HVCFReadException) {
	auto chromosomes_cache_it = chromosomes_cache.find(chromosome);
	if (chromosomes_cache_it == chromosomes_cache.end()) {
		return;
	}

	if (end_position < start_position) {
		return;
	}

	long long int start_position_offset = 0;
	long long int end_position_offset = 0;

	if ((start_position_offset = get_variant_offset_by_position_eq(chromosome, start_position)) < 0) {
		return;
	}

	if ((end_position_offset = get_variant_offset_by_position_eq(chromosome, end_position)) < 0) {
		return;
	}

	auto subsets_cache_it = samples_cache.subsets.find(subset);

	if (subsets_cache_it == samples_cache.subsets.end()) {
		return;
	}

	HDF5DataspaceIdentifier file_dataspace_id;
	HDF5DataspaceIdentifier memory_dataspace_id;

	hsize_t n_samples = subsets_cache_it->second.n_samples;
	hsize_t n_haplotypes = 2 * n_samples;

	hsize_t n_variants = end_position_offset - start_position_offset + 1;

	hsize_t file_offset_2D[2]{static_cast<hsize_t>(start_position_offset), 0};
	hsize_t counts_2D[2]{n_variants, 0};
	hsize_t mem_dims_2D[2]{n_variants, n_haplotypes};

	unique_ptr<double[]> haplotypes = unique_ptr<double[]>(new double[n_variants * n_haplotypes]);

	if ((file_dataspace_id = H5Dget_space(chromosomes_cache_it->second->haplotypes_id)) < 0) {
		throw HVCFReadException(__FILE__, __FUNCTION__, __LINE__, "Error while getting dataspace.");
	}

	if ((memory_dataspace_id = H5Screate_simple(2, mem_dims_2D, nullptr)) < 0) {
		throw HVCFReadException(__FILE__, __FUNCTION__, __LINE__, "Error while creating memory dataspace.");
	}

	if (H5Sselect_none(file_dataspace_id) < 0) {
		throw HVCFReadException(__FILE__, __FUNCTION__, __LINE__, "Error while making selection in dataspace.");
	}

	for (auto& chunk : subsets_cache_it->second.chunks) {
		file_offset_2D[1] = 2 * get<0>(chunk);
		counts_2D[1] = 2 * get<2>(chunk);

		if (H5Sselect_hyperslab(file_dataspace_id, H5S_SELECT_OR, file_offset_2D, NULL, counts_2D, NULL) < 0) {
			throw HVCFReadException(__FILE__, __FUNCTION__, __LINE__, "Error while making selection in dataspace.");
		}
	}

	if (H5Dread(chromosomes_cache_it->second->haplotypes_id, H5T_NATIVE_UCHAR, memory_dataspace_id, file_dataspace_id, H5P_DEFAULT, haplotypes.get()) < 0) {
		throw HVCFReadException(__FILE__, __FUNCTION__, __LINE__, "Error while reading from dataset");
	}

	file_dataspace_id.close();
	memory_dataspace_id.close();

	if (H5Tconvert(H5T_NATIVE_UCHAR, H5T_NATIVE_DOUBLE, n_variants * n_haplotypes, haplotypes.get(), nullptr, H5P_DEFAULT) < 0) {
		throw HVCFReadException(__FILE__, __FUNCTION__, __LINE__, "Error while converting datatypes.");
	}

	Mat<double> S(haplotypes.get(), n_haplotypes, n_variants, false, false); // this call doesn't copy matrix.
	Row<double> J(n_haplotypes, fill::ones);

	Mat<double> C1(J * S);
	Mat<double> C2(n_haplotypes - C1);
	Mat<double> M1(C1.t() * C1);

	Mat<double> R((n_haplotypes * S.t() * S - M1) / sqrt(M1 % (C2.t() * C2)));

//	cout << "R:" << endl;
//	R.raw_print();

//	cout << "Rsq:" << endl;
//	cout << square(R) << endl;

	hsize_t file_offset_1D[1]{static_cast<hsize_t>(start_position_offset)};
	hsize_t mem_dims_1D[1]{n_variants};

	variants_entry_type variants_buffer[n_variants];

	if ((file_dataspace_id = H5Dget_space(chromosomes_cache_it->second->variants_id)) < 0) {
		throw HVCFReadException(__FILE__, __FUNCTION__, __LINE__, "Error while getting dataspace.");
	}

	if ((memory_dataspace_id = H5Screate_simple(1, mem_dims_1D, nullptr)) < 0) {
		throw HVCFReadException(__FILE__, __FUNCTION__, __LINE__, "Error while creating memory dataspace.");
	}

	if (H5Sselect_hyperslab(file_dataspace_id, H5S_SELECT_SET, file_offset_1D, NULL, mem_dims_1D, NULL) < 0) {
		throw HVCFReadException(__FILE__, __FUNCTION__, __LINE__, "Error while making selection in dataspace.");
	}

	if (H5Dread(chromosomes_cache_it->second->variants_id, variants_entry_memory_datatype_id, memory_dataspace_id, file_dataspace_id, H5P_DEFAULT, variants_buffer) < 0) {
		throw HVCFReadException(__FILE__, __FUNCTION__, __LINE__, "Error while reading from dataset");
	}

	for (unsigned int i = 0u; i < n_variants; ++i) {
		for (unsigned int j = 0u; j < n_variants; ++j) {
			result.emplace_back(
				variants_buffer[i].name, variants_buffer[i].position,
				variants_buffer[j].name, variants_buffer[j].position,
				R(i, j), pow(R(i, j), 2.0));
		}
	}

	if (H5Dvlen_reclaim(variants_entry_memory_datatype_id, memory_dataspace_id, H5P_DEFAULT, variants_buffer) < 0) {
		throw HVCFReadException(__FILE__, __FUNCTION__, __LINE__, "Error while reclaiming HDF5 memory.");
	}
}

void HVCF::compute_ld(const string& chromosome, const string& subset, const string& lead_variant_name, unsigned long long int start_position, unsigned long long end_position, vector<variants_pair>& result) throw (HVCFReadException) {
	auto chromosomes_cache_it = chromosomes_cache.find(chromosome);
	if (chromosomes_cache_it == chromosomes_cache.end()) {
		return;
	}

	if (end_position < start_position) {
		return;
	}

	long long int lead_variant_offset = 0;
	long long int start_position_offset = 0;
	long long int end_position_offset = 0;

	if ((lead_variant_offset = get_variant_offset_by_name(chromosome, lead_variant_name)) < 0) {
		return;
	}

	if ((start_position_offset = get_variant_offset_by_position_ge(chromosome, start_position)) < 0) {
		return;
	}

	if ((end_position_offset = get_variant_offset_by_position_le(chromosome, end_position)) < 0) {
		return;
	}

	auto subsets_cache_it = samples_cache.subsets.find(subset);

	if (subsets_cache_it == samples_cache.subsets.end()) {
		return;
	}

	HDF5DataspaceIdentifier file_dataspace_id;
	HDF5DataspaceIdentifier memory_dataspace_id;

	hsize_t n_samples = subsets_cache_it->second.n_samples;
	hsize_t n_haplotypes = 2 * n_samples;

	hsize_t n_variants = 0;
	hsize_t lead_variant_local_offset = 0;

	if ((lead_variant_offset >= start_position_offset) && (lead_variant_offset <= end_position_offset)) {
		n_variants = end_position_offset - start_position_offset + 1;
		lead_variant_local_offset  = lead_variant_offset - start_position_offset;
	} else {
		n_variants = end_position_offset - start_position_offset + 2;
		if (lead_variant_offset < start_position_offset) {
			lead_variant_local_offset = 0;
		} else {
			lead_variant_local_offset = n_variants - 1;
		}
	}

	hsize_t file_offset1_2D[2]{static_cast<hsize_t>(lead_variant_offset), 0};
	hsize_t counts1_2D[2]{1, n_haplotypes};
	hsize_t file_offset2_2D[2]{static_cast<hsize_t>(start_position_offset), 0};
	hsize_t counts2_2D[2]{static_cast<hsize_t>(end_position_offset - start_position_offset + 1), n_haplotypes};
	hsize_t mem_dims_2D[2]{n_variants, n_haplotypes};

	unique_ptr<double[]> haplotypes = unique_ptr<double[]>(new double[n_variants * n_haplotypes]);

	if ((file_dataspace_id = H5Dget_space(chromosomes_cache_it->second->haplotypes_id)) < 0) {
		throw HVCFReadException(__FILE__, __FUNCTION__, __LINE__, "Error while getting dataspace.");
	}

	if ((memory_dataspace_id = H5Screate_simple(2, mem_dims_2D, nullptr)) < 0) {
		throw HVCFReadException(__FILE__, __FUNCTION__, __LINE__, "Error while creating memory dataspace.");
	}

	if (H5Sselect_none(file_dataspace_id) < 0) {
		throw HVCFReadException(__FILE__, __FUNCTION__, __LINE__, "Error while making selection in dataspace.");
	}

	for (auto& chunk : subsets_cache_it->second.chunks) {
		file_offset1_2D[1] = 2 * get<0>(chunk);
		counts1_2D[1] = 2 * get<2>(chunk);

		file_offset2_2D[1] = 2 * get<0>(chunk);
		counts2_2D[1] = 2 * get<2>(chunk);

		if (H5Sselect_hyperslab(file_dataspace_id, H5S_SELECT_OR, file_offset1_2D, NULL, counts1_2D, NULL) < 0) {
			throw HVCFReadException(__FILE__, __FUNCTION__, __LINE__, "Error while making selection in dataspace.");
		}

		if (H5Sselect_hyperslab(file_dataspace_id, H5S_SELECT_OR, file_offset2_2D, NULL, counts2_2D, NULL) < 0) {
			throw HVCFReadException(__FILE__, __FUNCTION__, __LINE__, "Error while making selection in dataspace.");
		}
	}

	if (H5Dread(chromosomes_cache_it->second->haplotypes_id, H5T_NATIVE_UCHAR, memory_dataspace_id, file_dataspace_id, H5P_DEFAULT, haplotypes.get()) < 0) {
		throw HVCFReadException(__FILE__, __FUNCTION__, __LINE__, "Error while reading from dataset");
	}

	file_dataspace_id.close();
	memory_dataspace_id.close();

	if (H5Tconvert(H5T_NATIVE_UCHAR, H5T_NATIVE_DOUBLE, n_variants * n_haplotypes, haplotypes.get(), nullptr, H5P_DEFAULT) < 0) {
		throw HVCFReadException(__FILE__, __FUNCTION__, __LINE__, "Error while converting datatypes.");
	}

	Mat<double> S(haplotypes.get(), n_haplotypes, n_variants, false, false); // this call doesn't copy matrix.
	Row<double> L(S.colptr(lead_variant_local_offset), n_haplotypes);
	Row<double> J(n_haplotypes, fill::ones);

	Mat<double> LC1(J * L.t());
	Mat<double> LC2(n_haplotypes - LC1);
	Mat<double> SC1(J * S);
	Mat<double> SC2(n_haplotypes - SC1);
	Mat<double> M1(LC1 * SC1);

	Mat<double> R((n_haplotypes * L * S - M1) / sqrt(M1 % (LC2 * SC2)));

//	cout << "R:" << endl;
//	R.raw_print();

//	cout << "Rsq:" << endl;
//	cout << square(R) << endl;

	hsize_t file_offset1_1D[1]{static_cast<hsize_t>(lead_variant_offset)};
	hsize_t counts1_1D[1]{1};
	hsize_t file_offset2_1D[1]{static_cast<hsize_t>(start_position_offset)};
	hsize_t counts2_1D[1]{static_cast<hsize_t>(end_position_offset - start_position_offset + 1)};
	hsize_t mem_dims_1D[1]{n_variants};

	variants_entry_type variants_buffer[n_variants];

	if ((file_dataspace_id = H5Dget_space(chromosomes_cache_it->second->variants_id)) < 0) {
		throw HVCFReadException(__FILE__, __FUNCTION__, __LINE__, "Error while getting dataspace.");
	}

	if ((memory_dataspace_id = H5Screate_simple(1, mem_dims_1D, nullptr)) < 0) {
		throw HVCFReadException(__FILE__, __FUNCTION__, __LINE__, "Error while creating memory dataspace.");
	}

	if (H5Sselect_hyperslab(file_dataspace_id, H5S_SELECT_SET, file_offset1_1D, NULL, counts1_1D, NULL) < 0) {
		throw HVCFReadException(__FILE__, __FUNCTION__, __LINE__, "Error while making selection in dataspace.");
	}

	if (H5Sselect_hyperslab(file_dataspace_id, H5S_SELECT_OR, file_offset2_1D, NULL, counts2_1D, NULL) < 0) {
		throw HVCFReadException(__FILE__, __FUNCTION__, __LINE__, "Error while making selection in dataspace.");
	}

	if (H5Dread(chromosomes_cache_it->second->variants_id, variants_entry_memory_datatype_id, memory_dataspace_id, file_dataspace_id, H5P_DEFAULT, variants_buffer) < 0) {
		throw HVCFReadException(__FILE__, __FUNCTION__, __LINE__, "Error while reading from dataset");
	}

	for (unsigned int i = 0; i < lead_variant_local_offset; ++i) {
		result.emplace_back(
				variants_buffer[lead_variant_local_offset].name, variants_buffer[lead_variant_local_offset].position,
				variants_buffer[i].name, variants_buffer[i].position,
				R.at(0, i), pow(R.at(0, i), 2.0));
	}

	for (unsigned int i = lead_variant_local_offset + 1; i < n_variants; ++i) {
		result.emplace_back(
				variants_buffer[lead_variant_local_offset].name, variants_buffer[lead_variant_local_offset].position,
				variants_buffer[i].name, variants_buffer[i].position,
				R.at(0, i), pow(R.at(0, i), 2.0));
	}

	if (H5Dvlen_reclaim(variants_entry_memory_datatype_id, memory_dataspace_id, H5P_DEFAULT, variants_buffer) < 0) {
		throw HVCFReadException(__FILE__, __FUNCTION__, __LINE__, "Error while reclaiming HDF5 memory.");
	}

	return;
}

void HVCF::extract_variants(const string& chromosome, unsigned long long int start_position, unsigned long long int end_position, vector<variant_info>& result) throw (HVCFReadException) {
	auto chromosomes_cache_it = chromosomes_cache.find(chromosome);
	if (chromosomes_cache_it == chromosomes_cache.end()) {
		return;
	}

	if (end_position < start_position) {
		return;
	}

	long long int start_position_offset = 0;
	long long int end_position_offset = 0;

	if ((start_position_offset = get_variant_offset_by_position_eq(chromosome, start_position)) < 0) {
		return;
	}

	if ((end_position_offset = get_variant_offset_by_position_eq(chromosome, end_position)) < 0) {
		return;
	}

	HDF5DataspaceIdentifier file_dataspace_id;
	HDF5DataspaceIdentifier memory_dataspace_id;

	hsize_t n_variants = end_position_offset - start_position_offset + 1;

	hsize_t file_offset[1]{static_cast<hsize_t>(start_position_offset)};
	hsize_t mem_dims[1]{n_variants};

	variants_entry_type variants_buffer[n_variants];

	if ((file_dataspace_id = H5Dget_space(chromosomes_cache_it->second->variants_id)) < 0) {
		throw HVCFReadException(__FILE__, __FUNCTION__, __LINE__, "Error while getting dataspace.");
	}

	if ((memory_dataspace_id = H5Screate_simple(1, mem_dims, nullptr)) < 0) {
		throw HVCFReadException(__FILE__, __FUNCTION__, __LINE__, "Error while creating memory dataspace.");
	}

	if (H5Sselect_hyperslab(file_dataspace_id, H5S_SELECT_SET, file_offset, NULL, mem_dims, NULL) < 0) {
		throw HVCFReadException(__FILE__, __FUNCTION__, __LINE__, "Error while making selection in dataspace.");
	}

	if (H5Dread(chromosomes_cache_it->second->variants_id, variants_entry_memory_datatype_id, memory_dataspace_id, file_dataspace_id, H5P_DEFAULT, variants_buffer) < 0) {
		throw HVCFReadException(__FILE__, __FUNCTION__, __LINE__, "Error while reading from dataset");
	}

	for (unsigned int i = 0u; i < n_variants; ++i) {
		result.emplace_back(variants_buffer[i].name, variants_buffer[i].ref, variants_buffer[i].alt, variants_buffer[i].position);
	}

	if (H5Dvlen_reclaim(variants_entry_memory_datatype_id, memory_dataspace_id, H5P_DEFAULT, variants_buffer) < 0) {
		throw HVCFReadException(__FILE__, __FUNCTION__, __LINE__, "Error while reclaiming HDF5 memory.");
	}
}

unsigned int HVCF::get_n_opened_objects() const {
	if (file_id >= 0) {
		return H5Fget_obj_count(file_id, H5F_OBJ_ALL);
	}
	return 0u;
}

unsigned int HVCF::get_n_all_opened_objects() {
	return H5Fget_obj_count(H5F_OBJ_ALL, H5F_OBJ_ALL);
}

void HVCF::compute_ld_test(const string& chromosome, const string& subset, const string& lead_variant_name, unsigned long long int start_position, unsigned long long end_position, vector<variants_pair>& result) throw (HVCFReadException) {
	std::chrono::time_point<std::chrono::system_clock> start, end;
	std::chrono::duration<double> elapsed_seconds;

	auto chromosomes_cache_it = chromosomes_cache.find(chromosome);
	if (chromosomes_cache_it == chromosomes_cache.end()) {
		return;
	}

	if (end_position < start_position) {
		return;
	}

	long long int lead_variant_offset = 0;
	long long int start_position_offset = 0;
	long long int end_position_offset = 0;

	start = std::chrono::system_clock::now();

	if ((lead_variant_offset = get_variant_offset_by_name(chromosome, lead_variant_name)) < 0) {
		return;
	}

	if ((start_position_offset = get_variant_offset_by_position_ge(chromosome, start_position)) < 0) {
		return;
	}

	if ((end_position_offset = get_variant_offset_by_position_le(chromosome, end_position)) < 0) {
		return;
	}

	end = std::chrono::system_clock::now();
	elapsed_seconds = end - start;
	cout << "Got offsets in " << elapsed_seconds.count() << " seconds" << endl;

	start = std::chrono::system_clock::now();

	auto subsets_cache_it = samples_cache.subsets.find(subset);

	if (subsets_cache_it == samples_cache.subsets.end()) {
		return;
	}

	HDF5DataspaceIdentifier file_dataspace_id;
	HDF5DataspaceIdentifier memory_dataspace_id;

	hsize_t n_samples = subsets_cache_it->second.n_samples;
	hsize_t n_haplotypes = 2 * n_samples;

	hsize_t n_variants = 0;
	hsize_t lead_variant_local_offset = 0;

	if ((lead_variant_offset >= start_position_offset) && (lead_variant_offset <= end_position_offset)) {
		n_variants = end_position_offset - start_position_offset + 1;
		lead_variant_local_offset  = lead_variant_offset - start_position_offset;
	} else {
		n_variants = end_position_offset - start_position_offset + 2;
		if (lead_variant_offset < start_position_offset) {
			lead_variant_local_offset = 0;
		} else {
			lead_variant_local_offset = n_variants - 1;
		}
	}

	hsize_t file_offset1_2D[2]{static_cast<hsize_t>(lead_variant_offset), 0};
	hsize_t counts1_2D[2]{1, n_haplotypes};
	hsize_t file_offset2_2D[2]{static_cast<hsize_t>(start_position_offset), 0};
	hsize_t counts2_2D[2]{static_cast<hsize_t>(end_position_offset - start_position_offset + 1), n_haplotypes};
	hsize_t mem_dims_2D[2]{n_variants, n_haplotypes};

	unique_ptr<double[]> haplotypes = unique_ptr<double[]>(new double[n_variants * n_haplotypes]);

	if ((file_dataspace_id = H5Dget_space(chromosomes_cache_it->second->haplotypes_id)) < 0) {
		throw HVCFReadException(__FILE__, __FUNCTION__, __LINE__, "Error while getting dataspace.");
	}

	if ((memory_dataspace_id = H5Screate_simple(2, mem_dims_2D, nullptr)) < 0) {
		throw HVCFReadException(__FILE__, __FUNCTION__, __LINE__, "Error while creating memory dataspace.");
	}

	if (H5Sselect_none(file_dataspace_id) < 0) {
		throw HVCFReadException(__FILE__, __FUNCTION__, __LINE__, "Error while making selection in dataspace.");
	}

	for (auto& chunk : subsets_cache_it->second.chunks) {
		file_offset1_2D[1] = 2 * get<0>(chunk);
		counts1_2D[1] = 2 * get<2>(chunk);

		file_offset2_2D[1] = 2 * get<0>(chunk);
		counts2_2D[1] = 2 * get<2>(chunk);

		if (H5Sselect_hyperslab(file_dataspace_id, H5S_SELECT_OR, file_offset1_2D, NULL, counts1_2D, NULL) < 0) {
			throw HVCFReadException(__FILE__, __FUNCTION__, __LINE__, "Error while making selection in dataspace.");
		}

		if (H5Sselect_hyperslab(file_dataspace_id, H5S_SELECT_OR, file_offset2_2D, NULL, counts2_2D, NULL) < 0) {
			throw HVCFReadException(__FILE__, __FUNCTION__, __LINE__, "Error while making selection in dataspace.");
		}
	}

	if (H5Dread(chromosomes_cache_it->second->haplotypes_id, H5T_NATIVE_UCHAR, memory_dataspace_id, file_dataspace_id, H5P_DEFAULT, haplotypes.get()) < 0) {
		throw HVCFReadException(__FILE__, __FUNCTION__, __LINE__, "Error while reading from dataset");
	}

	file_dataspace_id.close();
	memory_dataspace_id.close();

	end = std::chrono::system_clock::now();
	elapsed_seconds = end - start;
	cout << "Retrieved haplotypes in " << elapsed_seconds.count() << " seconds" << endl;

	start = std::chrono::system_clock::now();

	if (H5Tconvert(H5T_NATIVE_UCHAR, H5T_NATIVE_DOUBLE, n_variants * n_haplotypes, haplotypes.get(), nullptr, H5P_DEFAULT) < 0) {
		throw HVCFReadException(__FILE__, __FUNCTION__, __LINE__, "Error while converting datatypes.");
	}

	end = std::chrono::system_clock::now();
	elapsed_seconds = end - start;
	cout << "Type-casted haplotypes in " << elapsed_seconds.count() << " seconds" << endl;

	start = std::chrono::system_clock::now();
	Mat<double> S(haplotypes.get(), n_haplotypes, n_variants, false, false); // this call doesn't copy matrix.
	Row<double> L(S.colptr(lead_variant_local_offset), n_haplotypes);
	Row<double> J(n_haplotypes, fill::ones);

	Mat<double> LC1(J * L.t());
	Mat<double> LC2(n_haplotypes - LC1);
	Mat<double> SC1(J * S);
	Mat<double> SC2(n_haplotypes - SC1);
	Mat<double> M1(LC1 * SC1);

	Mat<double> R((n_haplotypes * L * S - M1) / sqrt(M1 % (LC2 * SC2)));

//	cout << "R:" << endl;
//	R.raw_print();

//	cout << "Rsq:" << endl;
//	cout << square(R) << endl;

	hsize_t file_offset1_1D[1]{static_cast<hsize_t>(lead_variant_offset)};
	hsize_t counts1_1D[1]{1};
	hsize_t file_offset2_1D[1]{static_cast<hsize_t>(start_position_offset)};
	hsize_t counts2_1D[1]{static_cast<hsize_t>(end_position_offset - start_position_offset + 1)};
	hsize_t mem_dims_1D[1]{n_variants};

	variants_entry_type variants_buffer[n_variants];

	end = std::chrono::system_clock::now();
	elapsed_seconds = end - start;
	cout << "Computed R in " << elapsed_seconds.count() << " seconds" << endl;

	start = std::chrono::system_clock::now();

	if ((file_dataspace_id = H5Dget_space(chromosomes_cache_it->second->variants_id)) < 0) {
		throw HVCFReadException(__FILE__, __FUNCTION__, __LINE__, "Error while getting dataspace.");
	}

	if ((memory_dataspace_id = H5Screate_simple(1, mem_dims_1D, nullptr)) < 0) {
		throw HVCFReadException(__FILE__, __FUNCTION__, __LINE__, "Error while creating memory dataspace.");
	}

	if (H5Sselect_hyperslab(file_dataspace_id, H5S_SELECT_SET, file_offset1_1D, NULL, counts1_1D, NULL) < 0) {
		throw HVCFReadException(__FILE__, __FUNCTION__, __LINE__, "Error while making selection in dataspace.");
	}

	if (H5Sselect_hyperslab(file_dataspace_id, H5S_SELECT_OR, file_offset2_1D, NULL, counts2_1D, NULL) < 0) {
		throw HVCFReadException(__FILE__, __FUNCTION__, __LINE__, "Error while making selection in dataspace.");
	}

	if (H5Dread(chromosomes_cache_it->second->variants_id, variants_entry_memory_datatype_id, memory_dataspace_id, file_dataspace_id, H5P_DEFAULT, variants_buffer) < 0) {
		throw HVCFReadException(__FILE__, __FUNCTION__, __LINE__, "Error while reading from dataset");
	}

	end = std::chrono::system_clock::now();
	elapsed_seconds = end - start;
	cout << "Read variant names and positions in " << elapsed_seconds.count() << " seconds" << endl;

	start = std::chrono::system_clock::now();

	for (unsigned int i = 0; i < lead_variant_local_offset; ++i) {
		result.emplace_back(
				variants_buffer[lead_variant_local_offset].name, variants_buffer[lead_variant_local_offset].position,
				variants_buffer[i].name, variants_buffer[i].position,
				R.at(0, i), pow(R.at(0, i), 2.0));
	}

	for (unsigned int i = lead_variant_local_offset + 1; i < n_variants; ++i) {
		result.emplace_back(
				variants_buffer[lead_variant_local_offset].name, variants_buffer[lead_variant_local_offset].position,
				variants_buffer[i].name, variants_buffer[i].position,
				R.at(0, i), pow(R.at(0, i), 2.0));
	}

	end = std::chrono::system_clock::now();
	elapsed_seconds = end - start;
	cout << "Formatted output in " << elapsed_seconds.count() << " seconds (" << n_variants << ")" << endl;

	if (H5Dvlen_reclaim(variants_entry_memory_datatype_id, memory_dataspace_id, H5P_DEFAULT, variants_buffer) < 0) {
		throw HVCFReadException(__FILE__, __FUNCTION__, __LINE__, "Error while reclaiming HDF5 memory.");
	}

	return;
}




}
