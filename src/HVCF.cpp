#include "include/HVCF.h"

namespace sph_umich_edu {

constexpr char HVCF::SAMPLES_GROUP[];
constexpr char HVCF::SAMPLES_ALL_DATASET[];
constexpr char HVCF::CHROMOSOMES_GROUP[];
constexpr char HVCF::VARIANTS_DATASET[];
constexpr char HVCF::HAPLOTYPES_DATASET[];
constexpr char HVCF::VARIABLE_LENGTH_STRING_TYPE[];
constexpr char HVCF::VARIANTS_ENTRY_TYPE[];
constexpr char HVCF::STRING_INDEX_ENTRY_TYPE[];
constexpr char HVCF::INTERVAL_INDEX_ENTRY_TYPE[];
constexpr char HVCF::HASH_INDEX_ENTRY_TYPE[];
constexpr char HVCF::ULL_INDEX_ENTRY_TYPE[];
constexpr char HVCF::VARIANT_NAMES_INDEX_GROUP[];
constexpr char HVCF::VARIANT_INTERVALS_INDEX_GROUP[];
constexpr char HVCF::VARIANT_INTERVALS_INDEX[];
constexpr char HVCF::VARIANT_HASH_INDEX[];
constexpr char HVCF::INDEX_BUCKETS[];


HVCF::HVCF() : HVCF(HVCFConfiguration()) {

}

HVCF::HVCF(const HVCFConfiguration& configuration) {
//	Disables automatic HDF5 error stack printing to stderr when function call returns negative value.
//	H5Eset_auto(H5E_DEFAULT, nullptr, nullptr);

	N_HASH_BUCKETS = configuration.n_hash_buckets;
	MAX_VARIANTS_IN_INTERVAL_BUCKET = configuration.max_variants_in_interval_bucket;
	VARIANTS_CHUNK_SIZE = configuration.variants_chunk_size;
	SAMPLES_CHUNK_SIZE = configuration.samples_chunk_size;
	COMPRESSION = configuration.compression;
	COMPRESSION_LEVEL = configuration.compression_level;

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

	if (H5Tinsert(memory_datatype_id, "position", HOFFSET(variants_entry_type, position), H5T_NATIVE_ULLONG) < 0) {
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

hid_t HVCF::create_strings_1D_dataset(const string& name, hid_t group_id, hsize_t chunk_size) throw (HVCFWriteException) {
	HDF5DataspaceIdentifier dataspace_id;
	HDF5DatasetIdentifier dataset_id;
	HDF5PropertyIdentifier dataset_property_id;
	HDF5DatatypeIdentifier native_string_datatype_id;

	hsize_t initial_dims[1]{0};
	hsize_t maximum_dims[1]{H5S_UNLIMITED};
	hsize_t chunk_dims[1]{chunk_size};

	if (name.length() == 0u) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Empty dataset name.");
	}

	if ((native_string_datatype_id = H5Topen(file_id, VARIABLE_LENGTH_STRING_TYPE, H5P_DEFAULT)) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while opening datatype.");
	}

	if ((dataspace_id = H5Screate_simple(1, initial_dims, maximum_dims)) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while creating dataspace.");
	}

	if ((dataset_property_id = H5Pcreate(H5P_DATASET_CREATE)) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while creating dataset property.");
	}

	if ((H5Pset_chunk(dataset_property_id, 1, chunk_dims) < 0) || (H5Pset_deflate(dataset_property_id, 9) < 0)) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while setting dataset properties.");
	}

	if ((dataset_id = H5Dcreate(group_id, name.c_str(), native_string_datatype_id, dataspace_id, H5P_DEFAULT, dataset_property_id, H5P_DEFAULT)) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while creating dataset.");
	}

	return dataset_id.release();
}

hid_t HVCF::create_hsize_1D_dataset(const string&name, hid_t group_id, hsize_t chunk_size) throw (HVCFWriteException) {
	HDF5DataspaceIdentifier dataspace_id;
	HDF5DatasetIdentifier dataset_id;
	HDF5PropertyIdentifier dataset_property_id;

	hsize_t initial_dims[1]{0};
	hsize_t maximum_dims[1]{H5S_UNLIMITED};
	hsize_t chunk_dims[1]{chunk_size};

	if (name.length() == 0u) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Empty dataset name.");
	}

	if ((dataspace_id = H5Screate_simple(1, initial_dims, maximum_dims)) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while creating dataspace.");
	}

	if ((dataset_property_id = H5Pcreate(H5P_DATASET_CREATE)) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while creating dataset property.");
	}

	if ((H5Pset_chunk(dataset_property_id, 1, chunk_dims) < 0) || (H5Pset_deflate(dataset_property_id, 9) < 0)) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while setting dataset properties.");
	}

	if ((dataset_id = H5Dcreate(group_id, name.c_str(), H5T_NATIVE_HSIZE, dataspace_id, H5P_DEFAULT, dataset_property_id, H5P_DEFAULT)) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while creating dataset.");
	}

	return dataset_id.release();
}

void HVCF::initialize_ull_index_buckets(hid_t chromosome_group_id, const char* index_group_name) throw (HVCFWriteException) {
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

	if ((index_group_id = H5Gopen(chromosome_group_id, index_group_name, H5P_DEFAULT)) < 0) {
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

void HVCF::initialize_string_index_buckets(hid_t chromosome_group_id, const char* index_group_name) throw (HVCFWriteException) {
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

	if ((index_group_id = H5Gopen(chromosome_group_id, index_group_name, H5P_DEFAULT)) < 0) {
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

void HVCF::write_intervals_index_bucket(hid_t chromosome_group_id, const vector<hsize_t>& offsets, interval_index_entry_type& interval_index_entry) throw (HVCFWriteException) {
	if (offsets.size() == 0u) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Empty list of offsets");
	}

	HDF5DatatypeIdentifier ull_index_entry_type_id;
	HDF5DatatypeIdentifier variants_entry_memory_datatype_id;
	HDF5DatatypeIdentifier ull_index_entry_memory_datatype_id;

	HDF5DatasetIdentifier dataset_id;
	HDF5DataspaceIdentifier file_dataspace_id;
	HDF5DataspaceIdentifier memory_dataspace_id;
	HDF5PropertyIdentifier dataset_property_id;
	HDF5GroupIdentifier index_group_id;

	hsize_t file_dims[1]{offsets.size()};
	hsize_t mem_dims[1]{offsets.size()};
	hsize_t offset[1]{0};

	variants_entry_type variants_buffer[offsets.size()];
	ull_index_entry_type index_entries_buffer[offsets.size()];

	// BEGIN: open comitted datatype.
	if ((ull_index_entry_type_id = H5Topen(file_id, ULL_INDEX_ENTRY_TYPE, H5P_DEFAULT)) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while opening datatype.");
	}
	// END: open comitted datatype.

	// BEGIN: create memory datatypes.
	try {
		variants_entry_memory_datatype_id = create_variants_entry_memory_datatype();
		ull_index_entry_memory_datatype_id = create_ull_index_entry_memory_datatype();
	} catch (HVCFCreateException &e) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while creating memory datatypes.");
	}
	// END: create memory datatypes.

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

	file_dataspace_id.close();
	dataset_id.close();
	// END: read variants.

	// BEGIN: sort bucket's entries.
	for (unsigned int i = 0; i < offsets.size(); ++i) {
		index_entries_buffer[i].ull_value = variants_buffer[i].position;
		index_entries_buffer[i].offset = offsets[i];
	}

	sort(index_entries_buffer, index_entries_buffer + offsets.size(),
			[] (const ull_index_entry_type& f, const ull_index_entry_type& s) -> bool {
				return (f.ull_value < s.ull_value);
	});
	// END: sort bucket's entries.

	// BEGIN: write bucket.
	if ((index_group_id = H5Gopen(chromosome_group_id, VARIANT_INTERVALS_INDEX_GROUP, H5P_DEFAULT)) < 0) {
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

	offset[0] = file_dims[0];
	file_dims[0] += mem_dims[0];

	if (H5Dset_extent(dataset_id, file_dims) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while setting dataset dimensions.");
	}

	if ((file_dataspace_id = H5Dget_space(dataset_id)) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while getting dataspace.");
	}

	if (H5Sselect_hyperslab(file_dataspace_id, H5S_SELECT_SET, offset, nullptr, mem_dims, nullptr) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while making selection in dataspace.");
	}

	if (H5Dwrite(dataset_id, ull_index_entry_memory_datatype_id, memory_dataspace_id, file_dataspace_id, H5P_DEFAULT, index_entries_buffer) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while writing to dataset.");
	}
	// END: write bucket.

	interval_index_entry.ull_value_1 = index_entries_buffer[0].ull_value;
	interval_index_entry.ull_value_2 = index_entries_buffer[offsets.size() - 1].ull_value;
	interval_index_entry.offset_1 = index_entries_buffer[0].offset;
	interval_index_entry.offset_2 = index_entries_buffer[offsets.size() - 1].offset;
	interval_index_entry.bucket_offset = offset[0];
	interval_index_entry.bucket_size = mem_dims[0];
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

	if ((index_group_id = H5Gopen(chromosome_group_id, VARIANT_INTERVALS_INDEX_GROUP, H5P_DEFAULT)) < 0) {
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

	if ((dataset_id = H5Dcreate(index_group_id, VARIANT_INTERVALS_INDEX, interval_index_entry_type_id, file_dataspace_id, H5P_DEFAULT, dataset_property_id, H5P_DEFAULT)) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while creating dataset.");
	}

	if (H5Dwrite(dataset_id, interval_index_entry_memory_datatype_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, interval_index_entries) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while writing to dataset.");
	}
}

void HVCF::write_names_index_bucket(hid_t chromosome_group_id, const vector<hsize_t>& offsets, hash_index_entry_type& hash_index_entry) throw (HVCFWriteException) {
	if (offsets.size() == 0u) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Empty list of offsets");
	}

	HDF5DatatypeIdentifier string_index_entry_type_id;
	HDF5DatatypeIdentifier variants_entry_memory_datatype_id;
	HDF5DatatypeIdentifier string_index_entry_memory_datatype_id;

	HDF5DatasetIdentifier dataset_id;
	HDF5DataspaceIdentifier file_dataspace_id;
	HDF5DataspaceIdentifier memory_dataspace_id;
	HDF5PropertyIdentifier dataset_property_id;
	HDF5GroupIdentifier index_group_id;

	hsize_t file_dims[1]{offsets.size()};
	hsize_t mem_dims[1]{offsets.size()};
	hsize_t offset[1]{0};

	variants_entry_type variants_buffer[offsets.size()];
	string_index_entry_type index_buffer[offsets.size()];

	// BEGIN: open comitted datatype.
	if ((string_index_entry_type_id = H5Topen(file_id, STRING_INDEX_ENTRY_TYPE, H5P_DEFAULT)) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while opening datatype.");
	}
	// END: open comitted datatype.

	// BEGIN: create memory datatypes.
	try {
		variants_entry_memory_datatype_id = create_variants_entry_memory_datatype();
		string_index_entry_memory_datatype_id = create_string_index_entry_memory_datatype();
	} catch (HVCFCreateException &e) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while creating memory datatypes.");
	}
	// END: create memory datatypes.

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

	file_dataspace_id.close();
	dataset_id.close();
	// END: read variants.

	// BEGIN: sort bucket's entries.
	for (unsigned int i = 0; i < offsets.size(); ++i) {
		index_buffer[i].string_value = variants_buffer[i].name;
		index_buffer[i].offset = offsets[i];
	}

	sort(index_buffer, index_buffer + offsets.size(),
			[] (const string_index_entry_type& f, const string_index_entry_type& s) -> bool {
				return (strcmp(f.string_value, s.string_value) < 0);
	});
	// END: sort bucket's entries.

	// BEGIN: write bucket.
	if ((index_group_id = H5Gopen(chromosome_group_id, VARIANT_NAMES_INDEX_GROUP, H5P_DEFAULT)) < 0) {
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

	offset[0] = file_dims[0];
	file_dims[0] += mem_dims[0];

	if (H5Dset_extent(dataset_id, file_dims) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while setting dataset dimensions.");
	}

	if ((file_dataspace_id = H5Dget_space(dataset_id)) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while getting dataspace.");
	}

	if (H5Sselect_hyperslab(file_dataspace_id, H5S_SELECT_SET, offset, nullptr, mem_dims, nullptr) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while making selection in dataspace.");
	}

	if (H5Dwrite(dataset_id, string_index_entry_memory_datatype_id, memory_dataspace_id, file_dataspace_id, H5P_DEFAULT, index_buffer) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while writing to dataset.");
	}

	if (H5Dvlen_reclaim(variants_entry_memory_datatype_id, memory_dataspace_id, H5P_DEFAULT, variants_buffer) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while reclaiming HDF5 memory.");
	}
	// END: write bucket.

	hash_index_entry.bucket_offset = offset[0];
	hash_index_entry.bucket_size = mem_dims[0];
}

void HVCF::write_names_index(hid_t chromosome_group_id, const hash_index_entry_type* hash_index_entries, unsigned int n_hash_index_entries) throw (HVCFWriteException) {
	if (n_hash_index_entries == 0u) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Empty list of offsets");
	}

	HDF5DatatypeIdentifier hash_index_entry_type_id;
	HDF5DatatypeIdentifier hash_index_entry_memory_datatype_id;

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

	// BEGIN: create memory datatype for index.
	try {
		hash_index_entry_memory_datatype_id = create_hash_index_entry_memory_datatype();
	} catch (HVCFCreateException &e) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while creating memory datatypes.");
	}
	// END: create memory datatype index.

	if ((index_group_id = H5Gopen(chromosome_group_id, VARIANT_NAMES_INDEX_GROUP, H5P_DEFAULT)) < 0) {
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

	if ((dataset_id = H5Dcreate(index_group_id, VARIANT_HASH_INDEX, hash_index_entry_type_id, file_dataspace_id, H5P_DEFAULT, dataset_property_id, H5P_DEFAULT)) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while creating dataset.");
	}

	if (H5Dwrite(dataset_id, hash_index_entry_memory_datatype_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, hash_index_entries) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while writing to dataset.");
	}
}

void HVCF::create_indices(hid_t chromosome_group_id) throw (HVCFWriteException) {
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
	if ((names_index_group_id = H5Gcreate(chromosome_group_id, VARIANT_NAMES_INDEX_GROUP, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT)) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while creating group.");
	}

	if ((intervals_index_group_id = H5Gcreate(chromosome_group_id, VARIANT_INTERVALS_INDEX_GROUP, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT)) < 0) {
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
		throw HVCFReadException(__FILE__, __FUNCTION__, __LINE__, "Error while getting dataspace dimensions.");
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
			hash_value = hash_function(buffer[i].name) % N_HASH_BUCKETS;
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
	initialize_ull_index_buckets(chromosome_group_id, VARIANT_INTERVALS_INDEX_GROUP);
	initialize_string_index_buckets(chromosome_group_id, VARIANT_NAMES_INDEX_GROUP);

	hash_index_entry_type hash_index_keys[N_HASH_BUCKETS];
	for (unsigned int i = 0u; i < N_HASH_BUCKETS; ++i) {
		hash_index_keys[i].bucket_offset = 0;
		hash_index_keys[i].bucket_size = 0;
	}

	for (auto&& bucket : names_index_buckets) {
		write_names_index_bucket(chromosome_group_id, bucket.second, hash_index_keys[bucket.first]);
	}
	write_names_index(chromosome_group_id, hash_index_keys, N_HASH_BUCKETS);

	interval_index_entry_type interval_index_keys[intervals.size()];
	for (unsigned int i = 0; i < intervals.size(); ++i) {
		write_intervals_index_bucket(chromosome_group_id, intervals[i], interval_index_keys[i]);
	}
	write_intervals_index(chromosome_group_id, interval_index_keys, intervals.size());
	// END: write indices on disk.
}

void HVCF::create(const string& name) throw (HVCFWriteException) {
	HDF5DatatypeIdentifier datatype_id;
	HDF5DatatypeIdentifier native_string_datatype_id;
	size_t native_string_datatype_size = 0;

	this->name = name;

	if ((file_id = H5Fcreate(this->name.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT)) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while creating file.");
	}

	// BEGIN: create and commit variable length string datatype.
	if (((datatype_id = H5Tcopy(H5T_C_S1)) < 0) || (H5Tset_size(datatype_id, H5T_VARIABLE) < 0)) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while creating variable length string datatype.");
	}

	if (H5Tcommit(file_id, VARIABLE_LENGTH_STRING_TYPE, datatype_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while committing variable length string datatype.");
	}

	if ((native_string_datatype_id = H5Topen(file_id, VARIABLE_LENGTH_STRING_TYPE, H5P_DEFAULT)) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while opening datatype.");
	}

	if ((native_string_datatype_size = H5Tget_size(native_string_datatype_id)) == 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while obtaining datatype size.");
	}
	datatype_id.close();
	// END: create and commit variable length string datatype.

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
	if ((datatype_id = H5Tcreate(H5T_COMPOUND, native_string_datatype_size + sizeof(unsigned long long int))) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while creating compound datatype.");
	}

	if (H5Tinsert(datatype_id, "name", 0, native_string_datatype_id) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while adding new member to compound datatype.");
	}

	if (H5Tinsert(datatype_id, "position", native_string_datatype_size, H5T_NATIVE_ULLONG) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while adding new member to compound datatype.");
	}

	if	(H5Tcommit(file_id, VARIANTS_ENTRY_TYPE, datatype_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT) < 0) {
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
}

void HVCF::set_samples(const vector<string>& samples) throw (HVCFWriteException) {
	HDF5DatasetIdentifier samples_all_dataset_id;
	HDF5DataspaceIdentifier file_dataspace_id;
	HDF5DataspaceIdentifier memory_dataspace_id;
	HDF5DatatypeIdentifier native_string_datatype_id;

	hsize_t mem_dims[1]{samples.size()};
	hsize_t file_dims[1]{0};
	hsize_t file_offset[1]{0};
	const char* buffer[samples.size()];

	for (unsigned int i = 0u; i < samples.size(); ++i) {
		buffer[i] = samples.at(i).c_str();
	}

	if ((native_string_datatype_id = H5Topen(file_id, VARIABLE_LENGTH_STRING_TYPE, H5P_DEFAULT)) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while opening datatype.");
	}

	if ((samples_all_dataset_id = create_strings_1D_dataset(SAMPLES_ALL_DATASET, samples_group_id, 1000)) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while creating dataset.");
	}

	if ((file_dataspace_id = H5Dget_space(samples_all_dataset_id)) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while getting dataspace.");
	}

	if (H5Sget_simple_extent_dims(file_dataspace_id, file_dims, nullptr) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while getting dataspace dimensions.");
	}

	file_dataspace_id.close();

	file_offset[0] = file_dims[0];
	file_dims[0] += mem_dims[0];

	if (H5Dset_extent(samples_all_dataset_id, file_dims) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while setting dataset dimensions.");
	}

	if ((file_dataspace_id = H5Dget_space(samples_all_dataset_id)) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while getting dataspace.");
	}

	if ((memory_dataspace_id = H5Screate_simple(1, mem_dims, nullptr)) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while creating memory dataspace.");
	}

	if (H5Sselect_hyperslab(file_dataspace_id, H5S_SELECT_SET, file_offset, nullptr, mem_dims, nullptr) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while making selection in dataspace.");
	}

	if (H5Dwrite(samples_all_dataset_id, native_string_datatype_id, memory_dataspace_id, file_dataspace_id, H5P_DEFAULT, buffer) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while writing to dataset.");
	}
}

void HVCF::set_population(const string& name, const std::vector<string>& samples) throw (HVCFWriteException) {
	HDF5DatasetIdentifier population_dataset_id;
	HDF5DataspaceIdentifier file_dataspace_id;
	HDF5DataspaceIdentifier memory_dataspace_id;

	map<string, unsigned int> all_samples;
	auto all_samples_it = all_samples.end();
	hsize_t buffer[samples.size()];

	hsize_t mem_dims[1]{samples.size()};
	hsize_t file_dims[1]{samples.size()};
	hsize_t file_offset[1]{0};

	for (auto&& sample : get_samples()) {
		all_samples.emplace(std::move(sample), all_samples.size());
	}

	for (unsigned int i = 0u; i < samples.size(); ++i) {
		all_samples_it = all_samples.find(samples[i]);
		if (all_samples_it == all_samples.end()) {
			throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Sample not found.");
		}
		buffer[i] = static_cast<hsize_t>(all_samples_it->second);
	}

	std::sort(buffer, buffer + samples.size(), std::less_equal<unsigned int>());

	population_dataset_id = create_hsize_1D_dataset(name, samples_group_id, 100);

	if (H5Dset_extent(population_dataset_id, file_dims) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while setting dataset dimensions.");
	}

	if ((file_dataspace_id = H5Dget_space(population_dataset_id)) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while getting dataspace.");
	}

	if ((memory_dataspace_id = H5Screate_simple(1, mem_dims, nullptr)) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while creating memory dataspace.");
	}

	if (H5Sselect_hyperslab(file_dataspace_id, H5S_SELECT_SET, file_offset, nullptr, mem_dims, nullptr) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while making selection in dataspace.");
	}

	if (H5Dwrite(population_dataset_id, H5T_NATIVE_HSIZE, memory_dataspace_id, file_dataspace_id, H5P_DEFAULT, buffer) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while writing to dataset.");
	}
}

void HVCF::write_variant(const Variant& variant) throw (HVCFWriteException) {
	const string& chromosome = variant.get_chrom().get_value();
	auto chromosomes_it = chromosomes.end();
	auto buffers_it = write_buffers.end();

	if (variant.get_alt().get_values().size() != 1) { // Support only bi-allelic (for computing LD it is file, but must be extended).
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
		write_haplotypes(chromosomes_it->second->get(), buffers_it->second->get_haplotypes_buffer(), buffers_it->second->get_n_variants(), buffers_it->second->get_n_haplotypes());
		write_variants(chromosomes_it->second->get(), buffers_it->second->get_variants_buffer(), buffers_it->second->get_n_variants());
		buffers_it->second->reset();
	}

	buffers_it->second->add_variant(variant);
}

void HVCF::flush_write_buffer() throw (HVCFWriteException) {
	auto buffers_it = write_buffers.end();

	for (auto&& entry : chromosomes) {
		buffers_it = write_buffers.find(entry.first);
		if (!buffers_it->second->is_empty()) {
			write_haplotypes(entry.second->get(), buffers_it->second->get_haplotypes_buffer(), buffers_it->second->get_n_variants(), buffers_it->second->get_n_haplotypes());
			write_variants(entry.second->get(), buffers_it->second->get_variants_buffer(), buffers_it->second->get_n_variants());
			buffers_it->second->reset();
		}
	}
}

void HVCF::create_indices() throw (HVCFWriteException) {
	for (auto&& chromosome : chromosomes) {
		create_indices(chromosome.second->get());
	}
}

hsize_t HVCF::get_n_samples() throw (HVCFReadException) {
	HDF5DatasetIdentifier samples_all_dataset_id;
	HDF5DataspaceIdentifier file_dataspace_id;
	hsize_t file_dims[1]{0};

	if ((samples_all_dataset_id = H5Dopen(samples_group_id, SAMPLES_ALL_DATASET, H5P_DEFAULT)) < 0) {
		throw HVCFOpenException(__FILE__, __FUNCTION__, __LINE__, "Error while opening dataset.");
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
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while opening datatype.");
	}

	if ((samples_all_dataset_id = H5Dopen(samples_group_id, SAMPLES_ALL_DATASET, H5P_DEFAULT)) < 0) {
		throw HVCFOpenException(__FILE__, __FUNCTION__, __LINE__, "Error while opening dataset.");
	}

	if ((file_dataspace_id = H5Dget_space(samples_all_dataset_id)) < 0) {
		throw HVCFReadException(__FILE__, __FUNCTION__, __LINE__, "Error while getting dataspace.");
	}

	if (H5Sget_simple_extent_dims(file_dataspace_id, file_dims, nullptr) < 0) {
		throw HVCFReadException(__FILE__, __FUNCTION__, __LINE__, "Error while getting dataspace dimensions.");
	}

	file_dataspace_id.close();

	char* buffer[file_dims[0]];
	vector<string> samples;

	if (H5Dread(samples_all_dataset_id, native_string_datatype_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, buffer) < 0) {
		throw HVCFReadException(__FILE__, __FUNCTION__, __LINE__, "Error while reading from dataset");
	}

	for (unsigned int i = 0; i < file_dims[0]; ++i) {
		samples.emplace_back(buffer[i]);
	}

	if ((file_dataspace_id = H5Dget_space(samples_all_dataset_id)) < 0) {
		throw HVCFReadException(__FILE__, __FUNCTION__, __LINE__, "Error while getting dataspace.");
	}

	if (H5Dvlen_reclaim(native_string_datatype_id, file_dataspace_id, H5P_DEFAULT, buffer) < 0) {
		throw HVCFReadException(__FILE__, __FUNCTION__, __LINE__, "Error while reclaiming HDF5 memory.");
	}

	return samples;
}

vector<string> HVCF::get_population(const string& name) throw (HVCFReadException) {
	HDF5DatasetIdentifier samples_all_dataset_id;
	HDF5DatasetIdentifier dataset_id;
	HDF5DataspaceIdentifier file_dataspace_id;
	HDF5DataspaceIdentifier mem_dataspace_id;
	HDF5DatatypeIdentifier native_string_datatype_id;

	vector<string> samples;
	hsize_t file_dims[1]{0};
	hsize_t mem_dims[1]{0};

	if ((native_string_datatype_id = H5Topen(file_id, VARIABLE_LENGTH_STRING_TYPE, H5P_DEFAULT)) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while opening datatype.");
	}

	if ((dataset_id = H5Dopen(samples_group_id, name.c_str(), H5P_DEFAULT)) < 0) {
		throw HVCFReadException(__FILE__, __FUNCTION__, __LINE__, "Error while opening dataset.");
	}

	if ((file_dataspace_id = H5Dget_space(dataset_id)) < 0) {
		throw HVCFReadException(__FILE__, __FUNCTION__, __LINE__, "Error while getting dataspace.");
	}

	if (H5Sget_simple_extent_dims(file_dataspace_id, file_dims, nullptr) < 0) {
		throw HVCFReadException(__FILE__, __FUNCTION__, __LINE__, "Error while getting dataspace extent.");
	}

	file_dataspace_id.close();

	mem_dims[0] = file_dims[0];

	hsize_t index_buffer[file_dims[0]];
	char* samples_buffer[file_dims[0]];

	if (H5Dread(dataset_id, H5T_NATIVE_HSIZE, H5S_ALL, H5S_ALL, H5P_DEFAULT, index_buffer) < 0) {
		throw HVCFReadException(__FILE__, __FUNCTION__, __LINE__, "Error while reading from dataset.");
	}

	dataset_id.close();

	if ((samples_all_dataset_id = H5Dopen(samples_group_id, SAMPLES_ALL_DATASET, H5P_DEFAULT)) < 0) {
		throw HVCFOpenException(__FILE__, __FUNCTION__, __LINE__, "Error while opening dataset.");
	}

	if ((file_dataspace_id = H5Dget_space(samples_all_dataset_id)) < 0) {
		throw HVCFReadException(__FILE__, __FUNCTION__, __LINE__, "Error while getting dataspace.");
	}

	if (H5Sselect_elements(file_dataspace_id, H5S_SELECT_SET, file_dims[0], index_buffer) < 0) {
		throw HVCFReadException(__FILE__, __FUNCTION__, __LINE__, "Error while making selection in dataspace.");
	}

	if ((mem_dataspace_id = H5Screate_simple(1, mem_dims, nullptr)) < 0) {
		throw HVCFReadException(__FILE__, __FUNCTION__, __LINE__, "Error while creating dataspace.");
	}

	if (H5Dread(samples_all_dataset_id, native_string_datatype_id, mem_dataspace_id, file_dataspace_id, H5P_DEFAULT, samples_buffer) < 0) {
		throw HVCFReadException(__FILE__, __FUNCTION__, __LINE__, "Error while reading from dataset.");
	}

	for (unsigned int i = 0; i < file_dims[0]; ++i) {
		samples.emplace_back(samples_buffer[i]);
	}

	if (H5Dvlen_reclaim(native_string_datatype_id, mem_dataspace_id, H5P_DEFAULT, samples_buffer) < 0) {
		throw HVCFReadException(__FILE__, __FUNCTION__, __LINE__, "Error while reclaiming HDF5 memory.");
	}

	return samples;
}

void HVCF::open(const string& name) throw (HVCFOpenException) {
	this->name = name;

	H5AC_cache_config_t config;
	config.version = H5AC__CURR_CACHE_CONFIG_VERSION;

	HDF5PropertyIdentifier file_access_property_id;
	file_access_property_id = H5Pcreate(H5P_FILE_ACCESS);

	if (file_access_property_id < 0) {
		cout << "blia!" << endl;
	}
	if (H5Pget_mdc_config(file_access_property_id, &config) < 0) {
		cout << "blia2!" << endl;
	}


	config.incr_mode = H5C_cache_incr_mode::H5C_incr__off;
	config.decr_mode = H5C_cache_decr_mode::H5C_decr__off;
	config.set_initial_size = true;
	config.initial_size = 32 * 1024 * 1024;


	if (H5Pset_mdc_config(file_access_property_id, &config) < 0) {
		cout << "blia!" << endl;
	}

	if (H5Pset_cache(file_access_property_id, 0, 1000, 16 * 1024 * 1024, 0.75) < 0) {
		cout << "blia!" << endl;
	}

	if ((file_id = H5Fopen(this->name.c_str(), H5F_ACC_RDONLY, file_access_property_id)) < 0) {
		throw HVCFOpenException(__FILE__, __FUNCTION__, __LINE__, "Error while opening file.");
	}

	HDF5PropertyIdentifier file_access_property_id2;
	file_access_property_id2 = H5Fget_access_plist(file_id);
	H5AC_cache_config_t cache_config;
	cache_config.version = H5AC__CURR_CACHE_CONFIG_VERSION;
	if (H5Pget_mdc_config(file_access_property_id2, &cache_config) < 0) {
		cout << "blia3!" << endl;
	}
	cout << "Opened file configuration:" << endl;
	cout << "\tInitital size: " << cache_config.initial_size << endl;

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
}

void HVCF::close() throw (HVCFCloseException) {
	for (auto&& entry : chromosomes) {
		entry.second->close();
	}

	chromosomes_group_id.close();
	samples_group_id.close();
	file_id.close();
	name.clear();
}

hsize_t HVCF::get_n_variants() throw (HVCFReadException) {
	hsize_t total = 0;

	for (auto&& entry : chromosomes) {
		total += get_n_variants(entry.first);
	}

	return total;
}

hsize_t HVCF::get_n_variants(const string& chromosome) throw (HVCFReadException) {
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
	auto chromosomes_it = chromosomes.find(chromosome);

	if (chromosomes_it == chromosomes.end()) {
		return -1;
	}

	HDF5GroupIdentifier index_group_id;
	HDF5DatasetIdentifier dataset_id;
	HDF5DataspaceIdentifier dataspace_id;
	HDF5DataspaceIdentifier memory_dataspace_id;
	HDF5DatatypeIdentifier interval_index_entry_memory_datatype_id;
	HDF5DatatypeIdentifier ull_index_key_memory_datatype_id;

	hsize_t file_dims[1]{0};

	try {
		ull_index_key_memory_datatype_id = create_ull_index_entry_memory_datatype();
		interval_index_entry_memory_datatype_id = create_interval_index_entry_memory_datatype();
	} catch (HVCFCreateException &e) {
		throw HVCFReadException(__FILE__, __FUNCTION__, __LINE__, "Error while creating memory datatypes.");
	}

	if ((index_group_id = H5Gopen(chromosomes_it->second->get(), VARIANT_INTERVALS_INDEX_GROUP, H5P_DEFAULT)) < 0) {
		throw HVCFReadException(__FILE__, __FUNCTION__, __LINE__, "Error while opening group.");
	}

	if ((dataset_id = H5Dopen(index_group_id, VARIANT_INTERVALS_INDEX, H5P_DEFAULT)) < 0) {
		throw HVCFReadException(__FILE__, __FUNCTION__, __LINE__, "Error while opening dataset.");
	}

	if ((dataspace_id = H5Dget_space(dataset_id)) < 0) {
		throw HVCFReadException(__FILE__, __FUNCTION__, __LINE__, "Error while getting dataspace.");
	}

	if (H5Sget_simple_extent_dims(dataspace_id, file_dims, NULL) < 0) {
		throw HVCFReadException(__FILE__, __FUNCTION__, __LINE__, "Error while getting dataspace dimensions.");
	}

	interval_index_entry_type interval_index_entries_buffer[file_dims[0]];

	if (H5Dread(dataset_id, interval_index_entry_memory_datatype_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, interval_index_entries_buffer) < 0) {
		throw HVCFReadException(__FILE__, __FUNCTION__, __LINE__, "Error while reading from dataset.");
	}

	dataspace_id.close();
	dataset_id.close();

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
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while creating memory dataspace.");
	}

	if ((dataset_id = H5Dopen(index_group_id, INDEX_BUCKETS, H5P_DEFAULT)) < 0) {
		throw HVCFReadException(__FILE__, __FUNCTION__, __LINE__, "Error while opening dataset.");
	}

	if ((dataspace_id = H5Dget_space(dataset_id)) < 0) {
		throw HVCFReadException(__FILE__, __FUNCTION__, __LINE__, "Error while getting dataspace.");
	}

	if (H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, offset, NULL, mem_dims, NULL) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while making selection in dataspace.");
	}

	if (H5Dread(dataset_id, ull_index_key_memory_datatype_id, memory_dataspace_id, dataspace_id, H5P_DEFAULT, ull_index_entries_buffer) < 0) {
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
	auto chromosomes_it = chromosomes.find(chromosome);

	if (chromosomes_it == chromosomes.end()) {
		return -1;
	}

	HDF5GroupIdentifier index_group_id;
	HDF5DatasetIdentifier dataset_id;
	HDF5DataspaceIdentifier dataspace_id;
	HDF5DataspaceIdentifier memory_dataspace_id;
	HDF5DatatypeIdentifier interval_index_entry_memory_datatype_id;
	HDF5DatatypeIdentifier ull_index_entry_memory_datatype_id;

	hsize_t file_dims[1]{0};

	try {
		interval_index_entry_memory_datatype_id = create_interval_index_entry_memory_datatype();
		ull_index_entry_memory_datatype_id = create_ull_index_entry_memory_datatype();
	} catch (HVCFCreateException &e) {
		throw HVCFReadException(__FILE__, __FUNCTION__, __LINE__, "Error while creating memory datatypes.");
	}

	if ((index_group_id = H5Gopen(chromosomes_it->second->get(), VARIANT_INTERVALS_INDEX_GROUP, H5P_DEFAULT)) < 0) {
		throw HVCFReadException(__FILE__, __FUNCTION__, __LINE__, "Error while opening group.");
	}

	if ((dataset_id = H5Dopen(index_group_id, VARIANT_INTERVALS_INDEX, H5P_DEFAULT)) < 0) {
		throw HVCFReadException(__FILE__, __FUNCTION__, __LINE__, "Error while opening dataset.");
	}

	if ((dataspace_id = H5Dget_space(dataset_id)) < 0) {
		throw HVCFReadException(__FILE__, __FUNCTION__, __LINE__, "Error while getting dataspace.");
	}

	if (H5Sget_simple_extent_dims(dataspace_id, file_dims, NULL) < 0) {
		throw HVCFReadException(__FILE__, __FUNCTION__, __LINE__, "Error while getting dataspace dimensions.");
	}

	interval_index_entry_type interval_index_entries_buffer[file_dims[0]];

	if (H5Dread(dataset_id, interval_index_entry_memory_datatype_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, interval_index_entries_buffer) < 0) {
		throw HVCFReadException(__FILE__, __FUNCTION__, __LINE__, "Error while reading from dataset.");
	}

	dataspace_id.close();
	dataset_id.close();

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
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while creating memory dataspace.");
	}

	if ((dataset_id = H5Dopen(index_group_id, INDEX_BUCKETS, H5P_DEFAULT)) < 0) {
		throw HVCFReadException(__FILE__, __FUNCTION__, __LINE__, "Error while opening dataset.");
	}

	if ((dataspace_id = H5Dget_space(dataset_id)) < 0) {
		throw HVCFReadException(__FILE__, __FUNCTION__, __LINE__, "Error while getting dataspace.");
	}

	if (H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, offset, NULL, mem_dims, NULL) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while making selection in dataspace.");
	}

	if (H5Dread(dataset_id, ull_index_entry_memory_datatype_id, memory_dataspace_id, dataspace_id, H5P_DEFAULT, ull_index_entries_buffer) < 0) {
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
	auto chromosomes_it = chromosomes.find(chromosome);

	if (chromosomes_it == chromosomes.end()) {
		return -1;
	}

	HDF5GroupIdentifier index_group_id;
	HDF5DatasetIdentifier dataset_id;
	HDF5DataspaceIdentifier dataspace_id;
	HDF5DataspaceIdentifier memory_dataspace_id;
	HDF5DatatypeIdentifier interval_index_entry_memory_datatype_id;
	HDF5DatatypeIdentifier ull_index_entry_memory_datatype_id;

	hsize_t file_dims[1]{0};

	try {
		ull_index_entry_memory_datatype_id = create_ull_index_entry_memory_datatype();
		interval_index_entry_memory_datatype_id = create_interval_index_entry_memory_datatype();
	} catch (HVCFCreateException &e) {
		throw HVCFReadException(__FILE__, __FUNCTION__, __LINE__, "Error while creating memory datatypes.");
	}

	if ((index_group_id = H5Gopen(chromosomes_it->second->get(), VARIANT_INTERVALS_INDEX_GROUP, H5P_DEFAULT)) < 0) {
		throw HVCFReadException(__FILE__, __FUNCTION__, __LINE__, "Error while opening group.");
	}

	if ((dataset_id = H5Dopen(index_group_id, VARIANT_INTERVALS_INDEX, H5P_DEFAULT)) < 0) {
		throw HVCFReadException(__FILE__, __FUNCTION__, __LINE__, "Error while opening dataset.");
	}

	if ((dataspace_id = H5Dget_space(dataset_id)) < 0) {
		throw HVCFReadException(__FILE__, __FUNCTION__, __LINE__, "Error while getting dataspace.");
	}

	if (H5Sget_simple_extent_dims(dataspace_id, file_dims, NULL) < 0) {
		throw HVCFReadException(__FILE__, __FUNCTION__, __LINE__, "Error while getting dataspace dimensions.");
	}

	interval_index_entry_type interval_index_entries_buffer[file_dims[0]];

	if (H5Dread(dataset_id, interval_index_entry_memory_datatype_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, interval_index_entries_buffer) < 0) {
		throw HVCFReadException(__FILE__, __FUNCTION__, __LINE__, "Error while reading from dataset.");
	}

	dataspace_id.close();
	dataset_id.close();

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
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while creating memory dataspace.");
	}

	if ((dataset_id = H5Dopen(index_group_id, INDEX_BUCKETS, H5P_DEFAULT)) < 0) {
		throw HVCFReadException(__FILE__, __FUNCTION__, __LINE__, "Error while opening dataset.");
	}

	if ((dataspace_id = H5Dget_space(dataset_id)) < 0) {
		throw HVCFReadException(__FILE__, __FUNCTION__, __LINE__, "Error while getting dataspace.");
	}

	if (H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, offset, NULL, mem_dims, NULL) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while making selection in dataspace.");
	}

	if (H5Dread(dataset_id, ull_index_entry_memory_datatype_id, memory_dataspace_id, dataspace_id, H5P_DEFAULT, ull_index_entries_buffer) < 0) {
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
	auto chromosomes_it = chromosomes.find(chromosome);

	if (chromosomes_it == chromosomes.end()) {
		return -1;
	}

	HDF5GroupIdentifier index_group_id;
	HDF5DatasetIdentifier dataset_id;
	HDF5DataspaceIdentifier dataspace_id;
	HDF5DataspaceIdentifier memory_dataspace_id;
	HDF5DatatypeIdentifier hash_index_entry_memory_datatype_id;
	HDF5DatatypeIdentifier string_index_entry_memory_datatype_id;

	std::hash<string> hash_function;

	hsize_t hash_value = hash_function(name) % N_HASH_BUCKETS;

	hsize_t offset[1]{hash_value};
	hsize_t mem_dims[1]{1};

	hash_index_entry_type hash_index_entry;

	try {
		hash_index_entry_memory_datatype_id = create_hash_index_entry_memory_datatype();
		string_index_entry_memory_datatype_id = create_string_index_entry_memory_datatype();
	} catch (HVCFCreateException &e) {
		throw HVCFReadException(__FILE__, __FUNCTION__, __LINE__, "Error while creating memory datatypes.");
	}

	if ((index_group_id = H5Gopen(chromosomes_it->second->get(), VARIANT_NAMES_INDEX_GROUP, H5P_DEFAULT)) < 0) {
		throw HVCFReadException(__FILE__, __FUNCTION__, __LINE__, "Error while opening group.");
	}

	if ((dataset_id = H5Dopen(index_group_id, VARIANT_HASH_INDEX, H5P_DEFAULT)) < 0) {
		throw HVCFReadException(__FILE__, __FUNCTION__, __LINE__, "Error while opening dataset.");
	}

	if ((dataspace_id = H5Dget_space(dataset_id)) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while getting dataspace.");
	}

	if ((memory_dataspace_id = H5Screate_simple(1, mem_dims, nullptr)) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while creating memory dataspace.");
	}

	if (H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, offset, NULL, mem_dims, NULL) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while making selection in dataspace.");
	}

	if (H5Dread(dataset_id, hash_index_entry_memory_datatype_id, memory_dataspace_id, dataspace_id, H5P_DEFAULT, &hash_index_entry) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while reading from dataset");
	}

	dataspace_id.close();
	memory_dataspace_id.close();
	dataset_id.close();

	if ((hash_index_entry.bucket_offset == 0) && (hash_index_entry.bucket_size == 0)) {
		return -1;
	}

	offset[0] = hash_index_entry.bucket_offset;
	mem_dims[0] = hash_index_entry.bucket_size;

	string_index_entry_type string_index_entries_buffer[hash_index_entry.bucket_size];

	if ((dataset_id = H5Dopen(index_group_id, INDEX_BUCKETS, H5P_DEFAULT)) < 0) {
		throw HVCFReadException(__FILE__, __FUNCTION__, __LINE__, "Error while opening dataset.");
	}

	if ((dataspace_id = H5Dget_space(dataset_id)) < 0) {
		throw HVCFReadException(__FILE__, __FUNCTION__, __LINE__, "Error while getting dataspace.");
	}

	if ((memory_dataspace_id = H5Screate_simple(1, mem_dims, nullptr)) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while creating memory dataspace.");
	}

	if (H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, offset, NULL, mem_dims, NULL) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while making selection in dataspace.");
	}

	if (H5Dread(dataset_id, string_index_entry_memory_datatype_id, memory_dataspace_id, dataspace_id, H5P_DEFAULT, string_index_entries_buffer) < 0) {
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

void HVCF::compute_ld(const string& chromosome, unsigned long long int start_position, unsigned long long end_position, vector<variants_pair>& result) throw (HVCFReadException) {
	auto chromosome_it = chromosomes.find(chromosome);

	if (chromosome_it == chromosomes.end()) {
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

	HDF5DatasetIdentifier dataset_id;
	HDF5DataspaceIdentifier file_dataspace_id;
	HDF5DataspaceIdentifier memory_dataspace_id;

	hsize_t n_samples = get_n_samples();
	hsize_t n_haplotypes = 2 * n_samples;

	hsize_t n_variants = end_position_offset - start_position_offset + 1;

	hsize_t file_offset_2D[2]{static_cast<hsize_t>(start_position_offset), 0};
	hsize_t mem_dims_2D[2]{n_variants, n_haplotypes};

	unique_ptr<unsigned char[]> haplotypes = unique_ptr<unsigned char[]>(new unsigned char[n_variants * n_haplotypes]);

	if ((dataset_id = H5Dopen(chromosome_it->second->get(), HAPLOTYPES_DATASET, H5P_DEFAULT)) < 0) {
		throw HVCFReadException(__FILE__, __FUNCTION__, __LINE__, "Error while opening dataset.");
	}

	if ((file_dataspace_id = H5Dget_space(dataset_id)) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while getting dataspace.");
	}

	if ((memory_dataspace_id = H5Screate_simple(2, mem_dims_2D, nullptr)) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while creating memory dataspace.");
	}

	if (H5Sselect_hyperslab(file_dataspace_id, H5S_SELECT_SET, file_offset_2D, NULL, mem_dims_2D, NULL) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while making selection in dataspace.");
	}

	if (H5Dread(dataset_id, H5T_NATIVE_UCHAR, memory_dataspace_id, file_dataspace_id, H5P_DEFAULT, haplotypes.get()) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while reading from dataset");
	}

	dataset_id.close();
	file_dataspace_id.close();
	memory_dataspace_id.close();

	unique_ptr<double[]> dbl_haplotypes = unique_ptr<double[]>(new double[n_variants * n_haplotypes]);
	for (unsigned long long int i = 0u; i < n_variants * n_haplotypes; ++i) {
		dbl_haplotypes[i] = (double)haplotypes[i];
	}

	Mat<double> S(dbl_haplotypes.get(), n_haplotypes, n_variants, false, false); // this call doesn't copy matrix.
	Row<double> J(n_haplotypes, fill::ones);

	Mat<double> C1(J * S);
	Mat<double> C2(n_haplotypes - C1);
	Mat<double> M1(C1.t() * C1);

	Mat<double> R((n_haplotypes * S.t() * S - M1) / sqrt(M1 % (C2.t() * C2)));

//	cout << "R:" << endl;
//	R.raw_print();

//	cout << "Rsq:" << endl;
//	cout << square(R) << endl;

	HDF5DatatypeIdentifier variants_entry_memory_datatype_id;

	hsize_t file_offset_1D[1]{static_cast<hsize_t>(start_position_offset)};
	hsize_t mem_dims_1D[1]{n_variants};

	variants_entry_type variants_buffer[n_variants];

	try {
		variants_entry_memory_datatype_id = create_variants_entry_memory_datatype();
	} catch (HVCFCreateException &e) {
		throw HVCFReadException(__FILE__, __FUNCTION__, __LINE__, "Error while creating memory datatypes.");
	}

	if ((dataset_id = H5Dopen(chromosome_it->second->get(), VARIANTS_DATASET, H5P_DEFAULT)) < 0) {
		throw HVCFReadException(__FILE__, __FUNCTION__, __LINE__, "Error while opening dataset.");
	}

	if ((file_dataspace_id = H5Dget_space(dataset_id)) < 0) {
		throw HVCFReadException(__FILE__, __FUNCTION__, __LINE__, "Error while getting dataspace.");
	}

	if ((memory_dataspace_id = H5Screate_simple(1, mem_dims_1D, nullptr)) < 0) {
		throw HVCFReadException(__FILE__, __FUNCTION__, __LINE__, "Error while creating memory dataspace.");
	}

	if (H5Sselect_hyperslab(file_dataspace_id, H5S_SELECT_SET, file_offset_1D, NULL, mem_dims_1D, NULL) < 0) {
		throw HVCFReadException(__FILE__, __FUNCTION__, __LINE__, "Error while making selection in dataspace.");
	}

	if (H5Dread(dataset_id, variants_entry_memory_datatype_id, memory_dataspace_id, file_dataspace_id, H5P_DEFAULT, variants_buffer) < 0) {
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

void HVCF::compute_ld(const string& chromosome, const string& lead_variant_name, unsigned long long int start_position, unsigned long long end_position, vector<variants_pair>& result) throw (HVCFReadException) {
	auto chromosome_it = chromosomes.find(chromosome);

	if (chromosome_it == chromosomes.end()) {
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

	HDF5DatasetIdentifier dataset_id;
	HDF5DataspaceIdentifier file_dataspace_id;
	HDF5DataspaceIdentifier memory_dataspace_id;

	hsize_t n_samples = get_n_samples();
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

	unique_ptr<unsigned char[]> haplotypes = unique_ptr<unsigned char[]>(new unsigned char[n_variants * n_haplotypes]); // why smart pointer???

	if ((dataset_id = H5Dopen(chromosome_it->second->get(), HAPLOTYPES_DATASET, H5P_DEFAULT)) < 0) {
		throw HVCFReadException(__FILE__, __FUNCTION__, __LINE__, "Error while opening dataset.");
	}

	if ((file_dataspace_id = H5Dget_space(dataset_id)) < 0) {
		throw HVCFReadException(__FILE__, __FUNCTION__, __LINE__, "Error while getting dataspace.");
	}

	if ((memory_dataspace_id = H5Screate_simple(2, mem_dims_2D, nullptr)) < 0) {
		throw HVCFReadException(__FILE__, __FUNCTION__, __LINE__, "Error while creating memory dataspace.");
	}

	if (H5Sselect_hyperslab(file_dataspace_id, H5S_SELECT_SET, file_offset1_2D, NULL, counts1_2D, NULL) < 0) {
		throw HVCFReadException(__FILE__, __FUNCTION__, __LINE__, "Error while making selection in dataspace.");
	}

	if (H5Sselect_hyperslab(file_dataspace_id, H5S_SELECT_OR, file_offset2_2D, NULL, counts2_2D, NULL) < 0) {
		throw HVCFReadException(__FILE__, __FUNCTION__, __LINE__, "Error while making selection in dataspace.");
	}

	if (H5Dread(dataset_id, H5T_NATIVE_UCHAR, memory_dataspace_id, file_dataspace_id, H5P_DEFAULT, haplotypes.get()) < 0) {
		throw HVCFReadException(__FILE__, __FUNCTION__, __LINE__, "Error while reading from dataset");
	}

	dataset_id.close();
	file_dataspace_id.close();
	memory_dataspace_id.close();

	unique_ptr<double[]> dbl_haplotypes = unique_ptr<double[]>(new double[n_variants * n_haplotypes]);
	for (unsigned long long int i = 0u; i < n_variants * n_haplotypes; ++i) {
		dbl_haplotypes[i] = (double)haplotypes[i];
	}

	Mat<double> S(dbl_haplotypes.get(), n_haplotypes, n_variants, false, false); // this call doesn't copy matrix.
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

	HDF5DatatypeIdentifier variants_entry_memory_datatype_id;

	hsize_t file_offset1_1D[1]{static_cast<hsize_t>(lead_variant_offset)};
	hsize_t counts1_1D[1]{1};
	hsize_t file_offset2_1D[1]{static_cast<hsize_t>(start_position_offset)};
	hsize_t counts2_1D[1]{static_cast<hsize_t>(end_position_offset - start_position_offset + 1)};
	hsize_t mem_dims_1D[1]{n_variants};

	variants_entry_type variants_buffer[n_variants];

	try {
		variants_entry_memory_datatype_id = create_variants_entry_memory_datatype();
	} catch (HVCFCreateException &e) {
		throw HVCFReadException(__FILE__, __FUNCTION__, __LINE__, "Error while creating memory datatypes.");
	}

	if ((dataset_id = H5Dopen(chromosome_it->second->get(), VARIANTS_DATASET, H5P_DEFAULT)) < 0) {
		throw HVCFReadException(__FILE__, __FUNCTION__, __LINE__, "Error while opening dataset.");
	}

	if ((file_dataspace_id = H5Dget_space(dataset_id)) < 0) {
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

	if (H5Dread(dataset_id, variants_entry_memory_datatype_id, memory_dataspace_id, file_dataspace_id, H5P_DEFAULT, variants_buffer) < 0) {
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

unsigned int HVCF::get_n_opened_objects() const {
	if (file_id >= 0) {
		return H5Fget_obj_count(file_id, H5F_OBJ_ALL);
	}
	return 0u;
}

unsigned int HVCF::get_n_all_opened_objects() {
	return H5Fget_obj_count(H5F_OBJ_ALL, H5F_OBJ_ALL);
}

void HVCF::chunk_read_test(const string& chromosome, unsigned long long int start_position, unsigned long long end_position) throw (HVCFReadException) {
	auto chromosome_it = chromosomes.find(chromosome); //TODO: check if chromosome exists

	long long int start_position_offset = get_variant_offset_by_position_eq(chromosome, start_position); // TODO: put 0 if not exists
	long long int end_position_offset = get_variant_offset_by_position_eq(chromosome, end_position); // TODO: put file_dims[0] - 1 offset if not exists

//	cout << lead_variant_offset << " " << start_position_offset << " " << end_position_offset << endl;

	HDF5DatasetIdentifier dataset_id;
	HDF5DataspaceIdentifier file_dataspace_id;
	HDF5DataspaceIdentifier memory_dataspace_id;

	hsize_t n_samples = get_n_samples(); // TODO: compute this only once on file opening
	hsize_t n_haplotypes = n_samples + n_samples; // TODO: compute this only once on file opening
	long long int n_variants = end_position_offset - start_position_offset + 1;

//	cout << n_variants << " " << n_haplotypes << endl;

	hsize_t file_offset[2]{start_position_offset, 0};
	hsize_t mem_dims[2]{n_variants, n_haplotypes};

	unique_ptr<unsigned char[]> haplotypes = unique_ptr<unsigned char[]>(new unsigned char[n_variants * n_haplotypes]);

	if ((dataset_id = H5Dopen(chromosome_it->second->get(), HAPLOTYPES_DATASET, H5P_DEFAULT)) < 0) {
		throw HVCFReadException(__FILE__, __FUNCTION__, __LINE__, "Error while opening dataset.");
	}

	if ((file_dataspace_id = H5Dget_space(dataset_id)) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while getting dataspace.");
	}

	if ((memory_dataspace_id = H5Screate_simple(2, mem_dims, nullptr)) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while creating memory dataspace.");
	}

	if (H5Sselect_hyperslab(file_dataspace_id, H5S_SELECT_SET, file_offset, NULL, mem_dims, NULL) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while making selection in dataspace.");
	}

	if (H5Dread(dataset_id, H5T_NATIVE_UCHAR, memory_dataspace_id, file_dataspace_id, H5P_DEFAULT, haplotypes.get()) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while reading from dataset");
	}

//	for (long long int i = 0; i < n_variants; ++i) {
//		cout << (start_position_offset + i) << "=";
//		for (hsize_t j = 0; j < n_haplotypes; ++j) {
//			cout << ((int)haplotypes[i * n_haplotypes + j]);
//		}
//		cout << endl;
//	}

}

}
