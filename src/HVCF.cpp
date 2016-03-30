#include "include/HVCF.h"

namespace sph_umich_edu {

constexpr char HVCF::SAMPLES_GROUP[];
constexpr char HVCF::SAMPLES_ALL_DATASET[];
constexpr char HVCF::CHROMOSOMES_GROUP[];
constexpr char HVCF::VARIANTS_DATASET[];
constexpr char HVCF::HAPLOTYPES_DATASET[];
constexpr char HVCF::VARIABLE_LENGTH_STRING_TYPE[];
constexpr char HVCF::VARIANTS_ENTRY_TYPE[];
constexpr char HVCF::STRING_INDEX_KEY_TYPE[];
constexpr char HVCF::ULL_INDEX_KEY_TYPE[];
constexpr char HVCF::VARIANT_NAMES_INDEX[];
constexpr char HVCF::VARIANT_POSITIONS_INDEX[];


HVCF::HVCF() {
//	Disables automatic HDF5 error stack printing to stderr when function call returns negative value.
//	H5Eset_auto(H5E_DEFAULT, nullptr, nullptr);
}

HVCF::~HVCF() {

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

hid_t HVCF::create_haplotypes_dataset(hid_t group_id, hsize_t variants_chunk_size, hsize_t samples_chunk_size) throw (HVCFWriteException) {
	HDF5DataspaceIdentifier dataspace_id;
	HDF5DatasetIdentifier dataset_id;
	HDF5PropertyIdentifier dataset_property_id;

	hsize_t n_samples = get_n_samples();
	hsize_t n_haplotypes = n_samples + n_samples;

	if (samples_chunk_size > n_haplotypes) {
		samples_chunk_size = n_haplotypes;
	}

	hsize_t initial_dims[2]{0, n_haplotypes};
	hsize_t maximum_dims[2]{H5S_UNLIMITED, n_haplotypes};
	hsize_t chunk_dims[2]{variants_chunk_size, samples_chunk_size};

	if ((dataspace_id = H5Screate_simple(2, initial_dims, maximum_dims)) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while creating dataspace.");
	}

	if ((dataset_property_id = H5Pcreate(H5P_DATASET_CREATE)) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while creating dataset property.");
	}

	if ((H5Pset_chunk(dataset_property_id, 2, chunk_dims) < 0) || (H5Pset_deflate(dataset_property_id, 9) < 0)) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while setting dataset properties.");
	}

//	if ((H5Pset_chunk(dataset_property_id, 2, chunk_dims) < 0) || (H5Pset_szip(dataset_property_id, H5_SZIP_NN_OPTION_MASK, 32) < 0)) {
//		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while setting dataset properties.");
//	}

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

	if ((H5Pset_chunk(dataset_property_id, 1, chunk_dims) < 0) || (H5Pset_deflate(dataset_property_id, 9) < 0)) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while setting dataset properties.");
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

	dataset_id = create_haplotypes_dataset(group_id, 2000, 4000);
	dataset_id.close();

	dataset_id = create_variants_dataset(group_id, 2000);
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
	HDF5DatatypeIdentifier native_string_datatype_id;
	HDF5DatatypeIdentifier memory_datatype_id;
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
	if ((native_string_datatype_id = H5Topen(file_id, VARIABLE_LENGTH_STRING_TYPE, H5P_DEFAULT)) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while opening datatype.");
	}

	if ((memory_datatype_id = H5Tcreate(H5T_COMPOUND, sizeof(variants_entry_type))) < 0) {
		throw HVCFReadException(__FILE__, __FUNCTION__, __LINE__, "Error while creating compound datatype.");
	}

	if (H5Tinsert(memory_datatype_id, "name", HOFFSET(variants_entry_type, name), native_string_datatype_id) < 0) {
		throw HVCFReadException(__FILE__, __FUNCTION__, __LINE__, "Error while adding new member to compound datatype.");
	}

	if (H5Tinsert(memory_datatype_id, "position", HOFFSET(variants_entry_type, position), H5T_NATIVE_ULLONG) < 0) {
		throw HVCFReadException(__FILE__, __FUNCTION__, __LINE__, "Error while adding new member to compound datatype.");
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

	if (H5Dwrite(dataset_id, memory_datatype_id, memory_dataspace_id, file_dataspace_id, H5P_DEFAULT, buffer) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while writing to dataset.");
	}
	// END: write.
}

void HVCF::write_positions_index_bucket(hid_t chromosome_group_id, const string& hash, const vector<hsize_t>& offsets) throw (HVCFWriteException) {
	if (hash.length() == 0u) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Empty group name.");
	}

	if (offsets.size() == 0u) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Empty list of offsets");
	}

	HDF5DatatypeIdentifier native_string_datatype_id;
	HDF5DatatypeIdentifier index_key_type_id;
	HDF5DatatypeIdentifier variants_entry_memory_datatype_id;
	HDF5DatatypeIdentifier ull_index_key_memory_datatype_id;

	HDF5DatasetIdentifier dataset_id;
	HDF5DataspaceIdentifier file_dataspace_id;
	HDF5DataspaceIdentifier memory_dataspace_id;
	HDF5PropertyIdentifier dataset_property_id;
	HDF5GroupIdentifier index_group_id;

	hsize_t dims[1]{offsets.size()};
	hsize_t mem_dims[1]{offsets.size()};

	variants_entry_type variants_buffer[offsets.size()];
	ull_index_key_type keys_buffer[offsets.size()];

	// BEGIN: open comitted datatype.
	if ((native_string_datatype_id = H5Topen(file_id, VARIABLE_LENGTH_STRING_TYPE, H5P_DEFAULT)) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while opening datatype.");
	}

	if ((index_key_type_id = H5Topen(file_id, ULL_INDEX_KEY_TYPE, H5P_DEFAULT)) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while opening datatype.");
	}
	// END: open comitted datatype.

	// BEGIN: create memory datatype for variants.
	if ((variants_entry_memory_datatype_id = H5Tcreate(H5T_COMPOUND, sizeof(variants_entry_type))) < 0) {
		throw HVCFReadException(__FILE__, __FUNCTION__, __LINE__, "Error while creating compound datatype.");
	}

	if (H5Tinsert(variants_entry_memory_datatype_id, "name", HOFFSET(variants_entry_type, name), native_string_datatype_id) < 0) {
		throw HVCFReadException(__FILE__, __FUNCTION__, __LINE__, "Error while adding new member to compound datatype.");
	}

	if (H5Tinsert(variants_entry_memory_datatype_id, "position", HOFFSET(variants_entry_type, position), H5T_NATIVE_ULLONG) < 0) {
		throw HVCFReadException(__FILE__, __FUNCTION__, __LINE__, "Error while adding new member to compound datatype.");
	}
	// END: create memory datatype for variants.

	// BEGIN: create memory datatype for index.
	if ((ull_index_key_memory_datatype_id = H5Tcreate(H5T_COMPOUND, sizeof(ull_index_key_type))) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while creating compound datatype.");
	}

	if (H5Tinsert(ull_index_key_memory_datatype_id, "ull_value", HOFFSET(ull_index_key_type, ull_value), H5T_NATIVE_ULLONG) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while adding new member to compound datatype.");
	}

	if (H5Tinsert(ull_index_key_memory_datatype_id, "offset", HOFFSET(ull_index_key_type, offset), H5T_NATIVE_HSIZE) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while adding new member to compound datatype.");
	}
	// END: create memory datatype index.

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

	file_dataspace_id.close();
	memory_dataspace_id.close();
	dataset_id.close();

	for (unsigned int i = 0; i < offsets.size(); ++i) {
		keys_buffer[i].ull_value = variants_buffer[i].position;
		keys_buffer[i].offset = offsets[i];
	}

	sort(keys_buffer, keys_buffer + offsets.size(),
			[] (const ull_index_key_type& f, const ull_index_key_type& s) -> bool {
				return (f.ull_value < s.ull_value);
	});

	// BEGIN: write bucket.
	if ((index_group_id = H5Gopen(chromosome_group_id, VARIANT_POSITIONS_INDEX, H5P_DEFAULT)) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while opening group.");
	}

	if ((file_dataspace_id = H5Screate_simple(1, dims, NULL)) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while creating dataspace.");
	}

	if ((dataset_property_id = H5Pcreate(H5P_DATASET_CREATE)) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while creating dataset property.");
	}

	if ((dataset_id = H5Dcreate(index_group_id, hash.c_str(), index_key_type_id, file_dataspace_id, H5P_DEFAULT, dataset_property_id, H5P_DEFAULT)) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while creating dataset.");
	}

	if (H5Dwrite(dataset_id, ull_index_key_memory_datatype_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, keys_buffer) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while writing to dataset.");
	}
	// END: write bucket.
}

void HVCF::write_names_index_bucket(hid_t chromosome_group_id, const string& hash, const vector<hsize_t>& offsets) throw (HVCFWriteException) {
	if (hash.length() == 0u) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Empty group name.");
	}

	if (offsets.size() == 0u) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Empty list of offsets");
	}

	HDF5DatatypeIdentifier native_string_datatype_id;
	HDF5DatatypeIdentifier index_key_type_id;
	HDF5DatatypeIdentifier variants_entry_memory_datatype_id;
	HDF5DatatypeIdentifier string_index_key_memory_datatype_id;

	HDF5DatasetIdentifier dataset_id;
	HDF5DataspaceIdentifier file_dataspace_id;
	HDF5DataspaceIdentifier memory_dataspace_id;
	HDF5PropertyIdentifier dataset_property_id;
	HDF5GroupIdentifier index_group_id;

	hsize_t dims[1]{offsets.size()};
	hsize_t mem_dims[1]{offsets.size()};

	variants_entry_type variants_buffer[offsets.size()];
	string_index_key_type index_buffer[offsets.size()];

	// BEGIN: open comitted datatype.
	if ((native_string_datatype_id = H5Topen(file_id, VARIABLE_LENGTH_STRING_TYPE, H5P_DEFAULT)) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while opening datatype.");
	}

	if ((index_key_type_id = H5Topen(file_id, STRING_INDEX_KEY_TYPE, H5P_DEFAULT)) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while opening datatype.");
	}
	// END: open comitted datatype.

	// BEGIN: create memory datatype for variants.
	if ((variants_entry_memory_datatype_id = H5Tcreate(H5T_COMPOUND, sizeof(variants_entry_type))) < 0) {
		throw HVCFReadException(__FILE__, __FUNCTION__, __LINE__, "Error while creating compound datatype.");
	}

	if (H5Tinsert(variants_entry_memory_datatype_id, "name", HOFFSET(variants_entry_type, name), native_string_datatype_id) < 0) {
		throw HVCFReadException(__FILE__, __FUNCTION__, __LINE__, "Error while adding new member to compound datatype.");
	}

	if (H5Tinsert(variants_entry_memory_datatype_id, "position", HOFFSET(variants_entry_type, position), H5T_NATIVE_ULLONG) < 0) {
		throw HVCFReadException(__FILE__, __FUNCTION__, __LINE__, "Error while adding new member to compound datatype.");
	}
	// END: create memory datatype for variants.

	// BEGIN: create memory datatype for index.
	if ((string_index_key_memory_datatype_id = H5Tcreate(H5T_COMPOUND, sizeof(string_index_key_type))) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while creating compound datatype.");
	}

	if (H5Tinsert(string_index_key_memory_datatype_id, "string_value", HOFFSET(string_index_key_type, string_value), native_string_datatype_id) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while adding new member to compound datatype.");
	}

	if (H5Tinsert(string_index_key_memory_datatype_id, "offset", HOFFSET(string_index_key_type, offset), H5T_NATIVE_HSIZE) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while adding new member to compound datatype.");
	}
	// END: create memory datatype for index.

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

	for (unsigned int i = 0; i < offsets.size(); ++i) {
		index_buffer[i].string_value = variants_buffer[i].name;
		index_buffer[i].offset = offsets[i];
	}

	file_dataspace_id.close();
	dataset_id.close();

	sort(index_buffer, index_buffer + offsets.size(),
			[] (const string_index_key_type& f, const string_index_key_type& s) -> bool {
				return (strcmp(f.string_value, s.string_value) < 0);
	});

	// BEGIN: write bucket.
	if ((index_group_id = H5Gopen(chromosome_group_id, VARIANT_NAMES_INDEX, H5P_DEFAULT)) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while opening group.");
	}

	if ((file_dataspace_id = H5Screate_simple(1, dims, NULL)) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while creating dataspace.");
	}

	if ((dataset_property_id = H5Pcreate(H5P_DATASET_CREATE)) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while creating dataset property.");
	}

	if ((dataset_id = H5Dcreate(index_group_id, hash.c_str(), index_key_type_id, file_dataspace_id, H5P_DEFAULT, dataset_property_id, H5P_DEFAULT)) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while creating dataset.");
	}

	if (H5Dwrite(dataset_id, string_index_key_memory_datatype_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, index_buffer) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while writing to dataset.");
	}
	// END: write bucket.

	if (H5Dvlen_reclaim(variants_entry_memory_datatype_id, memory_dataspace_id, H5P_DEFAULT, variants_buffer) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while reclaiming HDF5 memory.");
	}
}

void HVCF::create_indices(hid_t chromosome_group_id) throw (HVCFWriteException) {
	HDF5GroupIdentifier names_index_group_id;
	HDF5GroupIdentifier positions_index_group_id;
	HDF5DatasetIdentifier dataset_id;
	HDF5DataspaceIdentifier file_dataspace_id;
	HDF5DatatypeIdentifier native_string_datatype_id;
	HDF5DatatypeIdentifier memory_datatype_id;
	HDF5DataspaceIdentifier memory_dataspace_id;

	const unsigned int read_chunk_size = 100000;

	hsize_t file_dims[1]{0};
	hsize_t file_offset[1]{0};
	hsize_t mem_dims[1]{read_chunk_size};

	variants_entry_type buffer[read_chunk_size];

	unordered_map<unsigned int, vector<hsize_t>> names_index_buckets;
	auto names_index_buckets_it = names_index_buckets.end();

	unordered_map<unsigned int, vector<hsize_t>> positions_index_buckets;
	auto positions_index_buckets_it = positions_index_buckets.end();

	std::hash<string> hash_function;
	size_t hash_value;

	// BEGIN: create groups for indices.
	if ((names_index_group_id = H5Gcreate(chromosome_group_id, VARIANT_NAMES_INDEX, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT)) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while creating group.");
	}

	if ((positions_index_group_id = H5Gcreate(chromosome_group_id, VARIANT_POSITIONS_INDEX, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT)) < 0) {
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
	if ((native_string_datatype_id = H5Topen(file_id, VARIABLE_LENGTH_STRING_TYPE, H5P_DEFAULT)) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while opening datatype.");
	}

	if ((memory_datatype_id = H5Tcreate(H5T_COMPOUND, sizeof(variants_entry_type))) < 0) {
		throw HVCFReadException(__FILE__, __FUNCTION__, __LINE__, "Error while creating compound datatype.");
	}

	if (H5Tinsert(memory_datatype_id, "name", HOFFSET(variants_entry_type, name), native_string_datatype_id) < 0) {
		throw HVCFReadException(__FILE__, __FUNCTION__, __LINE__, "Error while adding new member to compound datatype.");
	}

	if (H5Tinsert(memory_datatype_id, "position", HOFFSET(variants_entry_type, position), H5T_NATIVE_ULLONG) < 0) {
		throw HVCFReadException(__FILE__, __FUNCTION__, __LINE__, "Error while adding new member to compound datatype.");
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

			hash_value = buffer[i].position % N_HASH_BUCKETS;
			positions_index_buckets_it = positions_index_buckets.find(hash_value);
			if (positions_index_buckets_it == positions_index_buckets.end()) {
				positions_index_buckets_it = positions_index_buckets.emplace(hash_value, vector<hsize_t>()).first;
			}
			positions_index_buckets_it->second.push_back(file_offset[0]);

			file_offset[0] += 1;
		}

		if (H5Dvlen_reclaim(memory_datatype_id, memory_dataspace_id, H5P_DEFAULT, buffer) < 0) {
			throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while reclaiming HDF5 memory.");
		}

		memory_dataspace_id.close();
	}
	// END: read data in chunks.

	// BEGIN: write indices on disk.
	for (auto&& bucket : names_index_buckets) {
		if (bucket.second.size() == 0) {
			continue;
		}
		write_names_index_bucket(chromosome_group_id, to_string(bucket.first), bucket.second);
	}

	for (auto&& bucket : positions_index_buckets) {
		if (bucket.second.size() == 0) {
			continue;
		}
		write_positions_index_bucket(chromosome_group_id, to_string(bucket.first), bucket.second);
	}
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

	if (H5Tcommit(file_id, ULL_INDEX_KEY_TYPE, datatype_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT) < 0) {
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

	if	(H5Tcommit(file_id, STRING_INDEX_KEY_TYPE, datatype_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT) < 0) {
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
		//flush!
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

	if ((file_id = H5Fopen(this->name.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT)) < 0) {
		throw HVCFOpenException(__FILE__, __FUNCTION__, __LINE__, "Error while opening file.");
	}

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

long long int HVCF::get_variant_offset_by_position(const string& chromosome, unsigned long long int position) throw (HVCFReadException) {
	auto chromosomes_it = chromosomes.find(chromosome);

	if (chromosomes_it == chromosomes.end()) {
		return -1;
	}

	HDF5GroupIdentifier index_group_id;
	HDF5DatasetIdentifier dataset_id;
	HDF5DataspaceIdentifier dataspace_id;

	if ((index_group_id = H5Gopen(chromosomes_it->second->get(), VARIANT_POSITIONS_INDEX, H5P_DEFAULT)) < 0) {
		throw HVCFReadException(__FILE__, __FUNCTION__, __LINE__, "Error while opening group.");
	}

	string hash_value = std::to_string(position % N_HASH_BUCKETS);

	hsize_t file_dims[1]{0};

	htri_t exists = H5Lexists(index_group_id, hash_value.c_str(), H5P_DEFAULT);
	if (exists < 0) {
		throw HVCFReadException(__FILE__, __FUNCTION__, __LINE__, "Error while hash bucket.");
	}

	if (exists != true) {
		return -1;
	}

	if ((dataset_id = H5Dopen(index_group_id, hash_value.c_str(), H5P_DEFAULT)) < 0) {
		throw HVCFReadException(__FILE__, __FUNCTION__, __LINE__, "Error while opening dataset.");
	}

	if ((dataspace_id = H5Dget_space(dataset_id)) < 0) {
		throw HVCFReadException(__FILE__, __FUNCTION__, __LINE__, "Error while getting dataspace.");
	}

	if (H5Sget_simple_extent_dims(dataspace_id, file_dims, NULL) < 0) {
		throw HVCFReadException(__FILE__, __FUNCTION__, __LINE__, "Error while getting dataspace dimensions.");
	}

	ull_index_key_type keys_buffer[file_dims[0]];

	ull_index_key_type search_key;
	search_key.ull_value = position;

	HDF5DatatypeIdentifier memory_index_key_type;

	// BEGIN: create memory datatype.
	if ((memory_index_key_type = H5Tcreate(H5T_COMPOUND, sizeof(ull_index_key_type))) < 0) {
		throw HVCFReadException(__FILE__, __FUNCTION__, __LINE__, "Error while creating compound datatype.");
	}

	if (H5Tinsert(memory_index_key_type, "ull_value", HOFFSET(ull_index_key_type, ull_value), H5T_NATIVE_ULLONG) < 0) {
		throw HVCFReadException(__FILE__, __FUNCTION__, __LINE__, "Error while adding new member to compound datatype.");
	}

	if (H5Tinsert(memory_index_key_type, "offset", HOFFSET(ull_index_key_type, offset), H5T_NATIVE_HSIZE) < 0) {
		throw HVCFReadException(__FILE__, __FUNCTION__, __LINE__, "Error while adding new member to compound datatype.");
	}
	// END: create memory datatype.

	if (H5Dread(dataset_id, memory_index_key_type, H5S_ALL, H5S_ALL, H5P_DEFAULT, keys_buffer) < 0) {
		throw HVCFReadException(__FILE__, __FUNCTION__, __LINE__, "Error while reading from dataset.");
	}

	auto result = lower_bound(keys_buffer, keys_buffer + file_dims[0], search_key,
			[] (const ull_index_key_type& f, const ull_index_key_type& s) -> bool {
				return (f.ull_value < s.ull_value);
			});

	if ((result == keys_buffer + file_dims[0]) || (result->ull_value != position)) {
		return -1;
	}

	return result->offset;
}

long long int HVCF::get_variant_offset_by_name(const string& chromosome, const string& name) throw (HVCFReadException) {
	auto chromosomes_it = chromosomes.find(chromosome);

	if (chromosomes_it == chromosomes.end()) {
		return -1;
	}

	HDF5GroupIdentifier index_group_id;
	HDF5DatasetIdentifier dataset_id;
	HDF5DataspaceIdentifier dataspace_id;
	HDF5DatatypeIdentifier native_string_datatype_id;

	std::hash<string> hash_function;

	if (chromosomes_it == chromosomes.end()) {
		return -1;
	}

	if ((native_string_datatype_id = H5Topen(file_id, VARIABLE_LENGTH_STRING_TYPE, H5P_DEFAULT)) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while opening datatype.");
	}


	if ((index_group_id = H5Gopen(chromosomes_it->second->get(), VARIANT_NAMES_INDEX, H5P_DEFAULT)) < 0) {
		throw HVCFReadException(__FILE__, __FUNCTION__, __LINE__, "Error while opening group.");
	}

	string hash_value = std::to_string(hash_function(name) % N_HASH_BUCKETS);

	htri_t exists = H5Lexists(index_group_id, hash_value.c_str(), H5P_DEFAULT);
	if (exists < 0) {
		throw HVCFReadException(__FILE__, __FUNCTION__, __LINE__, "Error while hash bucket.");
	}

	if (exists != true) {
		return -1;
	}

	hsize_t file_dims[1]{0};

	if ((dataset_id = H5Dopen(index_group_id, hash_value.c_str(), H5P_DEFAULT)) < 0) {
		throw HVCFReadException(__FILE__, __FUNCTION__, __LINE__, "Error while opening dataset.");
	}

	if ((dataspace_id = H5Dget_space(dataset_id)) < 0) {
		throw HVCFReadException(__FILE__, __FUNCTION__, __LINE__, "Error while getting dataspace.");
	}

	if (H5Sget_simple_extent_dims(dataspace_id, file_dims, NULL) < 0) {
		throw HVCFReadException(__FILE__, __FUNCTION__, __LINE__, "Error while getting dataspace dimensions.");
	}

	string_index_key_type keys_buffer[file_dims[0]];

	string_index_key_type search_key;
	search_key.string_value = (char*)name.c_str();

	HDF5DatatypeIdentifier memory_index_key_type;

	// BEGIN: create memory datatype.
	if ((memory_index_key_type = H5Tcreate(H5T_COMPOUND, sizeof(string_index_key_type))) < 0) {
		throw HVCFReadException(__FILE__, __FUNCTION__, __LINE__, "Error while creating compound datatype.");
	}

	if (H5Tinsert(memory_index_key_type, "string_value", HOFFSET(string_index_key_type, string_value), native_string_datatype_id) < 0) {
		throw HVCFReadException(__FILE__, __FUNCTION__, __LINE__, "Error while adding new member to compound datatype.");
	}

	if (H5Tinsert(memory_index_key_type, "offset", HOFFSET(string_index_key_type, offset), H5T_NATIVE_HSIZE) < 0) {
		throw HVCFReadException(__FILE__, __FUNCTION__, __LINE__, "Error while adding new member to compound datatype.");
	}
	// END: create memory datatype.

	if (H5Dread(dataset_id, memory_index_key_type, H5S_ALL, H5S_ALL, H5P_DEFAULT, keys_buffer) < 0) {
		throw HVCFReadException(__FILE__, __FUNCTION__, __LINE__, "Error while reading from dataset.");
	}

	auto result = lower_bound(keys_buffer, keys_buffer + file_dims[0], search_key,
			[] (const string_index_key_type& f, const string_index_key_type& s) -> bool {
				return (strcmp(f.string_value, s.string_value) < 0);
			});

	long long int offset = 0;

	if ((result == keys_buffer + file_dims[0]) || (strcmp(result->string_value, name.c_str())) != 0) {
		offset = -1;
	} else {
		offset = result->offset;
	}

	if (H5Dvlen_reclaim(memory_index_key_type, dataspace_id, H5P_DEFAULT, keys_buffer) < 0) {
		throw HVCFReadException(__FILE__, __FUNCTION__, __LINE__, "Error while reclaiming HDF5 memory.");
	}

	return offset;
}

void HVCF::chunk_read_test(const string& chromosome, unsigned long long int start_position, unsigned long long end_position) throw (HVCFReadException) {
	auto chromosome_it = chromosomes.find(chromosome); //TODO: check if chromosome exists

	long long int start_position_offset = get_variant_offset_by_position(chromosome, start_position); // TODO: put 0 if not exists
	long long int end_position_offset = get_variant_offset_by_position(chromosome, end_position); // TODO: put file_dims[0] - 1 offset if not exists

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

void HVCF::compute_ld(const string& chromosome, unsigned long long int start_position, unsigned long long end_position) throw (HVCFReadException) {
	auto chromosome_it = chromosomes.find(chromosome); //TODO: check if chromosome exists

	long long int start_position_offset = get_variant_offset_by_position(chromosome, start_position); // TODO: put 0 if not exists
	long long int end_position_offset = get_variant_offset_by_position(chromosome, end_position); // TODO: put file_dims[0] - 1 offset if not exists


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

//	std::chrono::time_point<std::chrono::system_clock> start, end;
//	std::chrono::duration<double> elapsed_seconds;
//
//	start = std::chrono::system_clock::now();
	unique_ptr<double[]> dbl_haplotypes = unique_ptr<double[]>(new double[n_variants * n_haplotypes]);
	for (long long int i = 0; i < n_variants * n_haplotypes; ++i) {
		dbl_haplotypes[i] = (double)haplotypes[i];
	}
//	end = std::chrono::system_clock::now();
//	elapsed_seconds = end - start;
//	cout << "1: " << elapsed_seconds.count() << endl;

//	Mat<double> S(dbl_haplotypes.get(), n_haplotypes, n_variants, false, false); // don't copy matrix.
//	Row<double> J(n_haplotypes, fill::ones);
//
//	Mat<double> C1(J * S);
//	Mat<double> C2(n_haplotypes - C1);
//	Mat<double> M1(C1.t() * C1);
//
//	Mat<double> R((n_haplotypes * S.t() * S - M1) / sqrt(M1 % (C2.t() * C2)));

//	start = std::chrono::system_clock::now();
	Mat<double> S(dbl_haplotypes.get(), n_haplotypes, n_variants, false, false); // don't copy matrix.
//	end = std::chrono::system_clock::now();
//	elapsed_seconds = end - start;
//	cout << "2: " << elapsed_seconds.count() << endl;

//	start = std::chrono::system_clock::now();
	Row<double> J(n_haplotypes, fill::ones);
//	end = std::chrono::system_clock::now();
//	elapsed_seconds = end - start;
//	cout << "3: " << elapsed_seconds.count() << endl;

//	start = std::chrono::system_clock::now();
	Mat<double> C1(J * S);
	Mat<double> C2(n_haplotypes - C1);
	Mat<double> M1(C1.t() * C1);
//	end = std::chrono::system_clock::now();
//	elapsed_seconds = end - start;
//	cout << "4: " << elapsed_seconds.count() << endl;

//	start = std::chrono::system_clock::now();
	Mat<double> R((n_haplotypes * S.t() * S - M1) / sqrt(M1 % (C2.t() * C2)));
//	end = std::chrono::system_clock::now();
//	elapsed_seconds = end - start;
//	cout << "5: " << elapsed_seconds.count() << endl;

//	cout << "R:" << endl;
//	R.raw_print();
//
//	cout << "Rsq:" << endl;
//	cout << square(R) << endl;

}

void HVCF::compute_ld(const string& chromosome, const string& lead_variant_name, unsigned long long int start_position, unsigned long long end_position) throw (HVCFReadException) {
	auto chromosome_it = chromosomes.find(chromosome); //TODO: check if chromosome exists

	long long int lead_variant_offset = get_variant_offset_by_name(chromosome, lead_variant_name); // TODO: check if exists
	long long int start_position_offset = get_variant_offset_by_position(chromosome, start_position); // TODO: put 0 if not exists
	long long int end_position_offset = get_variant_offset_by_position(chromosome, end_position); // TODO: put file_dims[0] - 1 offset if not exists

	HDF5DatasetIdentifier dataset_id;
	HDF5DataspaceIdentifier file_dataspace_id;
	HDF5DataspaceIdentifier memory_dataspace_id;

	hsize_t n_samples = get_n_samples(); // TODO: compute this only once on file opening
	hsize_t n_haplotypes = n_samples + n_samples; // TODO: compute this only once on file opening

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

//	cout << lead_variant_offset << " " << start_position_offset << " " << end_position_offset << endl;
//	cout << n_variants << " " << n_haplotypes << endl;

	hsize_t file_offset1[2]{lead_variant_offset, 0};
	hsize_t counts1[2]{1, n_haplotypes};
	hsize_t file_offset2[2]{start_position_offset, 0};
	hsize_t counts2[2]{end_position_offset - start_position_offset + 1, n_haplotypes};
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

	if (H5Sselect_hyperslab(file_dataspace_id, H5S_SELECT_SET, file_offset1, NULL, counts1, NULL) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while making selection in dataspace.");
	}

	if (H5Sselect_hyperslab(file_dataspace_id, H5S_SELECT_OR, file_offset2, NULL, counts2, NULL) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while making selection in dataspace.");
	}

	if (H5Dread(dataset_id, H5T_NATIVE_UCHAR, memory_dataspace_id, file_dataspace_id, H5P_DEFAULT, haplotypes.get()) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while reading from dataset");
	}

	unique_ptr<double[]> dbl_haplotypes = unique_ptr<double[]>(new double[n_variants * n_haplotypes]);
	for (long long int i = 0; i < n_variants * n_haplotypes; ++i) {
		dbl_haplotypes[i] = (double)haplotypes[i];
	}

	Mat<double> S(dbl_haplotypes.get(), n_haplotypes, n_variants, false, false); // don't copy matrix.
	Row<double> L(S.colptr(lead_variant_local_offset), n_haplotypes);
	Row<double> J(n_haplotypes, fill::ones);

//  also correct
//	Mat<double> C1(J * S);
//	Mat<double> C2(n_haplotypes - C1);
//	double l1 = C1(0, lead_variant_local_offset);
//	double l2 = n_haplotypes - l1;
//	Mat<double> M1(l1 * C1);
//
//	Mat<double> R((n_haplotypes * L * S - M1) / sqrt(M1 % (l2 * C2)));

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



}
