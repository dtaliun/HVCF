#include "include/HVCF.h"

namespace sph_umich_edu {

constexpr char HVCF::SAMPLES_GROUP[];
constexpr char HVCF::POPULATIONS_GROUP[];
constexpr char HVCF::VARIANTS_GROUP[];
constexpr char HVCF::HAPLOTYPES_GROUP[];
constexpr char HVCF::SAMPLE_NAMES_DATASET[];
constexpr char HVCF::VARIANT_NAMES_DATASET[];

HVCF::HVCF():
		file_id(numeric_limits<hid_t>::min()),
		samples_group_id(numeric_limits<hid_t>::min()),
		populations_group_id(numeric_limits<hid_t>::min()),
		variants_group_id(numeric_limits<hid_t>::min()),
		haplotypes_group_id(numeric_limits<hid_t>::min()),
		sample_names_dataset_id(numeric_limits<hid_t>::min()),
		variant_names_dataset_id(numeric_limits<hid_t>::min()),
		native_string_datatype_id(numeric_limits<hid_t>::min()) {

//	Disables automatic HDF5 error stack printing to stderr when function call returns negative value.
//	H5Eset_auto(H5E_DEFAULT, nullptr, nullptr);
}

HVCF::~HVCF() {

}

hid_t HVCF::create_strings_dataset(const string& name, hid_t group_id, hsize_t chunk_size) {
	hsize_t initial_dims[1]{0};
	hsize_t maximum_dims[1]{H5S_UNLIMITED};
	hsize_t chunk_dims[1]{chunk_size};

	hid_t dataspace_id = numeric_limits<hid_t>::min();
	hid_t dataset_creation_plist = numeric_limits<hid_t>::min();
	hid_t dataset_id = numeric_limits<hid_t>::min();

	if ((dataspace_id = H5Screate_simple(1, initial_dims, maximum_dims)) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, this->name.c_str());
	}

	if ((dataset_creation_plist = H5Pcreate(H5P_DATASET_CREATE)) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, this->name.c_str());
	}

	H5Pset_chunk(dataset_creation_plist, 1, chunk_dims);
	H5Pset_deflate(dataset_creation_plist, 9);

	dataset_id = H5Dcreate(group_id, name.c_str(), native_string_datatype_id, dataspace_id, H5P_DEFAULT, dataset_creation_plist, H5P_DEFAULT);

	H5Pclose(dataset_creation_plist);
	H5Sclose(dataspace_id);

	return dataset_id;
}

void HVCF::create(const string& name) throw (HVCFWriteException) {
	this->name = name;

	if ((file_id = H5Fcreate(this->name.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT)) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, name.c_str());
	}

	if ((samples_group_id = H5Gcreate(file_id, SAMPLES_GROUP, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT)) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, name.c_str());
	}

	if ((populations_group_id = H5Gcreate(file_id, POPULATIONS_GROUP, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT)) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, name.c_str());
	}

	if ((variants_group_id = H5Gcreate(file_id, VARIANTS_GROUP, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT)) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, name.c_str());
	}

	if ((haplotypes_group_id = H5Gcreate(file_id, HAPLOTYPES_GROUP, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT)) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, name.c_str());
	}

	if ((native_string_datatype_id = H5Tcopy(H5T_C_S1)) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, name.c_str());
	}
	H5Tset_size(native_string_datatype_id, H5T_VARIABLE);

	if ((sample_names_dataset_id = create_strings_dataset(SAMPLE_NAMES_DATASET, samples_group_id, 100)) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, name.c_str());
	}

	if ((variant_names_dataset_id = create_strings_dataset(VARIANT_NAMES_DATASET, variants_group_id, 10000)) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, name.c_str());
	}
}

void HVCF::write_samples(const string& name, const vector<string>& samples) throw (HVCFWriteException) {


}

void HVCF::open(const string& name) throw (HVCFOpenException) {
	this->name = name;

	if ((file_id = H5Fopen(this->name.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT)) < 0) {
		throw HVCFOpenException(__FILE__, __FUNCTION__, __LINE__, name.c_str());
	}

	if ((samples_group_id = H5Gopen(file_id, SAMPLES_GROUP, H5P_DEFAULT)) < 0) {
		throw HVCFOpenException(__FILE__, __FUNCTION__, __LINE__, name.c_str());
	}

	if ((populations_group_id = H5Gopen(file_id, POPULATIONS_GROUP, H5P_DEFAULT)) < 0) {
		throw HVCFOpenException(__FILE__, __FUNCTION__, __LINE__, name.c_str());
	}

	if ((variants_group_id = H5Gopen(file_id, VARIANTS_GROUP, H5P_DEFAULT)) < 0) {
		throw HVCFOpenException(__FILE__, __FUNCTION__, __LINE__, name.c_str());
	}

	if ((haplotypes_group_id = H5Gopen(file_id, HAPLOTYPES_GROUP, H5P_DEFAULT)) < 0) {
		throw HVCFOpenException(__FILE__, __FUNCTION__, __LINE__, name.c_str());
	}

	if ((sample_names_dataset_id = H5Dopen(samples_group_id, SAMPLE_NAMES_DATASET, H5P_DEFAULT)) < 0) {
		throw HVCFOpenException(__FILE__, __FUNCTION__, __LINE__, name.c_str());
	}

	if ((variant_names_dataset_id = H5Dopen(variants_group_id, VARIANT_NAMES_DATASET, H5P_DEFAULT)) < 0) {
		throw HVCFOpenException(__FILE__, __FUNCTION__, __LINE__, name.c_str());
	}

	if ((native_string_datatype_id = H5Tcopy(H5T_C_S1)) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, name.c_str());
	}
	H5Tset_size(native_string_datatype_id, H5T_VARIABLE);
}

void HVCF::close() throw (HVCFCloseException) {
	if (H5Tclose(native_string_datatype_id)  < 0) {
		throw HVCFCloseException(__FILE__, __FUNCTION__, __LINE__, name.c_str());
	}

	if (H5Dclose(sample_names_dataset_id) < 0) {
		throw HVCFCloseException(__FILE__, __FUNCTION__, __LINE__, name.c_str());
	}

	if (H5Dclose(variant_names_dataset_id) < 0) {
		throw HVCFCloseException(__FILE__, __FUNCTION__, __LINE__, name.c_str());
	}

	if (H5Gclose(haplotypes_group_id) < 0) {
		throw HVCFCloseException(__FILE__, __FUNCTION__, __LINE__, name.c_str());
	}

	if (H5Gclose(variants_group_id) < 0) {
		throw HVCFCloseException(__FILE__, __FUNCTION__, __LINE__, name.c_str());
	}

	if (H5Gclose(populations_group_id) < 0) {
		throw HVCFCloseException(__FILE__, __FUNCTION__, __LINE__, name.c_str());
	}

	if (H5Gclose(samples_group_id) < 0) {
		throw HVCFCloseException(__FILE__, __FUNCTION__, __LINE__, name.c_str());
	}

	if (H5Fclose(file_id) < 0) {
		throw HVCFCloseException(__FILE__, __FUNCTION__, __LINE__, name.c_str());
	}

	haplotypes_group_id = numeric_limits<hid_t>::min();
	variants_group_id = numeric_limits<hid_t>::min();
	populations_group_id = numeric_limits<hid_t>::min();
	samples_group_id = numeric_limits<hid_t>::min();
	file_id = numeric_limits<hid_t>::min();

	name.clear();
}

}
