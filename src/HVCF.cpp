#include "include/HVCF.h"

namespace sph_umich_edu {

constexpr char HVCF::SAMPLES[];
constexpr char HVCF::POPULATIONS[];
constexpr char HVCF::VARIANTS[];
constexpr char HVCF::HAPLOTYPES[];

HVCF::HVCF() :
		file_id(numeric_limits<hid_t>::min()),
		samples_group_id(numeric_limits<hid_t>::min()),
		populations_group_id(numeric_limits<hid_t>::min()),
		variants_group_id(numeric_limits<hid_t>::min()),
		haplotypes_group_id(numeric_limits<hid_t>::min()) {

}

HVCF::~HVCF() {

}

void HVCF::create(const string& name) throw (HVCFWriteException) {
	this->name = name;

	if ((file_id = H5Fcreate(this->name.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT)) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, name.c_str());
	}

	if ((samples_group_id = H5Gcreate(file_id, SAMPLES, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT)) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, name.c_str());
	}

	if ((populations_group_id = H5Gcreate(file_id, POPULATIONS, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT)) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, name.c_str());
	}

	if ((variants_group_id = H5Gcreate(file_id, VARIANTS, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT)) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, name.c_str());
	}

	if ((haplotypes_group_id = H5Gcreate(file_id, HAPLOTYPES, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT)) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, name.c_str());
	}
}

void HVCF::write_samples(const string& name, const vector<string>& samples) throw (HVCFWriteException) {
	hsize_t initial_dims[1]{0};
	hsize_t maximum_dims[1]{H5S_UNLIMITED};
	hsize_t chunk_dims[1]{100};

	hid_t names_dataspace_id = H5Screate_simple(1, initial_dims, maximum_dims);

	hid_t dataset_creation_plist = H5Pcreate(H5P_DATASET_CREATE);
	H5Pset_chunk(dataset_creation_plist, 1, chunk_dims);
	H5Pset_deflate(dataset_creation_plist, 9);

	hid_t native_string_datatype_id = H5Tcopy(H5T_C_S1);
	H5Tset_size(native_string_datatype_id, H5T_VARIABLE);

	hid_t names_dataset_id = H5Dcreate(samples_group_id, name.c_str(), native_string_datatype_id, names_dataspace_id, H5P_DEFAULT, dataset_creation_plist, H5P_DEFAULT);

	H5Pclose(dataset_creation_plist);
	H5Tclose(native_string_datatype_id);
	H5Sclose(names_dataspace_id);
	H5Dclose(names_dataset_id);
}

void HVCF::open(const string& name) throw (HVCFOpenException) {
	this->name = name;

	if ((file_id = H5Fopen(this->name.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT)) < 0) {
		throw HVCFOpenException(__FILE__, __FUNCTION__, __LINE__, name.c_str());
	}

	if ((samples_group_id = H5Gopen(file_id, SAMPLES, H5P_DEFAULT)) < 0) {
		throw HVCFOpenException(__FILE__, __FUNCTION__, __LINE__, name.c_str());
	}

	if ((populations_group_id = H5Gopen(file_id, SAMPLES, H5P_DEFAULT)) < 0) {
		throw HVCFOpenException(__FILE__, __FUNCTION__, __LINE__, name.c_str());
	}

	if ((variants_group_id = H5Gopen(file_id, SAMPLES, H5P_DEFAULT)) < 0) {
		throw HVCFOpenException(__FILE__, __FUNCTION__, __LINE__, name.c_str());
	}

	if ((haplotypes_group_id = H5Gopen(file_id, SAMPLES, H5P_DEFAULT)) < 0) {
		throw HVCFOpenException(__FILE__, __FUNCTION__, __LINE__, name.c_str());
	}
}

void HVCF::close() throw (HVCFCloseException) {
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
