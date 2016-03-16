#include "include/HVCF.h"

namespace sph_umich_edu {

constexpr char HVCF::SAMPLES_GROUP[];
constexpr char HVCF::VARIANTS_GROUP[];
constexpr char HVCF::HAPLOTYPES_GROUP[];
constexpr char HVCF::SAMPLES_ALL_DATASET[];
constexpr char HVCF::VARIANT_NAMES_DATASET[];

HVCF::HVCF():
//		file_id(numeric_limits<hid_t>::min()),
		samples_group_id(numeric_limits<hid_t>::min()),
		variants_group_id(numeric_limits<hid_t>::min()),
		haplotypes_group_id(numeric_limits<hid_t>::min()),
		samples_all_dataset_id(numeric_limits<hid_t>::min()),
		variant_names_dataset_id(numeric_limits<hid_t>::min()),
		native_string_datatype_id(numeric_limits<hid_t>::min()) {

//	Disables automatic HDF5 error stack printing to stderr when function call returns negative value.
//	H5Eset_auto(H5E_DEFAULT, nullptr, nullptr);
}

HVCF::~HVCF() {

}

hid_t HVCF::create_strings_1D_dataset(const string& name, hid_t group_id, hsize_t chunk_size) throw (HVCFWriteException) {
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

hid_t HVCF::create_hsize_1D_dataset(const string&name, hid_t group_id, hsize_t chunk_size) throw (HVCFWriteException) {
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

	dataset_id = H5Dcreate(group_id, name.c_str(), H5T_NATIVE_HSIZE, dataspace_id, H5P_DEFAULT, dataset_creation_plist, H5P_DEFAULT);

	H5Pclose(dataset_creation_plist);
	H5Sclose(dataspace_id);

	return dataset_id;
}

void HVCF::create(const string& name) throw (HVCFWriteException) {
	this->name = name;

//	if ((file_id = H5Fcreate(this->name.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT)) < 0) {
//		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, name.c_str());
//	}

	file_id.set(H5Fcreate(this->name.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT));
	if (file_id.get() < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, name.c_str());
	}

	if ((samples_group_id = H5Gcreate(file_id.get(), SAMPLES_GROUP, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT)) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, name.c_str());
	}

	if ((variants_group_id = H5Gcreate(file_id.get(), VARIANTS_GROUP, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT)) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, name.c_str());
	}

	if ((haplotypes_group_id = H5Gcreate(file_id.get(), HAPLOTYPES_GROUP, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT)) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, name.c_str());
	}

	if ((native_string_datatype_id = H5Tcopy(H5T_C_S1)) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, name.c_str());
	}
	H5Tset_size(native_string_datatype_id, H5T_VARIABLE);

	if ((samples_all_dataset_id = create_strings_1D_dataset(SAMPLES_ALL_DATASET, samples_group_id, 100)) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, name.c_str());
	}

	if ((variant_names_dataset_id = create_strings_1D_dataset(VARIANT_NAMES_DATASET, variants_group_id, 10000)) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, name.c_str());
	}
}

void HVCF::set_samples(const vector<string>& samples) throw (HVCFWriteException) {
	hid_t file_dataspace_id = numeric_limits<hid_t>::min();
	hid_t memory_dataspace_id = numeric_limits<hid_t>::min();
	hsize_t mem_dims[1]{samples.size()};
	hsize_t file_dims[1]{0};
	hsize_t file_offset[1]{0};
	const char* buffer[samples.size()];

	for (unsigned int i = 0u; i < samples.size(); ++i) {
		buffer[i] = samples.at(i).c_str();
	}

	file_dataspace_id = H5Dget_space(samples_all_dataset_id);
	H5Sget_simple_extent_dims(file_dataspace_id, file_dims, nullptr);
	H5Sclose(file_dataspace_id);

	file_offset[0] = file_dims[0];
	file_dims[0] += mem_dims[0];

	if (H5Dset_extent(samples_all_dataset_id, file_dims) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, name.c_str());
	}

	file_dataspace_id = H5Dget_space(samples_all_dataset_id);
	memory_dataspace_id = H5Screate_simple(1, mem_dims, nullptr);

	if (H5Sselect_hyperslab(file_dataspace_id, H5S_SELECT_SET, file_offset, nullptr, mem_dims, nullptr) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, name.c_str());
	}

	if (H5Dwrite(samples_all_dataset_id, native_string_datatype_id, memory_dataspace_id, file_dataspace_id, H5P_DEFAULT, buffer) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, name.c_str());
	}

	H5Sclose(memory_dataspace_id);
	H5Sclose(file_dataspace_id);
}

void HVCF::set_population(const string& name, const vector<string>& samples) throw (HVCFWriteException) {
	map<string, unsigned int> all_samples;
	auto all_samples_it = all_samples.end();
	hsize_t buffer[samples.size()];
	hid_t population_dataset_id = numeric_limits<hid_t>::min();
	hid_t file_dataspace_id = numeric_limits<hid_t>::min();
	hid_t memory_dataspace_id = numeric_limits<hid_t>::min();
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
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, name.c_str());
	}

	file_dataspace_id = H5Dget_space(population_dataset_id);
	memory_dataspace_id = H5Screate_simple(1, mem_dims, nullptr);

	if (H5Sselect_hyperslab(file_dataspace_id, H5S_SELECT_SET, file_offset, nullptr, mem_dims, nullptr) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, name.c_str());
	}

	if (H5Dwrite(population_dataset_id, H5T_NATIVE_HSIZE, memory_dataspace_id, file_dataspace_id, H5P_DEFAULT, buffer) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, name.c_str());
	}

	H5Sclose(file_dataspace_id);
	H5Sclose(memory_dataspace_id);
	H5Dclose(population_dataset_id);
}

vector<string> HVCF::get_samples() throw (HVCFReadException) {
	vector<string> samples;

	hsize_t file_dims[1]{0};

	hid_t file_dataspace_id = numeric_limits<hid_t>::min();

	file_dataspace_id = H5Dget_space(samples_all_dataset_id);
	H5Sget_simple_extent_dims(file_dataspace_id, file_dims, nullptr);
	H5Sclose(file_dataspace_id);

	char* buffer[file_dims[0]];

	if (H5Dread(samples_all_dataset_id, native_string_datatype_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, buffer) < 0) {
		throw HVCFReadException(__FILE__, __FUNCTION__, __LINE__, name.c_str());
	}

	for (unsigned int i = 0; i < file_dims[0]; ++i) {
		samples.emplace_back(buffer[i]);
	}

	file_dataspace_id = H5Dget_space(samples_all_dataset_id);
	H5Dvlen_reclaim(native_string_datatype_id, file_dataspace_id, H5P_DEFAULT, buffer);
	H5Sclose(file_dataspace_id);

	return samples;
}

vector<string> HVCF::get_population(const string& name) throw (HVCFReadException) {
	vector<string> samples;
	hsize_t file_dims[1]{0};
	hsize_t mem_dims[1]{0};
	hid_t dataset_id = numeric_limits<hid_t>::min();
	hid_t file_dataspace_id = numeric_limits<hid_t>::min();
	hid_t mem_dataspace_id = numeric_limits<hid_t>::min();

	if ((dataset_id = H5Dopen(samples_group_id, name.c_str(), H5P_DEFAULT)) < 0) {
		throw HVCFReadException(__FILE__, __FUNCTION__, __LINE__, "Error while opening dataset.");
	}

	file_dataspace_id = H5Dget_space(dataset_id);
	H5Sget_simple_extent_dims(file_dataspace_id, file_dims, nullptr);
	H5Sclose(file_dataspace_id);

	mem_dims[0] = file_dims[0];

	hsize_t index_buffer[file_dims[0]];
	char* samples_buffer[file_dims[0]];

	if (H5Dread(dataset_id, H5T_NATIVE_HSIZE, H5S_ALL, H5S_ALL, H5P_DEFAULT, index_buffer) < 0) {
		throw HVCFReadException(__FILE__, __FUNCTION__, __LINE__, this->name.c_str());
	}

	H5Dclose(dataset_id);

	file_dataspace_id = H5Dget_space(samples_all_dataset_id);
	if (H5Sselect_elements(file_dataspace_id, H5S_SELECT_SET, file_dims[0], index_buffer) < 0) {
		throw HVCFReadException(__FILE__, __FUNCTION__, __LINE__, this->name.c_str());
	}

	mem_dataspace_id = H5Screate_simple(1, mem_dims, nullptr);

	if (H5Dread(samples_all_dataset_id, native_string_datatype_id, mem_dataspace_id, file_dataspace_id, H5P_DEFAULT, samples_buffer) < 0) {
		throw HVCFReadException(__FILE__, __FUNCTION__, __LINE__, this->name.c_str());
	}

	for (unsigned int i = 0; i < file_dims[0]; ++i) {
		samples.emplace_back(samples_buffer[i]);
	}

	H5Dvlen_reclaim(native_string_datatype_id, mem_dataspace_id, H5P_DEFAULT, samples_buffer);

	H5Sclose(mem_dataspace_id);
	H5Sclose(file_dataspace_id);

	return samples;
}

void HVCF::open(const string& name) throw (HVCFOpenException) {
	this->name = name;

//	if ((file_id = H5Fopen(this->name.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT)) < 0) {
//		throw HVCFOpenException(__FILE__, __FUNCTION__, __LINE__, name.c_str());
//	}

	file_id.set(H5Fopen(this->name.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT));
	if (file_id.get() < 0) {
		throw HVCFOpenException(__FILE__, __FUNCTION__, __LINE__, name.c_str());
	}

	if ((samples_group_id = H5Gopen(file_id.get(), SAMPLES_GROUP, H5P_DEFAULT)) < 0) {
		throw HVCFOpenException(__FILE__, __FUNCTION__, __LINE__, name.c_str());
	}

	if ((variants_group_id = H5Gopen(file_id.get(), VARIANTS_GROUP, H5P_DEFAULT)) < 0) {
		throw HVCFOpenException(__FILE__, __FUNCTION__, __LINE__, name.c_str());
	}

	if ((haplotypes_group_id = H5Gopen(file_id.get(), HAPLOTYPES_GROUP, H5P_DEFAULT)) < 0) {
		throw HVCFOpenException(__FILE__, __FUNCTION__, __LINE__, name.c_str());
	}

	if ((samples_all_dataset_id = H5Dopen(samples_group_id, SAMPLES_ALL_DATASET, H5P_DEFAULT)) < 0) {
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

	if (H5Dclose(samples_all_dataset_id) < 0) {
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

	if (H5Gclose(samples_group_id) < 0) {
		throw HVCFCloseException(__FILE__, __FUNCTION__, __LINE__, name.c_str());
	}

//	if (H5Fclose(file_id) < 0) {
//		throw HVCFCloseException(__FILE__, __FUNCTION__, __LINE__, name.c_str());
//	}
	file_id.close();

	haplotypes_group_id = numeric_limits<hid_t>::min();
	variants_group_id = numeric_limits<hid_t>::min();
	samples_group_id = numeric_limits<hid_t>::min();
//	file_id = numeric_limits<hid_t>::min();

	name.clear();
}

}
