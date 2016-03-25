#include "include/HVCF.h"

namespace sph_umich_edu {

constexpr char HVCF::SAMPLES_GROUP[];
constexpr char HVCF::VARIANTS_GROUP[];
constexpr char HVCF::SAMPLES_ALL_DATASET[];
constexpr char HVCF::HAPLOTYPES_DATASET[];
constexpr char HVCF::VARIANT_NAMES_DATASET[];

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

hid_t HVCF::create_haplotypes_2D_dataset(const string& name, hid_t group_id, hsize_t variants_chunk_size, hsize_t samples_chunk_size) throw (HVCFWriteException) {
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

	if (name.length() == 0u) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Empty dataset name.");
	}

	if ((dataspace_id = H5Screate_simple(2, initial_dims, maximum_dims)) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while creating dataspace.");
	}

	if ((dataset_property_id = H5Pcreate(H5P_DATASET_CREATE)) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while creating dataset property.");
	}

	if ((H5Pset_chunk(dataset_property_id, 2, chunk_dims) < 0) || (H5Pset_deflate(dataset_property_id, 9) < 0)) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while setting dataset properties.");
	}

	if ((dataset_id = H5Dcreate(group_id, name.c_str(), H5T_NATIVE_UCHAR, dataspace_id, H5P_DEFAULT, dataset_property_id, H5P_DEFAULT)) < 0) {
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

	if ((group_id = H5Gcreate(variants_group_id, name.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT)) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while creating group.");
	}

	dataset_id = create_haplotypes_2D_dataset(HAPLOTYPES_DATASET, group_id, 100000, 10000);
	dataset_id.close();

	dataset_id = create_strings_1D_dataset(VARIANT_NAMES_DATASET, group_id, 100000);
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
		throw HVCFOpenException(__FILE__, __FUNCTION__, __LINE__, "Error while opening dataset.");
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

void HVCF::write_names(hid_t group_id, char* const* buffer, unsigned int n_variants) throw (HVCFWriteException) {
	HDF5DatasetIdentifier dataset_id;
	HDF5DataspaceIdentifier file_dataspace_id;
	HDF5DataspaceIdentifier memory_dataspace_id;

	hsize_t mem_dims[1]{n_variants};
	hsize_t file_dims[1]{0};
	hsize_t file_offset[1]{0};

	if ((dataset_id = H5Dopen(group_id, VARIANT_NAMES_DATASET, H5P_DEFAULT)) < 0) {
		throw HVCFOpenException(__FILE__, __FUNCTION__, __LINE__, "Error while opening dataset.");
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

	if (H5Dwrite(dataset_id, native_string_datatype_id, memory_dataspace_id, file_dataspace_id, H5P_DEFAULT, buffer) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while writing to dataset.");
	}
}

void HVCF::create(const string& name) throw (HVCFWriteException) {
	this->name = name;

	if ((file_id = H5Fcreate(this->name.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT)) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while creating file.");
	}

	if ((samples_group_id = H5Gcreate(file_id, SAMPLES_GROUP, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT)) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while creating group.");
	}

	if ((variants_group_id = H5Gcreate(file_id, VARIANTS_GROUP, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT)) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while creating group.");
	}

	if (((native_string_datatype_id = H5Tcopy(H5T_C_S1)) < 0) || (H5Tset_size(native_string_datatype_id, H5T_VARIABLE) < 0)) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while creating datatype.");
	}

	if ((samples_all_dataset_id = create_strings_1D_dataset(SAMPLES_ALL_DATASET, samples_group_id, 1000)) < 0) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while creating dataset.");
	}
}

void HVCF::set_samples(const vector<string>& samples) throw (HVCFWriteException) {
	HDF5DataspaceIdentifier file_dataspace_id;
	HDF5DataspaceIdentifier memory_dataspace_id;

	hsize_t mem_dims[1]{samples.size()};
	hsize_t file_dims[1]{0};
	hsize_t file_offset[1]{0};
	const char* buffer[samples.size()];

	for (unsigned int i = 0u; i < samples.size(); ++i) {
		buffer[i] = samples.at(i).c_str();
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

void HVCF::set_population(const string& name, const vector<string>& samples) throw (HVCFWriteException) {
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
	auto buffers_it = buffers.end();

	if (chromosomes.count(chromosome) == 0) {
		chromosomes_it = chromosomes.emplace(chromosome, std::move(unique_ptr<HDF5GroupIdentifier>(new HDF5GroupIdentifier()))).first;
		buffers_it = buffers.emplace(chromosome, std::move(unique_ptr<IOBuffer>(new IOBuffer(999, get_n_samples())))).first;
		chromosomes_it->second->set(create_chromosome_group(chromosome));
	} else {
		chromosomes_it = chromosomes.find(chromosome);
		buffers_it = buffers.find(chromosome);
	}

	if (buffers_it->second->is_full()) {
		//flush!
		write_haplotypes(chromosomes_it->second->get(), buffers_it->second->get_haplotypes_buffer(), buffers_it->second->get_n_variants(), buffers_it->second->get_n_haplotypes());
		write_names(chromosomes_it->second->get(), buffers_it->second->get_names_buffer(), buffers_it->second->get_n_variants());
		buffers_it->second->reset();
	}

	buffers_it->second->add_variant(variant);

//	cout << chromosome << " " << variant.get_pos().get_value() << endl;
}

void HVCF::flush_write_buffers() throw (HVCFWriteException) {
	auto buffer_it = buffers.end();

	for (auto&& entry : chromosomes) {
		buffer_it = buffers.find(entry.first);
		if (!buffer_it->second->is_empty()) {
			write_haplotypes(entry.second->get(), buffer_it->second->get_haplotypes_buffer(), buffer_it->second->get_n_variants(), buffer_it->second->get_n_haplotypes());
			write_names(entry.second->get(), buffer_it->second->get_names_buffer(), buffer_it->second->get_n_variants());
			buffer_it->second->reset();
		}
	}
}

hsize_t HVCF::get_n_samples() throw (HVCFReadException) {
	HDF5DataspaceIdentifier file_dataspace_id;
	hsize_t file_dims[1]{0};

	if ((file_dataspace_id = H5Dget_space(samples_all_dataset_id)) < 0) {
		throw HVCFReadException(__FILE__, __FUNCTION__, __LINE__, "Error while getting dataspace.");
	}

	if (H5Sget_simple_extent_dims(file_dataspace_id, file_dims, nullptr) < 0) {
		throw HVCFReadException(__FILE__, __FUNCTION__, __LINE__, "Error while getting dataspace dimensions.");
	}

	return file_dims[0];
}

vector<string> HVCF::get_samples() throw (HVCFReadException) {
	HDF5DataspaceIdentifier file_dataspace_id;
	hsize_t file_dims[1]{0};

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
	HDF5DatasetIdentifier dataset_id;
	HDF5DataspaceIdentifier file_dataspace_id;
	HDF5DataspaceIdentifier mem_dataspace_id;

	vector<string> samples;
	hsize_t file_dims[1]{0};
	hsize_t mem_dims[1]{0};

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

	if ((variants_group_id = H5Gopen(file_id, VARIANTS_GROUP, H5P_DEFAULT)) < 0) {
		throw HVCFOpenException(__FILE__, __FUNCTION__, __LINE__, "Error while opening group.");
	}

	if ((samples_all_dataset_id = H5Dopen(samples_group_id, SAMPLES_ALL_DATASET, H5P_DEFAULT)) < 0) {
		throw HVCFOpenException(__FILE__, __FUNCTION__, __LINE__, "Error while opening dataset.");
	}

	if (((native_string_datatype_id = H5Tcopy(H5T_C_S1)) < 0) || (H5Tset_size(native_string_datatype_id, H5T_VARIABLE) < 0)) {
		throw HVCFWriteException(__FILE__, __FUNCTION__, __LINE__, "Error while creating datatype.");
	}
}

void HVCF::close() throw (HVCFCloseException) {
	native_string_datatype_id.close();

	for (auto&& entry : chromosomes) {
		entry.second->close();
	}

	samples_all_dataset_id.close();
	variants_group_id.close();
	samples_group_id.close();
	file_id.close();
	name.clear();
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
