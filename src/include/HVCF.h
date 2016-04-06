#ifndef SRC_HVCF_H_
#define SRC_HVCF_H_

#include <iostream>
#include <limits>
#include <vector>
#include <map>
#include <set>
#include <unordered_map>
#include <memory>
#include <iterator>
#include <algorithm>
#include <cstdlib>
#include <chrono>

#define ARMA_NO_DEBUG
#include <armadillo>

#include "hdf5.h"

#include "Types.h"
#include "HVCFOpenException.h"
#include "HVCFCloseException.h"
#include "HVCFWriteException.h"
#include "HVCFReadException.h"
#include "HVCFCreateException.h"
#include "HDF5FileIdentifier.h"
#include "HDF5GroupIdentifier.h"
#include "HDF5DatasetIdentifier.h"
#include "HDF5DatatypeIdentifier.h"
#include "HDF5DataspaceIdentifier.h"
#include "HDF5PropertyIdentifier.h"
#include "HVCFConfiguration.h"
#include "../../../auxc/MiniVCF/src/include/VCFReader.h"
#include "WriteBuffer.h"
#include "../blosc/blosc_filter.h"

using namespace std;
using namespace arma;

namespace sph_umich_edu {

class HVCF {
private:
	string name;

	HDF5FileIdentifier file_id;
	HDF5GroupIdentifier samples_group_id;
	HDF5GroupIdentifier chromosomes_group_id;

	unsigned int N_HASH_BUCKETS;
	unsigned int MAX_VARIANTS_IN_INTERVAL_BUCKET;
	unsigned int VARIANTS_CHUNK_SIZE;
	unsigned int SAMPLES_CHUNK_SIZE;
	const char* COMPRESSION;
	unsigned int COMPRESSION_LEVEL;

	static constexpr char CHROMOSOMES_GROUP[] = "chromosomes";
	static constexpr char SAMPLES_GROUP[] = "samples";
	static constexpr char VARIANTS_DATASET[] = "variants";
	static constexpr char HAPLOTYPES_DATASET[] = "haplotypes";
	static constexpr char SAMPLE_NAMES_DATASET[] = "names";
	static constexpr char SAMPLE_SUBSETS_DATASET[] = "subsets";

	static constexpr char VARIABLE_LENGTH_STRING_TYPE[] = "variable_length_string_type";
	static constexpr char VARIANTS_ENTRY_TYPE[] = "variants_entry_type";
	static constexpr char SUBSETS_ENTRY_TYPE[] = "subsets_entry_type";
	static constexpr char STRING_INDEX_ENTRY_TYPE[] = "string_index_entry_type";
	static constexpr char INTERVAL_INDEX_ENTRY_TYPE[] = "interval_index_entry_type";
	static constexpr char HASH_INDEX_ENTRY_TYPE[] = "hash_index_entry_type";
	static constexpr char ULL_INDEX_ENTRY_TYPE[] = "ull_index_entry_type";
	static constexpr char VARIANT_NAMES_INDEX_GROUP[] = "names_index";
	static constexpr char VARIANT_INTERVALS_INDEX_GROUP[] = "intervals_index";
	static constexpr char VARIANT_INTERVALS_INDEX[] = "intervals";
	static constexpr char VARIANT_HASH_INDEX[] = "hashes";
	static constexpr char INDEX_BUCKETS[] = "buckets";

	unordered_map<string, unique_ptr<HDF5GroupIdentifier>> chromosomes;
	unordered_map<string, unique_ptr<WriteBuffer>> write_buffers;

	hid_t create_variants_entry_memory_datatype() throw (HVCFCreateException);
	hid_t create_subsets_entry_memory_datatype() throw (HVCFCreateException);
	hid_t create_ull_index_entry_memory_datatype() throw (HVCFCreateException);
	hid_t create_string_index_entry_memory_datatype() throw (HVCFCreateException);
	hid_t create_interval_index_entry_memory_datatype() throw (HVCFCreateException);
	hid_t create_hash_index_entry_memory_datatype() throw (HVCFCreateException);

	hid_t create_sample_names_dataset(hid_t group_id, hsize_t chunk_size) throw (HVCFWriteException);
	hid_t create_sample_subsets_dataset(hid_t group_id, hsize_t chunk_size) throw (HVCFWriteException);
	hid_t create_haplotypes_dataset(hid_t group_id, hsize_t variants_chunk_size, hsize_t samples_chunk_size) throw (HVCFWriteException);
	hid_t create_variants_dataset(hid_t group_id, hsize_t chunk_size) throw (HVCFWriteException);
	hid_t create_chromosome_group(const string& name) throw (HVCFWriteException);

	void initialize_ull_index_buckets(hid_t chromosome_group_id, const char* index_group_name) throw (HVCFWriteException);
	void initialize_string_index_buckets(hid_t chromosome_group_id, const char* index_group_name) throw (HVCFWriteException);

	void write_intervals_index_bucket(hid_t chromosome_group_id, const vector<hsize_t>& offsets, interval_index_entry_type& interval_index_entry) throw (HVCFWriteException);
	void write_intervals_index(hid_t chromosome_group_id, const interval_index_entry_type* interval_index_entries, unsigned int n_interval_index_entries) throw (HVCFWriteException);

	void write_names_index_bucket(hid_t chromosome_group_id, const vector<hsize_t>& offsets, hash_index_entry_type& hash_index_entry) throw (HVCFWriteException);
	void write_names_index(hid_t chromosome_group_id, const hash_index_entry_type* hash_index_entries, unsigned int n_hash_index_entries) throw (HVCFWriteException);

	void write_haplotypes(hid_t group_id, const unsigned char* buffer, unsigned int n_variants, unsigned int n_haplotypes) throw (HVCFWriteException);
	void write_variants(hid_t group_id, const variants_entry_type* buffer, unsigned int n_variants) throw (HVCFWriteException);

	void create_indices(hid_t chromosome_group_id) throw (HVCFWriteException);

public:
	HVCF();
	HVCF(const HVCFConfiguration& configuration);
	virtual ~HVCF() noexcept;

	void create(const string& name) throw (HVCFWriteException);

	void set_samples(const vector<string>& samples) throw (HVCFWriteException);
	void create_sample_subset(const string& name, const vector<string>& samples) throw (HVCFWriteException);

	void write_variant(const Variant& variant) throw (HVCFWriteException);
	void flush_write_buffer() throw (HVCFWriteException);
	void create_indices() throw (HVCFWriteException);

	void open(const string& name) throw (HVCFOpenException);
	void close() throw (HVCFCloseException);

	hsize_t get_n_samples() throw (HVCFReadException);
	vector<string> get_samples() throw (HVCFReadException);
	unsigned int get_n_sample_subsets() throw (HVCFReadException);
	vector<string> get_sample_subsets() throw (HVCFReadException);
	unsigned int get_n_samples_in_subset(const string& name) throw (HVCFReadException);
	vector<string> get_samples_in_subset(const string& name) throw (HVCFReadException);

	hsize_t get_n_variants() throw (HVCFReadException);
	hsize_t get_n_variants_in_chromosome(const string& chromosome) throw (HVCFReadException);

	long long int get_variant_offset_by_position_eq(const string& chromosome, unsigned long long int position) throw (HVCFReadException);
	long long int get_variant_offset_by_position_ge(const string& chromosome, unsigned long long int position) throw (HVCFReadException);
	long long int get_variant_offset_by_position_le(const string& chromosome, unsigned long long int position) throw (HVCFReadException);
	long long int get_variant_offset_by_name(const string& chromosome, const string& name) throw (HVCFReadException);

	void compute_ld(const string& chromosome, unsigned long long int start_position, unsigned long long end_position, vector<variants_pair>& result) throw (HVCFReadException);
	void compute_ld(const string& chromosome, const string& lead_variant_name, unsigned long long int start_position, unsigned long long end_position, vector<variants_pair>& result) throw (HVCFReadException);

	unsigned int get_n_opened_objects() const;
	static unsigned int get_n_all_opened_objects();

	void chunk_read_test(const string& chromosome, unsigned long long int start_position, unsigned long long end_position) throw (HVCFReadException);

	void chunk_read_test2() throw (HVCFReadException);
};

}

#endif
