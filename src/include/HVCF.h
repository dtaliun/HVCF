#ifndef SRC_HVCF_H_
#define SRC_HVCF_H_

#include <iostream>
#include <limits>
#include <vector>
#include <map>
#include <memory>
#include <iterator>
#include <algorithm>
#include "hdf5.h"
#include "HVCFOpenException.h"
#include "HVCFCloseException.h"
#include "HVCFWriteException.h"
#include "HVCFReadException.h"

using namespace std;

namespace sph_umich_edu {

class HVCF {
private:
	string name;
	hid_t file_id;
	hid_t samples_group_id;
	hid_t variants_group_id;
	hid_t haplotypes_group_id;
	hid_t samples_all_dataset_id;
	hid_t variant_names_dataset_id;
	hid_t native_string_datatype_id;

	static constexpr char SAMPLES_GROUP[] = "samples";
	static constexpr char VARIANTS_GROUP[] = "variants";
	static constexpr char HAPLOTYPES_GROUP[] = "haplotypes";
	static constexpr char SAMPLES_ALL_DATASET[] = "ALL";
	static constexpr char VARIANT_NAMES_DATASET[] = "names";

	hid_t create_strings_1D_dataset(const string& name, hid_t group_id, hsize_t chunk_size) throw (HVCFWriteException);
	hid_t create_ulong_1D_dataset(const string&name, hid_t group_id, hsize_t chunk_size) throw (HVCFWriteException);

public:
	HVCF();
	virtual ~HVCF();

	void create(const string& name) throw (HVCFWriteException);
	void set_samples(const vector<string>& samples) throw (HVCFWriteException);
	void set_population(const string& name, const vector<string>& samples) throw (HVCFWriteException);

	vector<string> get_samples() throw (HVCFReadException);

	void open(const string& name) throw (HVCFOpenException);


	void close() throw (HVCFCloseException);
};

}

#endif
