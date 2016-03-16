#ifndef SRC_HVCF_H_
#define SRC_HVCF_H_

#include <iostream>
#include <limits>
#include <vector>
#include <map>
#include <memory>
#include <iterator>
#include <algorithm>
#include <chrono>
#include "hdf5.h"
#include "HVCFOpenException.h"
#include "HVCFCloseException.h"
#include "HVCFWriteException.h"
#include "HVCFReadException.h"
#include "HDF5FileIdentifier.h"
#include "HDF5GroupIdentifier.h"
#include "HDF5DatasetIdentifier.h"
#include "HDF5DatatypeIdentifier.h"
#include "HDF5DataspaceIdentifier.h"
#include "HDF5PropertyIdentifier.h"

using namespace std;

namespace sph_umich_edu {

class HVCF {
private:
	string name;

	HDF5FileIdentifier file_id;
	HDF5GroupIdentifier samples_group_id;
	HDF5GroupIdentifier variants_group_id;
	HDF5GroupIdentifier haplotypes_group_id;
	HDF5DatasetIdentifier samples_all_dataset_id;
	HDF5DatasetIdentifier variant_names_dataset_id;
	HDF5DatatypeIdentifier native_string_datatype_id;

	static constexpr char SAMPLES_GROUP[] = "samples";
	static constexpr char VARIANTS_GROUP[] = "variants";
	static constexpr char HAPLOTYPES_GROUP[] = "haplotypes";
	static constexpr char SAMPLES_ALL_DATASET[] = "ALL";
	static constexpr char VARIANT_NAMES_DATASET[] = "names";

	hid_t create_strings_1D_dataset(const string& name, hid_t group_id, hsize_t chunk_size) throw (HVCFWriteException);
	hid_t create_hsize_1D_dataset(const string&name, hid_t group_id, hsize_t chunk_size) throw (HVCFWriteException);

public:
	HVCF();
	virtual ~HVCF();

	void create(const string& name) throw (HVCFWriteException);
	void set_samples(const vector<string>& samples) throw (HVCFWriteException);
	void set_population(const string& name, const vector<string>& samples) throw (HVCFWriteException);

	vector<string> get_samples() throw (HVCFReadException);
	vector<string> get_population(const string& name) throw (HVCFReadException);

	void open(const string& name) throw (HVCFOpenException);

	void close() throw (HVCFCloseException);

	unsigned int get_n_opened_objects() const;
	static unsigned int get_n_all_opened_objects();
};

}

#endif
