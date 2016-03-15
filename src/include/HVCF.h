#ifndef SRC_HVCF_H_
#define SRC_HVCF_H_

#include <iostream>
#include <limits>
#include <vector>
#include <memory>
#include "hdf5.h"
#include "HVCFOpenException.h"
#include "HVCFCloseException.h"
#include "HVCFWriteException.h"

using namespace std;

namespace sph_umich_edu {

class HVCF {
private:
	string name;
	hid_t file_id;
	hid_t samples_group_id;
	hid_t populations_group_id;
	hid_t variants_group_id;
	hid_t haplotypes_group_id;
	hid_t sample_names_dataset_id;
	hid_t variant_names_dataset_id;
	hid_t native_string_datatype_id;

	static constexpr char SAMPLES_GROUP[] = "samples";
	static constexpr char POPULATIONS_GROUP[] = "populations";
	static constexpr char VARIANTS_GROUP[] = "variants";
	static constexpr char HAPLOTYPES_GROUP[] = "haplotypes";
	static constexpr char SAMPLE_NAMES_DATASET[] = "names";
	static constexpr char VARIANT_NAMES_DATASET[] = "names";

	hid_t create_strings_dataset(const string& name, hid_t group_id, hsize_t chunk_size);

public:
	HVCF();
	virtual ~HVCF();

	void create(const string& name) throw (HVCFWriteException);
	void write_samples(const string& name, const vector<string>& samples) throw (HVCFWriteException);

	void open(const string& name) throw (HVCFOpenException);


	void close() throw (HVCFCloseException);
};

}

#endif
