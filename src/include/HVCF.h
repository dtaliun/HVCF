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

	void open_datatypes() throw (std::exception);
	void close_datatypes() throw (std::exception);

	static constexpr char SAMPLES[] = "samples";
	static constexpr char POPULATIONS[] = "populations";
	static constexpr char VARIANTS[] = "variants";
	static constexpr char HAPLOTYPES[] = "haplotypes";

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
