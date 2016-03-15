#include <array>
#include <gtest/gtest.h>
#include <cmath>
#include "../src/include/HVCF.h"

class HVCFTest : public::testing::Test {
protected:
	virtual ~HVCFTest() {
	}

	virtual void SetUp() {
	}

	virtual void TearDown() {
	}
};

TEST_F(HVCFTest, Create) {
	sph_umich_edu::HVCF hvcf;

	vector<string> samples{"sample1", "sample2", "sample3", "sample4", "sample5"};

	hvcf.create("test.h5");
	hvcf.set_samples(samples);
	hvcf.close();

	hvcf.open("test.h5");
	hvcf.close();
}
