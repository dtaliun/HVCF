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

	vector<string> in_samples{"sample1", "sample2", "sample3", "sample4", "sample5"};
	vector<string> in_pop1_samples{"sample5", "sample2", "sample3"};

	hvcf.create("test.h5");
	hvcf.set_samples(in_samples);
	hvcf.set_population("POP1", in_pop1_samples);
	hvcf.close();

	hvcf.open("test.h5");
	vector<string> out_samples = std::move(hvcf.get_samples());
	hvcf.close();

	ASSERT_EQ(in_samples.size(), out_samples.size());
	for (unsigned int i = 0u; i < in_samples.size(); ++i) {
		ASSERT_EQ(in_samples.at(i), out_samples.at(i));
	}

}
