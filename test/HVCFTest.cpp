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
	vector<string> in_samples{"sample1", "sample2", "sample3", "sample4", "sample5"};
	vector<string> in_pop1_samples{"sample5", "sample2", "sample3"};
	vector<string> in_pop1_samples_ordered{"sample2", "sample3", "sample5"};

	// test if all handlers are closed on exceptions
	bool exception = false;
	try {
		sph_umich_edu::HVCF hvcf;
		ASSERT_EQ(0u, hvcf.get_n_opened_objects());
		ASSERT_EQ(0u, sph_umich_edu::HVCF::get_n_all_opened_objects());
		hvcf.create("test.h5");
		ASSERT_EQ(6u, hvcf.get_n_opened_objects());
		ASSERT_EQ(7u, sph_umich_edu::HVCF::get_n_all_opened_objects());
		hvcf.set_samples(in_samples);
		ASSERT_EQ(6u, hvcf.get_n_opened_objects());
		ASSERT_EQ(7u, sph_umich_edu::HVCF::get_n_all_opened_objects());
		hvcf.set_population("", in_pop1_samples);
	} catch (sph_umich_edu::HVCFWriteException &e) {
		exception = true;
	}
	ASSERT_TRUE(exception);
	ASSERT_EQ(0u, sph_umich_edu::HVCF::get_n_all_opened_objects());

	sph_umich_edu::HVCF hvcf;

	ASSERT_EQ(0u, hvcf.get_n_opened_objects());
	ASSERT_EQ(0u, sph_umich_edu::HVCF::get_n_all_opened_objects());
	hvcf.create("test.h5");
	ASSERT_EQ(6u, hvcf.get_n_opened_objects());
	ASSERT_EQ(7u, sph_umich_edu::HVCF::get_n_all_opened_objects());
	hvcf.set_samples(in_samples);
	ASSERT_EQ(5u, hvcf.get_n_samples());
	ASSERT_EQ(6u, hvcf.get_n_opened_objects());
	hvcf.set_population("POP1", in_pop1_samples);
	ASSERT_EQ(6u, hvcf.get_n_opened_objects());
	hvcf.close();
	ASSERT_EQ(0u, hvcf.get_n_opened_objects());
	ASSERT_EQ(0u, sph_umich_edu::HVCF::get_n_all_opened_objects());

	hvcf.open("test.h5");
	ASSERT_EQ(6u, hvcf.get_n_opened_objects());
	ASSERT_EQ(7u, sph_umich_edu::HVCF::get_n_all_opened_objects());
	vector<string> out_samples = std::move(hvcf.get_samples());
	ASSERT_EQ(6u, hvcf.get_n_opened_objects());
	vector<string> out_pop1_samples_ordered = std::move(hvcf.get_population("POP1"));
	ASSERT_EQ(6u, hvcf.get_n_opened_objects());
	hvcf.close();
	ASSERT_EQ(0u, hvcf.get_n_opened_objects());
	ASSERT_EQ(0u, sph_umich_edu::HVCF::get_n_all_opened_objects());

	ASSERT_EQ(in_samples.size(), out_samples.size());
	for (unsigned int i = 0u; i < in_samples.size(); ++i) {
		ASSERT_EQ(in_samples[i], out_samples[i]);
	}

	ASSERT_EQ(in_pop1_samples_ordered.size(), out_pop1_samples_ordered.size());
	for (unsigned int i = 0u; i < in_pop1_samples_ordered.size(); ++i) {
		ASSERT_EQ(in_pop1_samples_ordered[i], out_pop1_samples_ordered[i]);
	}

}

TEST_F(HVCFTest, WriteVCF) {
	{ // 'dummy' scope to check if HVCF object closes every opened HDF5 identifier on its destruction
		sph_umich_edu::HVCF hvcf;
		sph_umich_edu::VCFReader vcf;

		ASSERT_EQ(0u, hvcf.get_n_opened_objects());
		ASSERT_EQ(0u, sph_umich_edu::HVCF::get_n_all_opened_objects());

		vcf.open("1000G_phase3.EUR.chr20-22.10K.vcf.gz");
		hvcf.create("test.h5");

		hvcf.set_samples(std::move(vcf.get_variant().get_samples()));

		ASSERT_EQ(503u, hvcf.get_n_samples());

		ASSERT_EQ(6u, hvcf.get_n_opened_objects());
		ASSERT_EQ(7u, sph_umich_edu::HVCF::get_n_all_opened_objects());

		while (vcf.read_next_variant()) {
			hvcf.write_variant(vcf.get_variant());
		}
		vcf.close();

		ASSERT_EQ(9u, hvcf.get_n_opened_objects());
		ASSERT_EQ(10u, sph_umich_edu::HVCF::get_n_all_opened_objects());
	}

	ASSERT_EQ(0u, sph_umich_edu::HVCF::get_n_all_opened_objects());
}
