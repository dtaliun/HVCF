#include <array>
#include <gtest/gtest.h>
#include <cmath>
#include <chrono>
#include "../src/include/HVCF.h"

class HVCFTestReadWrite : public::testing::Test {
protected:
	virtual ~HVCFTestReadWrite() {
	}

	virtual void SetUp() {

	}

	virtual void TearDown() {
	}
};

TEST_F(HVCFTestReadWrite, Create) {
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
		ASSERT_EQ(3u, hvcf.get_n_opened_objects());
		ASSERT_EQ(3u, sph_umich_edu::HVCF::get_n_all_opened_objects());
		hvcf.set_samples(in_samples);
		ASSERT_EQ(3u, hvcf.get_n_opened_objects());
		ASSERT_EQ(3u, sph_umich_edu::HVCF::get_n_all_opened_objects());
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
	ASSERT_EQ(3u, hvcf.get_n_opened_objects());
	ASSERT_EQ(3u, sph_umich_edu::HVCF::get_n_all_opened_objects());
	hvcf.set_samples(in_samples);
	ASSERT_EQ(5u, hvcf.get_n_samples());
	ASSERT_EQ(3u, hvcf.get_n_opened_objects());
	hvcf.set_population("POP1", in_pop1_samples);
	ASSERT_EQ(3u, hvcf.get_n_opened_objects());
	hvcf.close();
	ASSERT_EQ(0u, hvcf.get_n_opened_objects());
	ASSERT_EQ(0u, sph_umich_edu::HVCF::get_n_all_opened_objects());

	hvcf.open("test.h5");
	ASSERT_EQ(3u, hvcf.get_n_opened_objects());
	ASSERT_EQ(3u, sph_umich_edu::HVCF::get_n_all_opened_objects());
	vector<string> out_samples = std::move(hvcf.get_samples());
	ASSERT_EQ(3u, hvcf.get_n_opened_objects());
	vector<string> out_pop1_samples_ordered = std::move(hvcf.get_population("POP1"));
	ASSERT_EQ(3u, hvcf.get_n_opened_objects());
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

TEST_F(HVCFTestReadWrite, WriteVCF) {
	{ // 'dummy' scope to check if HVCF object closes every opened HDF5 identifier on its destruction
		sph_umich_edu::HVCF hvcf;
		sph_umich_edu::VCFReader vcf;

		ASSERT_EQ(0u, hvcf.get_n_opened_objects());
		ASSERT_EQ(0u, sph_umich_edu::HVCF::get_n_all_opened_objects());

		vcf.open("1000G_phase3.EUR.chr20-22.10K.vcf.gz");
		hvcf.create("test.h5");

		hvcf.set_samples(std::move(vcf.get_variant().get_samples()));

		ASSERT_EQ(503u, hvcf.get_n_samples());

		ASSERT_EQ(3u, hvcf.get_n_opened_objects());
		ASSERT_EQ(3u, sph_umich_edu::HVCF::get_n_all_opened_objects());

		while (vcf.read_next_variant()) {
			hvcf.write_variant(vcf.get_variant());
		}
		hvcf.flush_write_buffer();

		hvcf.create_indices();
		ASSERT_EQ(6u, hvcf.get_n_opened_objects());

		vcf.close();

		ASSERT_EQ(6u, hvcf.get_n_opened_objects());
		ASSERT_EQ(6u, sph_umich_edu::HVCF::get_n_all_opened_objects());
	}

	ASSERT_EQ(0u, sph_umich_edu::HVCF::get_n_all_opened_objects());

	sph_umich_edu::HVCF hvcf;

	ASSERT_EQ(0u, hvcf.get_n_opened_objects());
	ASSERT_EQ(0u, sph_umich_edu::HVCF::get_n_all_opened_objects());

	hvcf.open("test.h5");
	ASSERT_EQ(6u, hvcf.get_n_opened_objects());
	ASSERT_EQ(6u, sph_umich_edu::HVCF::get_n_all_opened_objects());

	ASSERT_EQ(0u, hvcf.get_n_variants("19"));
	ASSERT_EQ(6u, hvcf.get_n_opened_objects());

	ASSERT_EQ(9930u, hvcf.get_n_variants("20"));
	ASSERT_EQ(6u, hvcf.get_n_opened_objects());

	ASSERT_EQ(9941u, hvcf.get_n_variants("21"));
	ASSERT_EQ(6u, hvcf.get_n_opened_objects());

	ASSERT_EQ(9959u, hvcf.get_n_variants("22"));
	ASSERT_EQ(6u, hvcf.get_n_opened_objects());

	ASSERT_EQ(29830u, hvcf.get_n_variants());
	ASSERT_EQ(6u, hvcf.get_n_opened_objects());

	hvcf.close();
	ASSERT_EQ(0u, hvcf.get_n_opened_objects());
	ASSERT_EQ(0u, sph_umich_edu::HVCF::get_n_all_opened_objects());
}

TEST_F(HVCFTestReadWrite, VariantLookupByPosition) {
	std::chrono::time_point<std::chrono::system_clock> start, end;
	std::chrono::duration<double> elapsed_seconds;
	sph_umich_edu::HVCF hvcf;

	ASSERT_EQ(0u, hvcf.get_n_opened_objects());
	ASSERT_EQ(0u, sph_umich_edu::HVCF::get_n_all_opened_objects());

	hvcf.open("test.h5");
	ASSERT_EQ(6u, hvcf.get_n_opened_objects());
	ASSERT_EQ(6u, sph_umich_edu::HVCF::get_n_all_opened_objects());

	// --
	ASSERT_EQ(-1, hvcf.get_variant_offset_by_position_eq("XYZ", 60343ul));
	ASSERT_EQ(6u, hvcf.get_n_opened_objects());

	ASSERT_EQ(-1, hvcf.get_variant_offset_by_position_ge("XYZ", 60343ul));
	ASSERT_EQ(6u, hvcf.get_n_opened_objects());

	ASSERT_EQ(-1, hvcf.get_variant_offset_by_position_le("XYZ", 60343ul));
	ASSERT_EQ(6u, hvcf.get_n_opened_objects());

	// --
	start = std::chrono::system_clock::now();
	ASSERT_LE(-1, hvcf.get_variant_offset_by_position_eq("20", 0ul));
	end = std::chrono::system_clock::now();
	elapsed_seconds = end - start;
	GTEST_LOG_(INFO) << "Lookup by single position (eq) = " << elapsed_seconds.count() << " sec";
	ASSERT_EQ(6u, hvcf.get_n_opened_objects());

	start = std::chrono::system_clock::now();
	ASSERT_LE(0, hvcf.get_variant_offset_by_position_ge("20", 0ul));
	end = std::chrono::system_clock::now();
	elapsed_seconds = end - start;
	GTEST_LOG_(INFO) << "Lookup by single position (ge) = " << elapsed_seconds.count() << " sec";
	ASSERT_EQ(6u, hvcf.get_n_opened_objects());

	start = std::chrono::system_clock::now();
	ASSERT_LE(-1, hvcf.get_variant_offset_by_position_le("20", 0ul));
	end = std::chrono::system_clock::now();
	elapsed_seconds = end - start;
	GTEST_LOG_(INFO) << "Lookup by single position (le) = " << elapsed_seconds.count() << " sec";
	ASSERT_EQ(6u, hvcf.get_n_opened_objects());

	// --
	start = std::chrono::system_clock::now();
	ASSERT_EQ(0, hvcf.get_variant_offset_by_position_ge("20", 60000ul));
	end = std::chrono::system_clock::now();
	elapsed_seconds = end - start;
	GTEST_LOG_(INFO) << "Lookup by single position (ge) = " << elapsed_seconds.count() << " sec";
	ASSERT_EQ(6u, hvcf.get_n_opened_objects());

	start = std::chrono::system_clock::now();
	ASSERT_EQ(-1, hvcf.get_variant_offset_by_position_le("20", 60000ul));
	end = std::chrono::system_clock::now();
	elapsed_seconds = end - start;
	GTEST_LOG_(INFO) << "Lookup by single position (le) = " << elapsed_seconds.count() << " sec";
	ASSERT_EQ(6u, hvcf.get_n_opened_objects());

	// --
	start = std::chrono::system_clock::now();
	ASSERT_EQ(0, hvcf.get_variant_offset_by_position_eq("20", 60343ul));
	end = std::chrono::system_clock::now();
	elapsed_seconds = end - start;
	GTEST_LOG_(INFO) << "Lookup by single position (eq) = " << elapsed_seconds.count() << " sec";
	ASSERT_EQ(6u, hvcf.get_n_opened_objects());

	start = std::chrono::system_clock::now();
	ASSERT_EQ(0, hvcf.get_variant_offset_by_position_ge("20", 60343ul));
	end = std::chrono::system_clock::now();
	elapsed_seconds = end - start;
	GTEST_LOG_(INFO) << "Lookup by single position (ge) = " << elapsed_seconds.count() << " sec";
	ASSERT_EQ(6u, hvcf.get_n_opened_objects());

	start = std::chrono::system_clock::now();
	ASSERT_EQ(0, hvcf.get_variant_offset_by_position_le("20", 60343ul));
	end = std::chrono::system_clock::now();
	elapsed_seconds = end - start;
	GTEST_LOG_(INFO) << "Lookup by single position (le) = " << elapsed_seconds.count() << " sec";
	ASSERT_EQ(6u, hvcf.get_n_opened_objects());

	// --
	start = std::chrono::system_clock::now();
	ASSERT_EQ(6202, hvcf.get_variant_offset_by_position_eq("20", 263529ul));
	end = std::chrono::system_clock::now();
	elapsed_seconds = end - start;
	GTEST_LOG_(INFO) << "Lookup by single position (eq) = " << elapsed_seconds.count() << " sec";
	ASSERT_EQ(6u, hvcf.get_n_opened_objects());

	start = std::chrono::system_clock::now();
	ASSERT_EQ(6202, hvcf.get_variant_offset_by_position_ge("20", 263529ul));
	end = std::chrono::system_clock::now();
	elapsed_seconds = end - start;
	GTEST_LOG_(INFO) << "Lookup by single position (ge) = " << elapsed_seconds.count() << " sec";
	ASSERT_EQ(6u, hvcf.get_n_opened_objects());

	start = std::chrono::system_clock::now();
	ASSERT_EQ(6202, hvcf.get_variant_offset_by_position_le("20", 263529ul));
	end = std::chrono::system_clock::now();
	elapsed_seconds = end - start;
	GTEST_LOG_(INFO) << "Lookup by single position (le) = " << elapsed_seconds.count() << " sec";
	ASSERT_EQ(6u, hvcf.get_n_opened_objects());

	// --
	start = std::chrono::system_clock::now();
	ASSERT_EQ(9929, hvcf.get_variant_offset_by_position_eq("20", 372328ul));
	end = std::chrono::system_clock::now();
	elapsed_seconds = end - start;
	GTEST_LOG_(INFO) << "Lookup by single position (eq) = " << elapsed_seconds.count() << " sec";
	ASSERT_EQ(6u, hvcf.get_n_opened_objects());

	start = std::chrono::system_clock::now();
	ASSERT_EQ(9929, hvcf.get_variant_offset_by_position_ge("20", 372328ul));
	end = std::chrono::system_clock::now();
	elapsed_seconds = end - start;
	GTEST_LOG_(INFO) << "Lookup by single position (ge) = " << elapsed_seconds.count() << " sec";
	ASSERT_EQ(6u, hvcf.get_n_opened_objects());

	start = std::chrono::system_clock::now();
	ASSERT_EQ(9929, hvcf.get_variant_offset_by_position_le("20", 372328ul));
	end = std::chrono::system_clock::now();
	elapsed_seconds = end - start;
	GTEST_LOG_(INFO) << "Lookup by single position (le) = " << elapsed_seconds.count() << " sec";
	ASSERT_EQ(6u, hvcf.get_n_opened_objects());

	// --
	start = std::chrono::system_clock::now();
	ASSERT_EQ(-1, hvcf.get_variant_offset_by_position_eq("20", 372004ul));
	end = std::chrono::system_clock::now();
	elapsed_seconds = end - start;
	GTEST_LOG_(INFO) << "Lookup by single position (eq) = " << elapsed_seconds.count() << " sec";
	ASSERT_EQ(6u, hvcf.get_n_opened_objects());

	start = std::chrono::system_clock::now();
	ASSERT_EQ(9915, hvcf.get_variant_offset_by_position_ge("20", 372004ul));
	end = std::chrono::system_clock::now();
	elapsed_seconds = end - start;
	GTEST_LOG_(INFO) << "Lookup by single position (ge) = " << elapsed_seconds.count() << " sec";
	ASSERT_EQ(6u, hvcf.get_n_opened_objects());

	start = std::chrono::system_clock::now();
	ASSERT_EQ(9914, hvcf.get_variant_offset_by_position_le("20", 372004ul));
	end = std::chrono::system_clock::now();
	elapsed_seconds = end - start;
	GTEST_LOG_(INFO) << "Lookup by single position (le) = " << elapsed_seconds.count() << " sec";
	ASSERT_EQ(6u, hvcf.get_n_opened_objects());

	// --
	start = std::chrono::system_clock::now();
	ASSERT_EQ(-1, hvcf.get_variant_offset_by_position_eq("20", 348210ul));
	end = std::chrono::system_clock::now();
	elapsed_seconds = end - start;
	GTEST_LOG_(INFO) << "Lookup by single position (eq) = " << elapsed_seconds.count() << " sec";
	ASSERT_EQ(6u, hvcf.get_n_opened_objects());

	start = std::chrono::system_clock::now();
	ASSERT_EQ(9000, hvcf.get_variant_offset_by_position_ge("20", 348210ul));
	end = std::chrono::system_clock::now();
	elapsed_seconds = end - start;
	GTEST_LOG_(INFO) << "Lookup by single position (ge) = " << elapsed_seconds.count() << " sec";
	ASSERT_EQ(6u, hvcf.get_n_opened_objects());

	start = std::chrono::system_clock::now();
	ASSERT_EQ(8999, hvcf.get_variant_offset_by_position_le("20", 348210ul));
	end = std::chrono::system_clock::now();
	elapsed_seconds = end - start;
	GTEST_LOG_(INFO) << "Lookup by single position (le) = " << elapsed_seconds.count() << " sec";
	ASSERT_EQ(6u, hvcf.get_n_opened_objects());

	// --
	start = std::chrono::system_clock::now();
	ASSERT_EQ(-1, hvcf.get_variant_offset_by_position_eq("20", 3372004ul));
	end = std::chrono::system_clock::now();
	elapsed_seconds = end - start;
	GTEST_LOG_(INFO) << "Lookup by single position (eq) = " << elapsed_seconds.count() << " sec";
	ASSERT_EQ(6u, hvcf.get_n_opened_objects());

	start = std::chrono::system_clock::now();
	ASSERT_EQ(-1, hvcf.get_variant_offset_by_position_ge("20", 3372004ul));
	end = std::chrono::system_clock::now();
	elapsed_seconds = end - start;
	GTEST_LOG_(INFO) << "Lookup by single position (ge) = " << elapsed_seconds.count() << " sec";
	ASSERT_EQ(6u, hvcf.get_n_opened_objects());

	start = std::chrono::system_clock::now();
	ASSERT_EQ(9929, hvcf.get_variant_offset_by_position_le("20", 3372004ul));
	end = std::chrono::system_clock::now();
	elapsed_seconds = end - start;
	GTEST_LOG_(INFO) << "Lookup by single position (le) = " << elapsed_seconds.count() << " sec";
	ASSERT_EQ(6u, hvcf.get_n_opened_objects());

	hvcf.close();
	ASSERT_EQ(0u, hvcf.get_n_opened_objects());
	ASSERT_EQ(0u, sph_umich_edu::HVCF::get_n_all_opened_objects());
}

TEST_F(HVCFTestReadWrite, VariantLookupByName) {
	std::chrono::time_point<std::chrono::system_clock> start, end;
	std::chrono::duration<double> elapsed_seconds;
	sph_umich_edu::HVCF hvcf;

	ASSERT_EQ(0u, hvcf.get_n_opened_objects());
	ASSERT_EQ(0u, sph_umich_edu::HVCF::get_n_all_opened_objects());

	hvcf.open("test.h5");
	ASSERT_EQ(6u, hvcf.get_n_opened_objects());
	ASSERT_EQ(6u, sph_umich_edu::HVCF::get_n_all_opened_objects());

	ASSERT_EQ(-1, hvcf.get_variant_offset_by_name("XYZ", "20:60343_G/A"));
	ASSERT_EQ(6u, hvcf.get_n_opened_objects());

	start = std::chrono::system_clock::now();
	ASSERT_EQ(-1, hvcf.get_variant_offset_by_name("20", "20:282263_T/A"));
	end = std::chrono::system_clock::now();
	elapsed_seconds = end - start;
	GTEST_LOG_(INFO) << "Lookup by single name = " << elapsed_seconds.count() << " sec";
	ASSERT_EQ(6u, hvcf.get_n_opened_objects());

	start = std::chrono::system_clock::now();
	ASSERT_EQ(0, hvcf.get_variant_offset_by_name("20", "20:60343_G/A"));
	end = std::chrono::system_clock::now();
	elapsed_seconds = end - start;
	GTEST_LOG_(INFO) << "Lookup by single name = " << elapsed_seconds.count() << " sec";
	ASSERT_EQ(6u, hvcf.get_n_opened_objects());

	start = std::chrono::system_clock::now();
	ASSERT_EQ(6828, hvcf.get_variant_offset_by_name("20", "20:282218_C/CATGCAAGGCCCT"));
	end = std::chrono::system_clock::now();
	elapsed_seconds = end - start;
	GTEST_LOG_(INFO) << "Lookup by single name = " << elapsed_seconds.count() << " sec";
	ASSERT_EQ(6u, hvcf.get_n_opened_objects());

	start = std::chrono::system_clock::now();
	ASSERT_EQ(9929, hvcf.get_variant_offset_by_name("20", "20:372328_G/A"));
	end = std::chrono::system_clock::now();
	elapsed_seconds = end - start;
	GTEST_LOG_(INFO) << "Lookup by single name = " << elapsed_seconds.count() << " sec";
	ASSERT_EQ(6u, hvcf.get_n_opened_objects());

	hvcf.close();
	ASSERT_EQ(0u, hvcf.get_n_opened_objects());
	ASSERT_EQ(0u, sph_umich_edu::HVCF::get_n_all_opened_objects());
}

TEST_F(HVCFTestReadWrite, LargeFileTest) {
	std::chrono::time_point<std::chrono::system_clock> start, end;
	std::chrono::duration<double> elapsed_seconds;
	sph_umich_edu::HVCF hvcf;
	sph_umich_edu::VCFReader vcf;

	ASSERT_EQ(0u, hvcf.get_n_opened_objects());
	ASSERT_EQ(0u, sph_umich_edu::HVCF::get_n_all_opened_objects());

	vcf.open("1000G_phase3.EUR.chr20.vcf.gz");
	hvcf.create("test_large.h5");

	hvcf.set_samples(std::move(vcf.get_variant().get_samples()));

	ASSERT_EQ(503u, hvcf.get_n_samples());

	ASSERT_EQ(3u, hvcf.get_n_opened_objects());
	ASSERT_EQ(3u, sph_umich_edu::HVCF::get_n_all_opened_objects());

	while (vcf.read_next_variant()) {
		hvcf.write_variant(vcf.get_variant());
	}
	hvcf.flush_write_buffer();

	hvcf.create_indices();
	ASSERT_EQ(4u, hvcf.get_n_opened_objects());

	vcf.close();

	ASSERT_EQ(4u, hvcf.get_n_opened_objects());
	ASSERT_EQ(4u, sph_umich_edu::HVCF::get_n_all_opened_objects());

	hvcf.close();
	ASSERT_EQ(0u, hvcf.get_n_opened_objects());
	ASSERT_EQ(0u, sph_umich_edu::HVCF::get_n_all_opened_objects());

	hvcf.open("test_large.h5");
	ASSERT_EQ(4u, hvcf.get_n_opened_objects());
	ASSERT_EQ(4u, sph_umich_edu::HVCF::get_n_all_opened_objects());

	// --
	ASSERT_EQ(-1, hvcf.get_variant_offset_by_position_eq("XYZ", 60343ul));
	ASSERT_EQ(4u, hvcf.get_n_opened_objects());

	ASSERT_EQ(-1, hvcf.get_variant_offset_by_position_ge("XYZ", 60343ul));
	ASSERT_EQ(4u, hvcf.get_n_opened_objects());

	ASSERT_EQ(-1, hvcf.get_variant_offset_by_position_le("XYZ", 60343ul));
	ASSERT_EQ(4u, hvcf.get_n_opened_objects());

	// --
	start = std::chrono::system_clock::now();
	ASSERT_LE(-1, hvcf.get_variant_offset_by_position_eq("20", 0ul));
	end = std::chrono::system_clock::now();
	elapsed_seconds = end - start;
	GTEST_LOG_(INFO) << "Lookup by single position (eq) = " << elapsed_seconds.count() << " sec";
	ASSERT_EQ(4u, hvcf.get_n_opened_objects());

	start = std::chrono::system_clock::now();
	ASSERT_LE(0, hvcf.get_variant_offset_by_position_ge("20", 0ul));
	end = std::chrono::system_clock::now();
	elapsed_seconds = end - start;
	GTEST_LOG_(INFO) << "Lookup by single position (ge) = " << elapsed_seconds.count() << " sec";
	ASSERT_EQ(4u, hvcf.get_n_opened_objects());

	start = std::chrono::system_clock::now();
	ASSERT_LE(-1, hvcf.get_variant_offset_by_position_le("20", 0ul));
	end = std::chrono::system_clock::now();
	elapsed_seconds = end - start;
	GTEST_LOG_(INFO) << "Lookup by single position (le) = " << elapsed_seconds.count() << " sec";
	ASSERT_EQ(4u, hvcf.get_n_opened_objects());

	// --
	start = std::chrono::system_clock::now();
	ASSERT_EQ(0, hvcf.get_variant_offset_by_position_eq("20", 60795ul));
	end = std::chrono::system_clock::now();
	elapsed_seconds = end - start;
	GTEST_LOG_(INFO) << "Lookup by single position (eq) = " << elapsed_seconds.count() << " sec";
	ASSERT_EQ(4u, hvcf.get_n_opened_objects());

	start = std::chrono::system_clock::now();
	ASSERT_EQ(0, hvcf.get_variant_offset_by_position_ge("20", 60795ul));
	end = std::chrono::system_clock::now();
	elapsed_seconds = end - start;
	GTEST_LOG_(INFO) << "Lookup by single position (ge) = " << elapsed_seconds.count() << " sec";
	ASSERT_EQ(4u, hvcf.get_n_opened_objects());

	start = std::chrono::system_clock::now();
	ASSERT_EQ(0, hvcf.get_variant_offset_by_position_le("20", 60795ul));
	end = std::chrono::system_clock::now();
	elapsed_seconds = end - start;
	GTEST_LOG_(INFO) << "Lookup by single position (le) = " << elapsed_seconds.count() << " sec";
	ASSERT_EQ(4u, hvcf.get_n_opened_objects());

	// --
	start = std::chrono::system_clock::now();
	ASSERT_EQ(256157, hvcf.get_variant_offset_by_position_eq("20", 32959324ul));
	end = std::chrono::system_clock::now();
	elapsed_seconds = end - start;
	GTEST_LOG_(INFO) << "Lookup by single position (eq) = " << elapsed_seconds.count() << " sec";
	ASSERT_EQ(4u, hvcf.get_n_opened_objects());

	start = std::chrono::system_clock::now();
	ASSERT_EQ(256157, hvcf.get_variant_offset_by_position_ge("20", 32959324ul));
	end = std::chrono::system_clock::now();
	elapsed_seconds = end - start;
	GTEST_LOG_(INFO) << "Lookup by single position (eq) = " << elapsed_seconds.count() << " sec";
	ASSERT_EQ(4u, hvcf.get_n_opened_objects());

	start = std::chrono::system_clock::now();
	ASSERT_EQ(256157, hvcf.get_variant_offset_by_position_le("20", 32959324ul));
	end = std::chrono::system_clock::now();
	elapsed_seconds = end - start;
	GTEST_LOG_(INFO) << "Lookup by single position (eq) = " << elapsed_seconds.count() << " sec";
	ASSERT_EQ(4u, hvcf.get_n_opened_objects());

	// --
	start = std::chrono::system_clock::now();
	ASSERT_EQ(520853, hvcf.get_variant_offset_by_position_eq("20", 62963628ul));
	end = std::chrono::system_clock::now();
	elapsed_seconds = end - start;
	GTEST_LOG_(INFO) << "Lookup by single position (eq) = " << elapsed_seconds.count() << " sec";
	ASSERT_EQ(4u, hvcf.get_n_opened_objects());

	start = std::chrono::system_clock::now();
	ASSERT_EQ(520853, hvcf.get_variant_offset_by_position_ge("20", 62963628ul));
	end = std::chrono::system_clock::now();
	elapsed_seconds = end - start;
	GTEST_LOG_(INFO) << "Lookup by single position (ge) = " << elapsed_seconds.count() << " sec";
	ASSERT_EQ(4u, hvcf.get_n_opened_objects());

	start = std::chrono::system_clock::now();
	ASSERT_EQ(520853, hvcf.get_variant_offset_by_position_le("20", 62963628ul));
	end = std::chrono::system_clock::now();
	elapsed_seconds = end - start;
	GTEST_LOG_(INFO) << "Lookup by single position (le) = " << elapsed_seconds.count() << " sec";
	ASSERT_EQ(4u, hvcf.get_n_opened_objects());

	// --
	start = std::chrono::system_clock::now();
	ASSERT_EQ(-1, hvcf.get_variant_offset_by_position_eq("20", 167311ul));
	end = std::chrono::system_clock::now();
	elapsed_seconds = end - start;
	GTEST_LOG_(INFO) << "Lookup by single position (eq) = " << elapsed_seconds.count() << " sec";
	ASSERT_EQ(4u, hvcf.get_n_opened_objects());

	start = std::chrono::system_clock::now();
	ASSERT_EQ(1000, hvcf.get_variant_offset_by_position_ge("20", 167311ul));
	end = std::chrono::system_clock::now();
	elapsed_seconds = end - start;
	GTEST_LOG_(INFO) << "Lookup by single position (ge) = " << elapsed_seconds.count() << " sec";
	ASSERT_EQ(4u, hvcf.get_n_opened_objects());

	start = std::chrono::system_clock::now();
	ASSERT_EQ(999, hvcf.get_variant_offset_by_position_le("20", 167311ul));
	end = std::chrono::system_clock::now();
	elapsed_seconds = end - start;
	GTEST_LOG_(INFO) << "Lookup by single position (ge) = " << elapsed_seconds.count() << " sec";
	ASSERT_EQ(4u, hvcf.get_n_opened_objects());

	//--
	start = std::chrono::system_clock::now();
	ASSERT_EQ(-1, hvcf.get_variant_offset_by_position_eq("20", 3372004ul));
	end = std::chrono::system_clock::now();
	elapsed_seconds = end - start;
	GTEST_LOG_(INFO) << "Lookup by single position (eq) = " << elapsed_seconds.count() << " sec";
	ASSERT_EQ(4u, hvcf.get_n_opened_objects());

	start = std::chrono::system_clock::now();
	ASSERT_EQ(30449, hvcf.get_variant_offset_by_position_ge("20", 3372004ul));
	end = std::chrono::system_clock::now();
	elapsed_seconds = end - start;
	GTEST_LOG_(INFO) << "Lookup by single position (ge) = " << elapsed_seconds.count() << " sec";
	ASSERT_EQ(4u, hvcf.get_n_opened_objects());

	start = std::chrono::system_clock::now();
	ASSERT_EQ(30448, hvcf.get_variant_offset_by_position_le("20", 3372004ul));
	end = std::chrono::system_clock::now();
	elapsed_seconds = end - start;
	GTEST_LOG_(INFO) << "Lookup by single position (le) = " << elapsed_seconds.count() << " sec";
	ASSERT_EQ(4u, hvcf.get_n_opened_objects());

	// --
	ASSERT_EQ(-1, hvcf.get_variant_offset_by_name("XYZ", "20:60343_G/A"));
	ASSERT_EQ(4u, hvcf.get_n_opened_objects());

	start = std::chrono::system_clock::now();
	ASSERT_EQ(-1, hvcf.get_variant_offset_by_name("20", "20:282263_T/A"));
	end = std::chrono::system_clock::now();
	elapsed_seconds = end - start;
	GTEST_LOG_(INFO) << "Lookup by single name = " << elapsed_seconds.count() << " sec";
	ASSERT_EQ(4u, hvcf.get_n_opened_objects());

	start = std::chrono::system_clock::now();
	ASSERT_EQ(0, hvcf.get_variant_offset_by_name("20", "20:60795_G/C"));
	end = std::chrono::system_clock::now();
	elapsed_seconds = end - start;
	GTEST_LOG_(INFO) << "Lookup by single name = " << elapsed_seconds.count() << " sec";
	ASSERT_EQ(4u, hvcf.get_n_opened_objects());

	start = std::chrono::system_clock::now();
	ASSERT_EQ(264834, hvcf.get_variant_offset_by_name("20", "20:34244739_A/ACTTTTAAATATATGGGCTTTTAAATATAAGCCCATCTCTACTAAAATTTATGGG"));
	end = std::chrono::system_clock::now();
	elapsed_seconds = end - start;
	GTEST_LOG_(INFO) << "Lookup by single name = " << elapsed_seconds.count() << " sec";
	ASSERT_EQ(4u, hvcf.get_n_opened_objects());

	start = std::chrono::system_clock::now();
	ASSERT_EQ(520853, hvcf.get_variant_offset_by_name("20", "20:62963628_C/G"));
	end = std::chrono::system_clock::now();
	elapsed_seconds = end - start;
	GTEST_LOG_(INFO) << "Lookup by single name = " << elapsed_seconds.count() << " sec";
	ASSERT_EQ(4u, hvcf.get_n_opened_objects());

	hvcf.close();
	ASSERT_EQ(0u, hvcf.get_n_opened_objects());
	ASSERT_EQ(0u, sph_umich_edu::HVCF::get_n_all_opened_objects());
}

TEST_F(HVCFTestReadWrite, ChunkReadTest) {
	std::chrono::time_point<std::chrono::system_clock> start, end;
	std::chrono::duration<double> elapsed_seconds;
	sph_umich_edu::HVCF hvcf;

	hvcf.open("test_large.h5");
	ASSERT_EQ(4u, hvcf.get_n_opened_objects());
	ASSERT_EQ(4u, sph_umich_edu::HVCF::get_n_all_opened_objects());

	start = std::chrono::system_clock::now();
	hvcf.chunk_read_test("20", 61955ul, 167900ul);
	end = std::chrono::system_clock::now();
	elapsed_seconds = end - start;
	GTEST_LOG_(INFO) << "1000 variants = " << elapsed_seconds.count() << " sec";
	ASSERT_EQ(4u, hvcf.get_n_opened_objects());

	start = std::chrono::system_clock::now();
	hvcf.chunk_read_test("20", 3322215ul, 4397713ul);
	end = std::chrono::system_clock::now();
	elapsed_seconds = end - start;
	GTEST_LOG_(INFO) << "10000 variants = " << elapsed_seconds.count() << " sec";
	ASSERT_EQ(4u, hvcf.get_n_opened_objects());

	start = std::chrono::system_clock::now();
	hvcf.chunk_read_test("20", 3322215ul, 14879114ul);
	end = std::chrono::system_clock::now();
	elapsed_seconds = end - start;
	GTEST_LOG_(INFO) << "100000 variants = " << elapsed_seconds.count() << " sec";
	ASSERT_EQ(4u, hvcf.get_n_opened_objects());

	hvcf.close();
	ASSERT_EQ(0u, hvcf.get_n_opened_objects());
	ASSERT_EQ(0u, sph_umich_edu::HVCF::get_n_all_opened_objects());

}

