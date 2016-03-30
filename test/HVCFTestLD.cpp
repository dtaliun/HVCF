#include <array>
#include <gtest/gtest.h>
#include <cmath>
#include <chrono>
#include "../src/include/HVCF.h"

class HVCFTestLD : public::testing::Test {
protected:
	unordered_map<unsigned long long int, unordered_map<unsigned long long int, double>> precomputed_ld;

	virtual ~HVCFTestLD() {
	}

	virtual void SetUp() {
		sph_umich_edu::GzipReader reader;
		regex header_regex("^CHR[[:space:]]+POS1[[:space:]]+POS2[[:space:]]+N_CHR[[:space:]]+R\\^2[[:space:]]+D[[:space:]]+Dprime$");
		regex separator_regex("[[:space:]]+");
		const cregex_token_iterator end;
		unsigned int i = 0u;
		unsigned long long int position1 = 0ul;
		unsigned long long int position2 = 0ul;
		double rsquare = 0.0;

		auto precomputed_ld_it = precomputed_ld.end();

		reader.set_file_name("vcftools_ld.tbl.gz");
		reader.open();

		char* line = reader.get_line();

		if (reader.read_line() >= 0) {
			if (!regex_match(line, header_regex)) {
				return;
			}
		}

		while (reader.read_line() >= 0) {
			cregex_token_iterator fields_iter(line, line + strlen(line), separator_regex, -1);
			i = 0u;
			position1 = position2 = numeric_limits<unsigned long long int>::min();
			rsquare = numeric_limits<double>::quiet_NaN();
			while (fields_iter != end) {
				switch (i) {
				case 1:
					position1 = stoull(fields_iter->str(), nullptr, 10);
					break;
				case 2:
					position2 = stoull(fields_iter->str(), nullptr, 10);
					break;
				case 4:
					rsquare = stod(fields_iter->str(), nullptr);
					break;
				default:
					break;
				}
				++fields_iter;
				++i;
			}

			precomputed_ld_it = precomputed_ld.emplace(position1, unordered_map<unsigned long long int, double>()).first;
			precomputed_ld_it->second.emplace(position2, rsquare);

			precomputed_ld_it = precomputed_ld.emplace(position2, unordered_map<unsigned long long int, double>()).first;
			precomputed_ld_it->second.emplace(position1, rsquare);
		}

		reader.close();
	}

	virtual void TearDown() {
	}
};

TEST_F(HVCFTestLD, LD) {
	std::chrono::time_point<std::chrono::system_clock> start, end;
	std::chrono::duration<double> elapsed_seconds;
	sph_umich_edu::HVCF hvcf;
	sph_umich_edu::VCFReader vcf;

	// BEGIN: create test HVCF file.
	ASSERT_EQ(0u, hvcf.get_n_opened_objects());
	ASSERT_EQ(0u, sph_umich_edu::HVCF::get_n_all_opened_objects());

	vcf.open("1000G_phase3.EUR.chr20.LD_test.vcf.gz");
	hvcf.create("test_ld.h5");

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

	hvcf.close();
	ASSERT_EQ(0u, hvcf.get_n_opened_objects());
	ASSERT_EQ(0u, sph_umich_edu::HVCF::get_n_all_opened_objects());
	// END: create test HVCF file.

	hvcf.open("test_ld.h5");
	ASSERT_EQ(4u, hvcf.get_n_opened_objects());
	ASSERT_EQ(4u, sph_umich_edu::HVCF::get_n_all_opened_objects());

	start = std::chrono::system_clock::now();
	hvcf.compute_ld("20", 11650214ul, 60759931ul);
	end = std::chrono::system_clock::now();
	elapsed_seconds = end - start;
	GTEST_LOG_(INFO) << "Elapsed time = " << elapsed_seconds.count() << " sec";
	ASSERT_EQ(4u, hvcf.get_n_opened_objects());

	start = std::chrono::system_clock::now();
	hvcf.compute_ld("20", "20:11650214_G/A", 14403183ul, 55378791ul);
	end = std::chrono::system_clock::now();
	elapsed_seconds = end - start;
	GTEST_LOG_(INFO) << "Elapsed time = " << elapsed_seconds.count() << " sec";
	ASSERT_EQ(4u, hvcf.get_n_opened_objects());

	start = std::chrono::system_clock::now();
	hvcf.compute_ld("20", "20:60759931_C/T", 14403183ul, 55378791ul);
	end = std::chrono::system_clock::now();
	elapsed_seconds = end - start;
	GTEST_LOG_(INFO) << "Elapsed time = " << elapsed_seconds.count() << " sec";
	ASSERT_EQ(4u, hvcf.get_n_opened_objects());

	start = std::chrono::system_clock::now();
	hvcf.compute_ld("20", "20:46211051_A/G", 14403183ul, 55378791ul);
	end = std::chrono::system_clock::now();
	elapsed_seconds = end - start;
	GTEST_LOG_(INFO) << "Elapsed time = " << elapsed_seconds.count() << " sec";
	ASSERT_EQ(4u, hvcf.get_n_opened_objects());

//	cout << "precomputed ld:" << precomputed_ld.at(11650214ul).at(55378791ul) << endl;

	hvcf.close();
	ASSERT_EQ(0u, hvcf.get_n_opened_objects());
	ASSERT_EQ(0u, sph_umich_edu::HVCF::get_n_all_opened_objects());

}
