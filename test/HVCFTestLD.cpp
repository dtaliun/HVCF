#include <array>
#include <gtest/gtest.h>
#include <cmath>
#include <chrono>
#include "../src/include/HVCF.h"

using namespace std;

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

		for (auto&& entry : precomputed_ld) {
			bool all_nan = true;
			for (auto&& entry_entry : entry.second) {
				if (!std::isnan(entry_entry.second)) {
					all_nan = false;
					break;
				}
			}
			if (all_nan) {
				entry.second.emplace(entry.first, numeric_limits<double>::quiet_NaN());
			} else {
				entry.second.emplace(entry.first, 1.0);
			}
		}
	}

	virtual void TearDown() {
	}
};

TEST_F(HVCFTestLD, LD) {
	std::chrono::time_point<std::chrono::system_clock> start, end;
	std::chrono::duration<double> elapsed_seconds;
	sph_umich_edu::HVCF hvcf;
	vector<sph_umich_edu::variants_pair> result;

	// BEGIN: create test HVCF file.
	ASSERT_EQ(0u, hvcf.get_n_opened_objects());
	ASSERT_EQ(0u, sph_umich_edu::HVCF::get_n_all_opened_objects());

	hvcf.create("test_ld.h5");
	hvcf.import_vcf("1000G_phase3.EUR.chr20.LD_test.vcf.gz");

	ASSERT_EQ(503u, hvcf.get_n_samples());
	ASSERT_EQ(4u, hvcf.get_n_opened_objects());
	ASSERT_EQ(4u, sph_umich_edu::HVCF::get_n_all_opened_objects());

	hvcf.close();
	ASSERT_EQ(0u, hvcf.get_n_opened_objects());
	ASSERT_EQ(0u, sph_umich_edu::HVCF::get_n_all_opened_objects());
	// END: create test HVCF file.

	hvcf.open("test_ld.h5");
	ASSERT_EQ(4u, hvcf.get_n_opened_objects());
	ASSERT_EQ(4u, sph_umich_edu::HVCF::get_n_all_opened_objects());

	start = std::chrono::system_clock::now();
	result.clear();
	hvcf.compute_ld("20", "ALL", 11650214ul, 60759931ul, result);
	end = std::chrono::system_clock::now();
	elapsed_seconds = end - start;
	GTEST_LOG_(INFO) << "Elapsed time = " << elapsed_seconds.count() << " sec";
	ASSERT_EQ(81u, result.size());
	for (auto&& pair : result) {
		if (std::isnan(precomputed_ld.at(pair.position1).at(pair.position2))) {
			ASSERT_TRUE(pair.rsquare);
		} else {
			ASSERT_NEAR(pair.rsquare, precomputed_ld.at(pair.position1).at(pair.position2), 0.00000001);
		}
	}
	ASSERT_EQ(4u, hvcf.get_n_opened_objects());

	result.clear();
	start = std::chrono::system_clock::now();
	hvcf.compute_ld("20", "ALL", 11650214ul, 11650214ul, result);
	end = std::chrono::system_clock::now();
	elapsed_seconds = end - start;
	GTEST_LOG_(INFO) << "Elapsed time = " << elapsed_seconds.count() << " sec";
	ASSERT_EQ(1u, result.size());
	for (auto&& pair : result) {
		if (std::isnan(precomputed_ld.at(pair.position1).at(pair.position2))) {
			ASSERT_TRUE(pair.rsquare);
		} else {
			ASSERT_NEAR(pair.rsquare, precomputed_ld.at(pair.position1).at(pair.position2), 0.00000001);
		}
	}
	ASSERT_EQ(4u, hvcf.get_n_opened_objects());

	result.clear();
	start = std::chrono::system_clock::now();
	hvcf.compute_ld("20", "ALL", 160759931ul, 260759931ul, result);
	end = std::chrono::system_clock::now();
	elapsed_seconds = end - start;
	GTEST_LOG_(INFO) << "Elapsed time = " << elapsed_seconds.count() << " sec";
	ASSERT_EQ(0u, result.size());
	for (auto&& pair : result) {
		if (std::isnan(precomputed_ld.at(pair.position1).at(pair.position2))) {
			ASSERT_TRUE(pair.rsquare);
		} else {
			ASSERT_NEAR(pair.rsquare, precomputed_ld.at(pair.position1).at(pair.position2), 0.00000001);
		}
	}
	ASSERT_EQ(4u, hvcf.get_n_opened_objects());

	result.clear();
	start = std::chrono::system_clock::now();
	hvcf.compute_ld("20", "ALL", "20:11650214_G/A", 14403183ul, 55378791ul, result);
	end = std::chrono::system_clock::now();
	elapsed_seconds = end - start;
	GTEST_LOG_(INFO) << "Elapsed time = " << elapsed_seconds.count() << " sec";
	ASSERT_EQ(7u, result.size());
	for (auto&& pair : result) {
		if (std::isnan(precomputed_ld.at(pair.position1).at(pair.position2))) {
			ASSERT_TRUE(pair.rsquare);
		} else {
			ASSERT_NEAR(pair.rsquare, precomputed_ld.at(pair.position1).at(pair.position2), 0.00000001);
		}
	}
	ASSERT_EQ(4u, hvcf.get_n_opened_objects());

	result.clear();
	start = std::chrono::system_clock::now();
	hvcf.compute_ld("20", "ALL", "20:60759931_C/T", 14403183ul, 55378791ul, result);
	end = std::chrono::system_clock::now();
	elapsed_seconds = end - start;
	GTEST_LOG_(INFO) << "Elapsed time = " << elapsed_seconds.count() << " sec";
	ASSERT_EQ(7u, result.size());
	for (auto&& pair : result) {
		if (std::isnan(precomputed_ld.at(pair.position1).at(pair.position2))) {
			ASSERT_TRUE(pair.rsquare);
		} else {
			ASSERT_NEAR(pair.rsquare, precomputed_ld.at(pair.position1).at(pair.position2), 0.00000001);
		}
	}
	ASSERT_EQ(4u, hvcf.get_n_opened_objects());

	result.clear();
	start = std::chrono::system_clock::now();
	hvcf.compute_ld("20", "ALL", "20:46211051_A/G", 14403183ul, 55378791ul, result);
	end = std::chrono::system_clock::now();
	elapsed_seconds = end - start;
	GTEST_LOG_(INFO) << "Elapsed time = " << elapsed_seconds.count() << " sec";
	ASSERT_EQ(6u, result.size());
	for (auto&& pair : result) {
		if (std::isnan(precomputed_ld.at(pair.position1).at(pair.position2))) {
			ASSERT_TRUE(pair.rsquare);
		} else {
			ASSERT_NEAR(pair.rsquare, precomputed_ld.at(pair.position1).at(pair.position2), 0.00000001);
		}
	}
	ASSERT_EQ(4u, hvcf.get_n_opened_objects());

	result.clear();
	start = std::chrono::system_clock::now();
	hvcf.compute_ld("20", "ALL", "20:46211051_A/G", 19485821ul, 19485821ul, result);
	end = std::chrono::system_clock::now();
	elapsed_seconds = end - start;
	GTEST_LOG_(INFO) << "Elapsed time = " << elapsed_seconds.count() << " sec";
	ASSERT_EQ(1u, result.size());
	for (auto&& pair : result) {
		if (std::isnan(precomputed_ld.at(pair.position1).at(pair.position2))) {
			ASSERT_TRUE(pair.rsquare);
		} else {
			ASSERT_NEAR(pair.rsquare, precomputed_ld.at(pair.position1).at(pair.position2), 0.00000001);
		}
	}
	ASSERT_EQ(4u, hvcf.get_n_opened_objects());

	result.clear();
	start = std::chrono::system_clock::now();
	hvcf.compute_ld("20", "ALL", "20:19485821_A/G", 19485821ul, 19485821ul, result);
	end = std::chrono::system_clock::now();
	elapsed_seconds = end - start;
	GTEST_LOG_(INFO) << "Elapsed time = " << elapsed_seconds.count() << " sec";
	ASSERT_EQ(0u, result.size());
	for (auto&& pair : result) {
		if (std::isnan(precomputed_ld.at(pair.position1).at(pair.position2))) {
			ASSERT_TRUE(pair.rsquare);
		} else {
			ASSERT_NEAR(pair.rsquare, precomputed_ld.at(pair.position1).at(pair.position2), 0.00000001);
		}
	}
	ASSERT_EQ(4u, hvcf.get_n_opened_objects());

	result.clear();
	start = std::chrono::system_clock::now();
	hvcf.compute_ld("20", "ALL", "20:46211051_A/G", 16600000ul, 50000000ul, result);
	end = std::chrono::system_clock::now();
	elapsed_seconds = end - start;
	GTEST_LOG_(INFO) << "Elapsed time = " << elapsed_seconds.count() << " sec";
	ASSERT_EQ(3u, result.size());
	for (auto&& pair : result) {
		if (std::isnan(precomputed_ld.at(pair.position1).at(pair.position2))) {
			ASSERT_TRUE(pair.rsquare);
		} else {
			ASSERT_NEAR(pair.rsquare, precomputed_ld.at(pair.position1).at(pair.position2), 0.00000001);
		}
	}
	ASSERT_EQ(4u, hvcf.get_n_opened_objects());

	result.clear();
	start = std::chrono::system_clock::now();
	hvcf.compute_ld("20", "ALL", "20:46211051_A/G", 100ul, 16600000ul, result);
	end = std::chrono::system_clock::now();
	elapsed_seconds = end - start;
	GTEST_LOG_(INFO) << "Elapsed time = " << elapsed_seconds.count() << " sec";
	ASSERT_EQ(2u, result.size());
	for (auto&& pair : result) {
		if (std::isnan(precomputed_ld.at(pair.position1).at(pair.position2))) {
			ASSERT_TRUE(pair.rsquare);
		} else {
			ASSERT_NEAR(pair.rsquare, precomputed_ld.at(pair.position1).at(pair.position2), 0.00000001);
		}
	}
	ASSERT_EQ(4u, hvcf.get_n_opened_objects());

	result.clear();
	start = std::chrono::system_clock::now();
	hvcf.compute_ld("20", "ALL", "20:46211051_A/G", 50000000ul, 166000000ul, result);
	end = std::chrono::system_clock::now();
	elapsed_seconds = end - start;
	GTEST_LOG_(INFO) << "Elapsed time = " << elapsed_seconds.count() << " sec";
	ASSERT_EQ(3u, result.size());
	for (auto&& pair : result) {
		if (std::isnan(precomputed_ld.at(pair.position1).at(pair.position2))) {
			ASSERT_TRUE(pair.rsquare);
		} else {
			ASSERT_NEAR(pair.rsquare, precomputed_ld.at(pair.position1).at(pair.position2), 0.00000001);
		}
	}
	ASSERT_EQ(4u, hvcf.get_n_opened_objects());

	result.clear();
	start = std::chrono::system_clock::now();
	hvcf.compute_ld("20", "ALL", "20:46211051_A/G", 166000000ul, 166000000ul, result);
	end = std::chrono::system_clock::now();
	elapsed_seconds = end - start;
	GTEST_LOG_(INFO) << "Elapsed time = " << elapsed_seconds.count() << " sec";
	ASSERT_EQ(0u, result.size());
	for (auto&& pair : result) {
		if (std::isnan(precomputed_ld.at(pair.position1).at(pair.position2))) {
			ASSERT_TRUE(pair.rsquare);
		} else {
			ASSERT_NEAR(pair.rsquare, precomputed_ld.at(pair.position1).at(pair.position2), 0.00000001);
		}
	}
	ASSERT_EQ(4u, hvcf.get_n_opened_objects());

	hvcf.close();
	ASSERT_EQ(0u, hvcf.get_n_opened_objects());
	ASSERT_EQ(0u, sph_umich_edu::HVCF::get_n_all_opened_objects());

}

TEST_F(HVCFTestLD, LargeLD) {
	std::chrono::time_point<std::chrono::system_clock> start, end;
	std::chrono::duration<double> elapsed_seconds;
	sph_umich_edu::HVCF hvcf;
	vector<sph_umich_edu::variants_pair> result;

	ASSERT_EQ(0u, hvcf.get_n_opened_objects());
	ASSERT_EQ(0u, sph_umich_edu::HVCF::get_n_all_opened_objects());

	hvcf.open("test_large.h5");
	ASSERT_EQ(4u, hvcf.get_n_opened_objects());
	ASSERT_EQ(4u, sph_umich_edu::HVCF::get_n_all_opened_objects());

	result.clear();
	start = std::chrono::system_clock::now();
	hvcf.compute_ld("20", "ALL", 61955ul, 167870ul, result);
	end = std::chrono::system_clock::now();
	elapsed_seconds = end - start;
	GTEST_LOG_(INFO) << "1000 variants pairwise LD = " << elapsed_seconds.count() << " sec";
	ASSERT_EQ(1000u * 1000u, result.size());
	ASSERT_EQ(4u, hvcf.get_n_opened_objects());

	result.clear();
	start = std::chrono::system_clock::now();
	hvcf.compute_ld("20", "ALL", "20:61955_C/T", 61955ul, 167900ul, result);
	end = std::chrono::system_clock::now();
	elapsed_seconds = end - start;
	GTEST_LOG_(INFO) << "1000 variants vs 1 LD = " << elapsed_seconds.count() << " sec";
	ASSERT_EQ(1000u, result.size());
	ASSERT_EQ(4u, hvcf.get_n_opened_objects());

	result.clear();
	start = std::chrono::system_clock::now();
	hvcf.compute_ld("20", "ALL", 3322215ul, 4397641ul, result);
	end = std::chrono::system_clock::now();
	elapsed_seconds = end - start;
	GTEST_LOG_(INFO) << "10000 variants pairwise LD = " << elapsed_seconds.count() << " sec";
	ASSERT_EQ(10000u * 10000u, result.size());
	ASSERT_EQ(4u, hvcf.get_n_opened_objects());

	result.clear();
	start = std::chrono::system_clock::now();
	hvcf.compute_ld("20", "ALL", "20:3514537_T/C", 3322215ul, 4397713ul, result);
	end = std::chrono::system_clock::now();
	elapsed_seconds = end - start;
	GTEST_LOG_(INFO) << "10000 variants vs 1 LD = " << elapsed_seconds.count() << " sec";
	ASSERT_EQ(10000u, result.size());
	ASSERT_EQ(4u, hvcf.get_n_opened_objects());

	result.clear();
	start = std::chrono::system_clock::now();
	hvcf.compute_ld("20", "ALL", "20:3514537_T/C", 3322215ul, 14879114ul, result);
	end = std::chrono::system_clock::now();
	elapsed_seconds = end - start;
	GTEST_LOG_(INFO) << "100000 variants vs 1 LD = " << elapsed_seconds.count() << " sec";
	ASSERT_EQ(100000u, result.size());
	ASSERT_EQ(4u, hvcf.get_n_opened_objects());

	hvcf.close();
	ASSERT_EQ(0u, hvcf.get_n_opened_objects());
	ASSERT_EQ(0u, sph_umich_edu::HVCF::get_n_all_opened_objects());
}
