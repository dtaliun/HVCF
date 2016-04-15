#include <boost/python.hpp>
#include <boost/python/def.hpp>
#include <boost/python/return_by_value.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#include <boost/python/exception_translator.hpp>
#include <python2.7/Python.h>

#include "../src/include/HVCF.h"

using namespace sph_umich_edu;
using namespace boost::python;

void translator(const HVCFException& e) {
	PyErr_SetString(PyExc_UserWarning, e.what());
}

void (HVCF::*compute_region_ld)(const string& chromosome, const string& subset, unsigned long long int start_position, unsigned long long end_position, vector<ld_query_result>& result) = &HVCF::compute_ld;
void (HVCF::*compute_lead_ld)(const string& chromosome, const string& subset, const string& lead_variant_name, unsigned long long int start_position, unsigned long long end_position, vector<ld_query_result>& result) = &HVCF::compute_ld;
void (HVCF::*extract_haplotypes_for_variant)(const string& subset, const string& chromosome, const string& variant_name, vector<variant_haplotypes_query_result>& result) = &HVCF::extract_haplotypes;
void (HVCF::*extract_haplotypes_for_sample)(const string& sample, const string& chromosome, unsigned long long int start, unsigned long long int end, vector<sample_haplotypes_query_result>& result) = &HVCF::extract_haplotypes;

BOOST_PYTHON_MODULE(PyHVCF)
{
	register_exception_translator<HVCFException>(&translator);

	def("get_n_all_opened_objects", HVCF::get_n_all_opened_objects);

	class_<variant_query_result>("VariantQueryResult", init<const char*, const char*, const char*, unsigned long long int>())
			.def_readonly("name", &variant_query_result::name)
			.def_readonly("ref", &variant_query_result::ref)
			.def_readonly("alt", &variant_query_result::alt)
			.def_readonly("position", &variant_query_result::position)
		;

	class_<frequency_query_result>("FrequencyQueryResult", init<const char*, const char*, const char*, unsigned long long int, double, double>())
			.def_readonly("name", &frequency_query_result::name)
			.def_readonly("ref", &frequency_query_result::ref)
			.def_readonly("alt", &frequency_query_result::alt)
			.def_readonly("position", &frequency_query_result::position)
			.def_readonly("ref_af", &frequency_query_result::ref_af)
			.def_readonly("alt_af", &frequency_query_result::alt_af)
		;

	class_<ld_query_result>("LDQueryResult", init<const char*, unsigned long int, const char*, unsigned long int , double, double>())
			.def_readonly("name1", &LDQueryResult::name1)
			.def_readonly("position1", &LDQueryResult::position1)
			.def_readonly("name2", &LDQueryResult::name2)
			.def_readonly("position2", &LDQueryResult::position2)
			.def_readonly("r", &LDQueryResult::r)
			.def_readonly("rsquare", &LDQueryResult::rsquare)
		;

	class_<variant_haplotypes_query_result>("VariantHaplotypesQueryResult")
			.def_readonly("sample", &VariantHaplotypesQueryResult::sample)
			.def_readonly("allele1", &VariantHaplotypesQueryResult::allele1)
			.def_readonly("allele2", &VariantHaplotypesQueryResult::allele2)
		;

	class_<sample_haplotypes_query_result>("SampleHaplotypesQueryResult")
			.def_readonly("name", &SampleHaplotypesQueryResult::name)
			.def_readonly("position", &SampleHaplotypesQueryResult::position)
			.def_readonly("allele1", &SampleHaplotypesQueryResult::allele1)
			.def_readonly("allele2", &SampleHaplotypesQueryResult::allele2)
		;

	class_<vector<string>>("NamesVector")
			.def(vector_indexing_suite<std::vector<std::string>>())
		;

	class_<vector<double>>("DoublesVector")
			.def(vector_indexing_suite<std::vector<double>>())
		;

	class_<vector<variant_query_result>>("Variants")
			.def(vector_indexing_suite<std::vector<variant_query_result>>())
		;

	class_<vector<frequency_query_result>>("Frequencies")
			.def(vector_indexing_suite<std::vector<frequency_query_result>>())
		;

	class_<vector<ld_query_result>>("Pairs")
			.def(vector_indexing_suite<std::vector<ld_query_result>>())
		;

	class_<vector<variant_haplotypes_query_result>>("VariantHaplotypes")
			.def(vector_indexing_suite<std::vector<variant_haplotypes_query_result>>())
		;

	class_<vector<sample_haplotypes_query_result>>("SampleHaplotypes")
			.def(vector_indexing_suite<std::vector<sample_haplotypes_query_result>>())
		;

	class_<HVCF, boost::noncopyable>("HVCF")
			.def("create", &HVCF::create)
			.def("open", &HVCF::open)
			.def("close", &HVCF::close)
			.def("import_vcf", &HVCF::import_vcf)
			.def("create_sample_subset", &HVCF::create_sample_subset)
			.def("get_n_samples", &HVCF::get_n_samples)
			.def("get_samples", &HVCF::get_samples, return_value_policy<return_by_value>())
			.def("get_n_sample_subsets", &HVCF::get_n_sample_subsets)
			.def("get_sample_subsets", &HVCF::get_sample_subsets, return_value_policy<return_by_value>())
			.def("get_n_samples_in_subset", &HVCF::get_n_samples_in_subset)
			.def("get_samples_in_subset", &HVCF::get_samples_in_subset, return_value_policy<return_by_value>())
			.def("get_chromosomes", &HVCF::get_chromosomes, return_value_policy<return_by_value>())
			.def("has_chromosome", &HVCF::has_chromosome)
			.def("get_chromosome_start", &HVCF::get_chromosome_start)
			.def("get_chromosome_end", &HVCF::get_chromosome_end)
			.def("get_n_variants", &HVCF::get_n_variants)
			.def("get_n_variants_in_chromosome", &HVCF::get_n_variants_in_chromosome)
			.def("compute_ld", compute_region_ld)
			.def("compute_ld", compute_lead_ld)
			.def("compute_frequencies", &HVCF::compute_frequencies)
			.def("extract_variants", &HVCF::extract_variants)
			.def("extract_haplotypes", extract_haplotypes_for_variant)
			.def("extract_haplotypes", extract_haplotypes_for_sample)
			.def("get_n_opened_objects", &HVCF::get_n_opened_objects)
		;
}
