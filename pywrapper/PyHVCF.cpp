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

void (HVCF::*compute_region_ld)(const string& chromosome, const string& subset, unsigned long long int start_position, unsigned long long end_position, vector<variants_pair>& result) = &HVCF::compute_ld;
void (HVCF::*compute_lead_ld)(const string& chromosome, const string& subset, const string& lead_variant_name, unsigned long long int start_position, unsigned long long end_position, vector<variants_pair>& result) = &HVCF::compute_ld;
void (HVCF::*extract_haplotypes_for_variant)(const string& chromosome, const string& subset, const string& variant_name, vector<unsigned char>& result) = &HVCF::extract_haplotypes;
void (HVCF::*extract_haplotypes_for_sample)(const string& chromosome, const string& sample, unsigned long long int start, unsigned long long int end, vector<unsigned char>& result) = &HVCF::extract_haplotypes;

BOOST_PYTHON_MODULE(PyHVCF)
{
	register_exception_translator<HVCFException>(&translator);

	def("get_n_all_opened_objects", HVCF::get_n_all_opened_objects);

	class_<variant_info>("Variant", init<const char*, const char*, const char*, unsigned long long int>())
			.def_readonly("name", &variant_info::name)
			.def_readonly("ref", &variant_info::ref)
			.def_readonly("alt", &variant_info::alt)
			.def_readonly("position", &variant_info::position)
		;

	class_<variants_pair>("VariantsPair", init<const char*, unsigned long int, const char*, unsigned long int , double, double>())
			.def_readonly("name1", &VariantsPair::name1)
			.def_readonly("position1", &VariantsPair::position1)
			.def_readonly("name2", &VariantsPair::name2)
			.def_readonly("position2", &VariantsPair::position2)
			.def_readonly("r", &VariantsPair::r)
			.def_readonly("rsquare", &VariantsPair::rsquare)
		;

	class_<vector<string>>("NamesVector")
			.def(vector_indexing_suite<std::vector<std::string>>())
		;

	class_<vector<unsigned char>>("UCharVector")
			.def(vector_indexing_suite<std::vector<unsigned char>>())
		;

	class_<vector<variant_info>>("VariantsVector")
			.def(vector_indexing_suite<std::vector<variant_info>>())
		;

	class_<vector<variants_pair>>("PairsVector")
			.def(vector_indexing_suite<std::vector<variants_pair>>())
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
			.def("extract_variants", &HVCF::extract_variants)
			.def("extract_haplotypes", extract_haplotypes_for_variant)
			.def("extract_haplotypes", extract_haplotypes_for_sample)
			.def("get_n_opened_objects", &HVCF::get_n_opened_objects)
		;
}
