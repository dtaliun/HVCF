#include <boost/python.hpp>
#include <boost/python/def.hpp>
#include <boost/python/return_by_value.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>

#include "../src/include/HVCF.h"

using namespace sph_umich_edu;
using namespace boost::python;

void (HVCF::*compute_region_ld)(const string& chromosome, unsigned long long int start_position, unsigned long long end_position, vector<variants_pair>& result) = &HVCF::compute_ld;
void (HVCF::*compute_lead_ld)(const string& chromosome, const string& lead_variant_name, unsigned long long int start_position, unsigned long long end_position, vector<variants_pair>& result) = &HVCF::compute_ld;

BOOST_PYTHON_MODULE(PyHVCF)
{
	def("get_n_all_opened_objects", HVCF::get_n_all_opened_objects);

	class_<variants_pair>("VariantsPair", init<const char*, unsigned long int, const char*, unsigned long int , double, double>())
			.def_readonly("name1", &VariantsPair::name1)
			.def_readonly("position1", &VariantsPair::position1)
			.def_readonly("name2", &VariantsPair::name2)
			.def_readonly("position2", &VariantsPair::position2)
			.def_readonly("r", &VariantsPair::r)
			.def_readonly("rsquare", &VariantsPair::rsquare)
		;

	class_<vector<string>>("Vector")
			.def(vector_indexing_suite<std::vector<std::string>>())
		;

	class_<vector<variants_pair>>("PairsVector")
			.def(vector_indexing_suite<std::vector<variants_pair>>())
		;

	class_<HVCF, boost::noncopyable>("HVCF")
			.def("open", &HVCF::open)
			.def("close", &HVCF::close)
			.def("get_n_samples", &HVCF::get_n_samples)
			.def("get_n_variants", &HVCF::get_n_variants)
			.def("get_n_variants_in_chromosome", &HVCF::get_n_variants_in_chromosome)
			.def("get_samples", &HVCF::get_samples, return_value_policy<return_by_value>())
			.def("get_population", &HVCF::get_population, return_value_policy<return_by_value>())
			.def("get_n_opened_objects", &HVCF::get_n_opened_objects)
			.def("compute_ld", compute_region_ld)
			.def("compute_ld", compute_lead_ld)
		;
}
