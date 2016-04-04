#include <boost/python.hpp>
#include <boost/python/def.hpp>
#include <boost/python/return_by_value.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>

#include "../src/include/HVCF.h"

using namespace sph_umich_edu;
using namespace boost::python;

BOOST_PYTHON_MODULE(PyHVCF)
{
	def("get_n_all_opened_objects", HVCF::get_n_all_opened_objects);

	class_<vector<string>>("vector")
			.def(vector_indexing_suite<std::vector<std::string>>())
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
			;
}
