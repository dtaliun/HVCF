#include <boost/python.hpp>
#include <boost/python/def.hpp>

#include "../src/include/HVCF.h"

using namespace sph_umich_edu;
using namespace boost::python;

BOOST_PYTHON_MODULE(PyHVCF)
{
	class_<HVCF, boost::noncopyable>("HVCF")
			.def("open", &HVCF::open)
			.def("close", &HVCF::close)
			.def("get_n_samples", &HVCF::get_n_samples)
			.def("get_n_variants", &HVCF::get_n_variants)
			.def("get_n_variants_in_chromosome", &HVCF::get_n_variants_in_chromosome)
			;
}
