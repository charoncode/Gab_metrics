#include <boost/test/unit_test.hpp>

#include <srvf/interp.h>
#include <srvf/util.h>

BOOST_AUTO_TEST_SUITE(interp_tests)

BOOST_AUTO_TEST_CASE(lookup_test1)
{
  double table_data[]={ 0.0, 0.00001, 0.5, 0.55, 1.0 };
  double tv[]={ -5.0, -0.0001, 0.0, 0.00001, 0.49999, 0.5, 0.9999, 1.0, 23.0 };
  size_t expv[]={ 0, 0, 0, 1, 1, 2, 3, 4, 4 };
  size_t ntable=sizeof(table_data)/sizeof(double);
  size_t ncases=sizeof(expv)/sizeof(size_t);
  std::vector<double> table(&table_data[0],&table_data[ntable]);
  for (size_t i=0; i<ncases; ++i)
  {
    size_t idx=srvf::interp::lookup(table,tv[i]);
    BOOST_CHECK_EQUAL(idx,expv[i]);
  }
}

BOOST_AUTO_TEST_CASE(interp_const_test1)
{
  double samps_data[] = {0.25, -0.5, 1.0, -1.5};
  double tv_data[] =
  {
    0.0, 0.1,   0.24999, 0.25, 0.3, 0.499, 
    0.5, 0.749, 0.75,    0.8,  0.99, 1.0
  };
  double exp_data[] =
  {
    0.25, 0.25, 0.25, -0.5, -0.5, -0.5, 
    1.0,  1.0,  -1.5, -1.5, -1.5, -1.5
  };
  size_t ntv=sizeof(tv_data)/sizeof(double);

  srvf::Pointset samps(1,4,samps_data);
  std::vector<double> params=srvf::util::linspace(0.0,1.0,5);
  std::vector<double> tv(&tv_data[0],&tv_data[ntv]);

  srvf::Pointset result = srvf::interp::interp_const(samps,params,tv);
  BOOST_REQUIRE_EQUAL(result.dim(), samps.dim());
  BOOST_REQUIRE_EQUAL(result.npts(), tv.size());
  for (size_t i=0; i<ntv; ++i)
  {
    BOOST_CHECK_CLOSE(result[i][0],exp_data[i],1e-9);
  }
}

BOOST_AUTO_TEST_SUITE_END()
