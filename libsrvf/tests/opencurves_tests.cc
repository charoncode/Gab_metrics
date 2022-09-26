#include <fstream>
#include <cmath>
#include <boost/test/unit_test.hpp>

#include <srvf/matrix.h>
#include <srvf/pointset.h>
#include <srvf/srvf.h>
#include <srvf/qmap.h>
#include <srvf/opencurves.h>


BOOST_AUTO_TEST_SUITE(opencurves_tests)

BOOST_AUTO_TEST_CASE(shooting_vector_test1)
{
  double samps1_data[] = {1.0};
  double samps2_data[] = {1.0, -1.0};
  double exp_data[] = {M_PI_2, -M_PI_2};
  double exp_params_data[] = {0.0, 0.5, 1.0};
  size_t exp_ncp = sizeof(exp_params_data) / sizeof(double);

  srvf::Pointset samps1(1, 1, samps1_data);
  srvf::Pointset samps2(1, 2, samps2_data);
  srvf::Srvf Q1(samps1);
  srvf::Srvf Q2(samps2);

  srvf::Srvf Sv = srvf::opencurves::shooting_vector(Q1, Q2);
  BOOST_CHECK_EQUAL(Sv.dim(), 1);
  BOOST_CHECK_EQUAL(Sv.ncp(), exp_ncp);
  for (size_t i=0; i<exp_ncp; ++i)
  {
    BOOST_CHECK_EQUAL(Sv.params()[i], exp_params_data[i]);
  }
  for (size_t i=0; i+1<exp_ncp; ++i)
  {
    BOOST_CHECK_EQUAL(Sv.samps()[i][0], exp_data[i]);
  }
}


BOOST_AUTO_TEST_SUITE_END()
