#include <boost/test/unit_test.hpp>

#ifdef USE_GSL

#include <srvf/srvf.h>
#include <srvf/rotate.h>
#include <srvf/matrix.h>


BOOST_AUTO_TEST_SUITE(rotate_tests)

BOOST_AUTO_TEST_CASE(optimal_rotation_test1)
{
  double samps1_data[] = 
  {
    0.0, -1.0,  0.0, 
    1.0,  0.0, -1.0
  };
  double samps2_data[] = 
  {
    1.0, 0.0, -1.0,
    0.0, 1.0,  0.0
  };
  double expected[2][2] = 
  {
    {0.0, -1.0},
    {1.0,  0.0}
  };
  srvf::Pointset samps1(2, 3, samps1_data, srvf::Pointset::POINT_PER_COLUMN);
  srvf::Pointset samps2(2, 3, samps2_data, srvf::Pointset::POINT_PER_COLUMN);
  srvf::Srvf Q1(samps1);
  srvf::Srvf Q2(samps2);
  srvf::Matrix R = srvf::optimal_rotation(Q1, Q2);
  
  BOOST_REQUIRE_EQUAL(R.rows(),2);
  BOOST_REQUIRE_EQUAL(R.cols(),2);
  for (size_t i=0; i<2; ++i)
  {
    for (size_t j=0; j<2; ++j)
    {
      BOOST_CHECK_CLOSE(R(i,j), expected[i][j], 1e-4);
    }
  }
}

BOOST_AUTO_TEST_SUITE_END()

#endif // USE_GSL
