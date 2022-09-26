#include <boost/test/unit_test.hpp>

#include <srvf/numeric.h>


BOOST_AUTO_TEST_SUITE(numeric_tests)


BOOST_AUTO_TEST_CASE(almost_equal_test1)
{
  BOOST_CHECK_EQUAL(srvf::numeric::almost_equal(0.0, 0.0000005), true);
  BOOST_CHECK_EQUAL(srvf::numeric::almost_equal(0.0, 0.001), false);
  BOOST_CHECK_EQUAL(srvf::numeric::almost_equal(4.9999995, 5.0), true);
  BOOST_CHECK_EQUAL(srvf::numeric::almost_equal(5.0251, 5.0), false);
}


BOOST_AUTO_TEST_SUITE_END()
