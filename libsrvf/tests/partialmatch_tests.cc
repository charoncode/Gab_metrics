#include <boost/test/unit_test.hpp>

#include <srvf/plf.h>
#include <srvf/srvf.h>
#include <srvf/qmap.h>
#include <srvf/numeric.h>
#include <srvf/partialmatch.h>

#define MY_CHECK_CLOSE(a,b) \
 BOOST_CHECK_EQUAL(srvf::numeric::almost_equal((a),(b)), true)

#define MY_REQUIRE_CLOSE(a,b) \
 BOOST_REQUIRE_EQUAL(srvf::numeric::almost_equal((a),(b)), true)


BOOST_AUTO_TEST_SUITE(partialmatch_tests)

BOOST_AUTO_TEST_CASE(edge_weight_test1)
{
  double samps1_data[] = { 1.0, -1.0, 1.0, -1.0 };
  double samps2_data[] = { -1.0, 1.0, -1.0 };

  double params1_data[] = { 0.0, .25, .5, .75, 1.0 };
  double params2_data[] = { 0.0, 1.0/3.0, 2.0/3.0, 1.0 };

  size_t ncp1 = sizeof(params1_data) / sizeof(double);
  size_t ncp2 = sizeof(params2_data) / sizeof(double);

  srvf::Pointset (1, ncp2, samps2_data);

  srvf::Srvf Q1( srvf::Pointset (1, ncp1-1, samps1_data), 
    std::vector<double>(&params1_data[0], &params1_data[ncp1]) );
  srvf::Srvf Q2( srvf::Pointset (1, ncp2-1, samps2_data), 
    std::vector<double>(&params2_data[0], &params2_data[ncp2]) );

  double tv1_data[] = { 0.0, 1.0/3.0, 2.0/3.0, 1.0 };
  double tv2_data[] = { 0.0, 0.5, 1.0 };

  std::vector<double> tv1(&tv1_data[0], &tv1_data[4]);
  std::vector<double> tv2(&tv2_data[0], &tv2_data[4]);
  
  srvf::pmatch::MatchingGraph G = 
    srvf::pmatch::calculate_edge_weights(Q1, Q2, tv1, tv2);

  srvf::pmatch::MatchingGraph Gexp(4,3);
  Gexp(0,0,1,1) = 1.513747151;
  Gexp(1,1,3,2) = 1.262891712;

  MY_CHECK_CLOSE(G(0,0,1,1), Gexp(0,0,1,1));
  MY_CHECK_CLOSE(G(1,1,3,2), Gexp(1,1,3,2));
}

BOOST_AUTO_TEST_SUITE_END()
