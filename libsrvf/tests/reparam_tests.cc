#include <boost/test/unit_test.hpp>

#include <srvf/reparam.h>
#include <srvf/srvf.h>
#include <srvf/plf.h>
#include <srvf/pointset.h>
#include <srvf/util.h>

BOOST_AUTO_TEST_SUITE(reparam_tests)

BOOST_AUTO_TEST_CASE(match_cost_test1)
{
  double samps1_data[]={1.0, -1.0};
  double samps2_data[]={1.0,  1.0};
  srvf::Pointset samps1(1, 2, samps1_data);
  srvf::Pointset samps2(1, 2, samps2_data);
  std::vector<double> tv1 = srvf::util::linspace(0.0, 1.0, 3);
  std::vector<double> tv2 = srvf::util::linspace(0.0, 1.0, 3);
  srvf::Srvf Q1(samps1, tv1);
  srvf::Srvf Q2(samps2, tv1);
  size_t Q1_start_idx[] = { 0, 0, 0, 0, 1 };
  size_t Q2_start_idx[] = { 0, 1, 0, 0, 0 };
  double exp_costs[]=
  {0.08578643763, 1.5, 2.0, 1.5, 2.914213562};
  size_t idx[5][4]={
    {0, 1, 0, 2},
    {0, 2, 1, 2},
    {0, 2, 0, 2},
    {0, 2, 0, 1},
    {1, 2, 0, 2}
  };

  for (size_t i=0; i<5; ++i)
  {
    double c = srvf::opencurves::edge_weight(
        Q1, Q2, tv1, tv2, 
        idx[i][0], idx[i][1], idx[i][2], idx[i][3],
        Q1_start_idx[i], Q2_start_idx[i]);
    BOOST_CHECK_CLOSE(c, exp_costs[i], 1e-4);
  }
}

BOOST_AUTO_TEST_CASE(match_cost_test2)
{
  double samps1_data[] =
  {
    1.0/3.0, -2.0, 1.0,
    2.0,     -1.0, -1.0
  };
  double samps2_data[] =
  {
    0.5, -0.5, 2.0,  0.9,
    0.2,  0.1, 0.3, -1.1
  };
  srvf::Pointset samps1(2, 3, samps1_data, srvf::Pointset::POINT_PER_COLUMN);
  srvf::Pointset samps2(2, 4, samps2_data, srvf::Pointset::POINT_PER_COLUMN);
  std::vector<double> tv1 = srvf::util::linspace(0.0, 1.0, 4);
  std::vector<double> tv2 = srvf::util::linspace(0.0, 1.0, 5);
  srvf::Srvf Q1(samps1, tv1);
  srvf::Srvf Q2(samps2, tv2);
  
  double c = srvf::opencurves::edge_weight(
    Q1, Q2, tv1, tv2, 1, 3, 1, 4, 1, 1);
  BOOST_CHECK_CLOSE(c, 3.1716, 1e-3);
}

BOOST_AUTO_TEST_SUITE_END()
