#include <boost/test/unit_test.hpp>

#include <srvf/matrix.h>
#include <srvf/pointset.h>
#include <srvf/plf.h>
#include <srvf/srvf.h>
#include <srvf/qmap.h>
#include <srvf/util.h>

BOOST_AUTO_TEST_SUITE(qmap_tests)

BOOST_AUTO_TEST_CASE(plf_to_srvf_test1)
{
  double samps_data[]={0.0, 1.0/3.0, 0.0, 1.0/3.0};
  double exp_vals[]={1.0, -1.0, 1.0};
  std::vector<double> params=srvf::util::linspace(0.0,1.0,4);
  srvf::Pointset samps(1,4,samps_data);
  srvf::Plf F(samps,params);
  srvf::Srvf Q=srvf::plf_to_srvf(F);

  BOOST_REQUIRE_EQUAL(Q.dim(),1);
  BOOST_REQUIRE_EQUAL(Q.ncp(),4);

  // The Srvf should have same parameters as the Plf
  for (int i=0; i<4; ++i)
  {
    BOOST_CHECK_CLOSE(Q.params()[i],params[i],1e-9);
  }
  // Check samples
  for (int i=0; i<3; ++i)
  {
    BOOST_CHECK_CLOSE(Q.samps()[i][0],exp_vals[i],1e-9);
  }
}

BOOST_AUTO_TEST_CASE(plf_to_srvf_test2)
{
  double samps_data[]={
    0.0, 4.0, 4.0, 3.0, 2.0,
    0.0, 0.0, 2.0, 2.0, 1.0
  };
  double exp_vals[]={
    4.0, 0.0,      -2.0, -1.6817928,
    0.0, 2.828427,  0.0, -1.6817928
  };
  std::vector<double> params=srvf::util::linspace(0.0,1.0,5);
  srvf::Pointset samps(2,5,samps_data,srvf::Pointset::POINT_PER_COLUMN);
  srvf::Plf F(samps,params);
  srvf::Srvf Q=srvf::plf_to_srvf(F);

  BOOST_REQUIRE_EQUAL(Q.dim(),2);
  BOOST_REQUIRE_EQUAL(Q.ncp(),5);

  // The Srvf should have same parameters as the Plf
  for (size_t i=0; i<5; ++i)
  {
    BOOST_CHECK_CLOSE(Q.params()[i],params[i],1e-4);
  }
  // Check samples
  for (size_t i=0; i<4; ++i)
  {
    for (size_t j=0; j<2; ++j)
    {
      BOOST_CHECK_CLOSE(Q.samps()[i][j],exp_vals[j*4+i],1e-4);
    }
  }
}

BOOST_AUTO_TEST_CASE(srvf_to_plf_test1)
{
  double samps_data[] = {
    1.0, -1.0, 1.0, -1.0, 
    -1.0, 1.0, -1.0, 1.0
  };
  double exp_data[] = {
    0.0, M_SQRT2/4.0, 0.0, M_SQRT2/4.0, 0.0, 
    0.0, -M_SQRT2/4.0, 0.0, -M_SQRT2/4.0, 0.0
  };
  std::vector<double> params = srvf::util::linspace(0.0, 1.0, 5);
  srvf::Pointset samps(2, 4, samps_data, srvf::Pointset::POINT_PER_COLUMN);
  srvf::Srvf Q(samps, params);
  srvf::Plf F = srvf_to_plf(Q);
  for (size_t i=0; i<F.samps().dim(); ++i)
  {
    for (size_t j=0; j<F.samps().npts(); ++j)
    {
      BOOST_CHECK_SMALL(F.samps()[j][i] - exp_data[5*i+j], 1e-3);
    }
  }
}

BOOST_AUTO_TEST_SUITE_END()
