#include <boost/test/unit_test.hpp>

#include <srvf/srvf.h>
#include <srvf/plf.h>
#include <srvf/matrix.h>
#include <srvf/util.h>

BOOST_AUTO_TEST_SUITE(srvf_tests)

BOOST_AUTO_TEST_CASE(ctor_test1)
{
  srvf::Srvf z1(0.0, 1.0, std::vector<double>(2,0.25));
  BOOST_CHECK_EQUAL(z1.dim(), 2);
  BOOST_CHECK_EQUAL(z1.ncp(), 2);
  BOOST_CHECK_EQUAL(z1.domain_lb(), 0.0);
  BOOST_CHECK_EQUAL(z1.domain_ub(), 1.0);

  srvf::Point samps = z1.evaluate(0.5);
  BOOST_REQUIRE_EQUAL(samps.dim(), z1.dim());
  BOOST_CHECK_EQUAL(samps[0], 0.25);
  BOOST_CHECK_EQUAL(samps[1], 0.25);
}

BOOST_AUTO_TEST_CASE(ctor_test2)
{
  double samps_data[] = 
  {
    1.2, -0.5, 0.98, 100.01, -98.03, 
    0.0,  3.8, -0.5, -55.33, 100.00
  };
  double params_data[] = {0.0, 0.1, 0.12, 0.65, 0.8, 1.0};
  size_t ncp=sizeof(params_data)/sizeof(double);
  srvf::Pointset samps(2,ncp-1,samps_data,srvf::Pointset::POINT_PER_COLUMN);
  std::vector<double> params(&params_data[0],&params_data[ncp]);

  srvf::Srvf Q(samps,params);
  srvf::Srvf Qc1(Q);
  srvf::Srvf Qc2 = Q;

  BOOST_CHECK_EQUAL(Qc1.dim(), Q.dim());
  BOOST_CHECK_EQUAL(Qc1.ncp(), Q.ncp());
  BOOST_CHECK_EQUAL(Qc2.dim(), Q.dim());
  BOOST_CHECK_EQUAL(Qc2.ncp(), Q.ncp());
}

BOOST_AUTO_TEST_CASE(evaluate_test1)
{
  double samps_data[] = 
  {
    1.2, -0.5, 0.98, 100.01, -98.03, 
    0.0,  3.8, -0.5, -55.33, 100.00
  };
  double params_data[] = {0.0, 0.1, 0.12, 0.65, 0.8, 1.0};
  double tv_data[] = {
    0.0, 0.01, 0.0999, 0.1, 0.11, 0.11999, 
    0.12, 0.649, 0.65, 0.799, 0.8, 0.99, 1.0
  };
  double exp[] = 
  {
    1.2, 1.2, 1.2, -0.5, -0.5, -0.5, 0.98, 0.98, 
    100.01, 100.01, -98.03, -98.03, -98.03, 
    0.0, 0.0, 0.0, 3.8, 3.8, 3.8, -0.5, -0.5, 
    -55.33, -55.33, 100.00, 100.00, 100.00
  };
  size_t ncp=sizeof(params_data)/sizeof(double);
  size_t ntv=sizeof(tv_data)/sizeof(double);
  srvf::Pointset samps(2,ncp-1,samps_data,srvf::Pointset::POINT_PER_COLUMN);
  std::vector<double> params(&params_data[0],&params_data[ncp]);
  std::vector<double> tv(&tv_data[0],&tv_data[ntv]);

  srvf::Srvf Q(samps,params);
  srvf::Pointset result = Q.evaluate(tv);
  BOOST_REQUIRE_EQUAL(result.dim(), Q.dim());
  BOOST_REQUIRE_EQUAL(result.npts(), tv.size());
  for (size_t i=0; i<result.npts(); ++i)
  {
    for (size_t j=0; j<result.dim(); ++j)
    {
      BOOST_CHECK_CLOSE(result[i][j],exp[j*ntv+i],1e-9);
    }
  }
}

BOOST_AUTO_TEST_CASE(l2_norm_test1)
{
  double samps_data[] = {1.0, -1.0, 1.0, -1.0, -1.0};
  std::vector<double> params=srvf::util::linspace(0.0, 1.0, 6);
  srvf::Pointset samps(1, 5, samps_data);
  srvf::Srvf Q(samps, params);
  double exp=1.0;
  double nrm=srvf::l2_norm(Q);
  BOOST_CHECK_CLOSE(exp,nrm,1e-9);
}

BOOST_AUTO_TEST_CASE(l2_product_test1)
{
  double params1_data[] = {0.0, 0.25, 0.8, 1.0};
  double params2_data[] = {0.0, 0.5, 1.0};
  double samps1_data[] = {1.0, -1.0, 1.0};
  double samps2_data[] = {1.0, -1.0};
  double expected=0.1;

  size_t ncp1=sizeof(params1_data)/sizeof(double);
  size_t ncp2=sizeof(params2_data)/sizeof(double);

  std::vector<double> params1(&params1_data[0],&params1_data[ncp1]);
  std::vector<double> params2(&params2_data[0],&params2_data[ncp2]);
  srvf::Pointset samps1(1,ncp1-1,samps1_data);
  srvf::Pointset samps2(1,ncp2-1,samps2_data);
  srvf::Srvf Q1(samps1,params1);
  srvf::Srvf Q2(samps2,params2);

  double ip=srvf::l2_product(Q1,Q2);
  BOOST_CHECK_CLOSE(ip,expected,1e-9);
}

BOOST_AUTO_TEST_CASE(l2_product_test2)
{
  double params1_data[] = {0.0, 0.2500000, 0.8, 0.99995, 1.0};
  double params2_data[] = {0.0, 0.0000005, 0.25, 0.80001, 1.0};
  double samps1_data[] = {1.0, -1.0, 1.0, -1.0};
  double samps2_data[] = {1.0, -1.0, 1.0, -1.0};
  double expected=-0.999879;

  size_t ncp1=sizeof(params1_data)/sizeof(double);
  size_t ncp2=sizeof(params2_data)/sizeof(double);

  std::vector<double> params1(&params1_data[0],&params1_data[ncp1]);
  std::vector<double> params2(&params2_data[0],&params2_data[ncp2]);
  srvf::Pointset samps1(1,ncp1-1,samps1_data);
  srvf::Pointset samps2(1,ncp2-1,samps2_data);
  srvf::Srvf Q1(samps1,params1);
  srvf::Srvf Q2(samps2,params2);

  double ip=srvf::l2_product(Q1,Q2);
  BOOST_CHECK_CLOSE(ip,expected,1e-9);
}

BOOST_AUTO_TEST_CASE(l2_distance_test1)
{
  double params1_data[] = {0.0, 0.25, 0.8, 1.0};
  double params2_data[] = {0.0, 0.5, 1.0};
  double samps1_data[] = {1.0, -1.0, 1.0};
  double samps2_data[] = {1.0, -1.0};
  double expected=1.341640786;

  size_t ncp1=sizeof(params1_data)/sizeof(double);
  size_t ncp2=sizeof(params2_data)/sizeof(double);

  std::vector<double> params1(&params1_data[0],&params1_data[ncp1]);
  std::vector<double> params2(&params2_data[0],&params2_data[ncp2]);
  srvf::Pointset samps1(1,ncp1-1,samps1_data);
  srvf::Pointset samps2(1,ncp2-1,samps2_data);
  srvf::Srvf Q1(samps1,params1);
  srvf::Srvf Q2(samps2,params2);

  double d=srvf::l2_distance(Q1,Q2);
  BOOST_CHECK_CLOSE(d,expected,1e-4);
}

BOOST_AUTO_TEST_CASE(l2_distance_test2)
{
  double params1_data[] = {0.0, 0.2500000, 0.8, 0.99995, 1.0};
  double params2_data[] = {0.0, 0.0000005, 0.25, 0.80001, 1.0};
  double samps1_data[] = {1.0, -1.0, 1.0, -1.0};
  double samps2_data[] = {1.0, -1.0, 1.0, -1.0};
  double expected=1.999939499;

  size_t ncp1=sizeof(params1_data)/sizeof(double);
  size_t ncp2=sizeof(params2_data)/sizeof(double);

  std::vector<double> params1(&params1_data[0],&params1_data[ncp1]);
  std::vector<double> params2(&params2_data[0],&params2_data[ncp2]);
  srvf::Pointset samps1(1,ncp1-1,samps1_data);
  srvf::Pointset samps2(1,ncp2-1,samps2_data);
  srvf::Srvf Q1(samps1,params1);
  srvf::Srvf Q2(samps2,params2);

  double d=srvf::l2_distance(Q1,Q2);
  BOOST_CHECK_CLOSE(d,expected,1e-4);
}

BOOST_AUTO_TEST_CASE(rotate_test1)
{
  double samps_data[] = 
  {
    0.2, 0.15, -1.2, 0.99, 5.0, 
    -1.0, 0.0, 1.0, 1.0, 0.0
  };
  size_t ncp=6; // number of sample points +1
  size_t dim=2;
  double R_data[] = 
  {
    0.8660254038, -0.5,
    0.5, 0.8660254038
  };
  double exp_data[] = {
     0.673205, 0.129904, -1.539230, 0.357365, 4.330127,
    -0.766025, 0.075000,  0.266025, 1.361025, 2.500000
  };
  srvf::Pointset samps(dim,ncp-1,samps_data,srvf::Pointset::POINT_PER_COLUMN);
  srvf::Matrix R(dim,dim,R_data);

  srvf::Srvf Q(samps);
  Q.rotate(R);

  BOOST_REQUIRE_EQUAL(Q.dim(),dim);
  BOOST_REQUIRE_EQUAL(Q.ncp(),ncp);
  for (size_t i=0; i<Q.samps().npts(); ++i)
  {
    for (size_t j=0; j<Q.samps().dim(); ++j)
    {
      BOOST_CHECK_CLOSE(Q.samps()[i][j],exp_data[j*(ncp-1)+i],1e-3);
    }
  }
}

BOOST_AUTO_TEST_CASE(linear_combination_test1)
{
  double params1_data[] = {0.0, 0.25, 0.8, 1.0};
  double params2_data[] = {0.0, 0.5, 1.0};
  double samps1_data[] = {1.0, -1.0, 1.0};
  double samps2_data[] = {1.0, -1.0};
  double w1=0.25;
  double w2=0.5;
  double exp_params[] = {0.0, 0.25, 0.5, 0.8, 1.0};
  double exp_samps[] = {0.75, 0.25, -0.75, -0.25};

  size_t ncp1=sizeof(params1_data)/sizeof(double);
  size_t ncp2=sizeof(params2_data)/sizeof(double);
  size_t exp_ncp=sizeof(exp_params)/sizeof(double);

  std::vector<double> params1(&params1_data[0],&params1_data[ncp1]);
  std::vector<double> params2(&params2_data[0],&params2_data[ncp2]);
  srvf::Pointset samps1(1,ncp1-1,samps1_data);
  srvf::Pointset samps2(1,ncp2-1,samps2_data);
  srvf::Srvf Q1(samps1,params1);
  srvf::Srvf Q2(samps2,params2);
  srvf::Srvf Q=srvf::linear_combination(Q1,Q2,w1,w2);

  BOOST_REQUIRE_EQUAL(Q.dim(),Q1.dim());
  BOOST_REQUIRE_EQUAL(Q.ncp(),exp_ncp);
  for (size_t i=0; i<exp_ncp; ++i)
  {
    BOOST_CHECK_CLOSE(Q.params()[i],exp_params[i],1e-9);
  }
  for (size_t i=0; i<exp_ncp-1; ++i)
  {
    BOOST_CHECK_CLOSE(Q.samps()[i][0],exp_samps[i],1e-9);
  }
}

BOOST_AUTO_TEST_CASE(refinement_test1)
{
  double params_data[] = {0.0, 0.25, 0.8, 1.0};
  double tv_data[] = {0.0, 0.5, 0.99};
  double samps_data[] = {1.0, -1.0, 1.0};
  double exp_params[] = {0.0, 0.25, 0.5, 0.8, 0.99, 1.0};
  double exp_samps[] = {1.0, -1.0, -1.0, 1.0, 1.0};

  size_t ncp=sizeof(params_data)/sizeof(double);
  size_t ntv=sizeof(tv_data)/sizeof(double);
  size_t exp_ncp=sizeof(exp_params)/sizeof(double);

  std::vector<double> params(&params_data[0],&params_data[ncp]);
  srvf::Pointset samps(1,ncp-1,samps_data);
  std::vector<double> tv(&tv_data[0],&tv_data[ntv]);
  srvf::Srvf Q(samps,params);
  srvf::Srvf Qr=srvf::refinement(Q,tv);

  BOOST_REQUIRE_EQUAL(Qr.dim(),Q.dim());
  BOOST_REQUIRE_EQUAL(Qr.ncp(),exp_ncp);
  for (size_t i=0; i<exp_ncp; ++i)
  {
    BOOST_CHECK_CLOSE(Qr.params()[i],exp_params[i],1e-9);
  }
  for (size_t i=0; i<exp_ncp-1; ++i)
  {
    BOOST_CHECK_CLOSE(Qr.samps()[i][0],exp_samps[i],1e-9);
  }
}

BOOST_AUTO_TEST_CASE(gamma_action_test1)
{
  double Q_samps_data[] = {1.0};
  double Q_params_data[] = {0.0, 1.0};
  double gamma_samps_data[] = {0.0, 0.5, 1.0};
  double gamma_params_data[] = {0.0, 1.0/3.0, 1.0};
  double exp_params[] = {0.0, 1.0/3.0, 1.0};
  double exp_samps[] = {1.224744871, 0.8660254038};

  size_t Q_dim=1;
  size_t Q_ncp=sizeof(Q_params_data)/sizeof(double);
  size_t gamma_ncp=sizeof(gamma_samps_data)/sizeof(double);
  size_t exp_ncp=sizeof(exp_params)/sizeof(double);
  
  srvf::Pointset Q_samps(Q_dim,Q_ncp-1,Q_samps_data);
  std::vector<double> Q_params(&Q_params_data[0],&Q_params_data[Q_ncp]);
  srvf::Pointset gamma_samps(1,gamma_ncp,gamma_samps_data);
  std::vector<double> gamma_params(&gamma_params_data[0],
                                   &gamma_params_data[gamma_ncp]);

  srvf::Srvf Q(Q_samps,Q_params);
  srvf::Plf gamma(gamma_samps,gamma_params);
  srvf::Srvf Qgamma=srvf::gamma_action(Q,gamma);

  BOOST_REQUIRE_EQUAL(Qgamma.dim(),1);
  BOOST_REQUIRE_EQUAL(Qgamma.ncp(),exp_ncp);
  for (size_t i=0; i<exp_ncp; ++i)
  {
    BOOST_CHECK_CLOSE(Qgamma.params()[i],exp_params[i],1e-3);
  }
  for (size_t i=0; i<exp_ncp-1; ++i)
  {
    BOOST_CHECK_CLOSE(Qgamma.samps()[i][0],exp_samps[i],1e-3);
  }
}

BOOST_AUTO_TEST_CASE(constant_speed_test1)
{
  double samps_data[] = { M_SQRT2, -0.81649658 };
  double params_data[] = { 0.0, 0.25, 1.0 };
  double exp_samps[] = { 1.0, -1.0 };
  double exp_params[] = { 0.0, 0.5, 1.0 };

  srvf::Pointset samps(1, 2, samps_data);
  std::vector<double> params(&(params_data[0]), &(params_data[3]));
  srvf::Srvf Q(samps, params);
  srvf::Srvf Qcs = srvf::constant_speed_param(Q);

  BOOST_REQUIRE_EQUAL(Qcs.dim(), Q.dim());
  BOOST_REQUIRE_EQUAL(Qcs.ncp(), Q.ncp());
  
  for (size_t i=0; i<Qcs.params().size(); ++i)
  {
    BOOST_CHECK_SMALL(Qcs.params()[i] - exp_params[i], 1e-4);
  }
  for (size_t i=0; i<Qcs.samps().npts(); ++i)
  {
    BOOST_CHECK_SMALL(Qcs.samps()[i][0] - exp_samps[i], 1e-4);
  }
}

BOOST_AUTO_TEST_CASE(constant_speed_test2)
{
  double samps_data[] = { 1.0, 0.0, 1.0 };
  double exp_samps[] = { sqrt(2.0/3.0) };
  double exp_params[] = { 0.0, 1.0 };

  srvf::Pointset samps(1, 3, samps_data);
  std::vector<double> params = srvf::util::linspace(0.0, 1.0, 4);
  srvf::Srvf Q(samps, params);
  srvf::Srvf Qcs = srvf::constant_speed_param(Q);

  BOOST_REQUIRE_EQUAL(Qcs.dim(), Q.dim());
  BOOST_REQUIRE_EQUAL(Qcs.ncp(), 2);
  
  for (size_t i=0; i<Qcs.params().size(); ++i)
  {
    BOOST_CHECK_SMALL(Qcs.params()[i] - exp_params[i], 1e-4);
  }
  for (size_t i=0; i<Qcs.samps().npts(); ++i)
  {
    BOOST_CHECK_SMALL(Qcs.samps()[i][0] - exp_samps[i], 1e-4);
  }
}

BOOST_AUTO_TEST_CASE(constant_speed_test3)
{
  double samps_data[] = { 1.0, -1.0, 0.0, 1.0, 2.0 };
  double rsf = sqrt(7.0/5.0);
  double exp_samps[] = { rsf, -rsf, rsf };
  double exp_params[] = { 0.0, 1.0/7.0, 2.0/7.0, 1.0 };

  srvf::Pointset samps(1, 5, samps_data);
  std::vector<double> params = srvf::util::linspace(0.0, 1.0, 6);
  srvf::Srvf Q(samps, params);
  srvf::Srvf Qcs = srvf::constant_speed_param(Q);

  BOOST_REQUIRE_EQUAL(Qcs.dim(), Q.dim());
  BOOST_REQUIRE_EQUAL(Qcs.ncp(), 4);
  
  for (size_t i=0; i<Qcs.params().size(); ++i)
  {
    BOOST_CHECK_SMALL(Qcs.params()[i] - exp_params[i], 1e-4);
  }
  for (size_t i=0; i<Qcs.samps().npts(); ++i)
  {
    BOOST_CHECK_SMALL(Qcs.samps()[i][0] - exp_samps[i], 1e-4);
  }
}

BOOST_AUTO_TEST_SUITE_END()
