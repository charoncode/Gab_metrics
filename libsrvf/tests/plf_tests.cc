#include <boost/test/unit_test.hpp>

#include <srvf/plf.h>
#include <srvf/pointset.h>
#include <srvf/matrix.h>
#include <srvf/util.h>

BOOST_AUTO_TEST_SUITE(plf_tests)

BOOST_AUTO_TEST_CASE(evaluate_test1)
{
  std::vector<double> samps_data=srvf::util::linspace(0.0,1.0,2);
  srvf::Pointset samps(1,2,samps_data);
  srvf::Plf F(samps);
  std::vector<double> tv=srvf::util::linspace(0.0,1.0,5);
  double expdata[]={0.0, 0.25, 0.5, 0.75, 1.0};
  srvf::Pointset Ftv1 = F.evaluate(tv);
  BOOST_REQUIRE_EQUAL(Ftv1.dim(), F.dim());
  BOOST_REQUIRE_EQUAL(Ftv1.npts(), tv.size());
  for (size_t i=0; i<5; ++i)
  {
    srvf::Point Ftv2 = F.evaluate(tv[i]);
    BOOST_REQUIRE_EQUAL(Ftv2.dim(), F.dim());
    BOOST_CHECK_EQUAL(Ftv1[i][0],expdata[i]);
    BOOST_CHECK_EQUAL(Ftv2[0],expdata[i]);
  }
}

BOOST_AUTO_TEST_CASE(evaluate_test2)
{
  std::vector<double> uv=srvf::util::linspace(0.0,1.0,500);
  std::vector<double> uvl=srvf::util::linspace(0.0,1.0,999);
  srvf::Pointset X(3,500);
  for (size_t i=0; i<3; ++i)
  {
    for (size_t j=0; j<500; ++j)
    {
      X[j][i]=(double)j;
    }
  }
  srvf::Plf F(X,uv);
  srvf::Pointset res = F.evaluate(uvl);
  BOOST_REQUIRE_EQUAL(res.dim(), F.dim());
  BOOST_REQUIRE_EQUAL(res.npts(), uvl.size());
  for (size_t i=0; i<3; ++i)
  {
    for (size_t j=0; j<999; ++j)
    {
      double ev=(double)(j) * 499.0 / 998.0;
      BOOST_CHECK_CLOSE(ev,res[j][i],1e-5);
    }
  }
}

BOOST_AUTO_TEST_CASE(preimages_test1)
{
  std::vector<double> samps_data=srvf::util::linspace(0.0,1.0,2);
  std::vector<double> uv=srvf::util::linspace(0.0,1.0,13);
  srvf::Pointset samps(1,2,samps_data);
  srvf::Plf F(samps);
  std::vector<double> Fiuv = F.preimages(uv);
  for (size_t i=0; i<uv.size(); ++i)
  {
    BOOST_CHECK_CLOSE(uv[i],Fiuv[i],1e-5);
  }
}

BOOST_AUTO_TEST_CASE(preimages_test2)
{
  double samps_data[]={0.0, 0.5, 0.5, 1.0};
  srvf::Pointset samps(1,4,samps_data);
  std::vector<double> tv=srvf::util::linspace(0.0,1.0,5);
  double exp_data[]={0.0, 1.0/6.0, 2.0/3.0, 5.0/6.0, 1.0};
  srvf::Plf F(samps);
  std::vector<double> Fitv = F.preimages(tv);
  for (size_t i=0; i<tv.size(); ++i)
  {
    BOOST_CHECK_CLOSE(exp_data[i],Fitv[i],1e-5);
  }
}

BOOST_AUTO_TEST_CASE(arc_length_test1)
{
  double samps_data[]=
  {
    0.0, 1.0, 1.0, 0.0,
    0.0, 0.0, 1.0, 1.0
  };
  srvf::Plf F(srvf::Pointset(2,4,samps_data,srvf::Pointset::POINT_PER_COLUMN));
  double ev=3.0;
  double av=F.arc_length();
  BOOST_CHECK_CLOSE(ev,av,1e-5);
}

BOOST_AUTO_TEST_CASE(translate_test1)
{
  srvf::Pointset samps(3,20,0.0);
  srvf::Plf F(samps);
  double v_data[]={1.0,2.0,3.0};
  srvf::Point v(&v_data[0],&v_data[3]);
  F.translate(v);
  for (size_t i=0; i<F.samps().npts(); ++i)
  {
    for (size_t j=0; j<F.samps().dim(); ++j)
    {
      BOOST_CHECK_EQUAL(F.samps()[i][j],v[j]);
    }
  }
}

BOOST_AUTO_TEST_CASE(rotate_test1)
{
  // Rotation by 90 degrees counter-clockwise
  double R_data[]={
    0.0, -1.0,
    1.0, 0.0
  };
  srvf::Matrix R(2,2,R_data);
  srvf::Pointset samps(2,3,1.0);
  srvf::Plf F(samps);
  F.rotate(R);
  for (size_t i=0; i<F.samps().npts(); ++i)
  {
    BOOST_CHECK_EQUAL(F.samps()[i][0],-1.0);
    BOOST_CHECK_EQUAL(F.samps()[i][1],1.0);
  }
}

BOOST_AUTO_TEST_CASE(scale_test1)
{
  srvf::Matrix samps_data(3,57,1.0);
  srvf::Pointset samps(samps_data);
  srvf::Plf F(samps);
  F.scale(2.0);
  for (size_t i=0; i<F.samps().npts(); ++i)
  {
    for (size_t j=0; j<F.samps().dim(); ++j)
    {
      BOOST_CHECK_EQUAL(F.samps()[i][j],2.0);
    }
  }
}

BOOST_AUTO_TEST_CASE(linear_combination_test1)
{
  double params1_data[]= { 0.0, 0.1, 0.2, 0.9, 1.0 };
  double params2_data[]= { 0.0, 0.3, 0.4, 0.6, 1.0 };
  size_t ncp1=sizeof(params1_data)/sizeof(double);
  size_t ncp2=sizeof(params2_data)/sizeof(double);
  srvf::Pointset samps1(1,5,1.0);
  srvf::Pointset samps2(1,5,-1.0);
  std::vector<double> params1(&params1_data[0],&params1_data[ncp1]);
  std::vector<double> params2(&params2_data[0],&params2_data[ncp2]);
  double exp_params_data[]={ 0.0, 0.1, 0.2, 0.3, 0.4, 0.6, 0.9, 1.0 };
  srvf::Plf F1(samps1,params1);
  srvf::Plf F2(samps2,params2);
  srvf::Plf Fr=srvf::linear_combination(F1,F2,0.75,0.25);

  BOOST_REQUIRE_EQUAL(Fr.samps().npts(),8);
  BOOST_REQUIRE_EQUAL(Fr.ncp(),8);
  for (size_t i=0; i<8; ++i)
  {
    BOOST_CHECK_EQUAL(Fr.samps()[i][0],0.5);
    BOOST_CHECK_EQUAL(Fr.params()[i],exp_params_data[i]);
  }
}

BOOST_AUTO_TEST_CASE(composition_test1)
{
  double params1_data[]={0.0, 0.25, 0.5, 1.0};
  double params2_data[]={0.0, 1.0/3.0, 2.0/3.0, 1.0};
  double exp_params[]={0.0, 1.0/6.0, 1.0/3.0, 2.0/3.0, 1.0};
  double samps1_data[]={0.0, 0.5, 0.0, 0.5};
  double samps2_data[]={0.0, 0.5, 0.5, 1.0};
  double exp_samps[]={0.0, 0.5, 0.0, 0.0, 0.5};
  size_t n1=sizeof(params1_data)/sizeof(double);
  size_t n2=sizeof(params2_data)/sizeof(double);
  size_t n3=sizeof(exp_params)/sizeof(double);
  std::vector<double> params1(&params1_data[0],&params1_data[n1]);
  std::vector<double> params2(&params2_data[0],&params2_data[n2]);
  srvf::Pointset samps1(1,n1,samps1_data);
  srvf::Pointset samps2(1,n2,samps2_data);
  srvf::Plf F1(samps1,params1);
  srvf::Plf F2(samps2,params2);
  srvf::Plf F12=srvf::composition(F1,F2);
  
  BOOST_REQUIRE_EQUAL(F12.ncp(),n3);
  for (size_t i=0; i<n3; ++i)
  {
    BOOST_CHECK_EQUAL(F12.params()[i],exp_params[i]);
    BOOST_CHECK_EQUAL(F12.samps()[i][0],exp_samps[i]);
  }
}

BOOST_AUTO_TEST_CASE(inverse_test1)
{
  double params_data[]={0.0, 0.25, 0.5, 0.75, 1.0};
  double samps_data[]={0.0, 1.0/3.0, 1.0/3.0, 2.0/3.0, 1.0};
  size_t n1=sizeof(params_data)/sizeof(double);
  size_t nexp=n1;
  std::vector<double> params(&params_data[0],&params_data[n1]);
  srvf::Pointset samps(1,n1,samps_data);
  srvf::Plf F(samps,params);
  srvf::Plf Fi=srvf::inverse(F);
  
  BOOST_REQUIRE_EQUAL(Fi.ncp(),nexp);
  for (size_t i=0; i<nexp; ++i)
  {
    BOOST_CHECK_EQUAL(Fi.params()[i],samps_data[i]);
    BOOST_CHECK_EQUAL(Fi.samps()[i][0],params_data[i]);
  }
}

BOOST_AUTO_TEST_SUITE_END()
