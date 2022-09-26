#include <boost/test/unit_test.hpp>

#include <srvf/util.h>
#include <srvf/matrix.h>

BOOST_AUTO_TEST_SUITE(util_tests)

BOOST_AUTO_TEST_CASE(linspace_test1)
{
  double expdata[]={ 0.0, 0.25, 0.5, 0.75, 1.0 };
  std::vector<double> A=srvf::util::linspace(0.0,1.0,5);
  BOOST_CHECK_EQUAL(A.size(),5);
  for (size_t i=0; i<5; ++i)
  {
    BOOST_CHECK_EQUAL(A[i],expdata[i]);
  }
}

BOOST_AUTO_TEST_CASE(unique_test1)
{
  double v1_data[]={ 0.67, 0.23, 0.9, 0.42 };
  double v2_data[]={ 1.0, 0.0, 0.23, 0.89999, 0.4200001, 0.67  };
  double exp_data[]={ 0.0, 0.23, 0.42, 0.4200001, 0.67, 0.89999, 0.9, 1.0 };
  size_t nv1=sizeof(v1_data)/sizeof(double);
  size_t nv2=sizeof(v2_data)/sizeof(double);
  size_t nexp=sizeof(exp_data)/sizeof(double);
  std::vector<double> v1(&v1_data[0], &v1_data[nv1]);
  std::vector<double> v2(&v2_data[0], &v2_data[nv2]);
  std::vector<double> v=srvf::util::unique(v1, v2, 1e-9);
  BOOST_REQUIRE_EQUAL(v.size(),nexp);
  for (size_t i=0; i<v.size(); ++i)
  {
    BOOST_CHECK_CLOSE(v[i], exp_data[i],1e-9);
  }
}

BOOST_AUTO_TEST_CASE(diff_test1)
{
  srvf::Matrix Xdata(10,2,1.523);
  srvf::Pointset X(Xdata);
  srvf::Pointset dX=srvf::util::diff(X);
  BOOST_REQUIRE_EQUAL(dX.dim(),2);
  BOOST_REQUIRE_EQUAL(dX.npts(),9);
  for (size_t i=0; i<dX.npts(); ++i)
  {
    for (size_t j=0; j<dX.dim(); ++j)
    {
      BOOST_CHECK_CLOSE(dX[i][j],0.0,1e-9);
    }
  }
}

BOOST_AUTO_TEST_CASE(diff_test2)
{
  size_t dim=3, npts=11;
  srvf::Pointset X(dim,npts);
  std::vector<double> tv=srvf::util::linspace(0.0,1.0,npts);
  for (size_t i=0; i<npts; ++i)
  {
    for (size_t j=0; j<dim; ++j)
    {
      X[i][j]=(double)(i+j);
    }
  }
  srvf::Pointset dX=srvf::util::diff(X,tv);
  BOOST_REQUIRE_EQUAL(dX.dim(),3);
  BOOST_REQUIRE_EQUAL(dX.npts(),10);
  for (size_t i=0; i<dX.npts(); ++i)
  {
    for (size_t j=0; j<dX.dim(); ++j)
    {
      BOOST_CHECK_CLOSE(dX[i][j],10.0,1e-6);
    }
  }
}

BOOST_AUTO_TEST_SUITE_END()
