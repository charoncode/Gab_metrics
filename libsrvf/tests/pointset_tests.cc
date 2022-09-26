#include <boost/test/unit_test.hpp>

#include <srvf/point.h>
#include <srvf/pointset.h>

#include <algorithm>
#include <iterator>


#define MY_CHECK_CLOSE(a,b) \
 BOOST_CHECK_EQUAL(srvf::numeric::almost_equal((a),(b)), true)

#define MY_REQUIRE_CLOSE(a,b) \
 BOOST_REQUIRE_EQUAL(srvf::numeric::almost_equal((a),(b)), true)


BOOST_AUTO_TEST_SUITE(pointset_tests)

BOOST_AUTO_TEST_CASE(ctor_test1)
{
  srvf::Pointset S1(2,5);
  BOOST_REQUIRE_EQUAL(S1.dim(), 2);
  BOOST_REQUIRE_EQUAL(S1.npts(), 5);

  for (size_t i=0; i<5; ++i)
    for (size_t j=0; j<2; ++j)
      MY_CHECK_CLOSE(S1[i][j], 0.0);
}

BOOST_AUTO_TEST_CASE(ctor_test2)
{
  srvf::Pointset S1(6,500,-42.0);
  BOOST_REQUIRE_EQUAL(S1.dim(), 6);
  BOOST_REQUIRE_EQUAL(S1.npts(), 500);

  for (size_t i=0; i<500; ++i)
    for (size_t j=0; j<6; ++j)
      MY_CHECK_CLOSE(S1[i][j], -42.0);
}

BOOST_AUTO_TEST_CASE(ctor_test3)
{
  double points_data[] = {
    0.0, 1.0, 2.0, 
    3.0, 4.0, 5.0
  };
  srvf::Matrix points(2,3,points_data);
  srvf::Pointset X1(points,srvf::Pointset::POINT_PER_ROW);
  srvf::Pointset X2(points,srvf::Pointset::POINT_PER_COLUMN);

  BOOST_REQUIRE_EQUAL(X1.dim(),3);
  BOOST_REQUIRE_EQUAL(X1.npts(),2);
  BOOST_REQUIRE_EQUAL(X2.dim(),2);
  BOOST_REQUIRE_EQUAL(X2.npts(),3);

  for (int i=0; i<2; ++i)
  {
    for (int j=0; j<3; ++j)
    {
      MY_CHECK_CLOSE(X1[i][j],points_data[3*i+j]);
      MY_CHECK_CLOSE(X2[j][i],points_data[3*i+j]);
    }
  }
}

BOOST_AUTO_TEST_CASE(ctor_test4)
{
  size_t dim=3;
  size_t npts=2;
  double samps_data[] = { 
    0.0, 1.0, 2.0, 
    3.0, 4.0, 5.0 
  };
  srvf::Pointset S(dim, npts, 
    std::vector<double>(&samps_data[0], &samps_data[dim*npts]));

  BOOST_REQUIRE_EQUAL(S.dim(), dim);
  BOOST_REQUIRE_EQUAL(S.npts(), npts);

  for (size_t i=0; i<npts; ++i)
    for (size_t j=0; j<dim; ++j)
      MY_CHECK_CLOSE(S[i][j], samps_data[3*i+j]);
}

BOOST_AUTO_TEST_CASE(push_test1)
{
  srvf::Point P(3,0.2);
  srvf::Pointset S;

  for (int i=1; i<=5; ++i)
  {
    S.push_back(P);
    BOOST_CHECK_EQUAL(S.dim(),3);
    BOOST_CHECK_EQUAL(S.npts(),i);
    for (size_t j=0; j<P.dim(); ++j)
    {
      MY_CHECK_CLOSE(S[S.npts()-1][j], P[j]);
    }
  }
}

BOOST_AUTO_TEST_CASE(distance_and_dot_product_test1)
{
  double points_data[] = 
  {
    0.0, 1.0, 2.0, 
    3.0, 4.0, 5.0
  };
  double dists1[2][2] = 
  {
    {0.0, 5.196152423},
    {5.196152423, 0.0}
  };
  double dists2[3][3] = 
  {
    {0.0, 1.414213562, 2.828427125},
    {1.414213562, 0.0, 1.414213562},
    {2.828427125, 1.414213562, 0.0}
  };
  double dots1[2][2] = 
  {
    {5.0, 14.0},
    {14.0, 50.0}
  };
  double dots2[3][3] = 
  {
    {9.0, 12.0, 15.0},
    {12.0, 17.0, 22.0},
    {15.0, 22.0, 29.0}
  };
  srvf::Matrix data(2,3,points_data);
  srvf::Pointset X1(data,srvf::Pointset::POINT_PER_ROW);
  srvf::Pointset X2(data,srvf::Pointset::POINT_PER_COLUMN);

  for (int i=0; i<2; ++i)
  {
    for (int j=0; j<2; ++j)
    {
      MY_CHECK_CLOSE(X1.distance(i,j),dists1[i][j]);
      MY_CHECK_CLOSE(X1.dot_product(i,j),dots1[i][j]);
    }
  }

  for (int i=0; i<3; ++i)
  {
    for (int j=0; j<3; ++j)
    {
      MY_CHECK_CLOSE(X2.distance(i,j),dists2[i][j]);
      MY_CHECK_CLOSE(X2.dot_product(i,j),dots2[i][j]);
    }
  }
}

BOOST_AUTO_TEST_SUITE_END()
