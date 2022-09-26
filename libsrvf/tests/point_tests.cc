#include <boost/test/unit_test.hpp>

#include <srvf/numeric.h>
#include <srvf/point.h>


#define MY_CHECK_CLOSE(a,b) \
 BOOST_CHECK_EQUAL(srvf::numeric::almost_equal((a),(b)), true)


BOOST_AUTO_TEST_SUITE(point_tests)

BOOST_AUTO_TEST_CASE(ctor_test1)
{
  srvf::Point P;
  BOOST_CHECK_EQUAL(P.dim(), 0);
}

BOOST_AUTO_TEST_CASE(ctor_test2)
{
  srvf::Point P(2);
  BOOST_CHECK_EQUAL(P.dim(), 2);
  BOOST_CHECK_EQUAL(P[0], 0.0);
  BOOST_CHECK_EQUAL(P[1], 0.0);
}

BOOST_AUTO_TEST_CASE(ctor_test3)
{
  double data[] = { 0.0, 1.0, 2.0, 3.0, 4.0 };
  size_t ndata = sizeof(data) / sizeof(double);

  srvf::Point P(&data[0], &data[ndata]);
  BOOST_REQUIRE_EQUAL(P.dim(), ndata);
  for (size_t i=0; i<P.dim(); ++i)
  {
    BOOST_CHECK_EQUAL(P[i], data[i]);
  }
}

BOOST_AUTO_TEST_CASE(ctor_test4)
{
  double data[] = { 0.0, 1.0, 2.0, 3.0, 4.0 };
  size_t ndata = sizeof(data) / sizeof(double);

  srvf::Point P(&data[0], &data[ndata]);
  BOOST_REQUIRE_EQUAL(P.dim(), ndata);
  for (size_t i=0; i<P.dim(); ++i)
  {
    BOOST_CHECK_EQUAL(P[i], data[i]);
  }
}

BOOST_AUTO_TEST_CASE(plus_minus_test1)
{
  double data1[] = { 0.0f, 1.0f, 2.0f, 3.0f, 4.0f };
  size_t ndata1 = sizeof(data1) / sizeof(double);
  double data2[] = { 0.0f, 1.0f, 2.0f, 3.0f, 4.0f };
  size_t ndata2 = sizeof(data2) / sizeof(double);
  double expdata[] = { 0.0f, 2.0f, 4.0f, 6.0f, 8.0f };
  size_t nexpdata = sizeof(expdata) / sizeof(double);

  srvf::Point P1(&data1[0], &data1[ndata1]);
  srvf::Point P2(&data2[0], &data2[ndata2]);
  srvf::Point P3 = P1 + P2;
  srvf::Point P4(P1);
  P4 += P2;

  BOOST_REQUIRE_EQUAL(P3.dim(), nexpdata);
  BOOST_REQUIRE_EQUAL(P4.dim(), nexpdata);

  for (size_t i=0; i<P3.dim(); ++i)
  {
    BOOST_CHECK_EQUAL(P3[i], expdata[i]);
    BOOST_CHECK_EQUAL(P4[i], expdata[i]);
  }

  P4 = P1 - P2;
  srvf::Point P5(P1);
  P5 -= P2;

  BOOST_REQUIRE_EQUAL(P3.dim(), nexpdata);
  BOOST_REQUIRE_EQUAL(P4.dim(), nexpdata);

  for (size_t i=0; i<P3.dim(); ++i)
  {
    BOOST_CHECK_EQUAL(P4[i], 0.0);
    BOOST_CHECK_EQUAL(P5[i], 0.0);
  }
}

BOOST_AUTO_TEST_CASE(times_test1)
{
  double data[] = { -0.1, 0.2, 10.5, -100.0, 0.0 };
  size_t ndata = sizeof(data) / sizeof(double);
  double exp_data[] = { 0.01, -0.02, -1.05, 10.0, -0.0 };
  double s = -0.1;

  srvf::Point P1(&data[0], &data[ndata]);
  srvf::Point P2 = P1 * s;
  P1 *= s;

  BOOST_REQUIRE_EQUAL(P1.dim(), ndata);
  BOOST_REQUIRE_EQUAL(P2.dim(), ndata);

  for (size_t i=0; i<P1.dim(); ++i)
  {
    MY_CHECK_CLOSE(P1[i], exp_data[i]);
    MY_CHECK_CLOSE(P2[i], exp_data[i]);
  }
}

BOOST_AUTO_TEST_CASE(norm_test1)
{
  double data[] = { 1.0, 1.0, 1.0, 1.0 };
  size_t ndata = sizeof(data) / sizeof(double);

  srvf::Point P(&data[0], &data[ndata]);
  double nrm = P.norm();
  MY_CHECK_CLOSE(nrm, 2.0);
}

BOOST_AUTO_TEST_CASE(dot_product_test1)
{
  double data1[] = { 1.0, 1.0, 1.0, 1.0 };
  size_t ndata1 = sizeof(data1) / sizeof(double);
  double data2[] = { -1.0, -0.5, -2.0, -0.25 };
  size_t ndata2 = sizeof(data2) / sizeof(double);

  srvf::Point P1(&data1[0], &data1[ndata1]);
  srvf::Point P2(&data2[0], &data2[ndata2]);

  double ip = P1.dot_product(P2);
  MY_CHECK_CLOSE(ip, -3.75);
}

BOOST_AUTO_TEST_CASE(distance_to_test1)
{
  double data1[] = { 1.0, 1.0, 1.0, 1.0 };
  size_t ndata1 = sizeof(data1) / sizeof(double);
  double data2[] = { -1.0, 0.0, -1.0, 1.0 };
  size_t ndata2 = sizeof(data2) / sizeof(double);

  srvf::Point P1(&data1[0], &data1[ndata1]);
  srvf::Point P2(&data2[0], &data2[ndata2]);

  double d = P1.distance_to(P2);
  MY_CHECK_CLOSE(d, 3.0);
}

BOOST_AUTO_TEST_CASE(equals_test1)
{
  double data1[] = { 0.0, 1.0, 1.0, 1000.0 };
  size_t ndata1 = sizeof(data1) / sizeof(double);
  double data2[] = { -0.0, 1.0000001, 1.0, 1000.5 };
  size_t ndata2 = sizeof(data2) / sizeof(double);

  srvf::Point P1(&data1[0], &data1[ndata1]);
  srvf::Point P2(&data2[0], &data2[ndata2]);

  BOOST_CHECK_EQUAL( (P1 == P2), true );
}

BOOST_AUTO_TEST_CASE(equals_test2)
{
  double data1[] = { 0.0, 1.0, 1.0, 1000.0 };
  size_t ndata1 = sizeof(data1) / sizeof(double);
  double data2[] = { 0.0, 1.0011, 1.0, 1000.0 };
  size_t ndata2 = sizeof(data2) / sizeof(double);

  srvf::Point P1(&data1[0], &data1[ndata1]);
  srvf::Point P2(&data2[0], &data2[ndata2]);
  srvf::Point P3(&data1[0], &data1[ndata1-1]);

  BOOST_CHECK_EQUAL( (P1 == P2), false );
  BOOST_CHECK_EQUAL( (P1 == P3), false );
}

BOOST_AUTO_TEST_CASE(less_test1)
{
  double data1[] = { 0.0, 1.0, 1.0, 1000.0, 10.0 };
  size_t ndata1 = sizeof(data1) / sizeof(double);
  double data2[] = { 0.0, 1.0, 1.0, 1000.1, 1.0 };
  size_t ndata2 = sizeof(data2) / sizeof(double);

  srvf::Point P1(&data1[0], &data1[ndata1]);
  srvf::Point P2(&data2[0], &data2[ndata2]);
  srvf::Point P3(&data2[0], &data2[3]);

  BOOST_CHECK_EQUAL( (P1 < P1), false );
  BOOST_CHECK_EQUAL( (P1 < P2), true );
  BOOST_CHECK_EQUAL( (P1 < P3), false );
  BOOST_CHECK_EQUAL( (P2 < P1), false );
  BOOST_CHECK_EQUAL( (P2 < P3), false );
  BOOST_CHECK_EQUAL( (P3 < P1), true );
  BOOST_CHECK_EQUAL( (P3 < P2), true );
}

BOOST_AUTO_TEST_CASE(is_multiple_test1)
{
  double data1[] = { 0.0, 1.0, -2.0, 3.0, 4.0 };
  size_t ndata1 = sizeof(data1) / sizeof(double);
  double data2[] = { 0.0, 0.5, -1.0, 1.5, 2.0 };
  size_t ndata2 = sizeof(data2) / sizeof(double);

  srvf::Point P1(&data1[0], &data1[ndata1]);
  srvf::Point P2(&data2[0], &data2[ndata2]);

  BOOST_CHECK_EQUAL(P1.on_same_ray(P2), true);
}

BOOST_AUTO_TEST_SUITE_END()
