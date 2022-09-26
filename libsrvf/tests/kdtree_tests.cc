#include <boost/test/unit_test.hpp>

#include <srvf/matrix.h>
#include <srvf/kdtree.h>
#include <srvf/numeric.h>

#include <cstdlib>
#include <vector>
#include <algorithm>

#define MY_CHECK_CLOSE(a,b) \
 BOOST_CHECK_EQUAL(srvf::numeric::almost_equal((a),(b)), true)

#define MY_REQUIRE_CLOSE(a,b) \
 BOOST_REQUIRE_EQUAL(srvf::numeric::almost_equal((a),(b)), true)


static srvf::Matrix randmat_(size_t r, size_t c, int modulus)
{
  srvf::Matrix result(r,c);
  for (size_t i=0; i<result.size(); ++i)
    result[i] = (double)(rand() % modulus);
  return result;
}


BOOST_AUTO_TEST_SUITE(kdtree_tests)

BOOST_AUTO_TEST_CASE(basic_test1)
{
  srvf::KdTree<srvf::Matrix> tree;
  std::vector<srvf::Matrix> mats;
  for (size_t i=0; i<5000; ++i)
  {
    mats.push_back(randmat_(1,8,2000));
    tree.insert(mats.back());
  }

  for (size_t i=0; i<100; ++i)
  {
    srvf::Matrix query = randmat_(1,8,2000);
    double expected_min_dist=1e9;
    for (size_t j=0; j<mats.size(); ++j)
    {
      double cur_dist = query.distance_to(mats[j]);
      expected_min_dist = std::min(expected_min_dist, cur_dist);
    }

    double actual_min_dist=1e9;
    srvf::Matrix actual_nbr=tree.get_nearest_neighbor(query, &actual_min_dist);
    MY_CHECK_CLOSE(actual_min_dist, expected_min_dist);
  }
}

BOOST_AUTO_TEST_CASE(basic_test2)
{
  srvf::KdTree<srvf::Matrix> tree;
  std::vector<srvf::Matrix> mats;
  std::vector<bool> intree;
  for (size_t i=0; i<100; ++i)
  {
    mats.push_back(randmat_(1,3,100));
    bool inserted = tree.insert_cond(mats.back(), 10.0);
    intree.push_back(inserted);

    double min_dist=1e9;
    for (size_t j=0; j+1<i; ++j)
    {
      if (intree[j])
        min_dist = std::min(min_dist, mats[i].distance_to(mats[j]));
    }

    BOOST_CHECK_EQUAL( (min_dist >= 10.0), inserted );
  }
}

BOOST_AUTO_TEST_CASE(timing_test1)
{
  srvf::KdTree<srvf::Matrix> tree;

  for (size_t i=0; i<10000; ++i)
  {
    tree.insert_cond(randmat_(1,3,10), 0.05);
  }
}

BOOST_AUTO_TEST_SUITE_END()
