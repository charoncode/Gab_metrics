#include <boost/test/unit_test.hpp>

#include <srvf/minmaxheap.h>

#include <cstdlib>
#include <vector>


template <class T>
static void knuth_shuffle_(std::vector<T> v)
{
  for (size_t i=0; i+1<v.size(); ++i)
  {
    size_t j = i + ((size_t)rand()) % (v.size()-i);
    if (j != i) { T tmp = v[i]; v[i] = v[j]; v[j] = tmp; }
  }
}


template <class T>
static std::vector<T> linspace_(T a, T b)
{
  std::vector<T> result;
  for (T t=a; t<=b; result.push_back(t++));
  return result;
}


BOOST_AUTO_TEST_SUITE(minmaxheap_tests)

BOOST_AUTO_TEST_CASE(basic_test1)
{
  std::vector<int> data = linspace_(1,100);
  knuth_shuffle_(data);
  srvf::MinMaxHeap<int> H1(data);
  srvf::MinMaxHeap<int> H2;
  for (size_t i=0; i<data.size(); H2.insert(data[i++]));

  for (size_t i=1; i<=100; ++i)
  {
    BOOST_CHECK_EQUAL(H1.min(), i);
    BOOST_CHECK_EQUAL(H2.min(), i);

    H1.remove_min();
    H2.remove_min();
  }
}

BOOST_AUTO_TEST_CASE(basic_test2)
{
  std::vector<int> data = linspace_(1,100);
  knuth_shuffle_(data);
  srvf::MinMaxHeap<int> H1(data);
  srvf::MinMaxHeap<int> H2;
  for (size_t i=0; i<data.size(); H2.insert(data[i++]));

  for (size_t i=100; i>=1; --i)
  {
    BOOST_CHECK_EQUAL(H1.max(), i);
    BOOST_CHECK_EQUAL(H2.max(), i);

    H1.remove_max();
    H2.remove_max();
  }
}

BOOST_AUTO_TEST_SUITE_END()
