#include <cmath>
#include <boost/test/unit_test.hpp>

#include <srvf/matrix.h>
#include <srvf/pointset.h>
#include <srvf/plf.h>
#include <srvf/srvf.h>
#include <srvf/qmap.h>
#include <srvf/functions.h>
using srvf::functions::match_vertex_t;
#include <srvf/util.h>

#include <vector>
#include <deque>
#include <iostream>
#include <iterator>
#include <time.h>
#include <cstdlib>


#define MY_CHECK_CLOSE(a,b) \
 BOOST_CHECK_EQUAL(srvf::numeric::almost_equal((a),(b)), true)

#define MY_REQUIRE_CLOSE(a,b) \
 BOOST_REQUIRE_EQUAL(srvf::numeric::almost_equal((a),(b)), true)


BOOST_AUTO_TEST_SUITE(functions_tests)

BOOST_AUTO_TEST_CASE(edge_variation_test1)
{
  double samps_data[] = { 1.0, -1.0, 1.0, -1.0 };
  std::vector<double> params = srvf::util::linspace(0.0, 1.0, 5);
  size_t i1 = 1;  // params[1] is a peak
  size_t i2 = 4;  // params[4] is a valley
  double expected = -0.5;

  srvf::Pointset samps(1, 4, samps_data);
  srvf::Srvf Q(samps, params);

  double actual = srvf::functions::TestAccess::edge_variation(Q, i1, i2);
  BOOST_CHECK_SMALL(actual-expected,1e-6);
}

BOOST_AUTO_TEST_CASE(edge_score_test1)
{
  double samps1_data[] = { 1.0, -1.0, 1.0, -1.0 };
  double samps2_data[] = { -1.0, 1.0 };
  std::vector<double> params1 = srvf::util::linspace(0.0, 1.0, 5);
  std::vector<double> params2 = srvf::util::linspace(0.0, 1.0, 3);
  size_t sc=1, sr=0;
  size_t tc=4, tr=1;
  double expected = 0.5;

  srvf::Pointset samps1(1, 4, samps1_data);
  srvf::Pointset samps2(1, 2, samps2_data);
  srvf::Srvf Q1(samps1, params1);
  srvf::Srvf Q2(samps2, params2);

  double actual = srvf::functions::TestAccess::edge_score(
    Q1, Q2, sc, sr, tc, tr);
  BOOST_CHECK_SMALL(actual-expected,1e-6);
}

BOOST_AUTO_TEST_CASE(build_gamma_segment_test1)
{
  double samps1_data[] = { 1.0, -1.0, 1.0, -1.0, 1.0, -1.0, 1.0, -1.0 };
  double samps2_data[] = { 1.0, -1.0, 1.0 };
  double G1exp[] = { 0.0, 0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875 };
  double G2exp[] = { 0.0, 0.166667, 0.166667, 0.333333, 0.666667, 
                     0.833333, 0.833333, 1.0 };
  size_t Glen = sizeof(G1exp) / sizeof(double);
  size_t ncp1 = sizeof(samps1_data) / sizeof(double) + 1;
  size_t ncp2 = sizeof(samps2_data) / sizeof(double) + 1;
  std::vector<double> params1 = srvf::util::linspace(0.0, 1.0, ncp1);
  std::vector<double> params2 = srvf::util::linspace(0.0, 1.0, ncp2);
  size_t sc=0, sr=0;
  size_t tc=7, tr=3;

  srvf::Pointset samps1(1, ncp1-1, samps1_data);
  srvf::Pointset samps2(1, ncp2-1, samps2_data);
  srvf::Srvf Q1(samps1, params1);
  srvf::Srvf Q2(samps2, params2);
  std::vector<double> G1samps(1,0.0);
  std::vector<double> G2samps(1,0.0);

  srvf::functions::TestAccess::build_gamma_segment(
    Q1, Q2, sc, sr, tc, tr, G1samps, G2samps);

  BOOST_REQUIRE_EQUAL(G1samps.size(), Glen);
  BOOST_REQUIRE_EQUAL(G2samps.size(), Glen);

  for (size_t i=0; i<Glen; ++i)
  {
    BOOST_CHECK_SMALL(G1exp[i] - G1samps[i], 1e-3);
    BOOST_CHECK_SMALL(G2exp[i] - G2samps[i], 1e-3);
  }
}


BOOST_AUTO_TEST_CASE(calculate_scores_test1)
{
  double samps1_data[] = { 0.0, 0.2, 0.0, 0.2, 0.0, 0.2 };
  double samps2_data[] = { 1.0/3.0, 0.0, 1.0/3.0, 0.0 };
  size_t nsamps1 = sizeof(samps1_data) / sizeof(double);
  size_t nsamps2 = sizeof(samps2_data) / sizeof(double);
  srvf::Plf F1(srvf::Pointset(1, nsamps1, samps1_data), 
               srvf::util::linspace(0.0, 1.0, nsamps1) );
  srvf::Plf F2(srvf::Pointset(1, nsamps2, samps2_data), 
               srvf::util::linspace(0.0, 1.0, nsamps2) );
  srvf::Srvf Q1 = srvf::plf_to_srvf(F1);
  srvf::Srvf Q2 = srvf::plf_to_srvf(F2);
  std::map<match_vertex_t,double> score;
  std::map<match_vertex_t,match_vertex_t> preds;
  match_vertex_t start_vertex(nsamps1,nsamps2);
  score[start_vertex] = 0.0;
  score[match_vertex_t(0,1)] = 0.0;
  score[match_vertex_t(1,0)] = 0.0;
  preds[match_vertex_t(0,1)] = start_vertex;
  preds[match_vertex_t(1,0)] = start_vertex;
  srvf::functions::TestAccess::calculate_scores(Q1,Q2,score,preds);

}


BOOST_AUTO_TEST_CASE(segment_translate_test1)
{
  double samps1_data[]={ 1.0, -1.0, 1.0, -1.0 };
  double samps2_data[]={ 1.0, -1.0, 1.0, -1.0, 1.0 };
  size_t nsamps1 = sizeof(samps1_data) / sizeof(double);
  size_t nsamps2 = sizeof(samps2_data) / sizeof(double);
  srvf::Srvf Q1(srvf::Pointset(1,nsamps1,samps1_data));
  srvf::Srvf Q2(srvf::Pointset(1,nsamps2,samps2_data));

  std::deque<match_vertex_t> path;
  path.push_back(match_vertex_t(0,0));
  path.push_back(match_vertex_t(3,1));
  path.push_back(match_vertex_t(4,4));

  // Cases for translation from Q1 to Q2
  size_t segs12[] = { 0, 0, 0, 0, 0, 0,
                      1, 1, 1, 1 };
  double ts12[]= { 0.0, 0.125, 0.25, 0.5, 0.625, 0.75, 
                   0.75, 10.0/12.0, 11.0/12.0, 1.0 };
  double exps12[] = { 0.0, 0.05, 0.1, 0.1, 0.15, 0.2, 
                      0.2, 1.0/3.0, 2.0/3.0, 0.8 };
  size_t ncases12 = sizeof(segs12) / sizeof(size_t);

  for (size_t i=0; i<ncases12; ++i)
  {
    double s = srvf::functions::TestAccess::segment_translate(
      Q1, Q2, path[segs12[i]].first, path[segs12[i]].second, 
      path[segs12[i]+1].first, path[segs12[i]+1].second, ts12[i]);
    MY_CHECK_CLOSE(s, exps12[i]);
  }
}


BOOST_AUTO_TEST_SUITE_END()
