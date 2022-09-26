/*
 * LibSRVF - a shape analysis library using the square root velocity framework.
 *
 * Copyright (C) 2012   FSU Statistical Shape Analysis and Modeling Group
 * 
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>
 */
#include <srvf/functions.h>

#include <srvf/srvf.h>
#include <srvf/plf.h>

#include <stack>
#include <algorithm>
#include <iterator>
#include <iostream>
#include <stdexcept>


namespace srvf
{

namespace functions
{


static inline bool
is_peak_(const Srvf &Q, size_t idx)
{
  return ( (idx>0 && Q.samps()[idx-1][0] > 0.0) ||
           (idx+1<Q.ncp() && Q.samps()[idx][0] < 0.0) );
}


static inline double
time_up_(const Srvf &Q, size_t i1, size_t i2)
{
  double v=0.0;
  for (size_t i=i1; i<i2; i+=2)
  {
    v += Q.params()[i+1] - Q.params()[i];
  }
  return v;
}


static inline double
edge_variation_(const Srvf &Q, size_t i1, size_t i2)
{
  double v=0.0;
  for (size_t i=i1; i<i2; i+=2)
  {
    double Qi = Q.samps()[i][0];
    v += Qi * fabs(Qi) * (Q.params()[i+1]-Q.params()[i]);
  }
  return v;
}


static inline double 
edge_score_(const Srvf &Q1, const Srvf &Q2, 
  size_t sc, size_t sr, size_t tc, size_t tr)
{
  //double v1 = edge_variation_(Q1,sc,tc);
  //double v2 = edge_variation_(Q2,sr,tr);
  double v1 = time_up_(Q1,sc,tc);
  double v2 = time_up_(Q2,sr,tr);
  return sqrt(v1*v2);
}


static void calculate_scores_(const Srvf &Q1, const Srvf &Q2, 
  std::map<match_vertex_t,double> &score, 
  std::map<match_vertex_t,match_vertex_t> &pred)
{
  size_t n1 = Q1.ncp();
  size_t n2 = Q2.ncp();

  // Dynamic programming step: flood-fill the score matrix.
  for (size_t tr=1; tr<n2; ++tr)
  {
    size_t tc=1;
    if (is_peak_(Q1,tc) != is_peak_(Q2,tr)) ++tc;

    for (/* NOOP */; tc<n1; tc+=2)
    {
      match_vertex_t tv(tc,tr);
      double best_score = -1.0;
      match_vertex_t best_pred;
      
      for (size_t sr=1-(tr%2); sr<tr; sr+=2)
      {
        for (size_t sc=1-(tc%2); sc<tc; sc+=2)
        {
          match_vertex_t sv(sc,sr);
          double w = edge_score_(Q1,Q2,sc,sr,tc,tr);
          double s = score[sv] + w;
          if (s > best_score)
          {
            best_score = s;
            best_pred = sv;
          }
        }
      }

      if (best_score > -1.0)
      {
        score[tv] = best_score;
        pred[tv] = best_pred;
      }
    }
  }

  //for (size_t r=n2-1; r<n2; --r)
  //{
  //  for (size_t c=0; c<n1; ++c)
  //  {
  //    if (score.find(match_vertex_t(c,r)) != score.end())
  //    {
  //      std::cout << "[" << score[match_vertex_t(c,r)] << " : " << 
  //        "(" << pred[match_vertex_t(c,r)].first << "," <<
  //        pred[match_vertex_t(c,r)].second << ")-->" << 
  //        "(" << c << "," << r << ")] ";
  //    }
  //    else
  //    {
  //      std::cout << "[(" << c << "," << r << ") xxxxxxxxxxxxxxxxxx] ";
  //    }
  //  }
  //  std::cout << std::endl;
  //}
  //std::cout << std::endl;
}


// Finds the optimal path through the matching graph for Q1 and Q2.
static std::deque<match_vertex_t>
build_matching_path_ (const Srvf &Q1, const Srvf &Q2)
{
  std::map<match_vertex_t,double> score;
  std::map<match_vertex_t,match_vertex_t> pred;

  size_t n1 = Q1.ncp();
  size_t n2 = Q2.ncp();

  // All paths start at start_vertex.  If both functions start with a peak, 
  // or both functions start with a valley, then the start vertex is (0,0).  
  // Otherwise, we make a dummy vertex, and add zero-score edges from it 
  // to (1,0) and (0,1).
  match_vertex_t start_vertex;
  if (is_peak_(Q1,0) == is_peak_(Q2,0))
  {
    start_vertex = match_vertex_t(0,0);
    score[start_vertex] = 0.0;
  }
  else
  {
    start_vertex = match_vertex_t(n1,n2);  // dummy vertex: invalid indices
    score[start_vertex] = 0.0;
    score[match_vertex_t(0,1)] = 0.0;
    pred[match_vertex_t(0,1)] = start_vertex;
    score[match_vertex_t(1,0)] = 0.0;
    pred[match_vertex_t(1,0)] = start_vertex;
  }

  // Call the main DP routine
  calculate_scores_(Q1, Q2, score, pred);

  // Select the end vertex.  If both functions end with a peak or both 
  // functions end with a valley, the end vertex is (n1-1,n2-1).  Otherwise, 
  // the end vertex is (n1-2,n2-1) or (n1-1,n2-2), depending on which one 
  // has the higher score.
  match_vertex_t end_vertex(n1-1,n2-1);
  if (is_peak_(Q1,n1-1) != is_peak_(Q2,n2-1))
  {
    match_vertex_t end_cand1(n1-2,n2-1);
    match_vertex_t end_cand2(n1-1,n2-2);

    if (score[end_cand1] > score[end_cand2])
      end_vertex = end_cand1;
    else
      end_vertex = end_cand2;
  }

  // Reconstruct the optimal path from start_vertex to end_vertex.
  std::deque<match_vertex_t> path;
  path.push_front(end_vertex);
  while (path.front() != start_vertex) path.push_front(pred[path.front()]);

  return path;
}


// Builds the pieces of gamma1 and gamma2 corresponding to 
// the edge (sc,sr)-->(tc,tr) in the matching graph.  The 
// output is appended to G1samps and G2samps.
static inline void 
build_gamma_segment_(const Srvf &Q1, const Srvf &Q2, 
  size_t sc, size_t sr, size_t tc, size_t tr, 
  std::vector<double> &G1samps, std::vector<double> &G2samps)
{
  double v1 = time_up_(Q1, sc, tc);
  double v2 = time_up_(Q2, sr, tr);
  double R = v2 / v1;
  double u1 = 0.0;
  double u2 = 0.0;
  double t1 = Q1.params()[sc];
  double t2 = Q2.params()[sr];
  size_t c = sc;
  size_t r = sr;

  while (c<tc && r<tr)
  {
    if (c+1 < tc || r+1 < tr)
    {
      double p1 = (u1 + (Q1.params()[c+1] - t1)) / v1;
      double p2 = (u2 + (Q2.params()[r+1] - t2)) / v2;

      if (fabs(p1-p2) > 1e-4)
      {
        if (p1 < p2)
        {
          double lambda = t2 + R*(Q1.params()[c+1]-t1);
          G1samps.push_back(Q1.params()[c+1]);
          G1samps.push_back(Q1.params()[c+2]);
          G2samps.push_back(lambda);
          G2samps.push_back(lambda);
          u1 += (Q1.params()[c+1]-t1);
          u2 += (lambda-t2);
          t1 = Q1.params()[c+2];
          t2 = lambda;
          c += 2;
        }
        else
        {
          double lambda = t1 + (Q2.params()[r+1]-t2) / R;
          G1samps.push_back(lambda);
          G1samps.push_back(lambda);
          G2samps.push_back(Q2.params()[r+1]);
          G2samps.push_back(Q2.params()[r+2]);
          u1 += (lambda-t1);
          u2 += (Q2.params()[r+1]-t2);
          t1 = lambda;
          t2 = Q2.params()[r+2];
          r += 2;
        }
      }
      else
      {
        G1samps.push_back(Q1.params()[c+1]);
        G1samps.push_back(Q1.params()[c+2]);
        G2samps.push_back(Q2.params()[r+1]);
        G2samps.push_back(Q2.params()[r+2]);
        u1 += (Q1.params()[c+1]-t1);
        u2 += (Q2.params()[r+1]-t2);
        t1 = Q1.params()[c+2];
        t2 = Q2.params()[r+2];
        c += 2;
        r += 2;
      }
    }
    else
    {
      G1samps.push_back(Q1.params()[tc]);
      G2samps.push_back(Q2.params()[tr]);
      break;
    } 
  }
}


/**
 * Given a number t between \c Q1.params()[cs] and 
 * \c Q1.params()[ct], returns the smallest number s between 
 * \c Q2.params()[rs] and \c Q2.params()[rt] such that the percent of the 
 * total rise (or fall) up to s (relative to the total rise of F2 for the 
 * whole matching segment) is equal to the percent of the total rise (or fall) 
 * up to t (relative to the total rise of F1 for the whole matching segment).
 */
static double segment_translate_(const Srvf &Q1, const Srvf &Q2, 
  size_t cs, size_t rs, size_t ct, size_t rt, double t)
{
  // Check for special cases, where the source vertex is the dummy vertex
  if (cs >= Q1.ncp())
  {
    // target vertex is either (0,1) or (1,0)
    if (ct == 0)
      return (t < 1e-6 ? 0.0 : 1e9);
    else
      return (t < Q1.params()[1]+1e-6 ? 0.0 : 1e9);
  }

  // Check for easy cases, where t is at an end of the segment
  if (t < Q1.params()[cs] + 1e-6) 
    return (Q2.params()[rs]);
  if (t > Q1.params()[ct] - 1e-6) 
    return (Q2.params()[rt]);

  // Get rise of Q1 up to t, and for the whole segment
  double total_rise1 = 0.0;
  double rise_to_t = 0.0;
  for (size_t i=cs; i<ct; i+=2)
  {
    total_rise1 += (Q1.params()[i+1]-Q1.params()[i]);
    if (t >= Q1.params()[i+1])
      rise_to_t += (Q1.params()[i+1]-Q1.params()[i]);
    else if (t >= Q1.params()[i])
      rise_to_t += (t-Q1.params()[i]);
  }
  double pct_rise_to_t = rise_to_t / total_rise1;

  // Get total rise of Q2 for the whole segment
  double total_rise2 = 0.0;
  for (size_t i=rs; i<rt; i+=2)
  {
    total_rise2 += (Q2.params()[i+1]-Q2.params()[i]);
  }

  // Now search for s
  double needed_rise_to_s = total_rise2 * pct_rise_to_t;
  for (size_t i=rs; i<rt; i+=2)
  {
    // Notice the i+=2 in the for loop above -- we're only looping through 
    // the rise intervals
    double interval_rise = (Q2.params()[i+1] - Q2.params()[i]);
    if (interval_rise <= needed_rise_to_s)
    {
      needed_rise_to_s -= interval_rise;
    } else {
      return (Q2.params()[i] + needed_rise_to_s);
    }
  }

  return Q2.params()[rt];
}



// IMPORTANT:  Mu and all elements of Qs must already be normalized, so 
// that they only take on values in {1,-1}.
std::vector<Plf> optimal_reparam(const Srvf &Q1, const Srvf &Q2)
{
  std::deque<match_vertex_t> path = build_matching_path_(Q1, Q2);
  return build_gammas(Q1, Q2, path);
}


std::vector<Plf> build_gammas(const Srvf&Q1, const Srvf &Q2, 
  const std::deque<match_vertex_t> &path)
{
  size_t n1 = Q1.ncp();
  size_t n2 = Q2.ncp();
  size_t path_idx=0;

  // Translate the path into a pair of piecewise-linear reparametrization
  // functions gamma1 and gamma2.
  std::vector<double> G1samps(1,0.0);
  std::vector<double> G2samps(1,0.0);

  if ( path[path_idx].first == n1 )
  {
    // Start vertex is a dummy vertex, so we consume another vertex from 
    // our matching path.
    ++path_idx;
    G1samps.push_back(Q1.params()[path[path_idx].first]);
    G2samps.push_back(Q2.params()[path[path_idx].second]);
  }

  size_t sc = path[path_idx].first; 
  size_t sr = path[path_idx].second;
  ++path_idx;

  while (path_idx < path.size())
  {
    size_t tc = path[path_idx].first;
    size_t tr = path[path_idx].second;

    build_gamma_segment_(Q1,Q2,sc,sr,tc,tr,G1samps,G2samps);

    sr = tr;
    sc = tc;
    ++path_idx;
  }

  // Handle case where one function ends with a peak and the other ends 
  // with a valley.
  if (sc < n1-1 || sr < n2-1)
  {
    G1samps.push_back(Q1.domain_ub());
    G2samps.push_back(Q2.domain_ub());
  }

  // The parameters for gamma1 and gamma2 are arbitrary, as long as we use 
  // the same parameters for both functions.  We use uniformly-spaced 
  // parameters in [0,1] for both gamma1 and gamma2.
  std::vector<double> G1params = 
    util::linspace(0.0, 1.0, G1samps.size());
  std::vector<double> G2params = 
    util::linspace(0.0, 1.0, G2samps.size());

  std::vector<Plf> res;
  res.push_back(Plf(Pointset(1,G1samps.size(),G1samps), G1params));
  res.push_back(Plf(Pointset(1,G2samps.size(),G2samps), G2params));
  return res;
}


std::vector<Plf> 
groupwise_optimal_reparam(const Srvf &Mu, const std::vector<Srvf> &Qs)
{
  std::vector<std::deque<match_vertex_t> > paths(Qs.size());
  for (size_t i=0; i<Qs.size(); ++i)
  {
    paths[i] = build_matching_path_(Mu, Qs[i]);
  }
  return groupwise_build_gammas(Mu, Qs, paths);
}


std::vector<Plf> 
groupwise_build_gammas(const Srvf &Mu, const std::vector<Srvf> &Qs,
  const std::vector<std::deque<match_vertex_t> > &paths)
{
  size_t nfuncs = Qs.size();
  std::vector<size_t> alpha(nfuncs,0);  // index of last changepoint output
  size_t alpha_mu = 0;
  std::vector<size_t> beta(nfuncs,0);   // matching path index

  std::vector<Pointset> Gsamps(nfuncs,Pointset(1,1,0.0));
  Pointset Esamps(1,1,0.0);

  while (alpha_mu+1 < Mu.ncp())
  {
    // Find the parameter of Mu corresponding to the next peak or valley
    std::vector<double> cands(nfuncs);
    double next_param = Mu.params()[alpha_mu+1];
    for (size_t i=0; i<nfuncs; ++i)
    {
      if (beta[i]+1 < paths[i].size()){
        cands[i] = segment_translate_(Qs[i], Mu, 
          paths[i][beta[i]].second, paths[i][beta[i]].first, 
          paths[i][beta[i]+1].second, paths[i][beta[i]+1].first, 
          Qs[i].params()[alpha[i]+1]);
      }
      else
        cands[i] = 1.0;

      if (cands[i] < next_param)
      {
        next_param = cands[i];
      }
    }

    if (srvf::numeric::almost_equal(next_param, Mu.params()[alpha_mu+1]))
    {
      double tau_mu = next_param;
      Esamps.push_back(Point(1,next_param));
      ++alpha_mu;

      for (size_t i=0; i<nfuncs; ++i)
      {
        double tau_i;
        if (beta[i]+1 >= paths[i].size())
        {
          tau_i = 1.0;
        }
        else if (alpha_mu == paths[i][beta[i]+1].first)
        {
          tau_i = Qs[i].params()[paths[i][beta[i]+1].second];
          alpha[i] = paths[i][beta[i]+1].second;
          ++beta[i];
        }
        else
        {
          tau_i = segment_translate_(Mu, Qs[i], 
            paths[i][beta[i]].first, paths[i][beta[i]].second, 
            paths[i][beta[i]+1].first, paths[i][beta[i]+1].second, 
            tau_mu );
        }

        Gsamps[i].push_back(Point(1,tau_i));
      }
    }
    else
    {
      double tau_mu = next_param;
      Esamps.push_back(Point(1,tau_mu));
      for (size_t i=0; i<nfuncs; ++i)
      {
        double tau_i;
        if (beta[i] == paths[i].size())
        {
          tau_i = Qs[i].params()[alpha[i]];
        }
        else if (srvf::numeric::almost_equal(cands[i],next_param))
        {
          tau_i = Qs[i].params()[alpha[i]+1];
          ++alpha[i];
          if (alpha[i] == paths[i][beta[i]+1].second)
          {
            // This only happens during the first iteration, when the first 
            // segment of Qs[i] gets crushed to a point on Mu
            ++beta[i];
          }
        }
        else
        {
          tau_i = segment_translate_(Mu, Qs[i], 
            paths[i][beta[i]].first, paths[i][beta[i]].second, 
            paths[i][beta[i]+1].first, paths[i][beta[i]+1].second, 
            tau_mu);
        }

        Gsamps[i].push_back(Point(1,tau_i));
      }
    }
  }

  // Check for trailing half-zero columns
  bool found_trailing = false;
  for (size_t i=0; i<nfuncs; ++i)
  {
    if (alpha[i]+1 < Qs[i].ncp())
    {
      found_trailing = true;
      break;
    }
  }

  if (found_trailing)
  {
    Esamps.push_back(Point(1,1.0));
    for (size_t i=0; i<nfuncs; ++i)
    {
      Gsamps[i].push_back(Point(1,1.0));
    }
  }

  // Create the PLFs
  std::vector<Plf> result;
  for (size_t i=0; i<nfuncs; ++i)
  {
    result.push_back( 
      Plf(Gsamps[i], srvf::util::linspace(0.0, 1.0, Gsamps[i].npts())) );
  }
  result.push_back( 
    Plf(Esamps, srvf::util::linspace(0.0, 1.0, Esamps.npts())) );

  return result;
}


Srvf karcher_mean(const std::vector<Srvf> &Qs, double tol, size_t max_iters)
{
  // Make sure Qs all have the same domain and the same L^2 norm
  double Qnorm = l2_norm(Qs[0]);
  double domain_lb = Qs[0].domain_lb();
  double domain_ub = Qs[0].domain_ub();
  for (size_t i=1; i<Qs.size(); ++i)
  {
    if ( !srvf::numeric::almost_equal(Qs[i].domain_lb(), domain_lb) )
      throw std::invalid_argument("SRVFs in Qs have different domains");
    if ( !srvf::numeric::almost_equal(Qs[i].domain_ub(), domain_ub) )
      throw std::invalid_argument("SRVFs in Qs have different domains");
    if ( !srvf::numeric::almost_equal(l2_norm(Qs[i]), Qnorm) )
      throw std::invalid_argument("points in Qs lie on different spheres");
  }

  // Initialize the mean estimate.  This can be done in several 
  // ways; we simply make it a copy of Qs[0].
  Srvf Mu = Qs[0];

  for (size_t iter=0; iter < max_iters; ++iter)
  {
    // Replace mean with its orbit's unit-speed representative
    Mu = constant_speed_param(Mu);

    // Compute unit-speed orbit representatives for each Qs[i].
    std::vector<Srvf> Qsl;
    for (size_t i=0; i<Qs.size(); ++i)
    {
      Qsl.push_back(constant_speed_param(Qs[i]));
    }


    // Reparametrize to optimally register Mu to each of the Qsl[i]
    std::vector<Plf> GE = groupwise_optimal_reparam(Mu, Qsl);
    for (size_t i=0; i<Qsl.size(); ++i)
    {
      Qsl[i] = gamma_action(Qsl[i], GE[i]);
    }
    Mu = gamma_action(Mu, GE[Qsl.size()]);

    // Compute inner products and great-circle distances between 
    // reparametrized Qsl[i] and Mu
    std::vector<double> ips(Qsl.size());
    std::vector<double> dists(Qsl.size());
    for (size_t i=0; i<Qsl.size(); ++i)
    {
      ips[i] = l2_product(Mu, Qsl[i]);
      dists[i] = sphere_distance(Mu, Qsl[i]);
    }

    // Compute the negative of the gradient of the energy function.  
    // This is the average of the shooting vectors.
    // NGrad is initialized to the zero function.
    Srvf NGrad(domain_lb, domain_ub, std::vector<double>(1,0.0));
    for (size_t i=0; i<Qsl.size(); ++i)
    {
      Srvf W = linear_combination(Qsl[i], Mu, 1.0, -ips[i] / (Qnorm * Qnorm) );
      double Wnorm = l2_norm(W);
      if (Wnorm > 1e-4)
        NGrad = linear_combination(NGrad, W, 1.0, dists[i] / l2_norm(W) );
    }
    NGrad.scale( 1.0 / Qsl.size() );

    // Check for termination condition: is gradient norm below 
    // tol % of the sphere radius?
    double NGrad_norm = l2_norm(NGrad);
    std::cout << "Karcher mean:  gradient norm = " << NGrad_norm << std::endl;
    if ( (NGrad_norm / Qnorm) < tol )
    {
      break;
    }

    // Update
    Mu = linear_combination(Mu, NGrad, 1.0, 0.25);
    Mu.scale( Qnorm / l2_norm(Mu) );
  }

  return Mu;
}


////////////////////////////////////////////////////////////////////////
// Helper functions for unit testing
////////////////////////////////////////////////////////////////////////


double TestAccess::edge_variation(const Srvf &Q, size_t i1, size_t i2)
{
  return edge_variation_(Q, i1, i2);
}

double TestAccess::edge_score(const Srvf &Q1, const Srvf &Q2, 
  size_t sc, size_t sr, size_t tc, size_t tr)
{
  return edge_score_(Q1, Q2, sc, sr, tc, tr);
}

void TestAccess::calculate_scores(const Srvf &Q1, const Srvf &Q2, 
  std::map<match_vertex_t,double> &score, 
  std::map<match_vertex_t,match_vertex_t> &pred)
{
  calculate_scores_(Q1, Q2, score, pred);
}

void TestAccess::build_gamma_segment(const Srvf &Q1, const Srvf &Q2, 
  size_t sc, size_t sr, size_t tc, size_t tr, 
  std::vector<double> &G1samps, std::vector<double> &G2samps)
{
  build_gamma_segment_(Q1, Q2, sc, sr, tc, tr, G1samps, G2samps);
}

double TestAccess::segment_translate(const Srvf &Q1, const Srvf &Q2, 
  size_t cs, size_t rs, size_t ct, size_t rt, double t)
{
  return segment_translate_(Q1, Q2, cs, rs, ct, rt, t);
}

} // namespace srvf::functions

} // namespace srvf
