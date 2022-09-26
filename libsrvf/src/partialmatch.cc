/*
 * LibSRVF - a shape analysis library using the square root velocity framework.
 *
 * Copyright (C) 2012  FSU Statistical Shape Analysis and Modeling Group
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

#include <srvf/partialmatch.h>
#include <srvf/paretoset.h>
#include <srvf/rotselect.h>

#include <srvf/srvf.h>
#include <srvf/reparam.h>
#include <srvf/dpnbhd.h>
#include <srvf/util.h>

#include <cstdlib>
#include <cmath>
#include <vector>
#include <map>

namespace srvf
{

namespace pmatch
{


/**
 * Calculates the weights of all edges in the matching graph.
 */
MatchingGraph calculate_edge_weights (
  const srvf::Srvf &Q1, const srvf::Srvf &Q2, 
  const std::vector<double> &tv1, const std::vector<double> &tv2 )
{
  std::map<size_t,size_t> tv1_idx_to_Q1_idx = 
    srvf::util::build_lookup_map(tv1,Q1.params());
  std::map<size_t,size_t> tv2_idx_to_Q2_idx = 
    srvf::util::build_lookup_map(tv2,Q2.params());

  MatchingGraph result(tv1.size(), tv2.size());

  for (size_t ct=1; ct<tv1.size(); ++ct)
  {
    for (size_t rt=1; rt<tv2.size(); ++rt)
    {
      for (size_t i=0; i<DP_NBHD_SIZE; ++i)
      {
        size_t cs = ct - srvf::dp_nbhd[i][0];
        size_t rs = rt - srvf::dp_nbhd[i][1];
        if (cs < tv1.size() && rs < tv2.size()){
          double w = srvf::opencurves::edge_weight (
            Q1, Q2, tv1, tv2, 
            cs, ct, rs, rt, 
            tv1_idx_to_Q1_idx[cs], 
            tv2_idx_to_Q2_idx[rs] );

          if (w < 0.0)
            w = 0.0;

          result(cs, rs, ct, rt) = w;
        }
      }
    }
  }

  return result;
}


/**
 * Use the Floyd-Warshall algorithm to solve the all-pairs-shortest 
 * path problem in \a G.
 */
void calculate_match_scores (MatchingGraph &G)
{
  for (size_t cm=1; cm+1<G.width(); ++cm)
  {
  for (size_t rm=1; rm+1<G.height(); ++rm)
  {

    for (size_t cs=0; cs<cm; ++cs)
    {
    for (size_t rs=0; rs<rm; ++rs)
    {

      for (size_t i=0; i<DP_NBHD_SIZE; ++i)
      {
        size_t ct = cm + srvf::dp_nbhd[i][0];
        size_t rt = rm + srvf::dp_nbhd[i][1];

        if (ct >= G.width() || rt >= G.height()) continue;

        double cand_w = G(cs, rs, cm, rm) + G(cm, rm, ct, rt);
        if (cand_w < G(cs, rs, ct, rt))
          G(cs, rs, ct, rt) = cand_w;
      }
    
    }
    }

  }
  }
}


ParetoSet find_matches (
  const srvf::Srvf &Q1, const srvf::Srvf &Q2, 
  bool do_rots, size_t grid_width, size_t grid_height, 
  size_t nbuckets, double bucket_thresh )
{
  std::vector<double> tv1;
  std::vector<double> tv2;

  // Determine grid columns
  if (grid_width == 0)
    tv1 = Q1.params();
  else
    tv1 = srvf::util::linspace(Q1.domain_lb(), Q1.domain_ub(), grid_width);

  // Determine grid rows
  if (grid_height == 0)
    tv2 = Q2.params();
  else
    tv2 = srvf::util::linspace(Q2.domain_lb(), Q2.domain_ub(), grid_height);

  // Determine number of buckets in the Pareto set
  if (nbuckets == 0)
  {
    // Default: nbuckets = number of possible match lengths.
    size_t min_chunks = 2;
    size_t max_chunks = (tv1.size()-1) + (tv2.size()-1);
    nbuckets = max_chunks - min_chunks + 1;
  }

  ParetoSet S(nbuckets, bucket_thresh);

  if (do_rots == false)
  {
    MatchingGraph G = calculate_edge_weights(Q1, Q2, tv1, tv2);
    calculate_match_scores(G);

    for (size_t ct=1; ct<grid_width; ++ct)
    {
    for (size_t rt=1; rt<grid_height; ++rt)
    {

      for (size_t cs=0; cs<ct; ++cs)
      {
      for (size_t rs=0; rs<rt; ++rs)
      {
        double cur_dist = G(cs, rs, ct, rt);
        if (cur_dist < 1e5)
        {
          S.insert( PartialMatch( 
            tv1[cs], tv1[ct], 
            tv2[rs], tv2[rt],
            cur_dist ) );
        }
      }
      }
    
    }
    }
  }
  else
  {
    std::vector<Matrix> rotset = select_rotations(Q1, Q2, tv1, tv2);
    std::cout << rotset.size() << " rotations" << std::endl;

    for (size_t i=0; i<rotset.size(); ++i)
    {
      Srvf Q2r(Q2);
      Q2r.rotate(rotset[i]);

      MatchingGraph G = calculate_edge_weights(Q1, Q2r, tv1, tv2);
      calculate_match_scores(G);

      for (size_t ct=1; ct<grid_width; ++ct)
      {
      for (size_t rt=1; rt<grid_height; ++rt)
      {

        for (size_t cs=0; cs<ct; ++cs)
        {
        for (size_t rs=0; rs<rt; ++rs)
        {
          double cur_dist = G(cs, rs, ct, rt);
          if (cur_dist < 1e5)
          {
            S.insert( PartialMatch( 
              tv1[cs], tv1[ct], 
              tv2[rs], tv2[rt],
              cur_dist ) );
          }
        }
        }
      
      }
      }
    }
  }

  return S;
}


} // namespace pmatch
} // namespace srvf
