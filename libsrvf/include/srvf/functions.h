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
#ifndef SRVF_FUNCTIONS_H
#define SRVF_FUNCTIONS_H 1

#include "srvf.h"
#include "plf.h"

#include <cstddef>
#include <utility>
#include <vector>
#include <deque>
#include <map>


namespace srvf
{

namespace functions
{

typedef std::pair<size_t,size_t> match_vertex_t;


/**
 * Returns reparametrizations which optimally register \a Q1 and \a Q2.
 *
 * IMPORTANT:
 * \a Q1 and \a Q2 MUST be 1-D, constant-speed SRVFs in order for this 
 * function to work properly.  To get a constant-speed SRVF from an SRVF 
 * with arbitrary parametrization, use \c Srvf::constant_speed_param().
 *
 * Unlike \c srvf::opencurves::optimal_reparam(), which uses the older 
 * dynamic programming algorithm, this function implements the new 
 * specialized 1-D function matching.
 *
 * \param Q1 the SRVF of the first function
 * \param Q2 the SRVF of the second function
 * \return a \c vector containing two \c Plf's representing non-decreasing 
 *   functions \f$ \gamma_i:[0,1]\to[0,1] \f$.  When the first is applied 
 *   to \c Q1 and the second is applied to \c Q2 (using \c gamma_action()), 
 *   the result is a pair of representatives of the closed-up orbits of 
 *   \c Q1 and \c Q2 which are at minimum distance.
 */
std::vector<Plf> optimal_reparam(const Srvf &Q1, const Srvf &Q2);


/**
 * Builds the piecewise-linear reparametrizations which realize the optimal 
 * matching between \a Q1 and \a Q2 determined by \a path.
 */
std::vector<Plf> build_gammas(const Srvf&Q1, const Srvf &Q2, 
  const std::deque<match_vertex_t> &path);


/**
 * Find reparametrizations for the SRVFs in \a Mu and \a Qs which optimally 
 * match each function in \a Qs to \a Mu.
 *
 * IMPORTANT: \a Mu and all SRVFs in \a Qs MUST be 1-D, constant-speed SRVFs 
 * in order for this function to work properly.  To get a constant-speed 
 * SRVF from an SRVF with arbitrary parametrization, use 
 * \c Srvf::constant_speed_param().
 *
 * This function is mainly used as a step in the algorithm for computing 
 * Karcher means.  In particular, it is called by 
 * \c srvf::functions::karcher_mean().
 *
 * \param Mu the SRVF to which all functions will be matched
 * \param Qs contains the SRVFs which will be matched to \a Mu
 * \return a \c vector of \c Plf's representing the reparametrizations.  The 
 * first \c Qs.size() elements are to be applied to the elements of \a Qs, in 
 * the same order.  The last element is to be applied to \a Mu.
 */
std::vector<Plf> 
groupwise_optimal_reparam(const Srvf &Mu, const std::vector<Srvf> &Qs);


/**
 * Builds the piecewise-linear reparametrizations which realize the 
 * simultaneous optimal matchings between \a Mu and each SRVF in \a Qs, 
 * determined by the matching paths in \a paths.
 */
std::vector<Plf> groupwise_build_gammas(
  const Srvf &Mu, const std::vector<Srvf> &Qs, 
  const std::vector<std::deque<match_vertex_t> > &paths);


/**
 * Computes the Karcher mean of a collecton of 1-D SRVFs.
 *
 * Unlike \c srvf::opencurves::karcher_mean(), this function uses the 
 * new 1-D function matching algorithms.
 */
Srvf karcher_mean(const std::vector<Srvf> &Qs, 
                  double tol=1e-3, size_t max_iters=0);


/**
 * Don't use this class -- it's just here to give the unit tests 
 * access to internal routines, etc.
 */
class TestAccess
{
public:
  
  static double edge_variation(const Srvf &Q, size_t i1, size_t i2);
  static double edge_score(const Srvf &Q1, const Srvf &Q2, 
    size_t sc, size_t sr, size_t tc, size_t tr);
  static void calculate_scores(const Srvf &Q1, const Srvf &Q2, 
    std::map<match_vertex_t,double> &score, 
    std::map<match_vertex_t,match_vertex_t> &pred);
  static void build_gamma_segment(const Srvf &Q1, const Srvf &Q2, 
    size_t sc, size_t sr, size_t tc, size_t tr, 
    std::vector<double> &G1samps, std::vector<double> &G2samps);
  static double segment_translate(const Srvf &Q1, const Srvf &Q2, 
    size_t cs, size_t rs, size_t ct, size_t rt, double t);

private:
  
  TestAccess(){ }  // no instances
};

} // namespace srvf::functions

} // namespace srvf

#endif // SRVF_FUNCTIONS_H
