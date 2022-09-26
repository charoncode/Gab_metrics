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
#ifndef PARETOSET_H
#define PARETOSET_H 1

#include <srvf/fileio.h>
#include <srvf/matrix.h>
#include <srvf/minmaxheap.h>

#include <algorithm>
#include <cmath>
#include <vector>
#include <limits>
#include <stdexcept>


#define DEFAULT_BUCKET_THRESH 0.05

namespace srvf
{

namespace pmatch
{

/**
 * Represents a partial match between two curves.
 *
 * Let \f$ \beta_1,\beta_2:[0,1]\to R^n \f$ be two curves parametrized 
 * by arc length.  An instance of \c PartialMatch specifies two subintervals 
 * \f$ [a,b] \f$ and \f$ [c,d] \f$, which together define a match between 
 * the partial curves \f$ \beta_1([a,b]) \f$ and \f$ \beta_2([c,d]) \f$.
 */
struct PartialMatch
{
  
  /**
   * Create a new \c PartialMatch representing a full match between 
   * two arc-length parametrized curves.
   */
  PartialMatch()
   : a(0.0), b(1.0), c(0.0), d(1.0), dist(0.0)
  { }

  /**
   * Create a new \c PartialMatch with the given interval endpoints.
   */
  PartialMatch(double a, double b, double c, double d, double dist)
   : a(a), b(b), c(c), d(d), dist(dist)
  { }

  /**
   * Returns the match length, defined as the sum of the lengths of 
   * the two intervals making up the partial match.
   */
  double length() const
  {
    return (b-a) + (d-c);
  }

  /**
   * Comparison operator, for ordered data structures.
   */
  bool operator< (const PartialMatch &M) const
  {
    return dist < M.dist;
  }

  double a;
  double b;
  double c;
  double d;
  double dist;
};


/**
 * Manages the collection of Pareto-optimal partial matches.
 */
class ParetoSet
{
public:

  /**
   * Creates a \c ParetoSet with the specified number of buckets.
   */
  ParetoSet(size_t nbuckets=1, double thresh=DEFAULT_BUCKET_THRESH)
   : buckets_(nbuckets), thresh_(thresh)
  { }

  /**
   * Inserts the given partial match into the Pareto set.
   *
   * The match will only be inserted if its shape distance field 
   * is within tolerance of the minimum distance for the corresponding 
   * match length.
   *
   * \return \c true if the match was inserted, \c false otherwise
   */
  bool insert(const PartialMatch &M)
  {
    size_t bi = length_to_bucket_idx(M.length());
    bool inserted = false;

    if (!buckets_[bi].empty())
    {
      double bucket_dist = buckets_[bi].min().dist;

      if (M.dist < bucket_dist + thresh_)
      {
        buckets_[bi].insert(M);

        // The new match may decrease the minimum distance for 
        // the bucket.  If so, we remove all matches in the bucket 
        // whose distances are above tolerance.
        bucket_dist = std::min(bucket_dist, M.dist);
        while (buckets_[bi].max().dist > bucket_dist + thresh_)
        {
          buckets_[bi].remove_max();
        }
        inserted = true;
      }
    }
    else
    {
      buckets_[bi].insert(M);
      inserted = true;
    }

    return inserted;
  }

  /**
   * Returns the number of buckets.
   */
  inline size_t nbuckets()
  {
    return buckets_.size();
  }

  /**
   * Returns the bucket index for the given match length.
   */
  inline size_t length_to_bucket_idx(double l)
  {
    return (size_t)(round(l * (buckets_.size()-1) / 2.0));
  }

  /**
   * Returns a reference to the \c vector containing the matches 
   * for the bucket at index \a i.
   */
  inline const std::vector<PartialMatch> &operator[] (size_t i)
  {
    if (i > buckets_.size())
      throw std::out_of_range("Index i is out of range.");

    return buckets_[i].data();
  }

  /**
   * Returns the total number of matches in the Pareto set.
   *
   * This is usually strictly larger than the number of buckets.
   */
  size_t total_size() const
  {
    size_t res=0;
    for (size_t i=0; i<buckets_.size(); ++i)
    {
      res += buckets_[i].size();
    }
    return res;
  }

  /**
   * Returns the weighted-norm Salukwadze distance of the Pareto set.
   *
   * The Salukwadze distance is defined as follows:
   *
   * \f[ d_s = \inf \lambda + w\epsilon \f]
   *
   * where the \f$ \inf \f$ runs over all possible values of the partiality 
   * measure \f$ \lambda \f$ and the shape distance \f$ \epsilon \f$.  
   * The partiality is defined as the maximum possible match length minus 
   * the length of the match.
   */
  double salukwadze_dist(double w) const
  {
    double res = std::numeric_limits<double>::max();
    for (size_t i=0; i<buckets_.size(); ++i)
    {
      if (buckets_[i].empty()) continue;

      double cur_lambda = 2.0 - buckets_[i][0].length();
      double cur_eps = buckets_[i].min().dist;
      double cur_nrm = cur_lambda + w*cur_eps;
      res = std::min(res, cur_nrm);
    }
    return res;
  }

  /**
   * Find the knee of the Pareto curve.
   *
   * \param thresh the maximum sum of squared residuals
   * \return the index of the last bucket for which the sum of squared 
   *   residuals was less than \a thresh
   */
  size_t find_knee(double thresh)
  {
    double nst = 0.0;
    double tdy = 0.0;

    for (size_t i=1; i<nbuckets(); ++i)
    {
      if (buckets_[i].empty()) continue;

      // Update dot product and squared norm
      double ti = buckets_[i][0].length();
      double yi = buckets_[i][0].dist;
      nst += ti*ti;
      tdy += ti*yi;
      if (nst < 1e-6) continue;

      // Get regression error up to and including bucket i
      double m = tdy / nst;
      double err = 0.0;
      for (size_t j=0; j<=i; ++j)
      {
        if (buckets_[j].empty()) continue;
        double tj = buckets_[j][0].length();
        double yj = buckets_[j][0].dist;
        double dy = yj - m*tj;
        err += dy*dy;
      }
      if (err > thresh)
      {
        return i-1;
      }
    }

    // If still here, then the regression error never got above thresh, 
    // so we return the index of the last bucket.
    return nbuckets()-1;
  }

  /**
   * Compute the linear regression errors for the Pareto curve.
   */
  std::vector<double> regression_errors()
  {
    std::vector<double> result(nbuckets(),0.0);
    double nst = 0.0;
    double tdy = 0.0;

    for (size_t i=1; i<nbuckets(); ++i)
    {
      if (buckets_[i].empty()) continue;

      // Update dot product and squared norm
      double ti = buckets_[i][0].length();
      double yi = buckets_[i][0].dist;
      nst += ti*ti;
      tdy += ti*yi;
      if (nst < 1e-6) continue;

      // Get regression error up to and including bucket i
      double m = tdy / nst;
      for (size_t j=0; j<=i; ++j)
      {
        if (buckets_[j].empty()) continue;
        double tj = buckets_[j][0].length();
        double yj = buckets_[j][0].dist;
        double dy = yj - m*tj;
        result[i] += dy*dy;
      }
    }
    
    return result;
  }

    
  /**
   * Save the matches to a .csv file.
   *
   * For each bucket (i.e. each possible range of match lengths), 
   * a matrix will be written to the output.  This matrix will have 
   * a row for each match.  Each row will have the form
   *
   *  \code[a b c d dist]
   *
   * where \c [a,b] and \c [c,d] are the intervals defining the partial 
   * match, and \c dist is the shape distance for the match.
   */
  void save_csv(std::ostream &os, char fieldsep=' ', char linesep='\n')
  {
    std::vector<srvf::Matrix> data;
    for (size_t i=0; i<buckets_.size(); ++i)
    {
      data.push_back(srvf::Matrix(buckets_[i].size(), 5));
      for (size_t j=0; j<buckets_[i].size(); ++j)
      {
        data.back()(j,0) = buckets_[i][j].a;
        data.back()(j,1) = buckets_[i][j].b;
        data.back()(j,2) = buckets_[i][j].c;
        data.back()(j,3) = buckets_[i][j].d;
        data.back()(j,4) = buckets_[i][j].dist;
      }
    }
    srvf::io::save_csv(os, data, fieldsep, linesep);
  }

private:

  std::vector<srvf::MinMaxHeap<PartialMatch> > buckets_;
  double thresh_;
};

} // namespace pmatch
} // namespace srvf

#endif // PARETOSET_H
