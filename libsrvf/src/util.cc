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
#include <srvf/util.h>

#include <cstdlib>
#include <time.h>
#include <algorithm>
#include <stdexcept>


namespace srvf
{

namespace util
{

std::vector<double> linspace(double a, double b, size_t n)
{
  std::vector<double> R(n);
  if (n>1)
  {
    R[0]=a;
    double dt=(b-a)/(double)(n-1);
    for (size_t i=1; i<n; ++i)
    {
      R[i]=R[i-1]+dt;
    }
  }
  else if (n>0)
  {
    R[0]=a;
  }
  return R;
}


std::vector<double> random_vector(size_t len, double first, int dir)
{
  //unsigned int seed = (unsigned int)time(NULL);
  //srandom(seed);
  std::vector<double> res(len);

  res[0] = first;
  for (size_t i=1; i<len; ++i)
  {
    if (dir > 0)
      res[i] = res[i-1] + ((double)rand()) / (RAND_MAX);
    else if (dir < 0)
      res[i] = res[i-1] - ((double)rand()) / (RAND_MAX);
    else
      res[i] = ((double)rand()) / (RAND_MAX);
  }
  return res;
}


/**
 * Returns a new \c vector containing unique elements of \a v1 and \a v2, 
 * sorted in ascending order.
 *
 * The result will be a \c vector with \c N elements, where \c N is the 
 * number of elements in the set-theoretic union of \c v1 and \c v2.
 *
 * For the purposes of this function, two \c doubles \c a and \c b are 
 * considered unique if \c fabs(a-b)>=thresh.
 *
 * \param v1
 * \param v2
 * \param thresh
 */
std::vector<double> unique (std::vector<double> v1, 
                            std::vector<double> v2, 
                            double thresh)
{
  std::vector<double> res;
  std::sort(v1.begin(),v1.end());
  std::sort(v2.begin(),v2.end());
  
  size_t i1=0;
  size_t i2=0;
  while (i1<v1.size() && i2<v2.size())
  {
    if (v1[i1] < v2[i2]-thresh)
    {
      res.push_back(v1[i1++]);
    }
    else if (v2[i2] < v1[i1]-thresh)
    {
      res.push_back(v2[i2++]);
    }
    else
    {
      res.push_back(v1[i1++]);
      ++i2;
    }
  }
  while (i1<v1.size()) res.push_back(v1[i1++]);
  while (i2<v2.size()) res.push_back(v2[i2++]);
  return res;
}

/**
 * Compute the first-order differences of the points in \a X.
 *
 * \param X a \c Pointset
 * \return a new \c Pointset containing the first-order differences of \a X
 */
Pointset diff(const Pointset &X)
{
  Pointset res(X.dim(),X.npts()-1);
  for (size_t i=0; i<res.npts(); ++i)
  {
    weighted_sum(X,X,i+1,i,1.0,-1.0,res,i);
  }
  return res;
}

/**
 * Forward difference approximation of first derivative.
 *
 * \param X a \c Pointset
 * \param tv abscissae corresponding to the points in \a X
 * \return a new \c Pointset containing the forward difference approximation
 */
Pointset diff(const Pointset &X, const std::vector<double> &tv)
{
  if (X.npts()!=tv.size())
    throw std::invalid_argument("X and tv must have the same length.");

  Pointset res(X.dim(),X.npts()-1);
  for (size_t i=0; i<res.npts(); ++i)
  {
    double dt=tv[i+1]-tv[i];
    weighted_sum(X,X,i+1,i,1.0,-1.0,res,i);
    res.scale(i,1.0/dt);
  }
  return res;
}


/**
 * Builds a map for looking up the elements of \a tv1 in \a tv2.
 *
 * \a tv1 and \tv2 must be sorted in increasing order, and all elements of 
 * \a tv1 must lie between the first and last elements of \a tv2.
 * If these conditions are not met, behavior is undefined.
 *
 * Given an index \c i, returns the index \c j such that 
 * \c tv2[j]<=tv1[i]<tv2[j+1] (the last inequality will be non-strict 
 * if \c tv1[i] is equal to the last element of \c tv2).
 */
std::map<size_t,size_t> build_lookup_map (
  const std::vector<double> &tv1, const std::vector<double> &tv2 )
{
  std::map<size_t,size_t> result;

  for (size_t tv1i=0, tv2i=0; tv1i+1<tv1.size(); ++tv1i)
  {
    result[tv1i] = tv2i;
    while (tv2i+2 < tv2.size() && tv1[tv1i+1] > (tv2[tv2i+1]-1e-5))
      ++tv2i;
  }
  result[tv1.size()-1] = tv2.size()-2;

  return result;
}

} // namespace srvf::util

} // namespace srvf
