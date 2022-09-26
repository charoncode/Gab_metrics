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
#ifndef SRVF_INTERP_H
#define SRVF_INTERP_H 1

#include "point.h"
#include "pointset.h"
#include "matrix.h"

#include <cstddef>
#include <vector>


namespace srvf
{

namespace interp
{


/**
 * 1-D table lookup.
 *
 * Given a table of numbers \c table[0]<=table[1]<=...<=table[n-1] and 
 * a number t, returns the smallest integer \c idx such that 
 * \c table[idx]<=t<table[idx+1].  If \c t<table[0], then the result is 
 * \c 0, and if \c t>=table[n-1], then the result is \c n-1.
 *
 * \param table a \c vector of \c n \c double's sorted in non-decreasing order
 * \param t the query
 * \result an index between \c -1 and \c n-1, inclusive
 */
inline int lookup(const std::vector<double> &table, double t)
{
  size_t n=table.size();

  if (n>=2)
  {
    if (t<table[0]) return 0;
    else if (t>=table[n-1]) return n-1;
    else
    {
      size_t i1=0, i3=n-2;
      size_t i2=(i1+i3)/2;
      while(i1<i3)
      {
        if (t<table[i2]) i3=i2;
        else if (t>=table[i2+1]) i1=i2+1;
        else break;

        i2=(i1+i3)/2;
      }
      return i2;
    }
  }
  else
  {
    return 0;
  }
}

/**
 * Lookup using 1-D \c Pointset as table.
 */
inline int lookup(const Pointset &table, double t)
{
  size_t n=table.npts();

  if (n>=2)
  {
    if (t<table[0][0]) return 0;
    else if (t>=table[n-1][0]) return n-1;
    else
    {
      size_t i1=0, i3=n-2;
      size_t i2=(i1+i3)/2;
      while(i1<i3)
      {
        if (t<table[i2][0]) i3=i2;
        else if (t>=table[i2+1][0]) i1=i2+1;
        else break;

        i2=(i1+i3)/2;
      }
      return i2;
    }
  }
  else
  {
    return 0;
  }
}


/**
 * Linear interpolation.
 *
 * If \c params(i)==params(i+1) for some \c i, the interpolant has a jump 
 * discontinuity.  In this case, we take the function to be right-continuous 
 * at that point.
 *
 * \param samps a \c Matrix containing the sample points, one point per column
 * \param params the parameter values corresponding to \a samps
 * \param tv parameters at which to interpolated.  Must be non-decreasing.
 * \return a \c Pointset containing the values at the specified abscissae
 */
inline Pointset 
interp_linear(const Pointset &samps, 
              const std::vector<double> &params, 
              const std::vector<double> &tv)
{
  Pointset result(samps.dim(),tv.size());

  size_t idx=0;
  if (params.size()>1 && tv[0]>params[1])
  {
    idx=lookup(params, tv[0]);
  }

  for (size_t i=0; i<tv.size(); ++i)
  {
    double tvi=tv[i];
    while (idx<params.size()-1 && tvi>=params[idx+1]) 
    {
      ++idx;
    }
    if (tvi<params[idx])
    {
      tvi=params[idx];
    }

    if (idx<params.size()-1)
    {
      double w1 = params[idx+1]-tvi;
      double w2 = tvi-params[idx];
      double w = w1+w2;

      if (w>1e-6)
      {
        w1 /= w;
        w2 /= w;
      }
      else
      {
        w1 = 0.0;
        w2 = 1.0;
      }

      weighted_sum(samps,samps,idx,idx+1,w1,w2,result,i);
    }
    else
    {
      weighted_sum(samps,samps,idx,idx,1.0,0.0,result,i);
    }
  }

  return result;
}

/**
 * Linear interpolation for \c Plf::preimages().
 */
inline std::vector<double> 
interp_linear (const std::vector<double> &samps, 
               const Pointset &params, 
               const std::vector<double> &tv)
{
  if (params.dim() != 1)
    throw std::invalid_argument("params must be a 1-D pointset");
  if (params.npts() != samps.size())
    throw std::invalid_argument("params and samps must have same size");

  std::vector<double> result(tv.size());

  size_t idx=0;
  if (params.npts()>1 && tv[0]>params[1][0])
  {
    idx=lookup(params, tv[0]);
  }

  for (size_t i=0; i<tv.size(); ++i)
  {
    double tvi=tv[i];
    while (idx<params.npts()-1 && tvi>=params[idx+1][0]) 
    {
      ++idx;
    }
    if (tvi<params[idx][0])
    {
      tvi=params[idx][0];
    }

    if (idx<params.npts()-1)
    {
      double w1 = params[idx+1][0]-tvi;
      double w2 = tvi-params[idx][0];
      double w = w1+w2;

      if (w>1e-6)
      {
        w1 /= w;
        w2 /= w;
      }
      else
      {
        w1 = 0.0;
        w2 = 1.0;
      }

      result[i] = w1*samps[idx] + w2*samps[idx+1];
    }
    else
    {
      result[i] = samps[idx];
    }
  }

  return result;
}

/**
 * Piecewise-constant interpolation.
 *
 * \param samps a \c Matrix containing the sample points, one point per column
 * \param params the parameter values corresponding to \a samps
 * \param tv parameters at which to interpolated.  Must be non-decreasing.
 * \param result [output] a \c Matrix to hold the result
 */
inline Pointset 
interp_const(const Pointset &samps, const std::vector<double> &params, 
             const std::vector<double> &tv)
{
  Pointset result(samps.dim(), tv.size());

  size_t idx=0;
  if (params.size()>1 && tv[0]>params[1])
  {
    idx=lookup(params, tv[0]);
  }

  for (size_t i=0; i<tv.size(); ++i)
  {
    double tvi=tv[i];
    while (idx<params.size()-1 && tvi>=params[idx+1]) 
    {
      ++idx;
    }
    if (tvi<params[idx])
    {
      tvi=params[idx];
    }
    if (idx>=params.size()-1)
    {
      --idx;
    }
    weighted_sum(samps,samps,idx,idx,1.0,0.0,result,i);
  }

  return result;
}

} // namespace srvf::interp

} // namespace srvf

#endif // SRVF_INTERP_H
