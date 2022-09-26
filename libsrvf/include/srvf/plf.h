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
#ifndef SRVF_PLF_H
#define SRVF_PLF_H 1

#include "pointset.h"
#include "matrix.h"
#include "util.h"

#include <cstddef>
#include <stdexcept>
#include <vector>


namespace srvf {


/**
 * Represents a piecewise-linear function.
 */
class Plf
{
public:

  /**
   * Creates an empty \c Plf.
   */
  Plf()
   : samps_(), params_(0)
  { }
  
  /**
   * Creates a new \c Plf with the given sample points.
   *
   * Corresponding parameter values will be uniformly spaced from 0 to 1.
   *
   * \param samples
   */
  Plf(const Pointset &samples)
    : samps_(samples)
  {
    params_=srvf::util::linspace(0.0,1.0,samples.npts());
  }

  /**
   * Creates a new \c Plf with the given sample points and parameters.
   *
   * \param samples
   * \param parameters
   */
  Plf(const Pointset &samples, const std::vector<double> &parameters)
    : samps_(samples), params_(parameters)
  {
    if (samples.npts()!=parameters.size())
      throw std::invalid_argument("number of samples != number of parameters");
  }

  /**
   * Returns a new \c Plf which is constant on the given interval.
   */
  static Plf create_constant(double a, double b, const Point &p)
  { 
    Plf res(Pointset(2,p), std::vector<double>(2));
    res.params_[0]=a; 
    res.params_[1]=b; 
    return res;
  }

  /**
   * Returns a new 1-D \c Plf which represents the identity function on 
   * the interval \f$ [a,b] \f$.
   */
  static Plf create_identity(double a, double b)
  {
    Plf res(Pointset(1,2), std::vector<double>(2));
    res.samps_[0][0] = a;
    res.samps_[1][0] = b;
    res.params_[0] = a;
    res.params_[1] = b;
    return res;
  }

  /** Returns the sample points. */
  Pointset &samps() { return samps_; }

  /** Returns the sample points. */
  const Pointset &samps() const { return samps_; }

  /** Returns the changepoint parameters. */
  std::vector<double> &params() { return params_; }

  /** Returns the changepoint parameters. */
  const std::vector<double> &params() const { return params_; }

  /** Returns the dimension of the ambient space. */
  size_t dim() const { return samps_.dim(); }

  /** Returns the number of changepoints. */
  size_t ncp() const { return samps_.npts(); }

  /** Does this \c Plf represent the empty map? */
  bool is_empty() const { return (samps_.npts() == 0); }

  /** Returns the left endpoint of the domain interval. */
  double domain_lb() const 
  { return (params_.size()>0 ? params_[0] : 0.0); };

  /** Returns the right endpoint of the domain interval. */
  double domain_ub() const 
  { return (params_.size()>0 ? params_[params_.size()-1] : 0.0); };
 
  Point evaluate(double t) const;
  Pointset evaluate(const std::vector<double> &tv) const;
  std::vector<double> preimages(const std::vector<double> &tv) const;

  double arc_length() const;
  Point centroid() const;
  std::vector<Point> bounding_box() const;

  void translate(const Point &v);
  void rotate(const Matrix &R);
  void scale(double s);

  void scale_to_unit_arc_length();
  void translate_to_origin();

  friend Plf linear_combination(const Plf &F1, const Plf &F2, 
                                double w1, double w2);
  friend Plf composition(const Plf &F1, const Plf &F2);
  friend Plf inverse(const Plf &F);
  friend Plf constant_speed_param(const Plf &F);
  friend Plf constant_speed_reparam(const Plf &F, double lb, double ub);

private:
  Pointset samps_;
  std::vector<double> params_;
};

Plf linear_combination(const Plf &F1, const Plf &F2, double w1, double w2);
Plf composition(const Plf &F1, const Plf &F2);
Plf inverse(const Plf &F);
Plf constant_speed_param(const Plf &F);
Plf constant_speed_reparam(const Plf &F, double lb, double ub);

} // namespace srvf

#endif // SRVF_PLF_H
