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
#ifndef SRVF_SRVF_H
#define SRVF_SRVF_H 1

#include "pointset.h"
#include "util.h"

#include <vector>
#include <cstddef>
#include <stdexcept>


namespace srvf {

class Matrix;
class Plf;

/**
 * Represents a piecewise-constant square-root velocity function.
 *
 * Any continuous piecewise-linear function \f$ f:[a,b]\to R^n \f$ will 
 * have a piecewise-constant square-root velocity function 
 * \f$ q\in L^2([a,b],R^n) \f$.  To represent such a function, we only 
 * have to keep track of the changepoint parameters 
 * \f[ t_0=a \le t_1 \le \ldots \le t_k=b \f]
 * and the \f$ k\f$ function values between the \f$ t_i \f$.
 *
 * Since these functions are piecewise-constant, things like norms, distances, 
 * and inner products in the \f$ L^2 \f$ space can be computed exactly.
 */
class Srvf
{
public:

  /**
   * Creates an empty \c Srvf.
   */
  Srvf() 
   : samps_(), params_(0)
  { }

  /**
   * Creates a new \c Srvf with the specified values.
   *
   * Changepoint parameters will be uniformly spaced from 0 to 1.
   * \param samples the sample points (will be deep-copied)
   */
  Srvf(const Pointset &samples)
   : samps_(samples)
  { 
    params_=srvf::util::linspace(0.0, 1.0, samples.npts()+1);
  }

  /**
   * Creates a new \c Srvf with the specified values and 
   * changepoint parameters.
   *
   * \param samples the sample points (will be deep-copied)
   * \param params the changepoint parameters (will be deep-copied).  Must 
   *        be a \c 1xN \c Matrix, where \c N=samples().cols()+1.
   */
  Srvf(const Pointset &samples, const std::vector<double> &parameters)
   : samps_(samples), params_(parameters)
  {
    if (parameters.size() != samples.npts()+1)
      throw std::invalid_argument("parameters has bad size");
  }

  /**
   * Creates a new \c Srvf which represents a constant function on the 
   * specified interval.
   *
   * \param a the left endpoint of the parameter interval
   * \param b the right endpoint of the parameter interval
   * \param val a \c vector<double> representing the value of the \c Srvf.  
   *   The dimension of the \c Srvf will be equal to \c val.size().
   */
  Srvf(double a, double b, std::vector<double> val)
   : samps_(val.size(), 1, val, Pointset::POINT_PER_ROW), params_(2)
  {
    params_[0]=a;
    params_[1]=b;
  }

  Point evaluate(double t) const;
  Pointset evaluate(const std::vector<double> &tv) const;

  /** Returns a reference to the samples. */
  Pointset &samps() { return samps_; }

  /** Returns a reference to the samples. */
  const Pointset &samps() const { return samps_; }

  /** Returns a reference to the parameters. */
  std::vector<double> &params() { return params_; }

  /** Returns a reference to the parameters. */
  const std::vector<double> &params() const { return params_; }

  /** Returns the dimension of the ambient space. */
  size_t dim() const { return samps_.dim(); }

  /** Returns the number of change points. */
  size_t ncp() const { return params_.size(); }

  /** Does this \c Srvf represent the empty map? */
  bool is_empty() const { return ncp()<2; }

  /** Returns the left endpoint of the domain interval. */
  double domain_lb() const 
  { return (params_.size()>0 ? params_[0] : 0.0); };

  /** Returns the right endpoint of the domain interval. */
  double domain_ub() const 
  { return (params_.size()>0 ? params_[params_.size()-1] : 0.0); };

  void rotate(const Matrix &R);
  void scale(double sf);
  void scale_to_unit_norm();

  friend double l2_norm(const Srvf &Q);
  friend double l2_product(const Srvf &Q1, const Srvf &Q2);
  friend double l2_distance(const Srvf &Q1, const Srvf &Q2);
  friend double sphere_distance(const Srvf &Q1, const Srvf &Q2);
  friend Srvf   linear_combination(const Srvf &Q1, const Srvf &Q2, 
                  double w1, double w2);
  friend Srvf   refinement(const Srvf &Q, const std::vector<double> &tv);
  friend Srvf   gamma_action(const Srvf &Q, const Plf &gamma);
  friend Srvf   constant_speed_param(const Srvf &Q);

private:
  Pointset samps_;
  std::vector<double> params_;
};

double l2_norm(const Srvf &Q);
double l2_product(const Srvf &Q1, const Srvf &Q2);
double l2_distance(const Srvf &Q1, const Srvf &Q2);
double sphere_distance(const Srvf &Q1, const Srvf &Q2);
Srvf   linear_combination(const Srvf &Q1, const Srvf &Q2, 
         double w1, double w2);
Srvf   refinement(const Srvf &Q, const std::vector<double> &tv);
Srvf   gamma_action(const Srvf &Q, const Plf &gamma);
Srvf   constant_speed_param(const Srvf &Q);

} // namespace srvf

#endif // SRVF_SRVF_H
