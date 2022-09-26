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
#ifndef SRVF_POINT_H
#define SRVF_POINT_H 1

#include "numeric.h"

#include <cstddef>
#include <cmath>
#include <vector>
#include <iterator>
#include <stdexcept>


namespace srvf
{

/**
 * A point class supporting arbitrary dimensions.
 */
class Point
{
public:
  Point(){ }

  /**
   * Create an \a n dimensional point.
   *
   * \param n the dimension of the point to be created
   * \param v all components of the point will be initialized to this value
   *   (default is 0)
   */
  Point(size_t n, double v = 0.0)
   : data_(n, v)
  { }

  /**
   * Create a point and initialize its components with the given 
   * range of data.
   *
   * Data from the range \c [first,last) is copied into the newly-created 
   * point.
   *
   * \param first an \c input_iterator defining the beginning of the range
   * \param last an \c input_iterator defining the end of the range
   */
  template<typename _InputIterator>
  Point(_InputIterator first, _InputIterator last)
   : data_(first, last) 
  { }

  /**
   * Return a reference to the nth component of this point.
   */
  inline double &operator[] (size_t n)
  { 
    if (n >= dim()) throw std::out_of_range("Component index out of range");
    return data_[n];
  }

  /**
   * Return a \c const reference to the nth component of this point.
   */
  const inline double &operator[] (size_t n) const
  { 
    if (n >= dim()) throw std::out_of_range("Component index out of range");
    return data_[n];
  }

  /** Return the dimension of this point. */
  inline size_t dim() const { return data_.size(); }

  /**
   * Adds the specified point \a P to this one.
   * 
   * \param P the point to be added to this one
   */
  inline Point &operator+= (const Point &P)
  {
    if (dim() != P.dim()) throw std::invalid_argument("Dimension mismatch");

    for (size_t i=0; i<dim(); ++i)
    {
      data_[i] += P.data_[i];
    }
    return *this;
  }

  /**
   * Returns the sum of this point and \a P.
   *
   * \param P the point to be added to this one
   */
  inline Point operator+ (const Point &P) const
  {
    if (dim() != P.dim()) throw std::invalid_argument("Dimension mismatch");

    Point res(*this);
    res += P;
    return res;
  }

  /**
   * Subtracts the specified point \a P from this one.
   *
   * \param P the point to be subtracted from this one
   */
  inline Point &operator-= (const Point &P)
  {
    if (dim() != P.dim()) throw std::invalid_argument("Dimension mismatch");

    for (size_t i=0; i<dim(); ++i)
    {
      data_[i] -= P.data_[i];
    }
    return *this;
  }

  /**
   * Returns the difference of this point and \a P.
   *
   * \param P the point to be subtracted from this one
   */
  inline Point operator- (const Point &P) const
  {
    if (dim() != P.dim()) throw std::invalid_argument("Dimension mismatch");

    Point res(*this);
    res -= P;
    return res;
  }

  /**
   * Multiplies this point by the given scalar.
   *
   * \param s the scalar by which this point will be multiplied
   */
  inline Point &operator*= (double s)
  {
    for (size_t i=0; i<dim(); ++i)
    {
      data_[i] *= s;
    }
    return *this;
  }

  /**
   * Returns the scalar product of this point with the given number.
   *
   * \param s the scalar by which this point will be multiplied
   */
  inline Point operator* (double s) const
  {
    Point res(*this);
    for (size_t i=0; i<dim(); ++i)
    {
      res.data_[i] *= s;
    }
    return res;
  }

  /**
   * Divides this point by the given scalar.
   *
   * \param s the scalar by which this point will be divided
   */
  inline Point &operator/= (double s)
  {
    for (size_t i=0; i<dim(); ++i)
    {
      data_[i] /= s;
    }
    return *this;
  }

  /**
   * Returns the scalar product of this point with the reciprocal 
   * of the given number.
   *
   * \param s the scalar by which this point will be divided
   */
  inline Point operator/ (double s) const
  {
    Point res(*this);
    for (size_t i=0; i<dim(); ++i)
    {
      res.data_[i] /= s;
    }
    return res;
  }

  /**
   * Returns the Euclidean norm of this point.
   */
  inline double norm() const
  {
    double res=0.0;
    for (size_t i=0; i<dim(); ++i)
    {
      res += data_[i] * data_[i];
    }
    return sqrt(res);
  }

  /**
   * Returns the dot product of this point with \a P.
   *
   * \param P the other point
   */
  inline double dot_product(const Point &P) const
  {
    if (dim() != P.dim()) throw std::invalid_argument("Dimension mismatch");

    double res=0.0;
    for (size_t i=0; i<dim(); ++i)
    {
      res += data_[i] * P.data_[i];
    }
    return res;
  }

  /**
   * Returns the Euclidean distance from this point to \a P.
   *
   * The points must have the same dimension.
   *
   * \param P the other point
   */
  inline double distance_to(const Point &P) const
  {
    if (dim() != P.dim()) throw std::invalid_argument("Dimension mismatch");

    double res=0.0;
    for (size_t i=0; i<dim(); ++i)
    {
      double dpi = (*this)[i] - P[i];
      res += dpi*dpi;
    }
    return sqrt(res);
  }

  /**
   * Equality operator.
   *
   * Returns true if and only if 
   *  1) both points have the same dimension, and 
   *  2) corresponding components are equal, as determined by calling 
   *     \c almost_equal() with default thresholds.
   *
   * \param P the point to which this one will be compared
   */
  inline bool operator== (const Point &P) const
  { 
    if (dim() != P.dim()) return false;

    for (size_t i=0; i<dim(); ++i)
      if ( !srvf::numeric::almost_equal(data_[i], P.data_[i]) ) return false;

    return true;
  }

  /**
   * Lexicographic (dictionary) ordering on points.
   *
   * \param P the point to which this one will be compared
   */
  inline bool operator< (const Point &P) const
  { 
    size_t i;
    for (i=0; i<dim(); ++i)
    {
      if (i >= P.dim()) return false;
      else if ((*this)[i] < P[i]) return true;
      else if ((*this)[i] > P[i]) return false;
    }
    return (i < P.dim());
  }

  /**
   * Determine whether or not \c P lies on the same ray based at the origin
   * as this point.
   *
   * If the points have different dimensions, the result will be \c false.
   * If either point is zero (i.e. has norm less than \a tol), the result 
   * will be \c true.  Otherwise, we project both points radially onto 
   * the unit sphere and then compute their distance.  We return true 
   * if and only if this distance is less than \a tol.
   */
  inline bool on_same_ray(const Point &P, double tol=1e-4) const
  {
    // Points that have different dimensions aren't on the same ray
    if (P.dim() != dim()) return false;

    double n1 = this->norm();
    double n2 = P.norm();

    // Zero is on the same ray as any point
    if (n1 < tol || n2 < tol) return true;

    // If neither point is zero, compute the Euclidean distance between the 
    // radial projections of the points on the unit sphere.
    double d = 0.0;
    for (size_t i=0; i<dim(); ++i)
    {
      double di = fabs((*this)[i]/n1 - P[i]/n2);
      d += di*di;
    }
    d = sqrt(d);
    return (d < tol);
  }

  friend inline Point linear_combination(const Point &A, const Point &B, 
    double w1, double w2);

private:
  std::vector<double> data_;
};


inline Point linear_combination(const Point &A, const Point &B, 
  double w1, double w2)
{
  Point res(A.dim());
  for (size_t i=0; i<A.dim(); ++i)
  {
    res[i] = w1*A[i] + w2*B[i];
  }
  return res;
}

} // namespace srvf

#endif // SRVF_POINT_H
