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
#ifndef SRVF_MATRIX_H
#define SRVF_MATRIX_H 1

#include "point.h"

#include <cstddef>
#include <cmath>
#include <vector>


namespace srvf
{

/**
 * A basic matrix class.
 */
class Matrix
{
public:

  /**
   * Indicates either row-major or column-major ordering.
   * This is used when converting between \c Matrix objects and linear 
   * representations.  In row-major ordering, rows are stored contiguously 
   * in memory, and in column-major ordering, columns are stored 
   * contiguously.
   */
  enum Majorness
  {
    ROW_MAJOR=0,
    COLUMN_MAJOR=1
  };
  
  Matrix ();
  Matrix (size_t rows, size_t cols);
  Matrix (size_t rows, size_t cols, double val);
  Matrix (size_t rows, size_t cols, 
          const double *data, Majorness layout=ROW_MAJOR);
  Matrix (size_t rows, size_t cols, 
          const std::vector<double> &data, Majorness layout=ROW_MAJOR);
  Matrix (const Matrix &A);
  Matrix &operator= (const Matrix &A);

  /** Returns a new \c Matrix representing the identity matrix. */
  static Matrix identity(size_t dim)
  {
    Matrix I(dim, dim, 0.0);
    for (size_t i=0; i<dim; ++i)
    {
      I(i,i) = 1.0;
    }
    return I;
  }

  /** Returns the number of rows in this \c Matrix. */
  size_t rows() const 
  { return rows_; }

  /** Returns the number of columns in this \c Matrix. */
  size_t cols() const 
  { return cols_; }

  /** Returns the total number of entries in this \c Matrix. */
  size_t size() const 
  { return rows_*cols_; }

  /** 
   * Equality test.
   *
   * Returns true if and only if \a X is a \c Matrix having the same 
   * dimensions as this \c Matrix, with all corresponding entries equal.
   *
   * Elements are considered equal if they differ by less than \a thresh.
   */
  bool equals(const Matrix &X, double thresh=1e-9) const 
  {
    if (X.rows() != rows() || X.cols() != cols()) return false;
    for (size_t i=0; i<rows(); ++i)
    {
      for (size_t j=0; j<cols(); ++j)
      {
        if ( fabs((*this)(i,j) - X(i,j)) >= thresh )
        {
          return false;
        }
      }
    }
    return true;
  }

  void clear();
  void resize(size_t rows, size_t cols);
  
  /** 
   * Returns the \f$ n^{th} \f$ entry in the matrix.
   * Row-major ordering is used.  For example, in a \c 2x3 \c Matrix, 
   * calling this method with \c n=3 references the entry in row 1, column 0. 
   */ 
  double& operator[] (size_t n)
  { return data_[n]; }

  /** 
   * Returns the \f$ n^{th} \f$ entry in the matrix.
   * Row-major ordering is used.  For example, in a \c 2x3 \c Matrix, 
   * calling this method with \c n=3 references the entry in row 1, column 0. 
   */ 
  const double& operator[] (size_t n) const
  { return data_[n]; };

  /** 
   * Returns the \f$ n^{th} \f$ entry in the matrix.
   * Row-major ordering is used.  For example, in a \c 2x3 \c Matrix, 
   * calling this method with \c n=3 references the entry in row 1, column 0. 
   */ 
  double& operator() (size_t n)
  { return data_[n]; }

  /** 
   * Returns the \f$ n^{th} \f$ entry in the matrix.
   * Row-major ordering is used.  For example, in a \c 2x3 \c Matrix, 
   * calling this method with \c n=3 references the entry in row 1, column 0. 
   */ 
  const double& operator() (size_t n) const
  { return data_[n]; };

  /** Returns the specified entry. */
  double& operator() (size_t r, size_t c)
  { return data_[r*cols_+c]; }

  /** Returns the specified entry. */
  const double& operator() (size_t r, size_t c) const
  { return data_[r*cols_+c]; }

  /** Returns a pointer to the raw data for this \c Matrix. */
  std::vector<double> &data()
  { return data_; }

  /** Returns a \c const pointer to the raw data for this \c Matrix. */
  const std::vector<double> &data() const 
  { return data_; }

  /** Returns the norm of the matrix. */
  double norm() const
  {
    double result = 0.0;

    for (size_t i=0; i<size(); ++i)
    {
      double xi = (*this)[i];
      result += xi*xi;
    }

    return sqrt(result);
  }

  /** Returns the Euclidean distance from this \c Matrix to \a A. */
  double distance_to(const Matrix &A) const
  {
    double result = 0.0;

    for (size_t i=0; i<size(); ++i){
      double dxi = (*this)[i] - A[i];
      result += dxi*dxi;
    }

    return sqrt(result);
  }

  // In-place elementwise operations
  Matrix& operator+= (const Matrix &A);
  Matrix& operator-= (const Matrix &A);
  Matrix& operator*= (const Matrix &A);
  Matrix& operator/= (const Matrix &A);

  Matrix& operator+= (double v);
  Matrix& operator-= (double v);
  Matrix& operator*= (double v);
  Matrix& operator/= (double v);

  // Elementwise operations
  friend Matrix operator+ (const Matrix &A, const Matrix &B);
  friend Matrix operator- (const Matrix &A, const Matrix &B);
  friend Matrix operator* (const Matrix &A, const Matrix &B);
  friend Matrix operator/ (const Matrix &A, const Matrix &B);
  friend Matrix operator+ (const Matrix &A, double v);
  friend Matrix operator- (const Matrix &A, double v);
  friend Matrix operator* (const Matrix &A, double v);
  friend Matrix operator/ (const Matrix &A, double v);

  // Matrix multiplication
  friend Matrix product(const Matrix &A, const Matrix &B);
  friend Point  product(const Matrix &A, const Point &v);
  friend Point  product(const Point &v, const Matrix &A);

  // Related matrices
  friend Matrix transpose(const Matrix &A);

private:
  std::vector<double> data_;
  size_t rows_;
  size_t cols_;
};

Matrix operator+ (const Matrix &A, const Matrix &B);
Matrix operator- (const Matrix &A, const Matrix &B);
Matrix operator* (const Matrix &A, const Matrix &B);
Matrix operator/ (const Matrix &A, const Matrix &B);
Matrix operator+ (const Matrix &A, double v);
Matrix operator- (const Matrix &A, double v);
Matrix operator* (const Matrix &A, double v);
Matrix operator/ (const Matrix &A, double v);
Matrix product (const Matrix &A, const Matrix &B);
Point  product (const Matrix &A, const Point &v);
Point  product (const Point &v, const Matrix &A);
Matrix transpose(const Matrix &A);

} // namespace srvf


#endif // SRVF_MATRIX_H
