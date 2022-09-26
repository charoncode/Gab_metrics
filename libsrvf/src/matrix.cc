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
#include <srvf/matrix.h>

#include <cstring>
#include <stdexcept>


static const size_t SIZE_MAX_ = (size_t)-1;

namespace srvf{

#define __ELEMENTWISE_OP(A,B,C,op) \
  size_t nc=(A).size(); \
  for (size_t i=0; i<nc; ++i) { \
    (A)(i) = (B)(i) op (C)(i); \
  }

#define __SCALAR_OP(A,B,v,op) \
  size_t nc=(A).size(); \
  for (size_t i=0; i<nc; ++i) { \
    (A)(i) = (B)(i) op (v); \
  }
  
#define __MATRIX_MUL(R,A,B) \
  for (size_t i=0; i<(A).rows(); ++i) { \
    for (size_t j=0; j<(B).cols(); ++j) { \
      (R)(i,j)=0.0; \
      for (size_t k=0; k<(A).cols(); ++k) { \
        (R)(i,j) += (A)(i,k)*(B)(k,j); \
  } } } \


/** 
 * Creates an empty \c Matrix. 
 */
Matrix::Matrix () 
 : data_(0), rows_(0), cols_(0) 
{ }

/**
 * Create a new \c Matrix with the specified size.
 *
 * \param rows
 * \param cols
 */
Matrix::Matrix (size_t rows, size_t cols)
  : rows_(rows), cols_(cols)
{
  if (rows>0 && cols>0)
  {
    if (rows > SIZE_MAX_/cols)
      throw std::overflow_error("rows*cols exceeds SIZE_MAX");

    data_.resize(rows*cols);
  }
  else
  {
    rows_=0;
    cols_=0;
  }
}

/**
 * Create a new \c Matrix with the specified size, with all entries 
 * set to the given value.
 *
 * \param rows
 * \param cols
 * \param val
 */
Matrix::Matrix (size_t rows, size_t cols, double val)
  : rows_(rows), cols_(cols)
{
  if (rows>0 && cols>0)
  {
    if (rows > SIZE_MAX_/cols)
      throw std::overflow_error("rows*cols exceeds SIZE_MAX");

    size_t alloc_size=rows*cols;
    data_.resize(alloc_size);
    for (size_t i=0; i<alloc_size; ++i)
      data_[i]=val;
  }
  else
  {
    rows_=0;
    cols_=0;
  }
}

/**
 * Create a new \c Matrix initialized with the given data.
 *
 * A deep copy of the data is made.
 *
 * \param rows
 * \param cols
 * \param data pointer to a \c (rows)x(cols) array of \c doubles
 * \param layout \c ROW_MAJOR if \a data is in row major order, otherwise 
 *        \c COLUMN_MAJOR.  Default is \c ROW_MAJOR.
 */
Matrix::Matrix (size_t rows, size_t cols, const double *data, Majorness layout)
  : rows_(rows), cols_(cols)
{
  if (rows>0 && cols>0)
  {
    if (!data) 
      throw std::invalid_argument("data is null");
    if (rows > SIZE_MAX_/cols)
      throw std::overflow_error("rows*cols exceeds max element count");

    size_t alloc_size=rows*cols;
    if (layout==ROW_MAJOR)
    {
      data_.insert(data_.begin(), &(data[0]), &(data[alloc_size]));
    }
    else
    {
      data_.resize(alloc_size);
      for (size_t i=0; i<rows; ++i)
      {
        for (size_t j=0; j<cols; ++j)
        {
          data_[i*cols+j]=data[j*rows+i];
        }
      }
    }
  }
  else
  {
    rows_=0;
    cols_=0;
  }
}

/**
 * Create a new \c Matrix initialized with the given data.
 *
 * A deep copy of the data is made.
 *
 * \param rows
 * \param cols
 * \param data a \c vector containing the data
 * \param layout \c ROW_MAJOR if \a data is in row major order, otherwise 
 *        \c COLUMN_MAJOR.  Default is \c ROW_MAJOR.
 */
Matrix::Matrix (size_t rows, size_t cols, const std::vector<double> &data, 
                Majorness layout)
  : rows_(rows), cols_(cols)
{
  if (rows>0 && cols>0)
  {
    if (rows > SIZE_MAX_/cols)
      throw std::overflow_error("rows*cols exceeds SIZE_MAX");

    if (layout==ROW_MAJOR)
    {
      data_.insert(data_.begin(), data.begin(), data.end());
    }
    else
    {
      data_.resize(rows*cols);
      for (size_t i=0; i<rows; ++i)
      {
        for (size_t j=0; j<cols; ++j)
        {
          data_[i*cols+j]=data[j*rows+i];
        }
      }
    }
  }
  else
  {
    rows_=0;
    cols_=0;
  }
}

/**
 * Copy constructor.
 *
 * Creates a deep copy of the given \c Matrix.
 *
 * \param A
 */
Matrix::Matrix (const Matrix &A)
  : data_(A.data_), rows_(A.rows_), cols_(A.cols_)
{ }

/**
 * Assignment operator.
 *
 * Sets the current \c Matrix to a deep copy of the given \c Matrix.
 *
 * \param A
 */
Matrix& Matrix::operator= (const Matrix &A)
{
  if (this != &A)
  {
    data_ = A.data_;
    rows_ = A.rows_;
    cols_ = A.cols_;
  }
  return *this;
}

/**
 * Sets this \c Matrix to the empty matrix, freeing its resources.
 */
void Matrix::clear()
{
  data_.clear();
  rows_=0;
  cols_=0;
}

/**
 * Resizes this \c Matrix to have the specified size.
 *
 * If the matrix is non-empty, then any entries in the new matrix that 
 * were present in the old matrix will still have the same value as before 
 * the call to \c resize().  New entries are not initialized.
 *
 * If either of \a rows or \a cols is zero, the Matrix is cleared.
 *
 * \param rows the new number of rows
 * \param cols the new number of columns
 */
void Matrix::resize(size_t rows, size_t cols)
{
  if (rows==rows_ && cols==cols_)
  {
    return;
  }
  if (rows==0 || cols==0)
  {
    clear();
    return;
  }

  if (rows > SIZE_MAX_/cols)
    throw std::overflow_error("rows*cols exceeds max element count");

  data_.resize(rows*cols);
  rows_ = rows;
  cols_ = cols;
}

/**
 * Add \a A to this \c Matrix.
 *
 * \param A
 * \return a reference to this \c Matrix
 */
Matrix& Matrix::operator+= (const Matrix &A)
{
  if (rows_!=A.rows_ || cols_!=A.cols_)
    throw std::invalid_argument("size mismatch");

  __ELEMENTWISE_OP(*this,*this,A,+);
  return *this;
}

/**
 * Subtract \a A from this \c Matrix.
 *
 * \param A
 * \return a reference to this \c Matrix
 */
Matrix& Matrix::operator-= (const Matrix &A)
{
  if (rows_!=A.rows_ || cols_!=A.cols_)
    throw std::invalid_argument("size mismatch");

  __ELEMENTWISE_OP(*this,*this,A,-);
  return *this;
}

/**
 * Multiply this \c Matrix elementwise by \a A.
 *
 * \param A
 * \return a reference to this \c Matrix
 */
Matrix& Matrix::operator*= (const Matrix &A)
{
  if (rows_!=A.rows_ || cols_!=A.cols_)
    throw std::invalid_argument("size mismatch");

  __ELEMENTWISE_OP(*this,*this,A,*);
  return *this;
}

/**
 * Divide this \c Matrix elementwise by \a A.
 *
 * \param A
 * \return a reference to this \c Matrix
 */
Matrix& Matrix::operator/= (const Matrix &A)
{
  if (rows_!=A.rows_ || cols_!=A.cols_)
    throw std::invalid_argument("size mismatch");

  __ELEMENTWISE_OP(*this,*this,A,/);
  return *this;
}

/**
 * Add \a v to all elements of this \c Matrix.
 *
 * \param v
 * \return a reference to this \c Matrix
 */
Matrix& Matrix::operator+= (double v)
{
  __SCALAR_OP(*this,*this,v,+);
  return *this;
}

/**
 * Subtract \a v from all elements of this \c Matrix.
 *
 * \param v
 * \return a reference to this \c Matrix
 */
Matrix& Matrix::operator-= (double v)
{
  __SCALAR_OP(*this,*this,v,-);
  return *this;
}

/**
 * Multiply all elements of this \c Matrix by \a v.
 *
 * \param v
 * \return a reference to this \c Matrix
 */
Matrix& Matrix::operator*= (double v)
{
  __SCALAR_OP(*this,*this,v,*);
  return *this;
}

/**
 * Divide all elements of this \c Matrix by \a v.
 *
 * \param v
 * \return a reference to this \c Matrix
 */
Matrix& Matrix::operator/= (double v)
{
  __SCALAR_OP(*this,*this,v,/);
  return *this;
}

/**
 * Elementwise sum of two matrices.
 *
 * \param A
 * \param B
 * \return a new \c Matrix representing \c A+B
 */
Matrix operator+ (const Matrix &A, const Matrix &B)
{
  if (A.rows()!=B.rows() || A.cols()!=B.cols())
    throw std::invalid_argument("size mismatch");

  Matrix R(A.rows(),A.cols());
  __ELEMENTWISE_OP(R,A,B,+);
  return R;
}

/**
 * Elementwise difference of two matrices.
 *
 * \param A
 * \param B
 * \return a new \c Matrix representing \c A-B
 */
Matrix operator- (const Matrix &A, const Matrix &B)
{
  if (A.rows()!=B.rows() || A.cols()!=B.cols())
    throw std::invalid_argument("size mismatch");

  Matrix R(A.rows(),A.cols());
  __ELEMENTWISE_OP(R,A,B,-);
  return R;
}

/**
 * Elementwise product of two matrices.
 *
 * \param A
 * \param B
 * \return a new \c Matrix representing \c A.*B
 */
Matrix operator* (const Matrix &A, const Matrix &B)
{
  if (A.rows()!=B.rows() || A.cols()!=B.cols())
    throw std::invalid_argument("size mismatch");

  Matrix R(A.rows(),A.cols());
  __ELEMENTWISE_OP(R,A,B,*);
  return R;
}

/**
 * Elementwise division of two matrices.
 *
 * \param A
 * \param B
 * \return a new \c Matrix representing \c A.*B
 */
Matrix operator/ (const Matrix &A, const Matrix &B)
{
  if (A.rows()!=B.rows() || A.cols()!=B.cols())
    throw std::invalid_argument("size mismatch");

  Matrix R(A.rows(),A.cols());
  __ELEMENTWISE_OP(R,A,B,/);
  return R;
}

/**
 * Scalar addition.
 *
 * \param A
 * \param v
 */
Matrix operator+ (const Matrix &A, double v)
{
  Matrix R(A);
  R+=v;
  return R;
}

/**
 * Scalar subtraction.
 *
 * \param A
 * \param v
 */
Matrix operator- (const Matrix &A, double v)
{
  return A+(-v);
}

/**
 * Scalar multiplication.
 *
 * \param A
 * \param v
 */
Matrix operator* (const Matrix &A, double v)
{
  Matrix R(A);
  R*=v;
  return R;
}

/**
 * Scalar division.
 *
 * \param A
 * \param v
 */
Matrix operator/ (const Matrix &A, double v)
{
  Matrix R(A);
  R/=v;
  return R;
}

/**
 * Matrix product of two matrices.
 *
 * \param A
 * \param At \c true to multiply by the transpose of \a A
 * \param B
 * \param Bt \c true to multiply by the transpose of \a B
 * \return a new \c Matrix representing the matrix product
 */
Matrix product(const Matrix &A, const Matrix &B)
{
  if (A.cols() != B.rows()) 
    throw std::invalid_argument("size mismatch");

  Matrix R(A.rows(),B.cols());
  __MATRIX_MUL(R,A,B);
  return R;
}


/**
 * Returns the product of \a A and \a v.
 */
Point product(const Matrix &A, const Point &v)
{
  if (A.cols() != v.dim()) throw std::invalid_argument("size mismatch");

  Point res(A.rows(), 0.0);
  for (size_t i=0; i<A.rows(); ++i)
    for (size_t j=0; j<A.cols(); ++j)
      res[i] += A(i,j) * v[j];

  return res;
}


/**
 * Returns the product of \a v and \a A.
 */
Point product(const Point &v, const Matrix &A)
{
  if (A.rows() != v.dim()) throw std::invalid_argument("size mismatch");

  Point res(A.cols(), 0.0);
  for (size_t i=0; i<A.cols(); ++i)
    for (size_t j=0; j<A.rows(); ++j)
      res[i] += v[j] * A(j,i);

  return res;
}


/**
 * Returns a new \c Matrix representing the transpose of the given \c Matrix.
 *
 * \param A
 * \return the transpose of A
 */
Matrix transpose(const Matrix &A)
{
  Matrix At(A.cols(),A.rows());
  for (size_t i=0; i<A.rows(); ++i)
  {
    for (size_t j=0; j<A.cols(); ++j)
    {
      At(j,i)=A(i,j);
    }
  }
  return At;
}

} // namespace srvf
