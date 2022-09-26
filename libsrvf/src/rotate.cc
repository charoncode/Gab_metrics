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
#include <srvf/rotate.h>
#include <srvf/interp.h>
#include <srvf/util.h>
#include <srvf/exceptions.h>

#ifdef USE_GSL
  #include <gsl/gsl_matrix.h>
  #include <gsl/gsl_linalg.h>
  #include <gsl/gsl_blas.h>
#endif

#include <algorithm>
#include <cstddef>
#include <cmath>
#include <stdexcept>


namespace srvf
{

#ifdef USE_GSL

// Computes the sign of the determinant of A
static int sgndet_(gsl_matrix *A)
{
  int signum;
  int status=0;
  int res=0;

  gsl_matrix *LU = gsl_matrix_alloc(A->size1, A->size2);
  gsl_permutation *perm = gsl_permutation_alloc(A->size1);
  if (!LU || !perm)
  {
    status=-1;
    goto cleanup;
  }

  gsl_matrix_memcpy(LU, A);
  gsl_linalg_LU_decomp(LU, perm, &signum);
  res = gsl_linalg_LU_sgndet(LU, signum);

cleanup:
  if (LU)   gsl_matrix_free(LU);
  if (perm) gsl_permutation_free(perm);
  if (status==-1)
  {
    throw std::bad_alloc();
  }

  return res;
}


// Returns a dim x dim matrix representing rotation by 
// theta radians about 
static Matrix perturbation_matrix_(double theta, size_t dim)
{
  double sinx, cosx;

  sinx = sin(theta);
  cosx = cos(theta);

  Matrix R = Matrix::identity(dim);
  R(0,0) = cosx; R(0,1) = -sinx;
  R(1,0) = sinx; R(1,1) = cosx;
  return R;
}


// Computes integral of Q1 * Q2^T
static void build_A_(const Srvf &Q1, const Srvf &Q2, gsl_matrix *A)
{
  std::vector<double> tv=srvf::util::unique(Q1.params(), Q2.params());
  size_t dim = Q1.dim();
  size_t npts = tv.size();

  Pointset Q1vals = Q1.evaluate(tv);
  Pointset Q2vals = Q2.evaluate(tv);

  gsl_matrix_set_zero(A);
  for (size_t i=0; i<npts-1; ++i)
  {
    double dt = tv[i+1] - tv[i];

    for (size_t j=0; j<dim; ++j)
    {
      for (size_t k=0; k<dim; ++k)
      {
        double vcur = gsl_matrix_get(A, j, k);
        vcur += Q1vals[i][j] * Q2vals[i][k] * dt;
        gsl_matrix_set(A, j, k, vcur);
      }
    }
  }
}


// Version of build_A_ for partial matches
static void build_A_(const Srvf &Q1, const Srvf &Q2, gsl_matrix *A, 
                     double a, double b, double c, double d, 
                     size_t Q1_start_idx, size_t Q2_start_idx)
{
  double m = (d-c) / (b-a);
  double rm = sqrt(m);
  
  size_t Q1_idx = Q1_start_idx;
  size_t Q2_idx = Q2_start_idx;
  double t1 = a;
  double t2 = c;

  gsl_matrix_set_zero(A);
  gsl_matrix *Ainc = gsl_matrix_alloc(A->size1, A->size2);
  if (!Ainc) throw std::bad_alloc();

  while (t1 < (b-1e-5) && t2 < (d-1e-5))
  {
    double dx1 = std::min(b, Q1.params()[Q1_idx+1]) - t1;
    double dy1 = m * dx1;
    double dy2 = std::min(d, Q2.params()[Q2_idx+1]) - t2;
    double dx2 = dy2 / m;
    
    gsl_matrix_set_zero(Ainc);
    for (size_t i=0; i<Q1.dim(); ++i)
    {
      for (size_t j=0; j<Q1.dim(); ++j)
      {
        double Qij = Q1.samps()[Q1_idx][i] * (rm * Q2.samps()[Q2_idx][j]);
        gsl_matrix_set(Ainc, i, j, Qij);
      }
    }

    if ( fabs(dx1 - dx2) < 1e-5 )
    {
      gsl_matrix_scale(Ainc, dx1);
      gsl_matrix_add(A, Ainc);
      t1 += dx1;
      t2 += dy1;
      ++Q1_idx;
      ++Q2_idx;
    }
    else if (dx1 < dx2)
    {
      gsl_matrix_scale(Ainc, dx1);
      gsl_matrix_add(A, Ainc);
      t1 += dx1;
      t2 += dy1;
      ++Q1_idx;
    }
    else
    {
      gsl_matrix_scale(Ainc, dx2);
      gsl_matrix_add(A, Ainc);
      t1 += dx2;
      t2 += dy2;
      ++Q2_idx;
    }
  }

  gsl_matrix_free(Ainc);
}


/**
 * Computes the optimal rotational alignment of \a Q2 to \a Q1.
 *
 * Returns a \c Matrix representing a rotation \f$ R \f$ which minimizes 
 * the L^2 norm of \f$ Q_1 - R Q_2 \f$.
 *
 * \a Q1 and \a Q2 must have the same dimension.  If the dimension 
 * is 1, the result will be the 1x1 identity matrix.
 */
Matrix optimal_rotation (const Srvf &Q1, const Srvf &Q2)
{
  int status=0;  // success
  int sgndet;
  size_t dim = Q1.dim();
  Matrix R;
  Matrix Q2pert;

  if (Q1.dim() != Q2.dim())
    throw std::invalid_argument("Q1 and Q2 must have same dimension.");

  if (Q1.dim() < 2)
  {
    return Matrix(1,1,1.0);
  }

  gsl_matrix *A    = gsl_matrix_alloc(dim, dim);
  gsl_matrix *V    = gsl_matrix_alloc(dim, dim);
  gsl_vector *S    = gsl_vector_alloc(dim);
  gsl_matrix *res  = gsl_matrix_alloc(dim, dim);
  gsl_vector *work = gsl_vector_alloc(dim);

  if ( !A || !V || !S || !res || !work )
  {
    status=-1; // GSL alloc error, so throw after cleanup
    goto cleanup;
  }

  // Compute A = \int_0^1 q1(t) q2(t)^T dt
  build_A_(Q1, Q2, A);

  // Compute sign of det(A)
  // If A is singular, apply a small rotation to Q2 and try again.
  sgndet = sgndet_(A);
  while(sgndet == 0){
    Matrix Rpert = perturbation_matrix_(0.001, dim);
    Srvf Q2pert(Q2);
    Q2pert.rotate(Rpert);
    build_A_(Q1, Q2pert, A);
    sgndet = sgndet_(A);
  }

  // Compute SVD of A
  if (gsl_linalg_SV_decomp(A, V, S, work) != GSL_SUCCESS)
  { 
    status=-3;  // SVD failed, so throw after cleanup
    goto cleanup;
  }
  
  // If det(A) < 0, change the sign of the last column of V
  if (sgndet < 0){
    for (size_t i=0; i<dim; i++)
    {
      double x = gsl_matrix_get(V, i, dim-1);
      gsl_matrix_set(V, i, dim-1, -x);
    }
  }

  // Rotation matrix is AV^T
  gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, A, V, 0.0, res);

  // Convert result to an srvf::Matrix.
  // gsl_matrix and srvf::Matrix both use row-major order.
  R=Matrix(dim, dim, res->data, Matrix::ROW_MAJOR);

cleanup:
  if ( A )    gsl_matrix_free(A);
  if ( V )    gsl_matrix_free(V);
  if ( S )    gsl_vector_free(S);
  if ( work ) gsl_vector_free(work);
  if ( res )  gsl_matrix_free(res);

  // status!=0 indicates an error
  if (status==-1)
    throw std::bad_alloc();
  else if (status==-2)
    throw AlgorithmFailure("A is singular.");
  else if (status==-3)
    throw AlgorithmFailure("SVD failed.");

  return R;
}


/**
 * Computes the optimal rotational alignment for a partial match 
 * of \a Q2 to \a Q1.
 *
 * Returns a \c Matrix representing a rotation \f$ R \f$ which minimizes 
 * the L^2 distance of \f$ Q_1([a,b]) \f$ from \f$ R Q_2([c,d]) \f$.
 *
 * \a Q1 and \a Q2 must have the same dimension.  If the dimension 
 * is 1, the result will be the 1x1 identity matrix.
 *
 * \param Q1 the first SRVF
 * \param Q2 the second SRVF
 * \param Q1_start_idx (optional) the index of the parameter subinterval 
 *        of \a Q1 which contains \a a.
 * \param Q2_start_idx (optional) the index of the parameter subinterval 
 *        of \a Q2 which contains \a c.
 */
Matrix optimal_rotation(const Srvf &Q1, const Srvf &Q2, 
  double a, double b, double c, double d, 
  size_t Q1_start_idx, size_t Q2_start_idx )
{
  if (Q1_start_idx >= Q1.ncp()){
    int i = srvf::interp::lookup(Q1.params(), a);
    Q1_start_idx = 
      (i>=0 ? ((size_t)i<Q1.ncp()-1 ? (size_t)i : Q1.ncp()-2) : 0);
  }
  if (Q2_start_idx >= Q2.ncp()){
    int i = srvf::interp::lookup(Q2.params(), c);
    Q2_start_idx = 
      (i>=0 ? ((size_t)i<Q2.ncp()-1 ? (size_t)i : Q2.ncp()-2) : 0);
  }

  // Copied and pasted from above (gulp)
  int status=0;  // success
  int sgndet;
  size_t dim = Q1.dim();
  Matrix R;
  Matrix Q2pert;

  if (Q1.dim() != Q2.dim())
    throw std::invalid_argument("Q1 and Q2 must have same dimension.");

  if (Q1.dim() < 2)
  {
    return Matrix(1,1,1.0);
  }

  gsl_matrix *A    = gsl_matrix_alloc(dim, dim);
  gsl_matrix *V    = gsl_matrix_alloc(dim, dim);
  gsl_vector *S    = gsl_vector_alloc(dim);
  gsl_matrix *res  = gsl_matrix_alloc(dim, dim);
  gsl_vector *work = gsl_vector_alloc(dim);

  if ( !A || !V || !S || !res || !work )
  {
    status=-1; // GSL alloc error, so throw after cleanup
    goto cleanup;
  }

  // Compute A = \int_a^b q1(t) ((q2*l)(t))^T dt
  // where l(t) is the linear function passing through the grid 
  // corners (a,c) and (b,d).
  build_A_(Q1, Q2, A, a, b, c, d, Q1_start_idx, Q2_start_idx);

  // Compute sign of det(A)
  // If A is singular, apply a small rotation to Q2 and try again.
  sgndet = sgndet_(A);
  while(sgndet == 0){
    Matrix Rpert = perturbation_matrix_(0.001, dim);
    Srvf Q2pert(Q2);
    Q2pert.rotate(Rpert);
    build_A_(Q1, Q2pert, A);
    sgndet = sgndet_(A);
  }

  // Compute SVD of A
  if (gsl_linalg_SV_decomp(A, V, S, work) != GSL_SUCCESS)
  { 
    status=-3;  // SVD failed, so throw after cleanup
    goto cleanup;
  }
  
  // If det(A) < 0, change the sign of the last column of V
  if (sgndet < 0){
    for (size_t i=0; i<dim; i++)
    {
      double x = gsl_matrix_get(V, i, dim-1);
      gsl_matrix_set(V, i, dim-1, -x);
    }
  }

  // Rotation matrix is AV^T
  gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, A, V, 0.0, res);

  // Convert result to an srvf::Matrix.
  // gsl_matrix and srvf::Matrix both use row-major order.
  R=Matrix(dim, dim, res->data, Matrix::ROW_MAJOR);

cleanup:
  if ( A )    gsl_matrix_free(A);
  if ( V )    gsl_matrix_free(V);
  if ( S )    gsl_vector_free(S);
  if ( work ) gsl_vector_free(work);
  if ( res )  gsl_matrix_free(res);

  // status!=0 indicates an error
  if (status==-1)
    throw std::bad_alloc();
  else if (status==-2)
    throw AlgorithmFailure("A is singular.");
  else if (status==-3)
    throw AlgorithmFailure("SVD failed.");

  return R;
}

#else // USE_GSL

Matrix optimal_rotation (const Srvf &Q1, const Srvf &Q2)
{
  throw UnsupportedOperation("optimal_rotation");
}

Matrix optimal_rotation(const Srvf &Q1, const Srvf &Q2, 
  double a, double b, double c, double d, 
  size_t Q1_start_idx, size_t Q2_start_idx )
{
  throw UnsupportedOperation("optimal_rotation");
}

#endif // USE_GSL

} // namespace srvf
