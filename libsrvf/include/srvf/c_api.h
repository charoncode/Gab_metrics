/*
 * LibSRVF - a shape analysis library using the square root velocity framework.
 *
 * Copyright (C) 2014   FSU Statistical Shape Analysis and Modeling Group
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
#ifndef SRVF_C_API_H
#define SRVF_C_API_H 1

#include <cstdlib>

#define LIBSRVF_EXPORT extern "C"

typedef struct
{
  size_t rows;
  size_t cols;
  double *data;
} libsrvf_matrix_t;

LIBSRVF_EXPORT
libsrvf_matrix_t libsrvf_matrix_alloc(size_t rows, size_t cols);

LIBSRVF_EXPORT
void libsrvf_matrix_free(libsrvf_matrix_t x);

LIBSRVF_EXPORT
void libsrvf_array_copy(double *dest, const double *source, size_t m, size_t n, int transpose);

/**
 * samps is an MxN matrix representing M points in R^N
 * (one point per row, one component per column)
 *
 * params is a 1xM matrix, where samps(i,:) contains the function value at params(i)
 */
typedef struct
{
  libsrvf_matrix_t samps;
  libsrvf_matrix_t params;
} libsrvf_plf_t;

LIBSRVF_EXPORT
libsrvf_plf_t libsrvf_plf_alloc(size_t dim, size_t ncp);

LIBSRVF_EXPORT
void libsrvf_plf_free(libsrvf_plf_t f);

#define LIBSRVF_PLF_T_DIM(t) t.samps.cols
#define LIBSRVF_PLF_T_NCP(t) t.params.cols
#define LIBSRVF_PLF_T_SAMP_I_COMP_J(t,i,j) t.samps.data[(i)*t.samps.cols + (j)]

/**
 * samps is an MxN matrix representing M points in R^N
 * (one point per row, one component per column)
 *
 * params is a 1x(M+1) matrix, where samps(i,:) contains the function value
 * on the interval between params(i) and params(i+1)
 */
typedef struct
{
  libsrvf_matrix_t samps;
  libsrvf_matrix_t params;
} libsrvf_srvf_t;

LIBSRVF_EXPORT
libsrvf_srvf_t libsrvf_srvf_alloc(size_t dim, size_t ncp);

LIBSRVF_EXPORT
void libsrvf_srvf_free(libsrvf_srvf_t q);

#define LIBSRVF_SRVF_T_DIM(t) t.samps.cols
#define LIBSRVF_SRVF_T_NCP(t) t.params.cols
#define LIBSRVF_SRVF_T_SAMP_I_COMP_J(t,i,j) t.samps.data[(i)*t.samps.cols + (j)]

/**
 * Given SRVFs q1 and q2, find reparametrizations g1 and g2 such that q1*g1 and q2*g2 are 
 * optimally aligned.
 *
 * This function makes three assumptions about the input:
 * 1. Q1 and Q2 have unit norm.
 * 2. Q1 and Q2 have constant-speed parametrizations (i.e. |Qi(t)| = 1 for all t)
 * 3. The values of Q1 and Q2 alternate between 1 and -1 on adjacent intervals.
 *
 * inputs: two SRVFs with the same dimension
 * return: an array containing g1 and g2, in that order
 */
LIBSRVF_EXPORT
libsrvf_plf_t *libsrvf_fa_optimal_reparam(libsrvf_srvf_t q1, libsrvf_srvf_t q2);

/**
 * Finds reparametrizations that optimally align each SRVF in qs to qm.
 *
 * This function makes three assumptions about the input:
 * 1. All SRVFs have unit norm.
 * 2. All SRVFs have constant-speed parametrizations (i.e. |Qi(t)| = 1 for all t)
 * 3. The values of all SRVFs alternate between 1 and -1 on adjacent intervals.
 *
 * inputs: qm = the SRVF to which all of the others will be aligned
 *         qs = an array of SRVFs to be aligned to qm
 * return: an array of (nfuncs + 1) PLFs.  The first nfuncs PLFs are to be applied to 
 *         the SRVFs in qs, in order.  The last PLF is to be applied to qm.
 */
LIBSRVF_EXPORT
libsrvf_plf_t *libsrvf_fa_groupwise_reparam(libsrvf_srvf_t qm, libsrvf_srvf_t *qs, size_t nfuncs);

LIBSRVF_EXPORT libsrvf_srvf_t foo();
LIBSRVF_EXPORT libsrvf_matrix_t bar();

#endif // SRVF_C_API_H
