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
#include <srvf/c_api.h>
#include <srvf/pointset.h>
#include <srvf/plf.h>
#include <srvf/srvf.h>
#include <srvf/functions.h>

#include <vector>


LIBSRVF_EXPORT
libsrvf_matrix_t libsrvf_matrix_alloc(size_t rows, size_t cols)
{
  libsrvf_matrix_t x;
  x.rows = rows;
  x.cols = cols;
  x.data = new double[rows * cols];
  return x;
}

LIBSRVF_EXPORT
void libsrvf_matrix_free(libsrvf_matrix_t x)
{
  delete[] x.data;
}

LIBSRVF_EXPORT
void libsrvf_array_copy(double *dest, const double *source, size_t m, size_t n, int transpose)
{
  for (size_t i=0; i<m; ++i)
  {
    for (size_t j=0; j<n; ++j)
    {
      if (transpose) dest[j*m + i] = source[i*n + j];
      else dest[i*n + j] = source[i*n + j];
    }
  }
}

LIBSRVF_EXPORT
libsrvf_plf_t libsrvf_plf_alloc(size_t dim, size_t ncp)
{
  libsrvf_plf_t f;
  f.samps = libsrvf_matrix_alloc(ncp, dim);
  f.params = libsrvf_matrix_alloc(1, ncp);
  return f;
}

LIBSRVF_EXPORT
void libsrvf_plf_free(libsrvf_plf_t f)
{
  libsrvf_matrix_free(f.samps);
  libsrvf_matrix_free(f.params);
}

static void copy_plf_to_c_api_struct_(libsrvf_plf_t &dest, const srvf::Plf &source)
{
  size_t dim = source.dim();
  size_t ncp = source.ncp();

  LIBSRVF_PLF_T_DIM(dest) = dim;
  LIBSRVF_PLF_T_NCP(dest) = ncp;

  for (size_t i=0; i<ncp; ++i)
  {
    for (size_t j=0; j<dim; ++j)
    {
      LIBSRVF_PLF_T_SAMP_I_COMP_J(dest, i, j) = source.samps()[i][j];
    }
  }

  for (size_t i=0; i<ncp; ++i)
  {
    dest.params.data[i] = source.params()[i];
  }
}

LIBSRVF_EXPORT
libsrvf_srvf_t libsrvf_srvf_alloc(size_t dim, size_t ncp)
{
  libsrvf_srvf_t q;
  q.samps = libsrvf_matrix_alloc(ncp - 1, dim);
  q.params = libsrvf_matrix_alloc(1, ncp);
  return q;
}

LIBSRVF_EXPORT
void libsrvf_srvf_free(libsrvf_srvf_t q)
{
  libsrvf_matrix_free(q.samps);
  libsrvf_matrix_free(q.params);
}

static srvf::Srvf convert_srvf_from_c_api_struct_(const libsrvf_srvf_t &source)
{
  size_t dim = LIBSRVF_SRVF_T_DIM(source);
  size_t ncp = LIBSRVF_SRVF_T_NCP(source);

  srvf::Pointset samps = srvf::Pointset(dim, ncp - 1, source.samps.data);
  std::vector<double> params = std::vector<double> (source.params.data, source.params.data + ncp);
  srvf::Srvf result = srvf::Srvf(samps, params);

  return result;
}

LIBSRVF_EXPORT
libsrvf_plf_t* libsrvf_fa_optimal_reparam(libsrvf_srvf_t q1, libsrvf_srvf_t q2)
{
  srvf::Srvf Q1 = convert_srvf_from_c_api_struct_(q1);
  srvf::Srvf Q2 = convert_srvf_from_c_api_struct_(q2);

  std::vector<srvf::Plf> Gs = srvf::functions::optimal_reparam(Q1, Q2);

  libsrvf_plf_t* result = new libsrvf_plf_t[2];
  for (size_t i=0; i<2; ++i)
  {
    result[i] = libsrvf_plf_alloc(Gs[i].dim(), Gs[i].ncp());
    copy_plf_to_c_api_struct_(result[i], Gs[i]);
  }

  return result;
}

LIBSRVF_EXPORT
libsrvf_plf_t *libsrvf_fa_groupwise_reparam(libsrvf_srvf_t qm, libsrvf_srvf_t *qs, size_t nfuncs)
{
  srvf::Srvf Qm = convert_srvf_from_c_api_struct_(qm);
  std::vector<srvf::Srvf> Qs;
  for (size_t i=0; i<nfuncs; ++i)
  {
    Qs.push_back(convert_srvf_from_c_api_struct_(qs[i]));
  }

  std::vector<srvf::Plf> Gs = srvf::functions::groupwise_optimal_reparam(Qm, Qs);

  libsrvf_plf_t* result = new libsrvf_plf_t[nfuncs+1];

  result[nfuncs] = libsrvf_plf_alloc(Gs[nfuncs].dim(), Gs[nfuncs].ncp());
  copy_plf_to_c_api_struct_(result[nfuncs], Gs[nfuncs]);

  for (size_t i=0; i<nfuncs; ++i)
  {
    result[i] = libsrvf_plf_alloc(Gs[i].dim(), Gs[i].ncp());
    copy_plf_to_c_api_struct_(result[i], Gs[i]);
  }

  return result;
}
