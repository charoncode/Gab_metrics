/*
 * libsrvf
 * =======
 *
 * A shape analysis library using the square root velocity framework.
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
#include <mex.h>
#include "operatornew/newdelete.h"

#include <vector>

#include <srvf/c_api.h>
#include "mex_helpers.h"

void do_usage()
{
  mexPrintf(
    "USAGE: [Gm,TGm,Gs,TGs] = %s(Qm,Tm,Qs,Ts)\n"
    "Inputs:\n"
    "\tQm, Tm : sample points and parameters of the mean SRVF\n"
    "\tQs, Ts : sample points and parameters of the other SRVFs\n"
    "Outputs:\n"
    "\tGm,TGm = the reparametrization for Qm\n",
    "\tGs,TGs = the reparametrizations for the Qs\n",
    mexFunctionName()
  );
}

static void teardown_matrix_(mxArray *matrix)
{
  if (matrix) mxDestroyArray(matrix);
}

static void teardown_cell_array_(mxArray *array)
{
  if (array)
  {
    size_t ncells = mxGetNumberOfElements(array);
    for (size_t i=0; i<ncells; ++i)
    {
      teardown_matrix_(mxGetCell(array, i));
    }
    mxDestroyArray(array);
  }
}

extern "C"
{
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  size_t nfuncs;
  mwSize Gs_dim;
  mxArray *Gi_matrix;
  mxArray *Ti_matrix;
  int success = 1;

  // Check arguments
  if (nrhs != 4 || nlhs != 4)
  {
    mexPrintf("error: incorrect number of arguments.\n");
    do_usage();
    return;
  }
  if (!mxIsCell(prhs[2]) || !mxIsCell(prhs[3]))
  {
    mexPrintf("error: Qs and Ts must be cell arrays.\n");
    do_usage();
    return;
  }
  nfuncs = mxGetNumberOfElements(prhs[2]);
  if (mxGetNumberOfElements(prhs[3]) != nfuncs)
  {
    mexPrintf("error: Qs and Ts must have the same number of elements.\n");
    do_usage();
    return;
  }

  // Check dimensions of all the sample point and parameter matrices
  die_if_not_valid_srvf_(prhs[0], prhs[1]);
  for (size_t i=0; i<nfuncs; ++i)
  {
    die_if_not_valid_srvf_(mxGetCell(prhs[2], i), mxGetCell(prhs[3], i));
  }

  libsrvf_srvf_t qm = mex_args_to_libsrvf_srvf_t_(prhs[0], prhs[1]);

  libsrvf_srvf_t *qs = new libsrvf_srvf_t[nfuncs];
  for (size_t i=0; i<nfuncs; ++i)
  {
    qs[i] = mex_args_to_libsrvf_srvf_t_(mxGetCell(prhs[2], i), mxGetCell(prhs[3], i));
  }

  // Compute the groupwise alignment using libsrvf
  libsrvf_plf_t *gs = libsrvf_fa_groupwise_reparam(qm, qs, nfuncs);

  // Allocate output variables
  Gs_dim = (mwSize)nfuncs;
  plhs[2] = mxCreateCellArray(1, &Gs_dim);
  plhs[3] = mxCreateCellArray(1, &Gs_dim);
  if (!plhs[2] || !plhs[3])
  {
    mexPrintf("error: mxCreateCellArray() failed.\n");
    success = 0;
    goto cleanup;
  }

  plhs[0] = mxCreateDoubleMatrix(1, LIBSRVF_PLF_T_NCP(gs[nfuncs]), mxREAL);
  plhs[1] = mxCreateDoubleMatrix(1, LIBSRVF_PLF_T_NCP(gs[nfuncs]), mxREAL);
  if (!plhs[0] || !plhs[1])
  {
    mexPrintf("error: mxCreateDoubleMatrix() failed.\n");
    success = 0;
    goto cleanup;
  }

  for (size_t i=0; i<nfuncs; ++i)
  {
    Gi_matrix = mxCreateDoubleMatrix(1, LIBSRVF_PLF_T_NCP(gs[i]), mxREAL);
    Ti_matrix = mxCreateDoubleMatrix(1, LIBSRVF_PLF_T_NCP(gs[i]), mxREAL);
    if (!Gi_matrix || !Ti_matrix)
    {
      mexPrintf("error: mxCreateDoubleMatrix() failed.\n");
      success = 0;
      goto cleanup;
    }
    mxSetCell(plhs[2], i, Gi_matrix);
    mxSetCell(plhs[3], i, Ti_matrix);
  }

  // Copy into plhs
  libsrvf_array_copy(mxGetPr(plhs[0]), gs[nfuncs].samps.data, 1, LIBSRVF_PLF_T_NCP(gs[nfuncs]), 0);
  libsrvf_array_copy(mxGetPr(plhs[1]), gs[nfuncs].params.data, 1, LIBSRVF_PLF_T_NCP(gs[nfuncs]), 0);

  for (size_t i=0; i<nfuncs; ++i)
  {
    libsrvf_array_copy(mxGetPr(mxGetCell(plhs[2], i)), gs[i].samps.data, 1, LIBSRVF_PLF_T_NCP(gs[i]), 0);
    libsrvf_array_copy(mxGetPr(mxGetCell(plhs[3], i)), gs[i].params.data, 1, LIBSRVF_PLF_T_NCP(gs[i]), 0);
  }

cleanup:

  if (!success)
  {
    teardown_matrix_(plhs[0]);
    teardown_matrix_(plhs[1]);
    teardown_cell_array_(plhs[2]);
    teardown_cell_array_(plhs[3]);
  }

  libsrvf_srvf_free(qm);
  if (qs)
  {
    for (size_t i=0; i<nfuncs; ++i)
    {
      libsrvf_srvf_free(qs[i]);
    }
    delete[] qs;
  }

  // We're responsible for freeing gs
  for (size_t i=0; i<nfuncs + 1; ++i)
  {
    libsrvf_plf_free(gs[i]);
  }
  delete[] gs;
}
} // extern "C"
