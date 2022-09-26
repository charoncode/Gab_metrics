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

#include <srvf/c_api.h>

void do_usage()
{
  mexPrintf(
    "USAGE: [G1,T1,G2,T2] = %s(Q1,T1,Q2,T2)\n"
    "Inputs:\n"
    "\tQ1, T1 : sample points and parameters of the first SRVF\n"
    "\tQ2, T2 : sample points and parameters of the second SRVF\n"
    "\tQ1 and Q2 are matrices with one ROW for each component function and\n"
    "\t\tone COLUMN for each sample point.\n"
    "\tIMPORTANT: Both SRVFs must have unit norm, and must have constant-speed parametrizations.\n"
    "Outputs:\n"
    "\tG1,T1 : the reparametrization for Q1\n"
    "\tG2,T2 : the reparametrization for Q2\n",
    mexFunctionName()
  );
}


extern "C"
{
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  size_t dim;
  size_t ncp1, ncp2;

  // Check arguments
  if (nrhs != 4 || nlhs != 4)
  {
    mexPrintf("error: incorrect number of arguments.\n");
    do_usage();
    return;
  }


  // Get dimension and number of points
  dim = mxGetM(prhs[0]);
  if (mxGetM(prhs[2]) != dim)
  {
    mexPrintf("error: Q1 and Q2 have different dimensions.\n");
    do_usage();
    return;
  }
  ncp1 = mxGetN(prhs[1]);
  ncp2 = mxGetN(prhs[3]);
  if (mxGetN(prhs[0])+1 != ncp1)
  {
    mexPrintf("error: bad length on T1 (must be size(Q1,2) + 1).\n");
    do_usage();
    return;
  }
  if (mxGetN(prhs[2])+1 != ncp2)
  {
    mexPrintf("error: bad length on T2 (must be size(Q2,2) + 1).\n");
    do_usage();
    return;
  }

  // Each sample point is stored contiguously, both in the incoming Matlab
  // matrices, and in the libsrvf matrices.
  libsrvf_srvf_t q1 = libsrvf_srvf_alloc(dim, ncp1);
  libsrvf_array_copy(q1.samps.data, mxGetPr(prhs[0]), dim, ncp1 - 1, 0);
  libsrvf_array_copy(q1.params.data, mxGetPr(prhs[1]), 1, ncp1, 0);

  libsrvf_srvf_t q2 = libsrvf_srvf_alloc(dim, ncp2);
  libsrvf_array_copy(q2.samps.data, mxGetPr(prhs[2]), dim, ncp2 - 1, 0);
  libsrvf_array_copy(q2.params.data, mxGetPr(prhs[3]), 1, ncp2, 0);


  // Get the reparametrizations for Q1 and Q2 using libsrvf
  libsrvf_plf_t *Gs = libsrvf_fa_optimal_reparam(q1, q2);


  size_t g1_ncp = LIBSRVF_PLF_T_NCP(Gs[0]);
  size_t g2_ncp = LIBSRVF_PLF_T_NCP(Gs[1]);


  // Allocate output variables and copy Gs[0] and Gs[1] into plhs
  plhs[0] = mxCreateDoubleMatrix(1, g1_ncp, mxREAL);
  plhs[1] = mxCreateDoubleMatrix(1, g1_ncp, mxREAL);
  plhs[2] = mxCreateDoubleMatrix(1, g2_ncp, mxREAL);
  plhs[3] = mxCreateDoubleMatrix(1, g2_ncp, mxREAL);
  if (!plhs[0] || !plhs[1] || !plhs[2] || !plhs[3])
  {
    mexPrintf("error: mxCreateDoubleMatrix() failed.\n");

    if (plhs[0]) mxDestroyArray(plhs[0]);
    if (plhs[1]) mxDestroyArray(plhs[1]);
    if (plhs[2]) mxDestroyArray(plhs[2]);
    if (plhs[3]) mxDestroyArray(plhs[3]);

    goto cleanup;
  }

  libsrvf_array_copy(mxGetPr(plhs[0]), Gs[0].samps.data, 1, g1_ncp, 0);
  libsrvf_array_copy(mxGetPr(plhs[1]), Gs[0].params.data, 1, g1_ncp, 0);

  libsrvf_array_copy(mxGetPr(plhs[2]), Gs[1].samps.data, 1, g2_ncp, 0);
  libsrvf_array_copy(mxGetPr(plhs[3]), Gs[1].params.data, 1, g2_ncp, 0);

cleanup:

  libsrvf_srvf_free(q1);
  libsrvf_srvf_free(q2);

  // We're responsible for freeing Gs
  libsrvf_plf_free(Gs[0]);
  libsrvf_plf_free(Gs[1]);
  delete[] Gs;
}
} // extern "C"
