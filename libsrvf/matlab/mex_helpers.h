/*
 * libsrvf
 * =======
 *
 * A shape analysis library using the square root velocity framework.
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
#ifndef MEX_HELPERS_H
#define MEX_HELPERS_H 1

void die_if_not_valid_srvf_(const mxArray *samps, const mxArray *params)
{
  if (mxGetM(samps) != 1 || mxGetM(params) != 1)
  {
    mexErrMsgTxt("Sample and parameter matrices must have 1 row.\n");
  }
  if (mxGetN(samps)+1 != mxGetN(params))
  {
    mexErrMsgTxt("Parameter matrix must have 1 more column than sample matrix.\n");
  }
}

libsrvf_srvf_t mex_args_to_libsrvf_srvf_t_(const mxArray *samps, const mxArray *params)
{
  size_t dim = 1;
  size_t ncp = mxGetN(params);

  // Each sample point is stored contiguously, both in the incoming Matlab
  // matrices, and in the libsrvf matrices.
  libsrvf_srvf_t result = libsrvf_srvf_alloc(dim, ncp);
  libsrvf_array_copy(result.samps.data, mxGetPr(samps), dim, ncp - 1, 0);
  libsrvf_array_copy(result.params.data, mxGetPr(params), dim, ncp, 0);
  return result;
}

#endif // MEX_HELPERS_H
