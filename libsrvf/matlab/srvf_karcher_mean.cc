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
#include <srvf/srvf.h>
#include <srvf/reparam.h>
#include <srvf/rotate.h>
#include <srvf/opencurves.h>

#include <mex.h>
#include <vector>


void do_usage()
{
  mexPrintf(
    "USAGE: [Qm,Tm] = %s(Qs,Ts,tol,max_iters,do_rots,do_reparams)\n"
    "Inputs:\n"
    "\tQs, Ts : sample points and parameters of the SRVFs\n"
    "\ttol : stop when gradient norm is smaller than tol (default 1e-3)\n"
    "\tmax_iters : maximum number of gradient descent iterations (default is infinity)\n"
    "\tdo_rots : nonzero to rotationally align the SRVFs at every iteration (default is 1)\n"
    "\tdo_reparams : nonzero to reparameterize the SRVFs to match the mean at every iteration (default 1)\n"
    "Outputs:\n"
    "\tQm,Tm = sample points and parameters of the mean SRVF\n",
    mexFunctionName()
  );
}


extern "C"
{
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  size_t nfuncs;
  double tol = 1e-3;
  size_t max_iters = 0;
  bool do_rots = true;
  bool do_reparams = true;

  std::vector<srvf::Srvf> Qs;
  mxArray *sampsi_data;
  mxArray *paramsi_data;
  srvf::Pointset sampsi;
  std::vector<double> paramsi;
  size_t dim;

  srvf::Srvf Qm;
  double *Qm_data;
  double *Tm_data;


  // Check arguments
  if (nrhs < 2 || nrhs > 6 || nlhs != 2)
  {
    mexPrintf("error: incorrect number of arguments.\n");
    do_usage();
    return;
  }
  if (!mxIsCell(prhs[0]) || !mxIsCell(prhs[1]))
  {
    mexPrintf("error: Qs and Ts must be cell arrays.\n");
    do_usage();
    return;
  }
  nfuncs = mxGetNumberOfElements(prhs[0]);
  if (mxGetNumberOfElements(prhs[1]) != nfuncs)
  {
    mexPrintf("error: Qs and Ts must have the same number of elements.\n");
    do_usage();
    return;
  }
  if (nrhs > 2) tol = mxGetScalar(prhs[2]);
  if (nrhs > 3) max_iters = (size_t)mxGetScalar(prhs[3]);
  if (nrhs > 4) do_rots = (size_t)mxGetScalar(prhs[4]);
  if (nrhs > 5) do_reparams = (size_t)mxGetScalar(prhs[5]);


  // Empty collection?
  if (nfuncs == 0)
  {
    plhs[0] = mxCreateDoubleMatrix(0,0,mxREAL);
    plhs[1] = mxCreateDoubleMatrix(0,0,mxREAL);
    return;
  }


  // Get dimension
  dim = mxGetM(prhs[0]);


  // Create the SRVFs
  for (size_t i=0; i<nfuncs; ++i)
  {
    sampsi_data = mxGetCell(prhs[0],i);
    paramsi_data = mxGetCell(prhs[1],i);

    // More input checking
    if (mxGetM(sampsi_data) != dim)
    {
      mexPrintf("error: the SRVFs in Qs have different dimensions.\n");
      return;
    }
    if (mxGetN(paramsi_data) != mxGetN(sampsi_data)+1)
    {
      mexPrintf("error: Ts(%d) must have length Qs(%d)+1.\n", i, i);
      return;
    }

    sampsi = srvf::Pointset(dim, mxGetN(sampsi_data), mxGetPr(sampsi_data));
    paramsi = std::vector<double>(mxGetPr(paramsi_data), 
      mxGetPr(paramsi_data)+mxGetN(paramsi_data));

    Qs.push_back(srvf::Srvf(sampsi, paramsi));
  }


  // Compute the Karcher mean using libsrvf
  Qm = srvf::opencurves::karcher_mean(Qs, tol, max_iters, do_rots, do_reparams);


  // Allocate output variables
  plhs[0] = mxCreateDoubleMatrix(dim, Qm.samps().npts(), mxREAL);
  plhs[1] = mxCreateDoubleMatrix(1, Qm.ncp(), mxREAL);
  if (!plhs[0] || !plhs[1])
  {
    mexPrintf("error: mxCreateDoubleMatrix() failed.\n");
    if (plhs[0]) mxDestroyArray(plhs[0]);
    if (plhs[1]) mxDestroyArray(plhs[1]);
    return;
  }

  
  // Copy samples and parameters of Qm into plhs[0] and plhs[1]
  Qm_data = mxGetPr(plhs[0]);
  Tm_data = mxGetPr(plhs[1]);
  for (size_t i=0; i<Qm.samps().npts(); ++i)
  {
    for (size_t j=0; j<dim; ++j)
    {
      Qm_data[i*dim+j] = Qm.samps()[i][j];
    }
  }
  for (size_t i=0; i<Qm.params().size(); ++i)
  {
    Tm_data[i] = Qm.params()[i];
  }

}
} // extern "C"
