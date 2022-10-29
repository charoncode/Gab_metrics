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
#include <srvf/plf.h>
#include <srvf/srvf.h>
#include <srvf/qmap.h>
#include <srvf/pointset.h>
#include <srvf/partialmatch.h>
#include <srvf/util.h>
#include <srvf/interp.h>

#include <mex.h>
#include <vector>


void do_usage()
{
  mexPrintf(
    "USAGE: P = %s(X1,T1,Q2,T2,grid_cols,grid_rows,do_rots)\n"
    "Inputs:\n"
    "\tX1, T1 = sample points and parameters of first curve\n"
    "\tX2, T2 = sample points and parameters of second curve\n"
    "\tgrid_cols = number of possible breakpoints on first curve\n"
    "\tgrid_rows = number of possible breakpoints on second curve\n"
    "\tdo_rots : nonzero to optimize over rotations (default=0)\n"
    "Outputs:\n"
    "\tP = a Kx5 matrix representing the Pareto frontier.  Each row"
    " represents a partial match, and has the form\n"
    " [a b c d dist]\n"
    "where a and b are the indices of X1 where the match begins and ends, "
    "c and d are the indices of X2 where the match begins and ends, "
    "and dist gives the shape distance of the match\n",
    mexFunctionName()
  );
}


// Returns the index of the element of v which is closest to t.
static size_t find_nearest_(double t, const std::vector<double> &v)
{
  int idx = srvf::interp::lookup(v, t);
  if (idx < 0) return 0;

  size_t res = (size_t)idx;
  if (res+1 < v.size() && (fabs(t - v[res]) > fabs(v[res+1] - t)))
  {
    ++res;
  }
  return res;
}


extern "C"
{
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  // Check number of arguments
  if (nlhs != 1 || nrhs < 6 || nrhs > 7)
  { 
    mexPrintf("error: incorrect number of arguments\n");
    do_usage(); 
    return; 
  }

  // Check for incompatible SRVF dimensions
  size_t dim = mxGetM(prhs[0]);
  if (mxGetM(prhs[2]) != dim)
  { 
    mexPrintf("error: X1 and X2 have different dimensions\n"); 
    do_usage();
    return;
  }

  size_t npts1 = mxGetNumberOfElements(prhs[1]);
  size_t npts2 = mxGetNumberOfElements(prhs[3]);

  // Check for correct length on T1 and T2
  if ((npts1 != mxGetN(prhs[0])) || (npts2 != mxGetN(prhs[2])))
  {
    mexPrintf("error: incorrect length for T1 and/or T2\n");
    do_usage();
    return;
  }

  // Create sample pointsets.
  // Caller gives us the sample points in column-major, point-per-column 
  // ordering, which is the same thing as row-major, point-per-row ordering.
  srvf::Pointset samps1(dim, npts1, mxGetPr(prhs[0]), 
    srvf::Pointset::POINT_PER_ROW);
  srvf::Pointset samps2(dim, npts2, mxGetPr(prhs[2]), 
    srvf::Pointset::POINT_PER_ROW);

  // Create the parameter vectors T1 and T2
  double *T1_data = mxGetPr(prhs[1]);
  double *T2_data = mxGetPr(prhs[3]);
  std::vector<double> T1(&T1_data[0], &T1_data[npts1]);
  std::vector<double> T2(&T2_data[0], &T2_data[npts2]);

  // Get matching grid dimensions
  size_t grid_rows = mxGetScalar(prhs[4]);
  size_t grid_cols = mxGetScalar(prhs[5]);

  // Optimize over rotations?
  bool do_rots = false;
  if (nrhs > 6 && mxGetScalar(prhs[6]) != 0.0)
  {
    do_rots = true;
  }

  // Create Plfs
  srvf::Plf F1(samps1, T1);
  srvf::Plf F2(samps2, T2);

  // Scale F1 and F2 down by a common factor
  double L1 = F1.arc_length();
  double L2 = F2.arc_length();
  double scale_factor = std::min(L1,L2);
  if (scale_factor < 1e-6) scale_factor = 1.0;
  F1.scale( 1.0 / scale_factor );
  F2.scale( 1.0 / scale_factor );
  
  // Find the Pareto set
  srvf::Srvf Q1 = srvf::plf_to_srvf(F1);
  srvf::Srvf Q2 = srvf::plf_to_srvf(F2);
  srvf::pmatch::ParetoSet P = srvf::pmatch::find_matches (
    Q1, Q2, do_rots, grid_cols, grid_rows, (grid_rows+grid_cols)/2, 0.0005);
  size_t nmatches = P.total_size();

  // Allocate output variables
  plhs[0] = mxCreateDoubleMatrix(nmatches, 5, mxREAL);
  if (!plhs[0])
  {
    mexErrMsgIdAndTxt(mexFunctionName(), "mxCreateDoubleMatrix() failed.");
  }

  double *P_data = mxGetPr(plhs[0]);
  size_t Prow=0;
  for (size_t i=0; i<P.nbuckets(); ++i)
  {
    for (size_t j=0; j<P[i].size(); ++j)
    {
      // Get 1-based indices for match boundaries
      double ai = 1.0 + find_nearest_(P[i][j].a, T1);
      double bi = 1.0 + find_nearest_(P[i][j].b, T1);
      double ci = 1.0 + find_nearest_(P[i][j].c, T2);
      double di = 1.0 + find_nearest_(P[i][j].d, T2);

      P_data[Prow] = ai;
      P_data[nmatches+Prow] = bi;
      P_data[2*nmatches+Prow] = ci;
      P_data[3*nmatches+Prow] = di;
      P_data[4*nmatches+Prow] = P[i][j].dist;
      ++Prow;
    }
  }

}
} // extern "C"
