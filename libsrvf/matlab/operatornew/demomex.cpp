//
// new and delete C++ operators for MATLAB.
//
// DEMO FILE
//
// Petter Strandmark 2009
//

#include "mex.h"
#include "newdelete.h"



void mexFunction(int			nlhs, 		/* number of expected outputs */
				 mxArray		*plhs[],	/* mxArray output pointer array */
				 int			nrhs, 		/* number of inputs */
				 const mxArray	*prhs[]		/* mxArray input pointer array */)
{
	// input checks
	if (nrhs != 1 || nlhs != 1)
	{
		mexErrMsgTxt ("One output and one input, please.");
	}
	const mxArray *A = prhs[0];
	if (mxIsComplex(A))
	{
		mexErrMsgTxt ("Complex entries are not supported!");
	}
	
	// fetch its dimensions
	// actually, we must have m=n
	unsigned m = mxGetM(A);
	unsigned n = mxGetN(A);
	unsigned nzmax = mxGetNzmax(A);
	if (n < 1 || m < 1)
	{
		mexErrMsgTxt ("Matrix should not be empty");
	}
	
	double* data = mxGetPr(A);
	
	//Allocat data
	
	double* work = new double[m*n];
	
	for (int i = 0; i < m*n; i++) {
		work[i] = data[i];
	}
	
	
	
	for (int i = 0; i < m*n; i++) {
		work[i] *= 2;
	}
	
	

	//Create output
	plhs[0] = mxCreateDoubleMatrix(m, n, mxREAL);
	double* output = (double*)mxGetData(plhs[0]);
	for (int i = 0; i < m*n; i++) {
		output[i] = work[i];
	}
	
	delete[] work;

}

