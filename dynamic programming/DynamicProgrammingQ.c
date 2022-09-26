#include <math.h>
#include <stdlib.h>
#include "mex.h"

#define NNBRS   75

const int Nbrs[NNBRS][2] = { 
	{ 1, 1 }, 
	{ 1, 2 },
	{ 2, 1 },
	{ 2, 3 },
	{ 3, 2 },
	{ 1, 3 },
	{ 3, 1 },
	{ 1, 4 },
	{ 3, 4 },
	{ 4, 3 },
	{ 4, 1 },
	{ 1, 5 },
	{ 2, 5 },
	{ 3, 5 },
	{ 4, 5 },
	{ 5, 4 },
	{ 5, 3 },
	{ 5, 2 },
	{ 5, 1 },
	{ 1, 6 },
    { 2, 6 },
	{ 3, 6 },
	{ 4, 6 },
	{ 5, 6 },
	{ 6, 5 },
    { 6, 4 },
	{ 6, 3 },
	{ 6, 2 },
	{ 6, 1 },
	{ 7, 1 },
    { 7, 5 },
    { 7, 4 },
	{ 7, 3 },
	{ 7, 2 },
	{ 7, 6 },
	{ 1, 7 },
    { 2, 7 },
	{ 3, 7 },
	{ 4, 7 },
	{ 5, 7 },
	{ 6, 7 },
    { 1, 8 },
    { 2, 8 },
	{ 3, 8 },
	{ 4, 8 },
	{ 5, 8 },
	{ 6, 8 },
	{ 7, 8 },
	{ 8, 7 },
    { 8, 6 },
    { 8, 5 },
	{ 8, 4 },
	{ 8, 3 },
	{ 8, 2 },
	{ 8, 1 },
    { 1, 9 },
    { 2, 9 },
	{ 3, 9 },
	{ 4, 9 },
	{ 5, 9 },
	{ 6, 9 },
	{ 7, 9 },
    { 8, 9 },
    { 9, 8 },
	{ 9, 7 },
    { 9, 6 },
    { 9, 5 },
	{ 9, 4 },
	{ 9, 3 },
	{ 9, 2 },
	{ 9, 1 },
    { 1, 10 },
    { 2, 10 },
	{ 3, 10 },
	{ 4, 10 },
	{ 5, 10 },
	{ 6, 10 },
	{ 7, 10 },
    { 8, 10 },
    { 9, 10 },
    { 10, 9 }, 
    { 10, 8 },
	{ 10, 7 },
    { 10, 6 },
    { 10, 5 },
	{ 10, 4 },
	{ 10, 3 },
	{ 10, 2 },
	{ 10, 1 },
    { 1, 9 },
    { 2, 9 },
	{ 3, 9 },
	{ 4, 9 },
	{ 5, 9 },
	{ 6, 9 },
	{ 7, 9 },
    { 8, 9 },
    { 9, 8 },
	{ 9, 7 },
    { 9, 6 },
    { 9, 5 },
	{ 9, 4 },
	{ 9, 3 },
	{ 9, 2 },
	{ 9, 1 },
        { 1, 9 },
    { 2, 9 },
	{ 3, 9 },
	{ 4, 9 },
	{ 5, 9 },
	{ 6, 9 },
	{ 7, 9 },
    { 8, 9 },
    { 9, 8 },
	{ 9, 7 },
    { 9, 6 },
    { 9, 5 },
	{ 9, 4 },
	{ 9, 3 },
	{ 9, 2 },
	{ 9, 1 },
        { 1, 9 },
    { 2, 9 },
	{ 3, 9 },
	{ 4, 9 },
	{ 5, 9 },
	{ 6, 9 },
	{ 7, 9 },
    { 8, 9 },
    { 9, 8 },
	{ 9, 7 },
    { 9, 6 },
    { 9, 5 },
	{ 9, 4 },
	{ 9, 3 },
	{ 9, 2 },
	{ 9, 1 },
};

double CostFn2(const double *q1, const double *q2, const double *q2L, int k, int l, int i, int j, int n, int N, int M, double lam, double acons);
void thomas(double *x, const double *a, const double *b, double *c, int n);
void spline(double *D, const double *y, int n);
void lookupspline(double *t, int *k, double dist, double len, int n);
double evalspline(double t, const double D[2], const double y[2]);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
	int i, j, k, l, n, M, N, Eidx, Fidx, Ftmp, Fmin, Num, *Path, *x, *y, cnt;
	const double *q1, *q2;
	double *q2L, *yy, *D, *tmp, *E, Etmp, Emin, t, a, b, acons, lam = 0;

	if (nrhs != 5)
		mexErrMsgTxt("usage: [gam] = DynamicProgrammingQ(q1,q2,acons,(defunct)lam,(defunct)Disp)");

	if (!mxIsDouble(prhs[0]) || !mxIsDouble(prhs[1]) || !mxIsDouble(prhs[3]))
		mexErrMsgTxt("Expected double precision arguments.");

	if (mxGetNumberOfDimensions(prhs[0]) != 2 || mxGetNumberOfDimensions(prhs[1]) != 2)
		mexErrMsgTxt("First two arguments expected to be two dimensional.");

	n = mxGetM(prhs[0]);
	N = mxGetN(prhs[0]);
	
	if (n != mxGetM(prhs[1]) || N != mxGetN(prhs[1]))
		mexErrMsgTxt("Dimension mismatch between first and second argument.");

	if (nlhs > 1)
		mexErrMsgTxt("Expected one return.");

	plhs[0] = mxCreateDoubleMatrix(1, N, mxREAL);

	yy = mxGetPr(plhs[0]);
	q1 = mxGetPr(prhs[0]);
	q2 = mxGetPr(prhs[1]);
	acons = mxGetScalar(prhs[2]);


	
	M = 5*N;
	q2L = malloc(n*M*sizeof(double));

	D = malloc(2*N*sizeof(double));
	tmp = D + N;

	a = 1.0/N;
	b = 1.0;

	for (i = 0; i < n; ++i) {
		for (j = 0; j < N; ++j) {
			tmp[j] = q2[n*j + i];
		}

		spline(D, tmp, N);

		for (j = 0; j < M; ++j) {
			/* XXX: Extrapolation 1/M < 1/N */
			lookupspline(&t, &k, (j+1.0)/M - a, b - a, N);
			q2L[n*j + i] = evalspline(t, D+k, tmp+k);
		}
	}

	free(D);

	E = calloc(N*N, sizeof(double));
	Path = malloc(2*N*N*sizeof(int));

	for (i = 0; i < N; ++i) {
		E[N*i + 0] = 1;
		E[N*0 + i] = 1;
		Path[N*(N*0 + i) + 0] = -1;
		Path[N*(N*0 + 0) + i] = -1;
		Path[N*(N*1 + i) + 0] = -1;
		Path[N*(N*1 + 0) + i] = -1;
	}
	E[N*0 + 0] = 0;

	for (j = 1; j < N; ++j) {
		for (i = 1; i < N; ++i) {

			Emin = 100000;
			Eidx = 0;

			for (Num = 0; Num < NNBRS; ++Num) {
				k = i - Nbrs[Num][0];
				l = j - Nbrs[Num][1];

				if (k >= 0 && l >= 0) {
					Etmp = E[N*l + k] + CostFn2(q1,q2,q2L,k,l,i,j,n,N,M,lam,acons);
					if (Num == 0 || Etmp < Emin) {
                        Emin = Etmp;
						Eidx = Num;
					}
				}
			}
            
			E[N*j + i] = Emin;
			Path[N*(N*0 + j) + i] = i - Nbrs[Eidx][0];
			Path[N*(N*1 + j) + i] = j - Nbrs[Eidx][1];
		}
	}

	free(E);
	free(q2L);

	/* XXX: x, y assumed to be at most length N */
	x = malloc(2*N*sizeof(int));
	y = x + N;

	x[0] = N-1;
	y[0] = N-1;

	cnt = 1;
	while (x[cnt-1] > 0) {
		y[cnt] = Path[N*(N*0 + x[cnt-1]) + y[cnt-1]];
		x[cnt] = Path[N*(N*1 + x[cnt-1]) + y[cnt-1]];
		++cnt;
	}

	free(Path);
	for (i = 0, j = cnt-1; i < j; ++i, --j) {
		k = x[i];
		x[i] = x[j];
		x[j] = k;

		k = y[i];
		y[i] = y[j];
		y[j] = k;
	}

	for (i = 0; i < N; ++i) {

		Fmin = 100000;
		Fidx = 0;

		for (j = 0; j < cnt; ++j) {
			Ftmp = (int)fabs(i - x[j]);
			if (j == 0 || Ftmp < Fmin) {
				Fmin = Ftmp;
				Fidx = j;
			}
		}

		if (x[Fidx] == i) {
			yy[i] = (y[Fidx]+1);
		}
		else {
			if (x[Fidx] > i) {
				a = x[Fidx] - i;
				b = i - x[Fidx-1];
				yy[i] = (a*(y[Fidx-1]+1) + b*(y[Fidx]+1))/(a+b);
			}
			else {
				a = i - x[Fidx];
				b = x[Fidx+1] - i;
				yy[i] = (a*(y[Fidx+1]+1) + b*(y[Fidx]+1))/(a+b);
			}
		}

		yy[i] /= N;
	}

	free(x);
}

double CostFn2(const double *q1, const double *q2, const double *q2L, int k, int l, int i, int j, int n, int N, int M, double lam, const double acons) {
	double m = (j-l)/(double)(i-k), sqrtm = sqrt(m), E = 0, y, tmp, ip, fp, w, ww, normq1, normq2, scal, theta, calctmp;
	int x, idx, d;

	for (x = k; x <= i; ++x) {
		y = (x-k)*m + l + 1;
		fp = modf(y*M/N, &ip);
		idx = (int)(ip + (fp >= 0.5)) - 1;
        
        normq1 = 0;
        normq2 = 0;
        scal = 0;

        for (d = 0; d < n; ++d) {
            normq1 += q1[n*x + d]*q1[n*x + d];
            normq2 += q2L[n*idx + d]*q2L[n*idx + d];
            scal += q1[n*x + d]*sqrtm*q2L[n*idx + d];
        }
            
        normq1 = sqrt(normq1);
        normq2 = sqrtm*sqrt(normq2);

	    theta = acos(scal/(normq1*normq2));
        if (acons*theta<=M_PI/2) {
	        w = (normq1*normq2) * cos(acons*theta);
        }
        else{
            w = 0;
        }

        E += normq1*normq1 + normq2*normq2 - 2*w;
		

//         for (d = 0; d < n; ++d) {
// 			theta = acos((q1[n*x + d]*q2L[n*idx + d])/(sqrt(q1[n*x + d]*q1[n*x + d])*sqrt(q2L[n*idx + d]*q2L[n*idx + d])));
//             if (acons*theta<=M_PI/2) {
// 			w = q1[n*x + d]*sqrtm*q2L[n*idx + d] * cos(acons*theta);
//             }
//             else{
//             w = 0;
//             }
//             ww = q1[n*x + d]*sqrtm*q2L[n*idx + d];          
// 			E += q1[n*x + d]*q1[n*x + d]+sqrtm*q2L[n*idx + d]*sqrtm*q2L[n*idx + d] - 2*w;
// 		}
	}

	return E/N;
}

void thomas(double *x, const double *a, const double *b, double *c, int n) {
	double tmp;
	int i;

	c[0] /= b[0];
	x[0] /= b[0];

	for (i = 1; i < n; ++i) {
		tmp = 1/(b[i] - c[i-1] * a[i]);
		c[i] *= tmp;
		x[i] = (x[i] - x[i-1] * a[i])*tmp;
	}

	for (i = n-2; i >= 0; --i) {
		x[i] -= c[i]*x[i+1];
	}
}

void spline(double *D, const double *y, int n) {
	int i;
	double *a, *b, *c;

	a = malloc(3*n*sizeof(double));
	b = a + n;
	c = b + n;

	if (n < 4) {
		a[0] = 0;
		b[0] = 2;
		c[0] = 1;
		D[0] = 3*(y[1]-y[0]);

		a[n-1] = 1;
		b[n-1] = 2;
		c[n-1] = 0;
		D[n-1] = 3*(y[n-1]-y[n-2]);
	}
	else {
		a[0] = 0;
		b[0] = 2;
		c[0] = 4;
		D[0] = -5*y[0] + 4*y[1] + y[2];

		a[n-1] = 4;
		b[n-1] = 2;
		c[n-1] = 0;
		D[n-1] = 5*y[n-1] - 4*y[n-2] - y[n-3];
	}

	for (i = 1; i < n-1; ++i) {
		a[i] = 1;
		b[i] = 4;
		c[i] = 1;
		D[i] = 3*(y[i+1]-y[i-1]);
	}

	thomas(D, a, b, c, n);

	free(a);
}

void lookupspline(double *t, int *k, double dist, double len, int n) {
	*t = (n-1)*dist/len;
	*k = (int)floor(*t);

	*k = (*k > 0)*(*k);
	*k += (*k > n-2)*(n-2-*k);

	*t -= *k;
}

double evalspline(double t, const double D[2], const double y[2]) {
	double c[4];

	c[0] = y[0];
	c[1] = D[0];
	c[2] = 3*(y[1]-y[0])-2*D[0]-D[1];
	c[3] = 2*(y[0]-y[1])+D[0]+D[1];

	return t*(t*(t*c[3] + c[2]) + c[1]) + c[0];
}
