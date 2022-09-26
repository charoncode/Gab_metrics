#include <stdlib.h>
#include <time.h>
#include <CUnit/CUnit.h>

#include "dp_grid.h"


#if RAND_MAX < 100
#error "RAND_MAX < 100!!!"
#endif

static void random_partition( double *T, int n )
{
  double s;
  int i;

  T[0] = 0.0;
  for ( i=1; i<n; ++i )
  {
    T[i] = T[i-1] + (double)(rand() % 100);
  }

  s = T[n-1];
  for ( i=1; i<n; ++i )
  {
    T[i] /= s;
  }
}


static void test_dp_lookup_basic1()
{
  double T[] = { 0.0, 0.1, 0.4, 0.6, 1.0 };
  double params[] = { 
    0.0, 0.05, 0.09, 0.0999, 0.1, 0.11, 0.35, 0.3999, 0.4, 
    0.5, 0.6, 0.7, 0.75, 0.89, 0.9999, 1.0 };
  int indexes[] = {
    0, 0, 0, 0, 1, 1, 1, 1, 2, 
    2, 3, 3, 3, 3, 3, 3 };
  int nsamps, nparams;
  int i, idx;

  nsamps = sizeof(T) / sizeof(double);
  nparams = sizeof(params) / sizeof(double);

  for ( i=0; i<nparams; ++i )
  {
    idx = dp_lookup( T, nsamps, params[i] );
    CU_ASSERT_EQUAL( idx, indexes[i] );

    /*
    if ( idx < nsamps-2 )
      printf( "%f in [%f,%f)\n", params[i], T[idx], T[idx+1] );
    else
      printf( "%f in [%f,%f]\n", params[i], T[idx], T[idx+1] );
    */
  }
}


static void test_dp_lookup_random1()
{
  double T[10];
  int nsamps;
  int ntrials = 1000;
  double t;
  int i, idx;

  srand( time(NULL) );

  nsamps = sizeof(T) / sizeof(double);
  random_partition( T, nsamps );

  for ( i=0; i<ntrials; ++i )
  {
    t = (rand() % 100) / 99.0;
    idx = dp_lookup( T, nsamps, t );
    CU_ASSERT( (t>0.9999 && idx==nsamps-2) || (t >= T[idx] && t < T[idx+1]) );

    /*
    if ( idx < nsamps-2 )
      printf( "%f in [%f,%f)\n", t, T[idx], T[idx+1] );
    else
      printf( "%f in [%f,%f]\n", t, T[idx], T[idx+1] );
    */
  }
}

static void test_dp_lookup_random2()
{
  double T[256];
  int nsamps;
  int ntrials = 10000;
  double t;
  int i, idx;

  srand( time(NULL) );

  nsamps = sizeof(T) / sizeof(double);
  random_partition( T, nsamps );

  for ( i=0; i<ntrials; ++i )
  {
    t = (rand() % 100) / 99.0;
    idx = dp_lookup( T, nsamps, t );
    CU_ASSERT(idx<nsamps-1);
    CU_ASSERT( (t>0.9999 && idx==nsamps-2) || (t >= T[idx] && t < T[idx+1]) );

    /*
    if ( idx < nsamps-2 )
      printf( "%f in [%f,%f)\n", t, T[idx], T[idx+1] );
    else
      printf( "%f in [%f,%f]\n", t, T[idx], T[idx+1] );
    */
  }
}

static void test_dp_edge_weight_basic1()
{
  double T1[] = { 0.0, 1.0/6.0, 0.5, 4.0/6.0, 5.0/6.0, 1.0 };
  double Q1[] = { 1.0, -1.0, 1.0, -1.0, 1.0 };
  double T2[] = { 0.0, 2.0/6.0, 5.0/6.0, 1.0 };
  double Q2[] = { -1.0, 1.0, -1.0 };
  double tv1[] = { 0.0, 0.25, 0.5, 0.75, 1.0 };
  double tv2[] = { 0.0, 0.25, 0.5, 0.75, 1.0 };
  int nsamps1 = sizeof(T1) / sizeof(double);
  int nsamps2 = sizeof(T2) / sizeof(double);
  int ntv1 = sizeof(tv1) / sizeof(double);
  int ntv2 = sizeof(tv2) / sizeof(double);
  int dim = 1;

  double W[ntv2*ntv1*ntv2*ntv1]; /* edge weight matrix */
  double E[ntv2*ntv1];  /* DP cost matrix */
  int P[ntv2*ntv1];     /* DP predecessor matrix */
  int idxv1b[ntv1];
  int idxv2b[ntv2];
  double G[ntv1];   /* gamma function values */
  double T[ntv1];   /* gamma parameter values */
  int npts;

  double a[] = { 0.25, 0.75, 0.0, 0.0, 0.25, 1.0/6.0 };
  double b[] = { 0.75, 1.0, 0.75, 0.25, 0.5, 0.5 };
  double c[] = { 0.25, 0.25, 0.0, 0.0, 0.0, 1.0/6.0 };
  double d[] = { 0.5, 1.0, 0.25, 0.25, 0.25, 5.0/6.0 };
  double expected[] = { 0.51430, 0.90377, 0.90377, 0.66667, 0.0, 1.4714 };
  int ncases = sizeof(expected) / sizeof(double);
  int idxv1a, idxv2a;

  double actual;
  int i;

  for ( i=0; i<ncases; ++i )
  {
    /* Can't use dp_all_indexes() for this because a and c are 
     * not sorted in ascending order */
    idxv1a = dp_lookup(T1,nsamps1,a[i]);
    idxv2a = dp_lookup(T2,nsamps2,c[i]);

    actual = dp_edge_weight ( 
      Q1, T1, nsamps1, 
      Q2, T2, nsamps2, 
      dim, a[i], b[i], c[i], d[i], idxv1a, idxv2a );
    CU_ASSERT_DOUBLE_EQUAL( actual, expected[i], 0.0001 );
    if ( fabs(actual - expected[i]) > 0.0001 )
    {
      printf( "\nCase %d: expected %0.5f but got %0.5f\n", i, expected[i], actual );
    }
  }

  dp_all_indexes(T1,nsamps1,tv1,ntv1,idxv1b);
  dp_all_indexes(T2,nsamps2,tv2,ntv2,idxv2b);

  dp_all_edge_weights(Q1,T1,6,Q2,T2,4,1,tv1,idxv1b,ntv1,tv2,idxv2b,ntv2,W);
  dp_costs(Q1,T1,6,Q2,T2,4,1,tv1,idxv1b,ntv1,tv2,idxv2b,ntv2,E,P);
  npts = dp_build_gamma( P, tv1, ntv1, tv2, ntv2, G, T );

  printf( "gamma:\n" );
  for ( i=0; i<npts; ++i )
  {
    printf( "(%0.2f,%0.2f)\n", T[i], G[i] );
  }
}

static void test_dp_edge_weight_basic2()
{
  /* Q1 and Q2 are column-major arrays! */
  double Q1[2*3] = 
  {
    1,1,
    1,1,
    1,1
  };
  double T1[]={0,1.0/3.0,2.0/3.0,1};
  double Q2[2*4] = 
  {
    1,-1,
    -1,1,
    1,-1,
    -1,1
  };
  double T2[]={0,.25,.5,.75,1};
  double tv1[]={0,1.0/3.0,2.0/3.0,1};
  double tv2[]={0,0.5,1};
  int idxv1[4];
  int idxv2[3];
  double W[144]; /* edge weight matrix */
  double E[12];  /* DP cost matrix */
  int P[12];     /* DP predecessor matrix */
  double G[4];   /* gamma function values */
  double T[4];   /* gamma parameter values */
  double retval;

  double Eexp[] = { 
    0.0, 1e9, 1e9, 1e9, 
    1e9, 5.0/3.0, 7.0/3.0, 3, 
    1e9, 8.0/3.0, 10.0/3.0, 4.0
  };

  int nsamps1 = sizeof(T1) / sizeof(double);
  int nsamps2 = sizeof(T2) / sizeof(double);
  int ntv1 = sizeof(tv1) / sizeof(double);
  int ntv2 = sizeof(tv2) / sizeof(double);
  int dim = 2;
  int Gnsamps;
  int i;

  dp_all_indexes(T1,nsamps1,tv1,ntv1,idxv1);
  dp_all_indexes(T2,nsamps2,tv2,ntv2,idxv2);

  retval = dp_costs( Q1, T1, nsamps1, Q2, T2, nsamps2,
    dim, tv1, idxv1, ntv1, tv2, idxv2, ntv2, E, P );
  CU_ASSERT_DOUBLE_EQUAL( retval, E[ntv1*ntv2-1], 1e-3 );

  for ( i=0; i<ntv1*ntv2; ++i )
  {
    CU_ASSERT_DOUBLE_EQUAL( E[i], Eexp[i], 1e-3 );
  }

  dp_all_edge_weights (
    Q1, T1, nsamps1, Q2, T2, nsamps2, 
    dim, tv1, idxv1, ntv1, tv2, idxv2, ntv2, W );

  Gnsamps = dp_build_gamma(P,tv1,ntv1,tv2,ntv2,G,T);
  printf( "Gamma: " );
  for ( i=0; i<Gnsamps; ++i )
  {
    printf( "(%0.3f,%0.3f)", T[i], G[i] );
  }
  printf( "\n" );
}

static void test_dp_edge_weight_basic3()
{
  /* Q1 and Q2 are column-major arrays! */
  double Q1[2*7] = 
  {
    -0.12093, -3.29569,
    -3.99025,  1.89471,
     0.41784, -5.22274,
    -3.41510, -3.55044,
    -4.53925, -2.88102,
    -3.91785, -3.61769,
    -2.95956, -4.41621
  };
  double T1[8] =
  {
    0.00000, 0.14286, 0.28571, 0.42857, 
    0.57143, 0.71429, 0.85714, 1.00000
  };
  double Q2[2*9] = 
  {
    -3.94476, -4.47761,
    -5.49738, -1.25924,
    -6.12657,  0.84389,
    -6.08871, -0.48200,
    -5.89130, -1.34482,
    -6.17836,  0.40485,
    -6.01889,  1.04033,
    -6.09084,  0.86223,
    -5.45081, -2.80984
  };
  double T2[10] = 
  {
    0.00000, 0.11111, 0.22222, 0.33333, 0.44444, 
    0.55556, 0.66667, 0.77778, 0.88889, 1.00000
  };
  int nsamps1 = sizeof(T1) / sizeof(double);
  int nsamps2 = sizeof(T2) / sizeof(double);
  int idxv1[nsamps1];
  int idxv2[nsamps2];
  int dim = 2;
  double E[nsamps1*nsamps2];
  int P[nsamps1*nsamps2];
  double retval;

  dp_all_indexes(T1,nsamps1,T1,nsamps1,idxv1);
  dp_all_indexes(T2,nsamps2,T2,nsamps2,idxv2);

  retval = dp_costs( Q1, T1, nsamps1, Q2, T2, nsamps2,
    dim, T1, idxv1, nsamps1, T2, idxv2, nsamps2, E, P );
  CU_ASSERT_DOUBLE_EQUAL( retval, 19.629, 1e-3 );
}

static void test_dp_edge_weight_timing1()
{
  double *Q1=0, *T1=0;
  double *Q2=0, *T2=0;
  int *idxv1=0;
  int *idxv2=0;
  double *E=0;
  int *P=0;
  int nsamps1 = 200;
  int nsamps2 = 250;
  
  Q1 = (double*)malloc( (nsamps1-1)*sizeof(double) );
  Q2 = (double*)malloc( (nsamps2-1)*sizeof(double) );
  T1 = (double*)malloc( nsamps1*sizeof(double) );
  T2 = (double*)malloc( nsamps2*sizeof(double) );
  idxv1 = (int*)malloc( nsamps1*sizeof(int) );
  idxv2 = (int*)malloc( nsamps2*sizeof(int) );
  E = (double*)malloc( nsamps1*nsamps2*sizeof(double) );
  P = (int*)malloc( nsamps1*nsamps2*sizeof(int) );

  CU_ASSERT_PTR_NOT_NULL_FATAL(Q1);
  CU_ASSERT_PTR_NOT_NULL_FATAL(Q2);
  CU_ASSERT_PTR_NOT_NULL_FATAL(T1);
  CU_ASSERT_PTR_NOT_NULL_FATAL(T2);
  CU_ASSERT_PTR_NOT_NULL_FATAL(E);
  CU_ASSERT_PTR_NOT_NULL_FATAL(P);

  random_partition( T1, nsamps1 );
  random_partition( T2, nsamps2 );

  dp_all_indexes(T1,nsamps1,T1,nsamps1,idxv1);
  dp_all_indexes(T2,nsamps2,T2,nsamps2,idxv2);

  dp_costs( Q1, T1, nsamps1, Q2, T2, nsamps2, 1, 
    T1, idxv1, nsamps1, T2, idxv2, nsamps2, E, P );

  if ( Q1 ) free(Q1);
  if ( Q2 ) free(Q2);
  if ( T1 ) free(T1);
  if ( T2 ) free(T2);
  if ( idxv1) free(idxv1);
  if ( idxv2) free(idxv2);
  if ( E ) free(E);
  if ( P ) free(P);
}


static void test_dp_all_indexes_basic1()
{
  double p[]={0.0,0.2,0.33333,0.75,1.0};
  double tv[]={0.0,0.1,0.1999,0.2,0.33332,0.33333,0.74,0.999,1.0};
  int exp[]={0,0,0,1,1,2,2,3,3};
  int np=sizeof(p)/sizeof(double);
  int ntv=sizeof(tv)/sizeof(double);
  int idxv[ntv];
  int i;

  dp_all_indexes(p,np,tv,ntv,idxv);
  for ( i=0; i<ntv; ++i )
  {
    CU_ASSERT_EQUAL(idxv[i],exp[i]);
  }
}


CU_ErrorCode grid_tests_suite()
{
  CU_pSuite suite = CU_add_suite( "grid_tests", NULL, NULL );
  if ( suite == NULL )
  {
    return CU_get_error();
  }
  
  if ( CU_add_test( suite, "dp_lookup() basic test 1", 
                    (CU_TestFunc)test_dp_lookup_basic1 ) == NULL )
  {
    return CU_get_error();
  }

  if ( CU_add_test( suite, "dp_lookup() random test 1", 
                    (CU_TestFunc)test_dp_lookup_random1 ) == NULL )
  {
    return CU_get_error();
  }

  if ( CU_add_test( suite, "dp_lookup() random test 2", 
                    (CU_TestFunc)test_dp_lookup_random2 ) == NULL )
  {
    return CU_get_error();
  }

  if ( CU_add_test( suite, "dp_edge_weight() basic test 1", 
                    (CU_TestFunc)test_dp_edge_weight_basic1 ) == NULL )
  {
    return CU_get_error();
  }

  if ( CU_add_test( suite, "dp_edge_weight() basic test 2", 
                    (CU_TestFunc)test_dp_edge_weight_basic2 ) == NULL )
  {
    return CU_get_error();
  }

  if ( CU_add_test( suite, "dp_edge_weight() basic test 3", 
                    (CU_TestFunc)test_dp_edge_weight_basic3 ) == NULL )
  {
    return CU_get_error();
  }

  if ( CU_add_test( suite, "dp_edge_weight() timing test 1", 
                    (CU_TestFunc)test_dp_edge_weight_timing1 ) == NULL )
  {
    return CU_get_error();
  }

  if ( CU_add_test( suite, "dp_all_indexes() basic test 1", 
                    (CU_TestFunc)test_dp_all_indexes_basic1 ) == NULL )
  {
    return CU_get_error();
  }

  return CUE_SUCCESS;
}
