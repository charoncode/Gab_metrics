#include <stdio.h>
#include <CUnit/CUnit.h>
#include <CUnit/Basic.h>

#include "all_tests.h"


int main( int argc, char **argv )
{
  CU_ErrorCode ec;
  int i;

  /* Initialize the CUnit test registry */
  ec = CU_initialize_registry( );
  if ( ec != CUE_SUCCESS )
  {
    puts( "CU_initialize_registry() failed; exiting..." );
    return -1;
  }

  /* Call the setup routines for selected test suites */
  if ( argc == 1 ){
    grid_tests_suite();
  } else {
    for ( i=1; i<argc; ++i )
    {
      if ( strcmp( "grid", argv[i] ) == 0 )
        grid_tests_suite();
    }
  }

  /* Run the tests */
  CU_basic_set_mode(CU_BRM_VERBOSE);
  CU_basic_run_tests();
  CU_cleanup_registry();
  return CU_get_error();
}
