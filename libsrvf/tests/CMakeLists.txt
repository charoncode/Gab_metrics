include_directories(${LIBSRVF_SOURCE_DIR}/include)
link_directories(${LIBSRVF_BINARY_DIR})

add_executable(run_tests fileio_tests.cc function_tests.cc interp_tests.cc kdtree_tests.cc main.cc matrix_tests.cc minmaxheap_tests.cc numeric_tests.cc opencurves_tests.cc partialmatch_tests.cc plf_tests.cc pointset_tests.cc point_tests.cc qmap_tests.cc render_tests.cc reparam_tests.cc rotate_tests.cc srvf_tests.cc util_tests.cc)

target_link_libraries(run_tests srvf boost_unit_test_framework boost_test_exec_monitor)

if(GSL)
  set_target_properties(run_tests PROPERTIES COMPILE_FLAGS "-DUSE_GSL")
endif()

file(COPY data DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
