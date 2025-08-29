set(
  fppFiles
  "test_intrinsics.fypp"
  "test_intrinsics_pathological.fypp"
  "test_intrinsics_accuracy.fypp"
)
fypp_f90("${fyppFlags}" "${fppFiles}" outFiles)

# Standard test registration
ADDTEST(intrinsics)
ADDTEST(intrinsics_pathological)
ADDTEST(intrinsics_accuracy)

# Cross-optimization testing (disabled due to CMake dependency issues)
# TODO: Implement cross-optimization testing once CMake file dependencies are resolved
# The issue is that add_executable is called before fypp preprocessing generates the .F90 files
# if(BUILD_OPTIMIZATION_TESTS)
#   foreach(OPT_LEVEL IN LISTS OPTIMIZATION_LEVELS)
#     add_executable(test_intrinsics_${OPT_LEVEL} ${CMAKE_CURRENT_BINARY_DIR}/test_intrinsics.F90)
#     target_link_libraries(test_intrinsics_${OPT_LEVEL} "${PROJECT_NAME}" "test-drive::test-drive")
#     target_compile_options(test_intrinsics_${OPT_LEVEL} PRIVATE "-${OPT_LEVEL}")
#     add_test(NAME intrinsics_${OPT_LEVEL}
#              COMMAND $<TARGET_FILE:test_intrinsics_${OPT_LEVEL}>
#              WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
#   endforeach()
# endif()

# Conditional precision type testing
if(WITH_XDP)
  set_tests_properties(intrinsics intrinsics_pathological intrinsics_accuracy PROPERTIES LABELS "extended_precision")
endif()
if(WITH_QP)
  set_tests_properties(intrinsics intrinsics_pathological intrinsics_accuracy PROPERTIES LABELS "quadruple_precision")
endif()

# Set resource limits for CI environments
set_tests_properties(intrinsics intrinsics_pathological intrinsics_accuracy 
                     PROPERTIES 
                     TIMEOUT 300
                     RESOURCE_LOCK memory_intensive)
