set(
  fppFiles
  "test_intrinsics.fypp"
  "test_intrinsics_pathological.fypp"
  "test_intrinsics_accuracy.fypp"
)
fypp_f90("${fyppFlags}" "${fppFiles}" outFiles)

ADDTEST(intrinsics)
ADDTEST(intrinsics_pathological)
ADDTEST(intrinsics_accuracy)

add_executable(test_intrinsics_performance test_intrinsics_performance.f90)
target_link_libraries(test_intrinsics_performance PRIVATE fortran_stdlib)
add_test(NAME intrinsics_performance COMMAND test_intrinsics_performance)
