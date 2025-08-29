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
