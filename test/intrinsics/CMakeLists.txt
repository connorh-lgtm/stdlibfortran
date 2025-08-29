set(
  fppFiles
  "test_intrinsics.fypp"
  "test_intrinsics_pathological.fypp"
)
fypp_f90("${fyppFlags}" "${fppFiles}" outFiles)

ADDTEST(intrinsics)
