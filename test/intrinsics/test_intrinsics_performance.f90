program test_intrinsics_performance

    use stdlib_kinds, only: dp, sp
    use stdlib_intrinsics, only: stdlib_sum, stdlib_sum_kahan, stdlib_sum_adaptive
    implicit none

    integer, parameter :: &
        array_sizes(*) = [1000, 10000, 100000, 1000000]
    integer, parameter :: repeat = 1000
    integer :: i, j, k, lun
    real(dp) :: t1, t2, tdiff
    real(dp), allocatable :: test_array(:)
    real(dp) :: result_intrinsic, result_stdlib, result_kahan, result_adaptive
    real(dp) :: time_intrinsic, time_stdlib, time_kahan, time_adaptive

    open(newunit=lun, file="intrinsics_performance_log.txt", &
         access="sequential", action="write", form="formatted", &
         position="rewind")

    write(lun, '("| Algorithm     | Array Size | Time (s) | Relative Performance |")')
    write(lun, '("|---------------|------------|----------|---------------------|")')

    do i = 1, size(array_sizes)
        allocate(test_array(array_sizes(i)))
        call random_number(test_array)
        test_array = test_array * 1000.0_dp - 500.0_dp

        call cpu_time(t1)
        do j = 1, repeat
            result_intrinsic = sum(test_array)
        end do
        call cpu_time(t2)
        time_intrinsic = t2 - t1

        call cpu_time(t1)
        do j = 1, repeat
            result_stdlib = stdlib_sum(test_array)
        end do
        call cpu_time(t2)
        time_stdlib = t2 - t1

        call cpu_time(t1)
        do j = 1, repeat
            result_kahan = stdlib_sum_kahan(test_array)
        end do
        call cpu_time(t2)
        time_kahan = t2 - t1

        call cpu_time(t1)
        do j = 1, repeat
            result_adaptive = stdlib_sum_adaptive(test_array)
        end do
        call cpu_time(t2)
        time_adaptive = t2 - t1

        write(lun, '("|", a13, 2x, "|", i10, 2x, "|", f8.5, 2x, "|", f19.2, 2x, "|")') &
            'Intrinsic', array_sizes(i), time_intrinsic, 1.0_dp
        write(lun, '("|", a13, 2x, "|", i10, 2x, "|", f8.5, 2x, "|", f19.2, 2x, "|")') &
            'stdlib_sum', array_sizes(i), time_stdlib, time_stdlib/time_intrinsic
        write(lun, '("|", a13, 2x, "|", i10, 2x, "|", f8.5, 2x, "|", f19.2, 2x, "|")') &
            'Kahan', array_sizes(i), time_kahan, time_kahan/time_intrinsic
        write(lun, '("|", a13, 2x, "|", i10, 2x, "|", f8.5, 2x, "|", f19.2, 2x, "|")') &
            'Adaptive', array_sizes(i), time_adaptive, time_adaptive/time_intrinsic

        deallocate(test_array)
    end do

    close(lun)
    print *, "Performance benchmark completed. Results written to intrinsics_performance_log.txt"

end program test_intrinsics_performance
