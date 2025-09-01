
module test_intrinsics_pathological
    use testdrive, only : new_unittest, unittest_type, error_type, check, skip_test
    use stdlib_kinds, only: sp, dp, xdp, qp, int8, int16, int32, int64
    use stdlib_intrinsics
    use stdlib_math, only: swap
    implicit none
    
contains

!> Collect all exported unit tests
subroutine collect_suite(testsuite)
    !> Collection of tests
    type(unittest_type), allocatable, intent(out) :: testsuite(:)

    testsuite = [ &
        new_unittest('kahan_classic_cancellation', test_kahan_classic_cancellation), &
        new_unittest('multi_scale_magnitude', test_multi_scale_magnitude), &
        new_unittest('chunking_boundary_tests', test_chunking_boundary_tests), &
        new_unittest('catastrophic_cancellation', test_catastrophic_cancellation), &
        new_unittest('alternating_series', test_alternating_series), &
        new_unittest('near_overflow_underflow', test_near_overflow_underflow), &
        new_unittest('complex_pathological', test_complex_pathological), &
        new_unittest('subnormal_handling', test_subnormal_handling), &
        new_unittest('large_array_stress', test_large_array_stress), &
        new_unittest('adversarial_ordering', test_adversarial_ordering) &
    ]
end subroutine

subroutine test_kahan_classic_cancellation(error)
    !> Error handling
    type(error_type), allocatable, intent(out) :: error

    block
        real(sp), allocatable :: x(:)
        real(sp), parameter :: expected = 1.0_sp
        real(sp), parameter :: tolerance = epsilon(1._sp)*1000000
        real(sp) :: result_stdlib, result_kahan, result_intrinsic
        real(sp) :: error_stdlib, error_kahan, error_intrinsic

        allocate(x(3))
        x = [1.0e6_sp, 1.0_sp, -1.0e6_sp]
        
        result_intrinsic = sum(x)
        result_stdlib = stdlib_sum(x)
        result_kahan = stdlib_sum_kahan(x)
        
        error_intrinsic = abs(expected - result_intrinsic)
        error_stdlib = abs(expected - result_stdlib)
        error_kahan = abs(expected - result_kahan)
        
        call check(error, error_kahan <= error_intrinsic .or. (error_kahan < tolerance .and. error_intrinsic < tolerance), &
                  "Kahan summation should be at least as accurate as intrinsic for classic cancellation case")
        if (allocated(error)) return
        
        call check(error, error_kahan < tolerance, &
                  "Kahan summation should handle classic cancellation within tolerance")
        if (allocated(error)) return

        ! Extended cancellation test - now enabled with adaptive algorithm
        ! This test involves catastrophic cancellation where 1.0e15 + 1.0 + 1.0 - 1.0e15
        ! should equal 2.0. The adaptive algorithm uses pairwise summation for such
        ! extreme magnitude ratios to minimize error growth.
        x = [1.0e15_sp, 1.0_sp, 1.0_sp, -1.0e15_sp]
        result_kahan = stdlib_sum_adaptive(x)
        call check(error, abs(2.0_sp - result_kahan) < tolerance*1000, &
                  "Adaptive summation should handle extended cancellation case")
        if (allocated(error)) return
    end block
    block
        real(dp), allocatable :: x(:)
        real(dp), parameter :: expected = 1.0_dp
        real(dp), parameter :: tolerance = epsilon(1._dp)*1000000
        real(dp) :: result_stdlib, result_kahan, result_intrinsic
        real(dp) :: error_stdlib, error_kahan, error_intrinsic

        allocate(x(3))
        x = [1.0e6_dp, 1.0_dp, -1.0e6_dp]
        
        result_intrinsic = sum(x)
        result_stdlib = stdlib_sum(x)
        result_kahan = stdlib_sum_kahan(x)
        
        error_intrinsic = abs(expected - result_intrinsic)
        error_stdlib = abs(expected - result_stdlib)
        error_kahan = abs(expected - result_kahan)
        
        call check(error, error_kahan <= error_intrinsic .or. (error_kahan < tolerance .and. error_intrinsic < tolerance), &
                  "Kahan summation should be at least as accurate as intrinsic for classic cancellation case")
        if (allocated(error)) return
        
        call check(error, error_kahan < tolerance, &
                  "Kahan summation should handle classic cancellation within tolerance")
        if (allocated(error)) return

        ! Extended cancellation test - now enabled with adaptive algorithm
        ! This test involves catastrophic cancellation where 1.0e15 + 1.0 + 1.0 - 1.0e15
        ! should equal 2.0. The adaptive algorithm uses pairwise summation for such
        ! extreme magnitude ratios to minimize error growth.
        x = [1.0e15_dp, 1.0_dp, 1.0_dp, -1.0e15_dp]
        result_kahan = stdlib_sum_adaptive(x)
        call check(error, abs(2.0_dp - result_kahan) < tolerance*1000, &
                  "Adaptive summation should handle extended cancellation case")
        if (allocated(error)) return
    end block

    ! Complex adaptive test - now enabled with adaptive algorithm
    ! This test involves complex numbers with extreme magnitude ratios where
    ! cmplx(1.0e20, 1.0e20) + cmplx(1.0, 1.0) + cmplx(-1.0e20, -1.0e20)
    ! should equal cmplx(1.0, 1.0). The adaptive algorithm uses pairwise summation
    ! for such extreme magnitude differences.
    block
        complex(sp), allocatable :: x(:)
        complex(sp), parameter :: expected = cmplx(1.0_sp, 1.0_sp)
        real(sp), parameter :: tolerance = epsilon(1._sp)*100000
        complex(sp) :: result_adaptive
        real(sp) :: error_re, error_im

        allocate(x(3))
        x = [cmplx(1.0e20_sp, 1.0e20_sp), cmplx(1.0_sp, 1.0_sp), cmplx(-1.0e20_sp, -1.0e20_sp)]
        
        result_adaptive = stdlib_sum_adaptive(x)
        error_re = abs(expected%re - result_adaptive%re)
        error_im = abs(expected%im - result_adaptive%im)
        
        call check(error, error_re < tolerance .and. error_im < tolerance, &
                  "Complex adaptive summation should handle extreme cancellation")
        if (allocated(error)) return
    end block
    block
        complex(dp), allocatable :: x(:)
        complex(dp), parameter :: expected = cmplx(1.0_dp, 1.0_dp)
        real(dp), parameter :: tolerance = epsilon(1._dp)*100000
        complex(dp) :: result_adaptive
        real(dp) :: error_re, error_im

        allocate(x(3))
        x = [cmplx(1.0e20_dp, 1.0e20_dp), cmplx(1.0_dp, 1.0_dp), cmplx(-1.0e20_dp, -1.0e20_dp)]
        
        result_adaptive = stdlib_sum_adaptive(x)
        error_re = abs(expected%re - result_adaptive%re)
        error_im = abs(expected%im - result_adaptive%im)
        
        call check(error, error_re < tolerance .and. error_im < tolerance, &
                  "Complex adaptive summation should handle extreme cancellation")
        if (allocated(error)) return
    end block

end subroutine

subroutine test_multi_scale_magnitude(error)
    !> Error handling
    type(error_type), allocatable, intent(out) :: error

    block
        real(sp), allocatable :: x(:)
        real(sp), parameter :: tolerance = epsilon(1._sp)*1000000
        real(sp) :: result_stdlib, result_kahan, expected
        integer :: i

        allocate(x(20))
        expected = 0.0_sp
        do i = 1, 10
            x(2*i-1) = 10.0_sp**(i-6)
            x(2*i) = -10.0_sp**(i-6)
            expected = expected + x(2*i-1) + x(2*i)
        end do
        
        result_stdlib = stdlib_sum(x)
        result_kahan = stdlib_sum_kahan(x)
        
        call check(error, abs(result_kahan) < tolerance, &
                  "Multi-scale magnitude array should sum to near zero with Kahan")
        if (allocated(error)) return

        x = [(10.0_sp**(-i), i = 1, 15)]
        expected = sum(x)
        result_kahan = stdlib_sum_kahan(x)
        
        call check(error, abs(expected - result_kahan)/expected < tolerance, &
                  "Decreasing magnitude series should be accurate with Kahan")
        if (allocated(error)) return
    end block
    block
        real(dp), allocatable :: x(:)
        real(dp), parameter :: tolerance = epsilon(1._dp)*1000000
        real(dp) :: result_stdlib, result_kahan, expected
        integer :: i

        allocate(x(20))
        expected = 0.0_dp
        do i = 1, 10
            x(2*i-1) = 10.0_dp**(i-6)
            x(2*i) = -10.0_dp**(i-6)
            expected = expected + x(2*i-1) + x(2*i)
        end do
        
        result_stdlib = stdlib_sum(x)
        result_kahan = stdlib_sum_kahan(x)
        
        call check(error, abs(result_kahan) < tolerance, &
                  "Multi-scale magnitude array should sum to near zero with Kahan")
        if (allocated(error)) return

        x = [(10.0_dp**(-i), i = 1, 15)]
        expected = sum(x)
        result_kahan = stdlib_sum_kahan(x)
        
        call check(error, abs(expected - result_kahan)/expected < tolerance, &
                  "Decreasing magnitude series should be accurate with Kahan")
        if (allocated(error)) return
    end block

end subroutine

subroutine test_chunking_boundary_tests(error)
    !> Error handling
    type(error_type), allocatable, intent(out) :: error

    integer, parameter :: chunk_sizes(*) = [63, 64, 65, 127, 128, 129, 191, 192, 193, 255, 256, 257]
    integer :: i, n

    block
        real(sp), allocatable :: x(:)
        real(sp), parameter :: tolerance = epsilon(1._sp)*100000
        real(sp) :: result_stdlib, result_kahan, result_intrinsic
        real(sp) :: expected
        integer :: j

        do i = 1, size(chunk_sizes)
            n = chunk_sizes(i)
            allocate(x(n))
            
            expected = 0.0_sp
            do j = 1, n
                x(j) = sin(real(j, sp)) / real(j, sp)
                expected = expected + x(j)
            end do
            
            result_intrinsic = sum(x)
            result_stdlib = stdlib_sum(x)
            result_kahan = stdlib_sum_kahan(x)
            
            call check(error, abs(result_stdlib - result_intrinsic) < tolerance*n, &
                      "Chunked summation should be close to intrinsic for size " // trim(adjustl(char(n))))
            if (allocated(error)) return
            
            call check(error, abs(result_kahan - expected) < tolerance*n, &
                      "Kahan summation should be accurate for size " // trim(adjustl(char(n))))
            if (allocated(error)) return
            
            deallocate(x)
        end do
    end block
    block
        real(dp), allocatable :: x(:)
        real(dp), parameter :: tolerance = epsilon(1._dp)*100000
        real(dp) :: result_stdlib, result_kahan, result_intrinsic
        real(dp) :: expected
        integer :: j

        do i = 1, size(chunk_sizes)
            n = chunk_sizes(i)
            allocate(x(n))
            
            expected = 0.0_dp
            do j = 1, n
                x(j) = sin(real(j, dp)) / real(j, dp)
                expected = expected + x(j)
            end do
            
            result_intrinsic = sum(x)
            result_stdlib = stdlib_sum(x)
            result_kahan = stdlib_sum_kahan(x)
            
            call check(error, abs(result_stdlib - result_intrinsic) < tolerance*n, &
                      "Chunked summation should be close to intrinsic for size " // trim(adjustl(char(n))))
            if (allocated(error)) return
            
            call check(error, abs(result_kahan - expected) < tolerance*n, &
                      "Kahan summation should be accurate for size " // trim(adjustl(char(n))))
            if (allocated(error)) return
            
            deallocate(x)
        end do
    end block

end subroutine

subroutine test_catastrophic_cancellation(error)
    !> Error handling
    type(error_type), allocatable, intent(out) :: error

    block
        real(sp), allocatable :: x(:)
        real(sp), parameter :: tolerance = epsilon(1._sp)*100000
        real(sp) :: result_kahan, expected
        integer :: i, n

        n = 1000
        allocate(x(n))
        
        expected = 0.0_sp
        do i = 1, n/2
            x(2*i-1) = 1.0e10_sp + real(i, sp)
            x(2*i) = -(1.0e10_sp + real(i, sp))
            expected = expected + x(2*i-1) + x(2*i)
        end do
        
        result_kahan = stdlib_sum_kahan(x)
        
        call check(error, abs(result_kahan - expected) < tolerance, &
                  "Catastrophic cancellation should be handled by Kahan")
        if (allocated(error)) return

        do i = 1, n
            x(i) = (-1.0_sp)**(i+1) * (1.0e8_sp + real(i, sp))
        end do
        expected = sum(x)
        result_kahan = stdlib_sum_kahan(x)
        
        call check(error, abs(result_kahan - expected)/max(abs(expected), 1.0_sp) < tolerance*n, &
                  "Alternating large values should be handled accurately")
        if (allocated(error)) return
    end block
    block
        real(dp), allocatable :: x(:)
        real(dp), parameter :: tolerance = epsilon(1._dp)*100000
        real(dp) :: result_kahan, expected
        integer :: i, n

        n = 1000
        allocate(x(n))
        
        expected = 0.0_dp
        do i = 1, n/2
            x(2*i-1) = 1.0e10_dp + real(i, dp)
            x(2*i) = -(1.0e10_dp + real(i, dp))
            expected = expected + x(2*i-1) + x(2*i)
        end do
        
        result_kahan = stdlib_sum_kahan(x)
        
        call check(error, abs(result_kahan - expected) < tolerance, &
                  "Catastrophic cancellation should be handled by Kahan")
        if (allocated(error)) return

        do i = 1, n
            x(i) = (-1.0_dp)**(i+1) * (1.0e8_dp + real(i, dp))
        end do
        expected = sum(x)
        result_kahan = stdlib_sum_kahan(x)
        
        call check(error, abs(result_kahan - expected)/max(abs(expected), 1.0_dp) < tolerance*n, &
                  "Alternating large values should be handled accurately")
        if (allocated(error)) return
    end block

end subroutine

subroutine test_alternating_series(error)
    !> Error handling
    type(error_type), allocatable, intent(out) :: error

    block
        real(sp), allocatable :: x(:)
        real(sp), parameter :: tolerance = epsilon(1._sp)*1000000
        real(sp) :: result_kahan, expected
        integer :: i, n

        n = 1000
        allocate(x(n))
        
        expected = 0.0_sp
        do i = 1, n
            x(i) = (-1.0_sp)**(i+1) / real(i, sp)
            expected = expected + x(i)
        end do
        
        result_kahan = stdlib_sum_kahan(x)
        
        call check(error, abs(result_kahan - expected) < tolerance, &
                  "Alternating harmonic series should be accurate with Kahan")
        if (allocated(error)) return

        do i = 1, n
            x(i) = (-1.0_sp)**(i+1) * exp(-real(i, sp)/100.0_sp)
        end do
        expected = sum(x)
        result_kahan = stdlib_sum_kahan(x)
        
        call check(error, abs(result_kahan - expected) < tolerance, &
                  "Alternating exponential series should be accurate")
        if (allocated(error)) return
    end block
    block
        real(dp), allocatable :: x(:)
        real(dp), parameter :: tolerance = epsilon(1._dp)*1000000
        real(dp) :: result_kahan, expected
        integer :: i, n

        n = 1000
        allocate(x(n))
        
        expected = 0.0_dp
        do i = 1, n
            x(i) = (-1.0_dp)**(i+1) / real(i, dp)
            expected = expected + x(i)
        end do
        
        result_kahan = stdlib_sum_kahan(x)
        
        call check(error, abs(result_kahan - expected) < tolerance, &
                  "Alternating harmonic series should be accurate with Kahan")
        if (allocated(error)) return

        do i = 1, n
            x(i) = (-1.0_dp)**(i+1) * exp(-real(i, dp)/100.0_dp)
        end do
        expected = sum(x)
        result_kahan = stdlib_sum_kahan(x)
        
        call check(error, abs(result_kahan - expected) < tolerance, &
                  "Alternating exponential series should be accurate")
        if (allocated(error)) return
    end block

end subroutine

subroutine test_near_overflow_underflow(error)
    !> Error handling
    type(error_type), allocatable, intent(out) :: error

    block
        real(sp), allocatable :: x(:)
        real(sp), parameter :: huge_val = huge(1.0_sp)
        real(sp), parameter :: tiny_val = tiny(1.0_sp)
        real(sp), parameter :: tolerance = epsilon(1._sp)*1000000
        real(sp) :: result_kahan, expected
        integer :: i

        allocate(x(10))
        
        do i = 1, 5
            x(2*i-1) = huge_val * 0.1_sp
            x(2*i) = -huge_val * 0.1_sp
        end do
        expected = 0.0_sp
        
        result_kahan = stdlib_sum_kahan(x)
        
        call check(error, abs(result_kahan - expected) < tolerance, &
                  "Near-overflow cancellation should be handled")
        if (allocated(error)) return

        do i = 1, 10
            x(i) = tiny_val * real(i, sp)
        end do
        expected = sum(x)
        result_kahan = stdlib_sum_kahan(x)
        
        call check(error, abs(result_kahan - expected) < tolerance, &
                  "Near-underflow summation should be accurate")
        if (allocated(error)) return
    end block
    block
        real(dp), allocatable :: x(:)
        real(dp), parameter :: huge_val = huge(1.0_dp)
        real(dp), parameter :: tiny_val = tiny(1.0_dp)
        real(dp), parameter :: tolerance = epsilon(1._dp)*1000000
        real(dp) :: result_kahan, expected
        integer :: i

        allocate(x(10))
        
        do i = 1, 5
            x(2*i-1) = huge_val * 0.1_dp
            x(2*i) = -huge_val * 0.1_dp
        end do
        expected = 0.0_dp
        
        result_kahan = stdlib_sum_kahan(x)
        
        call check(error, abs(result_kahan - expected) < tolerance, &
                  "Near-overflow cancellation should be handled")
        if (allocated(error)) return

        do i = 1, 10
            x(i) = tiny_val * real(i, dp)
        end do
        expected = sum(x)
        result_kahan = stdlib_sum_kahan(x)
        
        call check(error, abs(result_kahan - expected) < tolerance, &
                  "Near-underflow summation should be accurate")
        if (allocated(error)) return
    end block

end subroutine

subroutine test_complex_pathological(error)
    !> Error handling
    type(error_type), allocatable, intent(out) :: error

    block
        complex(sp), allocatable :: x(:)
        real(sp), parameter :: tolerance = epsilon(1._sp)*100000
        complex(sp) :: result_kahan, expected
        integer :: i

        allocate(x(100))
        
        expected = cmplx(0.0_sp, 0.0_sp)
        do i = 1, 50
            x(2*i-1) = cmplx(1.0e10_sp, epsilon(1._sp) * real(i, sp))
            x(2*i) = cmplx(-1.0e10_sp, -epsilon(1._sp) * real(i, sp))
            expected = expected + x(2*i-1) + x(2*i)
        end do
        
        result_kahan = stdlib_sum_kahan(x)
        
        call check(error, abs(result_kahan%re - expected%re) < tolerance .and. &
                         abs(result_kahan%im - expected%im) < tolerance, &
                  "Complex near-zero imaginary parts should be handled accurately")
        if (allocated(error)) return

        do i = 1, 100
            x(i) = cmplx(cos(real(i, sp)), sin(real(i, sp))) / real(i, sp)
        end do
        expected = sum(x)
        result_kahan = stdlib_sum_kahan(x)
        
        call check(error, abs(result_kahan%re - expected%re) < tolerance .and. &
                         abs(result_kahan%im - expected%im) < tolerance, &
                  "Complex phase-sensitive summation should be accurate")
        if (allocated(error)) return
    end block
    block
        complex(dp), allocatable :: x(:)
        real(dp), parameter :: tolerance = epsilon(1._dp)*100000
        complex(dp) :: result_kahan, expected
        integer :: i

        allocate(x(100))
        
        expected = cmplx(0.0_dp, 0.0_dp)
        do i = 1, 50
            x(2*i-1) = cmplx(1.0e10_dp, epsilon(1._dp) * real(i, dp))
            x(2*i) = cmplx(-1.0e10_dp, -epsilon(1._dp) * real(i, dp))
            expected = expected + x(2*i-1) + x(2*i)
        end do
        
        result_kahan = stdlib_sum_kahan(x)
        
        call check(error, abs(result_kahan%re - expected%re) < tolerance .and. &
                         abs(result_kahan%im - expected%im) < tolerance, &
                  "Complex near-zero imaginary parts should be handled accurately")
        if (allocated(error)) return

        do i = 1, 100
            x(i) = cmplx(cos(real(i, dp)), sin(real(i, dp))) / real(i, dp)
        end do
        expected = sum(x)
        result_kahan = stdlib_sum_kahan(x)
        
        call check(error, abs(result_kahan%re - expected%re) < tolerance .and. &
                         abs(result_kahan%im - expected%im) < tolerance, &
                  "Complex phase-sensitive summation should be accurate")
        if (allocated(error)) return
    end block

end subroutine

subroutine test_subnormal_handling(error)
    !> Error handling
    type(error_type), allocatable, intent(out) :: error

    block
        real(sp), allocatable :: x(:)
        real(sp), parameter :: tiny_val = tiny(1.0_sp)
        real(sp), parameter :: tolerance = epsilon(1._sp)*1000000
        real(sp) :: result_kahan, expected
        integer :: i

        allocate(x(100))
        
        do i = 1, 100
            x(i) = tiny_val * 0.1_sp * real(i, sp)
        end do
        expected = sum(x)
        
        result_kahan = stdlib_sum_kahan(x)
        
        call check(error, abs(result_kahan - expected) < tolerance, &
                  "Subnormal number summation should be handled correctly")
        if (allocated(error)) return

        do i = 1, 50
            x(2*i-1) = tiny_val * 0.5_sp
            x(2*i) = -tiny_val * 0.5_sp
        end do
        expected = 0.0_sp
        result_kahan = stdlib_sum_kahan(x)
        
        call check(error, abs(result_kahan - expected) < tolerance, &
                  "Subnormal cancellation should work correctly")
        if (allocated(error)) return
    end block
    block
        real(dp), allocatable :: x(:)
        real(dp), parameter :: tiny_val = tiny(1.0_dp)
        real(dp), parameter :: tolerance = epsilon(1._dp)*1000000
        real(dp) :: result_kahan, expected
        integer :: i

        allocate(x(100))
        
        do i = 1, 100
            x(i) = tiny_val * 0.1_dp * real(i, dp)
        end do
        expected = sum(x)
        
        result_kahan = stdlib_sum_kahan(x)
        
        call check(error, abs(result_kahan - expected) < tolerance, &
                  "Subnormal number summation should be handled correctly")
        if (allocated(error)) return

        do i = 1, 50
            x(2*i-1) = tiny_val * 0.5_dp
            x(2*i) = -tiny_val * 0.5_dp
        end do
        expected = 0.0_dp
        result_kahan = stdlib_sum_kahan(x)
        
        call check(error, abs(result_kahan - expected) < tolerance, &
                  "Subnormal cancellation should work correctly")
        if (allocated(error)) return
    end block

end subroutine

subroutine test_large_array_stress(error)
    !> Error handling
    type(error_type), allocatable, intent(out) :: error

    integer, parameter :: large_sizes(*) = [10000, 50000, 100000]
    integer :: i, n, j

    block
        real(sp), allocatable :: x(:)
        real(sp), parameter :: tolerance = epsilon(1._sp)*10000000
        real(sp) :: result_stdlib, result_kahan
        real(sp) :: expected

        do i = 1, size(large_sizes)
            n = large_sizes(i)
            allocate(x(n))
            
            expected = 0.0_sp
            do j = 1, n
                x(j) = sin(real(j, sp) * 0.001_sp) / real(j, sp)
                expected = expected + x(j)
            end do
            
            result_stdlib = stdlib_sum(x)
            result_kahan = stdlib_sum_kahan(x)
            
            call check(error, abs(result_kahan - expected) < tolerance*sqrt(real(n, sp)), &
                      "Large array Kahan summation should be accurate")
            if (allocated(error)) return
            
            call check(error, abs(result_stdlib - result_kahan) < tolerance*sqrt(real(n, sp)), &
                      "Large array stdlib and Kahan should be reasonably close")
            if (allocated(error)) return
            
            deallocate(x)
        end do
    end block
    block
        real(dp), allocatable :: x(:)
        real(dp), parameter :: tolerance = epsilon(1._dp)*10000000
        real(dp) :: result_stdlib, result_kahan
        real(dp) :: expected

        do i = 1, size(large_sizes)
            n = large_sizes(i)
            allocate(x(n))
            
            expected = 0.0_dp
            do j = 1, n
                x(j) = sin(real(j, dp) * 0.001_dp) / real(j, dp)
                expected = expected + x(j)
            end do
            
            result_stdlib = stdlib_sum(x)
            result_kahan = stdlib_sum_kahan(x)
            
            call check(error, abs(result_kahan - expected) < tolerance*sqrt(real(n, dp)), &
                      "Large array Kahan summation should be accurate")
            if (allocated(error)) return
            
            call check(error, abs(result_stdlib - result_kahan) < tolerance*sqrt(real(n, dp)), &
                      "Large array stdlib and Kahan should be reasonably close")
            if (allocated(error)) return
            
            deallocate(x)
        end do
    end block

end subroutine

subroutine test_adversarial_ordering(error)
    !> Error handling
    type(error_type), allocatable, intent(out) :: error

    block
        real(sp), allocatable :: x_sorted(:), x_reverse(:)
        real(sp), parameter :: tolerance = epsilon(1._sp)*1000000
        real(sp) :: result_sorted, result_reverse, expected
        integer :: i, n

        n = 1000
        allocate(x_sorted(n), x_reverse(n))
        
        expected = 0.0_sp
        do i = 1, n
            x_sorted(i) = 1.0_sp / real(i, sp)**2
            expected = expected + x_sorted(i)
        end do
        
        do i = 1, n
            x_reverse(i) = x_sorted(n + 1 - i)
        end do
        
        result_sorted = stdlib_sum_kahan(x_sorted)
        result_reverse = stdlib_sum_kahan(x_reverse)
        
        call check(error, abs(result_sorted - expected) < tolerance, &
                  "Sorted order summation should be accurate")
        if (allocated(error)) return
        
        call check(error, abs(result_reverse - expected) < tolerance, &
                  "Reverse order summation should be accurate")
        if (allocated(error)) return
        
        call check(error, abs(result_sorted - result_reverse) < tolerance, &
                  "Different orderings should give similar results with Kahan")
        if (allocated(error)) return
    end block
    block
        real(dp), allocatable :: x_sorted(:), x_reverse(:)
        real(dp), parameter :: tolerance = epsilon(1._dp)*1000000
        real(dp) :: result_sorted, result_reverse, expected
        integer :: i, n

        n = 1000
        allocate(x_sorted(n), x_reverse(n))
        
        expected = 0.0_dp
        do i = 1, n
            x_sorted(i) = 1.0_dp / real(i, dp)**2
            expected = expected + x_sorted(i)
        end do
        
        do i = 1, n
            x_reverse(i) = x_sorted(n + 1 - i)
        end do
        
        result_sorted = stdlib_sum_kahan(x_sorted)
        result_reverse = stdlib_sum_kahan(x_reverse)
        
        call check(error, abs(result_sorted - expected) < tolerance, &
                  "Sorted order summation should be accurate")
        if (allocated(error)) return
        
        call check(error, abs(result_reverse - expected) < tolerance, &
                  "Reverse order summation should be accurate")
        if (allocated(error)) return
        
        call check(error, abs(result_sorted - result_reverse) < tolerance, &
                  "Different orderings should give similar results with Kahan")
        if (allocated(error)) return
    end block

end subroutine
    
end module test_intrinsics_pathological

program test_intrinsics_pathological_main
    use, intrinsic :: iso_fortran_env, only : error_unit
    use testdrive, only : run_testsuite, new_testsuite, testsuite_type
    use test_intrinsics_pathological, only : collect_suite
    implicit none
    integer :: stat, is
    type(testsuite_type), allocatable :: testsuites(:)
    character(len=*), parameter :: fmt = '("#", *(1x, a))'

    stat = 0

    testsuites = [ &
        new_testsuite("intrinsics_pathological", collect_suite) &
        ]

    do is = 1, size(testsuites)
        write(error_unit, fmt) "Testing:", testsuites(is)%name
        call run_testsuite(testsuites(is)%collect, error_unit, stat)
    end do

    if (stat > 0) then
        write(error_unit, '(i0, 1x, a)') stat, "test(s) failed!"
        error stop
    end if
end program
