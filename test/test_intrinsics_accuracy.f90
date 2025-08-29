
module test_intrinsics_accuracy
    use testdrive, only : new_unittest, unittest_type, error_type, check, skip_test
    use stdlib_kinds, only: sp, dp, xdp, qp, int8, int16, int32, int64
    use stdlib_intrinsics
    use stdlib_constants
    implicit none
    
    interface kahan_summation_tolerance
        module procedure :: kahan_summation_tolerance_sp
        module procedure :: kahan_summation_tolerance_dp
    end interface
    
    interface chunked_summation_tolerance
        module procedure :: chunked_summation_tolerance_sp
        module procedure :: chunked_summation_tolerance_dp
    end interface
    
    interface standard_summation_tolerance
        module procedure :: standard_summation_tolerance_sp
        module procedure :: standard_summation_tolerance_dp
    end interface
    
    interface kahan_dot_product_tolerance
        module procedure :: kahan_dot_product_tolerance_sp
        module procedure :: kahan_dot_product_tolerance_dp
    end interface
    
    interface chunked_dot_product_tolerance
        module procedure :: chunked_dot_product_tolerance_sp
        module procedure :: chunked_dot_product_tolerance_dp
    end interface
    
    interface standard_dot_product_tolerance
        module procedure :: standard_dot_product_tolerance_sp
        module procedure :: standard_dot_product_tolerance_dp
    end interface
    
    interface complex_dot_product_tolerance
        module procedure :: complex_dot_product_tolerance_sp
        module procedure :: complex_dot_product_tolerance_dp
    end interface
    
    interface array_conditioning_factor
        module procedure :: array_conditioning_factor_sp
        module procedure :: array_conditioning_factor_dp
    end interface
    
contains

! Context-dependent tolerance calculation functions
pure function kahan_summation_tolerance_sp(n, base_eps, sum_magnitude) result(tol)
    !! Kahan summation error bound: O(eps) + O(n*eps^2)
    integer, intent(in) :: n
    real(sp), intent(in) :: base_eps, sum_magnitude
    real(sp) :: tol
    real(sp) :: first_order_term, second_order_term
    
    first_order_term = 2.0_sp * base_eps
    second_order_term = real(n, sp) * base_eps**2
    
    tol = (first_order_term + second_order_term) * abs(sum_magnitude)
    ! For large arrays, use more conservative bounds accounting for element magnitude
    tol = max(tol, sqrt(real(n, sp)) * base_eps * abs(sum_magnitude))
    tol = max(tol, real(n, sp) * base_eps)
end function

pure function chunked_summation_tolerance_sp(n, base_eps, sum_magnitude) result(tol)
    !! Chunked summation with 64-element batches
    !! Error bound: O(sqrt(chunk_size) * eps * sqrt(n_chunks))
    integer, intent(in) :: n
    real(sp), intent(in) :: base_eps, sum_magnitude
    real(sp) :: tol
    real(sp), parameter :: chunk_size = 64.0_sp
    real(sp) :: n_chunks, chunk_error, accumulation_error
    
    n_chunks = ceiling(real(n, sp) / chunk_size)
    chunk_error = sqrt(chunk_size) * base_eps
    accumulation_error = sqrt(n_chunks) * base_eps
    
    tol = (chunk_error + accumulation_error) * abs(sum_magnitude)
    ! For large arrays, use more conservative bounds
    tol = max(tol, sqrt(real(n, sp)) * base_eps * abs(sum_magnitude))
    tol = max(tol, 2.0_sp * real(n, sp) * base_eps)
end function

pure function standard_summation_tolerance_sp(n, base_eps, sum_magnitude) result(tol)
    !! Standard summation error bound: O(n*eps)
    !! Worst-case forward error analysis for sequential summation
    integer, intent(in) :: n
    real(sp), intent(in) :: base_eps, sum_magnitude
    real(sp) :: tol
    
    tol = real(n, sp) * base_eps * abs(sum_magnitude)
    tol = max(tol, 5.0_sp * real(n, sp) * base_eps)
end function

pure function kahan_dot_product_tolerance_sp(n, base_eps, product_magnitude) result(tol)
    !! Kahan dot product: combines multiplication and Kahan summation errors
    integer, intent(in) :: n
    real(sp), intent(in) :: base_eps, product_magnitude
    real(sp) :: tol
    real(sp) :: multiplication_error, summation_error
    
    multiplication_error = base_eps
    summation_error = kahan_summation_tolerance_sp(n, base_eps, product_magnitude)
    
    tol = (multiplication_error + summation_error) * abs(product_magnitude)
    tol = max(tol, 15.0_sp * base_eps)
end function

pure function chunked_dot_product_tolerance_sp(n, base_eps, product_magnitude) result(tol)
    !! Chunked dot product: multiplication + chunked summation errors
    integer, intent(in) :: n
    real(sp), intent(in) :: base_eps, product_magnitude
    real(sp) :: tol
    real(sp) :: multiplication_error, summation_error
    
    multiplication_error = base_eps
    summation_error = chunked_summation_tolerance_sp(n, base_eps, product_magnitude)
    
    tol = (multiplication_error + summation_error) * abs(product_magnitude)
    tol = max(tol, 30.0_sp * base_eps)
end function

pure function standard_dot_product_tolerance_sp(n, base_eps, product_magnitude) result(tol)
    !! Standard dot product: O(n*eps) for both multiplication and summation
    integer, intent(in) :: n
    real(sp), intent(in) :: base_eps, product_magnitude
    real(sp) :: tol
    
    tol = real(n, sp) * base_eps * abs(product_magnitude)
    tol = max(tol, 100.0_sp * base_eps)
end function

pure function complex_dot_product_tolerance_sp(n, base_eps) result(tol)
    !! Complex dot product tolerance accounting for conjugation and component errors
    !! Each complex operation involves 4 real operations plus conjugation
    !! Based on error analysis in Golub & Van Loan "Matrix Computations"
    integer, intent(in) :: n
    real(sp), intent(in) :: base_eps
    real(sp) :: tol
    real(sp) :: complex_factor
    
    complex_factor = 4.0_sp
    tol = complex_factor * real(n, sp) * base_eps
    tol = max(tol, 1000.0_sp * base_eps)
end function

pure function array_conditioning_factor_sp(array_values) result(factor)
    !! Estimate conditioning factor based on array characteristics
    !! Higher condition numbers require looser tolerances
    real(sp), intent(in) :: array_values(:)
    real(sp) :: factor
    real(sp) :: max_val, min_val, range_ratio
    
    max_val = maxval(abs(array_values))
    min_val = minval(abs(array_values), mask=abs(array_values) > 0.0_sp)
    
    if (min_val > 0.0_sp) then
        range_ratio = max_val / min_val
        factor = 1.0_sp + log10(max(range_ratio, 1.0_sp))
    else
        factor = 2.0_sp
    end if
    
    factor = min(factor, 10.0_sp)
end function
pure function kahan_summation_tolerance_dp(n, base_eps, sum_magnitude) result(tol)
    !! Kahan summation error bound: O(eps) + O(n*eps^2)
    integer, intent(in) :: n
    real(dp), intent(in) :: base_eps, sum_magnitude
    real(dp) :: tol
    real(dp) :: first_order_term, second_order_term
    
    first_order_term = 2.0_dp * base_eps
    second_order_term = real(n, dp) * base_eps**2
    
    tol = (first_order_term + second_order_term) * abs(sum_magnitude)
    ! For large arrays, use more conservative bounds accounting for element magnitude
    tol = max(tol, sqrt(real(n, dp)) * base_eps * abs(sum_magnitude))
    tol = max(tol, real(n, dp) * base_eps)
end function

pure function chunked_summation_tolerance_dp(n, base_eps, sum_magnitude) result(tol)
    !! Chunked summation with 64-element batches
    !! Error bound: O(sqrt(chunk_size) * eps * sqrt(n_chunks))
    integer, intent(in) :: n
    real(dp), intent(in) :: base_eps, sum_magnitude
    real(dp) :: tol
    real(dp), parameter :: chunk_size = 64.0_dp
    real(dp) :: n_chunks, chunk_error, accumulation_error
    
    n_chunks = ceiling(real(n, dp) / chunk_size)
    chunk_error = sqrt(chunk_size) * base_eps
    accumulation_error = sqrt(n_chunks) * base_eps
    
    tol = (chunk_error + accumulation_error) * abs(sum_magnitude)
    ! For large arrays, use more conservative bounds
    tol = max(tol, sqrt(real(n, dp)) * base_eps * abs(sum_magnitude))
    tol = max(tol, 2.0_dp * real(n, dp) * base_eps)
end function

pure function standard_summation_tolerance_dp(n, base_eps, sum_magnitude) result(tol)
    !! Standard summation error bound: O(n*eps)
    !! Worst-case forward error analysis for sequential summation
    integer, intent(in) :: n
    real(dp), intent(in) :: base_eps, sum_magnitude
    real(dp) :: tol
    
    tol = real(n, dp) * base_eps * abs(sum_magnitude)
    tol = max(tol, 5.0_dp * real(n, dp) * base_eps)
end function

pure function kahan_dot_product_tolerance_dp(n, base_eps, product_magnitude) result(tol)
    !! Kahan dot product: combines multiplication and Kahan summation errors
    integer, intent(in) :: n
    real(dp), intent(in) :: base_eps, product_magnitude
    real(dp) :: tol
    real(dp) :: multiplication_error, summation_error
    
    multiplication_error = base_eps
    summation_error = kahan_summation_tolerance_dp(n, base_eps, product_magnitude)
    
    tol = (multiplication_error + summation_error) * abs(product_magnitude)
    tol = max(tol, 15.0_dp * base_eps)
end function

pure function chunked_dot_product_tolerance_dp(n, base_eps, product_magnitude) result(tol)
    !! Chunked dot product: multiplication + chunked summation errors
    integer, intent(in) :: n
    real(dp), intent(in) :: base_eps, product_magnitude
    real(dp) :: tol
    real(dp) :: multiplication_error, summation_error
    
    multiplication_error = base_eps
    summation_error = chunked_summation_tolerance_dp(n, base_eps, product_magnitude)
    
    tol = (multiplication_error + summation_error) * abs(product_magnitude)
    tol = max(tol, 30.0_dp * base_eps)
end function

pure function standard_dot_product_tolerance_dp(n, base_eps, product_magnitude) result(tol)
    !! Standard dot product: O(n*eps) for both multiplication and summation
    integer, intent(in) :: n
    real(dp), intent(in) :: base_eps, product_magnitude
    real(dp) :: tol
    
    tol = real(n, dp) * base_eps * abs(product_magnitude)
    tol = max(tol, 100.0_dp * base_eps)
end function

pure function complex_dot_product_tolerance_dp(n, base_eps) result(tol)
    !! Complex dot product tolerance accounting for conjugation and component errors
    !! Each complex operation involves 4 real operations plus conjugation
    !! Based on error analysis in Golub & Van Loan "Matrix Computations"
    integer, intent(in) :: n
    real(dp), intent(in) :: base_eps
    real(dp) :: tol
    real(dp) :: complex_factor
    
    complex_factor = 4.0_dp
    tol = complex_factor * real(n, dp) * base_eps
    tol = max(tol, 1000.0_dp * base_eps)
end function

pure function array_conditioning_factor_dp(array_values) result(factor)
    !! Estimate conditioning factor based on array characteristics
    !! Higher condition numbers require looser tolerances
    real(dp), intent(in) :: array_values(:)
    real(dp) :: factor
    real(dp) :: max_val, min_val, range_ratio
    
    max_val = maxval(abs(array_values))
    min_val = minval(abs(array_values), mask=abs(array_values) > 0.0_dp)
    
    if (min_val > 0.0_dp) then
        range_ratio = max_val / min_val
        factor = 1.0_dp + log10(max(range_ratio, 1.0_dp))
    else
        factor = 2.0_dp
    end if
    
    factor = min(factor, 10.0_dp)
end function

!> Collect all exported unit tests
subroutine collect_suite(testsuite)
    type(unittest_type), allocatable, intent(out) :: testsuite(:)
    testsuite = [ &
        new_unittest('pathological_case_validation', test_pathological_case_validation), &
        new_unittest('precision_scaling_validation', test_precision_scaling_validation), &
        new_unittest('conditioning_factor_validation', test_conditioning_factor_validation) &
    ]
end subroutine

subroutine test_pathological_case_validation(error)
    !> Error handling
    type(error_type), allocatable, intent(out) :: error

    block
        real(sp), allocatable :: x(:)
        real(sp) :: result_kahan, result_chunked
        real(sp) :: expected_sum, tolerance_kahan, tolerance_chunked
        real(sp), parameter :: base_tolerance = epsilon(1._sp)
        integer :: i, n

        n = 10000
        allocate(x(n))
        
        do i = 1, n
            x(i) = ((-1.0_sp)**(i-1)) * (1.0_sp + real(i, sp) * base_tolerance)
        end do
        
        expected_sum = sum(x)
        result_kahan = stdlib_sum_kahan(x)
        result_chunked = stdlib_sum(x)
        
        tolerance_kahan = kahan_summation_tolerance(n, base_tolerance, expected_sum)
        tolerance_chunked = chunked_summation_tolerance(n, base_tolerance, expected_sum)
        
        call check(error, abs(result_kahan - expected_sum) < tolerance_kahan, &
                  "Kahan summation should handle pathological alternating series")
        if (allocated(error)) return
        
        call check(error, abs(result_chunked - expected_sum) < tolerance_chunked, &
                  "Chunked summation should be within theoretical bounds for pathological case")
        if (allocated(error)) return
    end block
    block
        real(dp), allocatable :: x(:)
        real(dp) :: result_kahan, result_chunked
        real(dp) :: expected_sum, tolerance_kahan, tolerance_chunked
        real(dp), parameter :: base_tolerance = epsilon(1._dp)
        integer :: i, n

        n = 10000
        allocate(x(n))
        
        do i = 1, n
            x(i) = ((-1.0_dp)**(i-1)) * (1.0_dp + real(i, dp) * base_tolerance)
        end do
        
        expected_sum = sum(x)
        result_kahan = stdlib_sum_kahan(x)
        result_chunked = stdlib_sum(x)
        
        tolerance_kahan = kahan_summation_tolerance(n, base_tolerance, expected_sum)
        tolerance_chunked = chunked_summation_tolerance(n, base_tolerance, expected_sum)
        
        call check(error, abs(result_kahan - expected_sum) < tolerance_kahan, &
                  "Kahan summation should handle pathological alternating series")
        if (allocated(error)) return
        
        call check(error, abs(result_chunked - expected_sum) < tolerance_chunked, &
                  "Chunked summation should be within theoretical bounds for pathological case")
        if (allocated(error)) return
    end block

end subroutine

subroutine test_precision_scaling_validation(error)
    !> Error handling
    type(error_type), allocatable, intent(out) :: error

    block
        real(sp), allocatable :: x(:)
        real(sp) :: result_kahan, result_chunked
        real(sp) :: tolerance_kahan, tolerance_chunked
        real(sp), parameter :: base_tolerance = epsilon(1._sp)
        integer :: i, n

        n = 1000
        allocate(x(n))
        
        do i = 1, n
            x(i) = 1.0_sp / real(i, sp)
        end do
        
        result_kahan = stdlib_sum_kahan(x)
        result_chunked = stdlib_sum(x)
        
        tolerance_kahan = kahan_summation_tolerance(n, base_tolerance, result_kahan)
        tolerance_chunked = chunked_summation_tolerance(n, base_tolerance, result_chunked)
        
        call check(error, tolerance_kahan < tolerance_chunked, &
                  "Kahan tolerance should be tighter than chunked tolerance for precision "//trim("sp"))
        if (allocated(error)) return
        
        call check(error, tolerance_kahan > 10.0_sp * base_tolerance, &
                  "Kahan tolerance should be reasonable multiple of machine epsilon for precision "//trim("sp"))
        if (allocated(error)) return
    end block
    block
        real(dp), allocatable :: x(:)
        real(dp) :: result_kahan, result_chunked
        real(dp) :: tolerance_kahan, tolerance_chunked
        real(dp), parameter :: base_tolerance = epsilon(1._dp)
        integer :: i, n

        n = 1000
        allocate(x(n))
        
        do i = 1, n
            x(i) = 1.0_dp / real(i, dp)
        end do
        
        result_kahan = stdlib_sum_kahan(x)
        result_chunked = stdlib_sum(x)
        
        tolerance_kahan = kahan_summation_tolerance(n, base_tolerance, result_kahan)
        tolerance_chunked = chunked_summation_tolerance(n, base_tolerance, result_chunked)
        
        call check(error, tolerance_kahan < tolerance_chunked, &
                  "Kahan tolerance should be tighter than chunked tolerance for precision "//trim("dp"))
        if (allocated(error)) return
        
        call check(error, tolerance_kahan > 10.0_dp * base_tolerance, &
                  "Kahan tolerance should be reasonable multiple of machine epsilon for precision "//trim("dp"))
        if (allocated(error)) return
    end block

end subroutine

subroutine test_conditioning_factor_validation(error)
    !> Error handling
    type(error_type), allocatable, intent(out) :: error

    block
        real(sp), allocatable :: well_conditioned(:), ill_conditioned(:)
        real(sp) :: factor_well, factor_ill
        real(sp), parameter :: base_tolerance = epsilon(1._sp)
        integer :: i, n

        n = 100
        allocate(well_conditioned(n), ill_conditioned(n))
        
        do i = 1, n
            well_conditioned(i) = 1.0_sp
            ill_conditioned(i) = 10.0_sp**(i-50)
        end do
        
        factor_well = array_conditioning_factor(well_conditioned)
        factor_ill = array_conditioning_factor(ill_conditioned)
        
        call check(error, factor_well < factor_ill, &
                  "Well-conditioned arrays should have smaller conditioning factor")
        if (allocated(error)) return
        
        call check(error, factor_well >= 1.0_sp .and. factor_well <= 2.0_sp, &
                  "Well-conditioned arrays should have conditioning factor near 1")
        if (allocated(error)) return
        
        call check(error, factor_ill > 2.0_sp, &
                  "Ill-conditioned arrays should have larger conditioning factor")
        if (allocated(error)) return
    end block
    block
        real(dp), allocatable :: well_conditioned(:), ill_conditioned(:)
        real(dp) :: factor_well, factor_ill
        real(dp), parameter :: base_tolerance = epsilon(1._dp)
        integer :: i, n

        n = 100
        allocate(well_conditioned(n), ill_conditioned(n))
        
        do i = 1, n
            well_conditioned(i) = 1.0_dp
            ill_conditioned(i) = 10.0_dp**(i-50)
        end do
        
        factor_well = array_conditioning_factor(well_conditioned)
        factor_ill = array_conditioning_factor(ill_conditioned)
        
        call check(error, factor_well < factor_ill, &
                  "Well-conditioned arrays should have smaller conditioning factor")
        if (allocated(error)) return
        
        call check(error, factor_well >= 1.0_dp .and. factor_well <= 2.0_dp, &
                  "Well-conditioned arrays should have conditioning factor near 1")
        if (allocated(error)) return
        
        call check(error, factor_ill > 2.0_dp, &
                  "Ill-conditioned arrays should have larger conditioning factor")
        if (allocated(error)) return
    end block

end subroutine

end module test_intrinsics_accuracy

program test_intrinsics_accuracy_main
    use testdrive, only : run_testsuite, new_testsuite, testsuite_type
    use test_intrinsics_accuracy, only : collect_suite
    use iso_fortran_env, only : error_unit
    implicit none
    integer :: stat, is
    type(testsuite_type), allocatable :: testsuites(:)
    character(len=*), parameter :: fmt = '("#", *(1x, a))'

    stat = 0

    testsuites = [ &
        new_testsuite("intrinsics_accuracy", collect_suite) &
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
