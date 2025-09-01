
module test_intrinsics
    use testdrive, only : new_unittest, unittest_type, error_type, check, skip_test
    use stdlib_kinds, only: sp, dp, xdp, qp, int8, int16, int32, int64
    use stdlib_intrinsics
    use stdlib_math, only: swap
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
    
contains

! Context-dependent tolerance calculation functions
pure function kahan_summation_tolerance_sp(n, base_eps, sum_magnitude) result(tol)
    !! Kahan summation error bound: O(eps) + O(n*eps^2)
    !! Based on Higham "Accuracy and Stability of Numerical Algorithms" Theorem 4.8
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
    tol = max(tol, 500.0_sp * base_eps)
end function

pure function kahan_dot_product_tolerance_sp(n, base_eps, product_magnitude) result(tol)
    !! Kahan dot product: combines multiplication and Kahan summation errors
    !! Multiplication: O(eps), Summation: O(eps) for Kahan
    integer, intent(in) :: n
    real(sp), intent(in) :: base_eps, product_magnitude
    real(sp) :: tol
    real(sp) :: multiplication_error, summation_error
    
    multiplication_error = base_eps
    summation_error = kahan_summation_tolerance_sp(n, base_eps, product_magnitude)
    
    tol = (multiplication_error + summation_error) * abs(product_magnitude)
    tol = max(tol, 150.0_sp * base_eps)
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
pure function kahan_summation_tolerance_dp(n, base_eps, sum_magnitude) result(tol)
    !! Kahan summation error bound: O(eps) + O(n*eps^2)
    !! Based on Higham "Accuracy and Stability of Numerical Algorithms" Theorem 4.8
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
    tol = max(tol, 500.0_dp * base_eps)
end function

pure function kahan_dot_product_tolerance_dp(n, base_eps, product_magnitude) result(tol)
    !! Kahan dot product: combines multiplication and Kahan summation errors
    !! Multiplication: O(eps), Summation: O(eps) for Kahan
    integer, intent(in) :: n
    real(dp), intent(in) :: base_eps, product_magnitude
    real(dp) :: tol
    real(dp) :: multiplication_error, summation_error
    
    multiplication_error = base_eps
    summation_error = kahan_summation_tolerance_dp(n, base_eps, product_magnitude)
    
    tol = (multiplication_error + summation_error) * abs(product_magnitude)
    tol = max(tol, 150.0_dp * base_eps)
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

!> Collect all exported unit tests
subroutine collect_suite(testsuite)
    !> Collection of tests
    type(unittest_type), allocatable, intent(out) :: testsuite(:)

    testsuite = [ &
        new_unittest('sum', test_sum), &
        new_unittest('dot_product', test_dot_product), &
        new_unittest('reference_summation_validation', test_reference_summation_validation), &
        new_unittest('dot_product_reference_validation', test_dot_product_reference_validation), &
        new_unittest('kahan_superiority_demonstration', test_kahan_superiority_demonstration), &
        new_unittest('theoretical_error_bounds', test_theoretical_error_bounds), &
        new_unittest('complex_conjugation_accuracy', test_complex_conjugation_accuracy) &
    ]
end subroutine

subroutine test_sum(error)
    !> Error handling
    type(error_type), allocatable, intent(out) :: error

    !> Internal parameters and variables
    integer, parameter :: n = 1e3, ncalc = 3
    real(sp) :: u
    integer :: iter, i, j
    !====================================================================================
    block
        integer(int32), allocatable :: x(:)
        integer(int32), parameter :: total_sum = 0_int32
        integer(int32) :: xsum(ncalc), err(ncalc)
        logical, allocatable :: mask(:), nmask(:)

        allocate(x(n+1))
        do i = 1, n+1
            x(i) = i - n/2 - 1
        end do
        allocate(mask(n+1),source=.false.); mask(1:n+1:2) = .true.
        allocate(nmask(n+1)); nmask = .not.mask
        ! scramble array
        do i = 1, n+1
            call random_number(u) 
            j = 1 + floor(n*u)
            call swap( x(i), x(j) )
            call swap( mask(i), mask(j) )
            call swap( nmask(i), nmask(j) )
        end do
        
        xsum(1) = sum(x)        ! compiler intrinsic
        xsum(2) = stdlib_sum(x) ! chunked summation
        err(1:2) = abs(total_sum-xsum(1:2))
        
        call check(error, all(err(1:2)==0_int32) , "real sum is not accurate" )
        if (allocated(error)) return

        xsum(1) = sum(x,mask)+sum(x,nmask) ! compiler intrinsic
        xsum(2) = stdlib_sum(x,mask)+stdlib_sum(x,nmask) ! chunked summation
        err(1:2) = abs(total_sum-xsum(1:2))
        
        call check(error, all(err(1:2)==0_int32) , "masked real sum is not accurate" )
        if (allocated(error)) return
    end block
    block
        integer(int64), allocatable :: x(:)
        integer(int64), parameter :: total_sum = 0_int64
        integer(int64) :: xsum(ncalc), err(ncalc)
        logical, allocatable :: mask(:), nmask(:)

        allocate(x(n+1))
        do i = 1, n+1
            x(i) = i - n/2 - 1
        end do
        allocate(mask(n+1),source=.false.); mask(1:n+1:2) = .true.
        allocate(nmask(n+1)); nmask = .not.mask
        ! scramble array
        do i = 1, n+1
            call random_number(u) 
            j = 1 + floor(n*u)
            call swap( x(i), x(j) )
            call swap( mask(i), mask(j) )
            call swap( nmask(i), nmask(j) )
        end do
        
        xsum(1) = sum(x)        ! compiler intrinsic
        xsum(2) = stdlib_sum(x) ! chunked summation
        err(1:2) = abs(total_sum-xsum(1:2))
        
        call check(error, all(err(1:2)==0_int64) , "real sum is not accurate" )
        if (allocated(error)) return

        xsum(1) = sum(x,mask)+sum(x,nmask) ! compiler intrinsic
        xsum(2) = stdlib_sum(x,mask)+stdlib_sum(x,nmask) ! chunked summation
        err(1:2) = abs(total_sum-xsum(1:2))
        
        call check(error, all(err(1:2)==0_int64) , "masked real sum is not accurate" )
        if (allocated(error)) return
    end block

    block
        real(sp), allocatable :: x(:)
        real(sp), parameter :: total_sum = 4*atan(1._sp)
        real(sp) :: tolerance_kahan, tolerance_chunked, tolerance_standard
        real(sp) :: xsum(ncalc), err(ncalc)
        logical, allocatable :: mask(:), nmask(:)

        allocate(x(n))
        do i = 1, n 
            x(i) = 8*atan(1._sp)*(real(i,kind=sp)-0.5_sp)/real(n,kind=sp)**2
        end do
        allocate(mask(n),source=.false.); mask(1:n:2) = .true.
        allocate(nmask(n)); nmask = .not.mask
        ! scramble array
        do i = 1, n
            call random_number(u) 
            j = 1 + floor(n*u)
            call swap( x(i), x(j) )
            call swap( mask(i), mask(j) )
            call swap( nmask(i), nmask(j) )
        end do
        
        xsum(1) = sum(x)        ! compiler intrinsic
        xsum(2) = stdlib_sum_kahan(x) ! chunked Kahan summation
        xsum(3) = stdlib_sum(x)       ! chunked summation
        err(1:ncalc) = abs(1._sp-xsum(1:ncalc)/total_sum)
        
        tolerance_standard = standard_summation_tolerance(n, epsilon(1._sp), total_sum)
        tolerance_kahan = kahan_summation_tolerance(n, epsilon(1._sp), total_sum)
        tolerance_chunked = chunked_summation_tolerance(n, epsilon(1._sp), total_sum)
        
        call check(error, err(1) < tolerance_standard, "compiler intrinsic sum not within standard bounds")
        if (allocated(error)) return
        call check(error, err(2) < tolerance_kahan, "Kahan sum not within theoretical bounds")
        if (allocated(error)) return
        call check(error, err(3) < tolerance_chunked, "chunked sum not within theoretical bounds")
        if (allocated(error)) return

        xsum(1) = sum(x,mask)+sum(x,nmask) ! compiler intrinsic
        xsum(2) = stdlib_sum_kahan(x,mask)+stdlib_sum_kahan(x,nmask) ! chunked Kahan summation
        xsum(3) = stdlib_sum(x,mask)+stdlib_sum(x,nmask) ! chunked summation
        err(1:ncalc) = abs(1._sp-xsum(1:ncalc)/total_sum)
        
        call check(error, err(1) < tolerance_standard, "compiler intrinsic masked sum not within standard bounds")
        if (allocated(error)) return
        call check(error, err(2) < tolerance_kahan, "Kahan masked sum not within theoretical bounds")
        if (allocated(error)) return
        call check(error, err(3) < tolerance_chunked, "chunked masked sum not within theoretical bounds")
        if (allocated(error)) return
    end block
    block
        real(dp), allocatable :: x(:)
        real(dp), parameter :: total_sum = 4*atan(1._dp)
        real(dp) :: tolerance_kahan, tolerance_chunked, tolerance_standard
        real(dp) :: xsum(ncalc), err(ncalc)
        logical, allocatable :: mask(:), nmask(:)

        allocate(x(n))
        do i = 1, n 
            x(i) = 8*atan(1._dp)*(real(i,kind=dp)-0.5_dp)/real(n,kind=dp)**2
        end do
        allocate(mask(n),source=.false.); mask(1:n:2) = .true.
        allocate(nmask(n)); nmask = .not.mask
        ! scramble array
        do i = 1, n
            call random_number(u) 
            j = 1 + floor(n*u)
            call swap( x(i), x(j) )
            call swap( mask(i), mask(j) )
            call swap( nmask(i), nmask(j) )
        end do
        
        xsum(1) = sum(x)        ! compiler intrinsic
        xsum(2) = stdlib_sum_kahan(x) ! chunked Kahan summation
        xsum(3) = stdlib_sum(x)       ! chunked summation
        err(1:ncalc) = abs(1._dp-xsum(1:ncalc)/total_sum)
        
        tolerance_standard = standard_summation_tolerance(n, epsilon(1._dp), total_sum)
        tolerance_kahan = kahan_summation_tolerance(n, epsilon(1._dp), total_sum)
        tolerance_chunked = chunked_summation_tolerance(n, epsilon(1._dp), total_sum)
        
        call check(error, err(1) < tolerance_standard, "compiler intrinsic sum not within standard bounds")
        if (allocated(error)) return
        call check(error, err(2) < tolerance_kahan, "Kahan sum not within theoretical bounds")
        if (allocated(error)) return
        call check(error, err(3) < tolerance_chunked, "chunked sum not within theoretical bounds")
        if (allocated(error)) return

        xsum(1) = sum(x,mask)+sum(x,nmask) ! compiler intrinsic
        xsum(2) = stdlib_sum_kahan(x,mask)+stdlib_sum_kahan(x,nmask) ! chunked Kahan summation
        xsum(3) = stdlib_sum(x,mask)+stdlib_sum(x,nmask) ! chunked summation
        err(1:ncalc) = abs(1._dp-xsum(1:ncalc)/total_sum)
        
        call check(error, err(1) < tolerance_standard, "compiler intrinsic masked sum not within standard bounds")
        if (allocated(error)) return
        call check(error, err(2) < tolerance_kahan, "Kahan masked sum not within theoretical bounds")
        if (allocated(error)) return
        call check(error, err(3) < tolerance_chunked, "chunked masked sum not within theoretical bounds")
        if (allocated(error)) return
    end block

    block
        complex(sp), allocatable :: x(:)
        real(sp), parameter :: total_sum = 4*atan(1._sp)
        real(sp) :: tolerance_kahan, tolerance_chunked, tolerance_standard
        real(sp) :: err(ncalc)
        complex(sp) :: xsum(ncalc)
        logical, allocatable :: mask(:), nmask(:)

        allocate(x(n))
        do i = 1, n
            x(i) = (8*atan(1._sp)*(real(i,kind=sp)-0.5_sp)/n**2)*cmplx(1._sp,1._sp)
        end do
        
        allocate(mask(n),source=.false.); mask(1:n:2) = .true.
        allocate(nmask(n)); nmask = .not.mask
        ! scramble array
        do i = 1, n
            call random_number(u) 
            j = 1 + floor(n*u)
            call swap( x(i), x(j) )
            call swap( mask(i), mask(j) )
            call swap( nmask(i), nmask(j) )
        end do
        
        xsum(1) = sum(x)        ! compiler intrinsic
        xsum(2) = stdlib_sum_kahan(x) ! chunked Kahan summation
        xsum(3) = stdlib_sum(x)       ! chunked summation
        err(1:ncalc) = abs(1._sp-(xsum(1:ncalc)%re)/total_sum)
        
        tolerance_standard = standard_summation_tolerance(n, epsilon(1._sp), total_sum)
        tolerance_kahan = kahan_summation_tolerance(n, epsilon(1._sp), total_sum)
        tolerance_chunked = chunked_summation_tolerance(n, epsilon(1._sp), total_sum)
        
        call check(error, err(1) < tolerance_standard, "compiler intrinsic complex sum not within standard bounds")
        if (allocated(error)) return
        call check(error, err(2) < tolerance_kahan, "Kahan complex sum not within theoretical bounds")
        if (allocated(error)) return
        call check(error, err(3) < tolerance_chunked, "chunked complex sum not within theoretical bounds")
        if (allocated(error)) return

        xsum(1) = sum(x,mask)+sum(x,nmask) ! compiler intrinsic
        xsum(2) = stdlib_sum_kahan(x,mask)+stdlib_sum_kahan(x,nmask) ! chunked Kahan summation
        xsum(3) = stdlib_sum(x,mask)+stdlib_sum(x,nmask) ! chunked summation
        err(1:ncalc) = abs(1._sp-(xsum(1:ncalc)%re)/total_sum)
        
        call check(error, err(1) < tolerance_standard, "compiler intrinsic complex masked sum not within standard bounds")
        if (allocated(error)) return
        call check(error, err(2) < tolerance_kahan, "Kahan complex masked sum not within theoretical bounds")
        if (allocated(error)) return
        call check(error, err(3) < tolerance_chunked, "chunked complex masked sum not within theoretical bounds")
        if (allocated(error)) return
    end block
    block
        complex(dp), allocatable :: x(:)
        real(dp), parameter :: total_sum = 4*atan(1._dp)
        real(dp) :: tolerance_kahan, tolerance_chunked, tolerance_standard
        real(dp) :: err(ncalc)
        complex(dp) :: xsum(ncalc)
        logical, allocatable :: mask(:), nmask(:)

        allocate(x(n))
        do i = 1, n
            x(i) = (8*atan(1._dp)*(real(i,kind=dp)-0.5_dp)/n**2)*cmplx(1._dp,1._dp)
        end do
        
        allocate(mask(n),source=.false.); mask(1:n:2) = .true.
        allocate(nmask(n)); nmask = .not.mask
        ! scramble array
        do i = 1, n
            call random_number(u) 
            j = 1 + floor(n*u)
            call swap( x(i), x(j) )
            call swap( mask(i), mask(j) )
            call swap( nmask(i), nmask(j) )
        end do
        
        xsum(1) = sum(x)        ! compiler intrinsic
        xsum(2) = stdlib_sum_kahan(x) ! chunked Kahan summation
        xsum(3) = stdlib_sum(x)       ! chunked summation
        err(1:ncalc) = abs(1._dp-(xsum(1:ncalc)%re)/total_sum)
        
        tolerance_standard = standard_summation_tolerance(n, epsilon(1._dp), total_sum)
        tolerance_kahan = kahan_summation_tolerance(n, epsilon(1._dp), total_sum)
        tolerance_chunked = chunked_summation_tolerance(n, epsilon(1._dp), total_sum)
        
        call check(error, err(1) < tolerance_standard, "compiler intrinsic complex sum not within standard bounds")
        if (allocated(error)) return
        call check(error, err(2) < tolerance_kahan, "Kahan complex sum not within theoretical bounds")
        if (allocated(error)) return
        call check(error, err(3) < tolerance_chunked, "chunked complex sum not within theoretical bounds")
        if (allocated(error)) return

        xsum(1) = sum(x,mask)+sum(x,nmask) ! compiler intrinsic
        xsum(2) = stdlib_sum_kahan(x,mask)+stdlib_sum_kahan(x,nmask) ! chunked Kahan summation
        xsum(3) = stdlib_sum(x,mask)+stdlib_sum(x,nmask) ! chunked summation
        err(1:ncalc) = abs(1._dp-(xsum(1:ncalc)%re)/total_sum)
        
        call check(error, err(1) < tolerance_standard, "compiler intrinsic complex masked sum not within standard bounds")
        if (allocated(error)) return
        call check(error, err(2) < tolerance_kahan, "Kahan complex masked sum not within theoretical bounds")
        if (allocated(error)) return
        call check(error, err(3) < tolerance_chunked, "chunked complex masked sum not within theoretical bounds")
        if (allocated(error)) return
    end block

    ndarray : block
        use stdlib_strings, only: to_string
        real(sp), allocatable :: x(:,:,:)
        real(sp) :: tolerance_kahan, tolerance_chunked
        integer :: i

        allocate(x(100,100,10))
        call random_number(x)
        tolerance_chunked = chunked_summation_tolerance(size(x), epsilon(1._sp), sum(x))
        tolerance_kahan = kahan_summation_tolerance(size(x), epsilon(1._sp), sum(x))
        !> sum all elements
        call check(error, abs( sum(x) - stdlib_sum(x) ) < tolerance_chunked, "KO: full ndarray stdlib_sum" )
        if (allocated(error)) return

        call check(error, abs( sum(x) - stdlib_sum_kahan(x) ) < tolerance_kahan, "KO: full ndarray stdlib_sum_kahan" )
        if (allocated(error)) return

        !> sum over specific rank dim
        do i = 1, rank(x)
            call check(error, norm2( sum(x,dim=i) - stdlib_sum(x,dim=i) ) < tolerance_chunked ,&
                        "KO: ndarray stdlib_sum over dim "//to_string(i) )
            if (allocated(error)) return

            call check(error, norm2( sum(x,dim=i) - stdlib_sum_kahan(x,dim=i) ) < tolerance_kahan ,&
                        "KO: ndarray stdlib_sum_kahan over dim "//to_string(i) )
            if (allocated(error)) return
        end do
    end block ndarray

end subroutine

subroutine test_dot_product(error)
    !> Error handling
    type(error_type), allocatable, intent(out) :: error

    !> Internal parameters and variables
    integer, parameter :: n = 1e3, ncalc = 3
    real(sp) :: u
    integer :: iter, i, j
    !====================================================================================
    block
        real(sp), allocatable :: x(:)
        real(sp), parameter :: total_sum = 4*atan(1._sp)
        real(sp) :: tolerance_kahan, tolerance_chunked, tolerance_standard
        real(sp) :: xsum(ncalc), err(ncalc)

        allocate(x(n))
        do i = 1, n 
            x(i) = 2*sqrt( 2*atan(1._sp)*(real(i,kind=sp)-0.5_sp) )/n
        end do
        ! scramble array
        do i = 1, n
            call random_number(u) 
            j = 1 + floor(n*u)
            call swap( x(i), x(j) )
        end do
        
        xsum(1) = dot_product(x,x) ! compiler intrinsic
        xsum(2) = stdlib_dot_product_kahan(x,x) ! chunked Kahan summation
        xsum(3) = stdlib_dot_product(x,x)       ! chunked summation
        err(1:ncalc) = abs(1._sp-xsum(1:ncalc)/total_sum)
        
        tolerance_standard = standard_dot_product_tolerance(n, epsilon(1._sp), total_sum)
        tolerance_kahan = kahan_dot_product_tolerance(n, epsilon(1._sp), total_sum)
        tolerance_chunked = chunked_dot_product_tolerance(n, epsilon(1._sp), total_sum)
        
        call check(error, err(1) < tolerance_standard, "compiler intrinsic dot_product not within standard bounds")
        if (allocated(error)) return
        call check(error, err(2) < tolerance_kahan, "Kahan dot_product not within theoretical bounds")
        if (allocated(error)) return
        call check(error, err(3) < tolerance_chunked, "chunked dot_product not within theoretical bounds")
        if (allocated(error)) return
    end block
    block
        real(dp), allocatable :: x(:)
        real(dp), parameter :: total_sum = 4*atan(1._dp)
        real(dp) :: tolerance_kahan, tolerance_chunked, tolerance_standard
        real(dp) :: xsum(ncalc), err(ncalc)

        allocate(x(n))
        do i = 1, n 
            x(i) = 2*sqrt( 2*atan(1._dp)*(real(i,kind=dp)-0.5_dp) )/n
        end do
        ! scramble array
        do i = 1, n
            call random_number(u) 
            j = 1 + floor(n*u)
            call swap( x(i), x(j) )
        end do
        
        xsum(1) = dot_product(x,x) ! compiler intrinsic
        xsum(2) = stdlib_dot_product_kahan(x,x) ! chunked Kahan summation
        xsum(3) = stdlib_dot_product(x,x)       ! chunked summation
        err(1:ncalc) = abs(1._dp-xsum(1:ncalc)/total_sum)
        
        tolerance_standard = standard_dot_product_tolerance(n, epsilon(1._dp), total_sum)
        tolerance_kahan = kahan_dot_product_tolerance(n, epsilon(1._dp), total_sum)
        tolerance_chunked = chunked_dot_product_tolerance(n, epsilon(1._dp), total_sum)
        
        call check(error, err(1) < tolerance_standard, "compiler intrinsic dot_product not within standard bounds")
        if (allocated(error)) return
        call check(error, err(2) < tolerance_kahan, "Kahan dot_product not within theoretical bounds")
        if (allocated(error)) return
        call check(error, err(3) < tolerance_chunked, "chunked dot_product not within theoretical bounds")
        if (allocated(error)) return
    end block

    block
        complex(sp), allocatable :: x(:)
        real(sp), parameter :: total_sum = 4*atan(1._sp)
        real(sp) :: tolerance_kahan, tolerance_chunked, tolerance_standard
        real(sp) :: err(ncalc)
        complex(sp) :: xsum(ncalc)

        allocate(x(n))
        do i = 1, n
            x(i) = ( 2*sqrt( 2*atan(1._sp)*(real(i,kind=sp)-0.5_sp) ) / n )*cmplx(1._sp,1._sp)
        end do
        ! scramble array
        do i = 1, n
            call random_number(u) 
            j = 1 + floor(n*u)
            call swap( x(i), x(j) )
        end do
        
        xsum(1) = dot_product(x,x) ! compiler intrinsic
        xsum(2) = stdlib_dot_product_kahan(x,x) ! chunked Kahan dot_product
        xsum(3) = stdlib_dot_product(x,x)       ! chunked dot_product
        err(1:ncalc) = abs(1._sp-xsum(1:ncalc)%re/(2*total_sum))
        
        tolerance_standard = standard_dot_product_tolerance(n, epsilon(1._sp), 2*total_sum)
        tolerance_kahan = kahan_dot_product_tolerance(n, epsilon(1._sp), 2*total_sum)
        tolerance_chunked = chunked_dot_product_tolerance(n, epsilon(1._sp), 2*total_sum)
         
        call check(error, err(1) < tolerance_standard, "compiler intrinsic complex dot_product not within standard bounds")
        if (allocated(error)) return
        call check(error, err(2) < tolerance_kahan, "Kahan complex dot_product not within theoretical bounds")
        if (allocated(error)) return
        call check(error, err(3) < tolerance_chunked, "chunked complex dot_product not within theoretical bounds")
        if (allocated(error)) return
    end block

    block ! test for https://github.com/fortran-lang/stdlib/issues/1016
        complex(sp) :: x(128), y(128)
        real(sp) :: z(128,2)
        real(sp) :: tolerance
        real(sp) :: err(2)
        complex(sp) :: p(3)

        call random_number(z)
        x%re = z(:, 1); x%im = z(:, 2)
        call random_number(z)
        y%re = z(:, 1); y%im = z(:, 2)
        
        tolerance = complex_dot_product_tolerance(128, epsilon(1._sp))
        
        p(1) = dot_product(x,y) ! compiler intrinsic
        p(2) = stdlib_dot_product_kahan(x,y) ! chunked Kahan dot_product
        p(3) = stdlib_dot_product(x,y)       ! chunked dot_product
        err(1:2) = sqrt((p(2:3)%re - p(1)%re)**2 + (p(2:3)%im - p(1)%im)**2)
        
        call check(error, all(err(:)<tolerance) , "complex dot_product does not conform to the standard" )
        if (allocated(error)) return
    end block
    block
        complex(dp), allocatable :: x(:)
        real(dp), parameter :: total_sum = 4*atan(1._dp)
        real(dp) :: tolerance_kahan, tolerance_chunked, tolerance_standard
        real(dp) :: err(ncalc)
        complex(dp) :: xsum(ncalc)

        allocate(x(n))
        do i = 1, n
            x(i) = ( 2*sqrt( 2*atan(1._dp)*(real(i,kind=dp)-0.5_dp) ) / n )*cmplx(1._dp,1._dp)
        end do
        ! scramble array
        do i = 1, n
            call random_number(u) 
            j = 1 + floor(n*u)
            call swap( x(i), x(j) )
        end do
        
        xsum(1) = dot_product(x,x) ! compiler intrinsic
        xsum(2) = stdlib_dot_product_kahan(x,x) ! chunked Kahan dot_product
        xsum(3) = stdlib_dot_product(x,x)       ! chunked dot_product
        err(1:ncalc) = abs(1._dp-xsum(1:ncalc)%re/(2*total_sum))
        
        tolerance_standard = standard_dot_product_tolerance(n, epsilon(1._dp), 2*total_sum)
        tolerance_kahan = kahan_dot_product_tolerance(n, epsilon(1._dp), 2*total_sum)
        tolerance_chunked = chunked_dot_product_tolerance(n, epsilon(1._dp), 2*total_sum)
         
        call check(error, err(1) < tolerance_standard, "compiler intrinsic complex dot_product not within standard bounds")
        if (allocated(error)) return
        call check(error, err(2) < tolerance_kahan, "Kahan complex dot_product not within theoretical bounds")
        if (allocated(error)) return
        call check(error, err(3) < tolerance_chunked, "chunked complex dot_product not within theoretical bounds")
        if (allocated(error)) return
    end block

    block ! test for https://github.com/fortran-lang/stdlib/issues/1016
        complex(dp) :: x(128), y(128)
        real(dp) :: z(128,2)
        real(dp) :: tolerance
        real(dp) :: err(2)
        complex(dp) :: p(3)

        call random_number(z)
        x%re = z(:, 1); x%im = z(:, 2)
        call random_number(z)
        y%re = z(:, 1); y%im = z(:, 2)
        
        tolerance = complex_dot_product_tolerance(128, epsilon(1._dp))
        
        p(1) = dot_product(x,y) ! compiler intrinsic
        p(2) = stdlib_dot_product_kahan(x,y) ! chunked Kahan dot_product
        p(3) = stdlib_dot_product(x,y)       ! chunked dot_product
        err(1:2) = sqrt((p(2:3)%re - p(1)%re)**2 + (p(2:3)%im - p(1)%im)**2)
        
        call check(error, all(err(:)<tolerance) , "complex dot_product does not conform to the standard" )
        if (allocated(error)) return
    end block

end subroutine

subroutine test_reference_summation_validation(error)
    !> Error handling
    type(error_type), allocatable, intent(out) :: error

    block
        real(sp), allocatable :: x(:)
        real(sp) :: result_stdlib, result_kahan
        real(sp) :: reference_sum, relaxed_bound
        real(sp), parameter :: base_tolerance = epsilon(1._sp)
        integer :: i, n

        n = 1000
        allocate(x(n))
        
        do i = 1, n
            x(i) = 1.0_sp / real(i, sp)**2
        end do
        
        reference_sum = sum(x)
        
        result_stdlib = stdlib_sum(x)
        result_kahan = stdlib_sum_kahan(x)
        
        relaxed_bound = kahan_summation_tolerance(n, base_tolerance, reference_sum)
        
        call check(error, abs(result_kahan - reference_sum) < relaxed_bound, &
                  "Kahan summation should be within relaxed error bounds for harmonic series")
        if (allocated(error)) return
        
        call check(error, abs(result_stdlib - reference_sum) < relaxed_bound, &
                  "Stdlib summation should be within relaxed error bounds for harmonic series")
        if (allocated(error)) return
    end block
    block
        real(dp), allocatable :: x(:)
        real(dp) :: result_stdlib, result_kahan
        real(dp) :: reference_sum, relaxed_bound
        real(dp), parameter :: base_tolerance = epsilon(1._dp)
        integer :: i, n

        n = 1000
        allocate(x(n))
        
        do i = 1, n
            x(i) = 1.0_dp / real(i, dp)**2
        end do
        
        reference_sum = sum(x)
        
        result_stdlib = stdlib_sum(x)
        result_kahan = stdlib_sum_kahan(x)
        
        relaxed_bound = kahan_summation_tolerance(n, base_tolerance, reference_sum)
        
        call check(error, abs(result_kahan - reference_sum) < relaxed_bound, &
                  "Kahan summation should be within relaxed error bounds for harmonic series")
        if (allocated(error)) return
        
        call check(error, abs(result_stdlib - reference_sum) < relaxed_bound, &
                  "Stdlib summation should be within relaxed error bounds for harmonic series")
        if (allocated(error)) return
    end block

end subroutine

subroutine test_dot_product_reference_validation(error)
    !> Error handling
    type(error_type), allocatable, intent(out) :: error

    block
        real(sp), allocatable :: x(:), y(:)
        real(sp) :: result_stdlib, result_kahan
        real(sp) :: reference_dot, practical_bound
        real(sp), parameter :: base_tolerance = epsilon(1._sp)
        integer :: i, n

        n = 1000
        allocate(x(n), y(n))
        
        do i = 1, n
            x(i) = 1.0_sp / sqrt(real(i, sp))
            y(i) = 1.0_sp / sqrt(real(i, sp))
        end do
        
        reference_dot = dot_product(x, y)
        
        result_stdlib = stdlib_dot_product(x, y)
        result_kahan = stdlib_dot_product_kahan(x, y)
        
        practical_bound = kahan_dot_product_tolerance(n, base_tolerance, reference_dot)
        
        call check(error, abs(result_kahan - reference_dot) < practical_bound, &
                  "Kahan dot product should be within practical error bounds")
        if (allocated(error)) return
        
        call check(error, abs(result_stdlib - reference_dot) < practical_bound, &
                  "Stdlib dot product should be within practical error bounds")
        if (allocated(error)) return
    end block
    block
        real(dp), allocatable :: x(:), y(:)
        real(dp) :: result_stdlib, result_kahan
        real(dp) :: reference_dot, practical_bound
        real(dp), parameter :: base_tolerance = epsilon(1._dp)
        integer :: i, n

        n = 1000
        allocate(x(n), y(n))
        
        do i = 1, n
            x(i) = 1.0_dp / sqrt(real(i, dp))
            y(i) = 1.0_dp / sqrt(real(i, dp))
        end do
        
        reference_dot = dot_product(x, y)
        
        result_stdlib = stdlib_dot_product(x, y)
        result_kahan = stdlib_dot_product_kahan(x, y)
        
        practical_bound = kahan_dot_product_tolerance(n, base_tolerance, reference_dot)
        
        call check(error, abs(result_kahan - reference_dot) < practical_bound, &
                  "Kahan dot product should be within practical error bounds")
        if (allocated(error)) return
        
        call check(error, abs(result_stdlib - reference_dot) < practical_bound, &
                  "Stdlib dot product should be within practical error bounds")
        if (allocated(error)) return
    end block

end subroutine

subroutine test_kahan_superiority_demonstration(error)
    !> Error handling
    type(error_type), allocatable, intent(out) :: error

    block
        real(sp), allocatable :: x(:)
        real(sp) :: result_stdlib, result_kahan
        real(sp) :: expected, error_stdlib, error_kahan
        real(sp), parameter :: base_tolerance = epsilon(1._sp)
        integer :: i, n

        n = 10000
        allocate(x(n))
        
        do i = 1, n
            x(i) = ((-1.0_sp)**(i-1)) / real(i, sp)**3
        end do
        expected = sum(x)
        
        result_stdlib = stdlib_sum(x)
        result_kahan = stdlib_sum_kahan(x)
        
        error_stdlib = abs(expected - result_stdlib)
        error_kahan = abs(expected - result_kahan)
        
        call check(error, error_kahan <= error_stdlib, &
                  "Kahan summation should be at least as accurate as stdlib for alternating series")
        if (allocated(error)) return
        
        call check(error, error_kahan < kahan_summation_tolerance(n, base_tolerance, sum(abs(x))), &
                  "Kahan summation should handle alternating series within reasonable bounds")
        if (allocated(error)) return
    end block
    block
        real(dp), allocatable :: x(:)
        real(dp) :: result_stdlib, result_kahan
        real(dp) :: expected, error_stdlib, error_kahan
        real(dp), parameter :: base_tolerance = epsilon(1._dp)
        integer :: i, n

        n = 10000
        allocate(x(n))
        
        do i = 1, n
            x(i) = ((-1.0_dp)**(i-1)) / real(i, dp)**3
        end do
        expected = sum(x)
        
        result_stdlib = stdlib_sum(x)
        result_kahan = stdlib_sum_kahan(x)
        
        error_stdlib = abs(expected - result_stdlib)
        error_kahan = abs(expected - result_kahan)
        
        call check(error, error_kahan <= error_stdlib, &
                  "Kahan summation should be at least as accurate as stdlib for alternating series")
        if (allocated(error)) return
        
        call check(error, error_kahan < kahan_summation_tolerance(n, base_tolerance, sum(abs(x))), &
                  "Kahan summation should handle alternating series within reasonable bounds")
        if (allocated(error)) return
    end block

end subroutine

subroutine test_theoretical_error_bounds(error)
    !> Error handling
    type(error_type), allocatable, intent(out) :: error

    block
        real(sp), allocatable :: x(:)
        real(sp) :: result_kahan, exact_sum
        real(sp) :: practical_bound, actual_error
        real(sp), parameter :: base_tolerance = epsilon(1._sp)
        integer :: i, n

        n = 100
        allocate(x(n))
        
        do i = 1, n
            x(i) = 1.0_sp / real(i, sp)
        end do
        
        exact_sum = sum(x)
        result_kahan = stdlib_sum_kahan(x)
        actual_error = abs(result_kahan - exact_sum)
        
        practical_bound = kahan_summation_tolerance(n, base_tolerance, exact_sum)
        
        call check(error, actual_error <= practical_bound, &
                  "Kahan summation error should be within practical error bound")
        if (allocated(error)) return
        
        call check(error, result_kahan /= 0.0_sp, &
                  "Kahan summation should produce non-zero result for harmonic series")
        if (allocated(error)) return
    end block
    block
        real(dp), allocatable :: x(:)
        real(dp) :: result_kahan, exact_sum
        real(dp) :: practical_bound, actual_error
        real(dp), parameter :: base_tolerance = epsilon(1._dp)
        integer :: i, n

        n = 100
        allocate(x(n))
        
        do i = 1, n
            x(i) = 1.0_dp / real(i, dp)
        end do
        
        exact_sum = sum(x)
        result_kahan = stdlib_sum_kahan(x)
        actual_error = abs(result_kahan - exact_sum)
        
        practical_bound = kahan_summation_tolerance(n, base_tolerance, exact_sum)
        
        call check(error, actual_error <= practical_bound, &
                  "Kahan summation error should be within practical error bound")
        if (allocated(error)) return
        
        call check(error, result_kahan /= 0.0_dp, &
                  "Kahan summation should produce non-zero result for harmonic series")
        if (allocated(error)) return
    end block

end subroutine

subroutine test_complex_conjugation_accuracy(error)
    !> Error handling
    type(error_type), allocatable, intent(out) :: error

    block
        complex(sp), allocatable :: x(:), y(:)
        complex(sp) :: result_kahan, result_intrinsic
        complex(sp) :: expected_result
        real(sp), parameter :: base_tolerance = epsilon(1._sp)
        real(sp) :: error_magnitude, practical_bound
        integer :: i, n

        n = 10
        allocate(x(n), y(n))
        
        do i = 1, n
            x(i) = cmplx(real(i, sp), real(i, sp))
            y(i) = x(i)
        end do
        
        result_kahan = stdlib_dot_product_kahan(x, y)
        result_intrinsic = dot_product(x, y)
        expected_result = sum(conjg(x) * y)
        
        error_magnitude = abs(result_kahan - expected_result)
        practical_bound = complex_dot_product_tolerance(n, base_tolerance)
        
        call check(error, error_magnitude < practical_bound, &
                  "Complex dot product should match expected result")
        if (allocated(error)) return
        
        call check(error, abs(result_kahan - result_intrinsic) < practical_bound, &
                  "Kahan and intrinsic complex dot products should agree")
        if (allocated(error)) return
        
        call check(error, abs(expected_result) > 0.0_sp, &
                  "Complex dot product test should have non-zero expected result")
        if (allocated(error)) return
    end block
    block
        complex(dp), allocatable :: x(:), y(:)
        complex(dp) :: result_kahan, result_intrinsic
        complex(dp) :: expected_result
        real(dp), parameter :: base_tolerance = epsilon(1._dp)
        real(dp) :: error_magnitude, practical_bound
        integer :: i, n

        n = 10
        allocate(x(n), y(n))
        
        do i = 1, n
            x(i) = cmplx(real(i, dp), real(i, dp))
            y(i) = x(i)
        end do
        
        result_kahan = stdlib_dot_product_kahan(x, y)
        result_intrinsic = dot_product(x, y)
        expected_result = sum(conjg(x) * y)
        
        error_magnitude = abs(result_kahan - expected_result)
        practical_bound = complex_dot_product_tolerance(n, base_tolerance)
        
        call check(error, error_magnitude < practical_bound, &
                  "Complex dot product should match expected result")
        if (allocated(error)) return
        
        call check(error, abs(result_kahan - result_intrinsic) < practical_bound, &
                  "Kahan and intrinsic complex dot products should agree")
        if (allocated(error)) return
        
        call check(error, abs(expected_result) > 0.0_dp, &
                  "Complex dot product test should have non-zero expected result")
        if (allocated(error)) return
    end block

end subroutine
    
end module test_intrinsics

program tester
    use, intrinsic :: iso_fortran_env, only : error_unit
    use testdrive, only : run_testsuite, new_testsuite, testsuite_type
    use test_intrinsics, only : collect_suite
    implicit none
    integer :: stat, is
    type(testsuite_type), allocatable :: testsuites(:)
    character(len=*), parameter :: fmt = '("#", *(1x, a))'

    stat = 0

    testsuites = [ &
        new_testsuite("sparse", collect_suite) &
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
