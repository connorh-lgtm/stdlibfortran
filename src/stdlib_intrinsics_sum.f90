
submodule(stdlib_intrinsics) stdlib_intrinsics_sum
    !! ([Specification](../page/specs/stdlib_intrinsics.html))
    use stdlib_kinds
    use stdlib_constants
    implicit none

    integer, parameter :: ilp = int64
    integer(ilp), parameter :: optimal_chunk = 64_ilp
    
contains

!================= 1D Base implementations ============
! This implementation is based on https://github.com/jalvesz/fast_math
pure module function stdlib_sum_1d_int8(a) result(s)
    integer(int8), intent(in) :: a(:)
    integer(int8) :: s
    integer(int8) :: abatch(optimal_chunk)
    integer(ilp) :: i, n, r
    ! -----------------------------
    n = size(a,kind=ilp)
    r = mod(n,optimal_chunk)
    
    abatch(1:r)       = a(1:r)
    abatch(r+1:optimal_chunk) = zero_int8
    do i = r+1, n-r, optimal_chunk
     abatch(1:optimal_chunk) = abatch(1:optimal_chunk) + a(i:i+optimal_chunk-1)
    end do

    s = zero_int8
    do i = 1, optimal_chunk/2
      s = s + abatch(i)+abatch(optimal_chunk/2+i)
    end do
end function

pure module function stdlib_sum_1d_int8_mask(a,mask) result(s)
    integer(int8), intent(in) :: a(:)
    logical, intent(in) :: mask(:)
    integer(int8) :: s
    integer(int8) :: abatch(optimal_chunk)
    integer(ilp) :: i, n, r
    ! -----------------------------
    n = size(a,kind=ilp)
    r = mod(n,optimal_chunk)

    abatch(1:r)       = merge( zero_int8 , a(1:r) , mask(1:r) )
    abatch(r+1:optimal_chunk) = zero_int8
    do i = r+1, n-r, optimal_chunk
        abatch(1:optimal_chunk) = abatch(1:optimal_chunk) + merge( zero_int8 , a(i:i+optimal_chunk-1), mask(i:i+optimal_chunk-1) )
    end do
    
    s = zero_int8
    do i = 1, optimal_chunk/2
        s = s + abatch(i)+abatch(optimal_chunk/2+i)
    end do
end function
pure module function stdlib_sum_1d_int16(a) result(s)
    integer(int16), intent(in) :: a(:)
    integer(int16) :: s
    integer(int16) :: abatch(optimal_chunk)
    integer(ilp) :: i, n, r
    ! -----------------------------
    n = size(a,kind=ilp)
    r = mod(n,optimal_chunk)
    
    abatch(1:r)       = a(1:r)
    abatch(r+1:optimal_chunk) = zero_int16
    do i = r+1, n-r, optimal_chunk
     abatch(1:optimal_chunk) = abatch(1:optimal_chunk) + a(i:i+optimal_chunk-1)
    end do

    s = zero_int16
    do i = 1, optimal_chunk/2
      s = s + abatch(i)+abatch(optimal_chunk/2+i)
    end do
end function

pure module function stdlib_sum_1d_int16_mask(a,mask) result(s)
    integer(int16), intent(in) :: a(:)
    logical, intent(in) :: mask(:)
    integer(int16) :: s
    integer(int16) :: abatch(optimal_chunk)
    integer(ilp) :: i, n, r
    ! -----------------------------
    n = size(a,kind=ilp)
    r = mod(n,optimal_chunk)

    abatch(1:r)       = merge( zero_int16 , a(1:r) , mask(1:r) )
    abatch(r+1:optimal_chunk) = zero_int16
    do i = r+1, n-r, optimal_chunk
        abatch(1:optimal_chunk) = abatch(1:optimal_chunk) + merge( zero_int16 , a(i:i+optimal_chunk-1), mask(i:i+optimal_chunk-1) )
    end do
    
    s = zero_int16
    do i = 1, optimal_chunk/2
        s = s + abatch(i)+abatch(optimal_chunk/2+i)
    end do
end function
pure module function stdlib_sum_1d_int32(a) result(s)
    integer(int32), intent(in) :: a(:)
    integer(int32) :: s
    integer(int32) :: abatch(optimal_chunk)
    integer(ilp) :: i, n, r
    ! -----------------------------
    n = size(a,kind=ilp)
    r = mod(n,optimal_chunk)
    
    abatch(1:r)       = a(1:r)
    abatch(r+1:optimal_chunk) = zero_int32
    do i = r+1, n-r, optimal_chunk
     abatch(1:optimal_chunk) = abatch(1:optimal_chunk) + a(i:i+optimal_chunk-1)
    end do

    s = zero_int32
    do i = 1, optimal_chunk/2
      s = s + abatch(i)+abatch(optimal_chunk/2+i)
    end do
end function

pure module function stdlib_sum_1d_int32_mask(a,mask) result(s)
    integer(int32), intent(in) :: a(:)
    logical, intent(in) :: mask(:)
    integer(int32) :: s
    integer(int32) :: abatch(optimal_chunk)
    integer(ilp) :: i, n, r
    ! -----------------------------
    n = size(a,kind=ilp)
    r = mod(n,optimal_chunk)

    abatch(1:r)       = merge( zero_int32 , a(1:r) , mask(1:r) )
    abatch(r+1:optimal_chunk) = zero_int32
    do i = r+1, n-r, optimal_chunk
        abatch(1:optimal_chunk) = abatch(1:optimal_chunk) + merge( zero_int32 , a(i:i+optimal_chunk-1), mask(i:i+optimal_chunk-1) )
    end do
    
    s = zero_int32
    do i = 1, optimal_chunk/2
        s = s + abatch(i)+abatch(optimal_chunk/2+i)
    end do
end function
pure module function stdlib_sum_1d_int64(a) result(s)
    integer(int64), intent(in) :: a(:)
    integer(int64) :: s
    integer(int64) :: abatch(optimal_chunk)
    integer(ilp) :: i, n, r
    ! -----------------------------
    n = size(a,kind=ilp)
    r = mod(n,optimal_chunk)
    
    abatch(1:r)       = a(1:r)
    abatch(r+1:optimal_chunk) = zero_int64
    do i = r+1, n-r, optimal_chunk
     abatch(1:optimal_chunk) = abatch(1:optimal_chunk) + a(i:i+optimal_chunk-1)
    end do

    s = zero_int64
    do i = 1, optimal_chunk/2
      s = s + abatch(i)+abatch(optimal_chunk/2+i)
    end do
end function

pure module function stdlib_sum_1d_int64_mask(a,mask) result(s)
    integer(int64), intent(in) :: a(:)
    logical, intent(in) :: mask(:)
    integer(int64) :: s
    integer(int64) :: abatch(optimal_chunk)
    integer(ilp) :: i, n, r
    ! -----------------------------
    n = size(a,kind=ilp)
    r = mod(n,optimal_chunk)

    abatch(1:r)       = merge( zero_int64 , a(1:r) , mask(1:r) )
    abatch(r+1:optimal_chunk) = zero_int64
    do i = r+1, n-r, optimal_chunk
        abatch(1:optimal_chunk) = abatch(1:optimal_chunk) + merge( zero_int64 , a(i:i+optimal_chunk-1), mask(i:i+optimal_chunk-1) )
    end do
    
    s = zero_int64
    do i = 1, optimal_chunk/2
        s = s + abatch(i)+abatch(optimal_chunk/2+i)
    end do
end function
pure module function stdlib_sum_1d_sp(a) result(s)
    real(sp), intent(in) :: a(:)
    real(sp) :: s
    real(sp) :: abatch(optimal_chunk)
    integer(ilp) :: i, n, r
    ! -----------------------------
    n = size(a,kind=ilp)
    r = mod(n,optimal_chunk)
    
    abatch(1:r)       = a(1:r)
    abatch(r+1:optimal_chunk) = zero_sp
    do i = r+1, n-r, optimal_chunk
     abatch(1:optimal_chunk) = abatch(1:optimal_chunk) + a(i:i+optimal_chunk-1)
    end do

    s = zero_sp
    do i = 1, optimal_chunk/2
      s = s + abatch(i)+abatch(optimal_chunk/2+i)
    end do
end function

pure module function stdlib_sum_1d_sp_mask(a,mask) result(s)
    real(sp), intent(in) :: a(:)
    logical, intent(in) :: mask(:)
    real(sp) :: s
    real(sp) :: abatch(optimal_chunk)
    integer(ilp) :: i, n, r
    ! -----------------------------
    n = size(a,kind=ilp)
    r = mod(n,optimal_chunk)

    abatch(1:r)       = merge( zero_sp , a(1:r) , mask(1:r) )
    abatch(r+1:optimal_chunk) = zero_sp
    do i = r+1, n-r, optimal_chunk
        abatch(1:optimal_chunk) = abatch(1:optimal_chunk) + merge( zero_sp , a(i:i+optimal_chunk-1), mask(i:i+optimal_chunk-1) )
    end do
    
    s = zero_sp
    do i = 1, optimal_chunk/2
        s = s + abatch(i)+abatch(optimal_chunk/2+i)
    end do
end function
pure module function stdlib_sum_1d_dp(a) result(s)
    real(dp), intent(in) :: a(:)
    real(dp) :: s
    real(dp) :: abatch(optimal_chunk)
    integer(ilp) :: i, n, r
    ! -----------------------------
    n = size(a,kind=ilp)
    r = mod(n,optimal_chunk)
    
    abatch(1:r)       = a(1:r)
    abatch(r+1:optimal_chunk) = zero_dp
    do i = r+1, n-r, optimal_chunk
     abatch(1:optimal_chunk) = abatch(1:optimal_chunk) + a(i:i+optimal_chunk-1)
    end do

    s = zero_dp
    do i = 1, optimal_chunk/2
      s = s + abatch(i)+abatch(optimal_chunk/2+i)
    end do
end function

pure module function stdlib_sum_1d_dp_mask(a,mask) result(s)
    real(dp), intent(in) :: a(:)
    logical, intent(in) :: mask(:)
    real(dp) :: s
    real(dp) :: abatch(optimal_chunk)
    integer(ilp) :: i, n, r
    ! -----------------------------
    n = size(a,kind=ilp)
    r = mod(n,optimal_chunk)

    abatch(1:r)       = merge( zero_dp , a(1:r) , mask(1:r) )
    abatch(r+1:optimal_chunk) = zero_dp
    do i = r+1, n-r, optimal_chunk
        abatch(1:optimal_chunk) = abatch(1:optimal_chunk) + merge( zero_dp , a(i:i+optimal_chunk-1), mask(i:i+optimal_chunk-1) )
    end do
    
    s = zero_dp
    do i = 1, optimal_chunk/2
        s = s + abatch(i)+abatch(optimal_chunk/2+i)
    end do
end function
pure module function stdlib_sum_1d_csp(a) result(s)
    complex(sp), intent(in) :: a(:)
    complex(sp) :: s
    complex(sp) :: abatch(optimal_chunk)
    integer(ilp) :: i, n, r
    ! -----------------------------
    n = size(a,kind=ilp)
    r = mod(n,optimal_chunk)
    
    abatch(1:r)       = a(1:r)
    abatch(r+1:optimal_chunk) = zero_csp
    do i = r+1, n-r, optimal_chunk
     abatch(1:optimal_chunk) = abatch(1:optimal_chunk) + a(i:i+optimal_chunk-1)
    end do

    s = zero_csp
    do i = 1, optimal_chunk/2
      s = s + abatch(i)+abatch(optimal_chunk/2+i)
    end do
end function

pure module function stdlib_sum_1d_csp_mask(a,mask) result(s)
    complex(sp), intent(in) :: a(:)
    logical, intent(in) :: mask(:)
    complex(sp) :: s
    complex(sp) :: abatch(optimal_chunk)
    integer(ilp) :: i, n, r
    ! -----------------------------
    n = size(a,kind=ilp)
    r = mod(n,optimal_chunk)

    abatch(1:r)       = merge( zero_csp , a(1:r) , mask(1:r) )
    abatch(r+1:optimal_chunk) = zero_csp
    do i = r+1, n-r, optimal_chunk
        abatch(1:optimal_chunk) = abatch(1:optimal_chunk) + merge( zero_csp , a(i:i+optimal_chunk-1), mask(i:i+optimal_chunk-1) )
    end do
    
    s = zero_csp
    do i = 1, optimal_chunk/2
        s = s + abatch(i)+abatch(optimal_chunk/2+i)
    end do
end function
pure module function stdlib_sum_1d_cdp(a) result(s)
    complex(dp), intent(in) :: a(:)
    complex(dp) :: s
    complex(dp) :: abatch(optimal_chunk)
    integer(ilp) :: i, n, r
    ! -----------------------------
    n = size(a,kind=ilp)
    r = mod(n,optimal_chunk)
    
    abatch(1:r)       = a(1:r)
    abatch(r+1:optimal_chunk) = zero_cdp
    do i = r+1, n-r, optimal_chunk
     abatch(1:optimal_chunk) = abatch(1:optimal_chunk) + a(i:i+optimal_chunk-1)
    end do

    s = zero_cdp
    do i = 1, optimal_chunk/2
      s = s + abatch(i)+abatch(optimal_chunk/2+i)
    end do
end function

pure module function stdlib_sum_1d_cdp_mask(a,mask) result(s)
    complex(dp), intent(in) :: a(:)
    logical, intent(in) :: mask(:)
    complex(dp) :: s
    complex(dp) :: abatch(optimal_chunk)
    integer(ilp) :: i, n, r
    ! -----------------------------
    n = size(a,kind=ilp)
    r = mod(n,optimal_chunk)

    abatch(1:r)       = merge( zero_cdp , a(1:r) , mask(1:r) )
    abatch(r+1:optimal_chunk) = zero_cdp
    do i = r+1, n-r, optimal_chunk
        abatch(1:optimal_chunk) = abatch(1:optimal_chunk) + merge( zero_cdp , a(i:i+optimal_chunk-1), mask(i:i+optimal_chunk-1) )
    end do
    
    s = zero_cdp
    do i = 1, optimal_chunk/2
        s = s + abatch(i)+abatch(optimal_chunk/2+i)
    end do
end function

pure module function stdlib_sum_kahan_1d_sp(a) result(s)
    real(sp), intent(in) :: a(:)
    real(sp) :: s
    real(sp) :: sbatch(optimal_chunk)
    real(sp) :: cbatch(optimal_chunk)
    integer(ilp) :: i, n, r
    ! -----------------------------
    n = size(a,kind=ilp)
    r = mod(n,optimal_chunk)

    sbatch(1:r) = a(1:r)
    sbatch(r+1:optimal_chunk)  = zero_sp
    cbatch = zero_sp
    do i = r+1, n-r, optimal_chunk
      call kahan_kernel( a(i:i+optimal_chunk-1) , sbatch(1:optimal_chunk) , cbatch(1:optimal_chunk) )
    end do 

    s = zero_sp
    do i = 1,optimal_chunk
        call kahan_kernel( sbatch(i) , s , cbatch(i) )
    end do
end function

pure module function stdlib_sum_kahan_1d_sp_mask(a,mask) result(s)
    real(sp), intent(in) :: a(:)
    logical, intent(in) :: mask(:)
    real(sp) :: s
    real(sp) :: sbatch(optimal_chunk)
    real(sp) :: cbatch(optimal_chunk)
    integer(ilp) :: i, n, r
    ! -----------------------------
    n = size(a,kind=ilp)
    r = mod(n,optimal_chunk)
    
    sbatch(1:r) = merge( zero_sp , a(1:r) , mask(1:r) )
    sbatch(r+1:optimal_chunk)  = zero_sp
    cbatch = zero_sp
    do i = r+1, n-r, optimal_chunk
      call kahan_kernel( a(i:i+optimal_chunk-1) , sbatch(1:optimal_chunk) , cbatch(1:optimal_chunk), mask(i:i+optimal_chunk-1) )
    end do 

    s = zero_sp
    do i = 1,optimal_chunk
        call kahan_kernel( sbatch(i) , s , cbatch(i) )
    end do
end function
pure module function stdlib_sum_kahan_1d_dp(a) result(s)
    real(dp), intent(in) :: a(:)
    real(dp) :: s
    real(dp) :: sbatch(optimal_chunk)
    real(dp) :: cbatch(optimal_chunk)
    integer(ilp) :: i, n, r
    ! -----------------------------
    n = size(a,kind=ilp)
    r = mod(n,optimal_chunk)

    sbatch(1:r) = a(1:r)
    sbatch(r+1:optimal_chunk)  = zero_dp
    cbatch = zero_dp
    do i = r+1, n-r, optimal_chunk
      call kahan_kernel( a(i:i+optimal_chunk-1) , sbatch(1:optimal_chunk) , cbatch(1:optimal_chunk) )
    end do 

    s = zero_dp
    do i = 1,optimal_chunk
        call kahan_kernel( sbatch(i) , s , cbatch(i) )
    end do
end function

pure module function stdlib_sum_kahan_1d_dp_mask(a,mask) result(s)
    real(dp), intent(in) :: a(:)
    logical, intent(in) :: mask(:)
    real(dp) :: s
    real(dp) :: sbatch(optimal_chunk)
    real(dp) :: cbatch(optimal_chunk)
    integer(ilp) :: i, n, r
    ! -----------------------------
    n = size(a,kind=ilp)
    r = mod(n,optimal_chunk)
    
    sbatch(1:r) = merge( zero_dp , a(1:r) , mask(1:r) )
    sbatch(r+1:optimal_chunk)  = zero_dp
    cbatch = zero_dp
    do i = r+1, n-r, optimal_chunk
      call kahan_kernel( a(i:i+optimal_chunk-1) , sbatch(1:optimal_chunk) , cbatch(1:optimal_chunk), mask(i:i+optimal_chunk-1) )
    end do 

    s = zero_dp
    do i = 1,optimal_chunk
        call kahan_kernel( sbatch(i) , s , cbatch(i) )
    end do
end function
pure module function stdlib_sum_kahan_1d_csp(a) result(s)
    complex(sp), intent(in) :: a(:)
    complex(sp) :: s
    complex(sp) :: sbatch(optimal_chunk)
    complex(sp) :: cbatch(optimal_chunk)
    integer(ilp) :: i, n, r
    ! -----------------------------
    n = size(a,kind=ilp)
    r = mod(n,optimal_chunk)

    sbatch(1:r) = a(1:r)
    sbatch(r+1:optimal_chunk)  = zero_csp
    cbatch = zero_csp
    do i = r+1, n-r, optimal_chunk
      call kahan_kernel( a(i:i+optimal_chunk-1) , sbatch(1:optimal_chunk) , cbatch(1:optimal_chunk) )
    end do 

    s = zero_csp
    do i = 1,optimal_chunk
        call kahan_kernel( sbatch(i) , s , cbatch(i) )
    end do
end function

pure module function stdlib_sum_kahan_1d_csp_mask(a,mask) result(s)
    complex(sp), intent(in) :: a(:)
    logical, intent(in) :: mask(:)
    complex(sp) :: s
    complex(sp) :: sbatch(optimal_chunk)
    complex(sp) :: cbatch(optimal_chunk)
    integer(ilp) :: i, n, r
    ! -----------------------------
    n = size(a,kind=ilp)
    r = mod(n,optimal_chunk)
    
    sbatch(1:r) = merge( zero_csp , a(1:r) , mask(1:r) )
    sbatch(r+1:optimal_chunk)  = zero_csp
    cbatch = zero_csp
    do i = r+1, n-r, optimal_chunk
      call kahan_kernel( a(i:i+optimal_chunk-1) , sbatch(1:optimal_chunk) , cbatch(1:optimal_chunk), mask(i:i+optimal_chunk-1) )
    end do 

    s = zero_csp
    do i = 1,optimal_chunk
        call kahan_kernel( sbatch(i) , s , cbatch(i) )
    end do
end function
pure module function stdlib_sum_kahan_1d_cdp(a) result(s)
    complex(dp), intent(in) :: a(:)
    complex(dp) :: s
    complex(dp) :: sbatch(optimal_chunk)
    complex(dp) :: cbatch(optimal_chunk)
    integer(ilp) :: i, n, r
    ! -----------------------------
    n = size(a,kind=ilp)
    r = mod(n,optimal_chunk)

    sbatch(1:r) = a(1:r)
    sbatch(r+1:optimal_chunk)  = zero_cdp
    cbatch = zero_cdp
    do i = r+1, n-r, optimal_chunk
      call kahan_kernel( a(i:i+optimal_chunk-1) , sbatch(1:optimal_chunk) , cbatch(1:optimal_chunk) )
    end do 

    s = zero_cdp
    do i = 1,optimal_chunk
        call kahan_kernel( sbatch(i) , s , cbatch(i) )
    end do
end function

pure module function stdlib_sum_kahan_1d_cdp_mask(a,mask) result(s)
    complex(dp), intent(in) :: a(:)
    logical, intent(in) :: mask(:)
    complex(dp) :: s
    complex(dp) :: sbatch(optimal_chunk)
    complex(dp) :: cbatch(optimal_chunk)
    integer(ilp) :: i, n, r
    ! -----------------------------
    n = size(a,kind=ilp)
    r = mod(n,optimal_chunk)
    
    sbatch(1:r) = merge( zero_cdp , a(1:r) , mask(1:r) )
    sbatch(r+1:optimal_chunk)  = zero_cdp
    cbatch = zero_cdp
    do i = r+1, n-r, optimal_chunk
      call kahan_kernel( a(i:i+optimal_chunk-1) , sbatch(1:optimal_chunk) , cbatch(1:optimal_chunk), mask(i:i+optimal_chunk-1) )
    end do 

    s = zero_cdp
    do i = 1,optimal_chunk
        call kahan_kernel( sbatch(i) , s , cbatch(i) )
    end do
end function

!================= Pairwise summation for extreme cases ============
pure recursive function stdlib_sum_pairwise_1d_sp(a, start_idx, end_idx) result(s)
    real(sp), intent(in) :: a(:)
    integer(ilp), intent(in) :: start_idx, end_idx
    real(sp) :: s
    integer(ilp) :: mid, n
    
    n = end_idx - start_idx + 1
    if (n <= 4) then
        ! For small arrays, use specialized extreme cancellation handling
        s = stdlib_sum_extreme_cancellation_sp(a(start_idx:end_idx))
    else
        mid = start_idx + n / 2 - 1
        s = stdlib_sum_pairwise_1d_sp(a, start_idx, mid) + &
            stdlib_sum_pairwise_1d_sp(a, mid + 1, end_idx)
    end if
end function

pure function stdlib_sum_extreme_cancellation_sp(a) result(s)
    real(sp), intent(in) :: a(:)
    real(sp) :: s
    integer :: n, i, j
    real(sp), allocatable :: temp(:)
    logical, allocatable :: used(:)
    real(sp) :: abs_vals(size(a))
    integer :: sorted_indices(size(a))
    
    n = size(a)
    if (n == 0) then
        s = zero_sp
        return
    end if
    
    if (n == 1) then
        s = a(1)
        return
    end if
    
    ! For small arrays, try to pair canceling values first
    allocate(temp(n), used(n))
    temp = a
    used = .false.
    
    ! Calculate absolute values for sorting
    do i = 1, n
        abs_vals(i) = abs(a(i))
    end do
    
    ! Simple sorting by magnitude (bubble sort for small arrays)
    do i = 1, n
        sorted_indices(i) = i
    end do
    do i = 1, n-1
        do j = i+1, n
            if (abs_vals(sorted_indices(i)) < abs_vals(sorted_indices(j))) then
                ! Swap indices
                call swap_indices(sorted_indices(i), sorted_indices(j))
            end if
        end do
    end do
    
    ! Try to find canceling pairs among largest values first
    s = zero_sp
    do i = 1, n-1
        if (used(sorted_indices(i))) cycle
        do j = i+1, n
            if (used(sorted_indices(j))) cycle
            ! Check if values approximately cancel
            if (abs(temp(sorted_indices(i)) + temp(sorted_indices(j))) < &
                abs(temp(sorted_indices(i))) * 1.0e-10_sp) then
                s = s + (temp(sorted_indices(i)) + temp(sorted_indices(j)))
                used(sorted_indices(i)) = .true.
                used(sorted_indices(j)) = .true.
                exit
            end if
        end do
    end do
    
    ! Add remaining unpaired values
    do i = 1, n
        if (.not. used(sorted_indices(i))) then
            s = s + temp(sorted_indices(i))
        end if
    end do
    
contains
    pure subroutine swap_indices(a, b)
        integer, intent(inout) :: a, b
        integer :: temp_idx
        temp_idx = a
        a = b
        b = temp_idx
    end subroutine
end function

pure function detect_extreme_magnitude_ratio_sp(a) result(is_extreme)
    real(sp), intent(in) :: a(:)
    logical :: is_extreme
    real(sp) :: max_val, min_val, ratio
    real(sp), parameter :: extreme_threshold = 1.0e10_sp
    real(sp) :: abs_vals(size(a))
    integer :: i
    
    ! Calculate absolute values for both real and complex
    do i = 1, size(a)
        abs_vals(i) = abs(a(i))
    end do
    
    max_val = maxval(abs_vals)
    min_val = minval(abs_vals, mask=abs_vals > tiny(1.0_sp))
    if (min_val <= tiny(1.0_sp)) then
        min_val = tiny(1.0_sp)
    end if
    ratio = max_val / min_val
    is_extreme = ratio > extreme_threshold
end function

pure module function stdlib_sum_adaptive_1d_sp(a) result(s)
    real(sp), intent(in) :: a(:)
    real(sp) :: s
    
    if (detect_extreme_magnitude_ratio_sp(a)) then
        s = stdlib_sum_pairwise_1d_sp(a, 1_ilp, size(a, kind=ilp))
    else
        s = stdlib_sum_kahan(a)
    end if
end function

pure module function stdlib_sum_adaptive_1d_sp_mask(a, mask) result(s)
    real(sp), intent(in) :: a(:)
    logical, intent(in) :: mask(:)
    real(sp) :: s
    real(sp), allocatable :: masked_a(:)
    integer(ilp) :: i, count_true
    
    count_true = count(mask)
    if (count_true == 0) then
        s = zero_sp
        return
    end if
    
    allocate(masked_a(count_true))
    count_true = 0
    do i = 1, size(a)
        if (mask(i)) then
            count_true = count_true + 1
            masked_a(count_true) = a(i)
        end if
    end do
    
    s = stdlib_sum_adaptive_1d_sp(masked_a)
end function
pure recursive function stdlib_sum_pairwise_1d_dp(a, start_idx, end_idx) result(s)
    real(dp), intent(in) :: a(:)
    integer(ilp), intent(in) :: start_idx, end_idx
    real(dp) :: s
    integer(ilp) :: mid, n
    
    n = end_idx - start_idx + 1
    if (n <= 4) then
        ! For small arrays, use specialized extreme cancellation handling
        s = stdlib_sum_extreme_cancellation_dp(a(start_idx:end_idx))
    else
        mid = start_idx + n / 2 - 1
        s = stdlib_sum_pairwise_1d_dp(a, start_idx, mid) + &
            stdlib_sum_pairwise_1d_dp(a, mid + 1, end_idx)
    end if
end function

pure function stdlib_sum_extreme_cancellation_dp(a) result(s)
    real(dp), intent(in) :: a(:)
    real(dp) :: s
    integer :: n, i, j
    real(dp), allocatable :: temp(:)
    logical, allocatable :: used(:)
    real(dp) :: abs_vals(size(a))
    integer :: sorted_indices(size(a))
    
    n = size(a)
    if (n == 0) then
        s = zero_dp
        return
    end if
    
    if (n == 1) then
        s = a(1)
        return
    end if
    
    ! For small arrays, try to pair canceling values first
    allocate(temp(n), used(n))
    temp = a
    used = .false.
    
    ! Calculate absolute values for sorting
    do i = 1, n
        abs_vals(i) = abs(a(i))
    end do
    
    ! Simple sorting by magnitude (bubble sort for small arrays)
    do i = 1, n
        sorted_indices(i) = i
    end do
    do i = 1, n-1
        do j = i+1, n
            if (abs_vals(sorted_indices(i)) < abs_vals(sorted_indices(j))) then
                ! Swap indices
                call swap_indices(sorted_indices(i), sorted_indices(j))
            end if
        end do
    end do
    
    ! Try to find canceling pairs among largest values first
    s = zero_dp
    do i = 1, n-1
        if (used(sorted_indices(i))) cycle
        do j = i+1, n
            if (used(sorted_indices(j))) cycle
            ! Check if values approximately cancel
            if (abs(temp(sorted_indices(i)) + temp(sorted_indices(j))) < &
                abs(temp(sorted_indices(i))) * 1.0e-10_dp) then
                s = s + (temp(sorted_indices(i)) + temp(sorted_indices(j)))
                used(sorted_indices(i)) = .true.
                used(sorted_indices(j)) = .true.
                exit
            end if
        end do
    end do
    
    ! Add remaining unpaired values
    do i = 1, n
        if (.not. used(sorted_indices(i))) then
            s = s + temp(sorted_indices(i))
        end if
    end do
    
contains
    pure subroutine swap_indices(a, b)
        integer, intent(inout) :: a, b
        integer :: temp_idx
        temp_idx = a
        a = b
        b = temp_idx
    end subroutine
end function

pure function detect_extreme_magnitude_ratio_dp(a) result(is_extreme)
    real(dp), intent(in) :: a(:)
    logical :: is_extreme
    real(dp) :: max_val, min_val, ratio
    real(dp), parameter :: extreme_threshold = 1.0e10_dp
    real(dp) :: abs_vals(size(a))
    integer :: i
    
    ! Calculate absolute values for both real and complex
    do i = 1, size(a)
        abs_vals(i) = abs(a(i))
    end do
    
    max_val = maxval(abs_vals)
    min_val = minval(abs_vals, mask=abs_vals > tiny(1.0_dp))
    if (min_val <= tiny(1.0_dp)) then
        min_val = tiny(1.0_dp)
    end if
    ratio = max_val / min_val
    is_extreme = ratio > extreme_threshold
end function

pure module function stdlib_sum_adaptive_1d_dp(a) result(s)
    real(dp), intent(in) :: a(:)
    real(dp) :: s
    
    if (detect_extreme_magnitude_ratio_dp(a)) then
        s = stdlib_sum_pairwise_1d_dp(a, 1_ilp, size(a, kind=ilp))
    else
        s = stdlib_sum_kahan(a)
    end if
end function

pure module function stdlib_sum_adaptive_1d_dp_mask(a, mask) result(s)
    real(dp), intent(in) :: a(:)
    logical, intent(in) :: mask(:)
    real(dp) :: s
    real(dp), allocatable :: masked_a(:)
    integer(ilp) :: i, count_true
    
    count_true = count(mask)
    if (count_true == 0) then
        s = zero_dp
        return
    end if
    
    allocate(masked_a(count_true))
    count_true = 0
    do i = 1, size(a)
        if (mask(i)) then
            count_true = count_true + 1
            masked_a(count_true) = a(i)
        end if
    end do
    
    s = stdlib_sum_adaptive_1d_dp(masked_a)
end function
pure recursive function stdlib_sum_pairwise_1d_csp(a, start_idx, end_idx) result(s)
    complex(sp), intent(in) :: a(:)
    integer(ilp), intent(in) :: start_idx, end_idx
    complex(sp) :: s
    integer(ilp) :: mid, n
    
    n = end_idx - start_idx + 1
    if (n <= 4) then
        ! For small arrays, use specialized extreme cancellation handling
        s = stdlib_sum_extreme_cancellation_csp(a(start_idx:end_idx))
    else
        mid = start_idx + n / 2 - 1
        s = stdlib_sum_pairwise_1d_csp(a, start_idx, mid) + &
            stdlib_sum_pairwise_1d_csp(a, mid + 1, end_idx)
    end if
end function

pure function stdlib_sum_extreme_cancellation_csp(a) result(s)
    complex(sp), intent(in) :: a(:)
    complex(sp) :: s
    integer :: n, i, j
    complex(sp), allocatable :: temp(:)
    logical, allocatable :: used(:)
    real(sp) :: abs_vals(size(a))
    integer :: sorted_indices(size(a))
    
    n = size(a)
    if (n == 0) then
        s = zero_csp
        return
    end if
    
    if (n == 1) then
        s = a(1)
        return
    end if
    
    ! For small arrays, try to pair canceling values first
    allocate(temp(n), used(n))
    temp = a
    used = .false.
    
    ! Calculate absolute values for sorting
    do i = 1, n
        abs_vals(i) = abs(a(i))
    end do
    
    ! Simple sorting by magnitude (bubble sort for small arrays)
    do i = 1, n
        sorted_indices(i) = i
    end do
    do i = 1, n-1
        do j = i+1, n
            if (abs_vals(sorted_indices(i)) < abs_vals(sorted_indices(j))) then
                ! Swap indices
                call swap_indices(sorted_indices(i), sorted_indices(j))
            end if
        end do
    end do
    
    ! Try to find canceling pairs among largest values first
    s = zero_csp
    do i = 1, n-1
        if (used(sorted_indices(i))) cycle
        do j = i+1, n
            if (used(sorted_indices(j))) cycle
            ! Check if values approximately cancel
            if (abs(temp(sorted_indices(i)) + temp(sorted_indices(j))) < &
                abs(temp(sorted_indices(i))) * 1.0e-10_sp) then
                s = s + (temp(sorted_indices(i)) + temp(sorted_indices(j)))
                used(sorted_indices(i)) = .true.
                used(sorted_indices(j)) = .true.
                exit
            end if
        end do
    end do
    
    ! Add remaining unpaired values
    do i = 1, n
        if (.not. used(sorted_indices(i))) then
            s = s + temp(sorted_indices(i))
        end if
    end do
    
contains
    pure subroutine swap_indices(a, b)
        integer, intent(inout) :: a, b
        integer :: temp_idx
        temp_idx = a
        a = b
        b = temp_idx
    end subroutine
end function

pure function detect_extreme_magnitude_ratio_csp(a) result(is_extreme)
    complex(sp), intent(in) :: a(:)
    logical :: is_extreme
    real(sp) :: max_val, min_val, ratio
    real(sp), parameter :: extreme_threshold = 1.0e10_sp
    real(sp) :: abs_vals(size(a))
    integer :: i
    
    ! Calculate absolute values for both real and complex
    do i = 1, size(a)
        abs_vals(i) = abs(a(i))
    end do
    
    max_val = maxval(abs_vals)
    min_val = minval(abs_vals, mask=abs_vals > tiny(1.0_sp))
    if (min_val <= tiny(1.0_sp)) then
        min_val = tiny(1.0_sp)
    end if
    ratio = max_val / min_val
    is_extreme = ratio > extreme_threshold
end function

pure module function stdlib_sum_adaptive_1d_csp(a) result(s)
    complex(sp), intent(in) :: a(:)
    complex(sp) :: s
    
    if (detect_extreme_magnitude_ratio_csp(a)) then
        s = stdlib_sum_pairwise_1d_csp(a, 1_ilp, size(a, kind=ilp))
    else
        s = stdlib_sum_kahan(a)
    end if
end function

pure module function stdlib_sum_adaptive_1d_csp_mask(a, mask) result(s)
    complex(sp), intent(in) :: a(:)
    logical, intent(in) :: mask(:)
    complex(sp) :: s
    complex(sp), allocatable :: masked_a(:)
    integer(ilp) :: i, count_true
    
    count_true = count(mask)
    if (count_true == 0) then
        s = zero_csp
        return
    end if
    
    allocate(masked_a(count_true))
    count_true = 0
    do i = 1, size(a)
        if (mask(i)) then
            count_true = count_true + 1
            masked_a(count_true) = a(i)
        end if
    end do
    
    s = stdlib_sum_adaptive_1d_csp(masked_a)
end function
pure recursive function stdlib_sum_pairwise_1d_cdp(a, start_idx, end_idx) result(s)
    complex(dp), intent(in) :: a(:)
    integer(ilp), intent(in) :: start_idx, end_idx
    complex(dp) :: s
    integer(ilp) :: mid, n
    
    n = end_idx - start_idx + 1
    if (n <= 4) then
        ! For small arrays, use specialized extreme cancellation handling
        s = stdlib_sum_extreme_cancellation_cdp(a(start_idx:end_idx))
    else
        mid = start_idx + n / 2 - 1
        s = stdlib_sum_pairwise_1d_cdp(a, start_idx, mid) + &
            stdlib_sum_pairwise_1d_cdp(a, mid + 1, end_idx)
    end if
end function

pure function stdlib_sum_extreme_cancellation_cdp(a) result(s)
    complex(dp), intent(in) :: a(:)
    complex(dp) :: s
    integer :: n, i, j
    complex(dp), allocatable :: temp(:)
    logical, allocatable :: used(:)
    real(dp) :: abs_vals(size(a))
    integer :: sorted_indices(size(a))
    
    n = size(a)
    if (n == 0) then
        s = zero_cdp
        return
    end if
    
    if (n == 1) then
        s = a(1)
        return
    end if
    
    ! For small arrays, try to pair canceling values first
    allocate(temp(n), used(n))
    temp = a
    used = .false.
    
    ! Calculate absolute values for sorting
    do i = 1, n
        abs_vals(i) = abs(a(i))
    end do
    
    ! Simple sorting by magnitude (bubble sort for small arrays)
    do i = 1, n
        sorted_indices(i) = i
    end do
    do i = 1, n-1
        do j = i+1, n
            if (abs_vals(sorted_indices(i)) < abs_vals(sorted_indices(j))) then
                ! Swap indices
                call swap_indices(sorted_indices(i), sorted_indices(j))
            end if
        end do
    end do
    
    ! Try to find canceling pairs among largest values first
    s = zero_cdp
    do i = 1, n-1
        if (used(sorted_indices(i))) cycle
        do j = i+1, n
            if (used(sorted_indices(j))) cycle
            ! Check if values approximately cancel
            if (abs(temp(sorted_indices(i)) + temp(sorted_indices(j))) < &
                abs(temp(sorted_indices(i))) * 1.0e-10_dp) then
                s = s + (temp(sorted_indices(i)) + temp(sorted_indices(j)))
                used(sorted_indices(i)) = .true.
                used(sorted_indices(j)) = .true.
                exit
            end if
        end do
    end do
    
    ! Add remaining unpaired values
    do i = 1, n
        if (.not. used(sorted_indices(i))) then
            s = s + temp(sorted_indices(i))
        end if
    end do
    
contains
    pure subroutine swap_indices(a, b)
        integer, intent(inout) :: a, b
        integer :: temp_idx
        temp_idx = a
        a = b
        b = temp_idx
    end subroutine
end function

pure function detect_extreme_magnitude_ratio_cdp(a) result(is_extreme)
    complex(dp), intent(in) :: a(:)
    logical :: is_extreme
    real(dp) :: max_val, min_val, ratio
    real(dp), parameter :: extreme_threshold = 1.0e10_dp
    real(dp) :: abs_vals(size(a))
    integer :: i
    
    ! Calculate absolute values for both real and complex
    do i = 1, size(a)
        abs_vals(i) = abs(a(i))
    end do
    
    max_val = maxval(abs_vals)
    min_val = minval(abs_vals, mask=abs_vals > tiny(1.0_dp))
    if (min_val <= tiny(1.0_dp)) then
        min_val = tiny(1.0_dp)
    end if
    ratio = max_val / min_val
    is_extreme = ratio > extreme_threshold
end function

pure module function stdlib_sum_adaptive_1d_cdp(a) result(s)
    complex(dp), intent(in) :: a(:)
    complex(dp) :: s
    
    if (detect_extreme_magnitude_ratio_cdp(a)) then
        s = stdlib_sum_pairwise_1d_cdp(a, 1_ilp, size(a, kind=ilp))
    else
        s = stdlib_sum_kahan(a)
    end if
end function

pure module function stdlib_sum_adaptive_1d_cdp_mask(a, mask) result(s)
    complex(dp), intent(in) :: a(:)
    logical, intent(in) :: mask(:)
    complex(dp) :: s
    complex(dp), allocatable :: masked_a(:)
    integer(ilp) :: i, count_true
    
    count_true = count(mask)
    if (count_true == 0) then
        s = zero_cdp
        return
    end if
    
    allocate(masked_a(count_true))
    count_true = 0
    do i = 1, size(a)
        if (mask(i)) then
            count_true = count_true + 1
            masked_a(count_true) = a(i)
        end if
    end do
    
    s = stdlib_sum_adaptive_1d_cdp(masked_a)
end function

!================= N-D implementations ============
pure module function stdlib_sum_2d_int8( x , mask ) result( s )
    integer(int8), intent(in) :: x(:,:)
    logical, intent(in), optional :: mask(:,:)
    integer(int8) :: s
    if(.not.present(mask)) then
        s = sum_recast(x,size(x,kind=ilp))
    else
        s = sum_recast_mask(x,mask,size(x,kind=ilp))
    end if
contains
    pure integer(int8) function sum_recast(b,n)
      integer(ilp), intent(in) :: n
      integer(int8), intent(in) :: b(n)
      sum_recast = stdlib_sum(b)
    end function
    pure integer(int8) function sum_recast_mask(b,m,n)
      integer(ilp), intent(in) :: n
      integer(int8), intent(in) :: b(n)
      logical, intent(in) :: m(n)
      sum_recast_mask = stdlib_sum(b,m)
    end function
end function

pure module function stdlib_sum_2d_dim_int8( x , dim, mask ) result( s )
    integer(int8), intent(in) :: x(:,:)
    integer, intent(in):: dim
    logical, intent(in), optional :: mask(:,:)
    integer(int8) :: s(merge(size(x, 1), size(x, 2), mask=1<dim))
    integer :: j 

    if(.not.present(mask)) then
        if(dim<2)then
            do j = 1, size(x,dim=2)
                s(j) = stdlib_sum( x(:, j) )
            end do
        else
            do j = 1, size(x,dim=1)
                s(j) = stdlib_sum( x(j, :) )
            end do
        end if
    else 
        if(dim<2)then
            do j = 1, size(x,dim=2)
                s(j) = stdlib_sum( x(:, j), mask=mask(:, j) )
            end do
        else
            do j = 1, size(x,dim=1)
                s(j) = stdlib_sum( x(j, :), mask=mask(j, :) )
            end do
        end if
    end if

end function
pure module function stdlib_sum_3d_int8( x , mask ) result( s )
    integer(int8), intent(in) :: x(:,:,:)
    logical, intent(in), optional :: mask(:,:,:)
    integer(int8) :: s
    if(.not.present(mask)) then
        s = sum_recast(x,size(x,kind=ilp))
    else
        s = sum_recast_mask(x,mask,size(x,kind=ilp))
    end if
contains
    pure integer(int8) function sum_recast(b,n)
      integer(ilp), intent(in) :: n
      integer(int8), intent(in) :: b(n)
      sum_recast = stdlib_sum(b)
    end function
    pure integer(int8) function sum_recast_mask(b,m,n)
      integer(ilp), intent(in) :: n
      integer(int8), intent(in) :: b(n)
      logical, intent(in) :: m(n)
      sum_recast_mask = stdlib_sum(b,m)
    end function
end function

pure module function stdlib_sum_3d_dim_int8( x , dim, mask ) result( s )
    integer(int8), intent(in) :: x(:,:,:)
    integer, intent(in):: dim
    logical, intent(in), optional :: mask(:,:,:)
    integer(int8) :: s(merge(size(x, 1), size(x, 2), mask=1<dim), merge(size(x, 2), size(x, 3), mask=2<dim))
    integer :: j 

    if(.not.present(mask)) then
        if(dim<3)then
            do j = 1, size(x,dim=3)
                s(:, j) = stdlib_sum( x(:, :, j), dim=dim )
            end do
        else
            do j = 1, size(x,dim=1)
                s(j, :) = stdlib_sum( x(j, :, :), dim=2 )
            end do
        end if
    else 
        if(dim<3)then
            do j = 1, size(x,dim=3)
                s(:, j) = stdlib_sum( x(:, :, j), dim=dim, mask=mask(:, :, j) )
            end do
        else
            do j = 1, size(x,dim=1)
                s(j, :) = stdlib_sum( x(j, :, :), dim=2, mask=mask(j, :, :) )
            end do
        end if
    end if

end function
pure module function stdlib_sum_4d_int8( x , mask ) result( s )
    integer(int8), intent(in) :: x(:,:,:,:)
    logical, intent(in), optional :: mask(:,:,:,:)
    integer(int8) :: s
    if(.not.present(mask)) then
        s = sum_recast(x,size(x,kind=ilp))
    else
        s = sum_recast_mask(x,mask,size(x,kind=ilp))
    end if
contains
    pure integer(int8) function sum_recast(b,n)
      integer(ilp), intent(in) :: n
      integer(int8), intent(in) :: b(n)
      sum_recast = stdlib_sum(b)
    end function
    pure integer(int8) function sum_recast_mask(b,m,n)
      integer(ilp), intent(in) :: n
      integer(int8), intent(in) :: b(n)
      logical, intent(in) :: m(n)
      sum_recast_mask = stdlib_sum(b,m)
    end function
end function

pure module function stdlib_sum_4d_dim_int8( x , dim, mask ) result( s )
    integer(int8), intent(in) :: x(:,:,:,:)
    integer, intent(in):: dim
    logical, intent(in), optional :: mask(:,:,:,:)
    integer(int8) :: s(merge(size(x, 1), size(x, 2), mask=1<dim), merge(size(x, 2), size(x, 3), mask=2<dim), merge(size(x, 3),&
        & size(x, 4), mask=3<dim))
    integer :: j 

    if(.not.present(mask)) then
        if(dim<4)then
            do j = 1, size(x,dim=4)
                s(:, :, j) = stdlib_sum( x(:, :, :, j), dim=dim )
            end do
        else
            do j = 1, size(x,dim=1)
                s(j, :, :) = stdlib_sum( x(j, :, :, :), dim=3 )
            end do
        end if
    else 
        if(dim<4)then
            do j = 1, size(x,dim=4)
                s(:, :, j) = stdlib_sum( x(:, :, :, j), dim=dim, mask=mask(:, :, :, j) )
            end do
        else
            do j = 1, size(x,dim=1)
                s(j, :, :) = stdlib_sum( x(j, :, :, :), dim=3, mask=mask(j, :, :, :) )
            end do
        end if
    end if

end function
pure module function stdlib_sum_2d_int16( x , mask ) result( s )
    integer(int16), intent(in) :: x(:,:)
    logical, intent(in), optional :: mask(:,:)
    integer(int16) :: s
    if(.not.present(mask)) then
        s = sum_recast(x,size(x,kind=ilp))
    else
        s = sum_recast_mask(x,mask,size(x,kind=ilp))
    end if
contains
    pure integer(int16) function sum_recast(b,n)
      integer(ilp), intent(in) :: n
      integer(int16), intent(in) :: b(n)
      sum_recast = stdlib_sum(b)
    end function
    pure integer(int16) function sum_recast_mask(b,m,n)
      integer(ilp), intent(in) :: n
      integer(int16), intent(in) :: b(n)
      logical, intent(in) :: m(n)
      sum_recast_mask = stdlib_sum(b,m)
    end function
end function

pure module function stdlib_sum_2d_dim_int16( x , dim, mask ) result( s )
    integer(int16), intent(in) :: x(:,:)
    integer, intent(in):: dim
    logical, intent(in), optional :: mask(:,:)
    integer(int16) :: s(merge(size(x, 1), size(x, 2), mask=1<dim))
    integer :: j 

    if(.not.present(mask)) then
        if(dim<2)then
            do j = 1, size(x,dim=2)
                s(j) = stdlib_sum( x(:, j) )
            end do
        else
            do j = 1, size(x,dim=1)
                s(j) = stdlib_sum( x(j, :) )
            end do
        end if
    else 
        if(dim<2)then
            do j = 1, size(x,dim=2)
                s(j) = stdlib_sum( x(:, j), mask=mask(:, j) )
            end do
        else
            do j = 1, size(x,dim=1)
                s(j) = stdlib_sum( x(j, :), mask=mask(j, :) )
            end do
        end if
    end if

end function
pure module function stdlib_sum_3d_int16( x , mask ) result( s )
    integer(int16), intent(in) :: x(:,:,:)
    logical, intent(in), optional :: mask(:,:,:)
    integer(int16) :: s
    if(.not.present(mask)) then
        s = sum_recast(x,size(x,kind=ilp))
    else
        s = sum_recast_mask(x,mask,size(x,kind=ilp))
    end if
contains
    pure integer(int16) function sum_recast(b,n)
      integer(ilp), intent(in) :: n
      integer(int16), intent(in) :: b(n)
      sum_recast = stdlib_sum(b)
    end function
    pure integer(int16) function sum_recast_mask(b,m,n)
      integer(ilp), intent(in) :: n
      integer(int16), intent(in) :: b(n)
      logical, intent(in) :: m(n)
      sum_recast_mask = stdlib_sum(b,m)
    end function
end function

pure module function stdlib_sum_3d_dim_int16( x , dim, mask ) result( s )
    integer(int16), intent(in) :: x(:,:,:)
    integer, intent(in):: dim
    logical, intent(in), optional :: mask(:,:,:)
    integer(int16) :: s(merge(size(x, 1), size(x, 2), mask=1<dim), merge(size(x, 2), size(x, 3), mask=2<dim))
    integer :: j 

    if(.not.present(mask)) then
        if(dim<3)then
            do j = 1, size(x,dim=3)
                s(:, j) = stdlib_sum( x(:, :, j), dim=dim )
            end do
        else
            do j = 1, size(x,dim=1)
                s(j, :) = stdlib_sum( x(j, :, :), dim=2 )
            end do
        end if
    else 
        if(dim<3)then
            do j = 1, size(x,dim=3)
                s(:, j) = stdlib_sum( x(:, :, j), dim=dim, mask=mask(:, :, j) )
            end do
        else
            do j = 1, size(x,dim=1)
                s(j, :) = stdlib_sum( x(j, :, :), dim=2, mask=mask(j, :, :) )
            end do
        end if
    end if

end function
pure module function stdlib_sum_4d_int16( x , mask ) result( s )
    integer(int16), intent(in) :: x(:,:,:,:)
    logical, intent(in), optional :: mask(:,:,:,:)
    integer(int16) :: s
    if(.not.present(mask)) then
        s = sum_recast(x,size(x,kind=ilp))
    else
        s = sum_recast_mask(x,mask,size(x,kind=ilp))
    end if
contains
    pure integer(int16) function sum_recast(b,n)
      integer(ilp), intent(in) :: n
      integer(int16), intent(in) :: b(n)
      sum_recast = stdlib_sum(b)
    end function
    pure integer(int16) function sum_recast_mask(b,m,n)
      integer(ilp), intent(in) :: n
      integer(int16), intent(in) :: b(n)
      logical, intent(in) :: m(n)
      sum_recast_mask = stdlib_sum(b,m)
    end function
end function

pure module function stdlib_sum_4d_dim_int16( x , dim, mask ) result( s )
    integer(int16), intent(in) :: x(:,:,:,:)
    integer, intent(in):: dim
    logical, intent(in), optional :: mask(:,:,:,:)
    integer(int16) :: s(merge(size(x, 1), size(x, 2), mask=1<dim), merge(size(x, 2), size(x, 3), mask=2<dim), merge(size(x, 3),&
        & size(x, 4), mask=3<dim))
    integer :: j 

    if(.not.present(mask)) then
        if(dim<4)then
            do j = 1, size(x,dim=4)
                s(:, :, j) = stdlib_sum( x(:, :, :, j), dim=dim )
            end do
        else
            do j = 1, size(x,dim=1)
                s(j, :, :) = stdlib_sum( x(j, :, :, :), dim=3 )
            end do
        end if
    else 
        if(dim<4)then
            do j = 1, size(x,dim=4)
                s(:, :, j) = stdlib_sum( x(:, :, :, j), dim=dim, mask=mask(:, :, :, j) )
            end do
        else
            do j = 1, size(x,dim=1)
                s(j, :, :) = stdlib_sum( x(j, :, :, :), dim=3, mask=mask(j, :, :, :) )
            end do
        end if
    end if

end function
pure module function stdlib_sum_2d_int32( x , mask ) result( s )
    integer(int32), intent(in) :: x(:,:)
    logical, intent(in), optional :: mask(:,:)
    integer(int32) :: s
    if(.not.present(mask)) then
        s = sum_recast(x,size(x,kind=ilp))
    else
        s = sum_recast_mask(x,mask,size(x,kind=ilp))
    end if
contains
    pure integer(int32) function sum_recast(b,n)
      integer(ilp), intent(in) :: n
      integer(int32), intent(in) :: b(n)
      sum_recast = stdlib_sum(b)
    end function
    pure integer(int32) function sum_recast_mask(b,m,n)
      integer(ilp), intent(in) :: n
      integer(int32), intent(in) :: b(n)
      logical, intent(in) :: m(n)
      sum_recast_mask = stdlib_sum(b,m)
    end function
end function

pure module function stdlib_sum_2d_dim_int32( x , dim, mask ) result( s )
    integer(int32), intent(in) :: x(:,:)
    integer, intent(in):: dim
    logical, intent(in), optional :: mask(:,:)
    integer(int32) :: s(merge(size(x, 1), size(x, 2), mask=1<dim))
    integer :: j 

    if(.not.present(mask)) then
        if(dim<2)then
            do j = 1, size(x,dim=2)
                s(j) = stdlib_sum( x(:, j) )
            end do
        else
            do j = 1, size(x,dim=1)
                s(j) = stdlib_sum( x(j, :) )
            end do
        end if
    else 
        if(dim<2)then
            do j = 1, size(x,dim=2)
                s(j) = stdlib_sum( x(:, j), mask=mask(:, j) )
            end do
        else
            do j = 1, size(x,dim=1)
                s(j) = stdlib_sum( x(j, :), mask=mask(j, :) )
            end do
        end if
    end if

end function
pure module function stdlib_sum_3d_int32( x , mask ) result( s )
    integer(int32), intent(in) :: x(:,:,:)
    logical, intent(in), optional :: mask(:,:,:)
    integer(int32) :: s
    if(.not.present(mask)) then
        s = sum_recast(x,size(x,kind=ilp))
    else
        s = sum_recast_mask(x,mask,size(x,kind=ilp))
    end if
contains
    pure integer(int32) function sum_recast(b,n)
      integer(ilp), intent(in) :: n
      integer(int32), intent(in) :: b(n)
      sum_recast = stdlib_sum(b)
    end function
    pure integer(int32) function sum_recast_mask(b,m,n)
      integer(ilp), intent(in) :: n
      integer(int32), intent(in) :: b(n)
      logical, intent(in) :: m(n)
      sum_recast_mask = stdlib_sum(b,m)
    end function
end function

pure module function stdlib_sum_3d_dim_int32( x , dim, mask ) result( s )
    integer(int32), intent(in) :: x(:,:,:)
    integer, intent(in):: dim
    logical, intent(in), optional :: mask(:,:,:)
    integer(int32) :: s(merge(size(x, 1), size(x, 2), mask=1<dim), merge(size(x, 2), size(x, 3), mask=2<dim))
    integer :: j 

    if(.not.present(mask)) then
        if(dim<3)then
            do j = 1, size(x,dim=3)
                s(:, j) = stdlib_sum( x(:, :, j), dim=dim )
            end do
        else
            do j = 1, size(x,dim=1)
                s(j, :) = stdlib_sum( x(j, :, :), dim=2 )
            end do
        end if
    else 
        if(dim<3)then
            do j = 1, size(x,dim=3)
                s(:, j) = stdlib_sum( x(:, :, j), dim=dim, mask=mask(:, :, j) )
            end do
        else
            do j = 1, size(x,dim=1)
                s(j, :) = stdlib_sum( x(j, :, :), dim=2, mask=mask(j, :, :) )
            end do
        end if
    end if

end function
pure module function stdlib_sum_4d_int32( x , mask ) result( s )
    integer(int32), intent(in) :: x(:,:,:,:)
    logical, intent(in), optional :: mask(:,:,:,:)
    integer(int32) :: s
    if(.not.present(mask)) then
        s = sum_recast(x,size(x,kind=ilp))
    else
        s = sum_recast_mask(x,mask,size(x,kind=ilp))
    end if
contains
    pure integer(int32) function sum_recast(b,n)
      integer(ilp), intent(in) :: n
      integer(int32), intent(in) :: b(n)
      sum_recast = stdlib_sum(b)
    end function
    pure integer(int32) function sum_recast_mask(b,m,n)
      integer(ilp), intent(in) :: n
      integer(int32), intent(in) :: b(n)
      logical, intent(in) :: m(n)
      sum_recast_mask = stdlib_sum(b,m)
    end function
end function

pure module function stdlib_sum_4d_dim_int32( x , dim, mask ) result( s )
    integer(int32), intent(in) :: x(:,:,:,:)
    integer, intent(in):: dim
    logical, intent(in), optional :: mask(:,:,:,:)
    integer(int32) :: s(merge(size(x, 1), size(x, 2), mask=1<dim), merge(size(x, 2), size(x, 3), mask=2<dim), merge(size(x, 3),&
        & size(x, 4), mask=3<dim))
    integer :: j 

    if(.not.present(mask)) then
        if(dim<4)then
            do j = 1, size(x,dim=4)
                s(:, :, j) = stdlib_sum( x(:, :, :, j), dim=dim )
            end do
        else
            do j = 1, size(x,dim=1)
                s(j, :, :) = stdlib_sum( x(j, :, :, :), dim=3 )
            end do
        end if
    else 
        if(dim<4)then
            do j = 1, size(x,dim=4)
                s(:, :, j) = stdlib_sum( x(:, :, :, j), dim=dim, mask=mask(:, :, :, j) )
            end do
        else
            do j = 1, size(x,dim=1)
                s(j, :, :) = stdlib_sum( x(j, :, :, :), dim=3, mask=mask(j, :, :, :) )
            end do
        end if
    end if

end function
pure module function stdlib_sum_2d_int64( x , mask ) result( s )
    integer(int64), intent(in) :: x(:,:)
    logical, intent(in), optional :: mask(:,:)
    integer(int64) :: s
    if(.not.present(mask)) then
        s = sum_recast(x,size(x,kind=ilp))
    else
        s = sum_recast_mask(x,mask,size(x,kind=ilp))
    end if
contains
    pure integer(int64) function sum_recast(b,n)
      integer(ilp), intent(in) :: n
      integer(int64), intent(in) :: b(n)
      sum_recast = stdlib_sum(b)
    end function
    pure integer(int64) function sum_recast_mask(b,m,n)
      integer(ilp), intent(in) :: n
      integer(int64), intent(in) :: b(n)
      logical, intent(in) :: m(n)
      sum_recast_mask = stdlib_sum(b,m)
    end function
end function

pure module function stdlib_sum_2d_dim_int64( x , dim, mask ) result( s )
    integer(int64), intent(in) :: x(:,:)
    integer, intent(in):: dim
    logical, intent(in), optional :: mask(:,:)
    integer(int64) :: s(merge(size(x, 1), size(x, 2), mask=1<dim))
    integer :: j 

    if(.not.present(mask)) then
        if(dim<2)then
            do j = 1, size(x,dim=2)
                s(j) = stdlib_sum( x(:, j) )
            end do
        else
            do j = 1, size(x,dim=1)
                s(j) = stdlib_sum( x(j, :) )
            end do
        end if
    else 
        if(dim<2)then
            do j = 1, size(x,dim=2)
                s(j) = stdlib_sum( x(:, j), mask=mask(:, j) )
            end do
        else
            do j = 1, size(x,dim=1)
                s(j) = stdlib_sum( x(j, :), mask=mask(j, :) )
            end do
        end if
    end if

end function
pure module function stdlib_sum_3d_int64( x , mask ) result( s )
    integer(int64), intent(in) :: x(:,:,:)
    logical, intent(in), optional :: mask(:,:,:)
    integer(int64) :: s
    if(.not.present(mask)) then
        s = sum_recast(x,size(x,kind=ilp))
    else
        s = sum_recast_mask(x,mask,size(x,kind=ilp))
    end if
contains
    pure integer(int64) function sum_recast(b,n)
      integer(ilp), intent(in) :: n
      integer(int64), intent(in) :: b(n)
      sum_recast = stdlib_sum(b)
    end function
    pure integer(int64) function sum_recast_mask(b,m,n)
      integer(ilp), intent(in) :: n
      integer(int64), intent(in) :: b(n)
      logical, intent(in) :: m(n)
      sum_recast_mask = stdlib_sum(b,m)
    end function
end function

pure module function stdlib_sum_3d_dim_int64( x , dim, mask ) result( s )
    integer(int64), intent(in) :: x(:,:,:)
    integer, intent(in):: dim
    logical, intent(in), optional :: mask(:,:,:)
    integer(int64) :: s(merge(size(x, 1), size(x, 2), mask=1<dim), merge(size(x, 2), size(x, 3), mask=2<dim))
    integer :: j 

    if(.not.present(mask)) then
        if(dim<3)then
            do j = 1, size(x,dim=3)
                s(:, j) = stdlib_sum( x(:, :, j), dim=dim )
            end do
        else
            do j = 1, size(x,dim=1)
                s(j, :) = stdlib_sum( x(j, :, :), dim=2 )
            end do
        end if
    else 
        if(dim<3)then
            do j = 1, size(x,dim=3)
                s(:, j) = stdlib_sum( x(:, :, j), dim=dim, mask=mask(:, :, j) )
            end do
        else
            do j = 1, size(x,dim=1)
                s(j, :) = stdlib_sum( x(j, :, :), dim=2, mask=mask(j, :, :) )
            end do
        end if
    end if

end function
pure module function stdlib_sum_4d_int64( x , mask ) result( s )
    integer(int64), intent(in) :: x(:,:,:,:)
    logical, intent(in), optional :: mask(:,:,:,:)
    integer(int64) :: s
    if(.not.present(mask)) then
        s = sum_recast(x,size(x,kind=ilp))
    else
        s = sum_recast_mask(x,mask,size(x,kind=ilp))
    end if
contains
    pure integer(int64) function sum_recast(b,n)
      integer(ilp), intent(in) :: n
      integer(int64), intent(in) :: b(n)
      sum_recast = stdlib_sum(b)
    end function
    pure integer(int64) function sum_recast_mask(b,m,n)
      integer(ilp), intent(in) :: n
      integer(int64), intent(in) :: b(n)
      logical, intent(in) :: m(n)
      sum_recast_mask = stdlib_sum(b,m)
    end function
end function

pure module function stdlib_sum_4d_dim_int64( x , dim, mask ) result( s )
    integer(int64), intent(in) :: x(:,:,:,:)
    integer, intent(in):: dim
    logical, intent(in), optional :: mask(:,:,:,:)
    integer(int64) :: s(merge(size(x, 1), size(x, 2), mask=1<dim), merge(size(x, 2), size(x, 3), mask=2<dim), merge(size(x, 3),&
        & size(x, 4), mask=3<dim))
    integer :: j 

    if(.not.present(mask)) then
        if(dim<4)then
            do j = 1, size(x,dim=4)
                s(:, :, j) = stdlib_sum( x(:, :, :, j), dim=dim )
            end do
        else
            do j = 1, size(x,dim=1)
                s(j, :, :) = stdlib_sum( x(j, :, :, :), dim=3 )
            end do
        end if
    else 
        if(dim<4)then
            do j = 1, size(x,dim=4)
                s(:, :, j) = stdlib_sum( x(:, :, :, j), dim=dim, mask=mask(:, :, :, j) )
            end do
        else
            do j = 1, size(x,dim=1)
                s(j, :, :) = stdlib_sum( x(j, :, :, :), dim=3, mask=mask(j, :, :, :) )
            end do
        end if
    end if

end function
pure module function stdlib_sum_2d_sp( x , mask ) result( s )
    real(sp), intent(in) :: x(:,:)
    logical, intent(in), optional :: mask(:,:)
    real(sp) :: s
    if(.not.present(mask)) then
        s = sum_recast(x,size(x,kind=ilp))
    else
        s = sum_recast_mask(x,mask,size(x,kind=ilp))
    end if
contains
    pure real(sp) function sum_recast(b,n)
      integer(ilp), intent(in) :: n
      real(sp), intent(in) :: b(n)
      sum_recast = stdlib_sum(b)
    end function
    pure real(sp) function sum_recast_mask(b,m,n)
      integer(ilp), intent(in) :: n
      real(sp), intent(in) :: b(n)
      logical, intent(in) :: m(n)
      sum_recast_mask = stdlib_sum(b,m)
    end function
end function

pure module function stdlib_sum_2d_dim_sp( x , dim, mask ) result( s )
    real(sp), intent(in) :: x(:,:)
    integer, intent(in):: dim
    logical, intent(in), optional :: mask(:,:)
    real(sp) :: s(merge(size(x, 1), size(x, 2), mask=1<dim))
    integer :: j 

    if(.not.present(mask)) then
        if(dim<2)then
            do j = 1, size(x,dim=2)
                s(j) = stdlib_sum( x(:, j) )
            end do
        else
            do j = 1, size(x,dim=1)
                s(j) = stdlib_sum( x(j, :) )
            end do
        end if
    else 
        if(dim<2)then
            do j = 1, size(x,dim=2)
                s(j) = stdlib_sum( x(:, j), mask=mask(:, j) )
            end do
        else
            do j = 1, size(x,dim=1)
                s(j) = stdlib_sum( x(j, :), mask=mask(j, :) )
            end do
        end if
    end if

end function
pure module function stdlib_sum_3d_sp( x , mask ) result( s )
    real(sp), intent(in) :: x(:,:,:)
    logical, intent(in), optional :: mask(:,:,:)
    real(sp) :: s
    if(.not.present(mask)) then
        s = sum_recast(x,size(x,kind=ilp))
    else
        s = sum_recast_mask(x,mask,size(x,kind=ilp))
    end if
contains
    pure real(sp) function sum_recast(b,n)
      integer(ilp), intent(in) :: n
      real(sp), intent(in) :: b(n)
      sum_recast = stdlib_sum(b)
    end function
    pure real(sp) function sum_recast_mask(b,m,n)
      integer(ilp), intent(in) :: n
      real(sp), intent(in) :: b(n)
      logical, intent(in) :: m(n)
      sum_recast_mask = stdlib_sum(b,m)
    end function
end function

pure module function stdlib_sum_3d_dim_sp( x , dim, mask ) result( s )
    real(sp), intent(in) :: x(:,:,:)
    integer, intent(in):: dim
    logical, intent(in), optional :: mask(:,:,:)
    real(sp) :: s(merge(size(x, 1), size(x, 2), mask=1<dim), merge(size(x, 2), size(x, 3), mask=2<dim))
    integer :: j 

    if(.not.present(mask)) then
        if(dim<3)then
            do j = 1, size(x,dim=3)
                s(:, j) = stdlib_sum( x(:, :, j), dim=dim )
            end do
        else
            do j = 1, size(x,dim=1)
                s(j, :) = stdlib_sum( x(j, :, :), dim=2 )
            end do
        end if
    else 
        if(dim<3)then
            do j = 1, size(x,dim=3)
                s(:, j) = stdlib_sum( x(:, :, j), dim=dim, mask=mask(:, :, j) )
            end do
        else
            do j = 1, size(x,dim=1)
                s(j, :) = stdlib_sum( x(j, :, :), dim=2, mask=mask(j, :, :) )
            end do
        end if
    end if

end function
pure module function stdlib_sum_4d_sp( x , mask ) result( s )
    real(sp), intent(in) :: x(:,:,:,:)
    logical, intent(in), optional :: mask(:,:,:,:)
    real(sp) :: s
    if(.not.present(mask)) then
        s = sum_recast(x,size(x,kind=ilp))
    else
        s = sum_recast_mask(x,mask,size(x,kind=ilp))
    end if
contains
    pure real(sp) function sum_recast(b,n)
      integer(ilp), intent(in) :: n
      real(sp), intent(in) :: b(n)
      sum_recast = stdlib_sum(b)
    end function
    pure real(sp) function sum_recast_mask(b,m,n)
      integer(ilp), intent(in) :: n
      real(sp), intent(in) :: b(n)
      logical, intent(in) :: m(n)
      sum_recast_mask = stdlib_sum(b,m)
    end function
end function

pure module function stdlib_sum_4d_dim_sp( x , dim, mask ) result( s )
    real(sp), intent(in) :: x(:,:,:,:)
    integer, intent(in):: dim
    logical, intent(in), optional :: mask(:,:,:,:)
    real(sp) :: s(merge(size(x, 1), size(x, 2), mask=1<dim), merge(size(x, 2), size(x, 3), mask=2<dim), merge(size(x, 3), size(x,&
        & 4), mask=3<dim))
    integer :: j 

    if(.not.present(mask)) then
        if(dim<4)then
            do j = 1, size(x,dim=4)
                s(:, :, j) = stdlib_sum( x(:, :, :, j), dim=dim )
            end do
        else
            do j = 1, size(x,dim=1)
                s(j, :, :) = stdlib_sum( x(j, :, :, :), dim=3 )
            end do
        end if
    else 
        if(dim<4)then
            do j = 1, size(x,dim=4)
                s(:, :, j) = stdlib_sum( x(:, :, :, j), dim=dim, mask=mask(:, :, :, j) )
            end do
        else
            do j = 1, size(x,dim=1)
                s(j, :, :) = stdlib_sum( x(j, :, :, :), dim=3, mask=mask(j, :, :, :) )
            end do
        end if
    end if

end function
pure module function stdlib_sum_2d_dp( x , mask ) result( s )
    real(dp), intent(in) :: x(:,:)
    logical, intent(in), optional :: mask(:,:)
    real(dp) :: s
    if(.not.present(mask)) then
        s = sum_recast(x,size(x,kind=ilp))
    else
        s = sum_recast_mask(x,mask,size(x,kind=ilp))
    end if
contains
    pure real(dp) function sum_recast(b,n)
      integer(ilp), intent(in) :: n
      real(dp), intent(in) :: b(n)
      sum_recast = stdlib_sum(b)
    end function
    pure real(dp) function sum_recast_mask(b,m,n)
      integer(ilp), intent(in) :: n
      real(dp), intent(in) :: b(n)
      logical, intent(in) :: m(n)
      sum_recast_mask = stdlib_sum(b,m)
    end function
end function

pure module function stdlib_sum_2d_dim_dp( x , dim, mask ) result( s )
    real(dp), intent(in) :: x(:,:)
    integer, intent(in):: dim
    logical, intent(in), optional :: mask(:,:)
    real(dp) :: s(merge(size(x, 1), size(x, 2), mask=1<dim))
    integer :: j 

    if(.not.present(mask)) then
        if(dim<2)then
            do j = 1, size(x,dim=2)
                s(j) = stdlib_sum( x(:, j) )
            end do
        else
            do j = 1, size(x,dim=1)
                s(j) = stdlib_sum( x(j, :) )
            end do
        end if
    else 
        if(dim<2)then
            do j = 1, size(x,dim=2)
                s(j) = stdlib_sum( x(:, j), mask=mask(:, j) )
            end do
        else
            do j = 1, size(x,dim=1)
                s(j) = stdlib_sum( x(j, :), mask=mask(j, :) )
            end do
        end if
    end if

end function
pure module function stdlib_sum_3d_dp( x , mask ) result( s )
    real(dp), intent(in) :: x(:,:,:)
    logical, intent(in), optional :: mask(:,:,:)
    real(dp) :: s
    if(.not.present(mask)) then
        s = sum_recast(x,size(x,kind=ilp))
    else
        s = sum_recast_mask(x,mask,size(x,kind=ilp))
    end if
contains
    pure real(dp) function sum_recast(b,n)
      integer(ilp), intent(in) :: n
      real(dp), intent(in) :: b(n)
      sum_recast = stdlib_sum(b)
    end function
    pure real(dp) function sum_recast_mask(b,m,n)
      integer(ilp), intent(in) :: n
      real(dp), intent(in) :: b(n)
      logical, intent(in) :: m(n)
      sum_recast_mask = stdlib_sum(b,m)
    end function
end function

pure module function stdlib_sum_3d_dim_dp( x , dim, mask ) result( s )
    real(dp), intent(in) :: x(:,:,:)
    integer, intent(in):: dim
    logical, intent(in), optional :: mask(:,:,:)
    real(dp) :: s(merge(size(x, 1), size(x, 2), mask=1<dim), merge(size(x, 2), size(x, 3), mask=2<dim))
    integer :: j 

    if(.not.present(mask)) then
        if(dim<3)then
            do j = 1, size(x,dim=3)
                s(:, j) = stdlib_sum( x(:, :, j), dim=dim )
            end do
        else
            do j = 1, size(x,dim=1)
                s(j, :) = stdlib_sum( x(j, :, :), dim=2 )
            end do
        end if
    else 
        if(dim<3)then
            do j = 1, size(x,dim=3)
                s(:, j) = stdlib_sum( x(:, :, j), dim=dim, mask=mask(:, :, j) )
            end do
        else
            do j = 1, size(x,dim=1)
                s(j, :) = stdlib_sum( x(j, :, :), dim=2, mask=mask(j, :, :) )
            end do
        end if
    end if

end function
pure module function stdlib_sum_4d_dp( x , mask ) result( s )
    real(dp), intent(in) :: x(:,:,:,:)
    logical, intent(in), optional :: mask(:,:,:,:)
    real(dp) :: s
    if(.not.present(mask)) then
        s = sum_recast(x,size(x,kind=ilp))
    else
        s = sum_recast_mask(x,mask,size(x,kind=ilp))
    end if
contains
    pure real(dp) function sum_recast(b,n)
      integer(ilp), intent(in) :: n
      real(dp), intent(in) :: b(n)
      sum_recast = stdlib_sum(b)
    end function
    pure real(dp) function sum_recast_mask(b,m,n)
      integer(ilp), intent(in) :: n
      real(dp), intent(in) :: b(n)
      logical, intent(in) :: m(n)
      sum_recast_mask = stdlib_sum(b,m)
    end function
end function

pure module function stdlib_sum_4d_dim_dp( x , dim, mask ) result( s )
    real(dp), intent(in) :: x(:,:,:,:)
    integer, intent(in):: dim
    logical, intent(in), optional :: mask(:,:,:,:)
    real(dp) :: s(merge(size(x, 1), size(x, 2), mask=1<dim), merge(size(x, 2), size(x, 3), mask=2<dim), merge(size(x, 3), size(x,&
        & 4), mask=3<dim))
    integer :: j 

    if(.not.present(mask)) then
        if(dim<4)then
            do j = 1, size(x,dim=4)
                s(:, :, j) = stdlib_sum( x(:, :, :, j), dim=dim )
            end do
        else
            do j = 1, size(x,dim=1)
                s(j, :, :) = stdlib_sum( x(j, :, :, :), dim=3 )
            end do
        end if
    else 
        if(dim<4)then
            do j = 1, size(x,dim=4)
                s(:, :, j) = stdlib_sum( x(:, :, :, j), dim=dim, mask=mask(:, :, :, j) )
            end do
        else
            do j = 1, size(x,dim=1)
                s(j, :, :) = stdlib_sum( x(j, :, :, :), dim=3, mask=mask(j, :, :, :) )
            end do
        end if
    end if

end function
pure module function stdlib_sum_2d_csp( x , mask ) result( s )
    complex(sp), intent(in) :: x(:,:)
    logical, intent(in), optional :: mask(:,:)
    complex(sp) :: s
    if(.not.present(mask)) then
        s = sum_recast(x,size(x,kind=ilp))
    else
        s = sum_recast_mask(x,mask,size(x,kind=ilp))
    end if
contains
    pure complex(sp) function sum_recast(b,n)
      integer(ilp), intent(in) :: n
      complex(sp), intent(in) :: b(n)
      sum_recast = stdlib_sum(b)
    end function
    pure complex(sp) function sum_recast_mask(b,m,n)
      integer(ilp), intent(in) :: n
      complex(sp), intent(in) :: b(n)
      logical, intent(in) :: m(n)
      sum_recast_mask = stdlib_sum(b,m)
    end function
end function

pure module function stdlib_sum_2d_dim_csp( x , dim, mask ) result( s )
    complex(sp), intent(in) :: x(:,:)
    integer, intent(in):: dim
    logical, intent(in), optional :: mask(:,:)
    complex(sp) :: s(merge(size(x, 1), size(x, 2), mask=1<dim))
    integer :: j 

    if(.not.present(mask)) then
        if(dim<2)then
            do j = 1, size(x,dim=2)
                s(j) = stdlib_sum( x(:, j) )
            end do
        else
            do j = 1, size(x,dim=1)
                s(j) = stdlib_sum( x(j, :) )
            end do
        end if
    else 
        if(dim<2)then
            do j = 1, size(x,dim=2)
                s(j) = stdlib_sum( x(:, j), mask=mask(:, j) )
            end do
        else
            do j = 1, size(x,dim=1)
                s(j) = stdlib_sum( x(j, :), mask=mask(j, :) )
            end do
        end if
    end if

end function
pure module function stdlib_sum_3d_csp( x , mask ) result( s )
    complex(sp), intent(in) :: x(:,:,:)
    logical, intent(in), optional :: mask(:,:,:)
    complex(sp) :: s
    if(.not.present(mask)) then
        s = sum_recast(x,size(x,kind=ilp))
    else
        s = sum_recast_mask(x,mask,size(x,kind=ilp))
    end if
contains
    pure complex(sp) function sum_recast(b,n)
      integer(ilp), intent(in) :: n
      complex(sp), intent(in) :: b(n)
      sum_recast = stdlib_sum(b)
    end function
    pure complex(sp) function sum_recast_mask(b,m,n)
      integer(ilp), intent(in) :: n
      complex(sp), intent(in) :: b(n)
      logical, intent(in) :: m(n)
      sum_recast_mask = stdlib_sum(b,m)
    end function
end function

pure module function stdlib_sum_3d_dim_csp( x , dim, mask ) result( s )
    complex(sp), intent(in) :: x(:,:,:)
    integer, intent(in):: dim
    logical, intent(in), optional :: mask(:,:,:)
    complex(sp) :: s(merge(size(x, 1), size(x, 2), mask=1<dim), merge(size(x, 2), size(x, 3), mask=2<dim))
    integer :: j 

    if(.not.present(mask)) then
        if(dim<3)then
            do j = 1, size(x,dim=3)
                s(:, j) = stdlib_sum( x(:, :, j), dim=dim )
            end do
        else
            do j = 1, size(x,dim=1)
                s(j, :) = stdlib_sum( x(j, :, :), dim=2 )
            end do
        end if
    else 
        if(dim<3)then
            do j = 1, size(x,dim=3)
                s(:, j) = stdlib_sum( x(:, :, j), dim=dim, mask=mask(:, :, j) )
            end do
        else
            do j = 1, size(x,dim=1)
                s(j, :) = stdlib_sum( x(j, :, :), dim=2, mask=mask(j, :, :) )
            end do
        end if
    end if

end function
pure module function stdlib_sum_4d_csp( x , mask ) result( s )
    complex(sp), intent(in) :: x(:,:,:,:)
    logical, intent(in), optional :: mask(:,:,:,:)
    complex(sp) :: s
    if(.not.present(mask)) then
        s = sum_recast(x,size(x,kind=ilp))
    else
        s = sum_recast_mask(x,mask,size(x,kind=ilp))
    end if
contains
    pure complex(sp) function sum_recast(b,n)
      integer(ilp), intent(in) :: n
      complex(sp), intent(in) :: b(n)
      sum_recast = stdlib_sum(b)
    end function
    pure complex(sp) function sum_recast_mask(b,m,n)
      integer(ilp), intent(in) :: n
      complex(sp), intent(in) :: b(n)
      logical, intent(in) :: m(n)
      sum_recast_mask = stdlib_sum(b,m)
    end function
end function

pure module function stdlib_sum_4d_dim_csp( x , dim, mask ) result( s )
    complex(sp), intent(in) :: x(:,:,:,:)
    integer, intent(in):: dim
    logical, intent(in), optional :: mask(:,:,:,:)
    complex(sp) :: s(merge(size(x, 1), size(x, 2), mask=1<dim), merge(size(x, 2), size(x, 3), mask=2<dim), merge(size(x, 3),&
        & size(x, 4), mask=3<dim))
    integer :: j 

    if(.not.present(mask)) then
        if(dim<4)then
            do j = 1, size(x,dim=4)
                s(:, :, j) = stdlib_sum( x(:, :, :, j), dim=dim )
            end do
        else
            do j = 1, size(x,dim=1)
                s(j, :, :) = stdlib_sum( x(j, :, :, :), dim=3 )
            end do
        end if
    else 
        if(dim<4)then
            do j = 1, size(x,dim=4)
                s(:, :, j) = stdlib_sum( x(:, :, :, j), dim=dim, mask=mask(:, :, :, j) )
            end do
        else
            do j = 1, size(x,dim=1)
                s(j, :, :) = stdlib_sum( x(j, :, :, :), dim=3, mask=mask(j, :, :, :) )
            end do
        end if
    end if

end function
pure module function stdlib_sum_2d_cdp( x , mask ) result( s )
    complex(dp), intent(in) :: x(:,:)
    logical, intent(in), optional :: mask(:,:)
    complex(dp) :: s
    if(.not.present(mask)) then
        s = sum_recast(x,size(x,kind=ilp))
    else
        s = sum_recast_mask(x,mask,size(x,kind=ilp))
    end if
contains
    pure complex(dp) function sum_recast(b,n)
      integer(ilp), intent(in) :: n
      complex(dp), intent(in) :: b(n)
      sum_recast = stdlib_sum(b)
    end function
    pure complex(dp) function sum_recast_mask(b,m,n)
      integer(ilp), intent(in) :: n
      complex(dp), intent(in) :: b(n)
      logical, intent(in) :: m(n)
      sum_recast_mask = stdlib_sum(b,m)
    end function
end function

pure module function stdlib_sum_2d_dim_cdp( x , dim, mask ) result( s )
    complex(dp), intent(in) :: x(:,:)
    integer, intent(in):: dim
    logical, intent(in), optional :: mask(:,:)
    complex(dp) :: s(merge(size(x, 1), size(x, 2), mask=1<dim))
    integer :: j 

    if(.not.present(mask)) then
        if(dim<2)then
            do j = 1, size(x,dim=2)
                s(j) = stdlib_sum( x(:, j) )
            end do
        else
            do j = 1, size(x,dim=1)
                s(j) = stdlib_sum( x(j, :) )
            end do
        end if
    else 
        if(dim<2)then
            do j = 1, size(x,dim=2)
                s(j) = stdlib_sum( x(:, j), mask=mask(:, j) )
            end do
        else
            do j = 1, size(x,dim=1)
                s(j) = stdlib_sum( x(j, :), mask=mask(j, :) )
            end do
        end if
    end if

end function
pure module function stdlib_sum_3d_cdp( x , mask ) result( s )
    complex(dp), intent(in) :: x(:,:,:)
    logical, intent(in), optional :: mask(:,:,:)
    complex(dp) :: s
    if(.not.present(mask)) then
        s = sum_recast(x,size(x,kind=ilp))
    else
        s = sum_recast_mask(x,mask,size(x,kind=ilp))
    end if
contains
    pure complex(dp) function sum_recast(b,n)
      integer(ilp), intent(in) :: n
      complex(dp), intent(in) :: b(n)
      sum_recast = stdlib_sum(b)
    end function
    pure complex(dp) function sum_recast_mask(b,m,n)
      integer(ilp), intent(in) :: n
      complex(dp), intent(in) :: b(n)
      logical, intent(in) :: m(n)
      sum_recast_mask = stdlib_sum(b,m)
    end function
end function

pure module function stdlib_sum_3d_dim_cdp( x , dim, mask ) result( s )
    complex(dp), intent(in) :: x(:,:,:)
    integer, intent(in):: dim
    logical, intent(in), optional :: mask(:,:,:)
    complex(dp) :: s(merge(size(x, 1), size(x, 2), mask=1<dim), merge(size(x, 2), size(x, 3), mask=2<dim))
    integer :: j 

    if(.not.present(mask)) then
        if(dim<3)then
            do j = 1, size(x,dim=3)
                s(:, j) = stdlib_sum( x(:, :, j), dim=dim )
            end do
        else
            do j = 1, size(x,dim=1)
                s(j, :) = stdlib_sum( x(j, :, :), dim=2 )
            end do
        end if
    else 
        if(dim<3)then
            do j = 1, size(x,dim=3)
                s(:, j) = stdlib_sum( x(:, :, j), dim=dim, mask=mask(:, :, j) )
            end do
        else
            do j = 1, size(x,dim=1)
                s(j, :) = stdlib_sum( x(j, :, :), dim=2, mask=mask(j, :, :) )
            end do
        end if
    end if

end function
pure module function stdlib_sum_4d_cdp( x , mask ) result( s )
    complex(dp), intent(in) :: x(:,:,:,:)
    logical, intent(in), optional :: mask(:,:,:,:)
    complex(dp) :: s
    if(.not.present(mask)) then
        s = sum_recast(x,size(x,kind=ilp))
    else
        s = sum_recast_mask(x,mask,size(x,kind=ilp))
    end if
contains
    pure complex(dp) function sum_recast(b,n)
      integer(ilp), intent(in) :: n
      complex(dp), intent(in) :: b(n)
      sum_recast = stdlib_sum(b)
    end function
    pure complex(dp) function sum_recast_mask(b,m,n)
      integer(ilp), intent(in) :: n
      complex(dp), intent(in) :: b(n)
      logical, intent(in) :: m(n)
      sum_recast_mask = stdlib_sum(b,m)
    end function
end function

pure module function stdlib_sum_4d_dim_cdp( x , dim, mask ) result( s )
    complex(dp), intent(in) :: x(:,:,:,:)
    integer, intent(in):: dim
    logical, intent(in), optional :: mask(:,:,:,:)
    complex(dp) :: s(merge(size(x, 1), size(x, 2), mask=1<dim), merge(size(x, 2), size(x, 3), mask=2<dim), merge(size(x, 3),&
        & size(x, 4), mask=3<dim))
    integer :: j 

    if(.not.present(mask)) then
        if(dim<4)then
            do j = 1, size(x,dim=4)
                s(:, :, j) = stdlib_sum( x(:, :, :, j), dim=dim )
            end do
        else
            do j = 1, size(x,dim=1)
                s(j, :, :) = stdlib_sum( x(j, :, :, :), dim=3 )
            end do
        end if
    else 
        if(dim<4)then
            do j = 1, size(x,dim=4)
                s(:, :, j) = stdlib_sum( x(:, :, :, j), dim=dim, mask=mask(:, :, :, j) )
            end do
        else
            do j = 1, size(x,dim=1)
                s(j, :, :) = stdlib_sum( x(j, :, :, :), dim=3, mask=mask(j, :, :, :) )
            end do
        end if
    end if

end function

pure module function stdlib_sum_kahan_2d_sp( x , mask ) result( s )
    real(sp), intent(in) :: x(:,:)
    logical, intent(in), optional :: mask(:,:)
    real(sp) :: s
    if(.not.present(mask)) then
        s = sum_recast(x,size(x,kind=ilp))
    else
        s = sum_recast_mask(x,mask,size(x,kind=ilp))
    end if
contains
    pure real(sp) function sum_recast(b,n)
      integer(ilp), intent(in) :: n
      real(sp), intent(in) :: b(n)
      sum_recast = stdlib_sum_kahan(b)
    end function
    pure real(sp) function sum_recast_mask(b,m,n)
      integer(ilp), intent(in) :: n
      real(sp), intent(in) :: b(n)
      logical, intent(in) :: m(n)
      sum_recast_mask = stdlib_sum_kahan(b,m)
    end function
end function

pure module function stdlib_sum_kahan_2d_dim_sp( x , dim, mask ) result( s )
    real(sp), intent(in) :: x(:,:)
    integer, intent(in):: dim
    logical, intent(in), optional :: mask(:,:)
    real(sp) :: s(merge(size(x, 1), size(x, 2), mask=1<dim))
    integer :: j 

    if(.not.present(mask)) then
        if(dim<2)then
            do j = 1, size(x,dim=2)
                s(j) = stdlib_sum_kahan( x(:, j) )
            end do
        else
            do j = 1, size(x,dim=1)
                s(j) = stdlib_sum_kahan( x(j, :) )
            end do
        end if
    else 
        if(dim<2)then
            do j = 1, size(x,dim=2)
                s(j) = stdlib_sum_kahan( x(:, j), mask=mask(:, j) )
            end do
        else
            do j = 1, size(x,dim=1)
                s(j) = stdlib_sum_kahan( x(j, :), mask=mask(j, :) )
            end do
        end if
    end if

end function
pure module function stdlib_sum_kahan_3d_sp( x , mask ) result( s )
    real(sp), intent(in) :: x(:,:,:)
    logical, intent(in), optional :: mask(:,:,:)
    real(sp) :: s
    if(.not.present(mask)) then
        s = sum_recast(x,size(x,kind=ilp))
    else
        s = sum_recast_mask(x,mask,size(x,kind=ilp))
    end if
contains
    pure real(sp) function sum_recast(b,n)
      integer(ilp), intent(in) :: n
      real(sp), intent(in) :: b(n)
      sum_recast = stdlib_sum_kahan(b)
    end function
    pure real(sp) function sum_recast_mask(b,m,n)
      integer(ilp), intent(in) :: n
      real(sp), intent(in) :: b(n)
      logical, intent(in) :: m(n)
      sum_recast_mask = stdlib_sum_kahan(b,m)
    end function
end function

pure module function stdlib_sum_kahan_3d_dim_sp( x , dim, mask ) result( s )
    real(sp), intent(in) :: x(:,:,:)
    integer, intent(in):: dim
    logical, intent(in), optional :: mask(:,:,:)
    real(sp) :: s(merge(size(x, 1), size(x, 2), mask=1<dim), merge(size(x, 2), size(x, 3), mask=2<dim))
    integer :: j 

    if(.not.present(mask)) then
        if(dim<3)then
            do j = 1, size(x,dim=3)
                s(:, j) = stdlib_sum_kahan( x(:, :, j), dim=dim )
            end do
        else
            do j = 1, size(x,dim=1)
                s(j, :) = stdlib_sum_kahan( x(j, :, :), dim=2 )
            end do
        end if
    else 
        if(dim<3)then
            do j = 1, size(x,dim=3)
                s(:, j) = stdlib_sum_kahan( x(:, :, j), dim=dim, mask=mask(:, :, j) )
            end do
        else
            do j = 1, size(x,dim=1)
                s(j, :) = stdlib_sum_kahan( x(j, :, :), dim=2, mask=mask(j, :, :) )
            end do
        end if
    end if

end function
pure module function stdlib_sum_kahan_4d_sp( x , mask ) result( s )
    real(sp), intent(in) :: x(:,:,:,:)
    logical, intent(in), optional :: mask(:,:,:,:)
    real(sp) :: s
    if(.not.present(mask)) then
        s = sum_recast(x,size(x,kind=ilp))
    else
        s = sum_recast_mask(x,mask,size(x,kind=ilp))
    end if
contains
    pure real(sp) function sum_recast(b,n)
      integer(ilp), intent(in) :: n
      real(sp), intent(in) :: b(n)
      sum_recast = stdlib_sum_kahan(b)
    end function
    pure real(sp) function sum_recast_mask(b,m,n)
      integer(ilp), intent(in) :: n
      real(sp), intent(in) :: b(n)
      logical, intent(in) :: m(n)
      sum_recast_mask = stdlib_sum_kahan(b,m)
    end function
end function

pure module function stdlib_sum_kahan_4d_dim_sp( x , dim, mask ) result( s )
    real(sp), intent(in) :: x(:,:,:,:)
    integer, intent(in):: dim
    logical, intent(in), optional :: mask(:,:,:,:)
    real(sp) :: s(merge(size(x, 1), size(x, 2), mask=1<dim), merge(size(x, 2), size(x, 3), mask=2<dim), merge(size(x, 3), size(x,&
        & 4), mask=3<dim))
    integer :: j 

    if(.not.present(mask)) then
        if(dim<4)then
            do j = 1, size(x,dim=4)
                s(:, :, j) = stdlib_sum_kahan( x(:, :, :, j), dim=dim )
            end do
        else
            do j = 1, size(x,dim=1)
                s(j, :, :) = stdlib_sum_kahan( x(j, :, :, :), dim=3 )
            end do
        end if
    else 
        if(dim<4)then
            do j = 1, size(x,dim=4)
                s(:, :, j) = stdlib_sum_kahan( x(:, :, :, j), dim=dim, mask=mask(:, :, :, j) )
            end do
        else
            do j = 1, size(x,dim=1)
                s(j, :, :) = stdlib_sum_kahan( x(j, :, :, :), dim=3, mask=mask(j, :, :, :) )
            end do
        end if
    end if

end function
pure module function stdlib_sum_kahan_2d_dp( x , mask ) result( s )
    real(dp), intent(in) :: x(:,:)
    logical, intent(in), optional :: mask(:,:)
    real(dp) :: s
    if(.not.present(mask)) then
        s = sum_recast(x,size(x,kind=ilp))
    else
        s = sum_recast_mask(x,mask,size(x,kind=ilp))
    end if
contains
    pure real(dp) function sum_recast(b,n)
      integer(ilp), intent(in) :: n
      real(dp), intent(in) :: b(n)
      sum_recast = stdlib_sum_kahan(b)
    end function
    pure real(dp) function sum_recast_mask(b,m,n)
      integer(ilp), intent(in) :: n
      real(dp), intent(in) :: b(n)
      logical, intent(in) :: m(n)
      sum_recast_mask = stdlib_sum_kahan(b,m)
    end function
end function

pure module function stdlib_sum_kahan_2d_dim_dp( x , dim, mask ) result( s )
    real(dp), intent(in) :: x(:,:)
    integer, intent(in):: dim
    logical, intent(in), optional :: mask(:,:)
    real(dp) :: s(merge(size(x, 1), size(x, 2), mask=1<dim))
    integer :: j 

    if(.not.present(mask)) then
        if(dim<2)then
            do j = 1, size(x,dim=2)
                s(j) = stdlib_sum_kahan( x(:, j) )
            end do
        else
            do j = 1, size(x,dim=1)
                s(j) = stdlib_sum_kahan( x(j, :) )
            end do
        end if
    else 
        if(dim<2)then
            do j = 1, size(x,dim=2)
                s(j) = stdlib_sum_kahan( x(:, j), mask=mask(:, j) )
            end do
        else
            do j = 1, size(x,dim=1)
                s(j) = stdlib_sum_kahan( x(j, :), mask=mask(j, :) )
            end do
        end if
    end if

end function
pure module function stdlib_sum_kahan_3d_dp( x , mask ) result( s )
    real(dp), intent(in) :: x(:,:,:)
    logical, intent(in), optional :: mask(:,:,:)
    real(dp) :: s
    if(.not.present(mask)) then
        s = sum_recast(x,size(x,kind=ilp))
    else
        s = sum_recast_mask(x,mask,size(x,kind=ilp))
    end if
contains
    pure real(dp) function sum_recast(b,n)
      integer(ilp), intent(in) :: n
      real(dp), intent(in) :: b(n)
      sum_recast = stdlib_sum_kahan(b)
    end function
    pure real(dp) function sum_recast_mask(b,m,n)
      integer(ilp), intent(in) :: n
      real(dp), intent(in) :: b(n)
      logical, intent(in) :: m(n)
      sum_recast_mask = stdlib_sum_kahan(b,m)
    end function
end function

pure module function stdlib_sum_kahan_3d_dim_dp( x , dim, mask ) result( s )
    real(dp), intent(in) :: x(:,:,:)
    integer, intent(in):: dim
    logical, intent(in), optional :: mask(:,:,:)
    real(dp) :: s(merge(size(x, 1), size(x, 2), mask=1<dim), merge(size(x, 2), size(x, 3), mask=2<dim))
    integer :: j 

    if(.not.present(mask)) then
        if(dim<3)then
            do j = 1, size(x,dim=3)
                s(:, j) = stdlib_sum_kahan( x(:, :, j), dim=dim )
            end do
        else
            do j = 1, size(x,dim=1)
                s(j, :) = stdlib_sum_kahan( x(j, :, :), dim=2 )
            end do
        end if
    else 
        if(dim<3)then
            do j = 1, size(x,dim=3)
                s(:, j) = stdlib_sum_kahan( x(:, :, j), dim=dim, mask=mask(:, :, j) )
            end do
        else
            do j = 1, size(x,dim=1)
                s(j, :) = stdlib_sum_kahan( x(j, :, :), dim=2, mask=mask(j, :, :) )
            end do
        end if
    end if

end function
pure module function stdlib_sum_kahan_4d_dp( x , mask ) result( s )
    real(dp), intent(in) :: x(:,:,:,:)
    logical, intent(in), optional :: mask(:,:,:,:)
    real(dp) :: s
    if(.not.present(mask)) then
        s = sum_recast(x,size(x,kind=ilp))
    else
        s = sum_recast_mask(x,mask,size(x,kind=ilp))
    end if
contains
    pure real(dp) function sum_recast(b,n)
      integer(ilp), intent(in) :: n
      real(dp), intent(in) :: b(n)
      sum_recast = stdlib_sum_kahan(b)
    end function
    pure real(dp) function sum_recast_mask(b,m,n)
      integer(ilp), intent(in) :: n
      real(dp), intent(in) :: b(n)
      logical, intent(in) :: m(n)
      sum_recast_mask = stdlib_sum_kahan(b,m)
    end function
end function

pure module function stdlib_sum_kahan_4d_dim_dp( x , dim, mask ) result( s )
    real(dp), intent(in) :: x(:,:,:,:)
    integer, intent(in):: dim
    logical, intent(in), optional :: mask(:,:,:,:)
    real(dp) :: s(merge(size(x, 1), size(x, 2), mask=1<dim), merge(size(x, 2), size(x, 3), mask=2<dim), merge(size(x, 3), size(x,&
        & 4), mask=3<dim))
    integer :: j 

    if(.not.present(mask)) then
        if(dim<4)then
            do j = 1, size(x,dim=4)
                s(:, :, j) = stdlib_sum_kahan( x(:, :, :, j), dim=dim )
            end do
        else
            do j = 1, size(x,dim=1)
                s(j, :, :) = stdlib_sum_kahan( x(j, :, :, :), dim=3 )
            end do
        end if
    else 
        if(dim<4)then
            do j = 1, size(x,dim=4)
                s(:, :, j) = stdlib_sum_kahan( x(:, :, :, j), dim=dim, mask=mask(:, :, :, j) )
            end do
        else
            do j = 1, size(x,dim=1)
                s(j, :, :) = stdlib_sum_kahan( x(j, :, :, :), dim=3, mask=mask(j, :, :, :) )
            end do
        end if
    end if

end function
pure module function stdlib_sum_kahan_2d_csp( x , mask ) result( s )
    complex(sp), intent(in) :: x(:,:)
    logical, intent(in), optional :: mask(:,:)
    complex(sp) :: s
    if(.not.present(mask)) then
        s = sum_recast(x,size(x,kind=ilp))
    else
        s = sum_recast_mask(x,mask,size(x,kind=ilp))
    end if
contains
    pure complex(sp) function sum_recast(b,n)
      integer(ilp), intent(in) :: n
      complex(sp), intent(in) :: b(n)
      sum_recast = stdlib_sum_kahan(b)
    end function
    pure complex(sp) function sum_recast_mask(b,m,n)
      integer(ilp), intent(in) :: n
      complex(sp), intent(in) :: b(n)
      logical, intent(in) :: m(n)
      sum_recast_mask = stdlib_sum_kahan(b,m)
    end function
end function

pure module function stdlib_sum_kahan_2d_dim_csp( x , dim, mask ) result( s )
    complex(sp), intent(in) :: x(:,:)
    integer, intent(in):: dim
    logical, intent(in), optional :: mask(:,:)
    complex(sp) :: s(merge(size(x, 1), size(x, 2), mask=1<dim))
    integer :: j 

    if(.not.present(mask)) then
        if(dim<2)then
            do j = 1, size(x,dim=2)
                s(j) = stdlib_sum_kahan( x(:, j) )
            end do
        else
            do j = 1, size(x,dim=1)
                s(j) = stdlib_sum_kahan( x(j, :) )
            end do
        end if
    else 
        if(dim<2)then
            do j = 1, size(x,dim=2)
                s(j) = stdlib_sum_kahan( x(:, j), mask=mask(:, j) )
            end do
        else
            do j = 1, size(x,dim=1)
                s(j) = stdlib_sum_kahan( x(j, :), mask=mask(j, :) )
            end do
        end if
    end if

end function
pure module function stdlib_sum_kahan_3d_csp( x , mask ) result( s )
    complex(sp), intent(in) :: x(:,:,:)
    logical, intent(in), optional :: mask(:,:,:)
    complex(sp) :: s
    if(.not.present(mask)) then
        s = sum_recast(x,size(x,kind=ilp))
    else
        s = sum_recast_mask(x,mask,size(x,kind=ilp))
    end if
contains
    pure complex(sp) function sum_recast(b,n)
      integer(ilp), intent(in) :: n
      complex(sp), intent(in) :: b(n)
      sum_recast = stdlib_sum_kahan(b)
    end function
    pure complex(sp) function sum_recast_mask(b,m,n)
      integer(ilp), intent(in) :: n
      complex(sp), intent(in) :: b(n)
      logical, intent(in) :: m(n)
      sum_recast_mask = stdlib_sum_kahan(b,m)
    end function
end function

pure module function stdlib_sum_kahan_3d_dim_csp( x , dim, mask ) result( s )
    complex(sp), intent(in) :: x(:,:,:)
    integer, intent(in):: dim
    logical, intent(in), optional :: mask(:,:,:)
    complex(sp) :: s(merge(size(x, 1), size(x, 2), mask=1<dim), merge(size(x, 2), size(x, 3), mask=2<dim))
    integer :: j 

    if(.not.present(mask)) then
        if(dim<3)then
            do j = 1, size(x,dim=3)
                s(:, j) = stdlib_sum_kahan( x(:, :, j), dim=dim )
            end do
        else
            do j = 1, size(x,dim=1)
                s(j, :) = stdlib_sum_kahan( x(j, :, :), dim=2 )
            end do
        end if
    else 
        if(dim<3)then
            do j = 1, size(x,dim=3)
                s(:, j) = stdlib_sum_kahan( x(:, :, j), dim=dim, mask=mask(:, :, j) )
            end do
        else
            do j = 1, size(x,dim=1)
                s(j, :) = stdlib_sum_kahan( x(j, :, :), dim=2, mask=mask(j, :, :) )
            end do
        end if
    end if

end function
pure module function stdlib_sum_kahan_4d_csp( x , mask ) result( s )
    complex(sp), intent(in) :: x(:,:,:,:)
    logical, intent(in), optional :: mask(:,:,:,:)
    complex(sp) :: s
    if(.not.present(mask)) then
        s = sum_recast(x,size(x,kind=ilp))
    else
        s = sum_recast_mask(x,mask,size(x,kind=ilp))
    end if
contains
    pure complex(sp) function sum_recast(b,n)
      integer(ilp), intent(in) :: n
      complex(sp), intent(in) :: b(n)
      sum_recast = stdlib_sum_kahan(b)
    end function
    pure complex(sp) function sum_recast_mask(b,m,n)
      integer(ilp), intent(in) :: n
      complex(sp), intent(in) :: b(n)
      logical, intent(in) :: m(n)
      sum_recast_mask = stdlib_sum_kahan(b,m)
    end function
end function

pure module function stdlib_sum_kahan_4d_dim_csp( x , dim, mask ) result( s )
    complex(sp), intent(in) :: x(:,:,:,:)
    integer, intent(in):: dim
    logical, intent(in), optional :: mask(:,:,:,:)
    complex(sp) :: s(merge(size(x, 1), size(x, 2), mask=1<dim), merge(size(x, 2), size(x, 3), mask=2<dim), merge(size(x, 3),&
        & size(x, 4), mask=3<dim))
    integer :: j 

    if(.not.present(mask)) then
        if(dim<4)then
            do j = 1, size(x,dim=4)
                s(:, :, j) = stdlib_sum_kahan( x(:, :, :, j), dim=dim )
            end do
        else
            do j = 1, size(x,dim=1)
                s(j, :, :) = stdlib_sum_kahan( x(j, :, :, :), dim=3 )
            end do
        end if
    else 
        if(dim<4)then
            do j = 1, size(x,dim=4)
                s(:, :, j) = stdlib_sum_kahan( x(:, :, :, j), dim=dim, mask=mask(:, :, :, j) )
            end do
        else
            do j = 1, size(x,dim=1)
                s(j, :, :) = stdlib_sum_kahan( x(j, :, :, :), dim=3, mask=mask(j, :, :, :) )
            end do
        end if
    end if

end function
pure module function stdlib_sum_kahan_2d_cdp( x , mask ) result( s )
    complex(dp), intent(in) :: x(:,:)
    logical, intent(in), optional :: mask(:,:)
    complex(dp) :: s
    if(.not.present(mask)) then
        s = sum_recast(x,size(x,kind=ilp))
    else
        s = sum_recast_mask(x,mask,size(x,kind=ilp))
    end if
contains
    pure complex(dp) function sum_recast(b,n)
      integer(ilp), intent(in) :: n
      complex(dp), intent(in) :: b(n)
      sum_recast = stdlib_sum_kahan(b)
    end function
    pure complex(dp) function sum_recast_mask(b,m,n)
      integer(ilp), intent(in) :: n
      complex(dp), intent(in) :: b(n)
      logical, intent(in) :: m(n)
      sum_recast_mask = stdlib_sum_kahan(b,m)
    end function
end function

pure module function stdlib_sum_kahan_2d_dim_cdp( x , dim, mask ) result( s )
    complex(dp), intent(in) :: x(:,:)
    integer, intent(in):: dim
    logical, intent(in), optional :: mask(:,:)
    complex(dp) :: s(merge(size(x, 1), size(x, 2), mask=1<dim))
    integer :: j 

    if(.not.present(mask)) then
        if(dim<2)then
            do j = 1, size(x,dim=2)
                s(j) = stdlib_sum_kahan( x(:, j) )
            end do
        else
            do j = 1, size(x,dim=1)
                s(j) = stdlib_sum_kahan( x(j, :) )
            end do
        end if
    else 
        if(dim<2)then
            do j = 1, size(x,dim=2)
                s(j) = stdlib_sum_kahan( x(:, j), mask=mask(:, j) )
            end do
        else
            do j = 1, size(x,dim=1)
                s(j) = stdlib_sum_kahan( x(j, :), mask=mask(j, :) )
            end do
        end if
    end if

end function
pure module function stdlib_sum_kahan_3d_cdp( x , mask ) result( s )
    complex(dp), intent(in) :: x(:,:,:)
    logical, intent(in), optional :: mask(:,:,:)
    complex(dp) :: s
    if(.not.present(mask)) then
        s = sum_recast(x,size(x,kind=ilp))
    else
        s = sum_recast_mask(x,mask,size(x,kind=ilp))
    end if
contains
    pure complex(dp) function sum_recast(b,n)
      integer(ilp), intent(in) :: n
      complex(dp), intent(in) :: b(n)
      sum_recast = stdlib_sum_kahan(b)
    end function
    pure complex(dp) function sum_recast_mask(b,m,n)
      integer(ilp), intent(in) :: n
      complex(dp), intent(in) :: b(n)
      logical, intent(in) :: m(n)
      sum_recast_mask = stdlib_sum_kahan(b,m)
    end function
end function

pure module function stdlib_sum_kahan_3d_dim_cdp( x , dim, mask ) result( s )
    complex(dp), intent(in) :: x(:,:,:)
    integer, intent(in):: dim
    logical, intent(in), optional :: mask(:,:,:)
    complex(dp) :: s(merge(size(x, 1), size(x, 2), mask=1<dim), merge(size(x, 2), size(x, 3), mask=2<dim))
    integer :: j 

    if(.not.present(mask)) then
        if(dim<3)then
            do j = 1, size(x,dim=3)
                s(:, j) = stdlib_sum_kahan( x(:, :, j), dim=dim )
            end do
        else
            do j = 1, size(x,dim=1)
                s(j, :) = stdlib_sum_kahan( x(j, :, :), dim=2 )
            end do
        end if
    else 
        if(dim<3)then
            do j = 1, size(x,dim=3)
                s(:, j) = stdlib_sum_kahan( x(:, :, j), dim=dim, mask=mask(:, :, j) )
            end do
        else
            do j = 1, size(x,dim=1)
                s(j, :) = stdlib_sum_kahan( x(j, :, :), dim=2, mask=mask(j, :, :) )
            end do
        end if
    end if

end function
pure module function stdlib_sum_kahan_4d_cdp( x , mask ) result( s )
    complex(dp), intent(in) :: x(:,:,:,:)
    logical, intent(in), optional :: mask(:,:,:,:)
    complex(dp) :: s
    if(.not.present(mask)) then
        s = sum_recast(x,size(x,kind=ilp))
    else
        s = sum_recast_mask(x,mask,size(x,kind=ilp))
    end if
contains
    pure complex(dp) function sum_recast(b,n)
      integer(ilp), intent(in) :: n
      complex(dp), intent(in) :: b(n)
      sum_recast = stdlib_sum_kahan(b)
    end function
    pure complex(dp) function sum_recast_mask(b,m,n)
      integer(ilp), intent(in) :: n
      complex(dp), intent(in) :: b(n)
      logical, intent(in) :: m(n)
      sum_recast_mask = stdlib_sum_kahan(b,m)
    end function
end function

pure module function stdlib_sum_kahan_4d_dim_cdp( x , dim, mask ) result( s )
    complex(dp), intent(in) :: x(:,:,:,:)
    integer, intent(in):: dim
    logical, intent(in), optional :: mask(:,:,:,:)
    complex(dp) :: s(merge(size(x, 1), size(x, 2), mask=1<dim), merge(size(x, 2), size(x, 3), mask=2<dim), merge(size(x, 3),&
        & size(x, 4), mask=3<dim))
    integer :: j 

    if(.not.present(mask)) then
        if(dim<4)then
            do j = 1, size(x,dim=4)
                s(:, :, j) = stdlib_sum_kahan( x(:, :, :, j), dim=dim )
            end do
        else
            do j = 1, size(x,dim=1)
                s(j, :, :) = stdlib_sum_kahan( x(j, :, :, :), dim=3 )
            end do
        end if
    else 
        if(dim<4)then
            do j = 1, size(x,dim=4)
                s(:, :, j) = stdlib_sum_kahan( x(:, :, :, j), dim=dim, mask=mask(:, :, :, j) )
            end do
        else
            do j = 1, size(x,dim=1)
                s(j, :, :) = stdlib_sum_kahan( x(j, :, :, :), dim=3, mask=mask(j, :, :, :) )
            end do
        end if
    end if

end function

end submodule stdlib_intrinsics_sum
