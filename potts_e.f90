module big_uint_mod
    use iso_fortran_env, only: int64
    implicit none

    private
    public :: big_uint, bu_from_int, bu_is_zero, bu_add_inplace, bu_subtract
    public :: bu_compare, bu_pow_small, bu_to_string
    public :: bu_limbs

    ! 7 limbs with base 1e9 gives up to 63 decimal digits.
    integer, parameter :: bu_limbs = 7
    integer(int64), parameter :: bu_base = 1000000000_int64

    type :: big_uint
        integer(int64) :: limbs(bu_limbs) = 0_int64
    end type big_uint

contains

    function bu_from_int(v) result(x)
        integer, intent(in) :: v
        type(big_uint) :: x
        integer(int64) :: tmp

        tmp = int(v, int64)
        x%limbs = 0_int64
        x%limbs(1) = modulo(tmp, bu_base)
        tmp = tmp / bu_base
        if (tmp /= 0_int64) then
            x%limbs(2) = tmp
        end if
    end function bu_from_int

    logical function bu_is_zero(x)
        type(big_uint), intent(in) :: x
        bu_is_zero = all(x%limbs == 0_int64)
    end function bu_is_zero

    subroutine bu_add_inplace(x, y)
        type(big_uint), intent(inout) :: x
        type(big_uint), intent(in) :: y
        integer :: i
        integer(int64) :: carry, temp

        carry = 0_int64
        do i = 1, bu_limbs
            temp = x%limbs(i) + y%limbs(i) + carry
            if (temp >= bu_base) then
                x%limbs(i) = temp - bu_base
                carry = 1_int64
            else
                x%limbs(i) = temp
                carry = 0_int64
            end if
        end do
    end subroutine bu_add_inplace

    function bu_subtract(a, b) result(c)
        type(big_uint), intent(in) :: a, b
        type(big_uint) :: c
        integer :: i
        integer(int64) :: borrow, temp

        c%limbs = 0_int64
        borrow = 0_int64
        do i = 1, bu_limbs
            temp = a%limbs(i) - b%limbs(i) - borrow
            if (temp < 0_int64) then
                c%limbs(i) = temp + bu_base
                borrow = 1_int64
            else
                c%limbs(i) = temp
                borrow = 0_int64
            end if
        end do
    end function bu_subtract

    integer function bu_compare(a, b)
        type(big_uint), intent(in) :: a, b
        integer :: i

        bu_compare = 0
        do i = bu_limbs, 1, -1
            if (a%limbs(i) < b%limbs(i)) then
                bu_compare = -1
                return
            else if (a%limbs(i) > b%limbs(i)) then
                bu_compare = 1
                return
            end if
        end do
    end function bu_compare

    function bu_mul_small(a, m) result(c)
        type(big_uint), intent(in) :: a
        integer, intent(in) :: m
        type(big_uint) :: c
        integer :: i
        integer(int64) :: mm, carry, temp

        mm = int(m, int64)
        c%limbs = 0_int64
        carry = 0_int64
        do i = 1, bu_limbs
            temp = a%limbs(i) * mm + carry
            c%limbs(i) = modulo(temp, bu_base)
            carry = temp / bu_base
        end do
    end function bu_mul_small

    function bu_pow_small(base_val, exponent) result(res)
        integer, intent(in) :: base_val, exponent
        type(big_uint) :: res
        integer :: i

        res = bu_from_int(1)
        do i = 1, exponent
            res = bu_mul_small(res, base_val)
        end do
    end function bu_pow_small

    function bu_to_string(x) result(str)
        type(big_uint), intent(in) :: x
        character(len=:), allocatable :: str
        character(len=32) :: part
        integer :: i, hi

        hi = 1
        do i = bu_limbs, 1, -1
            if (x%limbs(i) /= 0_int64) then
                hi = i
                exit
            end if
        end do

        write(part, '(I0)') x%limbs(hi)
        str = trim(adjustl(part))
        do i = hi - 1, 1, -1
            write(part, '(I9.9)') x%limbs(i)
            str = str // part(1:9)
        end do
    end function bu_to_string

end module big_uint_mod

module potts_e
    use big_uint_mod
    use iso_fortran_env, only: real64
    implicit none

    private
    public :: calculate_coefficients, write_coefficients

    ! i in [1, l], j in [1, n], k in [1, q], b in [0, max_bonds - 1], index in [1, n_intra_states]
    integer, parameter :: l = 8, n = 8, q = 3
    integer, parameter :: max_bonds = l * n * 2 + 1
    integer, parameter :: n_intra_states = q ** l
    
    integer :: interaction(q, q)
    character :: boundary

    ! o(1+b, index) is the number of configurations with b bonds and the index-th intra-layer state,
    ! where a bond is defined as an edge connecting two spins with the *SAME* state
    type(big_uint), target, allocatable :: o_storage(:, :)
    type(big_uint), target, allocatable :: oo_storage(:, :)
    type(big_uint), pointer :: o(:,:) => null()
    type(big_uint), pointer :: oo(:,:) => null()
    type(big_uint), allocatable :: final_coefficients(:)
    integer, allocatable :: state_table(:, :)
    integer, allocatable :: intra_bonds_table(:)
    integer, allocatable :: transition_index(:, :, :)
    integer, allocatable :: transition_db(:, :, :)
    integer, allocatable :: pow_q(:)

contains

    subroutine initialize()
        integer :: index, k

        interaction(:, :) = 0
        do k = 1, q
            interaction(k, k) = 1
        end do

        call ensure_big_capacity()
        call allocate_work_arrays()
        call build_lookup_tables()

        o => o_storage
        oo => oo_storage

        o(:, :) = bu_from_int(0)
        do index = 1, n_intra_states
            o(1+intra_bonds_table(index), index) = bu_from_int(1)
        end do
    end subroutine initialize

    subroutine intra_layer(i)
        integer, intent(in) :: i
        integer :: k
        integer :: index, db, src_b, old_index
        oo(:, :) = bu_from_int(0)
        !$omp parallel do default(shared) private(index, k, db, old_index, src_b) schedule(static)
        do index = 1, n_intra_states
            do k = 1, q
                db = transition_db(i, k, index)
                old_index = transition_index(i, k, index)
                do src_b = 1, max_bonds - db
                    call bu_add_inplace(oo(src_b+db, index), o(src_b, old_index))
                end do
            end do
        end do
        !$omp end parallel do
        call swap_arrays()
    end subroutine intra_layer
    
    subroutine finalize_layer()
        integer :: index, db, src_b
        oo(:, :) = bu_from_int(0)
        !$omp parallel do default(shared) private(index, db, src_b) schedule(static)
        do index = 1, n_intra_states
            db = intra_bonds_table(index)
            do src_b = 1, max_bonds - db
                call bu_add_inplace(oo(src_b+db, index), o(src_b, index))
            end do
        end do
        !$omp end parallel do
        call swap_arrays()
    end subroutine finalize_layer

    subroutine calculate_coefficients(boundary_val)
        character, intent(in) :: boundary_val
        integer :: i, j
        boundary = boundary_val

        print *, 'calculating coefficients for q-state Potts model on l * n lattice with'
        print *, 'l = ', l, '(', boundary, ')'
        print *, 'n = ', n, '(o)'
        print *, 'q = ', q
        call initialize()
        do j = 2, n
            do i = 1, l
                call intra_layer(i)
            end do
            call finalize_layer()
        end do
        call sum_over_states(o, final_coefficients)
    end subroutine calculate_coefficients

    subroutine write_coefficients()
        integer :: b, unit, iostat
        character(len = 100) :: filename
        integer :: cmp
        type(big_uint) :: temp, theoretical_total, calculated_total, difference
        character(len=1) :: sign_prefix
        
        write(filename, '(A,I0,A,I0,A,I0,A,A,A)') &
            'e_l', l, '_n', n, '_q', q, '_', boundary, '.csv'
        open(newunit = unit, file = trim(filename), status = 'replace', action = 'write', iostat = iostat)
        if (iostat /= 0) then
            error stop 'failed to open file'
        end if
        write(unit, '(A)') '# bonds,count'
        
        do b = 0, max_bonds - 1
            temp = final_coefficients(1+b)
            if (.not. bu_is_zero(temp)) then
                write(unit, '(I0,A,A)') b, ',', bu_to_string(temp)
            end if
        end do
        
        close(unit)
        
        theoretical_total = bu_pow_small(q, l * n)
        calculated_total = sum_big_vector(final_coefficients)
        cmp = bu_compare(calculated_total, theoretical_total)
        sign_prefix = ''
        if (cmp >= 0) then
            difference = bu_subtract(calculated_total, theoretical_total)
        else
            difference = bu_subtract(theoretical_total, calculated_total)
            sign_prefix = '-'
        end if
        print '(A,A)', 'CSV file written: ', trim(filename)
        print '(A,A)', 'total (theo.): ', bu_to_string(theoretical_total)
        print '(A,A)', 'total (calc.): ', bu_to_string(calculated_total)
        print '(A,A,A)', 'difference:   ', sign_prefix, bu_to_string(difference)
    end subroutine write_coefficients

    subroutine swap_arrays()
        type(big_uint), pointer :: temp(:,:)
        temp => o
        o => oo
        oo => temp
    end subroutine swap_arrays

    subroutine allocate_work_arrays()
        if (.not. allocated(o_storage)) then
            allocate(o_storage(max_bonds, n_intra_states))
        end if
        if (.not. allocated(oo_storage)) then
            allocate(oo_storage(max_bonds, n_intra_states))
        end if
        if (.not. allocated(final_coefficients)) then
            allocate(final_coefficients(max_bonds))
        end if
        if (.not. allocated(state_table)) then
            allocate(state_table(l, n_intra_states))
        end if
        if (.not. allocated(intra_bonds_table)) then
            allocate(intra_bonds_table(n_intra_states))
        end if
        if (.not. allocated(transition_index)) then
            allocate(transition_index(l, q, n_intra_states))
        end if
        if (.not. allocated(transition_db)) then
            allocate(transition_db(l, q, n_intra_states))
        end if
        if (.not. allocated(pow_q)) then
            allocate(pow_q(l))
        end if
    end subroutine allocate_work_arrays

    subroutine ensure_big_capacity()
        integer :: required_digits, required_limbs
        real(real64) :: log10q

        log10q = log10(real(q, real64))
        required_digits = ceiling(real(l * n, real64) * log10q) + 1
        required_limbs = (required_digits + 8) / 9

        if (required_limbs > bu_limbs) then
            error stop 'big_uint capacity too small for current l, n, q'
        end if
    end subroutine ensure_big_capacity

    subroutine build_lookup_tables()
        integer :: index, i, k, idx0, max_i, b
        integer :: s_i

        pow_q(1) = 1
        do i = 2, l
            pow_q(i) = pow_q(i - 1) * q
        end do

        do index = 1, n_intra_states
            idx0 = index - 1
            do i = 1, l
                state_table(i, index) = mod(idx0 / pow_q(i), q) + 1
            end do
        end do

        do index = 1, n_intra_states
            do i = 1, l
                s_i = state_table(i, index)
                do k = 1, q
                    transition_db(i, k, index) = interaction(s_i, k)
                    transition_index(i, k, index) = index + (k - s_i) * pow_q(i)
                end do
            end do
        end do

        if (boundary == 'o') then
            max_i = l - 1
        else if (boundary == 'p') then
            max_i = l
        else
            error stop 'invalid boundary condition'
        end if

        do index = 1, n_intra_states
            b = 0
            do i = 1, max_i
                b = b + interaction(state_table(i, index), state_table(mod(i, l) + 1, index))
            end do
            intra_bonds_table(index) = b
        end do
    end subroutine build_lookup_tables

    subroutine sum_over_states(arr, colsum)
        type(big_uint), intent(in) :: arr(:, :)
        type(big_uint), intent(out) :: colsum(size(arr, 1))
        integer :: b, index

        colsum = bu_from_int(0)
        do index = 1, size(arr, 2)
            do b = 1, size(arr, 1)
                call bu_add_inplace(colsum(b), arr(b, index))
            end do
        end do
    end subroutine sum_over_states

    function sum_big_vector(vec) result(total)
        type(big_uint), intent(in) :: vec(:)
        type(big_uint) :: total
        integer :: i

        total = bu_from_int(0)
        do i = 1, size(vec)
            call bu_add_inplace(total, vec(i))
        end do
    end function sum_big_vector

end module potts_e