module potts_e
    implicit none

    private
    public :: calculate_coefficients, print_coefficients, write_coefficients

    ! i in [1, l], j in [1, n], k in [1, q], b in [0, max_bonds - 1], index in [1, n_intra_states]
    integer, parameter :: l = 8, n = 8, q = 3
    integer, parameter :: max_bonds = l * n * 2 + 1
    integer, parameter :: n_intra_states = q ** l
    
    integer :: interaction(q, q)
    character :: boundary

    ! o(1+b, index) is the number of configurations with b bonds and the index-th intra-layer state,
    ! where a bond is defined as an edge connecting two spins with the *SAME* state
    integer(kind = 16), target :: o_storage(max_bonds, n_intra_states)
    integer(kind = 16), target :: oo_storage(max_bonds, n_intra_states)
    integer(kind = 16), pointer :: o(:,:) => null()
    integer(kind = 16), pointer :: oo(:,:) => null()
    integer(kind = 16) :: final_coefficients(max_bonds)

contains

    subroutine initialize()
        integer :: index, k

        interaction(:, :) = 0
        do concurrent (k = 1:q)
            interaction(k, k) = 1
        end do

        o => o_storage
        oo => oo_storage

        o(:, :) = 0
        do concurrent (index = 1:n_intra_states)
            o(1+intra_bonds(index), index) = 1
        end do
    end subroutine initialize

    subroutine intra_layer(i)
        integer, intent(in) :: i
        integer :: k, old_states(l), new_states(l)
        integer :: index, db
        oo(:, :) = 0
        do concurrent (index = 1:n_intra_states)
            new_states = index_to_states(index)
            old_states = new_states
            do k = 1, q
                old_states(i) = k
                db = interaction(new_states(i), old_states(i))
                oo((1+db):, index) = oo((1+db):, index) + o(:(max_bonds-db), states_to_index(old_states))
            end do
        end do
        call swap_arrays()
    end subroutine intra_layer
    
    subroutine finalize_layer()
        integer :: index, db
        oo(:, :) = 0
        do concurrent (index = 1:n_intra_states)
            db = intra_bonds(index)
            oo((1+db):, index) = oo((1+db):, index) + o(:(max_bonds-db), index)
        end do
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
        final_coefficients = sum(o, 2)
    end subroutine calculate_coefficients

    subroutine print_coefficients()
        integer :: b
        integer(kind = 16) :: temp, theoretical_total

        print *, '      bonds                                    count'
        do b = 0, max_bonds - 1
            temp = final_coefficients(1+b)
            if (temp /= 0) then
                print *, b, temp
            end if
        end do

        theoretical_total = q
        theoretical_total = theoretical_total ** (l * n)
        print *, 'total (theo.):', theoretical_total
        print *, 'total (calc.):', sum(final_coefficients)
        print *, 'difference:   ', sum(final_coefficients) - theoretical_total
    end subroutine print_coefficients

    subroutine write_coefficients()
        integer :: b, unit, iostat
        character(len = 100) :: filename
        integer(kind = 16) :: temp, theoretical_total
        
        write(filename, '(A,I0,A,I0,A,I0,A,A,A)') &
            'e_l', l, '_n', n, '_q', q, '_', boundary, '.csv'
        open(newunit = unit, file = trim(filename), status = 'replace', action = 'write', iostat = iostat)
        if (iostat /= 0) then
            error stop 'failed to open file'
        end if
        write(unit, '(A)') '# bonds,count'
        
        do b = 0, max_bonds - 1
            temp = final_coefficients(1+b)
            if (temp /= 0) then
                write(unit, '(I0,A,I0)') b, ',', temp
            end if
        end do
        
        close(unit)
        
        theoretical_total = q
        theoretical_total = theoretical_total ** (l * n)
        print '(A,A)', 'CSV file written: ', trim(filename)
        print *, 'total (theo.):', theoretical_total
        print *, 'total (calc.):', sum(final_coefficients)
        print *, 'difference:   ', sum(final_coefficients) - theoretical_total
    end subroutine write_coefficients

    subroutine swap_arrays()
        integer(kind = 16), pointer :: temp(:,:)
        temp => o
        o => oo
        oo => temp
    end subroutine swap_arrays

    pure function states_to_index(states) result(index)
        ! states is a 1D array of length l, with values in [1, q]
        integer, intent(in) :: states(l)
        integer :: i
        integer :: index

        index = 0
        do concurrent (i = 1:l)
            index = index + (states(i) - 1) * q ** (i-1)
        end do
        index = index + 1
    end function states_to_index

    pure function index_to_states(index) result(states)
        integer, intent(in) :: index
        integer :: states(l), i

        do concurrent (i = 1:l)
            states(i) = mod((index - 1) / q ** (i-1), q) + 1
        end do
    end function index_to_states

    pure function intra_bonds(index) result(b)
        integer, intent(in) :: index
        integer :: b
        integer :: states(l), i, max_i

        states = index_to_states(index)
        if (boundary == 'o') then
            max_i = l - 1
        else if (boundary == 'p') then
            max_i = l
        else
            error stop 'invalid boundary condition'
        end if

        b = 0
        do concurrent (i = 1:max_i)
            b = b + interaction(states(i), states(mod(i, l) + 1))
        end do
    end function intra_bonds

end module potts_e