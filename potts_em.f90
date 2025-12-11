module potts_em
    implicit none

    private
    public :: calculate_coefficients, print_coefficients

    ! i in [1, l], j in [1, m], k in [1, q], b in [0, max_bonds - 1], index in [1, n_intra_states]
    integer, parameter :: l = 4, m = 4, q = 3
    integer, parameter :: max_bonds = l * m * 2 + 1
    integer, parameter :: n_intra_states = q ** l
    
    integer :: interaction(q, q)

    ! o(1+b, index) is the number of configurations with b bonds and the index-th intra-layer state,
    ! where a bond is defined as an edge connecting two spins with the *SAME* state
    integer(kind = 16), target :: o_storage(max_bonds, n_intra_states)
    integer(kind = 16), target :: oo_storage(max_bonds, n_intra_states)
    integer(kind = 16), pointer :: o(:,:) => null()
    integer(kind = 16), pointer :: oo(:,:) => null()
    integer(kind = 16) :: final_coefficients(max_bonds)

contains

    subroutine initialize(boundary)
        character(len = *), intent(in) :: boundary
        integer :: index, k

        interaction(:, :) = 0
        do concurrent (k = 1:q)
            interaction(k, k) = 1
        end do

        o => o_storage
        oo => oo_storage

        o(:, :) = 0
        do concurrent (index = 1:n_intra_states)
            o(1+intra_bonds(index, boundary), index) = 1
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
    
    subroutine finalize_layer(boundary)
        character(len = *), intent(in) :: boundary
        integer :: index, db
        oo(:, :) = 0
        do concurrent (index = 1:n_intra_states)
            db = intra_bonds(index, boundary)
            oo((1+db):, index) = oo((1+db):, index) + o(:(max_bonds-db), index)
        end do
        call swap_arrays()
    end subroutine finalize_layer

    subroutine calculate_coefficients(boundary)
        character(len = *), intent(in) :: boundary
        integer :: i, j

        print *, 'calculating coefficients for q-state Potts model on l * m lattice with'
        print *, 'l = ', l, '(', boundary, ' boundary)'
        print *, 'm = ', m, '(open boundary)'
        print *, 'q = ', q
        call initialize(boundary)
        do j = 2, m
            do i = 1, l
                call intra_layer(i)
            end do
            call finalize_layer(boundary)
        end do
        final_coefficients = sum(o, 2)
    end subroutine calculate_coefficients

    subroutine print_coefficients()
        integer :: b
        integer(kind = 16) :: theoretical_total

        print *, '      bonds                                    count'
        do b = 0, max_bonds - 1
            print *, b, final_coefficients(1+b)
        end do

        theoretical_total = q
        theoretical_total = theoretical_total ** (l * m)
        print *, 'total (theo.):', theoretical_total
        print *, 'total (calc.):', sum(final_coefficients)
        print *, 'difference:   ', sum(final_coefficients) - theoretical_total
    end subroutine print_coefficients

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

    pure function intra_bonds(index, boundary) result(b)
        character(len = *), intent(in) :: boundary
        integer, intent(in) :: index
        integer :: b
        integer :: states(l), i, max_i

        states = index_to_states(index)
        if (boundary == 'open') then
            max_i = l - 1
        else if (boundary == 'periodic') then
            max_i = l
        else
            error stop 'invalid boundary condition'
        end if

        b = 0
        do concurrent (i = 1:max_i)
            b = b + interaction(states(i), states(mod(i, l) + 1))
        end do
    end function intra_bonds

end module potts_em