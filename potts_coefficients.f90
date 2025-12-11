module potts_coefficients
    implicit none

    private
    public :: initialize, iterate_layers, print_coefficients

    ! i in [1, l], j in [1, m], k in [1, q], b in [0, max_bonds], index in [1, n_intra_states]
    integer(kind = 1), parameter :: l = 4, m = 4, q = 2
    integer, parameter :: max_bonds = l * m * 2 + 1
    integer, parameter :: n_intra_states = q ** l
    
    integer(kind = 1) :: interaction(q, q)

    ! o(1+b, index) is the number of configurations with b bonds and the index-th intra-layer state,
    ! where a bond is defined as an edge connecting two spins with the *SAME* state
    integer(kind = 16), target :: o_storage(max_bonds, n_intra_states)
    integer(kind = 16), target :: oo_storage(max_bonds, n_intra_states)
    integer(kind = 16), pointer :: o(:,:) => null()
    integer(kind = 16), pointer :: oo(:,:) => null()

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
        integer(kind = 1), intent(in) :: i
        integer(kind = 1) :: k, old_states(l), new_states(l)
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

    subroutine iterate_layers()
        integer(kind = 1) :: i, j
        do j = 2, m
            do i = 1, l
                call intra_layer(i)
            end do
            call finalize_layer()
        end do
    end subroutine iterate_layers

    subroutine print_coefficients()
        integer :: b
        integer(kind = 16) :: theoretical_total
        integer(kind = 16), allocatable :: total(:)
        allocate(total(max_bonds))
        total = sum(o, 2)

        print *, '      bonds                                    count'
        do b = 0, max_bonds
            print *, b, total(1+b)
        end do
        deallocate(total)

        theoretical_total = q
        theoretical_total = theoretical_total ** (l * m)
        print *, 'total configurations (theoretical):', theoretical_total
        print *, 'total configurations (calculation):', sum(o)
    end subroutine print_coefficients

    subroutine swap_arrays()
        integer(kind = 16), pointer :: temp(:,:)
        temp => o
        o => oo
        oo => temp
    end subroutine swap_arrays

    pure function states_to_index(states) result(index)
        ! states is a 1D array of length l, with values in [1, q]
        integer(kind = 1), intent(in) :: states(l)
        integer(kind = 1) :: i
        integer :: index

        index = 0
        do concurrent (i = 1:l)
            index = index + (states(i) - 1) * q ** (i-1)
        end do
        index = index + 1
    end function states_to_index

    pure function index_to_states(index) result(states)
        integer, intent(in) :: index
        integer(kind = 1) :: states(l), i

        do concurrent (i = 1:l)
            states(i) = mod((index - 1) / q ** (i-1), q) + 1
        end do
    end function index_to_states

    pure function intra_bonds(index) result(b)
        integer, intent(in) :: index
        integer :: b
        integer(kind = 1) :: states(l), i
        states = index_to_states(index)

        b = 0
        do concurrent (i = 1:l)
            b = b + interaction(states(i), states(mod(i, l) + 1))
        end do
    end function intra_bonds

end module potts_coefficients