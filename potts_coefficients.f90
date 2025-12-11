module potts_coefficients
    implicit none

    private
    public :: initialize, print_coefficients

    ! Fortran requires array dimensions to be determined at compile time
    integer(kind = 1), parameter :: l = 4, m = 4, q = 2
    integer :: n_bonds = l * m * 2
    integer(kind = 1), allocatable :: interaction(:,:)
    ! o(1+energy, index) is the number of configurations,
    ! where energy is equal to the number of bonds with the *SAME* state (antiferromagnetic)
    integer(kind = 16), allocatable :: o(:,:)
    integer(kind = 16), allocatable :: oo(:,:)

contains

    subroutine initialize()
        integer :: index, i

        allocate(interaction(q, q))
        interaction(:, :) = 0
        do concurrent (i = 1:q)
            interaction(i, i) = 1
        end do

        allocate(o(n_bonds, q ** l))
        allocate(oo(n_bonds, q ** l))
        o(:, :) = 0
        do concurrent (index = 1:q**l)
            o(1+layer_energy(index), index) = 1
        end do
    end subroutine initialize

    subroutine print_coefficients()
        integer :: i
        integer(kind = 16), allocatable :: total(:)
        allocate(total(n_bonds))
        total = sum(o, 2)

        print *, '     energy                                    count'
        do i = 1, n_bonds
            print *, i, total(i)
        end do
        deallocate(total)
    end subroutine print_coefficients

    subroutine swap_arrays()
        integer(kind = 16), allocatable :: temp(:,:)
        call move_alloc(from=o, to=temp)
        call move_alloc(from=oo, to=o)
        call move_alloc(from=temp, to=oo)
    end subroutine swap_arrays

    function states_to_index(states) result(index)
        ! states is a 1D array of length l, with values in [1, q]
        integer, intent(in) :: states(l)
        integer :: index, i
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

    pure function layer_energy(index) result(energy)
        integer, intent(in) :: index
        integer :: energy, states(l), i
        states = index_to_states(index)

        energy = 0
        do concurrent (i = 1:l)
            energy = energy + interaction(states(i), states(mod(i, l) + 1))
        end do
    end function layer_energy

end module potts_coefficients