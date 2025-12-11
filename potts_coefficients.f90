module potts_coefficients
    implicit none

    private
    public :: initialize, states_to_index, index_to_states

    ! Fortran requires array dimensions to be determined at compile time
    integer, parameter :: l = 4, m = 4, q = 2
    integer :: n_bonds = l * m * 2
    integer(kind = 16), allocatable :: o(:,:)
    integer(kind = 16), allocatable :: oo(:,:)

contains

    subroutine initialize()
        allocate(o(n_bonds, q ** l))
        allocate(oo(n_bonds, q ** l))
        o(:, :) = 0
        oo(:, :) = 0
    end subroutine initialize

    function states_to_index(states) result(index)
        ! states is a 1D array of length l, with values in [0, q-1]
        integer, intent(in) :: states(l)
        integer :: index, i
        index = 0
        do concurrent (i = 1:l)
            index = index + states(i) * q ** (i-1)
        end do
        index = index + 1
    end function states_to_index

    function index_to_states(index) result(states)
        integer, intent(in) :: index
        integer :: states(l), i
        do concurrent (i = 1:l)
            states(i) = mod((index - 1) / q ** (i-1), q)
        end do
    end function index_to_states

end module potts_coefficients