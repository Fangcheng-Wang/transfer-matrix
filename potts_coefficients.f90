module potts_coefficients
    implicit none

    private
    public :: initialize

    ! Fortran requires array dimensions to be determined at compile time
    integer, parameter :: l = 4, m = 4, q = 2
    integer :: n_bonds = l * m * 2
    integer(kind = 16), allocatable :: o(:,:)

contains

    subroutine initialize()
        allocate(o(n_bonds, q ** l))
        o(:, :) = 0
    end subroutine initialize

    function states_to_1d(states) result(index)
        ! states is a 1D array of length l, with values in [0, q-1]
        integer, intent(in) :: states(l)
        integer :: index, i
        index = states(1)
        do concurrent (i = 2:l)
            index = index * q + states(i)
        end do
        index = index + 1
    end function states_to_1d


end module potts_coefficients