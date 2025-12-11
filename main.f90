program test
    use potts_coefficients
    implicit none
    
    integer :: index

    index = states_to_index([0, 0, 0, 0])
    print *, index
    print *, index_to_states(index)
end program test