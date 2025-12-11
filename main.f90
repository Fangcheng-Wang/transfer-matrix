program test
    use potts_coefficients
    implicit none
    
    call initialize()
    call iterate_layers()
    call print_coefficients()
end program test