program test
    use potts_coefficients
    implicit none
    
    call calculate_coefficients('open')
    call print_coefficients()
end program test