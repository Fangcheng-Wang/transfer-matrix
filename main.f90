program test
    use potts2_em
    implicit none
    
    call calculate_coefficients('periodic')
    call print_coefficients()
end program test