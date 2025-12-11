program test
    use potts_em
    implicit none
    
    call calculate_coefficients('periodic')
    call print_coefficients()
end program test