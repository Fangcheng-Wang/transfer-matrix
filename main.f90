program test
    use potts_emm
    implicit none
    
    call calculate_coefficients('o')
    call write_coefficients()
end program test