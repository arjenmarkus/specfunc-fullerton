! test_complex.f90 --
!     Simple test program for the complex functions provided by Fullerton's library
!
program test_complex
    use fullerton_complex

    implicit none

    complex :: z, a, b
    integer :: i, j

    open( 10, file = 'test_complex_fullerton.out' )

    write( 10, '(30a16)' ) 'z (real)', 'z (imag)',  'fgamma (real)', 'fgamma (imag)',  'reciprocal_gamma', &
                           ' gamr (imag)', 'digamma (real)', 'digamma (imag)',  'log_gamma (real)', 'log_gamma (imag)'
    do j = -20,20
        do i = -20,20
            z = cmplx( 0.25*i, 0.35*j )

            write( 10, '(30g16.8)' ) z, fgamma(z), reciprocal_gamma(z), digamma(z), log_gamma(z)
        enddo
    enddo

    write( 10, '(30a16)' ) 'z (real)', 'z (imag)', 'lnrel (real)', 'lnrel (imag)', 'exprl (real)', &
                           'exprl (imag)', 'cbrt (real)', 'cbrt (imag)', 'log10 (real)', 'log10 (imag)'
    do j = -20,20
        do i = -20,20
            z = cmplx( 0.25*i, 0.35*j )

            write( 10, '(30g16.8)' ) z, lnrel(z), exprl(z), cbrt(z), log10(z)
        enddo
    enddo

    write( 10, '(30a16)' ) 'a (real)', 'a (imag)', 'b (real)', 'b (imag)', 'beta (real)', 'beta (imag)', &
                           'log_beta (real)', 'log_beta (imag)', 'atan2 (real)', 'atan2 (imag)'
    do j = 1,20
        do i = 1,20
            a = cmplx( 0.10*i, 0.10*j )
            b = cmplx( 0.15*j, 0.15*i )

            write( 10, '(30g16.8)' ) a, b, beta(a,b), log_beta(a,b), atan2(a,b)
        enddo
    enddo
end program test_complex
