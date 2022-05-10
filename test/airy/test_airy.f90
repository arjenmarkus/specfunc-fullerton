! test_airy.f90 --
!     Simple tests for the Airy functions
!
program test_airy
    use fullerton_airy

    implicit none

    real    :: dx, x, ya, yb, yad, ybd
    integer :: i

    open( 10, file = 'test_airy.out')

    dx = 0.1

    write( 10, '(a)' ) 'Regular range'
    do i = -100,100
        x   = i * dx
        ya  = airy_ai(x)
        yad = airy_aiprime(x)
        yb  = airy_bi(x)
        ybd = airy_biprime(x)
        write( 10, '(5g16.8)' ) x, ya, yad, yb, ybd
    enddo

    ! To be done:
    ! Ranges: -1. 1, 2, 4
    !write( 10, '(a)' ) 'Breakpoint x = sqrt(2.0*r1mach(1)) - besk1'
    !
    !xc = sqrt(4.0*r1mach(1))
    !xm = xc  * (1.0 - epsilon(1.0))  ; ym = besk1(xm)
    !x  = xc                          ; y  = besk1(x)
    !xp = xc  * (1.0 + epsilon(1.0))  ; yp = besk1(xp)
    !write( 10, '(7g16.8)' ) x, y, ym, yp, y-ym, yp-y, abs(y-ym)-abs(yp-y)

end program test_airy
