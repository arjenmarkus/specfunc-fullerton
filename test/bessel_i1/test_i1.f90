! test_i1.f90 --
!     Simple tests for I1
!
program test_i1
    use fullerton_bessel
    use fullerton_aux, only: r1mach

    implicit none

    real    :: xc, x, dx, xm, xp, y, ym, yp
    integer :: i

    open( 10, file = 'test_i1.out')

    dx = 0.1

    write( 10, '(a)' ) 'Regular range'
    do i = 1,101
        x = 0.00001 + (i-1) * dx
        write( 10, '(2g16.8)' ) x, bessel_i1(x)
    enddo

    write( 10, '(a)' ) 'Breakpoint x = sqrt(2.0*r1mach(1)) - bessel_i1/bessel_i1e'

    xc = sqrt(2.0*r1mach(1)) * (1.0 + 2.0 * epsilon(1.0))
    xm = xc  * (1.0 - epsilon(1.0))  ; ym = bessel_i1(xm)
    x  = xc                          ; y  = bessel_i1(x)
    xp = xc  * (1.0 + epsilon(1.0))  ; yp = bessel_i1(xp)
    write( 10, '(7g16.8)' ) x, y, ym, yp, y-ym, yp-y, abs(y-ym)-abs(yp-y)

    write( 10, '(a)' ) 'Breakpoint x = 3.0 - bessel_i1/bessel_i1e'

    xc = 3.0
    xm = xc  * (1.0 - epsilon(1.0))  ; ym = bessel_i1(xm)
    x  = xc                          ; y  = bessel_i1(x)
    xp = xc  * (1.0 + epsilon(1.0))  ; yp = bessel_i1(xp)
    write( 10, '(7g16.8)' ) x, y, ym, yp, y-ym, yp-y, abs(y-ym)-abs(yp-y)

    write( 10, '(a)' ) 'Breakpoint x = sqrt(8.0*r1mach(3) - bessel_i1/bessel_i1e'

    xc = sqrt(8.0*r1mach(3))
    xm = xc  * (1.0 - epsilon(1.0))  ; ym = bessel_i1(xm)
    x  = xc                          ; y  = bessel_i1(x)
    xp = xc  * (1.0 + epsilon(1.0))  ; yp = bessel_i1(xp)
    write( 10, '(7g16.8)' ) x, y, ym, yp, y-ym, yp-y, abs(y-ym)-abs(yp-y)

end program test_i1
