! test_i0.f90 --
!     Simple tests for I0
!
program test_i0
    use fullerton_bessel

    implicit none

    real    :: xc, x, dx, xm, xp, y, ym, yp
    integer :: i

    open( 10, file = 'test_i0.out')

    dx = 0.1

    write( 10, '(a)' ) 'Regular range'
    do i = 1,101
        x = 0.00001 + (i-1) * dx
        write( 10, '(2g16.8)' ) x, bessel_i0(x)
    enddo


    write( 10, '(a)' ) 'Breakpoint x = 3.0 - bessel_i0/bessel_i0e'

    xc = 3.0
    xm = xc  * (1.0 - epsilon(1.0))  ; ym = bessel_i0(xm)
    x  = xc                          ; y  = bessel_i0(x)
    xp = xc  * (1.0 + epsilon(1.0))  ; yp = bessel_i0(xp)
    write( 10, '(7g16.8)' ) x, y, ym, yp, y-ym, yp-y, abs(y-ym)-abs(yp-y)

    write( 10, '(a)' ) 'Breakpoint x = 8.0 - bessel_i0e'

    xc = 8.0
    xm = xc  * (1.0 - epsilon(1.0))  ; ym = bessel_i0(xm)
    x  = xc                          ; y  = bessel_i0(x)
    xp = xc  * (1.0 + epsilon(1.0))  ; yp = bessel_i0(xp)
    write( 10, '(7g16.8)' ) x, y, ym, yp, y-ym, yp-y, abs(y-ym)-abs(yp-y)

end program test_i0
