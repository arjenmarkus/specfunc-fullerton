! test_ei.f90 --
!     Simple tests for the ei function and related functions
!
program test_ei
    use fullerton_ei

    implicit none

    real    :: dx, x, yei, ye1, ysi, yci, yshi, ychi, yli, xc, xm, xp, ym, y, yp
    integer :: i

    open( 10, file = 'test_ei.out')

    dx = 0.1

    write( 10, '(a)' ) 'Regular range'
    write( 10, '(10a16)' ) 'x', 'ei', 'e1', 'si', 'ci', 'shi', 'chi', 'li'
    do i = 1,100
        x   = dx * (i - 0.5)
        yei  = integral_ei(x)
        ye1  = integral_e1(x)
        ysi  = integral_si(x)
        yci  = integral_ci(x)
        yshi = integral_shi(x)
        ychi = integral_chi(x)
        yli  = integral_li(x)
        write( 10, '(10g16.8)' ) x, yei, ye1, ysi, yci, yshi, ychi, yli
    enddo

    ! Ranges: 4.0
    write( 10, '(a)' ) 'Breakpoint x = 4.0 - si'

    xc = 4.0
    xm = xc  * (1.0 - epsilon(1.0))  ; ym = integral_si(xm)
    x  = xc                          ; y  = integral_si(x)
    xp = xc  * (1.0 + epsilon(1.0))  ; yp = integral_si(xp)
    write( 10, '(7g16.8)' ) x, y, ym, yp, y-ym, yp-y, abs(y-ym)-abs(yp-y)


    write( 10, '(a)' ) 'Breakpoint x = 4.0 - ci'

    xc = 4.0
    xm = xc  * (1.0 - epsilon(1.0))  ; ym = integral_ci(xm)
    x  = xc                          ; y  = integral_ci(x)
    xp = xc  * (1.0 + epsilon(1.0))  ; yp = integral_ci(xp)
    write( 10, '(7g16.8)' ) x, y, ym, yp, y-ym, yp-y, abs(y-ym)-abs(yp-y)


    write( 10, '(a)' ) 'Breakpoint x = 0.375 - shi'

    xc = 0.375
    xm = xc  * (1.0 - epsilon(1.0))  ; ym = integral_shi(xm)
    x  = xc                          ; y  = integral_shi(x)
    xp = xc  * (1.0 + epsilon(1.0))  ; yp = integral_shi(xp)
    write( 10, '(7g16.8)' ) x, y, ym, yp, y-ym, yp-y, abs(y-ym)-abs(yp-y)


    write( 10, '(a)' ) 'Breakpoint x = -10.0 - e1'

    xc = -10.0
    xm = xc  * (1.0 - epsilon(1.0))  ; ym = integral_e1(xm)
    x  = xc                          ; y  = integral_e1(x)
    xp = xc  * (1.0 + epsilon(1.0))  ; yp = integral_e1(xp)
    write( 10, '(7g16.8)' ) x, y, ym, yp, y-ym, yp-y, abs(y-ym)-abs(yp-y)

    write( 10, '(a)' ) 'Breakpoint x = -4.0 - e1'

    xc = -4.0
    xm = xc  * (1.0 - epsilon(1.0))  ; ym = integral_e1(xm)
    x  = xc                          ; y  = integral_e1(x)
    xp = xc  * (1.0 + epsilon(1.0))  ; yp = integral_e1(xp)
    write( 10, '(7g16.8)' ) x, y, ym, yp, y-ym, yp-y, abs(y-ym)-abs(yp-y)

    write( 10, '(a)' ) 'Breakpoint x = -1.0 - e1'

    xc = -1.0
    xm = xc  * (1.0 - epsilon(1.0))  ; ym = integral_e1(xm)
    x  = xc                          ; y  = integral_e1(x)
    xp = xc  * (1.0 + epsilon(1.0))  ; yp = integral_e1(xp)
    write( 10, '(7g16.8)' ) x, y, ym, yp, y-ym, yp-y, abs(y-ym)-abs(yp-y)

    write( 10, '(a)' ) 'Breakpoint x = 1.0 - e1'

    xc = 1.0
    xm = xc  * (1.0 - epsilon(1.0))  ; ym = integral_e1(xm)
    x  = xc                          ; y  = integral_e1(x)
    xp = xc  * (1.0 + epsilon(1.0))  ; yp = integral_e1(xp)
    write( 10, '(7g16.8)' ) x, y, ym, yp, y-ym, yp-y, abs(y-ym)-abs(yp-y)

    write( 10, '(a)' ) 'Breakpoint x = 4.0 - e1'

    xc = 4.0
    xm = xc  * (1.0 - epsilon(1.0))  ; ym = integral_e1(xm)
    x  = xc                          ; y  = integral_e1(x)
    xp = xc  * (1.0 + epsilon(1.0))  ; yp = integral_e1(xp)
    write( 10, '(7g16.8)' ) x, y, ym, yp, y-ym, yp-y, abs(y-ym)-abs(yp-y)

end program test_ei
