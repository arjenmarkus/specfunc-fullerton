! test_misc_exp.f90 --
!     Simple tests for the functions related to the exponential function
!
!     Note:
!     Dawson's function/integral F has a derivative 1 - 2x F(x)
!
!     Spence's function/integral may be incorrectly implemented - it should be a monotonously increasing function.
!     If you look at the value for negative x, then it does look like minus the function Li2 shown on Mathworld.
!     I am not sure what to make of it. Perhaps use the definition in a naïve way? Note, however, that the complex
!     extension has a branchcut at x = 1.
!
!     F(x) = - Li2(-x) = sum k=1 to inf x**k / k**2
!
program test_misc_exp
    use fullerton_misc

    implicit none

    real    :: dx, x, yalnrel, yexprel, ydaws, yspenc, xc, xm, xp, ym, y, yp
    integer :: i

    open( 10, file = 'test_misc_exp.out')

    dx = 0.1

    write( 10, '(a)' ) 'Regular range - complete misc_exp'
    write( 10, '(10a16)' ) 'x', 'alnrel', 'exprel', 'integral_dawson', 'integral_spence'
    do i = 1,100
        x      = dx * (i - 0.5)
        yalnrel = lnrel(x)
        yexprel = exprel(x)
        ydaws   = integral_dawson(x)
        yspenc  = integral_spence(x)
        write( 10, '(10g16.8)' ) x, yalnrel, yexprel, ydaws, yspenc
    enddo


    ! Ranges: 0.375
    write( 10, '(a)' ) 'Breakpoint x = 0.375 - lnrel'

    xc = 0.375
    xm = xc  * (1.0 - epsilon(1.0))  ; ym = lnrel(xm)
    x  = xc                          ; y  = lnrel(x)
    xp = xc  * (1.0 + epsilon(1.0))  ; yp = lnrel(xp)
    write( 10, '(7g16.8)' ) x, y, ym, yp, y-ym, yp-y, abs(y-ym)-abs(yp-y)


    ! Ranges: 0.5
    write( 10, '(a)' ) 'Breakpoint x = 0.5 - exprel'

    xc = 0.5
    xm = xc  * (1.0 - epsilon(1.0))  ; ym = exprel(xm)
    x  = xc                          ; y  = exprel(x)
    xp = xc  * (1.0 + epsilon(1.0))  ; yp = exprel(xp)
    write( 10, '(7g16.8)' ) x, y, ym, yp, y-ym, yp-y, abs(y-ym)-abs(yp-y)


    ! Ranges: 1.0, 4.0
    write( 10, '(a)' ) 'Breakpoint x = 1.0 - integral_dawson'

    xc = 1.0
    xm = xc  * (1.0 - epsilon(1.0))  ; ym = integral_dawson(xm)
    x  = xc                          ; y  = integral_dawson(x)
    xp = xc  * (1.0 + epsilon(1.0))  ; yp = integral_dawson(xp)
    write( 10, '(7g16.8)' ) x, y, ym, yp, y-ym, yp-y, abs(y-ym)-abs(yp-y)

    write( 10, '(a)' ) 'Breakpoint x = 4.0 - integral_dawson'

    xc = 4.0
    xm = xc  * (1.0 - epsilon(1.0))  ; ym = integral_dawson(xm)
    x  = xc                          ; y  = integral_dawson(x)
    xp = xc  * (1.0 + epsilon(1.0))  ; yp = integral_dawson(xp)
    write( 10, '(7g16.8)' ) x, y, ym, yp, y-ym, yp-y, abs(y-ym)-abs(yp-y)


    ! Ranges: -1.0, 0.0, 0.5, 1.0, 2.0
    write( 10, '(a)' ) 'Breakpoint x = -1.0 - integral_spence'

    xc = -1.0
    xm = xc  * (1.0 - epsilon(1.0))  ; ym = integral_spence(xm)
    x  = xc                          ; y  = integral_spence(x)
    xp = xc  * (1.0 + epsilon(1.0))  ; yp = integral_spence(xp)
    write( 10, '(7g16.8)' ) x, y, ym, yp, y-ym, yp-y, abs(y-ym)-abs(yp-y)

    write( 10, '(a)' ) 'Breakpoint x = 0.0 - integral_spenc'

    xc = 0.0
    xm = xc  * (1.0 - epsilon(1.0))  ; ym = integral_spence(xm)
    x  = xc                          ; y  = integral_spence(x)
    xp = xc  * (1.0 + epsilon(1.0))  ; yp = integral_spence(xp)
    write( 10, '(7g16.8)' ) x, y, ym, yp, y-ym, yp-y, abs(y-ym)-abs(yp-y)

    write( 10, '(a)' ) 'Breakpoint x = 1.0 - integral_spenc'

    xc = 1.0
    xm = xc  * (1.0 - epsilon(1.0))  ; ym = integral_spence(xm)
    x  = xc                          ; y  = integral_spence(x)
    xp = xc  * (1.0 + epsilon(1.0))  ; yp = integral_spence(xp)
    write( 10, '(7g16.8)' ) x, y, ym, yp, y-ym, yp-y, abs(y-ym)-abs(yp-y)

    write( 10, '(a)' ) 'Breakpoint x = 2.0 - integral_spenc'

    xc = 2.0
    xm = xc  * (1.0 - epsilon(1.0))  ; ym = integral_spence(xm)
    x  = xc                          ; y  = integral_spence(x)
    xp = xc  * (1.0 + epsilon(1.0))  ; yp = integral_spence(xp)
    write( 10, '(7g16.8)' ) x, y, ym, yp, y-ym, yp-y, abs(y-ym)-abs(yp-y)

    open( 11, file = 'test_spenc.out')
    do i = 1,100
        x      = -dx * (i - 0.5)
        yspenc  = integral_spence(x)
        write( 11, '(10g16.8)' ) x, yspenc
    enddo

end program test_misc_exp
