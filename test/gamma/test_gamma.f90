! test_gamma.f90 --
!     Simple tests for the functions related to the gamma function
!
!     TODO: negative arguments
!     TODO: negative integer arguments
!     TODO: more careful consideration of breakpoints
!
!     Override intrinsic gamma function?
!
program test_gamma
    use fullerton_gamma

    implicit none

    real    :: a, dx, x, yalngam, yalngama, ysgngam, ygami, ygamic, ygamit, ygamr, ygamma, ypsi, xc, xm, xp, ym, y, yp
    integer :: i, j

    open( 10, file = 'test_gamma.out')

    dx = 0.1

    write( 10, '(a)' ) 'Regular range - complete gamma'
    write( 10, '(10a16)' ) 'x', 'gamma', 'alngam', 'algams', 'sgngam', 'gamr', 'psi'
    do i = 1,150
        x      = dx * (i - 0.5)
        ygamma  = fgamma(x)
        yalngam = log_gamma(x)
        call sign_log_gamma( x, yalngama, ysgngam )
        ygamr   = reciprocal_gamma(x)
        ypsi    = digamma(x)
        write( 10, '(10g16.8)' ) x, ygamma, yalngam, yalngama, ysgngam, ygamr, ypsi
    enddo

    write( 10, '(a)' ) 'Regular range - incomplete gamma'
    write( 10, '(10a16)' ) 'a', 'x', 'gami', 'gamic', 'gamit'
    do j = 1,10
        a = 0.1 + (j-1)
        do i = 1,20
             x      = dx * (i - 0.5)
            ygami  = inc_gamma(a,x)
            ygamic = inc_gamma_compl(a,x)
            ygamit = gamma_tricomi(a,x)
            write( 10, '(10g16.8)' ) a, x, ygami, ygamic, ygamit
        enddo
    enddo

    ! Ranges: 10.0
    write( 10, '(a)' ) 'Breakpoint x = 10.0 - gamma'

    xc = 10.0
    xm = xc  * (1.0 - epsilon(1.0))  ; ym = fgamma(xm)
    x  = xc                          ; y  = fgamma(x)
    xp = xc  * (1.0 + epsilon(1.0))  ; yp = fgamma(xp)
    write( 10, '(7g16.8)' ) x, y, ym, yp, y-ym, yp-y, abs(y-ym)-abs(yp-y)


    ! Incomplete!
    write( 10, '(a)' ) 'Breakpoint x = 1.0 - gamit'

    a  = 1.5
    xc = 1.0
    xm = xc  * (1.0 - epsilon(1.0))  ; ym = gamma_tricomi(a,xm)
    x  = xc                          ; y  = gamma_tricomi(a,x)
    xp = xc  * (1.0 + epsilon(1.0))  ; yp = gamma_tricomi(a,xp)
    write( 10, '(7g16.8)' ) x, y, ym, yp, y-ym, yp-y, abs(y-ym)-abs(yp-y)

    ! Incomplete!
    write( 10, '(a)' ) 'Breakpoint x = 1.0 - gamic'

    a  = 1.5
    xc = 1.0
    xm = xc  * (1.0 - epsilon(1.0))  ; ym = inc_gamma_compl(a,xm)
    x  = xc                          ; y  = inc_gamma_compl(a,x)
    xp = xc  * (1.0 + epsilon(1.0))  ; yp = inc_gamma_compl(a,xp)
    write( 10, '(7g16.8)' ) x, y, ym, yp, y-ym, yp-y, abs(y-ym)-abs(yp-y)

    write( 10, '(a)' ) 'Breakpoint x = 2.0 - psi'

    xc = 2.0
    xm = xc  * (1.0 - epsilon(1.0))  ; ym = digamma(xm)
    x  = xc                          ; y  = digamma(x)
    xp = xc  * (1.0 + epsilon(1.0))  ; yp = digamma(xp)
    write( 10, '(7g16.8)' ) x, y, ym, yp, y-ym, yp-y, abs(y-ym)-abs(yp-y)

end program test_gamma
