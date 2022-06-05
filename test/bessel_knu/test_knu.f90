! test_knu.f90 --
!     Simple tests for Knu
!
program test_knu
    use fullerton_bessel

    implicit none

    real    :: bk(5), bkm(5), bkp(5)
    real    :: xnu, x, dx, xc, xm, xp
    integer :: i

    open( 10, file = 'test_knu.out')

    dx = 0.1

    write( 10, '(a)'   ) 'Regular range - nu = 0.6'
    write( 10, '(10a)' ) 'x               ', 'nu = 0.6       ', &
                         'nu = 1.6        ', 'nu = 2.6       ', &
                         'nu = 3.6        ', 'nu = 4.6       '
    do i = 1,101
        x = 0.00001 + (i-1) * dx
        xnu = 0.6
        call bessel_knu( xnu, x, bk )
        write( 10, '(20g16.8)' ) x, bk
    enddo

    write( 10, '(a)' ) 'Breakpoint x = 2.0'

    xc = 2.0
    xm = xc  * (1.0 - epsilon(1.0))  ; call bessel_knu( xnu, xm, bkm )
    x  = xc                          ; call bessel_knu( xnu, x,  bk  )
    xp = xc  * (1.0 + epsilon(1.0))  ; call bessel_knu( xnu, xp, bkp )
    do i = 1,size(bk)
        write( 10, '(7g16.8)' ) x, bkm(i), bk(i), bkp(i), abs(bkm(i)-bk(i)), abs(bkp(i)-bk(i))
    enddo

    write( 10, '(a)'   ) 'Regular range - nu = 0.4'
    write( 10, '(10a)' ) 'x               ', 'nu = 0.4       ', &
                         'nu = 1.4        ', 'nu = 2.4       ', &
                         'nu = 3.4        ', 'nu = 4.4       '
    do i = 1,101
        x = 0.00001 + (i-1) * dx
        xnu = 0.4
        call bessel_knu( xnu, x, bk )
        write( 10, '(20g16.8)' ) x, bk
    enddo

    write( 10, '(a)' ) 'Breakpoint x = 2.0'

    xc = 2.0
    xm = xc  * (1.0 - epsilon(1.0))  ; call bessel_knu( xnu, xm, bkm )
    x  = xc                          ; call bessel_knu( xnu, x,  bk  )
    xp = xc  * (1.0 + epsilon(1.0))  ; call bessel_knu( xnu, xp, bkp )
    do i = 1,size(bk)
        write( 10, '(7g16.8)' ) x, bkm(i), bk(i), bkp(i), abs(bkm(i)-bk(i)), abs(bkp(i)-bk(i))
    enddo

end program test_knu
