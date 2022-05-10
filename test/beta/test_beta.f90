! test_beta.f90 --
!     Simple tests for the beta function and variants
!
program test_beta
    use fullerton_beta

    implicit none

    real    :: a, b, dx, x, yb, ybi, ylogb, xc, xm, xp, ym, y, yp
    integer :: i, j, k

    open( 10, file = 'test_beta.out')

    dx = 0.2

    write( 10, '(a)' ) 'Regular range'
    do j = 1,10
        do i = 1,100
            a     = i * dx
            b     = dx + 10 * (j-1) * dx
            yb    = beta(a,b)
            ylogb = log_beta(a,b)
            write( 10, '(5g16.8)' ) a, b, yb, ylogb
        enddo
    enddo

    ! Ranges: 10.0
    write( 10, '(a)' ) 'Breakpoint a/b = 10.0 - albeta'

    a  =  1.0
    xc = 10.0
    xm = xc  * (1.0 - epsilon(1.0))  ; ym = log_beta(a,xm)
    x  = xc                          ; y  = log_beta(a,x)
    xp = xc  * (1.0 + epsilon(1.0))  ; yp = log_beta(a,xp)
    write( 10, '(7g16.8)' ) x, y, ym, yp, y-ym, yp-y, abs(y-ym)-abs(yp-y)


    write( 10, '(a)' ) 'Incomplete beta function'

    do j = 1,5
        do i = 1,5
            a = 0.1  + i-1
            b = 0.15 + j-1
            do k = 1,21
                x = (k-1) * 0.05
                ybi = beta_inc(x,a,b)
                write( 10, '(7g16.8)' ) x, a, b, ybi
            enddo
       enddo
   enddo

end program test_beta
