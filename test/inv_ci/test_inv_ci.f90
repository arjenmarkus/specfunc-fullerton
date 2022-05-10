! test_inv_ci.f90 --
!     Simple tests for the inverse cosine and cosh integrals
!
!     NOTE:
!     The breakpoint in inverse_chi reveals a problem - the function is not continuous there?
!
program test_inv_ci
    use fullerton_aux
    use fullerton_inv_ci

    implicit none

    real    :: dx, x, yinverse_ci, yinverse_chi, xc, xm, xp, ym, y, yp
    integer :: i

    open( 10, file = 'test_inv_ci.out')

    dx = 0.1

    write( 10, '(a)' ) 'Regular range'
    write( 10, '(10a16)' ) 'x', 'inverse_ci', 'inverse_chi'
    do i = 1,100
        x     = dx * (i - 0.5)
        yinverse_ci  = inverse_ci(x)
        yinverse_chi = inverse_chi(x)
        write( 10, '(10g16.8)' ) x, yinverse_ci, yinverse_chi
    enddo

    ! Ranges: 3.0
    write( 10, '(a)' ) 'Breakpoint x = 3.0 - inverse_chi'

    xc = 3.0
    xm = xc  * (1.0 - epsilon(1.0))  ; ym = inverse_chi(xm)
    x  = xc                          ; y  = inverse_chi(x)
    xp = xc  * (1.0 + epsilon(1.0))  ; yp = inverse_chi(xp)
    write( 10, '(7g16.8)' ) x, y, ym, yp, y-ym, yp-y, abs(y-ym)-abs(yp-y)


    write( 10, '(a)' ) 'Breakpoint x = 4.0 - inverse_ci'

    xc = 4.0
    xm = xc  * (1.0 - epsilon(1.0))  ; ym = inverse_ci(xm)
    x  = xc                          ; y  = inverse_ci(x)
    xp = xc  * (1.0 + epsilon(1.0))  ; yp = inverse_ci(xp)
    write( 10, '(7g16.8)' ) x, y, ym, yp, y-ym, yp-y, abs(y-ym)-abs(yp-y)

end program test_inv_ci
