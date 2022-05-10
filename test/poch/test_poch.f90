! test_poch.f90 --
!     Simple tests for the functions that evaluate the Pochhammer symbol and variant
!
program test_poch
    use fullerton_poch

    implicit none

    real    :: a, x
    integer :: i, j

    open( 10, file = 'test_poch.out' )

    write( 10, '(a)' )     'Regular range - pochhammer and pochhammer-1'
    write( 10, '(10a16)' ) 'a', 'x', 'poch', 'poch1'

    do i = 1,20
        a = (i-1) * 0.4
        do j= 1,30
            x = (j-1) * 0.3

            write( 10, '(10g16.8)' ) a, x, poch(a,x), poch1(a,x)
        enddo
    enddo

    write( 10, '(10a16)' ) 'a', 'x', 'poch1'

    do i = 1,10
        a = (i-1) * 0.4
        do j= 1,30
            x = (j-15) * 0.003

            write( 10, '(10g16.8)' ) a, x, poch(a,x), poch1(a,x)
        enddo
    enddo

    ! No provision for boundaries between the various ranges

end program test_poch
