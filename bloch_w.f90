!!! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!!! bloch_w
!!!
!!! Compute Bloch spectral weight by sorting PW coefficients 
!!! into appropriate groups.
!!!
!!! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
SUBROUTINE bloch_w(folds, nplane, igall, & ! <- args in 
                coeff, nnewk, & ! <- args in 
                w) ! -> args out
IMPLICIT NONE
INTEGER, INTENT(in) :: folds(3), nplane
INTEGER, INTENT(in) :: igall(3,nplane)
COMPLEX*8, INTENT(in) :: coeff(nplane)
INTEGER, INTENT(in) :: nnewk ! number of new k-points after unfolding
REAL(kind=8), INTENT(out) ::  w(nnewk) ! Bloch spectral weights
COMPLEX*8 :: TGroupC(folds(1),folds(2),folds(3),nplane)
REAL(kind=8) :: Sums(nnewk)
INTEGER :: counter(folds(1),folds(2),folds(3))
REAL(kind=8) :: sumtot
INTEGER :: remainder_x, remainder_y, remainder_z, j, k, l, p, el
INTEGER :: FX, FY, FZ

FX = folds(1); FY = folds(2); FZ = folds(3)
 
!! Initiates the counter and TGroupC elements at 0

DO j=1,FX
    DO k=1,FY
        DO l=1,FZ
            counter(j,k,l)=0
            DO p=1,nplane
                TGroupC(j,k,l,p)=0.0
            ENDDO
        ENDDO
    ENDDO
ENDDO

!! Sorts the PW coefficients

do j=1, nplane
    remainder_x=MODULO(igall(1,j), FX)
    remainder_y=MODULO(igall(2,j), FY)
    remainder_z=MODULO(igall(3,j), FZ)
    counter(remainder_x+1, remainder_y+1, remainder_z+1) = &
        counter(remainder_x+1, remainder_y+1, remainder_z+1)+1
    TGroupC(remainder_x+1, remainder_y+1, remainder_z+1, &
        counter(remainder_x+1, remainder_y+1, remainder_z+1))=coeff(j)
enddo

!! Sums the squares  of all coefficients per group
 
el=1
do j=1, FX
    do k=1, FY
        do l=1, FZ
            if (counter(j, k, l).gt.0) then
                do p=1, counter(j, k, l)
                    TGroupC(j, k, l,p)= &
                        TGroupC(j, k, l,p)*CONJG(TGroupC(j, k, l,p))
                enddo
                Sums(el)=SUM(TGroupC(j, k, l,1:counter(j, k, l)))
                el=el+1
            else 
                Sums(el)=0.0
                el=el+1
            endif
        enddo
    enddo
enddo
sumtot=SUM(Sums)
do j=1, nnewk
    w(j)=Sums(j)/sumtot
enddo
 
END SUBROUTINE bloch_w
