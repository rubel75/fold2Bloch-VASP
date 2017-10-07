!!! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!!! vcross
!!!
!!! This subroutine is a part of a WaveTrans package and computes 
!!! vector cross-product
!!!
!!! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

SUBROUTINE vcross(a, & ! -> args out 
            b,c) ! <- args in 
IMPLICIT REAL*8(a-h,o-z)
DIMENSION a(3),b(3),c(3)
  
a(1)=b(2)*c(3)-b(3)*c(2)
a(2)=b(3)*c(1)-b(1)*c(3)
a(3)=b(1)*c(2)-b(2)*c(1)
RETURN

END SUBROUTINE vcross