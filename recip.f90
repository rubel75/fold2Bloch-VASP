!!! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!!! recip
!!!
!!! This subroutine is a part of a WaveTrans package and computes 
!!! reciprocal properties.
!!!
!!! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

SUBROUTINE recip(ecut, a1, a2, a3, spinor, & ! <- args in 
            b1, b2, b3, nb1max, nb2max, nb3max, npmax) ! -> args out 
IMPLICIT NONE
REAL(kind=8), INTENT(in) :: ecut ! cut-off energy [eV]
REAL(kind=8), DIMENSION(1:3), &
    INTENT(in) :: a1, a2, a3 ! real space lattice vectors [Ang]
LOGICAL, INTENT(in) :: spinor ! WF are spinors?
REAL(kind=8), DIMENSION(1:3), INTENT(out) :: b1, b2, &
    b3 ! reciprocal lattice vectors [rad/Ang]
INTEGER, INTENT(out) :: nb1max, nb2max, nb3max ! max number of G vect.
INTEGER, INTENT(out) :: npmax
REAL(kind=8) :: Vcell ! volume of the unit cell
REAL(kind=8) :: b1mag, b2mag, b3mag ! magnitude of reciprocal 
                ! lattice vectors [rad/Ang]
REAL(kind=8) :: phi12, phi13, phi23 ! angles [rad]
REAL(kind=8) :: phi123, sinphi123
REAL(kind=8) :: vmag
REAL(kind=8), DIMENSION(1:3)  :: vtmp ! temp. variable
INTEGER :: nb1maxA, nb2maxA, nb3maxA, &
            nb1maxB, nb2maxB, nb3maxB, &
            nb1maxC, nb2maxC, nb3maxC
INTEGER :: npmaxA, npmaxB, npmaxC
INTEGER :: j ! dimensionality counter

!!$*   constant 'c' below is 2m/hbar**2 in units of 1/eV Ang^2 (value is
!!$*   adjusted in final decimal places to agree with VASP value; program
!!$*   checks for discrepancy of any results between this and VASP values)

REAL(kind=8), PARAMETER :: c=0.262465831
REAL(kind=8), PARAMETER :: pi=4.*atan(1.)

write(6,*) ' '
write(6,*) 'Computing reciprocal properties...'
call vcross(vtmp,a2,a3)
Vcell=a1(1)*vtmp(1)+a1(2)*vtmp(2)+a1(3)*vtmp(3)
write(6,*) 'volume unit cell [Ang3] =',sngl(Vcell)
call vcross(b1,a2,a3)
call vcross(b2,a3,a1)
call vcross(b3,a1,a2)
do j=1,3
   b1(j)=2.*pi*b1(j)/Vcell
   b2(j)=2.*pi*b2(j)/Vcell
   b3(j)=2.*pi*b3(j)/Vcell
enddo
b1mag=dsqrt(b1(1)**2+b1(2)**2+b1(3)**2)
b2mag=dsqrt(b2(1)**2+b2(2)**2+b2(3)**2)
b3mag=dsqrt(b3(1)**2+b3(2)**2+b3(3)**2)
write(6,*) 'reciprocal lattice vectors [rad/Ang]:'
write(6,*) 'b1 =',(sngl(b1(j)),j=1,3)
write(6,*) 'b2 =',(sngl(b2(j)),j=1,3)
write(6,*) 'b3 =',(sngl(b3(j)),j=1,3)
write(6,*) 'reciprocal lattice vector magnitudes [rad/Ang]:'
write(6,*) sngl(b1mag),sngl(b2mag),sngl(b3mag)

phi12=acos((b1(1)*b2(1)+b1(2)*b2(2)+b1(3)*b2(3))/(b1mag*b2mag))
call vcross(vtmp,b1,b2)
vmag=dsqrt(vtmp(1)**2+vtmp(2)**2+vtmp(3)**2)
sinphi123=(b3(1)*vtmp(1)+b3(2)*vtmp(2)+b3(3)*vtmp(3))/(vmag*b3mag)
nb1maxA=NINT(dsqrt(ecut*c)/(b1mag*abs(sin(phi12))))+1
nb2maxA=NINT(dsqrt(ecut*c)/(b2mag*abs(sin(phi12))))+1
nb3maxA=NINT(dsqrt(ecut*c)/(b3mag*abs(sinphi123)))+1
npmaxA=NINT(4.*pi*nb1maxA*nb2maxA*nb3maxA/3.)
      
phi13=acos((b1(1)*b3(1)+b1(2)*b3(2)+b1(3)*b3(3))/(b1mag*b3mag))
call vcross(vtmp,b1,b3)
vmag=dsqrt(vtmp(1)**2+vtmp(2)**2+vtmp(3)**2)
sinphi123=(b2(1)*vtmp(1)+b2(2)*vtmp(2)+b2(3)*vtmp(3))/(vmag*b2mag)
phi123=abs(asin(sinphi123))
nb1maxB=NINT(dsqrt(ecut*c)/(b1mag*abs(sin(phi13))))+1
nb2maxB=NINT(dsqrt(ecut*c)/(b2mag*abs(sinphi123)))+1
nb3maxB=NINT(dsqrt(ecut*c)/(b3mag*abs(sin(phi13))))+1
npmaxB=NINT(4.*pi*nb1maxB*nb2maxB*nb3maxB/3.)
      
phi23=acos((b2(1)*b3(1)+b2(2)*b3(2)+b2(3)*b3(3))/(b2mag*b3mag))
call vcross(vtmp,b2,b3)
vmag=dsqrt(vtmp(1)**2+vtmp(2)**2+vtmp(3)**2)
sinphi123=(b1(1)*vtmp(1)+b1(2)*vtmp(2)+b1(3)*vtmp(3))/(vmag*b1mag)
phi123=abs(asin(sinphi123))
nb1maxC=NINT(dsqrt(ecut*c)/(b1mag*abs(sinphi123)))+1
nb2maxC=NINT(dsqrt(ecut*c)/(b2mag*abs(sin(phi23))))+1
nb3maxC=NINT(dsqrt(ecut*c)/(b3mag*abs(sin(phi23))))+1 
npmaxC=NINT(4.*pi*nb1maxC*nb2maxC*nb3maxC/3.)

nb1max=max0(nb1maxA,nb1maxB,nb1maxC)
nb2max=max0(nb2maxA,nb2maxB,nb2maxC)
nb3max=max0(nb3maxA,nb3maxB,nb3maxC)

npmax=min0(npmaxA,npmaxB,npmaxC)
IF (spinor) THEN
    npmax=2*npmax
ENDIF

write(6,*) 'max. no. G values; 1,2,3 =',nb1max,nb2max,nb3max
write(6,*) 'estimated max. no. plane waves =',npmax
write(6,*) 'Computing reciprocal properties: done!'

END SUBROUTINE recip
