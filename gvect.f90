!!! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!!! gvect
!!!
!!! This subroutine calculates plane waves
!!!
!!! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

SUBROUTINE gvect(b1, b2, b3, nb1max, nb2max, nb3max, & ! <- args in
                nplane, npmax, wk, ecut, & ! <- args in
                spinor, igall) ! -> args out
IMPLICIT NONE
REAL(kind=8), DIMENSION(1:3), INTENT(in) :: b1, b2, &
    b3 ! reciprocal lattice vectors [rad/Ang]
INTEGER, INTENT(in) :: nb1max, nb2max, nb3max ! max number of G vect.
INTEGER, INTENT(in) :: npmax ! absolute max number of G vectors
INTEGER, INTENT(in) :: nplane ! no. of PWs
REAL(kind=8), INTENT(in) :: wk(3) ! k-point in units of b1, b2, b3
REAL(kind=8), INTENT(in) :: ecut ! cut-off energy [eV]
LOGICAL, INTENT(in) :: spinor ! WF spinor?
INTEGER, INTENT(out) :: igall(3,*) ! G vectors in units of b1, b2, b3
INTEGER :: ncnt ! PW counter
REAL(kind=8) :: sumkg(3) ! G-vector [rad/Ang]
REAL(kind=8) :: gtot ! magnitude |G| [rad/Ang]
REAL(kind=8) :: etot ! kinetic energy associated with |G|
INTEGER :: j, ig1, ig2, ig3, ig1p, ig2p, ig3p ! counters

!!$*   constant 'c' below is 2m/hbar**2 in units of 1/eV Ang^2 (value is
!!$*   adjusted in final decimal places to agree with VASP value; program
!!$*   checks for discrepancy of any results between this and VASP values)
! (Oleg Rubel Jul 19, 2018)
! in VASP manual HSQDTM = (plancks CONSTANT/(2*PI))**2/(2*ELECTRON MASS)
! HSQDTM = RYTOEV*AUTOA*AUTOA
! AUTOA=0.529177249_q,RYTOEV=13.605826_q
! thus HSQDTM = 13.605826*0.529177249*0.529177249 = 0.262465822502110

!REAL(kind=8), PARAMETER :: c=0.262465831 ! original from waveTrans, but there was a case when it did not work
!REAL(kind=8), PARAMETER :: c=0.262465822502110 ! according to constant.inc, but there was a case when it did not work
REAL(kind=8), PARAMETER :: c=0.26246585 ! that works error free so far...

ncnt=0
do ig3=0,2*nb3max
    ig3p=ig3
    if (ig3.gt.nb3max) ig3p=ig3-2*nb3max-1
    do ig2=0,2*nb2max
        ig2p=ig2
        if (ig2.gt.nb2max) ig2p=ig2-2*nb2max-1
        do ig1=0,2*nb1max
            ig1p=ig1
            if (ig1.gt.nb1max) ig1p=ig1-2*nb1max-1
            do j=1,3
                sumkg(j)=(wk(1)+ig1p)*b1(j)+ &
                    (wk(2)+ig2p)*b2(j)+(wk(3)+ig3p)*b3(j)
            enddo
            gtot=sqrt(sumkg(1)**2+sumkg(2)**2+sumkg(3)**2)
            etot=gtot**2/c
            if (etot.lt.ecut) then
                ncnt=ncnt+1
                igall(1,ncnt)=ig1p
                igall(2,ncnt)=ig2p
                igall(3,ncnt)=ig3p
            end if
        enddo
    enddo
enddo
IF (spinor) THEN
    if (2*ncnt.ne.nplane) then
        write(6,*) '*** error - computed no. of PWs', 2*ncnt, &
            ' != input no. of PWs', nplane
        stop
    endif
ELSE
    write(6,*) '*** error - computed no. of PWs', ncnt, &
            ' != input no. of PWs', nplane
    if (ncnt.ne.nplane) then
        write(6,*) '*** error - computed no. of PWs', ncnt, &
            ' != input no. of PWs', nplane
        if (ncnt*2.eq.nplane) then
            write(6,*) '2 x computed no. of PWs = input no. of PWs'
            write(6,*) 'Most likely you have a spinor WF. Try option -ncl'
            write(6,*) 'fold2Bloch WAVECAR FX:FY:FZ -ncl'
        endif
        stop
    endif
ENDIF
if (ncnt.gt.npmax) then
    write(6,*) '*** error - plane wave count', ncnt, &
        ' exceeds estimate', npmax
    stop
endif

END SUBROUTINE gvect
