!!! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!!! read_wavcar_head
!!!
!!! This subroutine will read WAVCAR first 2 records to determine
!!! the full record length, lattice vectors, number of k-points, etc
!!!
!!! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

SUBROUTINE read_wavcar_head(wffname, spinor, & ! <- args in 
            nrecl, a1, a2, a3, ecut, nwk, nband, nspin) ! -> args out 
IMPLICIT NONE
CHARACTER(len=*), INTENT(in) :: wffname ! name of the VASP WF file
LOGICAL, INTENT(in) :: spinor ! spinor WF
INTEGER, INTENT(out) :: nrecl ! record length WAVCAR
INTEGER, INTENT(out) :: nspin ! number of spins
INTEGER, INTENT(out) :: nwk ! number of k-points
INTEGER, INTENT(out) :: nband ! number of bands
REAL(kind=8), INTENT(out) :: ecut ! cut-off energy [eV]
REAL(kind=8), DIMENSION(1:3), &
    INTENT(out) :: a1, a2, a3 ! real space lattice vectors [Ang]
INTEGER :: nprec ! precision identifier
INTEGER :: iost ! I/O status
INTEGER :: j ! dimension counter
REAL(kind=8)  :: xnrecl ! record length WAVCAR
REAL(kind=8)  :: xnspin ! number of spins
REAL(kind=8)  :: xnprec ! precision identifier
REAL(kind=8)  :: xnwk ! number of k-points
REAL(kind=8)  :: xnband ! number of bands

!! Read 1st record

write(6,*) ' '
WRITE(6,*) 'Reading WF file head...'
INQUIRE (IOLENGTH=nrecl) xnrecl, xnspin, xnprec
WRITE (6,*) 'nrecl=', nrecl
OPEN(unit=12, file=TRIM(wffname), access='direct', FORM='UNFORMATTED', &
    recl=nrecl, iostat=iost, status='old')
WRITE(*,*) 'iost=', iost
IF (iost.ne.0) WRITE(6,*) 'open error - iostat =', iost            
READ(unit=12,rec=1) xnrecl, xnspin, xnprec
CLOSE(unit=12)
nrecl=nint(xnrecl) ! convert to integers
nspin=nint(xnspin)
nprec=nint(xnprec)
WRITE(6,*) 'record length  =',nrecl,' spins =',nspin, &
     ' prec flag ',nprec

!! Check spinor option activated correctly

IF((nspin.eq.2) .and. (spinor)) THEN
    WRITE(6,*) 'The WAVECAR has 2-spin wavefunction, while option [-ncl] is acticated'
    WRITE(6,*) 'Most likely the wavefunction is not spinor and [-ncl] option should be eliminated'
    STOP
ENDIF

!! Chech precision compatibility

IF(nprec.eq.45210) THEN
    WRITE(6,*) '*** error - WAVECAR_double requires complex*16'
    STOP
ENDIF

!! Reopen WAVCAR with a proper record length, read 2nd record

OPEN(unit=12, file=wffname, access='direct', FORM='UNFORMATTED', recl=nrecl, &
     iostat=iost, status='old')
WRITE(*,*) 'iost=', iost
IF (iost.ne.0) WRITE(6,*) 'open error - iostat =', iost
read(unit=12,rec=2) xnwk,xnband,ecut,(a1(j),j=1,3),(a2(j),j=1,3), &
     (a3(j),j=1,3)
nwk=nint(xnwk)
nband=nint(xnband)
write(6,*) 'no. k points =',nwk
write(6,*) 'no. bands =',nband
write(6,*) 'max. energy =',sngl(ecut)
write(6,*) 'real space lattice vectors [Ang]:'
write(6,*) 'a1 =',(sngl(a1(j)),j=1,3)
write(6,*) 'a2 =',(sngl(a2(j)),j=1,3)
write(6,*) 'a3 =',(sngl(a3(j)),j=1,3)
WRITE(6,*) 'Reading WF file head: done'


END SUBROUTINE read_wavcar_head
