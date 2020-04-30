!!! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! This code reads VASP WAVECAR file and computes Bloch spectral weights
! that are used to plot an effective band structure of supercells.
!
! Execution:
!   $ fold2Bloch WAVECAR FX:FY:FZ [-ncl]
!
! Options:
!   WAVECAR -- name of the input VASP wavefunction file
!   FX, FY, FZ -- multiplicity of the primitive cell that was used to 
!                 construct a supercell
!   -ncl -- the WAVECAR is produced by vasp_ncl code, which implies
!           that the wavefunctions are spinors
!           (defauls: assumption that WAVECAR comes from vasp_std or 
!           vasp_gam)
!
! Output:
!   Data are written to the file WAVECAR*.f2b in the following format
!   KX, KY, KZ, Eigenvalue (eV), Bloch weight
!
! Authors:
!   Oleg Rubel (McMaster University <oleg.v.rubel@gmail.com>)
!   Randall Feenstra (Carnegie Mellon University)
!   Michael Widom (Carnegie Mellon University)
!   Anton Bokhanchuk (Confederation College)
!
! Note: Reading wavefunctions is adapted from the WaveTrans code
!       http://www.andrew.cmu.edu/user/feenstra/wavetrans
!!! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

PROGRAM fold2Bloch

! Variables

implicit real*8 (a-h, o-z)
CHARACTER(len=256) :: arg1, arg2, arg3 ! command line imput arguments
INTEGER :: folds(3) ! number of folds X,Y,Z
CHARACTER(len=256) :: wffname ! wave function file name
INTEGER :: nrecl ! record length WAVCAR
INTEGER :: nspin ! number of spins
INTEGER :: nwk ! number of k-points
INTEGER :: nband ! number of bands
REAL(kind=8) :: ecut ! cut-off energy [eV]
REAL(kind=8), DIMENSION(1:3) :: a1, a2, a3 ! real space latt.vect.[Ang]
REAL(kind=8), DIMENSION(1:3) :: b1, b2, b3 ! recip. latt.vect.[rad/Ang]
INTEGER :: nb1max, nb2max, nb3max ! max number of G vect.
INTEGER :: npmax
COMPLEX*8, ALLOCATABLE :: coeff(:) ! PW coefficients
COMPLEX*16, ALLOCATABLE :: cener(:) ! eigenvalues
REAL(kind=8), ALLOCATABLE :: occ(:) ! occupancy of eigenstates
INTEGER, ALLOCATABLE :: igall(:,:) ! G vectors in units of b1, b2, b3
REAL(kind=8) :: wk(3) ! k-point
REAL(kind=8), ALLOCATABLE :: newk(:,:) ! new k-points in primitive BZ
INTEGER :: nnewk ! number of new k-points after unfolding
REAL(kind=8), ALLOCATABLE ::  w(:) ! Bloch spectral weights
INTEGER :: idf2b ! ID for the output file
CHARACTER(len=256) :: f2bfile ! name of the output file
LOGICAL :: spinor ! WF is spinor?

!!$*   constant 'c' below is 2m/hbar**2 in units of 1/eV Ang^2 (value is
!!$*   adjusted in final decimal places to agree with VASP value; program
!!$*   checks for discrepancy of any results between this and VASP values)

data c/0.262465831d0/ 
pi=4.*atan(1.)

!! Print version

WRITE(6,*) 'fold2Bloch for VASP version Oct 11, 2017'

!! Get command line imput arguments

CALL GETARG(1,arg1)
CALL GETARG(2,arg2)
CALL GETARG(3,arg3)
CALL check_in_args(COMMAND_ARGUMENT_COUNT(), arg1, arg2, arg3, & ! <- args in
    wffname, folds, spinor) ! -> args out 

write (*,*) "wffname=", TRIM(wffname)

!! Read head of the WAVCAR file

CALL read_wavcar_head(wffname, spinor, & ! <- args in 
    nrecl, a1, a2, a3, ecut, nwk, nband, nspin) ! -> args out 

!! Compute reciprocal properties

CALL recip(ecut, a1, a2, a3, spinor, & ! <- args in 
    b1, b2, b3, nb1max, nb2max, nb3max, npmax) ! -> args out 

!! Allocate arrays for occupancy, eigenvalues, unfolded k-points, and
!! Bloch spectral weight

allocate(occ(nband))
allocate(cener(nband))
nnewk = folds(1)*folds(2)*folds(3)
allocate(newk(3,nnewk))
allocate(w(nnewk))

!! Open output file for writing

!! Begin loops over spin, k-points and bands

irec=2
do isp=1,nspin
    write(*,*) ' '
    write(*,*) '******'
    write(*,*) 'reading spin ',isp, ' of ', nspin
    idf2b=13+isp
    CALL open_f2b(idf2b, isp, spinor, & ! <- args in
        f2bfile) ! -> args out
    
    !! Loop over k-points
    
    do iwk=1,nwk
        irec=irec+1
        read(unit=12,rec=irec) xnplane,(wk(i),i=1,3), &
            (cener(iband),occ(iband),iband=1,nband)
        nplane=nint(xnplane)
        write(6,*) 'k point #', iwk, ' of ', nwk
        write(6,*) 'k value =',(sngl(wk(j)),j=1,3)

        !! Allocate arrays for G-vectors and PW coefficients
        
        IF (spinor) THEN
            ALLOCATE(igall(3,nplane/2))
        ELSE
            ALLOCATE(igall(3,nplane))
        ENDIF
        ALLOCATE(coeff(nplane))

        !! Generate G-vectors for the specific k-point

        CALL gvect(b1, b2, b3, nb1max, nb2max, nb3max, & ! <- args in 
            nplane, npmax, wk, ecut, spinor, & ! <- args in 
            igall) ! -> args out
            
        !! Generated new unfolded k-points in primitive BZ
        
        CALL kgen_prim (wk, folds, nnewk, & ! <- args in
            newk) ! -> args out

        !! Loop over bands

        do iband=1,nband
            irec=irec+1
            read(unit=12,rec=irec) (coeff(iplane), &
                iplane=1,nplane)

            !! Compute Bloch spectral weight
            
            IF (spinor) THEN ! compute 1st spinor
                CALL bloch_w(folds, nplane/2, igall, & ! <- args in 
                    coeff(1:nplane/2), nnewk, & ! <- args in 
                    w) ! -> args out
            ELSE
                CALL bloch_w(folds, nplane, igall, & ! <- args in 
                    coeff, nnewk, & ! <- args in 
                    w) ! -> args out
            ENDIF
            
            !! Write unfolding result to *.f2b
            
            DO i=1, nnewk
                WRITE(idf2b,100) newk(1,i), newk(2,i), newk(3,i), &
                    REAL(cener(iband)), w(i)
100             format(f11.6, f11.6, f11.6, f11.6, f11.6)
            END DO
            
            !! Repeat for 2nd spinor if needed

            IF (spinor) THEN ! compute 2nd spinor
                CALL bloch_w(folds, nplane/2, igall, & ! <- args in 
                    coeff(nplane/2+1:nplane), nnewk, & ! <- args in 
                    w) ! -> args out
                DO i=1, nnewk
                    WRITE(idf2b+1,100) newk(1,i), newk(2,i), newk(3,i), &
                        REAL(cener(iband)), w(i)
                END DO
            ENDIF

        enddo ! bands        
        !! Deallocate arrays for G-vectors and PW coefficients
        deallocate (igall)
        deallocate (coeff)
    enddo ! k-point
    ! Summary per spin
    WRITE(6,*)
    WRITE(6,*) "     \/\/\/ UNFOLDING FINISHED SUCCESSFULLY \/\/\/"
    WRITE(6,*) "     Number of K points processed:", nwk
    WRITE(6,*) "     Data were written to: ", trim(f2bfile)
    WRITE(6,*) "     Data format: KX, KY, KZ, Eigenvalue(eV), Bloch weight"
    WRITE(6,*) "     Try matlab script ubs_dot.m for plotting."
    WRITE(6,*) "     You will need reciprocal lattice vectors b1, b2, b3 from the output above"
    close(unit=idf2b) ! *.f2b
    IF (spinor) close(unit=idf2b+1) ! *.f2b
enddo ! spin

close(unit=12) ! WAVECAR

END PROGRAM fold2Bloch
