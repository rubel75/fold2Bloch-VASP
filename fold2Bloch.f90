!!! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! This code reads VASP WAVECAR file and computes Bloch spectral weights
! that are used to plot an effective band structure of supercells.
!
! Execution:
!   $ fold2Bloch WAVECAR "P11 P12 P13:P21 P22 P23:P31 P32 P33" [-ncl]
!
! Options:
!   WAVECAR -- name of the input VASP wavefunction file
!   "P11 P12 P13:P21 P22 P23:P31 P32 P33" --  transformation matrix from 
!       primitive lattice vectors a_p to supercell lattice vectors a_s 
!       (same as in VESTA or the Bilbao Crystallographic Server):
!
!               a_s(i) = sum_j a_p(j)*P(j,i)      i,j = 1, 2, 3
!
!   -ncl -- the WAVECAR is produced by vasp_ncl code, which implies
!           that the wavefunctions are spinors
!           (default: assumption that WAVECAR comes from vasp_std or 
!           vasp_gam)
!
! Output:
!   Data are written to the file WAVECAR*.f2b in the following format
!   KX, KY, KZ, Eigenvalue (eV), Bloch weight
!
! Authors:
!   Oleg Rubel (McMaster University <rubelo@mcmaster.ca>)
!   Randall Feenstra (Carnegie Mellon University)
!   Michael Widom (Carnegie Mellon University)
!   Anton Bokhanchuk (Confederation College)
!
! Note: Reading wavefunctions is adapted from the WaveTrans code
!       http://www.andrew.cmu.edu/user/feenstra/wavetrans
!!! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

PROGRAM fold2Bloch

    ! Variables
    
    implicit none
    CHARACTER(len=256) :: arg1, arg2, arg3 ! command line input arguments
    CHARACTER(len=256) :: wffname ! wave function file name
    INTEGER :: &
        nrecl, &! record length WAVCAR
        nspin, &! number of spins
        nwk, &! number of k-points
        nband, &! number of bands
        nb1max, nb2max, nb3max, &! max number of G vect.
        npmax, &
        nnewk, &! number of new k-points after unfolding
        idf2b, &! ID for the output file
        irec, isp, iwk, i, iband, j, iplane, &! counters
        nplane, &! number of plane waves
        Dp2s(3,3) ! matrix to convert primitive cell to a supercell 
                  ! in direct (real) space
    INTEGER, ALLOCATABLE :: igall(:,:) ! G vectors in units of b1, b2, b3
    REAL(kind=8) :: &
        ecut, &! cut-off energy [eV]
        c, pi, &! parameters
        wk(3), &! k-point
        xnplane ! number of plane waves (real format)
    REAL(kind=8), DIMENSION(1:3) :: &
        a1, a2, a3, &! real space latt.vect.[Ang]
        b1, b2, b3 ! recip. latt.vect.[rad/Ang]
    REAL(kind=8), ALLOCATABLE :: &
        occ(:), &! occupancy of eigenstates
        newk(:,:), &! new k-points in primitive BZ
        w(:), w1(:), w2(:) ! Bloch spectral weights
    COMPLEX*8, ALLOCATABLE :: coeff(:) ! PW coefficients
    COMPLEX*16, ALLOCATABLE :: cener(:) ! eigenvalues
    CHARACTER(len=256) :: f2bfile ! name of the output file
    LOGICAL :: spinor ! WF is spinor?
    
    !!$*   constant 'c' below is 2m/hbar**2 in units of 1/eV Ang^2 (value is
    !!$*   adjusted in final decimal places to agree with VASP value; program
    !!$*   checks for discrepancy of any results between this and VASP values)
    
    data c/0.262465831d0/ 
    pi=4.*atan(1.)
    
    !! Print version
    
    WRITE(6,*) 'fold2Bloch for VASP version Aug 23, 2022'
    
    !! Get command line input arguments
    
    CALL GETARG(1,arg1)
    CALL GETARG(2,arg2)
    CALL GETARG(3,arg3)
    CALL check_in_args(COMMAND_ARGUMENT_COUNT(), arg1, arg2, arg3, & ! <- args in
        wffname, Dp2s, spinor) ! -> args out
    
    !! Read head of the WAVECAR file
    
    CALL read_wavcar_head(wffname, spinor, & ! <- args in 
        nrecl, a1, a2, a3, ecut, nwk, nband, nspin) ! -> args out 
    
    !! Compute reciprocal properties
    
    CALL recip(ecut, a1, a2, a3, spinor, & ! <- args in 
        b1, b2, b3, nb1max, nb2max, nb3max, npmax) ! -> args out 
    
    !! Allocate arrays for occupancy, eigenvalues, unfolded k-points, and
    !! Bloch spectral weight
    
    allocate(occ(nband))
    allocate(cener(nband))
    nnewk = Dp2s(1,1)*Dp2s(2,2)*Dp2s(3,3) + Dp2s(1,2)*Dp2s(2,3)*Dp2s(3,1) &
        +Dp2s(1,3)*Dp2s(2,1)*Dp2s(3,2) - Dp2s(3,1)*Dp2s(2,2)*Dp2s(1,3) &
        -Dp2s(1,1)*Dp2s(3,2)*Dp2s(2,3) - Dp2s(2,1)*Dp2s(1,2)*Dp2s(3,3) ! det(Dp2s)
    allocate(newk(3,nnewk))
    allocate(w(nnewk))
    IF (spinor) allocate(w1(nnewk), w2(nnewk))
    
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
            write(6,'(A,I0,A,I0)') 'k point # ', iwk, ' of ', nwk
            write(6,'(A,3F8.4)') 'k value =',(sngl(wk(j)),j=1,3)
    
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
            
            CALL kgen_prim(wk, Dp2s, nnewk, DBLE(1e-6), & ! <- args in
                newk) ! -> args out
            ! here 1e-6 is a tolerance (deltak) for finding unique k points
    
            !! Loop over bands
    
            do iband=1,nband
                irec=irec+1
                read(unit=12,rec=irec) (coeff(iplane), &
                    iplane=1,nplane)
    
                !! Compute Bloch spectral weights (not normalized)
                
                IF (spinor) THEN ! compute spinors
                    ! process alpha spinor component
                    CALL bloch_w(wk, newk, nnewk, Dp2s, DBLE(1e-6), &! <-- args in
                        igall, coeff(1:nplane/2), nplane/2, &! <-- args in
                        w1) ! --> args out
                    ! process beta spinor component
                    CALL bloch_w(wk, newk, nnewk, Dp2s, DBLE(1e-6), &! <-- args in
                        igall, coeff(nplane/2+1:nplane), nplane/2, &! <-- args in
                        w2) ! --> args out
                    w = w1 + w2 ! add weights for alpha and beta spinor components
                ELSE
                    CALL bloch_w(wk, newk, nnewk, Dp2s, DBLE(1e-6), &! <-- args in
                        igall, coeff, nplane, &! <-- args in
                        w) ! --> args out
                ENDIF
    
                !! normalize sum weights to 1
    
                w = w/SUM(w)
                
                !! Write unfolding result to *.f2b
                
                DO i=1, nnewk
                    WRITE(idf2b,100) newk(1,i), newk(2,i), newk(3,i), &
                        REAL(cener(iband)), w(i)
    100             format(f11.6, f11.6, f11.6, f11.6, f11.6)
                END DO
    
            enddo ! bands        
            !! Deallocate arrays for G-vectors and PW coefficients
            deallocate (igall)
            deallocate (coeff)
        enddo ! k-point
        ! Summary per spin
        WRITE(6,*)
        WRITE(6,'(A)') "\/\/\/ UNFOLDING FINISHED SUCCESSFULLY \/\/\/"
        WRITE(6,'(A,I0)') "Number of K points processed: ", nwk
        IF (spinor) THEN
            inquire(unit=idf2b, name=f2bfile) ! get file name from unit number
            WRITE(6,'(A)') "Data were written to: "//trim(f2bfile)
        ELSE
            inquire(unit=idf2b, name=f2bfile)
            WRITE(6,'(A)') "Data were written to: "//trim(f2bfile)
        ENDIF
        WRITE(6,'(A)') "Data format: KX, KY, KZ, Eigenvalue(eV), Bloch weight"
        WRITE(6,'(A)') "Try plotting tools available from "//&
            "https://github.com/rubel75/fold2Bloch-VASP/tree/master/utils"
        WRITE(6,'(A)') "You will need reciprocal lattice vectors "//&
            "b1, b2, b3 from the output above"
        close(unit=idf2b) ! *.f2b
    enddo ! spin
    
    close(unit=12) ! WAVECAR
    
    ! Epilogue
    
    WRITE(6,'(A)') ''
    WRITE(6,'(A)') 'If you find the results useful and publishable, '//&
        'we will appreciate citing the following papers:'
    WRITE(6,'(A)') '[1] O. Rubel, A. Bokhanchuk, S. J. Ahmed, and E. Assmann '//&
        '"Unfolding the band structure of disordered solids: from bound states '//&
        'to high-mobility Kane fermions", Phys. Rev. B 90, 115202 (2014).'
    WRITE(6,'(A)') '[2] L.-W. Wang, L. Bellaiche, S.-H. Wei, and A. Zunger '//&
        '"Majority representation of alloy electronic states", '//&
        'Phys. Rev. Lett. 80, 4725 (1998).'
    
    END PROGRAM fold2Bloch