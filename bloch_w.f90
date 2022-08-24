!!! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!!! bloch_w
!!!
!!! Compute Bloch spectral weight by sorting PW coefficients 
!!! into appropriate groups.
!!!
!!! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
SUBROUTINE bloch_w(ks, kp, vscale, Dp2s, toldk, G, &! <-- args in
    pwcoeffz, NV, &! <-- args in
    w) ! --> args out
implicit none
! external vars
integer, intent(in) :: &
    vscale, &! primit. -> superc. volume scale (real space)
    Dp2s(3,3), &! primit. -> superc. transform matrix (real space)
    NV, &! length of PW coefficient vector
    G(3,NV) ! matrix of PW lattice vectors
REAL(kind=8), intent(in) :: &
    ks(3), &! k point in supercell
    kp(3,vscale), &! ks point unfolded to primitive BZ
    toldk ! tolerance for finding unique k points
COMPLEX*8, intent(in) :: &
    pwcoeffz(NV) ! plane wave coefficients r/z (real/complex)
REAL(kind=8), intent(out) :: &
    w(vscale) ! weights after unfolding (0-1)
! internal vars
REAL(kind=8) :: &
    Ds2p(3,3), &! matrix to transf. direct supercell to primitive vectors
    ksG(3), &! ks + G
    kptmpi(3) ! intermediate (non unique new k point)
integer :: &
    i, j, &! counter
    countw(vscale) ! count entries into weight bins
logical :: &
    matchfound ! used to identify when a k point match found in a group 'kp'

!! construct Ds2p matrix
! Ds2p = inv(Dp2s)
    Ds2p(1,1) = (  Dp2s(2,2)*Dp2s(3,3) - Dp2s(2,3)*Dp2s(3,2) ) / REAL(vscale,8)
    Ds2p(2,1) = (- Dp2s(2,1)*Dp2s(3,3) + Dp2s(2,3)*Dp2s(3,1) ) / REAL(vscale,8)
    Ds2p(3,1) = (  Dp2s(2,1)*Dp2s(3,2) - Dp2s(2,2)*Dp2s(3,1) ) / REAL(vscale,8)
    Ds2p(1,2) = (- Dp2s(1,2)*Dp2s(3,3) + Dp2s(1,3)*Dp2s(3,2) ) / REAL(vscale,8)
    Ds2p(2,2) = (  Dp2s(1,1)*Dp2s(3,3) - Dp2s(1,3)*Dp2s(3,1) ) / REAL(vscale,8)
    Ds2p(3,2) = (- Dp2s(1,1)*Dp2s(3,2) + Dp2s(1,2)*Dp2s(3,1) ) / REAL(vscale,8)
    Ds2p(1,3) = (  Dp2s(1,2)*Dp2s(2,3) - Dp2s(1,3)*Dp2s(2,2) ) / REAL(vscale,8)
    Ds2p(2,3) = (- Dp2s(1,1)*Dp2s(2,3) + Dp2s(1,3)*Dp2s(2,1) ) / REAL(vscale,8)
    Ds2p(3,3) = (  Dp2s(1,1)*Dp2s(2,2) - Dp2s(1,2)*Dp2s(2,1) ) / REAL(vscale,8)

! initialize weights and its counter
w = DBLE(0)
countw = 0

! sorts the PW coefficients into bins of weights
do i=1, NV
    ksG = ks + (/G(1,i), G(2,i), G(3,i)/)
    kptmpi = MATMUL(ksG, Ds2p) ! convert supercell -> primitive basis
    ! bring k points into the range [0,1)
    do j=1,3
        kptmpi(j) = MODULO( kptmpi(j), DBLE(1))
        if (1-kptmpi(j) .lt. toldk) then ! 0.99999 -> 1 -> 0
            kptmpi(j) = DBLE(0)
        endif
    enddo ! j
    ! compare the kptmpi point with the list 'kp' and find which group 
    ! it belongs
    matchfound = .false.
    do j=1, vscale
        if ( ALL(ABS(kptmpi- kp(:,j)) < toldk) ) then
            ! point matched to one on the 'kp' list
            matchfound = .true.
            ! store absolute value squared of the PW coefficient in
            ! appropriate bin 'w'
            w(j) = w(j) + REAL(pwcoeffz(i)*CONJG(pwcoeffz(i)),8)
            countw(j) = countw(j) + 1 ! update counter
            exit ! loop
        endif
    enddo ! j (kp)
    if (.not.matchfound) then
        write(*,*) 'ERROR in subroutine Sort: unable to match k point (', &
            kptmpi, ') to a point in the group kp listed below'
        do j=1, vscale
            write(*,*) 'kp(', j, ') = ', kp(:,j)
        enddo
        ERROR STOP 1
    endif
enddo ! i (PW)

! check if all bins got weights
do i=1, vscale
    if (countw(j) .eq. 0) then
        write(*,*) 'ERROR in subroutine Sort: unable to populate k point (', &
            kp(:,i), ') from the group kp listed below. This is unlikely.'
        do j=1, vscale
            write(*,*) 'kp(', j, ') = ', kp(:,j)
        enddo
        ERROR STOP 1
    endif
enddo

! return _not_ normalized weights
 
END SUBROUTINE bloch_w
