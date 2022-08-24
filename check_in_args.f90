!!! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!!! check_in_args
!!!
!!! This subroutine will perform check of the input arguments that
!!! are passed in when fold2Bloch is executed. It is expected that
!!! there are only 2 or 3 arguments:
!!! WAVECAR (char) and 
!!! FX:FY:FZ (folds, integers) and 
!!! [-ncl] optional for spinor
!!!
!!! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

SUBROUTINE check_in_args(nargs, arg1, arg2, arg3, & ! <- args in
    wffname, Dp2s, spinor) ! -> args out 
IMPLICIT NONE
INTEGER, INTENT(in) :: nargs ! number of input arguments
CHARACTER(len=*), INTENT(in) :: arg1 ! name of the VASP WF file
CHARACTER(len=*), INTENT(in) :: arg2 ! folds for unfolding, e.g., 1:2:3
CHARACTER(len=*), INTENT(in) :: arg3 ! [-ncl]
CHARACTER(len=*), INTENT(out) :: wffname ! name of the VASP WF file
LOGICAL, INTENT(out) :: spinor ! spinor WF
INTEGER, INTENT(out) :: Dp2s(3,3) ! number of folds X,Y,Z
LOGICAL :: res ! used when checking the file presence
INTEGER :: &
    i, &! position of ":" in the folds input argument
    ios, &! I/O status
    det ! [P] matrix determinant
CHARACTER(len=256) :: folds ! cut folds [1:]2:3

!! Check if to print help

IF ((nargs .eq. 1).and.(arg1(1:2) .eq. '-h')) THEN
    WRITE(6,'(A)') 'Usage: fold2Bloch WAVECAR'//&
        ' "P11 P12 P13:P21 P22 P23:P31 P32 P33" [-ncl]'
    WRITE(*,'(A)') 'Notes:'
    WRITE(*,'(A)') '(1) [P] matrix (internally called Dp2s) is used to '//&
        'transform primitive a_p to supercell a_s lattice '//&
        'vectors (same as in VESTA):'
    WRITE(*,'(A)') '       a_s(i) = sum_j a_p(j)*P(j,i)      i,j = 1, 2, 3'
    WRITE(*,'(A)') '(2) Use quotations to input the [P] matrix _exactly_ as '//&
        'shown in this help'
    STOP
ENDIF

!! Check the number of input arguments

IF ((nargs .lt. 2) .or. (nargs .gt. 3)) THEN
    WRITE(6,'(A)') 'Number of input arguments: ', nargs, &
        'while 2 or 3 is expected'
    goto 912 ! show input options and exit with error
ENDIF

!! Check if WAVECAR exists

INQUIRE(file=arg1, exist=res)
if (.not. res) then
    WRITE(6,'(A)') 'The file with wavefunctions "', TRIM(arg1), &
        '" is not found' 
    goto 912 ! show input options and exit with error
endif
wffname=TRIM(arg1) ! name of the VASP WF file

!! Check "arg2" with folds to have a proper format FX:FY:FZ

i=INDEX(arg2, ':')
IF (i .eq. 0) THEN
    WRITE(6,'(A)') 'The second input field "', TRIM(arg2), &
            '" with the [P] matrix of folds is not recognized' 
    goto 912 ! show input options and exit with error
ENDIF

!! Extract number of [P] matrix (internally we call it [Dp2s])

folds = arg2
read(folds(1:i-1), *, iostat=ios) Dp2s(1,1), Dp2s(1,2), Dp2s(1,3)
if (ios.ne.0) then
    write (*,*) 'Unable to read the Dp2s(1,:) matrix'
    write (*,*) 'Relevant input line = ', trim(folds)
    write (*,*) 'Parsed input section = ', folds(1:i-1)
    write (*,*) 'while expected 3 numerical values separated by space.'
    GOTO 912 ! print error, usage options, and exit
endif
folds=folds(i+1:)
i=INDEX(folds, ':')
if (i.eq.0) then
    write(*,*) 'folds = ', TRIM(folds)
    write(*,*) 'Unknown number of folds. See below or type: "fold2Bloch -h" '//&
        'for more information.'
    GOTO 912 ! print error, usage options, and exit
endif
read (folds(1:i-1),*, iostat=ios) Dp2s(2,1), Dp2s(2,2), Dp2s(2,3)
if (ios.ne.0) then
    write (*,*) 'Unable to read the Dp2s(2,:) matrix'
    write (*,*) 'Parsed input section = ', folds(1:i-1)
    write (*,*) 'while expected 3 numerical values separated by space.'
    GOTO 912 ! print error, usage options, and exit
endif
read(folds(i+1:),*, iostat=ios) Dp2s(3,1), Dp2s(3,2), Dp2s(3,3)
if (ios.ne.0) then
    write (*,*) 'Unable to read the Dp2s(3,:) matrix'
    write (*,*) 'Parsed input section = ', folds(i+1:)
    write (*,*) 'while expected 3 numerical values separated by space.'
    GOTO 912 ! print error, usage options, and exit
endif

!! Check determinant of [Dp2s]

det = Dp2s(1,1)*Dp2s(2,2)*Dp2s(3,3) + Dp2s(1,2)*Dp2s(2,3)*Dp2s(3,1) &
    +Dp2s(1,3)*Dp2s(2,1)*Dp2s(3,2) - Dp2s(3,1)*Dp2s(2,2)*Dp2s(1,3) &
    -Dp2s(1,1)*Dp2s(3,2)*Dp2s(2,3) - Dp2s(2,1)*Dp2s(1,2)*Dp2s(3,3)
if (det.le.0) then
    write (*,*) 'The input matrix Dp2s that transforms the real space'
    write (*,*) 'primitive lattice vectors to a supercell leads to'
    write (*,*) 'the following volume change: ', det
    write (*,*) 'Parsed input section = ', folds(i+1:)
    write (*,*) 'FYI: Dp2s = ', Dp2s
    write (*,*) 'This means that the Dp2s matrix is not positively defined.'
    write (*,*) 'Also the volume change should be an integer number.'
    write (*,*) 'Please review your Dp2s input and rerun fold2Bloch.'
    GOTO 912 ! print error, usage options, and exit
else
    write (*,'(A)') 'Dp2s input matrix is successfully parsed.'
    write (*,'(A,3(I5,X),A)') '       | ', Dp2s(1,:), '|'
    write (*,'(A,3(I5,X),A)') 'Dp2s = | ', Dp2s(2,:), '|'
    write (*,'(A,3(I5,X),A)') '       | ', Dp2s(3,:), '|'
    write (*,'(A,I0)') 'The primitive to supercell volume scale is: ', &
        det
endif

!! Check "arg3" with [-ncl] option for spinors

IF (nargs .eq. 3) THEN
    IF (arg3 .eq. '-ncl') THEN
        spinor=.true.
        WRITE(6,'(A)') 'Wave functions are spinors'
    ELSE
        WRITE(6,*) 'Third option is not recognized'
        goto 912 ! show input options and exit with error
    ENDIF
ELSE
    spinor=.false.
ENDIF

write(6,'(A,I0)') 'Number of input arguments: ', nargs
write(6,'(A,A)') 'VASP WF file name: ', TRIM(wffname)
write(6,'(A,L)') 'Wavefunctions are spinors: ', spinor

return

! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!                           Errors
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

! show input options and exit with error
912 CONTINUE
WRITE(6,'(A)') 'Usage: fold2Bloch WAVECAR'//&
' "P11 P12 P13:P21 P22 P23:P31 P32 P33" [-ncl]'
WRITE(*,'(A)') 'Notes:'
WRITE(*,'(A)') '(1) [P] matrix (internally called Dp2s) is used to '//&
'transform primitive a_p to supercell a_s lattice '//&
'vectors (same as in VESTA):'
WRITE(*,'(A)') '       a_s(i) = sum_j a_p(j)*P(j,i)      i,j = 1, 2, 3'
WRITE(*,'(A)') '(2) Use quotations to input the [P] matrix _exactly_ as '//&
'shown in this help'
WRITE(6,*) 'Exit'
ERROR STOP

END SUBROUTINE check_in_args