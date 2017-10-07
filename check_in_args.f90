!!! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!!! check_in_args
!!!
!!! This subroutine will perform check of the imput arguments that
!!! are passed in when fold2Bloch is executed. It is expected that
!!! there are only 2 or 3 arguments:
!!! WAVECAR (char) and 
!!! FX:FY:FZ (folds, integers) and 
!!! [-ncl] optional for spinor
!!!
!!! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

SUBROUTINE check_in_args(nargs, arg1, arg2, arg3, & ! <- args in
    wffname, folds, spinor) ! -> args out 
IMPLICIT NONE
INTEGER, INTENT(in) :: nargs ! number of input arguments
CHARACTER(len=*), INTENT(in) :: arg1 ! name of the VASP WF file
CHARACTER(len=*), INTENT(in) :: arg2 ! folds for unfolding, e.g., 1:2:3
CHARACTER(len=*), INTENT(in) :: arg3 ! [-ncl]
CHARACTER(len=*), INTENT(out) :: wffname ! name of the VASP WF file
LOGICAL, INTENT(out) :: spinor ! spinor WF
INTEGER, INTENT(out) :: folds(3) ! number of folds X,Y,Z
LOGICAL :: res ! used when checking the file presence
INTEGER :: ifarg ! position of ":" in the folds input argument
INTEGER :: ios ! I/O status
CHARACTER(len=256) :: arg2cut ! cut folds [1:]2:3

!! Check if to print help

IF ((nargs .eq. 1).and.(arg1(1:2) .eq. '-h')) THEN
    WRITE(6,*) 'Usage: fold2Bloch WAVECAR FX:FY:FZ [-ncl]'
    WRITE(6,*) 'Exit'
    STOP
ENDIF

!! Check the number of input arguments

IF ((nargs .lt. 2) .or. (nargs .gt. 3)) THEN
    WRITE(6,*) 'Number of imput arguments: ', nargs, &
        'while 2 or 3 is expected'
    WRITE(6,*) 'Usage: fold2Bloch WAVECAR FX:FY:FZ [-ncl]'
    WRITE(6,*) 'Exit'
    STOP
ENDIF

!! Check if WAVECAR exists

INQUIRE(file=arg1, exist=res)
if (not(res)) then
    WRITE(6,*) 'The file with wavefunctions "', TRIM(arg1), &
        '" is not found' 
    WRITE(6,*) 'Usage: fold2Bloch WAVECAR FX:FY:FZ [-ncl]'
    WRITE(6,*) 'Exit'
    STOP
endif
wffname=TRIM(arg1) ! name of the VASP WF file

!! Check "arg2" with folds to have a proper format FX:FY:FZ

ifarg=INDEX(arg2, ':')
IF (ifarg .eq. 0) THEN
    WRITE(6,*) 'The second input field "', TRIM(arg2), &
            '" with the number of folds is not recognized' 
    WRITE(6,*) 'Usage: fold2Bloch WAVECAR FX:FY:FZ [-ncl]'
    WRITE(6,*) 'Exit'
    STOP
ENDIF

!! Extract number of folds

READ (arg2(1:ifarg-1), *, iostat=ios) folds(1)
IF ((ios.ne.0).or.(folds(1).le.0)) THEN
    WRITE(6,*) 'FX=', folds(1)
    WRITE(6,*) 'The mumber of folds has to be a positive integer greater than 0'
    WRITE(6,*) 'Usage: fold2Bloch WAVECAR FX:FY:FZ'
    WRITE(6,*) 'Exit'
    STOP
ENDIF
arg2cut=arg2(ifarg+1:) ! trim FX
ifarg=INDEX(arg2cut, ':')
IF (ifarg .eq. 0) THEN
    WRITE(6,*) 'The second input field "', TRIM(arg2), &
            '" with the number of folds is not recognized' 
    WRITE(6,*) 'Usage: fold2Bloch WAVECAR FX:FY:FZ [-ncl]'
    WRITE(6,*) 'Exit'
    STOP
ENDIF
READ (arg2cut(1:ifarg-1),*, iostat=ios) folds(2)
IF ((ios.ne.0).or.(folds(2).le.0)) then
    WRITE(6,*) 'FY=', folds(2)
    WRITE(6,*) 'The mumber of folds has to be a positive integer greater than 0'
    WRITE(6,*) 'Usage: fold2Bloch WAVECAR FX:FY:FZ  [-ncl]'
    WRITE(6,*) 'Exit'
    STOP
ENDIF
READ(arg2cut(ifarg+1:),*, iostat=ios) folds(3)
IF ((ios.ne.0).or.(folds(3).le.0)) then
    WRITE(6,*) 'FZ=', folds(3)
    WRITE(6,*) 'The mumber of folds has to be a positive integer greater than 0'
    WRITE(6,*) 'Usage: fold2Bloch WAVECAR FX:FY:FZ [-ncl]'
    WRITE(6,*) 'Exit'
    STOP
ENDIF

!! Check "arg3" with [-ncl] option for spinors

IF (nargs .eq. 3) THEN
    IF (arg3 .eq. '-ncl') THEN
        spinor=.true.
        WRITE(6,*) 'Wave functions are spinors'
    ELSE
        WRITE(6,*) 'Third option is not recognized'
        WRITE(6,*) 'Usage: fold2Bloch WAVECAR FX:FY:FZ [-ncl]'
        WRITE(6,*) 'Exit'
        STOP
    ENDIF
ELSE
    spinor=.false.
ENDIF

write(6,*) 'Number of input arguments: ', nargs
write(6,*) 'VASP WF file name: ', TRIM(wffname)
write(6,*) 'Number of folds (FX, FY, FZ):', folds
write(6,*) 'Wavefunctions are spinors:', spinor

END SUBROUTINE check_in_args
