!!! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!!! open_f2b
!!!
!!! It opens the file *.f2b to write output of unfolding and also 
!!! checks if this file already exists. In later case add '_#' to the 
!!! file name.
!!!
!!! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

SUBROUTINE open_f2b(idf2b, isp, spinor, & ! <- args in
                f2bfile) ! -> args out
IMPLICIT NONE
INTEGER, INTENT(inout) :: idf2b ! ID for the output file
INTEGER, INTENT(in) :: isp ! spin 1 or 2
LOGICAL, INTENT(in) :: spinor ! WF spinor?
CHARACTER(len=256), INTENT(out) :: f2bfile ! name of the output file
CHARACTER(len=256) :: f2bfile1 ! name of the output file
INTEGER*4 :: digit ! to make output file unique, e.g., WAVECAR.f2b_3
CHARACTER(len=4) :: cdigit ! name of the output file
LOGICAL :: res ! used when checking the file presence

write (cdigit,'(I0)') isp
IF (spinor) THEN
    f2bfile='WAVECAR_spinor.f2b'
ELSE
    f2bfile='WAVECAR_spin' // TRIM(cdigit) // '.f2b'
ENDIF

digit=0 ! initialize
res=.True. ! assume output file is already present
f2bfile1=f2bfile
DO WHILE (res)
    INQUIRE(file=f2bfile1, exist=res)
    IF (res) THEN ! file exists
        WRITE(6,*) 'The output file "', TRIM(f2bfile1), &
            '" exists, generate an alternative name'
        digit=digit+1
        IF (digit .gt. 100) THEN
            WRITE(6,*) 'You need to clean previous output files (>100)'
            WRITE(6,*) 'Exit'
            STOP
        ENDIF
        write (cdigit,'(I0)') digit
        f2bfile1=TRIM(f2bfile) // '_' // TRIM(cdigit)
        write(6,*) 'Will try the new output file name: ', TRIM(f2bfile1)
    ENDIF
END DO

f2bfile=f2bfile1
OPEN(idf2b, file=f2bfile, status='new', action='write')

END SUBROUTINE open_f2b
