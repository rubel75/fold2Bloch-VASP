!!! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!!! kgen_prim
!!!
!!!	Finds new K Values for each group, depending on the number of folds
!!! in x, y, and z directions.
!!!
!!! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

SUBROUTINE kgen_prim (wk, folds, nnewk, & ! <- args in
                newk) ! -> args out
 implicit none
 REAL(kind=8), intent(in) :: wk(3)
 INTEGER, intent(in) :: folds(3), nnewk
 REAL(kind=8), INTENT(out) :: newk (3,nnewk) ! new k-points in primitive BZ
 real :: field_x, field_y, field_z
 REAL(kind=8) :: TX, TY, TZ
 REAL(kind=8) :: X, Y, Z
 integer :: FX, FY, FZ
 integer :: loop, i, j, k
 
 X = wk(1); Y = wk(2); Z = wk(3)
 FX = folds(1); FY = folds(2); FZ = folds(3)

 !Brillouin zone range
 field_x=0.5*FX
 field_y=0.5*FY
 field_z=0.5*FZ
 loop=1 

do i=0, (FX-1)
    if ((X+i).gt.(field_x)) then
        TX=X-(field_x*2)+i
    else 
	TX=X+i
    endif
    do j=0, (FY-1)
        if ((Y+j).gt.(field_y)) then
            TY=Y-(field_y*2)+j
        else 
        TY=Y+j
	endif
        do k=0,(FZ-1)
            if ((Z+k).gt.(field_z)) then
                TZ=Z-(field_z*2)+k
            else 
                TZ=Z+k
            endif
            newk(1,loop)=TX
            newk(2,loop)=TY
            newk(3,loop)=TZ
            loop=loop+1
        enddo
    enddo
enddo

!! Normalize

newk(1,:)=newk(1,:)/FX
newk(2,:)=newk(2,:)/FY
newk(3,:)=newk(3,:)/FZ

END SUBROUTINE kgen_prim