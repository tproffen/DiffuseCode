MODULE kuplot_diffev_mod
!
!  Routines to plot fit results
!
!
INTEGER :: k_diff_gen
INTEGER :: k_diff_member
INTEGER :: k_diff_children
INTEGER :: k_diff_dimx
CHARACTER (LEN=1024) :: k_diff_trialfile
CHARACTER (LEN=1024) :: k_diff_resfile
CHARACTER (LEN=1024) :: k_diff_logfile
CHARACTER (LEN=1024) :: k_diff_sumfile
CHARACTER (LEN=1024) :: k_diff_curfile
CHARACTER (LEN=16), DIMENSION(:), ALLOCATABLE :: k_diff_name
!
CONTAINS
!
!*******************************************************************************
!
SUBROUTINE do_diffev_plot(line, length)
!
USE errlist_mod
USE get_params_mod
USE take_param_mod
!
IMPLICIT NONE
!
CHARACTER(LEN=*), INTENT(INOUT) :: line
INTEGER         , INTENT(INOUT) :: length
!
CHARACTER (LEN=10), PARAMETER :: genfile = 'GENERATION'
CHARACTER (LEN= 3), PARAMETER :: is_old  = 'old'
INTEGER, PARAMETER :: IRD = 77
INTEGER, PARAMETER :: MAXW = 4
!
CHARACTER (LEN=1024) :: string
CHARACTER (LEN=1024), DIMENSION(1:MAXW) :: cpara
INTEGER             , DIMENSION(1:MAXW) :: lpara
!
INTEGER              :: io_status
INTEGER              :: i
INTEGER              :: ianz
!
INTEGER, PARAMETER :: NOPTIONAL = 1
CHARACTER(LEN=1024), DIMENSION(NOPTIONAL) :: oname   !Optional parameter names
CHARACTER(LEN=1024), DIMENSION(NOPTIONAL) :: opara   !Optional parameter strings returned
INTEGER            , DIMENSION(NOPTIONAL) :: loname  !Lenght opt. para name
INTEGER            , DIMENSION(NOPTIONAL) :: lopara  !Lenght opt. para name returned
REAL               , DIMENSION(NOPTIONAL) :: owerte   ! Calculated values
INTEGER, PARAMETER                        :: ncalc = 1 ! Number of values to calculate 
INTEGER :: ipartial = 0
!
DATA oname  / 'partial'/
DATA loname /  7       /
opara  =  (/ '0.0000'  /)   ! Always provide fresh default values
lopara =  (/  6        /)
owerte =  (/  0.0      /)
!
!
!
CALL oeffne(IRD, genfile, is_old)
IF(ier_num /= 0) THEN
   CLOSE(IRD)
   RETURN
ENDIF
READ(IRD, '(a)', iostat=io_status) string      ! Generation header line
READ(IRD, *    , iostat=io_status) k_diff_gen, k_diff_member, k_diff_children, k_diff_dimx ! Generation values
READ(IRD, '(a)', iostat=io_status) string             ! trialfile  header line
READ(IRD, '(a)', iostat=io_status) k_diff_trialfile   ! trialfile  header line
READ(IRD, '(a)', iostat=io_status) string             ! result     header line
READ(IRD, '(a)', iostat=io_status) k_diff_resfile     ! result     header line
READ(IRD, '(a)', iostat=io_status) string             ! logfile    header line
READ(IRD, '(a)', iostat=io_status) k_diff_logfile     ! logfile    header line
READ(IRD, '(a)', iostat=io_status) string             ! Summary    header line
READ(IRD, '(a)', iostat=io_status) k_diff_sumfile     ! Summary    header line
READ(IRD, '(a)', iostat=io_status) string             ! Current    header line
READ(IRD, '(a)', iostat=io_status) k_diff_curfile     ! Current    header line
skip: DO
   READ(IRD, '(a)', iostat=io_status) string
   IF(string(1:11)=='# Parameter') EXIT skip
ENDDO skip
ALLOCATE(k_diff_name(-1:k_diff_dimx))
names: DO i=-1, k_diff_dimx
   READ(IRD, '(a)', iostat=io_status) string
   k_diff_name(i) = string(1:16)
ENDDO names
CLOSE(IRD)
!
CALL get_params (line, ianz, cpara, lpara, MAXW, length)
CALL get_optional(ianz, MAXW, cpara, lpara, NOPTIONAL,  ncalc, &
                  oname, loname, opara, lopara, owerte)
ipartial = NINT(owerte(1))
IF(ianz==1) THEN
   CALL kpara(MAXW, cpara, lpara, ianz, ipartial)
ELSEIF(ianz==3) THEN
   CALL kpar_par(MAXW, cpara, lpara, ianz, ipartial)
ELSE
   ier_num = -6
   ier_typ = ER_COMM
ENDIF
!
DEALLOCATE(k_diff_name)
!
END SUBROUTINE do_diffev_plot
!
!*******************************************************************************
!
SUBROUTINE kpara(MAXW, cpara, lpara, ianz, ipartial)
!
!  Create a plot of DIFFEV parameter versus Generation
!
USE kuplot_config
USE kuplot_mod
USE ber_params_mod
USE errlist_mod
USE get_params_mod
!
IMPLICIT NONE
!
INTEGER                            , INTENT(IN) :: MAXW
CHARACTER(LEN=*), DIMENSION(1:MAXW), INTENT(IN) :: cpara
INTEGER         , DIMENSION(1:MAXW), INTENT(IN) :: lpara
INTEGER                            , INTENT(IN) :: ianz
INTEGER                            , INTENT(IN) :: ipartial
!
CHARACTER (LEN=1024) :: infile
CHARACTER (LEN=1024) :: string
CHARACTER (LEN=   4) :: bef
INTEGER :: i
INTEGER              :: ipar1
INTEGER              :: success
INTEGER              :: length
LOGICAL              :: lname
LOGICAL              :: lexist
REAL   , DIMENSION(1:MAXW) :: werte
!
success = -1
lname   = .FALSE.
!
is_name: DO i=-1, k_diff_dimx
   IF(cpara(1)==k_diff_name(i)) THEN
      success = 0
      lname   = .TRUE.
      ipar1   = i
      EXIT is_name
   ENDIF
ENDDO is_name
IF(success/=0) THEN    ! Try to calculate parameter number
   CALL ber_params(ianz, cpara, lpara, werte, MAXW)
   IF(ier_num /=0) RETURN
   success = 0
   lname   = .FALSE.
   ipar1   = NINT(werte(1))
ENDIF
IF(ipar1<0 .OR. ipar1>k_diff_dimx) THEN   ! Parameter outside range
   ier_num = -63
   ier_typ = ER_APPL
   RETURN
ENDIF
!
!
IF(ipar1==0 .AND. ipartial>0) THEN
   infile = k_diff_sumfile(1:LEN_TRIM(k_diff_sumfile))//'.'//k_diff_name(ipar1)(1:LEN_TRIM(k_diff_name(ipar1)))
   WRITE(infile,'(a,''.'',i4.4)') infile(1:LEN_TRIM(infile)), ipartial
   INQUIRE(FILE=infile, EXIST=lexist)
   IF(.NOT.lexist) THEN      ! 
      ier_num = -64
      ier_typ = ER_APPL
      RETURN
   ENDIF
ELSE
infile = k_diff_sumfile(1:LEN_TRIM(k_diff_sumfile))//'.'//k_diff_name(ipar1)(1:LEN_TRIM(k_diff_name(ipar1)))
INQUIRE(FILE=infile, EXIST=lexist)
IF(.NOT.lexist) THEN      ! try par number as extension
   WRITE(infile,'(a,''.'',i4.4)')  k_diff_sumfile(1:LEN_TRIM(k_diff_sumfile)), ipar1
   INQUIRE(FILE=infile, EXIST=lexist)
   IF(.NOT.lexist) THEN      ! 
      ier_num = -64
      ier_typ = ER_APPL
      RETURN
   ENDIF
ENDIF
ENDIF
string = ' '
length = 1
CALL do_rese (string, length)
!
WRITE(string,'(a,a,a)') 'sc, ',infile(1:len_trim(infile)), ' , 1, 1, 2'
length = LEN_TRIM(string)
CALL do_load (string, length)
IF(ier_num/=0) RETURN
!
WRITE(string,'(a,a,a)') 'sc, ',infile(1:len_trim(infile)), ' , 1, 1, 3'
length = LEN_TRIM(string)
CALL do_load (string, length)
IF(ier_num/=0) RETURN
!
WRITE(string,'(a,a,a)') 'sc, ',infile(1:len_trim(infile)), ' , 1, 1, 4'
length = LEN_TRIM(string)
CALL do_load (string, length)
IF(ier_num/=0) RETURN
!
WRITE(string,'(a,a,a)') 'sc, ',infile(1:len_trim(infile)), ' , 1, 1, 5'
length = LEN_TRIM(string)
CALL do_load (string, length)
IF(ier_num/=0) RETURN
!
!
DO i=1,len(1)
   dy(i) = y(offxy (4-1)+i) 
ENDDO
!
IF(ymin(1)==ymax(1) .AND. ymin(2)==ymax(2) .AND. ymin(3)==ymax(3) .AND. &
   ABS(ymin(1)-ymin(2))<=ymin(1)*1.E-5 ) THEN
   WRITE(string,'(a,i4,'' '',a)') 'Fixed parameter No. ', ipar1, k_diff_name(ipar1)
   length = LEN_TRIM(string)
   call echo(string,length)
ELSE
   bef = 'mtyp'
   length = 4
   string = '1, 3'
   CALL para_seti (string, length, imarktyp, 1, maxkurvtot, bef, -3, 5000, .FALSE.)
   string = '2, 3'
   CALL para_seti (string, length, imarktyp, 1, maxkurvtot, bef, -3, 5000, .FALSE.)
   string = '3, 3'
   CALL para_seti (string, length, imarktyp, 1, maxkurvtot, bef, -3, 5000, .FALSE.)
!
   bef = 'mcol'
   string = '1, 3'
   CALL para_seti (string, length, imarkcol, 1, maxkurvtot, bef,  1,   15, .FALSE.)
   string = '2, 1'
   CALL para_seti (string, length, imarkcol, 1, maxkurvtot, bef,  1,   15, .FALSE.)
   string = '3, 1'
   CALL para_seti (string, length, imarkcol, 1, maxkurvtot, bef,  1,   15, .FALSE.)
!
   bef = 'ltyp'
   string = '1, 1'
   CALL para_seti (string, length, ilinetyp, 1, maxkurvtot, bef,  0,    5, .FALSE.)
   string = '2, 1'
   CALL para_seti (string, length, ilinetyp, 1, maxkurvtot, bef,  0,    5, .FALSE.)
   string = '3, 1'
   CALL para_seti (string, length, ilinetyp, 1, maxkurvtot, bef,  0,    5, .FALSE.)
!
   bef = 'lcol'
   string = '1, 3'
   CALL para_seti (string, length, ilinecol, 0, maxkurvtot, bef,  1,   15, .FALSE.)
   string = '2, 1'
   CALL para_seti (string, length, ilinecol, 0, maxkurvtot, bef,  1,   15, .FALSE.)
   string = '3, 1'
   CALL para_seti (string, length, ilinecol, 0, maxkurvtot, bef,  1,   15, .FALSE.)
!
   bef = 'etyp'
   string = '1, 2'
   CALL para_seti (string, length, ierr    , 1, maxkurvtot, bef,  0,    3, .FALSE.)
   string = '2, 0'
   CALL para_seti (string, length, ierr    , 1, maxkurvtot, bef,  0,    3, .FALSE.)
   string = '3, 0'
   CALL para_seti (string, length, ierr    , 1, maxkurvtot, bef,  0,    3, .FALSE.)
!
   length = 6
   string = '1, 0.2'
   CALL set_sizemark (string, length)
   string = '2, 0.2'
   CALL set_sizemark (string, length)
   string = '3, 0.2'
   CALL set_sizemark (string, length)
!
   string = '1, 1, 2, 3'
   length = 10
   CALL set_kfra (string, length)
   string = ' '
   length = 1
   CALL set_skal (string, length)
   string = ' '
   length = 1
   CALL set_mark (string, length)
!
   string ='DIFFEV Refinement'
   length = 17
   CALL para_set_title (string, length, titel (iwin, iframe, 1) )
!
   IF(ipar1==0 .AND. ipartial>0) THEN
      WRITE(string,'(a,i4,'' '',a,'' partial '', i4)')'Parameter No. ',ipar1, k_diff_name(ipar1), ipartial
   ELSE
      WRITE(string,'(a,i4,'' '',a)')'Parameter No. ',ipar1, k_diff_name(ipar1)
   ENDIF
   length = LEN_TRIM(string)
   CALL para_set_title (string, length, titel (iwin, iframe, 2) )
   CALL para_set_achse (string, length, achse (iwin, iframe, 2), lachse(iwin, iframe, 2) )
   string = 'Generation'
   length = 10
   CALL para_set_achse (string, length, achse (iwin, iframe, 1), lachse(iwin, iframe, 1) )
!
   ifname (iwin, iframe) = .FALSE.
!
   CALL do_plot (.false.)
ENDIF
!
END SUBROUTINE kpara
!
!*******************************************************************************
!
SUBROUTINE kpar_par(MAXW, cpara, lpara, ianz, ipartial)
!
USE kuplot_config
USE kuplot_mod
USE ber_params_mod
USE errlist_mod
USE get_params_mod
!
IMPLICIT NONE
!
INTEGER                            , INTENT(IN) :: MAXW
CHARACTER(LEN=*), DIMENSION(1:MAXW), INTENT(INOUT) :: cpara
INTEGER         , DIMENSION(1:MAXW), INTENT(INOUT) :: lpara
INTEGER                            , INTENT(INOUT) :: ianz
INTEGER                            , INTENT(IN) :: ipartial
!
CHARACTER (LEN=1024) :: infile0
CHARACTER (LEN=1024) :: infile1
CHARACTER (LEN=1024) :: infile2
CHARACTER (LEN=1024) :: parfile
CHARACTER (LEN=1024) :: string
CHARACTER (LEN=   4) :: bef
INTEGER :: i, iianz
INTEGER :: igen
INTEGER              :: ipar1
INTEGER              :: ipar2
INTEGER              :: success
INTEGER              :: length
INTEGER              :: rvalue_ind
LOGICAL              :: lname
LOGICAL              :: lexist
REAL                 :: rvalue_min
REAL   , DIMENSION(1:MAXW) :: werte
!
success = -1
lname   = .FALSE.
! Determine Generation to plot
!
IF(cpara(1)=='last'.OR.cpara(1)== '-1') THEN
   igen = k_diff_gen-1
ELSE
   iianz = 1
   CALL ber_params(iianz, cpara, lpara, werte, MAXW)
   IF(ier_num/=0) RETURN
   igen = NINT(werte(1))
ENDIF
CALL del_params (1, ianz, cpara, lpara, maxw)
!
! Determine first parameter to plot
!
success = -1
is_name1: DO i=-1, k_diff_dimx
   IF(cpara(1)==k_diff_name(i)) THEN
      success = 0
      lname   = .TRUE.
      ipar1   = i
      EXIT is_name1
   ENDIF
ENDDO is_name1
IF(success/=0) THEN    ! Try to calculate parameter number
   iianz = 1
   CALL ber_params(iianz, cpara, lpara, werte, MAXW)
   IF(ier_num /=0) RETURN
   success = 0
   lname   = .FALSE.
   ipar1   = NINT(werte(1))
ENDIF
CALL del_params (1, ianz, cpara, lpara, maxw)
!
! Determine second parameter to plot
!
success = -1
is_name2: DO i=-1, k_diff_dimx
   IF(cpara(1)==k_diff_name(i)) THEN
      success = 0
      lname   = .TRUE.
      ipar2   = i
      EXIT is_name2
   ENDIF
ENDDO is_name2
IF(success/=0) THEN    ! Try to calculate parameter number
   CALL ber_params(ianz, cpara, lpara, werte, MAXW)
   IF(ier_num /=0) RETURN
   success = 0
   lname   = .FALSE.
   ipar2   = NINT(werte(1))
ENDIF
!
!
IF(igen==k_diff_gen-1) THEN   ! In last Generation
   IF(k_diff_curfile /= ' ') THEN ! Current File should exist
      parfile = k_diff_curfile(1:LEN_TRIM(k_diff_curfile))
   ELSE
      parfile = k_diff_logfile(1:LEN_TRIM(k_diff_logfile))
   ENDIF
ELSE
   parfile = k_diff_logfile(1:LEN_TRIM(k_diff_logfile))
ENDIF
! Inquire Rvalue file
IF(ipartial>0) THEN
   infile0 = parfile(1:LEN_TRIM(parfile))//'.'//k_diff_name(0)(1:LEN_TRIM(k_diff_name(0)))
   WRITE(infile0,'(a,''.'',i4.4)') infile0(1:LEN_TRIM(infile0)), ipartial
   INQUIRE(FILE=infile0, EXIST=lexist)
   IF(.NOT.lexist) THEN      ! 
      ier_num = -66
      ier_typ = ER_APPL
      RETURN
   ENDIF
ELSE
   infile0 = parfile(1:LEN_TRIM(parfile))//'.'//k_diff_name(0)(1:LEN_TRIM(k_diff_name(0)))
   INQUIRE(FILE=infile0, EXIST=lexist)
   IF(.NOT.lexist) THEN      ! try par number as extension
      WRITE(infile0,'(a,''.'',i4.4)')  parfile(1:LEN_TRIM(parfile)), 0
      INQUIRE(FILE=infile0, EXIST=lexist)
      IF(.NOT.lexist) THEN      ! 
         ier_num = -66
         ier_typ = ER_APPL
         RETURN
      ENDIF
   ENDIF
ENDIF
! Inquire Current/Logfile for Parameter 1
IF(ipar1/=-1) THEN
   IF(ipar1==0 .AND. ipartial>0) THEN
      infile1 = parfile(1:LEN_TRIM(parfile))//'.'//k_diff_name(ipar1)(1:LEN_TRIM(k_diff_name(ipar1)))
      WRITE(infile1,'(a,''.'',I4.4)') infile1(1:LEN_TRIM(infile1)), ipartial
      IF(.NOT.lexist) THEN      ! 
         ier_num = -65
         ier_typ = ER_APPL
         RETURN
      ENDIF
   ELSE
      infile1 = parfile(1:LEN_TRIM(parfile))//'.'//k_diff_name(ipar1)(1:LEN_TRIM(k_diff_name(ipar1)))
      INQUIRE(FILE=infile1, EXIST=lexist)
      IF(.NOT.lexist) THEN      ! try par number as extension
         WRITE(infile1,'(a,''.'',i4.4)')  parfile(1:LEN_TRIM(parfile)), ipar1
         INQUIRE(FILE=infile1, EXIST=lexist)
         IF(.NOT.lexist) THEN      ! 
            ier_num = -65
            ier_typ = ER_APPL
            RETURN
         ENDIF
      ENDIF
   ENDIF
ELSE
   infile1 = infile0
ENDIF
! Inquire Current/Logfile for Parameter 2
IF(ipar2/=-1) THEN
   IF(ipar2==0 .AND. ipartial>0) THEN
      infile2 = parfile(1:LEN_TRIM(parfile))//'.'//k_diff_name(ipar2)(1:LEN_TRIM(k_diff_name(ipar2)))
      WRITE(infile2,'(a,''.'',I4.4)') infile2(1:LEN_TRIM(infile2)), ipartial
      IF(.NOT.lexist) THEN      ! 
         ier_num = -65
         ier_typ = ER_APPL
         RETURN
      ENDIF
   ELSE
      infile2 = parfile(1:LEN_TRIM(parfile))//'.'//k_diff_name(ipar2)(1:LEN_TRIM(k_diff_name(ipar2)))
      INQUIRE(FILE=infile2, EXIST=lexist)
      IF(.NOT.lexist) THEN      ! try par number as extension
         WRITE(infile2,'(a,''.'',i4.4)')  parfile(1:LEN_TRIM(parfile)), ipar2
         INQUIRE(FILE=infile2, EXIST=lexist)
         IF(.NOT.lexist) THEN      ! 
            ier_num = -65
            ier_typ = ER_APPL
            RETURN
         ENDIF
      ENDIF
   ENDIF
ELSE
   infile2 = infile0
ENDIF
string = ' '
length = 1
CALL do_rese (string, length)
!
WRITE(string,'(a,a,a,i5,a)') 'sc, ',infile0(1:len_trim(infile0)), ',', igen,' , 1, 2'
length = LEN_TRIM(string)
CALL do_load (string, length)
IF(ier_num/=0) RETURN
!
! Determine lowest R-value and its member number
rvalue_min = y(1)
rvalue_ind = 1
DO i=2,len(1)    ! len is the length of data sets
   IF(y(i)<rvalue_min) THEN
      rvalue_min = y(i)
      rvalue_ind = i
   ENDIF
ENDDO
!
string = ' '
length = 1
CALL do_rese (string, length)
!
! Load data sets for parameters 2 and 1
IF(ipar2 == -1) THEN
   WRITE(string,'(a,a,a,i5,a)') 'sc, ',infile2(1:len_trim(infile2)),',', igen,' , 1, 2'
ELSE
   WRITE(string,'(a,a,a,i5,a)') 'sc, ',infile2(1:len_trim(infile2)),',', igen,' , 3, 2'
ENDIF
length = LEN_TRIM(string)
CALL do_load (string, length)
IF(ier_num/=0) RETURN
!
IF(ipar1 == -1) THEN
   WRITE(string,'(a,a,a,i5,a)') 'sc, ',infile1(1:len_trim(infile1)),',', igen,' , 1, 2'
ELSE
   WRITE(string,'(a,a,a,i5,a)') 'sc, ',infile1(1:len_trim(infile1)),',', igen,' , 3, 2'
ENDIF
length = LEN_TRIM(string)
CALL do_load (string, length)
IF(ier_num/=0) RETURN
!
! Create the actual curves to be plotted
WRITE(string,'(a,i10)')  'points, ', len(1)
length = LEN_TRIM(string)
CALL do_allocate (string, length)
IF(ier_num/=0) RETURN
!
DO i=1,len(1)
   x(offxy(3-1)+i) = x(offxy (2-1)+i) 
   y(offxy(3-1)+i) = x(offxy (1-1)+i) 
ENDDO
string = 'minimum, 1'
length = 10
CALL do_allocate (string, length)
IF(ier_num/=0) RETURN
x(offxy(4-1)+1) = x(offxy(3-1) + rvalue_ind)
y(offxy(4-1)+1) = y(offxy(3-1) + rvalue_ind)
!
!
bef = 'mtyp'
length = 4
string = '3, 3'
CALL para_seti (string, length, imarktyp, 1, maxkurvtot, bef, -3, 5000, .FALSE.)
string = '4, 3'
CALL para_seti (string, length, imarktyp, 1, maxkurvtot, bef, -3, 5000, .FALSE.)
!
bef = 'mcol'
string = '3, 1'
CALL para_seti (string, length, imarkcol, 1, maxkurvtot, bef,  1,   15, .FALSE.)
string = '4, 3'
CALL para_seti (string, length, imarkcol, 1, maxkurvtot, bef,  1,   15, .FALSE.)
!
bef = 'ltyp'
string = '3, 0'
CALL para_seti (string, length, ilinetyp, 1, maxkurvtot, bef,  0,    5, .FALSE.)
string = '4, 0'
CALL para_seti (string, length, ilinetyp, 1, maxkurvtot, bef,  0,    5, .FALSE.)
!
bef = 'lcol'
string = '3, 1'
CALL para_seti (string, length, ilinecol, 0, maxkurvtot, bef,  1,   15, .FALSE.)
string = '4, 3'
CALL para_seti (string, length, ilinecol, 0, maxkurvtot, bef,  1,   15, .FALSE.)
!
bef = 'etyp'
string = '3, 0'
CALL para_seti (string, length, ierr    , 1, maxkurvtot, bef,  0,    3, .FALSE.)
string = '4, 0'
CALL para_seti (string, length, ierr    , 1, maxkurvtot, bef,  0,    3, .FALSE.)
!
length = 6
string = '4, 0.2'
CALL set_sizemark (string, length)
length = 9
IF(ipar1==0 .AND. ipar2/=0) THEN
   string = '3, 0.8, x'
ELSEIF(ipar1/=0 .AND. ipar2==0) THEN
   string = '3, 0.8, y'
ELSEIF(ipar1/=0 .AND. ipar2/=0) THEN
   string = '3, 0.8, 2'
ENDIF
CALL set_sizemark (string, length)
!
string = '1, 3, 4'
length =  7
CALL set_kfra (string, length)
string = 'xmin[3], xmax[3], ymin[3]-(ymax[3]-ymin[3])*0.05,ymax[3]+(ymax[3]-ymin[3])*0.05 '
length = LEN_TRIM(string)
CALL set_skal (string, length)
string = ' '
length = 1
CALL set_mark (string, length)
!
string ='DIFFEV Correlation between parameters'
length = 37
CALL para_set_title (string, length, titel (iwin, iframe, 1) )
!
IF(ipar2==0 .AND. ipartial>0) THEN
WRITE(string,'(a,i5,a,i4,'' '',a,i4,'' '',a,'' partial '',i4)') 'Generation ',igen,'; Parameter No. ',&
                ipar1, k_diff_name(ipar1)(1:LEN_TRIM(k_diff_name(ipar1))),  &
                ipar2, k_diff_name(ipar2)(1:LEN_TRIM(k_diff_name(ipar2))), ipartial
ELSEIF(ipar1==0 .AND. ipartial>0) THEN
WRITE(string,'(a,i5,a,i4,'' '',a,i4,'' partial '',i4,'' '',a)') 'Generation ',igen,'; Parameter No. ',&
                ipar1, k_diff_name(ipar1)(1:LEN_TRIM(k_diff_name(ipar1))), ipartial, &
                ipar2, k_diff_name(ipar2)(1:LEN_TRIM(k_diff_name(ipar2)))
ELSE
WRITE(string,'(a,i5,a,i4,'' '',a,i4,'' '',a)') 'Generation ',igen,'; Parameter No. ',&
                ipar1, k_diff_name(ipar1)(1:LEN_TRIM(k_diff_name(ipar1))),  &
                ipar2, k_diff_name(ipar2)(1:LEN_TRIM(k_diff_name(ipar2)))
ENDIF
!
length = LEN_TRIM(string)
CALL para_set_title (string, length, titel (iwin, iframe, 2) )
IF(ipar2==0 .AND. ipartial>0) THEN
   string = k_diff_name(ipar1)
   WRITE(string,'(a,'' '',i4)') string(1:LEN_TRIM(string)), ipartial
ELSEIF(ipar1==0 .AND. ipartial>0) THEN
   string = k_diff_name(ipar1)
   WRITE(string,'(a,'' '',i4)') string(1:LEN_TRIM(string)), ipartial
ELSE
   string = k_diff_name(ipar1)
ENDIF
length = LEN_TRIM(string)
CALL para_set_achse (string, length, achse (iwin, iframe, 1), lachse(iwin, iframe, 1) )
string = k_diff_name(ipar2)
length = LEN_TRIM(string)
CALL para_set_achse (string, length, achse (iwin, iframe, 2), lachse(iwin, iframe, 2) )
!
ifname (iwin, iframe) = .FALSE.
!
CALL do_plot (.false.)
string = 'data,4'
length = 6
CALL kuplot_do_show (string, length)
!
END SUBROUTINE kpar_par
!
END MODULE kuplot_diffev_mod
