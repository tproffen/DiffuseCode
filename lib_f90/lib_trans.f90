module lib_trans_mod
!
!  Low level routines to transform a 3D data set into a different coordinate system
!  No modules from DISCUS with metric tensor etc are assumed, all variables are
!  passed down
!  
!+
!
use precision_mod
!
private
public lib_trans_menu
!public get_length
!
contains
!
!*******************************************************************************
!
subroutine lib_trans_menu(idata, value, laver, outfile, old_inc, old_eck, old_vi, old_a0,    &
           old_win, qvalues, VAL_PDF, VAL_3DPDF,   &  
           new_outfile, new_inc, new_eck, new_vi, new_data)
!-
! main menu for a generalized transformation of 3D-data sets
! if idata == -1 ! DISCUS output data; write to disk
!    idata ==  0 ! KUPLOT load h5 data; copy into local storage
!    idata >   0 ! KUPLOT         data; replace in memory
!+
!
use calc_expr_mod
use class_macro_internal
use errlist_mod
use doact_mod
use kdo_all_mod
use lib_errlist_func
use lib_length
use lib_macro_func
use macro_mod
use precision_mod
use prompt_mod
use str_comp_mod
use sup_mod
!
implicit none
!
integer                             , intent(in)  :: idata        ! Old data set number 
integer                             , intent(in)  :: value        ! Intensity, 3DPDF etc
logical                             , intent(in)  :: laver        ! Average value ; currently not needed
character(len=*)                    , intent(in)  :: outfile      ! Output file name passed down to hdf5_write
integer           , dimension(3)    , intent(in)  :: old_inc      ! Number data points  original
real(kind=PREC_DP), dimension(3,4)  , intent(in)  :: old_eck      ! The Corners (ll, lr, ul, tl)
real(kind=PREC_DP), dimension(3,3)  , intent(in)  :: old_vi       ! Old vectors vi(hkl,1-3)
real(kind=PREC_DP), dimension(3)    , intent(in)  :: old_a0       ! Lattice parameters 
real(kind=PREC_DP), dimension(3)    , intent(in)  :: old_win      ! Lattice angles 
real(kind=PREC_DP), dimension(old_inc(1),old_inc(2),old_inc(3)), intent(in) :: qvalues   ! Old data 
integer                             , intent(in)  :: VAL_PDF      ! Number for a PDF as output
integer                             , intent(in)  :: VAL_3DPDF    ! Number for a 3DPDF as output
character(len=*)                    , intent(out) :: new_outfile  ! New Output file name 
integer           , dimension(3)    , intent(out) :: new_inc      ! New pixel numbers
real(kind=PREC_DP), dimension(3,3)  , intent(out) :: new_vi       ! New increment vectors
real(kind=PREC_DP), dimension(3,4)  , intent(out) :: new_eck      ! New Corners (ll, lr, ul, tl) 
real(kind=PREC_DP), dimension(:,:,:), allocatable, intent(out) :: new_data       ! New data set
!
character(len=PREC_STRING) :: zeile        ! Input line
character(len=PREC_STRING) :: line         ! Input line
character(len=6          ) :: befehl       ! The actual command
character(len=len(prompt)) :: orig_prompt
!
integer :: indxg
integer :: length            ! Input line length
integer :: lbef              ! Command length
integer :: lp                ! Line length
logical :: lend              ! if true terminate the menu
!
integer           , dimension(3)     :: new_icenter  ! Old center in pixels
logical           , dimension(3)     :: new_inc_user  = .false.  ! Try to determine inc's automatically?
logical           , dimension(4)     :: new_eck_user  = .false.  ! Try to determine eck's automatically?
logical           , dimension(3)     :: new_vi_user   = .false.  ! Try to determine vi's automatically?
logical                              :: new_zone_user = .false.  ! Try to determine zone automatically?
!
real(kind=PREC_DP), dimension(3)     :: new_zone     ! Zone axis in direct space
real(kind=PREC_DP), dimension(3)     :: new_absc     ! New abscissa in reciprocal space, normal to zone axis
real(kind=PREC_DP), dimension(3,3)   :: lib_tr_mat   ! Transformation matrix old => new
real(kind=PREC_DP), dimension(3,3)   :: lib_in_mat   ! Transformation matrix new => old
!logical :: success
!
call no_error
!
! Use default values 
call lib_trans_reset(outfile, lib_tr_mat, lib_in_mat, new_outfile,        &
     old_inc, old_eck, old_vi, new_inc, new_eck, new_vi, new_zone, new_icenter,     &
     new_inc_user, new_eck_user, new_vi_user, new_zone_user, new_data)
!
orig_prompt = prompt
prompt = prompt (1:len_str (prompt) ) //'/trans'
lend   = .FALSE.
!
loop_menu: do
   if(ier_num == 0) then
      CALL get_cmd (line, length, befehl, lbef, zeile, lp, prompt)  ! get new command
   endif
   if(lend) exit loop_menu
!
   IF(ier_num /= 0) THEN   ! Error in previous command or in get_cmd
      CALL errlist
      IF (ier_sta.ne.ER_S_LIVE) THEN
         IF (lmakro .OR. lmakro_error) THEN  ! Error within macro or termination errror
            IF(sprompt /= prompt ) THEN
               ier_num = -10
               ier_typ = ER_COMM
               ier_msg(1) = ' Error occured in data transformation menu'
               prompt_status = PROMPT_ON
               prompt = orig_prompt
               RETURN
            ELSE
               IF(lmacro_close) THEN
                  CALL macro_close
                  prompt_status = PROMPT_ON
               ENDIF
            ENDIF
         ENDIF
         IF (lblock) THEN
            ier_num = - 11
            ier_typ = ER_COMM
            prompt_status = PROMPT_ON
            prompt = orig_prompt
            RETURN
         ENDIF
         CALL no_error
         lmakro_error = .FALSE.
         sprompt = ' '
      ENDIF
      cycle loop_menu
   ENDIF
!
!  Normal menu commands
!
   IF (line == ' '      .or. line(1:1) == '#' .or. &
       line == char(13) .or. line(1:1) == '!'        ) cycle loop_menu
!
!     ----search for "="                                                
!                                                                       
   indxg = index (line, '=')
   IF (indxg.ne.0.AND..NOT. (str_comp (befehl, 'echo',   2, lbef, 4) ) &
                 .AND..NOT. (str_comp (befehl, 'system', 2, lbef, 6) )    &
                 .AND..NOT. (str_comp (befehl, 'help',   2, lbef, 4) .OR. &
                             str_comp (befehl, '?   ',   2, lbef, 4) )    &
                 .AND. INDEX(line,'==') == 0                            ) then
!                                                                       
!     ---evaluate an expression and assign the value to a variabble   
!                                                                       
      CALL do_math (line, indxg, length)
      cycle loop_menu
   endif
!
!     ----run transformation 'rese'                                            
!                                                                       
   if_rese: IF(str_comp (befehl, 'reset', 2, lbef, 5) ) then
      call lib_trans_reset(outfile, lib_tr_mat, lib_in_mat, new_outfile,        &
           old_inc, old_eck, old_vi, new_inc, new_eck, new_vi, new_zone, new_icenter,     &
           new_inc_user, new_eck_user, new_vi_user, new_zone_user, new_data)
      cycle loop_menu
   elseif(str_comp (befehl, 'exit', 2, lbef, 4) ) then  if_rese
      lend = .TRUE.
      exit loop_menu
   endif if_rese
!
!-- TRANSFORMATION commands
!
   if(str_comp(befehl, 'run', 3, lbef, 3)) then

      call lib_trans_main(idata, value, laver, old_inc, old_eck, old_vi, old_a0, old_win, &
           qvalues, VAL_PDF, VAL_3DPDF, ier_num, ier_typ, ER_IO, ER_APPL, &
           lib_tr_mat, lib_in_mat,                                              &
           new_outfile, new_inc, new_eck, new_vi, new_zone, new_icenter, &
           new_inc_user, new_eck_user, new_vi_user, new_zone_user,       &
           new_data)
   elseif(str_comp(befehl, 'set', 3, lbef, 3)) then
      call lib_trans_set(zeile, lp, new_outfile, new_inc, new_eck, new_vi,      &
                         new_zone, new_absc,                                    &
                         new_inc_user, new_eck_user, new_vi_user, new_zone_user)
!
! - Test common menu entries
!
   else
   CALL kdo_all(befehl, lbef, zeile, lp) 
      if(zeile == 'EXIT') then ! kdo_all detected "continue suite"
         lend = .TRUE.
      endif
   if(lend) exit loop_menu
!  if(success) cycle loop_menu
!  else
!     ier_num = -8
      ier_typ = ER_COMM
   endif
   if(lend) exit loop_menu
enddo loop_menu
!
prompt = orig_prompt
!
end subroutine lib_trans_menu
!
!*******************************************************************************
!
subroutine lib_trans_set(zeile, lp, new_outfile, new_inc, new_eck, new_vi,      &
                         new_zone, new_absc,                                    &
                         new_inc_user, new_eck_user, new_vi_user, new_zone_user)
!-
!  Set values for the transformation
!+
!
use errlist_mod
use get_params_mod
use str_comp_mod
use take_param_mod
!
implicit none
!
character(len=*)                  , intent(inout) :: zeile         ! input line
integer                           , intent(inout) :: lp            ! input length
character(len=*)                  , intent(out)   :: new_outfile   ! New output file
integer           , dimension(3)  , intent(out)   :: new_inc       ! New pixel numbers
real(kind=PREC_DP), dimension(3,4), intent(out)   :: new_eck       ! New Corners (ll, lr, ul, tl) 
real(kind=PREC_DP), dimension(3,3), intent(out)   :: new_vi        ! New vectors
real(kind=PREC_DP), dimension(3)  , intent(out)   :: new_zone      ! Zone axis in direct space
real(kind=PREC_DP), dimension(3)  , intent(out)   :: new_absc      ! New abscissa in reciprocal space, normal to zone axis
logical           , dimension(3)  , intent(out)   :: new_inc_user  ! Try to determine inc's automatically?
logical           , dimension(4)  , intent(out)   :: new_eck_user  ! Try to determine eck's automatically?
logical           , dimension(3)  , intent(out)   :: new_vi_user   ! Try to determine vi's automatically?
logical                           , intent(out)   :: new_zone_user ! Try to determine vi's automatically?
!
integer, parameter :: MAXP=4
!
integer :: ianz     ! Parameter number
character(len=PREC_STRING), dimension(MAXP) :: cpara   ! User parameters
integer                   , dimension(MAXP) :: lpara   ! Parameter length
!real(kind=PREC_DP)        , dimension(MAXP) :: werte   ! Numerical values
character(len=PREC_STRING)                  :: ccpara   ! User parameters
integer                                     :: llpara   ! Parameter length
real(kind=PREC_DP)        , dimension(MAXP) :: wwerte   ! Numerical values
!
integer, parameter :: NOPTIONAL = 10
integer, parameter :: O_ABS     =  1
integer, parameter :: O_ORD     =  2
integer, parameter :: O_TOP     =  3
integer, parameter :: O_OUTFILE =  4
integer, parameter :: O_NPOINTS =  5
integer, parameter :: O_LLB     =  6
integer, parameter :: O_RLB     =  7
integer, parameter :: O_LUB     =  8
integer, parameter :: O_LLT     =  9
integer, parameter :: O_ZONE    = 10
character(len=7)   , dimension(NOPTIONAL) :: oname   ! Optional parameter  name
character(len=MAX(PREC_STRING,LEN(zeile))), dimension(NOPTIONAL) :: opara   !Optional parameter strings returned
integer            , dimension(NOPTIONAL) :: loname  !Lenght opt. para name
integer            , dimension(NOPTIONAL) :: lopara  !Lenght opt. para name returned
logical            , dimension(NOPTIONAL) :: lpresent!opt. para is present
real(kind=PREC_DP) , dimension(NOPTIONAL) :: owerte   ! Calculated values
integer, parameter                        :: ncalc = 0 ! Number of values to calculate
data oname /'abs    ', 'ord    ', 'top    ', 'outfile', 'npoints',  &
            'll     ', 'lr     ', 'ul     ', 'tl     ', 'zone'    /
data loname/ 3       ,  3       ,  3       ,  7       ,  6       ,  &
             2       ,  2       ,  2       ,  2       ,  4        /
!
opara  = (/'auto     ', 'auto     ', 'auto     ', 'undefined', 'auto     ', &
           'auto     ', 'auto     ', 'auto     ', 'auto     ', 'auto     '  /)
lopara = (/ 4         ,  4         ,  4         ,  9         ,  4         , &
            4         ,  4         ,  4         ,  4         ,  4           /)
owerte = (/ 0.000D0   ,  0.000D0   ,  0.000D0   ,  0.000D0   ,  0.000D0   , &
            0.000D0   ,  0.000D0   ,  0.000D0   ,  0.000D0   ,  0.000D0     /)
!
call get_params(zeile, ianz, cpara, lpara, MAXP, lp)
if(ier_num/=0) return
!
call get_optional(ianz, MAXP, cpara, lpara, NOPTIONAL,  &
                  ncalc, oname, loname, opara, lopara, lpresent, owerte)
if(ier_num/=0) return
!
if(lpresent(O_NPOINTS)) then              ! User provided abscissa
   ccpara = opara(O_NPOINTS)
   llpara = lopara(O_NPOINTS)
   if(str_comp(ccpara, 'auto', 4, llpara, 4)) then
      new_inc_user(1) = .false.
   else
   call get_optional_multi(MAXP, ccpara, llpara, wwerte, ianz)
   if(ier_num/=0) return
      if(ianz==3) then
         new_inc(:) = nint(wwerte(1:3))
         new_inc_user(:) = .true.
      else
         ier_num = -6
         ier_typ = ER_COMM
         ier_msg(1) = 'Wrong number of data points'
         return
      endif
   endif
endif
!
if(lpresent(O_ABS)) then              ! User provided abscissa
   ccpara = opara(O_ABS)
   llpara = lopara(O_ABS)
   if(str_comp(ccpara, 'auto', 4, llpara, 4)) then
      new_vi_user(1) = .false.
   else
   call get_optional_multi(MAXP, ccpara, llpara, wwerte, ianz)
   if(ier_num/=0) return
      if(ianz==3) then
         new_vi(:,1) = wwerte(1:3)
         new_vi_user(1) = .true.
      else
         ier_num = -6
         ier_typ = ER_COMM
         ier_msg(1) = 'Wrong number of abscissa coordinates'
         return
      endif
   endif
endif
!
if(lpresent(O_ORD)) then              ! User provided abscissa
   ccpara = opara(O_ORD)
   llpara = lopara(O_ORD)
   if(str_comp(ccpara, 'auto', 4, llpara, 4)) then
      new_vi_user(2) = .false.
   else
   call get_optional_multi(MAXP, ccpara, llpara, wwerte, ianz)
   if(ier_num/=0) return
      if(ianz==3) then
         new_vi(:,2) = wwerte(1:3)
         new_vi_user(2) = .true.
      else
         ier_num = -6
         ier_typ = ER_COMM
         ier_msg(1) = 'Wrong number of ordinate coordinates'
         return
      endif
   endif
endif
!
if(lpresent(O_TOP)) then              ! User provided abscissa
   ccpara = opara(O_TOP)
   llpara = lopara(O_TOP)
   if(str_comp(ccpara, 'auto', 4, llpara, 4)) then
      new_vi_user(3) = .false.
   else
   call get_optional_multi(MAXP, ccpara, llpara, wwerte, ianz)
   if(ier_num/=0) return
      if(ianz==3) then
         new_vi(:,3) = wwerte(1:3)
         new_vi_user(3) = .true.
      else
         ier_num = -6
         ier_typ = ER_COMM
         ier_msg(1) = 'Wrong number of top axis coordinates'
         return
      endif
   endif
endif
!
if(lpresent(O_LLB)) then              ! User provided Left Lower Bottom
   ccpara = opara(O_LLB)
   llpara = lopara(O_LLB)
   if(str_comp(ccpara, 'auto', 4, llpara, 4)) then
      new_eck_user(1) = .false.
   else
      call get_optional_multi(MAXP, ccpara, llpara, wwerte, ianz)
      if(ier_num/=0) return
      if(ianz==3) then
         new_eck(:,1) = wwerte(1:3)
         new_eck_user(1) = .true.
      else
         ier_num = -6
         ier_typ = ER_COMM
         ier_msg(1) = 'Wrong number of LL corner coordinates'
         return
      endif
   endif
endif
!
if(lpresent(O_RLB)) then              ! User provided Right Lower Bottom = lr = lower right
   ccpara = opara(O_RLB)
   llpara = lopara(O_RLB)
   if(str_comp(ccpara, 'auto', 4, llpara, 4)) then
      new_eck_user(1) = .false.
   else
      call get_optional_multi(MAXP, ccpara, llpara, wwerte, ianz)
      if(ier_num/=0) return
      if(ianz==3) then
         new_eck(:,2) = wwerte(1:3)
         new_eck_user(2) = .true.
      else
         ier_num = -6
         ier_typ = ER_COMM
         ier_msg(1) = 'Wrong number of LR corner coordinates'
         return
      endif
   endif
endif
!
if(lpresent(O_LUB)) then              ! User provided Left Upper Bottom = ul = upper left 
   ccpara = opara(O_LUB)
   llpara = lopara(O_LUB)
   if(str_comp(ccpara, 'auto', 4, llpara, 4)) then
      new_eck_user(1) = .false.
   else
      call get_optional_multi(MAXP, ccpara, llpara, wwerte, ianz)
      if(ier_num/=0) return
      if(ianz==3) then
         new_eck(:,3) = wwerte(1:3)
         new_eck_user(3) = .true.
      else
         ier_num = -6
         ier_typ = ER_COMM
         ier_msg(1) = 'Wrong number of UL corner coordinates'
         return
      endif
   endif
endif
!
if(lpresent(O_LLT)) then              ! User provided Left Lower Top = tl = top left 
   ccpara = opara(O_LLT)
   llpara = lopara(O_LLT)
   if(str_comp(ccpara, 'auto', 4, llpara, 4)) then
      new_eck_user(1) = .false.
   else
      call get_optional_multi(MAXP, ccpara, llpara, wwerte, ianz)
      if(ier_num/=0) return
      if(ianz==3) then
         new_eck(:,4) = wwerte(1:3)
         new_eck_user(4) = .true.
      else
         ier_num = -6
         ier_typ = ER_COMM
         ier_msg(1) = 'Wrong number of TL corner coordinates'
         return
      endif
   endif
endif
!
!if(lpresent(O_INFILE)) then            ! User requested an input file
!   lib_tr_infile = opara(O_INFILE)
!endif
!
if(lpresent(O_OUTFILE)) then            ! User requested a new output file
   new_outfile = opara(O_OUTFILE)
endif
!
if(lpresent(O_ZONE)) then              ! User provided abscissa
   ccpara = opara(O_ZONE)
   llpara = lopara(O_ZONE)
   if(str_comp(ccpara, 'auto', 4, llpara, 4)) then
      new_zone_user = .false.
   else
   call get_optional_multi(MAXP, ccpara, llpara, wwerte, ianz)
   if(ier_num/=0) return
      if(ianz==3) then
         new_zone(:) = wwerte(1:3)
         new_zone_user = .true.
      else
         ier_num = -6
         ier_typ = ER_COMM
         ier_msg(1) = 'Wrong number of zone axis coordinates'
         return
      endif
   endif
endif
!
end subroutine lib_trans_set
!
!*******************************************************************************
!
subroutine lib_trans_main(idata, value, laver, old_inc, old_eck, old_vi, cr_a0, cr_win,     &
           qvalues, VAL_PDF, VAL_3DPDF, ier_num, ier_typ, ER_IO, ER_APPL,       &
           lib_tr_mat, lib_in_mat,                                              &
           new_outfile, new_inc, new_eck, new_vi, new_zone, new_icenter,        & 
           new_inc_user, new_eck_user, new_vi_user, new_zone_user,              &
           new_data)
!-
! Main interface to transform a 3D data set into a different orientation
!+
use gen_hdf_write_mod
use precision_mod
!
implicit none
!
integer                           , intent(in)  :: idata     ! Old data set number 
integer                           , intent(in)  :: value     ! Intensity, 3DPDF etc
logical                           , intent(in)  :: laver     ! Average value ; currently not needed
integer           , dimension(3)  , intent(in)  :: old_inc       ! Number data points  original
real(kind=PREC_DP), dimension(3,4), intent(in)  :: old_eck       ! The Corners (ll, lr, ul, tl)
real(kind=PREC_DP), dimension(3,3), intent(in)  :: old_vi        ! Old vectors vi(hkl,1-3)
real(kind=PREC_DP), dimension(3)  , intent(in)  :: cr_a0     ! Lattice parameters 
real(kind=PREC_DP), dimension(3)  , intent(in)  :: cr_win    ! Lattice angles 
real(kind=PREC_DP), dimension(old_inc(1),old_inc(2),old_inc(3)), intent(in) :: qvalues   ! Old data 
integer                           , intent(in)  :: VAL_PDF   ! Number for a PDF as output
integer                           , intent(in)  :: VAL_3DPDF ! Number for a 3DPDF as output
!real(kind=PREC_DP)                , intent(in)  :: valmax    ! Maximum value in original data
integer                           , intent(out) :: ier_num   ! Error number to report back
integer                           , intent(out) :: ier_typ   ! Error type   to report back
integer                           , intent(in)  :: ER_IO     ! Template to report an I/O error
integer                           , intent(in)  :: ER_APPL   ! Template to report an application error
!
real(kind=PREC_DP), dimension(3,3), intent(inout) :: lib_tr_mat   ! Transformation matrix old => new
real(kind=PREC_DP), dimension(3,3), intent(inout) :: lib_in_mat   ! Transformation matrix new => old
character(len=*)                  , intent(inout) :: new_outfile   ! Output file name passed down to hdf5_write
integer           , dimension(3)  , intent(out)   :: new_inc       ! New pixel numbers
real(kind=PREC_DP), dimension(3,3), intent(out)   :: new_vi        ! New increment vectors
real(kind=PREC_DP), dimension(3,4), intent(out)   :: new_eck       ! New Corners (ll, lr, ul, tl) 
real(kind=PREC_DP), dimension(3)  , intent(out)   :: new_zone      ! Zone axis in direct space
integer           , dimension(3)  , intent(out)   :: new_icenter   ! Number data points  original
logical           , dimension(3)  , intent(in)    :: new_inc_user  ! Try to determine inc's automatically?
logical           , dimension(4)  , intent(in)    :: new_eck_user  ! Try to determine eck's automatically?
logical           , dimension(3)  , intent(in)    :: new_vi_user   ! Try to determine vi's automatically?
logical                           , intent(in)    :: new_zone_user ! Try to determine vi's automatically?
real(kind=PREC_DP), dimension(:,:,:), allocatable, intent(out) :: new_data       ! New data set
!
real(kind=PREC_DP), dimension(3,3) :: temp_vi    ! Temporary copy
integer :: extr_abs, extr_ord, extr_top
!
!
integer           , dimension(3)     :: old_icenter            ! Old center in pixels
!integer           , dimension(3)     :: new_icenter            ! Old center in pixels
!
integer :: i,j,k   ! Indices in old data set
integer :: ii,jj,kk   ! Indices in old data set
!integer :: l,m,n   ! Indices in new data set
integer, dimension(3) :: lmn   ! Indices in new data set
real(kind=PREC_DP), dimension(3) :: vect   ! Pixel vector in old system
real(kind=PREC_DP), dimension(3) :: resl   ! Pixel vector in new system
real(kind=PREC_DP)               :: new_valmax   ! Pixel vector in new system
real(kind=PREC_DP), dimension(:,:,:), allocatable :: weight         ! Weights for pixels in new data
!
!new_vi = vz
!write(*,*) ' MAIN inc ', new_INC_user
!write(*,*) ' MAIN eck ', new_eck_user
!write(*,*) ' MAIN vi  ', new_vi_user
!write(*,*) ' MAIN zone', new_zone_user
call build_new(old_inc, old_eck, old_vi, old_icenter, cr_a0, cr_win, &
               lib_tr_mat, lib_in_mat,                 &
               new_inc, new_eck, new_vi, new_zone, new_icenter,  &
               new_inc_user, new_eck_user, new_vi_user, new_zone_user )
if(ier_num /= 0) return
!
allocate(new_data(new_inc(1), new_inc(2), new_inc(3)))
allocate(weight  (new_inc(1), new_inc(2), new_inc(3)))
new_data = 0.0
weight   = 0.0
!write(*,*) ' START TRANSFORMATION '
!write(*,*) ' old_inc ', old_inc
!write(*,*) ' old_icen', old_icenter
!write(*,*) ' Matrix  ', lib_tr_mat(:,1)
!write(*,*) ' Matrix  ', lib_tr_mat(:,2)
!write(*,*) ' Matrix  ', lib_tr_mat(:,3)
!write(*,*) ' new_inc ', new_inc
!write(*,*) ' new_icen', new_icenter
!write(*,*) ' new_LLB ', new_eck(:,1)
!write(*,*) ' new_RLB ', new_eck(:,2)
!write(*,*) ' new_LUB ', new_eck(:,3)
!write(*,*) ' new_LLT ', new_eck(:,4)
!write(*,*) ' new_abs ', new_vi (:,1)
!write(*,*) ' new_ord ', new_vi (:,2)
!write(*,*) ' new_top ', new_vi (:,3)
do kk=1, old_inc(3)
   k = kk - old_icenter(3)
   vect(3) = real(k,kind=PREC_DP)
   do jj=1, old_inc(2)
      j = jj - old_icenter(2)
      vect(2) = real(j,kind=PREC_DP)
      do ii=1, old_inc(1)
         i = ii - old_icenter(1)
         vect(1) = real(i,kind=PREC_DP)
         resl = matmul(lib_tr_mat, vect)
         lmn = int(resl(:)) + new_icenter(:)
         if(all(lmn>0) .and. all(lmn<=new_inc)) then
            new_data(lmn(1), lmn(2), lmn(3)) = new_data(lmn(1), lmn(2), lmn(3)) + &
                                               qvalues(ii,jj,kk)
            weight(lmn(1), lmn(2), lmn(3)) = weight(lmn(1), lmn(2), lmn(3)) + 1.0D0
         endif
      enddo
   enddo
enddo
where(weight>0.0D0)
   new_data = new_data/weight
end where
new_valmax = maxval(new_data)
!
!                        ! If idata >=0 do transformatioon in place, return 
!                        ! new_* to calling routine, everything will be deallocated 
!                        ! in the calling routine
if(idata < 0 ) then      ! Pure output data, write and deallocate
   temp_vi = new_vi
   extr_abs = maxloc(abs(temp_vi(:,1)), dim=1)
   temp_vi(extr_abs,2) = 0.0D0
   extr_ord = maxloc(abs(temp_vi(:,2)), dim=1)
   if(    extr_abs==1 .and. extr_ord==2) then
      extr_top = 3
   elseif(extr_abs==1 .and. extr_ord==3) then
      extr_top = 2
   elseif(extr_abs==2 .and. extr_ord==3) then
      extr_top = 1
   endif
CALL gen_hdf5_write (value, laver, new_outfile, new_inc, new_eck, new_vi, &
                     extr_abs, extr_ord, extr_top, &
                     cr_a0, cr_win, new_data,val_pdf, val_3Dpdf, new_valmax,  &
                     ier_num, ier_typ, ER_IO, ER_APPL)
!
deallocate(new_data)
endif
!
deallocate(weight)
!
end subroutine lib_trans_main
!
!*******************************************************************************
!
subroutine build_new(orig_old_inc, orig_old_eck, orig_old_vi, old_icenter, cr_a0, cr_win, &
                     lib_tr_mat, lib_in_mat,                &
                     new_inc, new_eck, new_vi, new_zone, new_icenter, &
                     orig_new_inc_user, orig_new_eck_user, orig_new_vi_user, orig_new_zone_user)
!-
! Calculate the new : Number data points
!                     Corners
!                     Increment vectors
!+
!
use errlist_mod
use lib_metric_mod
use matrix_mod
use precision_mod
!
implicit none
!
integer           , dimension(3)  , intent(in)    :: orig_old_inc  ! Old pixel numbers
real(kind=PREC_DP), dimension(3,4), intent(in)    :: orig_old_eck  ! Old Corners (ll, lr, ul, tl) 
real(kind=PREC_DP), dimension(3,3), intent(in)    :: orig_old_vi   ! Old increment vectors
integer           , dimension(3)  , intent(out)   :: old_icenter   ! Old center in pixels
real(kind=PREC_DP), dimension(3)  , intent(in)    :: cr_a0         ! Lattice parameters 
real(kind=PREC_DP), dimension(3)  , intent(in)    :: cr_win        ! Lattice angles 
real(kind=PREC_DP), dimension(3,3), intent(out)   :: lib_tr_mat    ! Transformation matrix old => new
real(kind=PREC_DP), dimension(3,3), intent(out)   :: lib_in_mat    ! Transformation matrix new => old
integer           , dimension(3)  , intent(inout) :: new_inc       ! New pixel numbers
real(kind=PREC_DP), dimension(3,3), intent(inout) :: new_vi        ! New increment vectors
real(kind=PREC_DP), dimension(3,4), intent(inout) :: new_eck       ! New Corners (ll, lr, ul, tl) 
real(kind=PREC_DP), dimension(3)  , intent(inout) :: new_zone      ! New zone axis
integer           , dimension(3)  , intent(out)   :: new_icenter   ! New center in pixels
logical           , dimension(3)  , intent(in)    :: orig_new_inc_user  ! Try to determine inc's automatically?
logical           , dimension(4)  , intent(in)    :: orig_new_eck_user  ! Try to determine eck's automatically?
logical           , dimension(3)  , intent(in)    :: orig_new_vi_user   ! Try to determine vi's automatically?
logical                           , intent(in)    :: orig_new_zone_user ! Try to determine vi's automatically?
!
real(kind=PREC_DP), parameter      :: TOL=1.0D-5   ! An uncertainty
!
character(len=3), dimension(3), parameter :: cname= (/ 'abs', 'ord', 'top'/)
integer               :: i
integer               :: ndims               ! Reciprocal space i 1-; 2-; 3-D
integer, dimension(3) :: points
!   Local copies of old_* This allows modification without changing the original
integer           , dimension(3)   :: old_inc       ! Old pixel numbers
real(kind=PREC_DP), dimension(3,3) :: old_vi        ! Old increment vectors
real(kind=PREC_DP), dimension(3,4) :: old_eck       ! Old Corners (ll, lr, ul, tl) 
logical           , dimension(3)   :: new_inc_user  ! Try to determine inc's automatically?
logical           , dimension(4)   :: new_eck_user  ! Try to determine eck's automatically?
logical           , dimension(3)   :: new_vi_user   ! Try to determine vi's automatically?
logical                            :: new_zone_user ! Try to determine vi's automatically?
!
real(kind=PREC_DP), dimension(3,3) :: gten   ! Direct metric tensor
real(kind=PREC_DP), dimension(3,3) :: rten   ! reciprocal metric tensor
real(kind=PREC_DP), dimension(3,3,3) :: eps  ! Direct space espilon tensor
real(kind=PREC_DP), dimension(3,3,3) :: reps  ! reciprocal space espilon tensor
real(kind=PREC_DP), dimension(3  ) :: new_zone_r   ! Direct metric tensor
real(kind=PREC_DP), dimension(3  ) :: new_zone_rn  ! Direct metric tensor
real(kind=PREC_DP)                 :: angle        ! An angle
!
!  Start by setting up a metric based on input data
!
call lib_tensor(gten, cr_a0, cr_win)
call matinv(gten, rten)
call lib_eps(gten, eps)
call lib_eps(rten, reps)
!
old_inc = orig_old_inc
old_eck = orig_old_eck
old_vi  = orig_old_vi
!
new_inc_user = orig_new_inc_user
new_eck_user = orig_new_eck_user
new_vi_user  = orig_new_vi_user
new_zone_user  = orig_new_zone_user
!
!write(*,*) ' INC U  ', new_inc_user
!write(*,*) ' ECK U  ', new_eck_user
!write(*,*) ' VI  U  ', new_vi_user
!write(*,*) ' ZON U  ', new_zone_user
ndims = 0                           ! Determine dimension
if(old_inc(1)>1) ndims=ndims+1
if(old_inc(2)>1) ndims=ndims+1
if(old_inc(3)>1) ndims=ndims+1
!write(*,*) ' DIMENSION ', ndims
call augment_dimension(ndims, old_inc, old_eck, old_vi,           &      ! Make dummy 3D
                              new_inc, new_eck, new_vi, new_zone, &  
                              gten, rten, eps, reps,              &
                       new_inc_user, new_eck_user, new_vi_user, new_zone_user)
!
!write(*,*) ' AUGMENT INC U  ', new_inc_user
!write(*,*) ' AUGMENT ECK U  ', new_eck_user
!write(*,*) ' AUGMENT VI  U  ', new_vi_user
!write(*,*) ' AUGMENT ZON U  ', new_zone_user
!
call error_check(ndims, new_inc_user, new_eck_user, new_vi_user, new_zone_user)
!
do i=1, 3
   if(mod(old_inc(i),2)==0) then
      old_icenter(i) = old_inc(i)/2
   else
      old_icenter(i) = (old_inc(i)+1)/2
   endif
enddo
!
if_zone:if(new_zone_user) then                  ! User provided zone axis
!write(*,*) ' MODE ZONE AXIS '
   call lib_d2r(gten, rten, new_zone, new_zone_r, new_zone_rn)  ! Transform Zone axis into reciprocal space
   call lib_vector_product(new_zone_rn, new_vi(:,1), new_vi(:,2), reps, gten)  ! Build direction of ordinate
   new_vi(:,2) = new_vi(:,2)/lib_blen(rten,new_vi(:,2))*lib_blen(rten,new_vi(:,1)) !Initially take, abs, ord, zone_rn
   new_vi(:,3) = new_zone_rn/lib_blen(rten,new_zone_rn)*lib_blen(rten,new_vi(:,1)) !Initially take, abs, ord, zone_rn
   call build_trans(old_vi, new_vi, lib_tr_mat, lib_in_mat) ! Build a temporary transformation matrix
!  Use build corners to determine temporary number of data points
   call build_corners(old_eck, old_inc, old_icenter, lib_tr_mat, lib_in_mat, &
                      points , new_eck, new_vi, new_icenter)
   if(all(new_inc_user)) then    ! User provided number of data points
      do i=1, 3                  ! Scale vectors by user increment points
!        if(new_inc(i)>old_inc(i)) then
            new_vi(:,i) = new_vi(:,i)*real(points (i),kind=PREC_DP)/real(new_inc(i),kind=PREC_DP)
!        endif
      enddo
   else                          ! Fully automatic
      do i=1, 3                  ! Scale vectors by old increment points
         new_vi(:,i) = new_vi(:,i)*real(points (i),kind=PREC_DP)/real(old_inc(i),kind=PREC_DP)
      enddo
   endif
   call build_trans(old_vi, new_vi, lib_tr_mat, lib_in_mat) ! Build a temporary transformation matrix
   call build_corners(old_eck, old_inc, old_icenter, lib_tr_mat, lib_in_mat, &
                      new_inc, new_eck, new_vi, new_icenter)
else  if_zone                                   ! No Zone axis mode
   if_corners:if(all(new_eck_user)) then         ! User provided all corners
!write(*,*) ' MODE Corners   '
      if(all(new_inc_user)) then      ! User provided all increment numbers
         if(new_inc(1)>1) new_vi(:,1) = (new_eck(:, 2) - new_eck(:,1))/real(new_inc(1)-1,kind=PREC_DP)   ! Abscissa vector
         if(new_inc(2)>1) new_vi(:,2) = (new_eck(:, 3) - new_eck(:,1))/real(new_inc(2)-1,kind=PREC_DP)   ! Ordinate vector
         if(new_inc(3)>1) new_vi(:,3) = (new_eck(:, 4) - new_eck(:,1))/real(new_inc(3)-1,kind=PREC_DP)   ! Top      vector
      elseif(all(new_vi_user)) then   ! User provided all vectors
!       Check if vectors are parallel and determine new_inc
         do i=1,3
            angle = lib_bang(rten, new_vi(:,i), (new_eck(:, i+1) - new_eck(:,1)))
            if(abs(angle)<TOL .or. abs(angle-180.0D0)<TOL) then
               points = 0
               where(new_vi(:,i) /= 0.0D0)
                  points(:) = nint(abs((new_eck(:, i+1) - new_eck(:,1))/new_vi(:,i)))
               end where
               new_inc(i) = maxval(points) + 1
               if(abs(angle-180.0D0)<TOL) then
                  new_vi(:,i) = -new_vi(:,i)
               endif
            else
               ier_num = -6
               ier_typ = ER_COMM
               ier_msg(1) = 'Increment vector '//cname(i)//' is not parallel to'
               ier_msg(3) = 'Vectors from  left lower bottom corner '
               return
            endif
         enddo
      else                            ! Use old increments as default
         new_vi(:,1) = (new_eck(:, 2) - new_eck(:,1))/real(old_inc(1)-1,kind=PREC_DP)   ! Abscissa vector
         new_vi(:,2) = (new_eck(:, 3) - new_eck(:,1))/real(old_inc(2)-1,kind=PREC_DP)   ! Ordinate vector
         new_vi(:,3) = (new_eck(:, 4) - new_eck(:,1))/real(old_inc(3)-1,kind=PREC_DP)   ! Top      vector
         new_inc = old_inc
      endif
      do i=1, 3
         if(mod(new_inc(i),2)==0) then
            new_icenter(i) = new_inc(i)/2
         else
            new_icenter(i) = (new_inc(i)+1)/2
         endif
      enddo
!
      call build_trans(old_vi, new_vi, lib_tr_mat, lib_in_mat)
   else if_corners                    ! User provided vectors
!      write(*,*) ' BUILD FROM vectors '
      call build_trans(old_vi, new_vi, lib_tr_mat, lib_in_mat)
      if(all(new_inc_user)) then      ! User provided number of data points as well, build corners
!write(*,*) ' CORN abs    ', new_vi(:,1), new_inc(1)
!write(*,*) ' CORN ord    ', new_vi(:,2), new_inc(2)
!write(*,*) ' CORN top    ', new_vi(:,3), new_inc(3)
!write(*,*)
!
         new_eck(:,1) = real(-1*(new_inc(1)-1)/2,kind=PREC_DP)*new_vi(:,1) +             &
                        real(-1*(new_inc(2)-1)/2,kind=PREC_DP)*new_vi(:,2) +             &
                        real(-1*(new_inc(3)-1)/2,kind=PREC_DP)*new_vi(:,3)
         new_eck(:,2) = real(   (new_inc(1)-1)/2,kind=PREC_DP)*new_vi(:,1) +             &
                        real(-1*(new_inc(2)-1)/2,kind=PREC_DP)*new_vi(:,2) +             &
                        real(-1*(new_inc(3)-1)/2,kind=PREC_DP)*new_vi(:,3)
         new_eck(:,3) = real(-1*(new_inc(1)-1)/2,kind=PREC_DP)*new_vi(:,1) +             &
                        real(   (new_inc(2)-1)/2,kind=PREC_DP)*new_vi(:,2) +             &
                        real(-1*(new_inc(3)-1)/2,kind=PREC_DP)*new_vi(:,3)
         new_eck(:,4) = real(-1*(new_inc(1)-1)/2,kind=PREC_DP)*new_vi(:,1) +             &
                        real(-1*(new_inc(2)-1)/2,kind=PREC_DP)*new_vi(:,2) +             &
                        real(   (new_inc(3)-1)/2,kind=PREC_DP)*new_vi(:,3)
         do i=1, 3
            if(mod(old_inc(i),2)==0) then
               new_icenter(i) = new_inc(i)/2
            else
               new_icenter(i) = (new_inc(i)+1)/2
            endif
         enddo
!
      elseif(.not.(all(new_inc_user) .or. all(new_eck_user))) then
         call build_corners(old_eck, old_inc, old_icenter, lib_tr_mat, lib_in_mat, &
                            new_inc, new_eck, new_vi, new_icenter)
      endif
   endif  if_corners
endif  if_zone
!
!write(*,*) ' INC U  ', new_inc_user
!write(*,*) ' ECK U  ', new_eck_user
!write(*,*) ' VI  U  ', new_vi_user
!write(*,*) ' ZON U  ', new_zone_user
!
!write(*,*) ' abs    ', new_vi(:,1)
!write(*,*) ' ord    ', new_vi(:,2)
!write(*,*) ' top    ', new_vi(:,3)
!write(*,*)
!
!write(*,*) ' LLB    ', new_eck(:,1)
!write(*,*) ' RLB    ', new_eck(:,2)
!write(*,*) ' LUB    ', new_eck(:,3)
!write(*,*) ' LLT    ', new_eck(:,4)
!write(*,*)
!
!write(*,*) ' ZONE   ', new_zone
!write(*,*)
!
!write(*,*) ' INC    ', new_inc
!write(*,*)
!write(*,*) ' TRANS  ', lib_tr_mat(:,1)
!write(*,*) ' TRANS  ', lib_tr_mat(:,2)
!write(*,*) ' TRANS  ', lib_tr_mat(:,3)
!write(*,*)
!
end subroutine build_new
!
!*******************************************************************************
!
subroutine augment_dimension(ndims, old_inc, old_eck, old_vi,                   &
                                    new_inc, new_eck, new_vi, new_zone,         &
                                    gten, rten, eps, reps,                      &
                       new_inc_user, new_eck_user, new_vi_user, new_zone_user)
!-
!  For 2-D and 1D data add missing coordinates 
!+
!
use lib_metric_mod
use precision_mod
!
implicit none
!
integer                           ,   intent(inout) :: ndims    ! Dimension
integer           , dimension(3)  ,   intent(inout) :: old_inc  ! Old pixel numbers
real(kind=PREC_DP), dimension(3,4),   intent(inout) :: old_eck  ! Old Corners (ll, lr, ul, tl) 
real(kind=PREC_DP), dimension(3,3),   intent(inout) :: old_vi   ! Old increment vectors
integer           , dimension(3)  ,   intent(inout) :: new_inc  ! New pixel numbers
real(kind=PREC_DP), dimension(3,4),   intent(inout) :: new_eck  ! New Corners (ll, lr, ul, tl) 
real(kind=PREC_DP), dimension(3,3),   intent(inout) :: new_vi   ! New increment vectors
real(kind=PREC_DP), dimension(3)  ,   intent(inout) :: new_zone ! New zone axis
real(kind=PREC_DP), dimension(3,3),   intent(in)    :: gten     ! Direct metric tensor
real(kind=PREC_DP), dimension(3,3),   intent(in)    :: rten     ! reciprocal metric tensor
real(kind=PREC_DP), dimension(3,3,3), intent(in)    :: eps      ! Direct metric tensor
real(kind=PREC_DP), dimension(3,3,3), intent(in)    :: reps     ! reciprocal metric tensor
logical           , dimension(3)    , intent(inout) :: new_inc_user  ! Try to determine inc's automatically?
logical           , dimension(4)    , intent(inout) :: new_eck_user  ! Try to determine eck's automatically?
logical           , dimension(3)    , intent(inout) :: new_vi_user   ! Try to determine vi's automatically?
logical                             , intent(inout) :: new_zone_user ! Try to determine vi's automatically?
!
real(kind=PREC_DP), dimension(3) :: u ! a vector
real(kind=PREC_DP), dimension(3) :: w ! a vector
real(kind=PREC_DP) :: r1 ! a random number
integer :: isflat
!
if(ndims==3) return          !  3D; nothing to do
!
if(ndims==2) then            ! 2D  
   isflat = minloc(old_inc, 1)
!   write(*,*) ' FLAT is ', isflat
   if(isflat==3) then
!      write(*,*) ' abs    ', old_vi(:,1)
!      write(*,*) ' ord    ', old_vi(:,2)
      call lib_vector_product(old_vi(:,1), old_vi(:,2), u, reps, gten)  ! Build direction of top axis
!      write(*,*) ' Normal ', u
      old_vi (:,3) = u
      old_eck(:,4) = old_eck(:,1) + u
      old_inc(3) = 1
      new_vi (:,3) = old_vi(:,3)
      new_eck(:,4) = old_eck(:,4)
      new_inc(3)   = 1
      new_zone     = u
!
      new_inc_user(3) = new_inc_user(2)   ! Copy user mode
      new_eck_user(4) = new_eck_user(3)   ! Copy user mode
      new_vi_user (3) = new_vi_user(2)    ! Copy user mode
   elseif(isflat==2) then
!      write(*,*) ' abs    ', old_vi(:,1)
!      write(*,*) ' top    ', old_vi(:,3)
      call lib_vector_product(old_vi(:,3), old_vi(:,1), u, reps, gten)  ! Build direction of ordinate
!      write(*,*) ' Normal ', u
      old_vi (:,2) = u
      old_eck(:,3) = old_eck(:,1) + u
      old_inc(2) = 1
      new_vi (:,2) = old_vi(:,2)
      new_eck(:,3) = old_eck(:,3)
      new_inc(2)   = 1
      new_zone     = u
!
      new_inc_user(2) = new_inc_user(3)   ! Copy user mode
      new_eck_user(3) = new_eck_user(4)   ! Copy user mode
      new_vi_user (2) = new_vi_user(3)    ! Copy user mode
   elseif(isflat==1) then
!      write(*,*) ' ord    ', old_vi(:,2)
!      write(*,*) ' top    ', old_vi(:,3)
      call lib_vector_product(old_vi(:,2), old_vi(:,3), u, reps, gten)  ! Build direction of ordinate
!      write(*,*) ' Normal ', u
      old_vi (:,1) = u
      old_eck(:,2) = old_eck(:,1) + u
      old_inc(1) = 1
      new_vi (:,1) = old_vi(:,1)
      new_eck(:,2) = old_eck(:,2)
      new_inc(1)   = 1
      new_zone     = u
!
      new_inc_user(1) = new_inc_user(2)   ! Copy user mode
      new_eck_user(2) = new_eck_user(3)   ! Copy user mode
      new_vi_user (1) = new_vi_user(2)    ! Copy user mode
   endif
elseif(ndims==1) then            ! 1D  
   isflat = maxloc(old_inc, 1)
!   write(*,*) ' FLAT is ', isflat
   if(isflat==1) then         ! Abscissa is non-flat
      call random_number(r1)
      w(1) = r1
      call random_number(r1)
      w(2) = r1
      call random_number(r1)
      w(3) = r1
      call lib_vector_product(old_vi(:,1), w, u, reps, gten)  ! Build a normal vector
      old_vi(:,2) = u
      old_eck(:,3) = old_eck(:,1) + u
      call lib_vector_product(old_vi(:,1), old_vi(:,2), u, reps, gten)  ! Build a normal vector
      old_vi(:,3) = u
      old_eck(:,4) = old_eck(:,1) + u
      old_inc(2) = 1
      old_inc(3) = 1
      new_inc(2) = 1
      new_inc(3) = 1
      new_eck(:,3) = old_eck(:,3)
      new_eck(:,4) = old_eck(:,4)
      new_vi (:,2) = old_vi (:,2)
      new_vi (:,3) = old_vi (:,3)
      new_inc_user(2:3) = new_inc_user(1)
      new_eck_user(3:4) = new_eck_user(1)
      new_vi_user (2:3) = new_vi_user(1)
      new_zone_user = .false.
   elseif(isflat==2) then     ! Ordinate is non-flat
      call random_number(r1)
      w(1) = r1
      call random_number(r1)
      w(2) = r1
      call random_number(r1)
      w(3) = r1
      call lib_vector_product(old_vi(:,2), w, u, reps, gten)  ! Build a normal vector
      old_vi(:,3) = u
      old_eck(:,4) = old_eck(:,1) + u
      call lib_vector_product(old_vi(:,2), old_vi(:,3), u, reps, gten)  ! Build a normal vector
      old_vi(:,1) = u
      old_eck(:,2) = old_eck(:,1) + u
      old_inc(3) = 1
      old_inc(1) = 1
      new_inc(3) = 1
      new_inc(1) = 1
      new_eck(:,4) = old_eck(:,4)
      new_eck(:,2) = old_eck(:,2)
      new_vi (:,3) = old_vi (:,3)
      new_vi (:,1) = old_vi (:,1)
      new_inc_user(1  ) = new_inc_user(2)
      new_inc_user(3  ) = new_inc_user(2)
      new_eck_user(2  ) = new_eck_user(3)
      new_eck_user(4  ) = new_eck_user(3)
      new_vi_user (1  ) = new_vi_user(2)
      new_vi_user (3  ) = new_vi_user(2)
      new_zone_user = .false.
   elseif(isflat==3) then     ! Top axis is non-flat
      call random_number(r1)
      w(1) = r1
      call random_number(r1)
      w(2) = r1
      call random_number(r1)
      w(3) = r1
      call lib_vector_product(old_vi(:,3), w, u, reps, gten)  ! Build a normal vector
      old_vi(:,1) = u
      old_eck(:,2) = old_eck(:,1) + u
      call lib_vector_product(old_vi(:,3), old_vi(:,1), u, reps, gten)  ! Build a normal vector
      old_vi(:,2) = u
      old_eck(:,3) = old_eck(:,1) + u
      old_inc(1) = 1
      old_inc(2) = 1
      new_inc(1) = 1
      new_inc(2) = 1
      new_eck(:,2) = old_eck(:,2)
      new_eck(:,3) = old_eck(:,3)
      new_vi (:,1) = old_vi (:,1)
      new_vi (:,2) = old_vi (:,2)
      new_inc_user(1  ) = new_inc_user(3)
      new_inc_user(2  ) = new_inc_user(3)
      new_eck_user(2  ) = new_eck_user(4)
      new_eck_user(3  ) = new_eck_user(4)
      new_vi_user (1  ) = new_vi_user(3)
      new_vi_user (2  ) = new_vi_user(3)
      new_zone_user = .false.
   endif
endif
!write(*,*) ' Augented  inc ', new_inc, new_inc_user
!write(*,*) ' Augented  eck ', new_eck(:,1), new_eck_user(1)
!write(*,*) ' Augented  eck ', new_eck(:,2), new_eck_user(2)
!write(*,*) ' Augented  eck ', new_eck(:,3), new_eck_user(3)
!write(*,*) ' Augented  eck ', new_eck(:,4), new_eck_user(4)
!write(*,*) ' Augented  inc ', new_vi, new_vi_user
!write(*,*) ' Augented  zon ', new_zone, new_zone_user
!
end subroutine augment_dimension
!
!*******************************************************************************
!
subroutine error_check(ndims, new_inc_user, new_eck_user, new_vi_user,   &
                       new_zone_user)
!-
! Perform error check on value completion and consistency
!+
!
use errlist_mod
!
implicit none
!
integer                           , intent(in)    :: ndims
logical           , dimension(3)  , intent(in)    :: new_inc_user  ! Try to determine inc's automatically?
logical           , dimension(4)  , intent(in)    :: new_eck_user  ! Try to determine eck's automatically?
logical           , dimension(3)  , intent(in)    :: new_vi_user   ! Try to determine vi's automatically?
logical                           , intent(in)    :: new_zone_user ! Try to determine vi's automatically?
!
if(.not.( all(new_eck_user) .or. (all(.not.new_eck_user)))) then
   ier_num = -6
   ier_typ = ER_COMM
   ier_msg(1) = 'Not all new corners were provided'
   return
endif
!
if(.not.((new_inc_user(1).eqv.new_inc_user(2)) .or. (new_inc_user(1).eqv.new_inc_user(3)) )) then
   ier_num = -6
   ier_typ = ER_COMM
   ier_msg(1) = 'Not all new number of points were provided'
   return
endif
!
if(.not.((new_vi_user(1).and. new_zone_user) .or.  &
         ((new_vi_user(1).eqv.new_vi_user(2)) .and. (new_vi_user(1).eqv.new_vi_user(3))) )) then
   ier_num = -6
   ier_typ = ER_COMM
   ier_msg(1) = 'Not all new increment vectors were provided'
   return
endif
!
if(new_zone_user .and. .not. new_vi_user(1)) then
   ier_num = -6
   ier_typ = ER_COMM
   ier_msg(1) = 'Zone axis mode, but abscissa is missing'
   return
endif
!
end subroutine error_check
!
!*******************************************************************************
!
subroutine build_trans(vi, vz, lib_tr_mat, lib_in_mat)
!-
!  Calculate the transformation matrix from old vi-vectors to new vectors
!
!       ( A, O, T )
!  vi = ( B, R, O )  Vectors are columns in matrix vi
!       ( S, D, P )
!
!  vz = TR*vi;   TR = zz * vi^-1
!+! 
!
use matrix_mod
use precision_mod
!
implicit none
!
real(kind=PREC_DP), dimension(3,3), intent(in)  :: vi   ! Old vectors vi(hkl,1-3)
real(kind=PREC_DP), dimension(3,3), intent(in)  :: vz   ! New vectors vi(hkl,1-3)
real(kind=PREC_DP), dimension(3,3), intent(out) :: lib_tr_mat   ! Transformation matrix old => new
real(kind=PREC_DP), dimension(3,3), intent(out) :: lib_in_mat   ! Transformation matrix new => old
!
real(kind=PREC_DP), dimension(3,3) :: temp   ! Temporary matrix 
!
call matinv(vz, temp)
lib_tr_mat = matmul(temp, vi)
call matinv(lib_tr_mat, lib_in_mat)
!
end subroutine build_trans
!
!*******************************************************************************
!
subroutine build_corners(eck, inc, old_icenter, lib_tr_mat, lib_in_mat, &
                         new_inc, new_eck, new_vi, new_icenter)
!-
! Build metric of the eight original corners 
! and their coordinates as multiples of the new axes
! Corners with left/right ; lower/upper ; bottom/top
!    bottom                 top 
! lub=3 -- rub=4       lut=7 -- rut=8
!  |        |           |        |
! llb=1 -- rlb=2       llt=5 -- rlt=6
!+
!
use precision_mod
!
implicit none
!
real(kind=PREC_DP), dimension(3,4), intent(in)   :: eck          ! The old Corners (ll, lr, ul, tl)
integer           , dimension(3)  , intent(in)   :: inc          ! Old Number data points  original
integer           , dimension(3)  , intent(out)  :: old_icenter  ! Old center in pixels
real(kind=PREC_DP), dimension(3,3), intent(in)   :: lib_tr_mat   ! Transformation matrix old => new
real(kind=PREC_DP), dimension(3,3), intent(in)   :: lib_in_mat   ! Transformation matrix new => old
integer           , dimension(3)  , intent(out)  :: new_inc      ! New Number data points  original
real(kind=PREC_DP), dimension(3,4), intent(out)  :: new_eck      ! The new corners (ll, lr, ul, tl)
real(kind=PREC_DP), dimension(3,3)  , intent(in) :: new_vi       ! New increment vectors
integer           , dimension(3)  , intent(out)  :: new_icenter  ! Old center in pixels
!
integer :: i, j   ! Dummy index
integer           , dimension(3,1:8) :: lib_iocorners ! Old corners in pixels away from center
integer           , dimension(3,1:8) :: lib_incorners ! Old corners in pixels away from center
real(kind=PREC_DP), dimension(3,1:8) :: lib_ocorners ! Old corners
real(kind=PREC_DP), dimension(3,1:8) :: lib_ncorners ! Old corners
integer           , dimension(3,1:8) :: new_icorners ! New corners in pixels away from center
real(kind=PREC_DP), dimension(3)     :: vect   ! Temporary vector 
real(kind=PREC_DP), dimension(3)     :: resl   ! Temporary vector 
!
lib_ocorners(:,1:3) =          eck(:,1:3)
lib_ocorners(:,4)   = lib_ocorners(:,2) + (eck(:,3) - eck(:,1))  ! right upper bottom
lib_ocorners(:,5)   =          eck(:,4)                          ! left  lower top
lib_ocorners(:,6)   = lib_ocorners(:,5) + (eck(:,2) - eck(:,1))  ! right lower top
lib_ocorners(:,7)   = lib_ocorners(:,5) + (eck(:,3) - eck(:,1))  ! left  upper top
lib_ocorners(:,8)   = lib_ocorners(:,4) + (eck(:,4) - eck(:,1))  ! right upper top
!
! Corners with respect to new coordinate system
do i=1, 8
  lib_ncorners(:,i) = matmul(lib_tr_mat, lib_ocorners(:,i))
enddo
!
do i=1, 3
   if(mod(inc(i),2)==0) then
      old_icenter(i) = inc(i)/2
   else
      old_icenter(i) = (inc(i)+1)/2
   endif
enddo
!
lib_iocorners(:, 1) =  (/1     -old_icenter(1), 1     -old_icenter(2), 1     -old_icenter(3)/)  ! left  lower bottom
lib_iocorners(:, 2) =  (/inc(1)-old_icenter(1), 1     -old_icenter(2), 1     -old_icenter(3)/)  ! right lower bottom
lib_iocorners(:, 3) =  (/1     -old_icenter(1), inc(2)-old_icenter(2), 1     -old_icenter(3)/)  ! left  upper bottom
lib_iocorners(:, 4) =  (/inc(1)-old_icenter(1), inc(2)-old_icenter(2), 1     -old_icenter(3)/)  ! right upper bottom
lib_iocorners(:, 5) =  (/1     -old_icenter(1), 1     -old_icenter(2), inc(3)-old_icenter(3)/)  ! left  lower top
lib_iocorners(:, 6) =  (/inc(1)-old_icenter(1), 1     -old_icenter(2), inc(3)-old_icenter(3)/)  ! right lower top
lib_iocorners(:, 7) =  (/1     -old_icenter(1), inc(2)-old_icenter(2), inc(3)-old_icenter(3)/)  ! left  upper top
lib_iocorners(:, 8) =  (/inc(1)-old_icenter(1), inc(2)-old_icenter(2), inc(3)-old_icenter(3)/)  ! right upper top
!
do i=1, 8
   vect = lib_iocorners(:,i)
   resl = matmul(lib_tr_mat, vect)
   lib_incorners(:,i) = resl
enddo
!
new_inc(1) = maxval(lib_incorners(1,:))-minval(lib_incorners(1,:)) + 1
new_inc(2) = maxval(lib_incorners(2,:))-minval(lib_incorners(2,:)) + 1
new_inc(3) = maxval(lib_incorners(3,:))-minval(lib_incorners(3,:)) + 1
!
do i=1, 3
   if(mod(new_inc(i),2)==0) then
      new_icenter(i) = new_inc(i)/2
   else
      new_icenter(i) = (new_inc(i)+1)/2
   endif
enddo
!
new_icorners(:, 1) =  (/1         - new_icenter(1), 1         - new_icenter(2), 1         - new_icenter(3)/)  ! left  lower bottom
new_icorners(:, 2) =  (/new_inc(1)- new_icenter(1), 1         - new_icenter(2), 1         - new_icenter(3)/)  ! right lower bottom
new_icorners(:, 3) =  (/1         - new_icenter(1), new_inc(2)- new_icenter(2), 1         - new_icenter(3)/)  ! left  upper bottom
new_icorners(:, 4) =  (/new_inc(1)- new_icenter(1), new_inc(2)- new_icenter(2), 1         - new_icenter(3)/)  ! right upper bottom
new_icorners(:, 5) =  (/1         - new_icenter(1), 1         - new_icenter(2), new_inc(3)- new_icenter(3)/)  ! left  lower top
new_icorners(:, 6) =  (/new_inc(1)- new_icenter(1), 1         - new_icenter(2), new_inc(3)- new_icenter(3)/)  ! right lower top
new_icorners(:, 7) =  (/1         - new_icenter(1), new_inc(2)- new_icenter(2), new_inc(3)- new_icenter(3)/)  ! left  upper top
new_icorners(:, 8) =  (/new_inc(1)- new_icenter(1), new_inc(2)- new_icenter(2), new_inc(3)- new_icenter(3)/)  ! right upper top
!
do j = 1, 3 ! First three corners are OK in this sequence
do i = 1, 3
   new_eck(i,j) = new_vi(i,1)*new_icorners(1,j) + new_vi(i,2)*new_icorners(2,j) + new_vi(i,3)*new_icorners(3,j)
enddo
enddo
do i = 1, 3
   new_eck(i,4) = new_vi(i,1)*new_icorners(1,5) + new_vi(i,2)*new_icorners(2,5) + new_vi(i,3)*new_icorners(3,5)  ! left lower top
enddo
!
end subroutine build_corners
!
!*******************************************************************************
!
pure function get_length(vec, gten) result(length)
!
use precision_mod
!
implicit none
!
real(kind=PREC_DP) :: length
!
real(kind=PREC_DP), dimension(3)  , intent(in) :: vec       ! A vector
real(kind=PREC_DP), dimension(3,3), intent(in) :: gten      ! Metric tensor
!
real(kind=PREC_DP), dimension(1,3) :: vec_t     ! vector "vec" transposed as row(1x3)
real(kind=PREC_DP), dimension(3,1) :: vec_u     ! vector "vec" as collumn (3x1)
real(kind=PREC_DP), dimension(1,1) :: point     ! A dummy matrix (1x1)
!
vec_t(1,:) = vec(:)
!
vec_u(:,1) = matmul(gten, vec)
point  =      matmul(vec_t, vec_u) 
length = sqrt(point(1,1))
!
end function get_length
!
!*******************************************************************************
!
subroutine lib_trans_reset (outfile, lib_tr_mat, lib_in_mat, new_outfile,       &
           old_inc, old_eck, old_vi, &
           new_inc, new_eck, new_vi, new_zone, new_icenter, new_inc_user, new_eck_user, new_vi_user, new_zone_user,   &
           new_data)
!-
! Reset all variables to system start
!+
!
use precision_mod
!
implicit none
!
character(len=*)                  , intent(in ) :: outfile      ! Output file name passed down to hdf5_write
integer           , dimension(3)  , intent(in ) :: old_inc      ! Old Number data points  original
real(kind=PREC_DP), dimension(3,4), intent(in ) :: old_eck      ! The old corners (ll, lr, ul, tl)
real(kind=PREC_DP), dimension(3,3), intent(in ) :: old_vi       ! Old increment vectors
real(kind=PREC_DP), dimension(3,3), intent(out) :: lib_tr_mat   ! Transformation matrix old => new
real(kind=PREC_DP), dimension(3,3), intent(out) :: lib_in_mat   ! Transformation matrix new => old
character(len=*)                  , intent(out) :: new_outfile  ! Output file name passed down to hdf5_write
integer           , dimension(3)  , intent(out) :: new_inc      ! New Number data points  original
real(kind=PREC_DP), dimension(3,4), intent(out) :: new_eck      ! The new corners (ll, lr, ul, tl)
real(kind=PREC_DP), dimension(3,3), intent(out) :: new_vi       ! New increment vectors
real(kind=PREC_DP), dimension(3)  , intent(out) :: new_zone      ! Zone axis in direct space
integer           , dimension(3)  , intent(out) :: new_icenter  ! Old center in pixels
logical           , dimension(3)  , intent(out) :: new_inc_user  ! Try to determine inc's automatically?
logical           , dimension(4)  , intent(out) :: new_eck_user  ! Try to determine eck's automatically?
logical           , dimension(3)  , intent(out) :: new_vi_user   ! Try to determine vi's automatically?
logical                           , intent(out) :: new_zone_user ! Try to determine vi's automatically?
real(kind=PREC_DP), dimension(:,:,:), allocatable, intent(inout) :: new_data       ! New data set
!
integer :: i
real(kind=PREC_DP), dimension(3,3), parameter :: unit_mat = reshape((/ 1.0D0, 0.0D0, 0.0D0, &
                                                                       0.0D0, 1.0D0, 0.0D0, &
                                                                       0.0D0, 0.0D0, 1.0D0 /), shape(unit_mat))
REAL(kind=PREC_DP), DIMENSION(1:3, 1:4) ::  eck      =      reshape((/ 0.0, 0.0,  0.0, &    ! (hkl, corner_number)
                                                                       5.0, 0.0,  0.0, &
                                                                       0.0, 5.0,  0.0, &
                                                                       0.0, 0.0,  0.0/),shape(eck))
!
lib_tr_mat = unit_mat
lib_in_mat = unit_mat
new_outfile = outfile
new_inc = old_inc
new_eck = old_eck
new_vi  = old_vi
new_zone = 0.0D0
do i=1, 3
   if(mod(new_inc(i),2)==0) then
      new_icenter(i) = new_inc(i)/2
   else
      new_icenter(i) = (new_inc(i)+1)/2
   endif
enddo
new_inc_user  = .false.
new_eck_user  = .false.
new_vi_user   = .false.
new_zone_user = .false.
if(allocated(new_data)) deallocate(new_data)
!
end subroutine lib_trans_reset
!
!*******************************************************************************
!
end module lib_trans_mod
