module kuplot_adt_mod
!-
! Routines to read the MAINZ ADT data
!+
use precision_mod
!
private
public adt_menu
!
integer, parameter         :: ADT_RODS  = 1
integer, parameter         :: ADT_PLANE = 2
integer, parameter         :: ADT_SPACE = 3
integer, parameter         :: ADT_BRAGG = 4
character(len=PREC_STRING) :: adt_data        ! Input MRC file
character(len=PREC_STRING) :: adt_cell        ! Input MRC cell information
character(len=PREC_STRING) :: adt_logfile     ! Output file for extraction
logical                            :: adt_log  = .false. ! Extract log on / off
integer                            :: adt_mode = 0! Extraction mode rods, plane, bulk
integer                            :: adt_node = 0! Node in global data structure
integer                            :: adt_ik   = 0! AData set in kuplot
real(kind=PREC_DP), dimension(3,3) :: adt_omat    ! Orientation matrix reciprocal space in pixels
real(kind=PREC_DP), dimension(3,3) :: adt_omat_t  ! Orientation matrix reciprocal space in pixels TRANSPOSED
real(kind=PREC_DP), dimension(3,3) :: adt_invmat  ! Inverse orientation matrix
real(kind=PREC_DP), dimension(3,3) :: adt_invmat_t! Inverse orientation matrix TRANSPOSED
real(kind=PREC_DP), dimension(3,3) :: adt_gten    ! Metric tensor
real(kind=PREC_DP), dimension(3,3) :: adt_rten    ! Reciprocal Metric tensor
real(kind=PREC_DP), dimension(6)   :: adt_px_cell ! Reciprocel cell in Pixels, degrees
real(kind=PREC_DP), dimension(6)   :: adt_di_cell ! Direct     cell in Angstr, degrees
real(kind=PREC_DP), dimension(6)   :: adt_re_cell ! Reciprocal cell in Angstr, degrees
real(kind=PREC_DP), dimension(3)   :: adt_limits  ! Upper limits in reciprocal space
real(kind=PREC_DP), dimension(3)   :: adt_rod     ! Rod directions to extract
real(kind=PREC_DP), dimension(3)   :: adt_sigma   ! Sigma in reciprocal units for extraction
real(kind=PREC_DP), dimension(3)   :: adt_steps   ! Steps in reciprocal units for extraction
real(kind=PREC_DP), dimension(3)   :: adt_zero    ! Pixel location of zero point
real(kind=PREC_DP)                 :: adt_resol   ! Resolution A^-1 pro pixel
!
contains
!
!*******************************************************************************
!
subroutine adt_menu
!-
!   Main menu for ADT
!+
!
use kuplot_kdo_common_mod, only:kuplot_kdo_common
!
use calc_expr_mod        , only:do_math
use class_macro_internal
use errlist_mod
use doact_mod
use lib_errlist_func
use lib_macro_func       , only:macro_close
use macro_mod
use prompt_mod
use str_comp_mod         , only:str_comp
use sup_mod              , only:get_cmd
!
implicit none
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
logical :: success
!
call no_error
!
orig_prompt = prompt
prompt = prompt (1:len_trim(prompt) ) //'/adt'
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
               ier_msg(1) = ' Error occured in adt menu'
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
   IF (indxg.ne.0.AND..NOT. (str_comp (befehl, 'echo', 2, lbef, 4) ) &
                 .AND..NOT. (str_comp (befehl, 'syst', 2, lbef, 4) )    &
                 .AND..NOT. (str_comp (befehl, 'help', 2, lbef, 4) .OR. &
                             str_comp (befehl, '?   ', 2, lbef, 4) )    &
                 .AND. INDEX(line,'==') == 0                            ) THEN
!                                                                       
!     ---evaluate an expression and assign the value to a variabble   
!                                                                       
      CALL do_math (line, indxg, length)
      cycle loop_menu
   endif
!
! - Test common menu entries
!
   CALL kuplot_kdo_common(befehl, lbef, line, length, zeile, lp, 'adt' , &
                          lend, success)
   if(lend) exit loop_menu
   if(success) cycle loop_menu

!                                                                       
!     ----run adt 'rese'                                            
!                                                                       
   if_rese: IF(str_comp (befehl, 'rese', 2, lbef, 4) ) THEN
      call adt_reset
      cycle loop_menu
   endif if_rese
!
!  if(str_comp (befehl, 'run', 2, lbef, 3) ) then
!     call adt_run
       if(str_comp (befehl, 'get', 2, lbef, 3) ) then
      call adt_get(zeile, lp)
   elseif(str_comp (befehl, 'show', 2, lbef, 4) ) then
      call adt_show
   elseif(str_comp (befehl, 'set', 2, lbef, 3) ) then
      call adt_set(zeile, lp)
   elseif(str_comp (befehl, 'trans', 2, lbef, 5) ) then
      call adt_trans(zeile, lp)
   elseif(str_comp (befehl, 'extract', 3, lbef, 7) ) then
      call adt_extract
   else
      ier_num = -8
      ier_typ = ER_COMM
   endif
   if(lend) exit loop_menu
!
enddo loop_menu
!
prompt = orig_prompt
!
end subroutine adt_menu
!
!*******************************************************************************
!
subroutine adt_show
!
use prompt_mod
!
implicit none
!
write(output_io,'(a)')   ' ADT SHOW '
write(output_io,'(a,a)') ' ADT DATA : ', adt_data(1:len_trim(adt_data))
write(output_io,'(a,a)') ' ADT CELL : ', adt_cell(1:len_trim(adt_cell))
write(output_io,'(a,15x,a,3x,a,13x,a)')    ' Orientation','Vectors in pixels','|','Inverse matrix'
write(output_io,'(a,2(2x,3(2x,f10.4)))') ' a*            :', adt_omat(1,:),adt_invmat(1,:)
write(output_io,'(a,2(2x,3(2x,f10.4)))') ' b*            :', adt_omat(2,:),adt_invmat(2,:)
write(output_io,'(a,2(2x,3(2x,f10.4)))') ' c*            :', adt_omat(3,:),adt_invmat(3,:)
write(output_io,'(a,15x,a,3x,a,13x,a)')    ' TRANSPOSED ','Vectors in pixels','|','Inverse matrix'
write(output_io,'(a,2(2x,3(2x,f10.4)))') ' a*            :', adt_omat_t(1,:),adt_invmat_t(1,:)
write(output_io,'(a,2(2x,3(2x,f10.4)))') ' b*            :', adt_omat_t(2,:),adt_invmat_t(2,:)
write(output_io,'(a,2(2x,3(2x,f10.4)))') ' c*            :', adt_omat_t(3,:),adt_invmat_t(3,:)
if(adt_mode==ADT_RODS) then
  write(output_io, *) 
  write(output_io, '(a)') ' Extracting rods in reciprocal space'
  write(output_io, '(a, 2x,3(2x,f10.4))') ' Rods parallel :', adt_rod
elseif(adt_mode==ADT_BRAGG) then
  write(output_io, '(a)') ' Extracting integrated Bragg data'
endif
  write(output_io, '(a, 2x,3(2x,f10.4))') ' Limits        :', adt_limits
  write(output_io, '(a, 2x,3(2x,f10.4))') ' steps         :', adt_steps
  write(output_io, '(a, 2x,3(2x,f10.4))') ' sigmas        :', adt_sigma
!
end subroutine adt_show
!
!*******************************************************************************
!
subroutine adt_get(zeile, length)
!
!
use ber_params_mod , only:ber_params
use errlist_mod
use get_params_mod , only:get_params
use precision_mod
use str_comp_mod   , only:str_comp
use take_param_mod
!
implicit none
!
character(len=*), intent(inout) :: zeile         ! input line
integer         , intent(inout) :: length        ! input length
!
!integer, parameter :: MIN_PARA = 3               ! Minimum parameter numbers
integer, parameter :: MAXW = 4                   ! Actual parameter number
!
character(len=PREC_STRING), dimension(MAXW) :: cpara
integer                   , dimension(MAXW) :: lpara
real(kind=PREC_DP)        , dimension(MAXW) :: werte
integer :: ianz          ! Number of parameters
!
integer, parameter :: NOPTIONAL = 2
integer, parameter :: O_DATA    = 1
integer, parameter :: O_CELL    = 2
character(LEN=   4), dimension(NOPTIONAL) :: oname   !Optional parameter names
character(LEN=PREC_STRING), dimension(NOPTIONAL) :: opara   !Optional parameter strings returned
integer            , dimension(NOPTIONAL) :: loname  !Lenght opt. para name
integer            , dimension(NOPTIONAL) :: lopara  !Lenght opt. para name returned
logical            , dimension(NOPTIONAL) :: lpresent!opt. para is present
real(kind=PREC_DP) , dimension(NOPTIONAL) :: owerte   ! Calculated values
integer, parameter                        :: ncalc = 0 ! Number of values to calculate 
!
data oname  / 'data', 'cell' /
data loname /  4    ,  4     /
opara  =  (/ 'unknown','unknown' /)   ! Always provide fresh default values
lopara =  (/  7,        7        /)
owerte =  (/  0.0,      0.0      /)
!
!
!
call get_params(zeile, ianz, cpara, lpara, maxw, length)
if(ier_num/=0) return
call get_optional(ianz, MAXW, cpara, lpara, NOPTIONAL,  ncalc, &
                  oname, loname, opara, lopara, lpresent, owerte)
if(ier_num/=0) return
!
if(lpresent(O_DATA)) then
   adt_data = opara(O_DATA)(1:lopara(O_DATA))
   if(adt_data /= 'unknown') call adt_read_data
endif
!
if(lpresent(O_CELL)) then
   adt_cell = opara(O_CELL)(1:lopara(O_CELL))
   if(adt_cell /= 'unknown') call adt_read_cell
endif
!
end subroutine adt_get
!
!*******************************************************************************
!
subroutine adt_set(zeile, length)
!
!
use ber_params_mod , only:ber_params
use errlist_mod
use get_params_mod , only:get_params
use precision_mod
use str_comp_mod   , only:str_comp
use take_param_mod
!
implicit none
!
character(len=*), intent(inout) :: zeile         ! input line
integer         , intent(inout) :: length        ! input length
!
!integer, parameter :: MIN_PARA = 3               ! Minimum parameter numbers
integer, parameter :: MAXW = 4                   ! Actual parameter number
!
character(len=PREC_STRING), dimension(MAXW) :: cpara
integer                   , dimension(MAXW) :: lpara
real(kind=PREC_DP)        , dimension(MAXW) :: werte
integer :: ianz          ! Number of parameters
!
integer, parameter :: NOPTIONAL = 7
integer, parameter :: O_MODE    = 1
integer, parameter :: O_LIMITS  = 2
integer, parameter :: O_STEPS   = 3
integer, parameter :: O_SIGMA   = 4
integer, parameter :: O_ZERO    = 5
integer, parameter :: O_RODS    = 6
integer, parameter :: O_LOGFILE = 7
character(LEN=   6), dimension(NOPTIONAL) :: oname   !Optional parameter names
character(LEN=PREC_STRING), dimension(NOPTIONAL) :: opara   !Optional parameter strings returned
integer            , dimension(NOPTIONAL) :: loname  !Lenght opt. para name
integer            , dimension(NOPTIONAL) :: lopara  !Lenght opt. para name returned
logical            , dimension(NOPTIONAL) :: lpresent!opt. para is present
real(kind=PREC_DP) , dimension(NOPTIONAL) :: owerte   ! Calculated values
integer, parameter                        :: ncalc = 0 ! Number of values to calculate 
!
data oname  / 'mode', 'limits', 'steps', 'sigma', 'zero', 'rod', 'log' /
data loname /  4    ,  6      ,  4     ,  5     ,  4    ,  3   ,  3 /
opara  =  (/ 'rods         ','[1.0,1.0,1.0]', '[1.0,1.0,0.1]', '[0.2,0.2,0.0]', &
             '[-1.,-1.,-1.]','[0.0,0.0,1.0]' ,'bragg.data   '/) !Provide fresh default values
lopara =  (/  4             , 13            ,  13            ,  13            , &
              13            , 13            , 10             /)
owerte =  (/  0.0           , 0.0           , 0.0           ,  0.0           , &
              0.0           , 0.0           , 0.0            /)
!
!
call get_params(zeile, ianz, cpara, lpara, maxw, length)
if(ier_num/=0) return
call get_optional(ianz, MAXW, cpara, lpara, NOPTIONAL,  ncalc, &
                  oname, loname, opara, lopara, lpresent, owerte)
if(ier_num/=0) return
!
if(lpresent(O_MODE)) then
   if(opara(O_MODE)=='rods') then
      adt_mode = ADT_RODS
   elseif(opara(O_MODE)=='bragg') then
      adt_mode = ADT_BRAGG
   endif
   if(lpresent(O_LOGFILE)) then
      adt_logfile = opara(O_LOGFILE)
      adt_log = .true.
   else
      adt_log = .false.
   endif
write(*,*) ' LOG ', adt_log, adt_logfile
endif
!
if(lpresent(O_LIMITS)) then
   call get_optional_multi(MAXW, opara(O_LIMITS), lopara(O_LIMITS), werte, ianz)
   adt_limits = werte(1:3)
endif
!
if(lpresent(O_SIGMA )) then
   call get_optional_multi(MAXW, opara(O_SIGMA ), lopara(O_SIGMA ), werte, ianz)
   adt_sigma  = werte(1:3)
endif
!
if(lpresent(O_STEPS )) then
   call get_optional_multi(MAXW, opara(O_STEPS ), lopara(O_STEPS ), werte, ianz)
   adt_steps  = werte(1:3)
endif
!
if(lpresent(O_ZERO )) then
   call get_optional_multi(MAXW, opara(O_ZERO ), lopara(O_ZERO ), werte, ianz)
   adt_zero   = werte(1:3)
endif
!
if(lpresent(O_RODS )) then
   call get_optional_multi(MAXW, opara(O_RODS ), lopara(O_RODS ), werte, ianz)
   adt_rod    = werte(1:3)
endif
!
if(ier_num/=0) return
!
end subroutine adt_set
!
!*******************************************************************************
!
subroutine adt_read_data
!
use kuplot_mod
use kuplot_load_mod, only:mrc_read_kuplot
!
use errlist_mod
use param_mod
!
implicit none
!
integer :: length
!
if(adt_data == ' ') then
   ier_num = -6
   ier_typ = ER_COMM
   ier_msg(1) = 'ADT data file name is blank'
   return
endif
!
if(adt_data(len_trim(adt_data)-3:len_trim(adt_data)) /= '.mrc') then
   adt_data(len_trim(adt_data)+1:len_trim(adt_data)+4) = '.mrc'
endif
!
length = len_trim(adt_data)
write(*,*) ' iz ', iz
call mrc_read_kuplot(adt_data, length, adt_node, adt_ik) 
write(*,*) ' iz ', iz
adt_resol = rpara(500)
!write(*,'(a,3(2x,f12.8))') ' res A^-1/px', adt_resol
!
end subroutine adt_read_data
!
!*******************************************************************************
!
subroutine adt_read_cell
!
use errlist_mod
use matrix_mod     , only:matinv
use support_mod    , only:oeffne
use trig_degree_mod, only: cosd, acosd
!
implicit none
!
integer, parameter :: IRD = 23
!
character(len=PREC_STRING) :: string
integer                    :: ios
!
if(adt_data == ' ') then
   ier_num = -6
   ier_typ = ER_COMM
   ier_msg(1) = 'ADT data file name is blank'
   return
endif
!
if(adt_cell(len_trim(adt_cell)-4:len_trim(adt_cell)) /= '.cell') then
   adt_cell(len_trim(adt_cell)+1:len_trim(adt_cell)+5) = '.cell'
endif

!
call oeffne(IRD, adt_cell, 'old')
if(ier_num/=0) then
   ier_msg(1) = 'Error opening ADT cell file'
   close(IRD)
   return
endif
!
read(IRD,'(a)', iostat=ios) string
if(string/='#ADTCELLVERSION 2') then
   ier_num = -6
   ier_typ = ER_COMM
   ier_msg(1) = 'ADT cell file is not recognized'
   ier_msg(2) = '1st line is not ''#ADTCELLVERSION 2'' '
   close(IRD)
   return
endif
read(IRD,'(a)', iostat=ios) string    !# reciprocal space cell vectors (pixels)
read(IRD,'(a)', iostat=ios) string
read(string,*) adt_omat(1,:)
read(IRD,'(a)', iostat=ios) string
read(string,*) adt_omat(2,:)
read(IRD,'(a)', iostat=ios) string
read(string,*) adt_omat(3,:)
adt_omat(:,2) = -1.0D0*adt_omat(:,2)
read(IRD,'(a)', iostat=ios) string    !# reciprocal space cell (pixels, deg)
read(IRD,'(a)', iostat=ios) string
read(string,*) adt_px_cell(1:3)       ! Length in pixels
read(IRD,'(a)', iostat=ios) string
read(string,*) adt_px_cell(4:6)       ! angles in degrees
read(IRD,'(a)', iostat=ios) string
read(IRD,'(a)', iostat=ios) string
read(string,*) adt_di_cell(1:3)       ! Direct space unit cell in Angstroem (?)
read(IRD,'(a)', iostat=ios) string
read(string,*) adt_di_cell(4:6)       ! angles in degrees
!
call matinv(adt_omat, adt_invmat)
adt_omat_t   = TRANSPOSE(adt_omat)
adt_invmat_t = TRANSPOSE(adt_invmat)
!
adt_gten(1,1) = adt_di_cell(1)**2
adt_gten(2,2) = adt_di_cell(2)**2
adt_gten(3,3) = adt_di_cell(3)**2
adt_gten(1,2) = adt_di_cell(1)*adt_di_cell(2)*cosd(adt_di_cell(6))
adt_gten(1,3) = adt_di_cell(1)*adt_di_cell(3)*cosd(adt_di_cell(5))
adt_gten(2,3) = adt_di_cell(2)*adt_di_cell(3)*cosd(adt_di_cell(4))
adt_gten(2,1) = adt_gten(1,2)
adt_gten(3,1) = adt_gten(1,3)
adt_gten(3,2) = adt_gten(2,3)
if(adt_resol/=0.0D0) then               ! Data were read, use resol from data
   adt_re_cell(1:3) = adt_px_cell(1:3)*adt_resol
   adt_re_cell(4:6) = adt_px_cell(4:6)       ! angles in degrees
   adt_rten(1,1) = adt_re_cell(1)**2
   adt_rten(2,2) = adt_re_cell(2)**2
   adt_rten(3,3) = adt_re_cell(3)**2
   adt_rten(1,2) = adt_re_cell(1)*adt_re_cell(2)*cosd(adt_re_cell(6))
   adt_rten(1,3) = adt_re_cell(1)*adt_re_cell(3)*cosd(adt_re_cell(5))
   adt_rten(2,3) = adt_re_cell(2)*adt_re_cell(3)*cosd(adt_re_cell(4))
   adt_rten(2,1) = adt_rten(1,2)
   adt_rten(3,1) = adt_rten(1,3)
   adt_rten(3,2) = adt_rten(2,3)
!  write(*,*)
!  write(*,'(a,3(2x,f12.4))') ' RECIPROCAL ', adt_re_cell(1:3)
!  write(*,*)
else                                   ! No data present, calc rec. lattice params
   call matinv(adt_gten, adt_rten)
   adt_re_cell(1) = sqrt(adt_rten(1,1))
   adt_re_cell(2) = sqrt(adt_rten(2,2))
   adt_re_cell(3) = sqrt(adt_rten(3,3))
   adt_re_cell(4) = acosd(adt_rten(2,3)/adt_re_cell(2)/adt_re_cell(3))
   adt_re_cell(5) = acosd(adt_rten(1,3)/adt_re_cell(1)/adt_re_cell(3))
   adt_re_cell(6) = acosd(adt_rten(1,2)/adt_re_cell(1)/adt_re_cell(2))
endif
write(*,*)
write(*,'(a,3(2x,f12.4))') ' Direct     ', adt_di_cell(1:3)
write(*,'(a,3(2x,f12.4))') ' Direct     ', adt_di_cell(4:6)
write(*,*)
write(*,'(a,3(2x,f12.4))') ' gten       ', adt_gten   (1,:)
write(*,'(a,3(2x,f12.4))') ' gten       ', adt_gten   (2,:)
write(*,'(a,3(2x,f12.4))') ' gten       ', adt_gten   (3,:)
write(*,*)
write(*,*)
write(*,'(a,3(2x,f12.4))') ' Reciprocal ', adt_re_cell(1:3)
write(*,'(a,3(2x,f12.4))') ' Reciprocal ', adt_re_cell(4:6)
write(*,*)
write(*,'(a,3(2x,f12.4))') ' rten       ', adt_rten   (1,:)
write(*,'(a,3(2x,f12.4))') ' rten       ', adt_rten   (2,:)
write(*,'(a,3(2x,f12.4))') ' rten       ', adt_rten   (3,:)
write(*,*)
write(*,*)
!adt_resol = (adt_re_cell(1)/adt_px_cell(1) + adt_re_cell(2)/adt_px_cell(2) + &
!             adt_re_cell(3)/adt_px_cell(3) ) / 3.0D0
!write(*,'(a,3(2x,f12.8))') ' res A^-1/px', adt_resol
!
close(IRD)
!
end subroutine adt_read_cell
!
!*******************************************************************************
!
subroutine adt_trans(zeile, length)
!
!
use ber_params_mod, only:ber_params
use errlist_mod
use get_params_mod, only:get_params
use precision_mod
use str_comp_mod  , only:str_comp
!
implicit none
!
character(len=*), intent(inout) :: zeile         ! input line
integer         , intent(inout) :: length        ! input length
!
integer, parameter :: MAXW = 6                   ! Actual parameter number
!
character(len=PREC_STRING), dimension(MAXW) :: cpara
integer                   , dimension(MAXW) :: lpara
real(kind=PREC_DP)        , dimension(MAXW) :: werte
integer :: ianz          ! Number of parameters
!
real(kind=PREC_DP)        , dimension(3) :: pixel
real(kind=PREC_DP)        , dimension(3) :: zero
real(kind=PREC_DP)        , dimension(3) :: hkl
!
call get_params(zeile, ianz, cpara, lpara, MAXW, length)
if(ier_num/=0) return
call ber_params(ianz, cpara, lpara, werte, MAXW)
!
pixel = werte(1:3)
zero  = werte(4:6)
!
hkl = matmul(adt_invmat_t,(pixel-zero))
!
write(*,'(a, 3(1x,f10.4))') ' HKL ', hkl
!
end subroutine adt_trans
!
!*******************************************************************************
!
subroutine adt_extract
!-
!  Extract the selected volume into the next data set
!+
!
use kuplot_mod
use kuplot_global
!
use errlist_mod
use lib_data_struc_type_mod
use lib_data_struc_h5 , only:data2local, local2data, dgl5_set_h5_is_ku, dgl5_set_ku_is_h5
use support_mod       , only:oeffne
!
implicit none
!
integer, parameter :: IWR = 95
character(len=PREC_STRING) :: string
integer                                    :: i,j,k      ! Loop indices
integer                                    :: ii,jj,kk   ! Loop indices
integer                                    :: ikk        ! Kuplot data set result
integer                                    :: nref       ! Reflection counter
integer                   , dimension(3)   :: idims      ! Deimensions of resulting grid
integer                   , dimension(3)   :: ipixel     ! Coordinates in ADT data
real(kind=PREC_DP)        , dimension(3)   :: pixel      ! Coordinates in ADT data
real(kind=PREC_DP)        , dimension(3)   :: hkl
real(kind=PREC_DP)        , dimension(3)   :: hkl_lim    ! HKL Limits including sigma
integer                   , dimension(2,3) :: pix_lim    ! PIXEL Limits including sigma
!
type(h5_data_struc) :: ik1
type(h5_data_struc) :: ik2
!
logical :: lout = .true.
real(kind=PREC_DP), dimension(3) :: com
real(kind=PREC_DP) :: weight
!
! Load the grid into a local data copy
!
call data2local(adt_ik  , ier_num, ier_typ, ik1%data_num   , ik1%infile,     &
     ik1%layer, ik1%is_direct, ik1%ndims, ik1%dims, ik1%is_grid,            &
     ik1%has_dxyz, ik1%has_dval, ik1%corners, ik1%vectors, ik1%cr_a0, ik1%cr_win,  &
     ik1%x, ik1%y, ik1%z, ik1%dx, ik1%dy, ik1%dz, ik1%datamap, ik1%sigma,       &
     ik1%llims, ik1%steps,  ik1%steps_full, ik1%minmaxval, ik1%minmaxcoor)
!
pix_lim(1,:) =  HUGE(1)
pix_lim(2,:) = -HUGE(1)
hkl_lim = adt_limits + adt_sigma     ! Set grand HKL limits
!
if(ALL(adt_zero<0.0D0)) then
   adt_zero = ik1%dims/2
endif
!write(*,*) ' ADT_zero at : ', adt_zero
!
! Transform the eight corners of +-hkl_lim to obtain Pixel limits
!
hkl = hkl_lim                                !(+++)
ipixel = int(matmul(adt_omat_t,hkl) + adt_zero)
pix_lim(1,:) = min(pix_lim(1,:), ipixel)
pix_lim(2,:) = max(pix_lim(2,:), ipixel)
!
hkl(1) = -hkl(1)                             !(-++)
ipixel = int(matmul(adt_omat_t,hkl) + adt_zero)
pix_lim(1,:) = min(pix_lim(1,:), ipixel)
pix_lim(2,:) = max(pix_lim(2,:), ipixel)
!
hkl(2) = -hkl(2)                             !(--+)
ipixel = int(matmul(adt_omat_t,hkl) + adt_zero)
pix_lim(1,:) = min(pix_lim(1,:), ipixel)
pix_lim(2,:) = max(pix_lim(2,:), ipixel)
!
hkl(1) = -hkl(1)                             !(+-+)
ipixel = int(matmul(adt_omat_t,hkl) + adt_zero)
pix_lim(1,:) = min(pix_lim(1,:), ipixel)
pix_lim(2,:) = max(pix_lim(2,:), ipixel)
!
hkl(3) = -hkl(3)                             !(+--)
ipixel = int(matmul(adt_omat_t,hkl) + adt_zero)
pix_lim(1,:) = min(pix_lim(1,:), ipixel)
pix_lim(2,:) = max(pix_lim(2,:), ipixel)
!
hkl(1) = -hkl(1)                             !(---)
ipixel = int(matmul(adt_omat_t,hkl) + adt_zero)
pix_lim(1,:) = min(pix_lim(1,:), ipixel)
pix_lim(2,:) = max(pix_lim(2,:), ipixel)
!
hkl(2) = -hkl(2)                             !(-+-)
ipixel = int(matmul(adt_omat_t,hkl) + adt_zero)
pix_lim(1,:) = min(pix_lim(1,:), ipixel)
pix_lim(2,:) = max(pix_lim(2,:), ipixel)
!
hkl(1) = -hkl(2)                             !(++-)
ipixel = int(matmul(adt_omat_t,hkl) + adt_zero)
pix_lim(1,:) = min(pix_lim(1,:), ipixel)
pix_lim(2,:) = max(pix_lim(2,:), ipixel)
!
pix_lim(1,:) = max(          1, pix_lim(1,:)  )
pix_lim(2,:) = min(ik1%dims(:), pix_lim(2,:)+1)
!
!write(*,'(a,3(2x,i6   ))') ' MIN PIXEL', pix_lim(1,:)
!write(*,'(a,3(2x,i6   ))') ' MAX PIXEL', pix_lim(2,:)
!
idims = nint(2*adt_limits/adt_steps)+1
!write(*,*) ' DIMENSIONS ', idims
!
if(adt_mode==ADT_RODS) then       ! Extract all rods
!
!write(*,*) ' DIMENSION ', ik1%dims, ' RODS ', idims 
!write(*,*) ' PIXEL LIM ', pix_lim(:,1)
!write(*,*) ' PIXEL LIM ', pix_lim(:,2)
!write(*,*) ' PIXEL LIM ', pix_lim(:,3)
!read(*,*) i
   allocate(ik2%datamap(idims(1), idims(2), idims(3)))
   ik2%datamap = 0.0D0
!write(*,*) ' ALLOCATED '
!read(*,*) i
   do k=pix_lim(1,3), pix_lim(2,3)
      pixel(3) = real(k,kind=PREC_DP)
      do j=pix_lim(1,2), pix_lim(2,2)
         pixel(2) = real(j,kind=PREC_DP)
         do i=pix_lim(1,1), pix_lim(2,1)
            pixel(1) = real(i,kind=PREC_DP)
            hkl = matmul(adt_invmat_t, (pixel-adt_zero))
            if(abs(hkl(1)-nint(hkl(1)))<=adt_sigma(1)) then
               if(real(nint(abs(hkl(1))),kind=PREC_DP)<=adt_limits(1)) then
               if(abs(hkl(2)-nint(hkl(2)))<=adt_sigma(2)) then
                  if(real(nint(abs(hkl(2))),kind=PREC_DP)<=adt_limits(2)) then
                  if(real(    (abs(hkl(3))),kind=PREC_DP)<=adt_limits(3)) then
                  ii = nint((hkl(1)+adt_limits(1))/adt_steps(1))+1
                  jj = nint((hkl(2)+adt_limits(2))/adt_steps(2))+1
                  kk = nint((hkl(3)+adt_limits(3))/adt_steps(3))+1
!write(*,'(a, 3i4,3f7.2,3i4,f10.2)') ' PIXEL, hkl ii', i,j,k, hkl, ii,jj,kk, ik1%datamap(i,j,k)
                  ik2%datamap(ii,jj,kk) = ik2%datamap(ii,jj,kk) + ik1%datamap(i,j,k)
               endif
               endif
               endif
               endif
            endif
         enddo
      enddo
   enddo
!
!write(*,*) ' DONE WITH LOOP '
!
   ikk = iz
   ik2%data_num = 0
   ik2%infile = 'extracted_rods'
   ik2%layer = idims(3)/2
   ik2%is_direct = .false.
   ik2%ndims     = 3
   ik2%dims      = idims
   ik2%is_grid   = .true.
   ik2%has_dxyz  = .false.
   ik2%has_dval  = .false.
   ik2%corners(:,1) = -adt_limits         ! Lower Left
   ik2%corners(:,2) = -adt_limits         ! Lower right
   ik2%corners(1,2) = -ik2%corners(1,2)
   ik2%corners(:,3) = -adt_limits         ! Upper left
   ik2%corners(2,3) = -ik2%corners(2,2)
   ik2%corners(:,4) = -adt_limits         ! Top left
   ik2%corners(3,4) = -ik2%corners(3,4)
   ik2%vectors(:,1) = ik2%corners(:,2) - ik2%corners(:,1)
   ik2%vectors(:,2) = ik2%corners(:,3) - ik2%corners(:,1)
   ik2%vectors(:,3) = ik2%corners(:,4) - ik2%corners(:,1)
   ik2%cr_a0          = adt_di_cell(1:3)
   ik2%cr_win         = adt_di_cell(4:6)
   ik2%llims          = -adt_limits
   ik2%steps(1)       = ik2%vectors(1,1)/(real(ik2%dims(1),kind=PREC_DP)-1.0D0)
   ik2%steps(2)       = ik2%vectors(2,2)/(real(ik2%dims(2),kind=PREC_DP)-1.0D0)
   ik2%steps(3)       = ik2%vectors(3,3)/(real(ik2%dims(3),kind=PREC_DP)-1.0D0)
   ik2%steps_full(:,1)= ik2%vectors(:,1)/(real(ik2%dims(1),kind=PREC_DP)-1.0D0)
   ik2%steps_full(:,2)= ik2%vectors(:,2)/(real(ik2%dims(2),kind=PREC_DP)-1.0D0)
   ik2%steps_full(:,3)= ik2%vectors(:,3)/(real(ik2%dims(3),kind=PREC_DP)-1.0D0)
!write(*,*) ' ABOUT TO SAVE ' 
!write(*,*) ' ik2 data_num  ', ik2%data_num
!write(*,*) ' ik2 nlayer    ', ik2%layer 
!write(*,*) ' ik2 ndim      ', ik2%ndims
!write(*,*) ' ik2  dim      ', ik2%dims
!write(*,*) ' ik2  x        ', allocated(ik2%x)
!write(*,*) ' ik2  ll       ', ik2%corners(:,1)
!write(*,*) ' ik2  lr       ', ik2%corners(:,2)
!write(*,*) ' ik2  ul       ', ik2%corners(:,3)
!write(*,*) ' ik2  tl       ', ik2%corners(:,4)
!write(*,*) ' ik2  vec x    ', ik2%vectors(:,1)
!write(*,*) ' ik2  vec y    ', ik2%vectors(:,2)
!write(*,*) ' ik2  vec z    ', ik2%vectors(:,3)
!write(*,*) ' ik2  steps    ', ik2%steps
!write(*,*) ' ik2  full     ', ik2%steps_full(1,1), ik2%steps_full(2,2) , ik2%steps_full(3,3)
!write(*,*) ' ik2  data     ', allocated(ik2%datamap)
!write(*,*) ' ik2  sigma    ', allocated(ik2%sigma)
!write(*,*) ' ik2  data     ', lbound(ik2%datamap), ubound(ik2%datamap)
!write(*,*) ' setp into local2data '
!read(*,*) i
write(*,*) ' ikk ' , ikk, iz
   call local2data(ikk, ier_num, ier_typ, ik2%data_num, ik2%infile, ik2%layer,  &
        ik2%is_direct, ik2%ndims, ik2%dims, ik2%is_grid, ik2%has_dxyz,             &
        ik2%has_dval, ik2%corners, ik2%vectors, ik2%cr_a0, ik2%cr_win, ik2%x, ik2%y,     &
        ik2%z, ik2%dx, ik2%dy, ik2%dz, ik2%datamap, ik2%sigma, ik2%llims, ik2%steps,  &
        ik2%steps_full)
write(*,*) ' ikk ' , ikk, iz, ik2%data_num
   call dgl5_set_h5_is_ku(iz, ik2%data_num)
   call dgl5_set_ku_is_h5(ik2%data_num, iz)
   ku_ndims(iz) = ik2%ndims
   string = 'adt_rods.inte'
   call data2kuplot(ikk, string, lout)
write(*,*) ' ikk ' , ikk, iz, ik2%data_num
!
!write(*,*) ' FINISHED extraction'
do j=1,idims(2)
do i=1,idims(1)
write(string,'(a,i3.3,a,i3.3,a)') 'RODS/rod_',abs(-4+i),'_',abs(-4+j),'_L.inte'
if((-4+i)<0) string(5:5)='-'
if((-4+j)<0) string(9:9)='-'

open(96,file=string(1:len_trim(string)), status='unknown')
do kk=1, idims(3)
  write(96,'(2f20.4)') -adt_limits(3)+(kk-1)*adt_steps(3), ik2%datamap(i,j,kk)
enddo
close(96)
enddo
enddo
!
elseif(adt_mode==ADT_BRAGG) then
   nref = 0
   if(adt_log) then
      call oeffne(IWR, adt_logfile, 'unknown')
      write(IWR,'(a, 43x,a,17x,a,17x,a,13x,a)') '#', 'Observed', 'Calculated', 'Calculated', 'Nominal'
      write(IWR,'(a6,1x,a,1x,a,8x,a,12x,a,21x,a,23x,a,16x,a)') '#   Nr',        &
      'Target', 'dummy', 'Weight', 'Pixel', 'Pixel', 'hkl', 'H   K   L'
   endif
   do k=-nint(adt_limits(3)), nint(adt_limits(3))
   do j=-nint(adt_limits(2)), nint(adt_limits(2))
   do i=-nint(adt_limits(1)), nint(adt_limits(1))
if(mod(i+j,2)==0) then
      hkl(3) = real(k, kind=PREC_DP)
      hkl(2) = real(j, kind=PREC_DP)
      hkl(1) = real(i, kind=PREC_DP)
      pixel = matmul(adt_omat_t, hkl) + adt_zero
      if(pixel(1)>5.0 .and. pixel(1)<ik1%dims(1)-5.0  .and. &
         pixel(2)>5.0 .and. pixel(2)<ik1%dims(2)-5.0  .and. &
         pixel(3)>5.0 .and. pixel(3)<ik1%dims(3)-5.0  ) then  
      com = 0.0D0
      weight = 0.0D0
      do kk = nint(pixel(3)-5), nint(pixel(3)+5)
      do jj = nint(pixel(2)-5), nint(pixel(2)+5)
      do ii = nint(pixel(1)-5), nint(pixel(1)+5)
         com(1) = com(1) + ii*ik1%datamap(ii,jj,kk)
         com(2) = com(2) + jj*ik1%datamap(ii,jj,kk)
         com(3) = com(3) + kk*ik1%datamap(ii,jj,kk)
         weight = weight + ik1%datamap(ii,jj,kk)
      enddo
      enddo
      enddo
      if(weight>1) then
         com = com/weight
         hkl = matmul(adt_invmat_t, (com-adt_zero))
         nref = nref + 1
         weight = weight/11.**3
         if(adt_log) then
         write(95,'(f6.0, 2(f5.1,2x),f15.6,3(3f8.2,2x),2x,3i4)') real(nref), 1.0, 0.0,  &
         1./weight, com, pixel, hkl, i,j,k
         endif
      endif
      endif
endif
   enddo
   enddo
   enddo
   close(IWR)
endif
!
call adt_clean(ik1)
call adt_clean(ik2)
!
end subroutine adt_extract
!
!*******************************************************************************
!
subroutine adt_clean(ik1)
!-
!  deallocate all local arrays
!+
use lib_data_struc_type_mod
!
type(h5_data_struc), intent(inout) :: ik1
!
if(allocated(ik1%x    )) deallocate(ik1%x    )
if(allocated(ik1%y    )) deallocate(ik1%y    )
if(allocated(ik1%z    )) deallocate(ik1%z    )
if(allocated(ik1%dx   )) deallocate(ik1%dx   )
if(allocated(ik1%dy   )) deallocate(ik1%dy   )
if(allocated(ik1%dz   )) deallocate(ik1%dz   )
if(allocated(ik1%datamap )) deallocate(ik1%datamap )
if(allocated(ik1%sigma)) deallocate(ik1%sigma)
!
!
end subroutine adt_clean
!
!*******************************************************************************
!
subroutine adt_reset
!
implicit none
!
adt_cell = ' '
adt_data = ' '
adt_mode = 0
!
end subroutine adt_reset
!
!*******************************************************************************
!
end module kuplot_adt_mod
