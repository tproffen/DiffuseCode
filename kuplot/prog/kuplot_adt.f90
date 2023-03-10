module kuplot_adt_mod
!-
! Routines to read the MAINZ ADT data
!+
use precision_mod
!
private
public adt_menu
!
integer, parameter         :: ADT_STAB   =-1
integer, parameter         :: ADT_RODS   = 1
integer, parameter         :: ADT_PLANES = 2
integer, parameter         :: ADT_SPACE  = 3
integer, parameter         :: ADT_BRAGG  = 4
integer, parameter         :: ADT_VOLUME = 5
integer, parameter         :: ADT_PRIMI  = 1
integer, parameter         :: ADT_CENT_A = 2
integer, parameter         :: ADT_CENT_B = 3
integer, parameter         :: ADT_CENT_C = 4
integer, parameter         :: ADT_CENT_F = 5
integer, parameter         :: ADT_CENT_I = 6
integer, parameter         :: ADT_CENT_R = 7
character(len=PREC_STRING) :: adt_data        ! Input MRC file
character(len=PREC_STRING) :: adt_cell        ! Input MRC cell information
character(len=PREC_STRING) :: adt_logfile     ! Output file for extraction
logical                            :: adt_init = .true.  ! ADT needs to be initialized
logical                            :: adt_log  = .false. ! Extract log on / off
integer                            :: adt_mode = 0! Extraction mode rods, plane, bulk
integer                            :: adt_node = 0! Node in global data structure
integer                            :: adt_ik   = 0! ADT Data set in kuplot
integer                            :: adt_lattice = ADT_PRIMI ! Lattice type
real(kind=PREC_DP), dimension(3,3) :: adt_omat    ! Orientation matrix reciprocal space in pixels
real(kind=PREC_DP), dimension(3,3) :: adt_omat_t  ! Orientation matrix reciprocal space in pixels TRANSPOSED
real(kind=PREC_DP), dimension(3,3) :: adt_invmat  ! Inverse orientation matrix
real(kind=PREC_DP), dimension(3,3) :: adt_invmat_t! Inverse orientation matrix TRANSPOSED
real(kind=PREC_DP), dimension(3,3) :: adt_gten    ! Metric tensor
real(kind=PREC_DP), dimension(3,3) :: adt_rten    ! Reciprocal Metric tensor
real(kind=PREC_DP), dimension(3,3) :: adt_transm  ! Matrix: hkl onto abs, od, interval
real(kind=PREC_DP), dimension(3,3) :: adt_transi  ! Matrix: abs, od, interval onto hkl
real(kind=PREC_DP), dimension(6)   :: adt_px_cell ! Reciprocel cell in Pixels, degrees
real(kind=PREC_DP), dimension(6)   :: adt_di_cell ! Direct     cell in Angstr, degrees
real(kind=PREC_DP), dimension(6)   :: adt_re_cell ! Reciprocal cell in Angstr, degrees
real(kind=PREC_DP), dimension(3)   :: adt_limits  ! Upper limits in reciprocal space
real(kind=PREC_DP), dimension(3)   :: adt_rod     ! Rod directions to extract
real(kind=PREC_DP), dimension(3)   :: adt_plane   ! Plane normal   to extract
real(kind=PREC_DP), dimension(3)   :: adt_absc    ! Abscissa for extracted data
real(kind=PREC_DP), dimension(3)   :: adt_ordi    ! Ordinate for extracted data
real(kind=PREC_DP), dimension(3)   :: adt_inter   ! Interval for extracted data
real(kind=PREC_DP), dimension(0:3) :: adt_labsc   ! LIMITS Abscissa for extracted data
real(kind=PREC_DP), dimension(0:3) :: adt_lordi   ! LIMITS Ordinate for extracted data
real(kind=PREC_DP), dimension(0:3) :: adt_linter  ! LIMITS Interval for extracted data
real(kind=PREC_DP), dimension(3)   :: adt_sigma   ! Sigma in reciprocal units for extraction
real(kind=PREC_DP), dimension(3)   :: adt_steps   ! Steps in reciprocal units for extraction
real(kind=PREC_DP), dimension(3)   :: adt_zero    = -1.0D0 !Pixel location of zero point
real(kind=PREC_DP)                 :: adt_resol   =  0.0D0 ! Resolution A^-1 pro pixel
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
character(len=9          ) :: befehl       ! The actual command
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
if(adt_init) then
   call adt_reset
   adt_init = .false.
endif
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
                 .AND..NOT. (str_comp (befehl, 'system', 2, lbef, 6) )    &
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
   if_rese: IF(str_comp (befehl, 'reset', 2, lbef, 5) ) THEN
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
   elseif(str_comp (befehl, 'transform', 2, lbef, 9) ) then
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
elseif(adt_mode==ADT_PLANES) then
  write(output_io, *) 
  write(output_io, '(a)') ' Extracting planes in reciprocal space'
  write(output_io, '(a, 2x,3(2x,f10.4))') ' Planes normal :', adt_plane
endif
  write(output_io, '(a, 2x,3(2x,f10.4))') ' Limits        :', adt_limits
  write(output_io, '(a, 2x,3(2x,f10.4))') ' steps         :', adt_steps
  write(output_io, '(a, 2x,3(2x,f10.4))') ' sigmas        :', adt_sigma
  write(output_io, '(a, 2x,3(2x,f10.4))') ' Zero at       :', adt_zero
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
logical :: lout
!
integer, parameter :: NOPTIONAL = 3
integer, parameter :: O_DATA    = 1
integer, parameter :: O_CELL    = 2
integer, parameter :: O_LOUT    = 3
character(LEN=   4), dimension(NOPTIONAL) :: oname   !Optional parameter names
character(LEN=PREC_STRING), dimension(NOPTIONAL) :: opara   !Optional parameter strings returned
integer            , dimension(NOPTIONAL) :: loname  !Lenght opt. para name
integer            , dimension(NOPTIONAL) :: lopara  !Lenght opt. para name returned
logical            , dimension(NOPTIONAL) :: lpresent!opt. para is present
real(kind=PREC_DP) , dimension(NOPTIONAL) :: owerte   ! Calculated values
integer, parameter                        :: ncalc = 0 ! Number of values to calculate 
!
data oname  / 'data', 'cell' ,'log ' /
data loname /  4    ,  4     , 3     /
opara  =  (/ 'unknown','unknown', 'off    ' /)   ! Always provide fresh default values
lopara =  (/  7,        7       ,  3        /)
owerte =  (/  0.0,      0.0     ,  0.0      /)
!
!
!
call get_params(zeile, ianz, cpara, lpara, maxw, length)
if(ier_num/=0) return
call get_optional(ianz, MAXW, cpara, lpara, NOPTIONAL,  ncalc, &
                  oname, loname, opara, lopara, lpresent, owerte)
if(ier_num/=0) return
!
lout = opara(O_LOUT)=='screen'
!
if(lpresent(O_DATA)) then
   adt_data = opara(O_DATA)(1:lopara(O_DATA))
   if(adt_data /= 'unknown') call adt_read_data
endif
!
if(lpresent(O_CELL)) then
   adt_cell = opara(O_CELL)(1:lopara(O_CELL))
   if(adt_cell /= 'unknown') call adt_read_cell(lout)
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
use take_param_mod
!
implicit none
!
character(len=*), intent(inout) :: zeile         ! input line
integer         , intent(inout) :: length        ! input length
!
!integer, parameter :: MIN_PARA = 3               ! Minimum parameter numbers
integer, parameter :: MAXW = 7                   ! Actual parameter number
!
character(len=PREC_STRING), dimension(MAXW) :: cpara
integer                   , dimension(MAXW) :: lpara
real(kind=PREC_DP)        , dimension(MAXW) :: werte
integer :: ianz          ! Number of parameters
logical :: lok
!
integer, parameter :: NOPTIONAL = 12
integer, parameter :: O_MODE    = 1
integer, parameter :: O_LIMITS  = 2
integer, parameter :: O_STEPS   = 3
integer, parameter :: O_SIGMA   = 4
integer, parameter :: O_ZERO    = 5
integer, parameter :: O_RODS    = 6
integer, parameter :: O_NORMAL  = 7
integer, parameter :: O_ABS     =  8
integer, parameter :: O_ORD     =  9
integer, parameter :: O_INTER   = 10
integer, parameter :: O_LOGFILE = 11
integer, parameter :: O_LATTICE = 12
character(LEN=   8), dimension(NOPTIONAL) :: oname   !Optional parameter names
character(LEN=PREC_STRING), dimension(NOPTIONAL) :: opara   !Optional parameter strings returned
integer            , dimension(NOPTIONAL) :: loname  !Lenght opt. para name
integer            , dimension(NOPTIONAL) :: lopara  !Lenght opt. para name returned
logical            , dimension(NOPTIONAL) :: lpresent!opt. para is present
real(kind=PREC_DP) , dimension(NOPTIONAL) :: owerte   ! Calculated values
integer, parameter                        :: ncalc = 0 ! Number of values to calculate 
!
data oname  / 'mode'        , 'limits'      , 'steps'       , 'sigma', &
              'zero'        , 'rod'         , 'normal'      , 'abs'  , &
              'ord'         , 'interval'    , 'output'      , 'lattice'      /
data loname /  4            ,  6            ,  4            ,  5     ,  &
               4            ,  3            ,  6            ,  3     ,  &
               3            ,  8            ,  6            ,  7            /
opara  =  (/ 'rods            ','[1.0,1.0,1.0]   ','[1.0,1.0,0.1]   ', '[0.2,0.2,0.0]   ', &
             '[-1.,-1.,-1.]   ','[0.0,0.0,1.0]   ','[0.0,0.0,1.0]   ', '[1.0,0.0,0.0, 5]', &
             '[0.0,1.0,0.0, 5]','[0.0,0.0,1.0, 5]','nolog.data      ', 'P               '  /) !Provide fresh default values
lopara =  (/  4             , 13            ,  13            ,  13            , &
              13            , 13            ,  13            ,  16            , &
              16            , 16            ,  13            , 1               /)
owerte =  (/  0.0           , 0.0           , 0.0            ,  0.0           , &
              0.0           , 0.0           , 0.0            ,  0.0           , &
              0.0           , 0.0           , 0.0            ,  1.0             /)
!
lok = .false.
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
   elseif(opara(O_MODE)=='stab') then
      adt_mode = ADT_STAB
   elseif(opara(O_MODE)=='plane') then
      adt_mode = ADT_PLANES
   elseif(opara(O_MODE)=='bragg') then
      adt_mode = ADT_BRAGG
   elseif(opara(O_MODE)=='volume') then
      adt_mode = ADT_VOLUME
   else
      ier_num = -6
      ier_typ = ER_COMM
      return
   endif
endif
!
if(lpresent(O_LOGFILE)) then
   adt_logfile = opara(O_LOGFILE)
   adt_log = .true.
   lok = .true.
else
   adt_log = adt_logfile /= 'nolog.data'
   lok = .true.
endif
!
if(lpresent(O_LIMITS)) then
   call get_optional_multi(MAXW, opara(O_LIMITS), lopara(O_LIMITS), werte, ianz)
   adt_limits = werte(1:3)
   lok = .true.
endif
!
if(lpresent(O_SIGMA )) then
   call get_optional_multi(MAXW, opara(O_SIGMA ), lopara(O_SIGMA ), werte, ianz)
   if(ianz==3) then
      adt_sigma  = werte(1:3)
   elseif(ianz==2) then
      adt_sigma(1:2) = werte(1:2)
      adt_sigma(3)   = 0.0D0
   elseif(ianz==1) then
      adt_sigma = werte(1)
   endif
   lok = .true.
endif
!
if(lpresent(O_STEPS )) then
   call get_optional_multi(MAXW, opara(O_STEPS ), lopara(O_STEPS ), werte, ianz)
   if(ianz==3) then
      adt_steps  = werte(1:3)
   elseif(ianz==2) then
      adt_steps(1:2) = werte(1:2)
      adt_steps(3)   = 0.0D0
   elseif(ianz==1) then
      adt_steps = werte(1)
   endif
   lok = .true.
endif
!
if(lpresent(O_ZERO )) then
   call get_optional_multi(MAXW, opara(O_ZERO ), lopara(O_ZERO ), werte, ianz)
   adt_zero   = werte(1:3)
   lok = .true.
endif
!
if(lpresent(O_RODS )) then
   call get_optional_multi(MAXW, opara(O_RODS ), lopara(O_RODS ), werte, ianz)
   adt_rod    = werte(1:3)
   if(ianz==4) then              ! DELTA[hkl]; maxsteps
      adt_inter  = werte(1:3)
      adt_linter = 1
      adt_linter    = (werte(4))
      adt_linter(0) = 1.0
   elseif(ianz==6) then          ! DELTA[hkl]; max[hkl]
      adt_inter  = werte(1:3)
      adt_linter = 3
      adt_linter    = (werte(4:7))
      adt_linter(0) = 1.0
   elseif(ianz==7) then          ! [hkl]; DELTA[hkl]; maxsteps
      adt_inter  = werte(4:6)
      adt_linter = 2
      adt_linter    = (werte(7))
      adt_linter(0) = 1.0
   endif
   lok = .true.
endif
!
if(lpresent(O_NORMAL)) then
   call get_optional_multi(MAXW, opara(O_NORMAL), lopara(O_NORMAL), werte, ianz)
   adt_plane  = werte(1:3)
   if(ianz==5) then
      adt_inter  = werte(4)
      adt_linter = 0
      adt_linter    = (werte(5))
      adt_linter(0) = 1.0
   elseif(ianz==4) then
      adt_inter = werte(4)
      adt_linter = 0
   endif
   lok = .true.
endif
!
if(lpresent(O_ABS )) then
   call get_optional_multi(MAXW, opara(O_ABS ), lopara(O_ABS ), werte, ianz)
   adt_absc   = werte(1:3)
   if(ianz==4) then
      adt_labsc    = 0
      adt_labsc(1) = (werte(4))
      adt_labsc(0) = 1
   elseif(ianz==6) then
      adt_labsc(1:3) = werte(4:6)
      adt_labsc(0) = 3
   endif
   lok = .true.
endif
!
if(lpresent(O_ORD )) then
   call get_optional_multi(MAXW, opara(O_ORD ), lopara(O_ORD ), werte, ianz)
   adt_ordi   = werte(1:3)
   if(ianz==4) then
      adt_lordi    = 0
      adt_lordi(1) = (werte(4))
      adt_lordi(0) = 1
   elseif(ianz==6) then
      adt_lordi(1:3) = werte(4:6)
      adt_lordi(0) = 3
   endif
   lok = .true.
endif
!
if(lpresent(O_INTER )) then
   call get_optional_multi(MAXW, opara(O_INTER ), lopara(O_INTER ), werte, ianz)
   adt_inter  = werte(1:3)
   if(ianz==4) then
      adt_linter(1) = nint(werte(4))
      adt_linter(0) = 1
   elseif(ianz==6) then
      adt_linter(1:3) = werte(4:6)
      adt_linter(0) = 3
   endif
   lok = .true.
endif
!
if(lpresent(O_LATTICE)) then
   if(opara(O_LATTICE)=='P') then
      adt_lattice = ADT_PRIMI
   elseif(opara(O_LATTICE)=='A') then
      adt_lattice = ADT_CENT_A
   elseif(opara(O_LATTICE)=='B') then
      adt_lattice = ADT_CENT_B
   elseif(opara(O_LATTICE)=='C') then
      adt_lattice = ADT_CENT_C
   elseif(opara(O_LATTICE)=='F') then
      adt_lattice = ADT_CENT_F
   elseif(opara(O_LATTICE)=='I') then
      adt_lattice = ADT_CENT_I
   elseif(opara(O_LATTICE)=='R') then
      adt_lattice = ADT_CENT_R
   else
      ier_num = -6
      ier_typ = ER_COMM
   endif
   lok = .true.
endif
!
if(.not.lok) then
   ier_num = -6
   ier_typ = ER_COMM
endif
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
call mrc_read_kuplot(adt_data, length, adt_node, adt_ik) 
adt_resol = rpara(500)
!write(*,'(a,3(2x,f12.8))') ' res A^-1/px', adt_resol
!
end subroutine adt_read_data
!
!*******************************************************************************
!
subroutine adt_read_cell(lout)
!
use errlist_mod
use matrix_mod     , only:matinv
use prompt_mod
use support_mod    , only:oeffne
use trig_degree_mod, only: cosd, acosd
!
implicit none
!
logical, intent(in) :: lout
integer, parameter :: IRD = 23
!
character(len=PREC_STRING) :: string
integer                    :: ios
!
if(adt_cell == ' ') then
   ier_num = -6
   ier_typ = ER_COMM
   ier_msg(1) = 'ADT cell file name is blank'
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
read(IRD,'(a)', iostat=ios) string    !# direct space cell
read(IRD,'(a)', iostat=ios) string
read(string,*) adt_di_cell(1:3)       ! Direct space unit cell in Angstroem (?)
read(IRD,'(a)', iostat=ios) string
read(string,*) adt_di_cell(4:6)       ! angles in degrees
!
call matinv(adt_omat, adt_invmat)
adt_omat_t   = TRANSPOSE(adt_omat)
adt_invmat_t = TRANSPOSE(adt_invmat)
!
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
   call matinv(adt_rten, adt_gten)
   adt_di_cell(1) = sqrt(adt_gten(1,1))
   adt_di_cell(2) = sqrt(adt_gten(2,2))
   adt_di_cell(3) = sqrt(adt_gten(3,3))
   adt_di_cell(4) = acosd(adt_gten(2,3)/adt_di_cell(2)/adt_di_cell(3))
   adt_di_cell(5) = acosd(adt_gten(1,3)/adt_di_cell(1)/adt_di_cell(3))
   adt_di_cell(6) = acosd(adt_gten(1,2)/adt_di_cell(1)/adt_di_cell(2))
!  write(*,*)
!  write(*,'(a,3(2x,f12.4))') ' RECIPROCAL ', adt_re_cell(1:3)
!  write(*,*)
else                                   ! No data present, calc rec. lattice params
   adt_gten(1,1) = adt_di_cell(1)**2
   adt_gten(2,2) = adt_di_cell(2)**2
   adt_gten(3,3) = adt_di_cell(3)**2
   adt_gten(1,2) = adt_di_cell(1)*adt_di_cell(2)*cosd(adt_di_cell(6))
   adt_gten(1,3) = adt_di_cell(1)*adt_di_cell(3)*cosd(adt_di_cell(5))
   adt_gten(2,3) = adt_di_cell(2)*adt_di_cell(3)*cosd(adt_di_cell(4))
   adt_gten(2,1) = adt_gten(1,2)
   adt_gten(3,1) = adt_gten(1,3)
   adt_gten(3,2) = adt_gten(2,3)
   call matinv(adt_gten, adt_rten)
   adt_re_cell(1) = sqrt(adt_rten(1,1))
   adt_re_cell(2) = sqrt(adt_rten(2,2))
   adt_re_cell(3) = sqrt(adt_rten(3,3))
   adt_re_cell(4) = acosd(adt_rten(2,3)/adt_re_cell(2)/adt_re_cell(3))
   adt_re_cell(5) = acosd(adt_rten(1,3)/adt_re_cell(1)/adt_re_cell(3))
   adt_re_cell(6) = acosd(adt_rten(1,2)/adt_re_cell(1)/adt_re_cell(2))
endif
if(lout) then
   write(output_io,*)
   write(output_io,'(a,3(2x,f12.4))') ' Direct     ', adt_di_cell(1:3)
   write(output_io,'(a,3(2x,f12.4))') ' Direct     ', adt_di_cell(4:6)
   write(output_io,*)
   write(output_io,'(a,3(2x,f12.4))') ' gten       ', adt_gten   (1,:)
   write(output_io,'(a,3(2x,f12.4))') ' gten       ', adt_gten   (2,:)
   write(output_io,'(a,3(2x,f12.4))') ' gten       ', adt_gten   (3,:)
   write(output_io,*)
   write(output_io,*)
   write(output_io,'(a,3(2x,f12.4))') ' Reciprocal ', adt_re_cell(1:3)
   write(output_io,'(a,3(2x,f12.4))') ' Reciprocal ', adt_re_cell(4:6)
   write(output_io,*)
   write(output_io,'(a,3(2x,f12.4))') ' rten       ', adt_rten   (1,:)
   write(output_io,'(a,3(2x,f12.4))') ' rten       ', adt_rten   (2,:)
   write(output_io,'(a,3(2x,f12.4))') ' rten       ', adt_rten   (3,:)
   write(output_io,*)
   write(output_io,*)
endif
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
use prompt_mod
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
write(output_io,'(a, 3(1x,f10.4))') ' HKL ', hkl
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
implicit none
!
if(adt_mode==ADT_RODS) then       ! Extract all rods
  call adt_extract_rods
elseif(adt_mode==ADT_PLANES) then
  call adt_extract_planes
elseif(adt_mode==ADT_VOLUME) then
  call adt_extract_volume
elseif(adt_mode==ADT_BRAGG) then
  call adt_extract_bragg
endif
!
end subroutine adt_extract
!
!*******************************************************************************
!
subroutine adt_extract_rods
!-
!  Extract the selected rods   into the next data set
!+
!
use kuplot_mod
use kuplot_global
!
use errlist_mod
use lib_data_struc_type_mod
use lib_data_struc_h5 , only:data2local, local2data, dgl5_set_h5_is_ku, dgl5_set_ku_is_h5
use matrix_mod        , only:matinv
use support_mod       , only:oeffne
!
implicit none
!
integer, parameter :: IWR = 95
character(len=PREC_STRING) :: string
integer                                    :: i,j,k      ! Loop indices
integer                                    :: ii,jj,kk   ! Loop indices
integer                                    :: ni, nj     ! Loop limits
integer                                    :: ikk        ! Kuplot data set result
integer                                    :: nref       ! Reflection counter
integer                   , dimension(3)   :: idims      ! Deimensions of resulting grid
integer                   , dimension(3)   :: ipixel     ! Coordinates in ADT data
real(kind=PREC_DP)        , dimension(3)   :: pixel      ! Coordinates in ADT data
real(kind=PREC_DP)        , dimension(3)   :: ppixel      ! Coordinates in ADT data
real(kind=PREC_DP)        , dimension(8)   :: pweight    ! Weight for the eight corners
real(kind=PREC_DP)        , dimension(3)   :: hkl
real(kind=PREC_DP)        , dimension(3)   :: hkl_lim    ! HKL Limits including sigma
real(kind=PREC_DP)                         :: lsteps     ! Max Steps along the rods
integer                   , dimension(2,3) :: pix_lim    ! PIXEL Limits including sigma
real(kind=PREC_DP), dimension(:,:,:), allocatable :: cweight
!
type bragg
   real(kind=PREC_DP)               :: weight
   real(kind=PREC_DP), dimension(3) :: com
   real(kind=PREC_DP), dimension(3) :: pixel
   real(kind=PREC_DP), dimension(3) :: hkl
   integer           , dimension(3) :: ihkl
end type bragg
type(bragg), dimension(:), allocatable :: bragg_list
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
if(ALL(adt_zero<0.0D0)) then
   adt_zero = ik1%dims/2
endif
!write(*,*) ' ADT_zero at : ', adt_zero
!
!  Copy temporary maximum steps into local storage
!
lsteps = adt_linter(1)
!
! Build matrix FinalPixel => hkl
!
adt_transi(:,1) = adt_absc
adt_transi(:,2) = adt_ordi
adt_transi(:,3) = adt_inter
call matinv(adt_transi, adt_transm)  ! Build matrix hkl => FinalPixel
!
! Determine maximum final dimensions
!
if(nint(adt_labsc(0))==1) then
   idims(1) = nint(2*adt_labsc(1)) + 1
else
   hkl = adt_labsc(1:3)
   ppixel = matmul(adt_transm, hkl)
   idims(1) = nint(2*ppixel(1)) + 1
endif
if(nint(adt_lordi(0))==1) then
   idims(2) = nint(2*adt_lordi(1)) + 1
else
   hkl = adt_lordi(1:3)
   ppixel = matmul(adt_transm, hkl)
   idims(2) = nint(2*ppixel(2)) + 1
endif
if(nint(adt_linter(0))==1) then       ! DELTA[hkl]; maxsteps
   idims(3) = nint(2*(adt_linter(1))) + 1
   adt_linter(1:3) = adt_inter*lsteps
elseif(nint(adt_linter(0))==2) then   ! [hkl], DELTA[hkl]; maxsteps
   idims(3) = nint(2*(adt_linter(1))) + 1
   adt_linter(1:3) = adt_inter*lsteps
elseif(nint(adt_linter(0))==3) then   ! [hkl], max  [hkl]
   hkl = adt_linter(1:3) 
   ppixel = matmul(adt_transm, hkl)
   idims(3) = nint(2*ppixel(3)) + 1
endif
!
!  Set hkl search limits to determine smallest volume in pixel space 
!
adt_limits = adt_labsc(1:3) + adt_lordi(1:3) + adt_linter(1:3)
!
!write(*,*)
!write(*,*) ' ADT_LABSC  ',adt_labsc
!write(*,*) ' ADT_LORDI  ',adt_lordi
!write(*,*) ' ADT_LINTER ',adt_linter
!write(*,*) ' ADT_LIMITS ',adt_limits
!write(*,*)
!
pix_lim(1,:) =  HUGE(1)
pix_lim(2,:) = -HUGE(1)
hkl_lim = adt_limits + adt_sigma     ! Set grand HKL limits
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
!write(*,*) ' ADT_LIMITS ', adt_limits
!write(*,*) ' ADT_NORMAL ', adt_plane
!write(*,*) ' ADT_STEPS  ', adt_steps
!write(*,*) ' ADT_ABS    ', adt_absc
!write(*,*) ' ADT_ORD    ', adt_ordi
!write(*,*) ' ADT_INTER  ', adt_inter
!
!write(*,*) ' ADT_LABS   ', adt_labsc
!write(*,*) ' ADT_LORD   ', adt_lordi
!write(*,*) ' ADT_LINTER ', adt_linter
!write(*,*)
!write(*,*)  ' TRANSI', adt_transi(1,:)
!write(*,*)  ' TRANSI', adt_transi(2,:)
!write(*,*)  ' TRANSI', adt_transi(3,:)
!!write(*,*)
!write(*,*)  ' TRANS ', adt_transm(1,:)
!write(*,*)  ' TRANS ', adt_transm(2,:)
!write(*,*)  ' TRANS ', adt_transm(3,:)
!!
!write(*,*) ' DIMENSIONS ', idims
!!
!!
!write(*,*) ' DIMENSION ', ik1%dims, ' RODS ', idims 
!write(*,*) ' PIXEL LIM ', pix_lim(:,1)
!write(*,*) ' PIXEL LIM ', pix_lim(:,2)
!write(*,*) ' PIXEL LIM ', pix_lim(:,3)
!read(*,*) i
!
allocate(ik2%datamap(idims(1), idims(2), idims(3)))
ik2%datamap = 0.0D0
!hkl = (/0,0,0/)
!ppixel = matmul(adt_transm,hkl)
!write(*,'(a,3(3(f7.2x,2x),2x),f7.2x)') ' HKL > ', hkl, ppixel, ppixel + (idims+1)/2, dot_product(adt_plane, hkl)
!hkl = (/0.0,0.0,0.5/)
!ppixel = matmul(adt_transm,hkl)
!write(*,'(a,3(3(f7.2x,2x),2x),f7.2x)') ' HKL > ', hkl, ppixel, ppixel + (idims+1)/2, dot_product(adt_plane, hkl)
!hkl = (/0.5,0.0,0.5/)
!ppixel = matmul(adt_transm,hkl)
!hkl = (/0.8,0.0,0.5/)
!ppixel = matmul(adt_transm,hkl)
!write(*,'(a,3(3(f7.2x,2x),2x),f7.2x)') ' HKL > ', hkl, ppixel, ppixel + (idims+1)/2, dot_product(adt_plane, hkl)
!
!write(*,*) ' ALLOCATED '
!read(*,*) i
do k=pix_lim(1,3), pix_lim(2,3)
   pixel(3) = real(k,kind=PREC_DP)
   do j=pix_lim(1,2), pix_lim(2,2)
      pixel(2) = real(j,kind=PREC_DP)
      do i=pix_lim(1,1), pix_lim(2,1)
         pixel(1) = real(i,kind=PREC_DP)
         hkl = matmul(adt_invmat_t, (pixel-adt_zero))
         ppixel = matmul(adt_transm,hkl)
         if(abs(ppixel(1)-nint(ppixel(1)))<=adt_sigma(1) .and. &
            abs(ppixel(1)-nint(ppixel(1)))<=adt_sigma(1)) then
            ii =  int(ppixel(1)) + (idims(1)+1)/2
            jj =  int(ppixel(2)) + (idims(2)+1)/2
            kk = nint(ppixel(3)) + (idims(3)+1)/2
            if(ii>1 .and. ii<=idims(1)  .and.   &
               jj>1 .and. jj<=idims(2)  .and.   &
               kk>1 .and. kk<=idims(3)) then
                  ik2%datamap(ii,jj,kk) = ik2%datamap(ii,jj,kk) + ik1%datamap(i,j,k)
            endif
         endif
      enddo
   enddo
enddo
!
!
!
ikk = iz
ik2%data_num = 0
if(mod(idims(3),2)==0) then
   ik2%layer = idims(3)/2
else
   ik2%layer = (idims(3)+1)/2
endif
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
ik2%vectors(:,1) = (ik2%corners(:,2) - ik2%corners(:,1))/(real(ik2%dims(1),kind=PREC_DP)-1.0D0)
ik2%vectors(:,2) = (ik2%corners(:,3) - ik2%corners(:,1))/(real(ik2%dims(2),kind=PREC_DP)-1.0D0)
ik2%vectors(:,3) = (ik2%corners(:,4) - ik2%corners(:,1))/(real(ik2%dims(3),kind=PREC_DP)-1.0D0)
ik2%cr_a0          = adt_di_cell(1:3)
ik2%cr_win         = adt_di_cell(4:6)
ik2%llims          = -adt_limits
ik2%steps(1)       = ik2%vectors(1,1)
ik2%steps(2)       = ik2%vectors(2,2)
ik2%steps(3)       = ik2%vectors(3,3)
ik2%steps_full(:,1)= ik2%vectors(:,1)
ik2%steps_full(:,2)= ik2%vectors(:,2)
ik2%steps_full(:,3)= ik2%vectors(:,3)
!write(*,*) ' IK2  dims ', ik2%dims, ' <> ', idims
!write(*,*) ' IK2 layer ', ik2%layer
!write(*,*) ' IK2 ll    ', ik2%corners(:,1)
!write(*,*) ' IK2 lr    ', ik2%corners(:,2)
!write(*,*) ' IK2 ul    ', ik2%corners(:,3)
!write(*,*) ' IK2 tl    ', ik2%corners(:,4)
!write(*,*) ' IK2 abs   ', ik2%vectors(:,1)
!write(*,*) ' IK2 ord   ', ik2%vectors(:,2)
!write(*,*) ' IK2 top   ', ik2%vectors(:,3)
!write(*,*) ' IK2 step  ', ik2%steps  (:  )
!write(*,*) ' IK2 step a', ik2%steps_full(:,1)
!write(*,*) ' IK2 step o', ik2%steps_full(:,2)
!write(*,*) ' IK2 step t', ik2%steps_full(:,3)
!write(*,*) ' IK2 llims ', ik2%llims
!write(*,*) ' IK2 cell  ', ik2%cr_a0, ik2%cr_win
!write(*,*)
!
ik2%infile = 'extracted_rods'
string = 'adt_rods.inte'
!
call local2data(ikk, ier_num, ier_typ, ik2%data_num, ik2%infile, ik2%layer,  &
     ik2%is_direct, ik2%ndims, ik2%dims, ik2%is_grid, ik2%has_dxyz,             &
     ik2%has_dval, ik2%corners, ik2%vectors, ik2%cr_a0, ik2%cr_win, ik2%x, ik2%y,     &
     ik2%z, ik2%dx, ik2%dy, ik2%dz, ik2%datamap, ik2%sigma, ik2%llims, ik2%steps,  &
     ik2%steps_full)
call dgl5_set_h5_is_ku(iz, ik2%data_num)
call dgl5_set_ku_is_h5(ik2%data_num, iz)
ku_ndims(iz) = ik2%ndims
call data2kuplot(ikk, string, lout)
!
call adt_clean(ik1)
call adt_clean(ik2)
!
end subroutine adt_extract_rods
!
!*******************************************************************************
!
subroutine adt_extract_planes
!-
!  Extract the selected planes into the next data set
!  The final hkl must fulfill the condition  [uvw] x [hkl] = N*steps with N integer
!+
!
use kuplot_mod
use kuplot_global
!
use errlist_mod
use lib_data_struc_type_mod
use lib_data_struc_h5 , only:data2local, local2data, dgl5_set_h5_is_ku, dgl5_set_ku_is_h5
use lib_metric_mod    , only:lib_d2r
use matrix_mod        , only:matinv
use support_mod       , only:oeffne
!
use lib_write_mod
!
implicit none
!
integer, parameter :: IWR = 95
character(len=PREC_STRING) :: string
integer                                    :: i,j,k      ! Loop indices
integer                                    :: ii,jj,kk   ! Loop indices
integer                                    :: ni, nj     ! Loop limits
integer                                    :: ikk        ! Kuplot data set result
integer                                    :: nref       ! Reflection counter
integer                   , dimension(3)   :: idims      ! Deimensions of resulting grid
integer                   , dimension(3)   :: ipixel     ! Coordinates in ADT data
real(kind=PREC_DP)        , dimension(3)   :: pixel      ! Coordinates in ADT data
real(kind=PREC_DP)        , dimension(3)   :: ppixel      ! Coordinates in ADT data
real(kind=PREC_DP)        , dimension(8)   :: pweight    ! Weight for the eight corners
real(kind=PREC_DP)        , dimension(3)   :: hkl
real(kind=PREC_DP)        , dimension(3)   :: hkl_lim    ! HKL Limits including sigma
real(kind=PREC_DP)                         :: rlayer     ! Current  [uvw] X [hkl] 
real(kind=PREC_DP)                         :: steps      ! Steps in [uvw] X [hkl] 
real(kind=PREC_DP)                         :: lsteps     ! Max Steps in [uvw] X [hkl] 
integer                   , dimension(2,3) :: pix_lim    ! PIXEL Limits including sigma
real(kind=PREC_DP), dimension(:,:,:), allocatable :: cweight
!
type bragg
   real(kind=PREC_DP)               :: weight
   real(kind=PREC_DP), dimension(3) :: com
   real(kind=PREC_DP), dimension(3) :: pixel
   real(kind=PREC_DP), dimension(3) :: hkl
   integer           , dimension(3) :: ihkl
end type bragg
type(bragg), dimension(:), allocatable :: bragg_list
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
if(ALL(adt_zero<0.0D0)) then
   adt_zero = ik1%dims/2
endif
!write(*,*) ' ADT_zero at : ', adt_zero
!
!  Copy temporary steps and maximum steps into local storage
!
steps  = adt_inter(1)     ! Steps in [uvw] x [hkl] 
lsteps = adt_linter(1)    ! Maxvalue of uvw] x [hkl]
!
! Transform [uvw] into [hkl]
!
call lib_d2r(adt_gten, adt_rten, adt_plane, hkl, adt_inter)
adt_inter = adt_inter*steps/dot_product(adt_plane, adt_inter)
!
!write(*,*) ' STEP, MAX ', steps, lsteps
!write(*,*) ' ADT_PLANE ', adt_plane
!write(*,*) ' ADT_INTER ', adt_inter
!write(*,*) ' UVWxHKL   ', dot_product(adt_plane, adt_inter)
!
! Build matrix FinalPixel => hkl
!
adt_transi(:,1) = adt_absc
adt_transi(:,2) = adt_ordi
adt_transi(:,3) = adt_inter
call matinv(adt_transi, adt_transm)  ! Build matrix hkl => FinalPixel
!
! Determine maximum final dimensions
!
if(nint(adt_labsc(0))==1) then
   idims(1) = nint(2*adt_labsc(1)) + 1
else
   hkl = adt_labsc(1:3)
   ppixel = matmul(adt_transm, hkl)
!   write(*,*) ' ABSCISASA  => ', ppixel
   idims(1) = nint(2*ppixel(1)) + 1
endif
if(nint(adt_lordi(0))==1) then
   idims(2) = nint(2*adt_lordi(1)) + 1
else
   hkl = adt_lordi(1:3)
   ppixel = matmul(adt_transm, hkl)
!   write(*,*) ' ORDINATE   => ', ppixel
   idims(2) = nint(2*ppixel(2)) + 1
endif
if(nint(adt_linter(0))==1) then
   idims(3) = nint(2*(lsteps/steps)) + 1
   adt_linter(1:3) = adt_inter*(lsteps/steps)/dot_product(adt_plane, adt_inter)
else
   hkl = adt_linter(1:3)
   ppixel = matmul(adt_transm, hkl)
!   write(*,*) ' INTERVAL   => ', ppixel
   idims(3) = nint(2*ppixel(3)) + 1
endif
!
!  Set hkl search limits to determine smallest volume in pixel space 
!
adt_limits(1) = max(adt_limits(1), maxval(adt_labsc), maxval(adt_lordi), maxval(adt_linter))
adt_limits(2) = max(adt_limits(2), maxval(adt_labsc), maxval(adt_lordi), maxval(adt_linter))
adt_limits(3) = max(adt_limits(3), maxval(adt_labsc), maxval(adt_lordi), maxval(adt_linter))
!
adt_limits = adt_labsc(1:3) + adt_lordi(1:3) + adt_linter(1:3)*steps
!write(*,*) 
!write(*,*) ' ADT_LABSC  ',adt_labsc 
!write(*,*) ' ADT_LORDI  ',adt_lordi 
!write(*,*) ' ADT_LINTER ',adt_linter
!write(*,*) ' ADT_LIMITS ',adt_limits
!write(*,*) 
!
pix_lim(1,:) =  HUGE(1)
pix_lim(2,:) = -HUGE(1)
hkl_lim = adt_limits + adt_sigma     ! Set grand HKL limits
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
!write(*,*) ' ADT_LIMITS ', adt_limits
!write(*,*) ' ADT_NORMAL ', adt_plane
!write(*,*) ' ADT_STEPS  ', adt_steps
!write(*,*) ' ADT_ABS    ', adt_absc
!write(*,*) ' ADT_ORD    ', adt_ordi
!write(*,*) ' ADT_INTER  ', adt_inter
!
!write(*,*) ' ADT_LABS   ', adt_labsc
!write(*,*) ' ADT_LORD   ', adt_lordi
!write(*,*) ' ADT_LINTER ', adt_linter
!write(*,*)
!write(*,*)  ' TRANSI', adt_transi(1,:)
!write(*,*)  ' TRANSI', adt_transi(2,:)
!write(*,*)  ' TRANSI', adt_transi(3,:)
!write(*,*)
!write(*,*)  ' TRANS ', adt_transm(1,:)
!write(*,*)  ' TRANS ', adt_transm(2,:)
!write(*,*)  ' TRANS ', adt_transm(3,:)
!hkl = (/0,1,0/)
!ppixel = matmul(adt_transm,hkl)
!write(*,'(a,3(3(f7.2x,2x),2x),f7.2x)') ' HKL > ', hkl, ppixel, ppixel  + (idims+1)/2, dot_product(adt_plane, hkl)
!hkl = (/0,0,0/)
!ppixel = matmul(adt_transm,hkl)
!write(*,'(a,3(3(f7.2x,2x),2x),f7.2x)') ' HKL > ', hkl, ppixel, ppixel + (idims+1)/2, dot_product(adt_plane, hkl)
!hkl = (/0.0,0.0,0.5/)
!ppixel = matmul(adt_transm,hkl)
!write(*,'(a,3(3(f7.2x,2x),2x),f7.2x)') ' HKL > ', hkl, ppixel, ppixel + (idims+1)/2, dot_product(adt_plane, hkl)
!hkl = (/0.0,0.0,1.0/)
!ppixel = matmul(adt_transm,hkl)
!write(*,'(a,3(3(f7.2x,2x),2x),f7.2x)') ' HKL > ', hkl, ppixel, ppixel + (idims+1)/2, dot_product(adt_plane, hkl)
!hkl = (/2.0,1.0,0.0/)
!ppixel = matmul(adt_transm,hkl)
!write(*,'(a,3(3(f7.2x,2x),2x),f7.2x)') ' HKL > ', hkl, ppixel, ppixel + (idims+1)/2, dot_product(adt_plane, hkl)
!
!write(*,*) ' DIMENSIONS ', idims
!
! Allocate final storage array
!
allocate(ik2%datamap(idims(1), idims(2), idims(3)))
ik2%datamap = 0.0D0
!
!write(*,*) ' ALLOCATED ', ubound(ik2%datamap)
!write(*,*) ' DIMENSION ', ik1%dims, ' PLANES', idims 
!write(*,*) ' PIXEL LIM ', pix_lim(:,1)
!write(*,*) ' PIXEL LIM ', pix_lim(:,2)
!write(*,*) ' PIXEL LIM ', pix_lim(:,3)
!write(*,*)
!write(*,*)  ' TRANSI', adt_invmat_t(1,:)
!write(*,*)  ' TRANSI', adt_invmat_t(2,:)
!write(*,*)  ' TRANSI', adt_invmat_t(3,:)
!write(*,*)
!!read(*,*) i
do k=pix_lim(1,3), pix_lim(2,3)
   pixel(3) = real(k,kind=PREC_DP)
   do j=pix_lim(1,2), pix_lim(2,2)
      pixel(2) = real(j,kind=PREC_DP)
      do i=pix_lim(1,1), pix_lim(2,1)
         pixel(1) = real(i,kind=PREC_DP)
         hkl = matmul(adt_invmat_t, (pixel-adt_zero))
         rlayer = dot_product(adt_plane, hkl)/steps
!if(abs(hkl(1)-0.0)<0.05 .and. abs(hkl(2)-0.0)<0.05) then
!if(abs(hkl(3)-0.0)<0.05                           ) then
!if(abs(hkl(3)/steps-nint(hkl(3)/steps))<0.02) then
ppixel = matmul(adt_transm, hkl)
            ii =  int(ppixel(1)) + (idims(1)+1)/2
            jj =  int(ppixel(2)) + (idims(2)+1)/2
            kk = nint(ppixel(3)) + (idims(3)+1)/2
!write(*,'(3i4,'':'',3(f7.2,2x),'':'',f7.2,2x, '':'',3(f7.2,2x),'':'',3i4)') &
!i,j,k, hkl, rlayer, ppixel, ii,jj,kk
!endif
!endif
         if(abs(rlayer-nint(rlayer))<=adt_sigma(1)) then
            ppixel = matmul(adt_transm, hkl)
            ii =  int(ppixel(1)) + (idims(1)+1)/2
            jj =  int(ppixel(2)) + (idims(2)+1)/2
            kk = nint(ppixel(3)) + (idims(3)+1)/2
            if(ii>1 .and. ii<=idims(1)  .and.   &
               jj>1 .and. jj<=idims(2)  .and.   &
               kk>1 .and. kk<=idims(3)) then
               ik2%datamap(ii,jj,kk) = ik2%datamap(ii,jj,kk) + ik1%datamap(i,j,k)
            endif
         endif
      enddo
   enddo
enddo
!
ikk = iz
ik2%data_num = 0
if(mod(idims(3),2)==0) then
   ik2%layer = idims(3)/2
else
   ik2%layer = (idims(3)+1)/2
endif
ik2%is_direct = .false.
ik2%ndims     = 3
ik2%dims      = idims
ik2%is_grid   = .true.
ik2%has_dxyz  = .false.
ik2%has_dval  = .false.
!
!
! Lower Left
ik2%corners(1,1) =  -adt_absc(1)*((idims(1)+1)/2-1)  - adt_ordi(1)*((idims(2)+1)/2-1) - adt_inter(1)*((idims(3)+1)/2-1)
ik2%corners(2,1) =  -adt_absc(2)*((idims(1)+1)/2-1)  - adt_ordi(2)*((idims(2)+1)/2-1) - adt_inter(2)*((idims(3)+1)/2-1)
ik2%corners(3,1) =  -adt_absc(3)*((idims(1)+1)/2-1)  - adt_ordi(3)*((idims(2)+1)/2-1) - adt_inter(3)*((idims(3)+1)/2-1)
! Lower Right
ik2%corners(1,2) =   adt_absc(1)*((idims(1)+1)/2-1)  - adt_ordi(1)*((idims(2)+1)/2-1) - adt_inter(1)*((idims(3)+1)/2-1)
ik2%corners(2,2) =   adt_absc(2)*((idims(1)+1)/2-1)  - adt_ordi(2)*((idims(2)+1)/2-1) - adt_inter(2)*((idims(3)+1)/2-1)
ik2%corners(3,3) =   adt_absc(3)*((idims(1)+1)/2-1)  - adt_ordi(3)*((idims(2)+1)/2-1) - adt_inter(3)*((idims(3)+1)/2-1)
! Upper Left
ik2%corners(1,3) =  -adt_absc(1)*((idims(1)+1)/2-1)  + adt_ordi(1)*((idims(2)+1)/2-1) - adt_inter(1)*((idims(3)+1)/2-1)
ik2%corners(2,3) =  -adt_absc(2)*((idims(1)+1)/2-1)  + adt_ordi(2)*((idims(2)+1)/2-1) - adt_inter(2)*((idims(3)+1)/2-1)
ik2%corners(3,3) =  -adt_absc(3)*((idims(1)+1)/2-1)  + adt_ordi(3)*((idims(2)+1)/2-1) - adt_inter(3)*((idims(3)+1)/2-1)
! Top   Left-1)
ik2%corners(1,4) =  -adt_absc(1)*((idims(1)+1)/2-1)  - adt_ordi(1)*((idims(2)+1)/2-1) + adt_inter(1)*((idims(3)+1)/2-1)
ik2%corners(2,4) =  -adt_absc(2)*((idims(1)+1)/2-1)  - adt_ordi(2)*((idims(2)+1)/2-1) + adt_inter(2)*((idims(3)+1)/2-1)
ik2%corners(3,4) =  -adt_absc(3)*((idims(1)+1)/2-1)  - adt_ordi(3)*((idims(2)+1)/2-1) + adt_inter(3)*((idims(3)+1)/2-1)
!
!
ik2%vectors(:,1) = adt_absc
ik2%vectors(:,2) = adt_ordi
ik2%vectors(:,3) = adt_inter
   
ik2%cr_a0          = adt_di_cell(1:3)
ik2%cr_win         = adt_di_cell(4:6)
ik2%llims          = -adt_limits
ik2%steps(1)       = maxval(abs(ik2%vectors(:,1)))
ik2%steps(2)       = maxval(abs(ik2%vectors(:,2)))
ik2%steps(3)       = maxval(abs(ik2%vectors(:,3)))
ik2%steps_full(:,1)= ik2%vectors(:,1)
ik2%steps_full(:,2)= ik2%vectors(:,2)
ik2%steps_full(:,3)= ik2%vectors(:,3)
!
!write(*,*) ' IK2  dims ', ik2%dims, ' <> ', idims
!write(*,*) ' IK2 layer ', ik2%layer
!write(*,*) ' IK2 ll    ', ik2%corners(:,1)
!write(*,*) ' IK2 lr    ', ik2%corners(:,2)
!write(*,*) ' IK2 ul    ', ik2%corners(:,3)
!write(*,*) ' IK2 tl    ', ik2%corners(:,4)
!write(*,*) ' IK2 abs   ', ik2%vectors(:,1)
!write(*,*) ' IK2 ord   ', ik2%vectors(:,2)
!write(*,*) ' IK2 top   ', ik2%vectors(:,3)
!write(*,*) ' IK2 step  ', ik2%steps  (:  )
!write(*,*) ' IK2 step a', ik2%steps_full(:,1)
!write(*,*) ' IK2 step o', ik2%steps_full(:,2)
!write(*,*) ' IK2 step t', ik2%steps_full(:,3)
!write(*,*) ' IK2 llims ', ik2%llims
!write(*,*) ' IK2 cell  ', ik2%cr_a0, ik2%cr_win
!write(*,*)
!
ik2%infile = 'extracted_planes'
string = 'adt_planes.inte'
call local2data(ikk, ier_num, ier_typ, ik2%data_num, ik2%infile, ik2%layer,  &
     ik2%is_direct, ik2%ndims, ik2%dims, ik2%is_grid, ik2%has_dxyz,             &
     ik2%has_dval, ik2%corners, ik2%vectors, ik2%cr_a0, ik2%cr_win, ik2%x, ik2%y,     &
     ik2%z, ik2%dx, ik2%dy, ik2%dz, ik2%datamap, ik2%sigma, ik2%llims, ik2%steps,  &
     ik2%steps_full)
call dgl5_set_h5_is_ku(iz, ik2%data_num)
call dgl5_set_ku_is_h5(ik2%data_num, iz)
ku_ndims(iz) = ik2%ndims
call data2kuplot(ikk, string, lout)
!
!
call adt_clean(ik1)
call adt_clean(ik2)
!
end subroutine adt_extract_planes
!
!*******************************************************************************
!
subroutine adt_extract_volume
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
integer                                    :: ni, nj     ! Loop limits
integer                                    :: ikk        ! Kuplot data set result
integer                                    :: nref       ! Reflection counter
integer                   , dimension(3)   :: idims      ! Deimensions of resulting grid
integer                   , dimension(3)   :: ipixel     ! Coordinates in ADT data
real(kind=PREC_DP)        , dimension(3)   :: pixel      ! Coordinates in ADT data
real(kind=PREC_DP)        , dimension(3)   :: ppixel      ! Coordinates in ADT data
real(kind=PREC_DP)        , dimension(8)   :: pweight    ! Weight for the eight corners
real(kind=PREC_DP)        , dimension(3)   :: hkl
real(kind=PREC_DP)        , dimension(3)   :: hkl_lim    ! HKL Limits including sigma
integer                   , dimension(2,3) :: pix_lim    ! PIXEL Limits including sigma
real(kind=PREC_DP), dimension(:,:,:), allocatable :: cweight
!
type bragg
   real(kind=PREC_DP)               :: weight
   real(kind=PREC_DP), dimension(3) :: com
   real(kind=PREC_DP), dimension(3) :: pixel
   real(kind=PREC_DP), dimension(3) :: hkl
   integer           , dimension(3) :: ihkl
end type bragg
type(bragg), dimension(:), allocatable :: bragg_list
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
!write(*,*) ' ADT_LIMITS ', adt_limits
!write(*,*) ' ADT_STEPS  ', adt_steps
idims = nint(2*adt_limits/adt_steps)+1
!write(*,*) ' DIMENSIONS ', idims
!
!
pixel = 0.25
      call adt_weight_dist(pixel, pweight)
!write(*,*) ' PWEIGHT ',pweight(1:4)
!write(*,*) ' PWEIGHT ',pweight(5:8)
!write(*,*) ' PWEIGHT ',sum(pweight)
   allocate(ik2%datamap(idims(1), idims(2), idims(3)))
   allocate(cweight     (idims(1), idims(2), idims(3)))
   ik2%datamap = 0.0D0
   cweight     = 0.0D0
!write(*,*) ' ALLOCATED ', ubound(ik2%datamap)
!write(*,*) ' DIMENSION ', ik1%dims, ' VOLU ', idims 
!write(*,*) ' PIXEL LIM ', pix_lim(:,1)
!write(*,*) ' PIXEL LIM ', pix_lim(:,2)
!write(*,*) ' PIXEL LIM ', pix_lim(:,3)
!read(*,*) i
do k=pix_lim(1,3), pix_lim(2,3)
   pixel(3) = real(k,kind=PREC_DP)
   do j=pix_lim(1,2), pix_lim(2,2)
      pixel(2) = real(j,kind=PREC_DP)
      do i=pix_lim(1,1), pix_lim(2,1)
         pixel(1) = real(i,kind=PREC_DP)
         hkl = matmul(adt_invmat_t, (pixel-adt_zero))
         if(real(nint(abs(hkl(3))),kind=PREC_DP)<=adt_limits(3)) then
            if(real(    (abs(hkl(2))),kind=PREC_DP)<=adt_limits(2)) then
               if(real(    (abs(hkl(1))),kind=PREC_DP)<=adt_limits(1)) then
                   ppixel = ((hkl   +adt_limits   )/adt_steps   )+1
                   call adt_weight_dist(ppixel, pweight)
                   ii =  int((hkl(1)+adt_limits(1))/adt_steps(1))+1
                   jj =  int((hkl(2)+adt_limits(2))/adt_steps(2))+1
                   kk =  int((hkl(3)+adt_limits(3))/adt_steps(3))+1
! if(k== 51 .and. j==101 .and. i==151) then
! write(*,*) ' HKL  ', hkl   
! write(*,*) 'Pixel ', ppixel
! write(*,*) ' PWEIGHT ',pweight(1:4)
! write(*,*) ' PWEIGHT ',pweight(5:8)
! write(*,*) ' PWEIGHT ',sum(pweight)
! write(*,*) ' II JJ KK', ii,jj,kk
! write(*,*) ' PREV    ', ik2%datamap(ii  ,jj  ,kk  ) , ' + ', ik1%datamap(i,j,k)
! endif
                  ik2%datamap(ii  ,jj  ,kk  ) = ik2%datamap(ii  ,jj  ,kk  ) + ik1%datamap(i,j,k) * pweight(1)
                  cweight    (ii  ,jj  ,kk  ) = cweight    (ii  ,jj  ,kk  ) +                      pweight(1)
                  ik2%datamap(ii+1,jj  ,kk  ) = ik2%datamap(ii+1,jj  ,kk  ) + ik1%datamap(i,j,k) * pweight(2)
                  cweight    (ii+1,jj  ,kk  ) = cweight    (ii+1,jj  ,kk  ) +                      pweight(2)
                  ik2%datamap(ii  ,jj+1,kk  ) = ik2%datamap(ii  ,jj+1,kk  ) + ik1%datamap(i,j,k) * pweight(3)
                  cweight    (ii  ,jj+1,kk  ) = cweight    (ii  ,jj+1,kk  ) +                      pweight(3)
                  ik2%datamap(ii+1,jj+1,kk  ) = ik2%datamap(ii+1,jj+1,kk  ) + ik1%datamap(i,j,k) * pweight(4)
                  cweight    (ii+1,jj+1,kk  ) = cweight    (ii+1,jj+1,kk  ) +                      pweight(4)
                  ik2%datamap(ii  ,jj  ,kk+1) = ik2%datamap(ii  ,jj  ,kk+1) + ik1%datamap(i,j,k) * pweight(5)
                  cweight    (ii  ,jj  ,kk+1) = cweight    (ii  ,jj  ,kk+1) +                      pweight(5)
                  ik2%datamap(ii+1,jj  ,kk+1) = ik2%datamap(ii+1,jj  ,kk+1) + ik1%datamap(i,j,k) * pweight(6)
                  cweight    (ii+1,jj  ,kk+1) = cweight    (ii+1,jj  ,kk+1) +                      pweight(6)
                  ik2%datamap(ii  ,jj+1,kk+1) = ik2%datamap(ii  ,jj+1,kk+1) + ik1%datamap(i,j,k) * pweight(7)
                  cweight    (ii  ,jj+1,kk+1) = cweight    (ii  ,jj+1,kk+1) +                      pweight(7)
                  ik2%datamap(ii+1,jj+1,kk+1) = ik2%datamap(ii+1,jj+1,kk+1) + ik1%datamap(i,j,k) * pweight(8)
                  cweight    (ii+1,jj+1,kk+1) = cweight    (ii+1,jj+1,kk+1) +                      pweight(8)
!                 ii = nint((hkl(1)+adt_limits(1))/adt_steps(1))+1
!                 jj = nint((hkl(2)+adt_limits(2))/adt_steps(2))+1
!                 kk = nint((hkl(3)+adt_limits(3))/adt_steps(3))+1
!write(*,'(a, 3i4,3f7.2,3i4,f10.2)') ' PIXEL, hkl ii', i,j,k, hkl, ii,jj,kk, ik1%datamap(i,j,k)
!                 ik2%datamap(ii,jj,kk) = ik2%datamap(ii,jj,kk) + ik1%datamap(i,j,k)
!if(k== 51 .and. j==101 .and. i==151) then
!write(*,*) ' AFTER   ', ik2%datamap(ii  ,jj  ,kk  ) !, cweight    (ii  ,jj  ,kk  )
!endif
               endif
            endif
         endif
      enddo
   enddo
enddo
!write(*,*) 'WEIGHT ', minval(cweight), maxval(cweight)
!write(*,*) 'DATA   ', minval(ik2%datamap), maxval(ik2%datamap)
!
where(cweight>0)
   ik2%datamap = ik2%datamap/cweight
end where
deallocate(cweight)

!
!write(*,*) ' DONE WITH LOOP '
!
!
ikk = iz
ik2%data_num = 0
if(mod(idims(3),2)==0) then
   ik2%layer = idims(3)/2
else
   ik2%layer = (idims(3)+1)/2
endif
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
ik2%vectors(:,1) = (ik2%corners(:,2) - ik2%corners(:,1))/(real(ik2%dims(1),kind=PREC_DP)-1.0D0)
ik2%vectors(:,2) = (ik2%corners(:,3) - ik2%corners(:,1))/(real(ik2%dims(2),kind=PREC_DP)-1.0D0)
ik2%vectors(:,3) = (ik2%corners(:,4) - ik2%corners(:,1))/(real(ik2%dims(3),kind=PREC_DP)-1.0D0)
ik2%cr_a0          = adt_di_cell(1:3)
ik2%cr_win         = adt_di_cell(4:6)
ik2%llims          = -adt_limits
ik2%steps(1)       = ik2%vectors(1,1)
ik2%steps(2)       = ik2%vectors(2,2)
ik2%steps(3)       = ik2%vectors(3,3)
ik2%steps_full(:,1)= ik2%vectors(:,1)
ik2%steps_full(:,2)= ik2%vectors(:,2)
ik2%steps_full(:,3)= ik2%vectors(:,3)
!write(*,*) ' IK2  dims ', ik2%dims, ' <> ', idims
!write(*,*) ' IK2 layer ', ik2%layer
!write(*,*) ' IK2 ll    ', ik2%corners(:,1)
!write(*,*) ' IK2 lr    ', ik2%corners(:,2)
!write(*,*) ' IK2 ul    ', ik2%corners(:,3)
!write(*,*) ' IK2 tl    ', ik2%corners(:,4)
!write(*,*) ' IK2 abs   ', ik2%vectors(:,1)
!write(*,*) ' IK2 ord   ', ik2%vectors(:,2)
!write(*,*) ' IK2 top   ', ik2%vectors(:,3)
!write(*,*) ' IK2 step  ', ik2%steps  (:  )
!write(*,*) ' IK2 step a', ik2%steps_full(:,1)
!write(*,*) ' IK2 step o', ik2%steps_full(:,2)
!write(*,*) ' IK2 step t', ik2%steps_full(:,3)
!write(*,*) ' IK2 llims ', ik2%llims
!!write(*,*) ' IK2 cell  ', ik2%cr_a0, ik2%cr_win
!write(*,*)
!
ik2%infile = 'extracted_volume'
string = 'adt_volume.inte'
!
call local2data(ikk, ier_num, ier_typ, ik2%data_num, ik2%infile, ik2%layer,  &
     ik2%is_direct, ik2%ndims, ik2%dims, ik2%is_grid, ik2%has_dxyz,             &
     ik2%has_dval, ik2%corners, ik2%vectors, ik2%cr_a0, ik2%cr_win, ik2%x, ik2%y,     &
     ik2%z, ik2%dx, ik2%dy, ik2%dz, ik2%datamap, ik2%sigma, ik2%llims, ik2%steps,  &
     ik2%steps_full)
call dgl5_set_h5_is_ku(iz, ik2%data_num)
call dgl5_set_ku_is_h5(ik2%data_num, iz)
ku_ndims(iz) = ik2%ndims
call data2kuplot(ikk, string, lout)
!
!
call adt_clean(ik1)
call adt_clean(ik2)
!
end subroutine adt_extract_volume
!
!*******************************************************************************
!
subroutine adt_extract_bragg
!-
!  Extract the integrated Bragg reflections into a file 
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
integer                                    :: ni, nj     ! Loop limits
integer                                    :: ikk        ! Kuplot data set result
integer                                    :: nref       ! Reflection counter
integer                   , dimension(3)   :: idims      ! Deimensions of resulting grid
integer                   , dimension(3)   :: ipixel     ! Coordinates in ADT data
real(kind=PREC_DP)        , dimension(3)   :: pixel      ! Coordinates in ADT data
real(kind=PREC_DP)        , dimension(3)   :: ppixel      ! Coordinates in ADT data
real(kind=PREC_DP)        , dimension(8)   :: pweight    ! Weight for the eight corners
real(kind=PREC_DP)        , dimension(3)   :: hkl
real(kind=PREC_DP)        , dimension(3)   :: hkl_lim    ! HKL Limits including sigma
integer                   , dimension(2,3) :: pix_lim    ! PIXEL Limits including sigma
real(kind=PREC_DP), dimension(:,:,:), allocatable :: cweight
!
type bragg
   real(kind=PREC_DP)               :: weight
   real(kind=PREC_DP), dimension(3) :: com
   real(kind=PREC_DP), dimension(3) :: pixel
   real(kind=PREC_DP), dimension(3) :: hkl
   integer           , dimension(3) :: ihkl
end type bragg
type(bragg), dimension(:), allocatable :: bragg_list
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
!
if(ALL(adt_zero<0.0D0)) then
   adt_zero = ik1%dims/2
endif
!write(*,*) ' ADT_zero at : ', adt_zero
!
!
nref = 0
if(adt_log) then
   call oeffne(IWR, adt_logfile, 'unknown')
   write(IWR,'(a, 43x,a,17x,a,17x,a,13x,a)') '#', 'Observed', 'Calculated', 'Calculated', 'Nominal'
   write(IWR,'(a6,1x,a,1x,a,8x,a,12x,a,21x,a,23x,a,16x,a)') '#   Nr',        &
   'Target', 'dummy', 'Weight', 'Pixel', 'Pixel', 'hkl', 'H   K   L'
endif
!  allocate(bragg_list((2*adt_limits(1)+1)*(2*adt_limits(2)+1)*(2*adt_limits(3)+1)))
!write(*,*) 'ADT LIMITS ', adt_limits
!write(*,*) 'ADT ZERO   ', adt_zero  
!write(*,*) ' ADT_STEPS ', adt_steps
!write(*,*) ' LOGFILE   ', adt_logfile(1:len_trim(adt_logfile))
!
do k=-nint(adt_limits(3)), nint(adt_limits(3))
   do j=-nint(adt_limits(2)), nint(adt_limits(2))
      do i=-nint(adt_limits(1)), nint(adt_limits(1))
         if(adt_lattice_cond(i,j,k)) then
            hkl(3) = real(k, kind=PREC_DP)
            hkl(2) = real(j, kind=PREC_DP)
            hkl(1) = real(i, kind=PREC_DP)
            pixel = matmul(adt_omat_t, hkl) + adt_zero
            if(pixel(1)>8.0 .and. pixel(1)<ik1%dims(1)-8.0  .and. &
               pixel(2)>8.0 .and. pixel(2)<ik1%dims(2)-8.0  .and. &
               pixel(3)>8.0 .and. pixel(3)<ik1%dims(3)-8.0  ) then  
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
!
call adt_clean(ik1)
!
end subroutine adt_extract_bragg
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
logical function adt_lattice_cond(h,k,l)
!
implicit none
!
integer, intent(in) :: h
integer, intent(in) :: k
integer, intent(in) :: l
!
adt_lattice_cond = .false.
!
if(adt_lattice == ADT_PRIMI) then
   adt_lattice_cond = .true.
elseif(adt_lattice == ADT_CENT_A) then
   adt_lattice_cond = mod(k+l,2)==0
elseif(adt_lattice == ADT_CENT_B) then
   adt_lattice_cond = mod(h+l,2)==0
elseif(adt_lattice == ADT_CENT_C) then
   adt_lattice_cond = mod(h+k,2)==0
elseif(adt_lattice == ADT_CENT_I) then
   adt_lattice_cond = mod(h+k+l,2)==0
elseif(adt_lattice == ADT_CENT_F) then
   adt_lattice_cond = mod(h+k,2)==0 .and. mod(h+l,2)==0 .and. mod(k+l,2)==0
endif
!
end function adt_lattice_cond
!
!*******************************************************************************
!
subroutine adt_weight_dist(pix, weight)
!-
!  Calculate a weight for the next eight corners
!+
!
use precision_mod
!
implicit none
!
real(kind=PREC_DP), dimension(3), intent(in)  :: pix
real(kind=PREC_DP), dimension(8), intent(out) :: weight
!
real(kind=PREC_DP), dimension(3) :: pix_fr
!
pix_fr = pix - real(int(pix), kind=PREC_DP)
!
weight(1) = 1.0D0 - pix_fr(1) + 1.0D0 - pix_fr(2) + 1.0 - pix_fr(3)
weight(2) =         pix_fr(1) + 1.0D0 - pix_fr(2) + 1.0 - pix_fr(3)
weight(3) = 1.0D0 - pix_fr(1)         + pix_fr(2) + 1.0 - pix_fr(3)
weight(4) =         pix_fr(1)         + pix_fr(2) + 1.0 - pix_fr(3)
weight(5) = 1.0D0 - pix_fr(1) + 1.0D0 - pix_fr(2) +       pix_fr(3)
weight(6) =         pix_fr(1) + 1.0D0 - pix_fr(2) +       pix_fr(3)
weight(7) = 1.0D0 - pix_fr(1)         + pix_fr(2) +       pix_fr(3)
weight(8) =         pix_fr(1)         + pix_fr(2) +       pix_fr(3)
weight = weight / 12.0D0
!
end subroutine adt_weight_dist
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
adt_zero = -1.0D0
adt_plane  = (/          0.0D0 , 0.0D0 , 1.0D0 /)
adt_absc   = (/          1.0D-2, 0.0D0 , 0.0D0 /)
adt_ordi   = (/          0.0D0 , 1.0D-2, 0.0D0 /)
adt_inter  = (/          0.0D0 , 0.0D0 , 1.0D0 /)
adt_labsc  = (/ 3.0D0,   3.0D0 , 0.0D0 , 0.0D0 /)
adt_lordi  = (/ 3.0D0,   0.0D0 , 3.0D0 , 0.0D0 /)
adt_linter = (/ 3.0D0,   0.0D0 , 0.0D0 , 3.0D0 /)
!
end subroutine adt_reset
!
!*******************************************************************************
!
end module kuplot_adt_mod
