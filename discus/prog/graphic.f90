MODULE output_menu
!
CONTAINS
!
!*****7*************************************************************************
!
SUBROUTINE do_niplps (linverse) 
!-                                                                      
!     This sublevel contains all routines used to write the output      
!     of the Fourier transform/Patterson to an output file in           
!     various formats.                                                  
!+                                                                      
USE crystal_mod
USE discus_config_mod 
USE diffuse_mod 
USE nexus_discus
!USE discus_mrc
USE discus_xplor
USE vtk_mod
USE output_mod 
use powder_mod           , only:pow_lperiod
USE powder_write_mod
use powder_pdf_hist_mod  , only:pow_pdf_hist_prep_period
USE chem_aver_mod
USE qval_mod
!
USE ber_params_mod
USE build_name_mod
USE calc_expr_mod
USE doact_mod 
USE do_eval_mod
USE do_wait_mod
USE errlist_mod 
use gen_hdf_write_mod
USE get_params_mod
USE learn_mod 
USE lib_help
USE lib_do_operating_mod
USE lib_echo
USE lib_errlist_func
USE lib_length
USE lib_macro_func
use lib_mrc_mod
USE class_macro_internal
USE take_param_mod
USE precision_mod
USE prompt_mod 
USE str_comp_mod
USE sup_mod
USE support_mod
!
use fourier_sup
!
use lib_trans_mod
!                                                                       
IMPLICIT none 
LOGICAL, INTENT(IN) :: linverse
!                                                                       
INTEGER, PARAMETER :: maxp = 11 
INTEGER, PARAMETER :: MAXFORM = 14
integer, parameter :: NEW = 0
integer, parameter :: OLD = 1
integer, parameter :: ADD = 2
!
integer, parameter :: PDF3D_NORMAL = 0
integer, parameter :: PDF3D_SHARP  = 1
!                                                                       
CHARACTER(LEN=8) :: befehl 
CHARACTER(LEN=LEN(prompt)) :: orig_prompt
!CHARACTER(LEN=14) :: cvalue (0:15) 
CHARACTER(LEN=22) :: cgraphik (0:MAXFORM) 
CHARACTER(LEN=PREC_STRING) :: infile 
CHARACTER(LEN=PREC_STRING) :: zeile 
CHARACTER(LEN=PREC_STRING) :: line, cpara (maxp) 
INTEGER, DIMENSION(MAXP) :: lpara  !(maxp) 
INTEGER :: ix, iy, ianz, value, lp, length, lbef 
integer :: i,j,k,l, ii     ! Dummy indices
integer :: pdf3d_mode  !  3D-PDF mode "normal" / "sharp"
INTEGER :: indxg 
LOGICAL :: laver, lread
LOGICAL :: l_val_limited=.FALSE.    ! for 3D PDF do we limit d-star
REAL(KIND=PREC_DP) :: xmin, ymin, xmax, ymax 
REAL(KIND=PREC_DP), DIMENSION(MAXP) :: werte
REAL(KIND=PREC_DP), DIMENSION(:,:,:), allocatable :: qvalues
REAL(KIND=PREC_DP)                  :: valmax
REAL(KIND=PREC_DP)                  ::  dsmax = 0.0
!                                                                          ! All new are dummy, not used in output menu
character(len=PREC_STRING)                                :: new_outfile   ! New file name  from transformed HDF5 data
integer                   , dimension(3)                  :: new_inc       ! New dimensions from transformed HDF5 data
real(kind=PREC_DP)        , dimension(3,4)                :: new_eck       ! New corners    from transformed HDF5 data
real(kind=PREC_DP)        , dimension(3,3)                :: new_vi        ! New vectors    from transformed HDF5 data
real(kind=PREC_DP)        , dimension(:,:,:), allocatable :: new_qvalues   ! New data       from transformed HDF5 data
! 
INTEGER, PARAMETER :: NOPTIONAL = 11
INTEGER, PARAMETER :: O_SCALE   = 1                  ! Scale factor for output values
INTEGER, PARAMETER :: O_MAXVAL  = 2                  ! Current SCALE for maxvalue
INTEGER, PARAMETER :: O_DSMAX   = 3                  ! Maximum d-star
INTEGER, PARAMETER :: O_QMAX    = 4                  ! Maximum Q
INTEGER, PARAMETER :: O_HKLMAX  = 5                  ! Maximum hkl vector
INTEGER, PARAMETER :: O_PATT    = 6                  ! Optional Patteron overlay for Vesta
INTEGER, PARAMETER :: O_SPATT   = 7                  ! Optional Patteron overlay for Vesta
INTEGER, PARAMETER :: O_DPATT   = 8                  ! Optional Patteron overlay for Vesta
INTEGER, PARAMETER :: O_MODE    = 9                  ! Mode if written into KUPLOT
INTEGER, PARAMETER :: O_SHARP   =10                  ! 3D-PDF type normal/sharpened
INTEGER, PARAMETER :: O_TRANS   = 11                 ! Transform 3D data into different orientation
CHARACTER(LEN=   6), DIMENSION(NOPTIONAL) :: oname   !Optional parameter names
character(len=8)                          :: cpatt   ! Optional patterson overlay
character(len=PREC_STRING)                :: spatt   ! Atoms selected for Patterson overlay
character(len=PREC_STRING)                :: dpatt   ! Atoms deselected for Patterson overlay
CHARACTER(LEN=MAX(PREC_STRING,LEN(line))), DIMENSION(NOPTIONAL) :: opara   !Optional parameter strings returned
INTEGER            , DIMENSION(NOPTIONAL) :: loname  !Lenght opt. para name
INTEGER            , DIMENSION(NOPTIONAL) :: lopara  !Lenght opt. para name returned
LOGICAL            , DIMENSION(NOPTIONAL) :: lpresent!opt. para is present
REAL(KIND=PREC_DP) , DIMENSION(NOPTIONAL) :: owerte   ! Calculated values
INTEGER, PARAMETER                        :: ncalc = 1 ! Number of values to calculate
integer, parameter, dimension(3) :: nxyzstart = (/0,0,0/)
!
!                                                                       
DATA oname  / 'scale', 'maxval', 'dsmax ', 'qmax  ', 'hklmax', 'patt'  ,'sel', 'des', 'mode', 'type','trans' /
DATA loname /  5,       6      ,  5      ,  4      ,  6      ,  4      , 3   ,  3   ,  4    ,  4    , 5     /
DATA cgraphik / 'Standard', 'Postscript', 'Pseudo Grey Map', 'Gnuplot', &
                'Portable Any Map', 'Powder Pattern', 'SHELX',          &
                'SHELXL List 5', 'SHELXL List 5 real HKL' ,             &
                '3d', 'nexus', 'vtk', 'MRC', 'HDF5', 'XPLOR' /                                    
!DATA cvalue / 'undefined     ', 'Intensity     ', 'Amplitude     ',&
!              'Phase angle   ', 'Real Part     ', 'Imaginary Part',&
!              'Random Phase  ', 'S(Q)          ', 'F(Q)          ',&
!              'f2aver = <f^2>', 'faver2 = <f>^2', 'faver = <f>   ',&
!              'Normal Inten  ', 'I(Q)          ', 'PDF           ',&
!              '3DPDF         '                                     &
!            /
!
DATA value / 1 / 
DATA laver / .false. / 
!
opara  = (/'1.000 ', 'data  ', '0.000 ', '0.000 ', '0.000 ', 'none  ', 'all   ', 'none  ', 'new   ', 'normal' , 'no    '/)
lopara = (/  5     ,  4      ,  5      ,  5      ,  5      ,  4      ,  3      ,  4      ,  3      ,  6       ,  2      /)
owerte = (/  1.000 , -1.000  ,  0.000  ,  0.000  ,  0.000  ,  0.00   ,  0.000  ,  0.000  ,  0.00   ,  0.00    ,  0.0    /)
out_scale = 1.0D0
zmin = ps_low * diffumax 
zmax = ps_high * diffumax 
orig_prompt = prompt
prompt = prompt (1:len_str (prompt) ) //'/output' 
10 CONTINUE 
!                                                                       
CALL no_error 
!                                                                       
CALL get_cmd (line, length, befehl, lbef, zeile, lp, prompt) 
main_if: IF (ier_num.eq.0) THEN 
   IF (line (1:1)  == ' '.or.line (1:1)  == '#' .or.   & 
       line == char(13) .or. line(1:1) == '!'  ) THEN
      IF(linteractive .or. lmakro) THEN
         GOTO 10
      ELSE
         RETURN
      ENDIF
   ENDIF
!                                                                       
!     search for "="                                                    
!                                                                       
   indxg = index (line, '=') 
   IF (indxg.ne.0.AND..NOT. (str_comp (befehl, 'echo',   2, lbef, 4) ) &
                 .AND..NOT. (str_comp (befehl, 'system', 2, lbef, 6) )    &
                 .AND..NOT. (str_comp (befehl, 'help',   2, lbef, 4) .OR. &
                             str_comp (befehl, '?   ',   2, lbef, 4) )    &
                 .AND. INDEX(line,'==') == 0                            ) THEN
!
!     --evaluate an expression and assign the value to a variabble       
!                                                                       
      CALL do_math (line, indxg, length) 
   ELSE 
!                                                                       
!------ execute a macro file                                            
!                                                                       
      IF (befehl (1:1) .eq.'@') THEN 
         IF (length.ge.2) THEN 
            line(1:length-1) = line(2:length)
            line(length:length) = ' '
            length = length - 1
            CALL file_kdo(line, length)
         ELSE 
            ier_num = - 13 
            ier_typ = ER_MAC 
         ENDIF 
!                                                                       
!     continues a macro 'continue'                                      
!                                                                       
      ELSEIF (str_comp (befehl, 'continue', 1, lbef, 8) ) THEN 
         CALL macro_continue (zeile, lp) 
!                                                                       
!------ Echo a string, just for interactive check in a macro 'echo'     
!                                                                       
      ELSEIF (str_comp (befehl, 'echo', 2, lbef, 4) ) THEN 
         CALL echo (zeile, lp) 
!                                                                       
!     Evaluate an expression, just for interactive check 'eval'         
!                                                                       
      ELSEIF (str_comp (befehl, 'evaluate', 2, lbef, 8) ) THEN 
         CALL do_eval (zeile, lp, .TRUE.) 
!                                                                       
!     Terminate output 'exit'                                           
!                                                                       
      ELSEIF (str_comp (befehl, 'exit', 2, lbef, 4) ) THEN 
         GOTO 9999 
!                                                                       
!     Determine format for output 'format'                              
!                                                                       
      ELSEIF (str_comp (befehl, 'format', 1, lbef, 6) ) THEN 
         CALL get_params (zeile, ianz, cpara, lpara, maxp, lp) 
         IF (ier_num.eq.0) THEN 
            CALL get_optional(ianz, MAXP, cpara, lpara, NOPTIONAL,  ncalc, &
                              oname, loname, opara, lopara, lpresent, owerte)
            IF (ianz.eq.1.or.ianz.eq.2.or.ianz==5) THEN 
!                                                                       
!     ------Switch output type to ASCII 3D  '3d'                      
!                                                                       
               IF(str_comp(cpara(1),'3d',2,lpara(1),2)) THEN                                        
                  ityp = 9 
!                                                                       
!     ------Switch output type to GNUPLOT 'gnup'                        
!                                                                       
               ELSEIF(str_comp(cpara(1),'gnup',1,lpara(1),4)) THEN                                             
                  ityp = 3 
!                                                                       
!     ------Switch output type to pgm 'pgm'                             
!                                                                       
               ELSEIF(str_comp(cpara(1),'pgm ',2,lpara(1),4)) THEN                                        
                  ityp = 2 
!                                                                       
!     ------Switch output type to postscript 'post'                     
!                                                                       
               ELSEIF(str_comp(cpara(1),'post',3,lpara(1),4)) THEN                                        
                  ityp = 1 
!                                                                       
!     ------Switch output type to powder pattern 'powd'                 
!                                                                       
               ELSEIF(str_comp(cpara(1),'powd',3,lpara(1),4)) THEN                                        
                  ityp = 5 
                  IF (ianz >= 2) THEN 
                     IF(str_comp(cpara(2), 'tth', 2, lpara(2), 3)) THEN
                        cpow_form = 'tth' 
                     ELSEIF(str_comp(cpara(2), 'q', 1, lpara(2), 1)) THEN
                        cpow_form = 'q  ' 
                     ELSEIF(str_comp(cpara(2), 'tof', 3, lpara(2), 3)) THEN
                        cpow_form = 'tof' 
!                    ELSEIF(str_comp(cpara(2), 'stl', 2, lpara(2), 3) ) THEN
!                       cpow_form = 'stl' 
!                    ELSEIF(str_comp(cpara(2), 'dst', 2, lpara(2), 3)) THEN
!                       cpow_form = 'dst' 
!                    ELSEIF(str_comp(cpara(2), 'lop', 2, lpara(2), 3) ) THEN
!                       cpow_form = 'lop' 
                     ELSE 
                        ier_num = - 6 
                        ier_typ = ER_COMM 
                        exit main_if
                     ENDIF 
                     out_user_limits      = .false.
                     IF(ianz==5) THEN
                        cpara(1:2)='0'
                        lpara(1:2)= 1 
                        CALL ber_params (ianz, cpara, lpara, werte, maxp)
                        IF(ier_num==0) THEN
                           out_user_inc         = 0
                           out_user_values(1:3) = nint(werte(3:5)*1.0D5)/1.0D5
                           out_user_inc(1)      = nint((out_user_values(2)-out_user_values(1))/out_user_values(3) ) + 1
                           out_user_limits      = .true.
                        ENDIF
                     ENDIF
                  ENDIF 
!                                                                       
!     ------Switch output type to PDF 'pdf'                 
!                                                                       
               ELSEIF(str_comp(cpara(1),'pdf',3,lpara(1),3)) THEN                                        
                  ityp = 5 
                  IF (ianz >= 2) THEN 
                     IF(str_comp(cpara(2), 'r', 1, lpara(2), 1)) THEN
                        cpow_form = 'r  ' 
                     ELSE 
                        ier_num = -6 
                        ier_typ = ER_COMM 
                     ENDIF 
                     out_user_limits      = .false.
                     IF(ianz==5) THEN
                        cpara(1:2)='0'
                        lpara(1:2)= 1 
                        CALL ber_params (ianz, cpara, lpara, werte, maxp)
                        IF(ier_num==0) THEN
                           out_user_values(1:3) = werte(3:5)
                           out_user_limits      = .true.
                        ENDIF
                     ENDIF
                  ENDIF 
!                                                              
!     ------Switch output type to ppm 'ppm'                             
!                                                                       
               ELSEIF(str_comp(cpara(1),'ppm ',2,lpara(1),4)) THEN                                        
                  ityp = 4 
!                                                                       
!     ------Switch output type to standard  'stan'                      
!                                                                       
               ELSEIF(str_comp(cpara(1),'stan',2,lpara(1),4)) THEN                                        
                  ityp = 0 
!                                                                       
!     ------Switch output type to Shelx 'shel', or 'hklf4'              
!                                                                       
               ELSEIF(str_comp(cpara(1),'shelx',2,lpara(1),5)) THEN                                        
                  ityp = 6 
                  CALL form_optional(lpresent(O_MAXVAL), O_MAXVAL, NOPTIONAL, opara, &
                       lopara, werte)
                  IF(ier_num == 0) THEN
                     valmax = werte(O_MAXVAL)
                  ENDIF
               ELSEIF(str_comp(cpara(1),'hklf4',2,lpara(1),5)) THEN                                        
                  ityp = 6 
                  CALL form_optional(lpresent(O_MAXVAL), O_MAXVAL, NOPTIONAL, opara, &
                       lopara, werte)
                  IF(ier_num == 0) THEN
                     valmax = werte(O_MAXVAL)
                  ENDIF
!                                                                       
!     ------Switch output type to Shelx LIST 5   'list5'                
!                                                                       
               ELSEIF(str_comp(cpara(1),'list5',2,lpara(1),5)) THEN                                        
                  ityp = 7 
                  CALL form_optional(lpresent(O_MAXVAL), O_MAXVAL, NOPTIONAL, opara, &
                       lopara, werte)
                  IF(ier_num == 0) THEN
                     valmax = werte(O_MAXVAL)
                  ENDIF
!                                                                       
!     ------Switch output type to Shelx LIST 5   'list9'                
!                                                                       
               ELSEIF(str_comp(cpara(1),'list9',2,lpara(1),5)) THEN                                        
                  ityp = 8 
!                                                                       
!     ------Switch output type to NeXus format   'nexus'                
!                                                                       
               ELSEIF(str_comp(cpara(1),'nexus',2,lpara(1),5)) THEN                                        
                  ityp = 10 
!                                                                       
!     ------Switch output type to VTK format   'vtk'
!                                                                       
               ELSEIF(str_comp(cpara(1),'vtk',2,lpara(1),3)) THEN
                  ityp = 11
!                                                                       
!     ------Switch output type to MRC   format   'mrc'                
!                                                                       
               ELSEIF (str_comp(cpara(1), 'mrc', 2, lpara(1), 3) ) THEN                                        
                  ityp = 12 
!                                                                       
!     ------Switch output type to HDF5  format   'hdf5'                
!                                                                       
               ELSEIF (str_comp(cpara(1), 'hdf5', 4, lpara(1), 4) ) THEN                                        
                  if(lpresent(O_SCALE)) then
                     if(lpresent(O_MAXVAL)) then
                        ier_num = -6
                        ier_typ = ER_FORT
                        ier_msg(1) = 'Only one of ''scale:'' and ''maxval'' is allowed'
                        exit main_if
                     else
                        valmax =  -1.0D0
                        out_scale = owerte(O_SCALE)
                        ityp = 13 
                     endif
                  else
                  CALL form_optional(lpresent(O_MAXVAL), O_MAXVAL, NOPTIONAL, opara, &
                       lopara, werte)
                  IF(ier_num == 0) THEN
                     if(nint(werte(O_MAXVAL))==-2) then
                        valmax = 10000.0D0
                     else
                        valmax = werte(O_MAXVAL)
                     endif
                     ityp = 13 
                  ENDIF
                  endif
!                                                                       
!     ------Switch output type to XPLOR format   'xplor'                
!                                                                       
               ELSEIF (str_comp(cpara(1), 'xplor', 2, lpara(1), 5) ) THEN                                        
                  ityp = 14 
!                                                                       
!     ------Switch output type to Vesta format   'vesta'                
!                                                                       
               ELSEIF (str_comp(cpara(1), 'vesta', 2, lpara(1), 5) ) THEN                                        
                  cpatt = opara(O_PATT)(1:min(len(cpatt),len_trim(opara(O_PATT))))
                  if(opara(O_SPATT)=='none' .or. opara(O_SPATT)=='all' .or. index(opara(O_SPATT),',')==0) then
                     spatt = opara(O_SPATT)                             ! Selected atoms  for Patterson overlay
                  else
                     if(opara(O_SPATT)(1:1)/='[' .or. opara(O_SPATT)(lopara(O_SPATT):lopara(O_SPATT))/=']') then
                     ier_num = -9
                     ier_typ = ER_FORT
                     ier_msg(1) = 'Multiple optional values must be '
                     ier_msg(2) = 'enclosed by []'
                     exit main_if
                     endif
                     spatt = opara(O_SPATT)(2:lopara(O_SPATT)-1)        ! Selected atoms  for Patterson overlay
                  endif
                  if(opara(O_DPATT)=='none' .or. opara(O_DPATT)=='all' .or. index(opara(O_DPATT),',')==0) then
                     dpatt = opara(O_DPATT)                             ! Selected atoms  for Patterson overlay
                  else
                  if(opara(O_DPATT)(1:1)/='[' .or. opara(O_DPATT)(lopara(O_DPATT):lopara(O_DPATT))/=']') then
                        ier_num = -9
                        ier_typ = ER_FORT
                        ier_msg(1) = 'Multiple optional values must be '
                        ier_msg(2) = 'enclosed by []'
                        exit main_if
                     endif
                     dpatt = opara(O_DPATT)(2:lopara(O_DPATT)-1)        ! Deselected atoms  for Patterson overlay
                  endif
                  ityp = 15 
               ELSE
                  ier_num = - 9 
                  ier_typ = ER_APPL 
               ENDIF 
            ELSE 
               ier_num = - 6 
               ier_typ = ER_COMM 
            ENDIF 
         ENDIF 
!                                                                       
!     help on output 'help'                                             
!                                                                       
      ELSEIF (str_comp (befehl, 'help', 2, lbef, 4) .or.str_comp (befehl&
        &, '?   ', 1, lbef, 4) ) THEN                                      
         IF (str_comp (zeile, 'errors', 2, lp, 6) ) THEN 
            lp = lp + 7 
            CALL do_hel ('discus '//zeile, lp) 
         ELSE 
            lp = lp + 14 
            CALL do_hel ('discus output '//zeile, lp) 
         ENDIF 
!                                                                       
!     read an old output file (only for standard file type' 'inpu'      
!                                                                       
      ELSEIF (str_comp (befehl, 'input', 1, lbef, 5) ) THEN 
         CALL get_params (zeile, ianz, cpara, lpara, maxp, lp) 
         IF (ier_num.eq.0) THEN 
            infile = cpara (1) 
            lread = .true. 
            CALL oeffne (1, infile, 'old') 
            IF (ier_num.eq.0) THEN 
               READ (1, * ) out_inc (1), out_inc (2) 
               READ (1, * ) xmin, xmax, ymin, ymax 
               READ (1, * ) zmax 
               zmin = zmax 
               BACKSPACE (1) 
!                                                                       
               DO iy = 1, out_inc (2) 
!              READ (1, * ) (dsi ( (ix - 1) * out_inc (2) + iy),  &
!              ix = 1, out_inc (1) )                              
               read(1,*) (dsi(ix,iy,1), ix = 1, out_inc (1))
!              DO ix = 1, out_inc (1) 
!                 zmax = max(zmax, REAL(dsi((ix - 1) * out_inc(2) + iy), kind=PREC_DP))
!                 zmin = min(zmin, REAL(dsi((ix - 1) * out_inc(2) + iy), kind=PREC_DP))
!              ENDDO 
               ENDDO 
               zmax = max(zmax,maxval(dsi))
               zmin = min(zmin,minval(dsi))
               WRITE (output_io, 1015, advance='no') zmin, zmax 
               READ ( *, *, end = 20) zmin, zmax 
   20                CONTINUE 
            ENDIF 
            CLOSE (1) 
         ENDIF 
!                                                                       
!     define name of output file 'outf'                                 
!                                                                       
      ELSEIF (str_comp (befehl, 'outfile', 1, lbef, 7) ) THEN 
         CALL get_params (zeile, ianz, cpara, lpara, maxp, lp) 
         IF (ier_num.eq.0) THEN 
            CALL get_optional(ianz, MAXP, cpara, lpara, NOPTIONAL,  ncalc, &
                              oname, loname, opara, lopara, lpresent, owerte)
            CALL do_build_name(ianz, cpara, lpara, werte, maxp, 1)
            IF (ier_num.eq.0) THEN 
               outfile = cpara(1)(1:lpara(1))
               if(opara(O_MODE)=='new') then
                  out_mode =  NEW
               elseif(opara(O_MODE)=='old') then
                  out_mode =  OLD
               elseif(opara(O_MODE)=='add') then
                  out_mode =  ADD
               endif
            ENDIF 
         ENDIF 
!                                                                       
!     Limit output ranget'
!                                                                       
      ELSEIF (str_comp (befehl, 'range', 2, lbef, 5) ) THEN 
         CALL output_range(zeile, lp)
!                                                                       
!     Reset output 'reset'
!                                                                       
      ELSEIF (str_comp (befehl, 'reset', 2, lbef, 5) ) THEN 
         CALL output_reset
!                                                                       
!     write output file 'run'                                           
!                                                                       
      ELSEIF (str_comp (befehl, 'run ', 2, lbef, 4) ) THEN 
         IF(four_was_run) THEN    ! A fourier has been calculated do output

            CALL chem_elem(.false.)
            CALL set_output (linverse) 
            IF(value==val_3DPDF) THEN
               if(.not. allocated(dsi3d)) then
                  ier_num = -180
                  ier_typ = ER_APPL
                  ier_msg(1) = 'No Fourier has been calculated '
                  ier_msg(2) = 'or 3D-PDF has been written since last calculation'
                  return
               endif
               dsi = dsi3d
               dsi3d = 0.0D0
               if(opara(O_SHARP)(1:5)=='sharp') then
                  do k=1,num(3)
                     do j=1,num(2)
                        do i=1,num(1)
!                          l =  (i - 1)*out_inc(3)*out_inc(2) +        &
!                               (j - 1)*out_inc(3)            + k
                           dsi(i,j,k) = dsi(i,j,k)/qval(i,j,k,10,  i, j, laver)
                           dsi3d(i,j,k) =          qval(i,j,k,10,  i, j, laver)
                        enddo
                     enddo
                  enddo
               endif
               CALL out_prep_3dpdf(laver, l_val_limited, dsmax)
               deallocate(dsi3d)
            elseif(value==val_3DBETA) then
               ii = 0
                        do i=1,num(1)
                     do j=1,num(2)
                  do k=1,num(3)
               ii = ii + 1
               dsi(i,j,k) = real(acsf(i,j,k) * conjg(acsf(i,j,k)), kind=PREC_DP)
                        enddo
                     enddo
                  enddo
               CALL out_prep_3dpdf(laver, l_val_limited, dsmax)
            ENDIF
            IF (ityp.eq.0) THEN 
               CALL do_output (value, laver, valmax) 
            ELSEIF (ityp.eq.1) THEN 
               CALL do_post (value, laver) 
            ELSEIF (ityp.eq.2) THEN 
               CALL do_pgm (value, laver) 
            ELSEIF (ityp.eq.3) THEN 
               CALL do_output (value, laver, valmax) 
            ELSEIF (ityp.eq.4) THEN 
               CALL do_ppm (value, laver) 
            ELSEIF (ityp.eq.5) THEN 
if(pow_lperiod .and. value==VAL_PDF) then
   call pow_pdf_hist_prep_period(.false.)   ! In case of periodic boundary conditions 
endif
               if(cpow_form == 'tof') then
                  call powder_out_tof(value, .false.)
               else
                  CALL powder_out (value, .false.)
               endif
            ELSEIF (ityp.eq.6) THEN 
               CALL do_output (value, laver, valmax) 
            ELSEIF (ityp.eq.7) THEN 
               CALL do_output (value, laver, valmax) 
            ELSEIF (ityp.eq.8) THEN 
               CALL do_output (value, laver, valmax) 
            ELSEIF (ityp.eq.9) THEN 
               CALL do_output (value, laver, valmax) 
            ELSEIF (ityp.eq.10) THEN 
               CALL nexus_write (value, laver) 
            ELSEIF (ityp.eq.11) THEN
               CALL vtk_write ()
            ELSEIF (ityp.eq.12) THEN
               if(allocated(qvalues)) deallocate(qvalues)
               allocate(qvalues(out_inc(1), out_inc(2), out_inc(3)))
               l = 0                                                       ! Copy proper "value"
               DO i = 1, out_inc(1)
                  DO j = 1, out_inc(2)
                     DO k = 1, out_inc(3)
                        l = l + 1
                        qvalues(i,j,k) = qval(i,j,k, value, i, j, laver)
                     ENDDO
                  ENDDO
               ENDDO
               call mrc_write(outfile,                                                   &
                    out_inc, 2, nxyzstart, cr_ar, cr_wrez, extr_abs, extr_ord,   &
                    extr_top, qvalues                                            &
                    )
            ELSEIF (ityp.eq.13) THEN                           ! HDF5 output
               if(allocated(qvalues)) deallocate(qvalues)
               allocate(qvalues(out_inc(1), out_inc(2), out_inc(3)))
               l = 0                                                       ! Copy proper "value"
               DO i = 1, out_inc(1)
                  DO j = 1, out_inc(2)
                     DO k = 1, out_inc(3)
                        l = l + 1
                        qvalues(i,j,k) = qval(i,j,k, value, i, j, laver)
                     ENDDO
                  ENDDO
               ENDDO
               qvalues = qvalues*out_scale
               CALL get_params (zeile, ianz, cpara, lpara, maxp, lp) 
               IF (ier_num.eq.0) THEN 
               CALL get_optional(ianz, MAXP, cpara, lpara, NOPTIONAL,  ncalc, &
                                 oname, loname, opara, lopara, lpresent, owerte)
               if(lpresent(O_TRANS) .and.                                       &
                  str_comp(opara(O_TRANS), 'yes', 3, lopara(O_trans),3)) then

                  call lib_trans_menu(-1, value, laver, outfile, out_inc, out_eck, out_vi,            &
                       cr_a0, cr_win, qvalues,VAL_PDF, VAL_3DPDF,       &
                       new_outfile, new_inc, new_eck, new_vi, new_qvalues)
!
               else
                  CALL gen_hdf5_write (value, laver, outfile, out_inc, out_eck, out_vi, &
                          out_extr_abs, out_extr_ord, out_extr_top,                     &
                          cr_a0, cr_win, qvalues,val_pdf, val_3Dpdf, valmax,            &
                          ier_num, ier_typ, ER_IO, ER_APPL)
               endif
               endif
               deallocate(qvalues)
!              CALL hdf5_write (value, laver, outfile, out_inc, out_eck, out_vi, &
!                      cr_a0, cr_win, qval,val_pdf, val_3Dpdf, valmax,           &
!                      ier_num, ier_typ, ER_IO, ER_APPL)
            ELSEIF (ityp.eq.14) THEN
               CALL xplor_write (value, laver)
            ELSEIF (ityp.eq.15) THEN
               CALL   grd_write (value, laver, cpatt, spatt, dpatt)
            ELSE 
               ier_num = - 9 
               ier_typ = ER_APPL 
            ENDIF 
            if(allocated(rpdf)) deallocate(rpdf)  !3DPDF is written, array no longer needed
         ELSE 
            ier_num = -118
            ier_typ = ER_APPL 
            ier_msg(1) = 'You need to calculate a Fourier / Patterson /'
            ier_msg(2) = 'Inverse Fourier / Powder / Fourier via Stack'
            ier_msg(3) = 'first, before an output can be written'
         ENDIF 
!                                                                       
!     Show current settings for output 'show'                           
!                                                                       
      ELSEIF (str_comp (befehl, 'show', 2, lbef, 4) ) THEN 
         WRITE (output_io, 3000) outfile 
         IF (ityp.lt.0.or.13.lt.ityp) THEN 
            WRITE (output_io, * ) 'ityp undefiniert ', ityp 
         ELSEIF (ityp.eq.5) THEN 
            WRITE (output_io, 3130) cgraphik (ityp), cpow_form 
         ELSE 
            WRITE (output_io, 3100) cgraphik (ityp) 
            IF (laver) THEN 
               WRITE (output_io, 3110) '<'//cvalue (value) //'>' 
            ELSE 
               WRITE (output_io, 3110) cvalue (value) 
            ENDIF 
         ENDIF 
         IF(value==val_3dpdf) THEN
            IF(l_val_limited) THEN
              WRITE(output_io, '(a,f4.1,a)') ' Rec. space limited to d-star <=: ',dsmax, 'A^-1'
            ENDIF
         ENDIF
         WRITE (output_io, 3060) braggmin, braggmax, diffumin,    &
         diffumax, diffuave, diffusig                             
         WRITE (output_io, 3080) 100.0 * ps_high, zmax 
         WRITE (output_io, 3090) 100.0 * ps_low, zmin 
!                                                                       
!-------Operating System Kommandos 'syst'                               
!                                                                       
      ELSEIF (str_comp (befehl, 'system', 2, lbef, 5) ) THEN 
         IF (zeile.ne.' ') THEN 
            CALL do_operating (zeile (1:lp), lp) 
         ELSE 
            ier_num = - 6 
            ier_typ = ER_COMM 
         ENDIF 
!                                                                       
!     Set threshold for intensity written to bitmaps 'thresh'           
!                                                                       
      ELSEIF (str_comp (befehl, 'threshold', 2, lbef, 9) ) THEN 
         CALL get_params (zeile, ianz, cpara, lpara, maxp, lp) 
         IF (ier_num.eq.0) THEN 
            IF (ianz.eq.2) THEN 
               IF (str_comp (cpara (1) , 'high', 1, lpara (1) , 4)&
               ) THEN                                             
                  CALL del_params (1, ianz, cpara, lpara, maxp) 
                  CALL ber_params (ianz, cpara, lpara, werte,     &
                  maxp)                                           
                  IF (ier_num.eq.0) THEN 
                     ps_high = werte (1) * 0.01 
                     zmax = diffumax * ps_high 
                  ENDIF 
               ELSEIF (str_comp (cpara (1) , 'low', 1, lpara (1) ,&
               3) ) THEN                                          
                  CALL del_params (1, ianz, cpara, lpara, maxp) 
                  CALL ber_params (ianz, cpara, lpara, werte,     &
                  maxp)                                           
                  IF (ier_num.eq.0) THEN 
                     ps_low = werte (1) * 0.01 
                     zmin = diffumax * ps_low 
                  ENDIF 
               ELSEIF (str_comp (cpara (1) , 'sigma', 1, lpara (1)&
               , 5) ) THEN                                        
                  CALL del_params (1, ianz, cpara, lpara, maxp) 
                  CALL ber_params (ianz, cpara, lpara, werte,     &
                  maxp)                                           
                  IF (ier_num.eq.0) THEN 
                     zmin = max(diffumin, diffuave-REAL(werte(1), kind=PREC_DP) * diffusig)
                     zmax = min(diffumax, diffuave+REAL(werte(1), kind=PREC_DP) * diffusig)
                     IF (diffumax.ne.0) THEN 
                        ps_high = zmax / diffumax 
                        ps_low = zmin / diffumax 
                     ELSE 
                        ps_high = 0.0 
                        ps_low = 0.0 
                     ENDIF 
                  ENDIF 
               ELSEIF(str_comp (cpara (1) , 'zmax', 3, lpara (1) &
               , 4) ) THEN                                        
                  CALL del_params (1, ianz, cpara, lpara, maxp) 
                  CALL ber_params (ianz, cpara, lpara, werte,     &
                  maxp)                                           
                  IF (ier_num.eq.0) THEN 
                     zmax = werte (1) 
                     IF (diffumax.ne.0) THEN 
                        ps_high = zmax / diffumax 
                     ELSE 
                        ps_high = 0.0 
                     ENDIF 
                  ENDIF 
               ELSEIF (str_comp (cpara (1) , 'zmin', 3, lpara (1) &
                  , 4) ) THEN                                        
                  CALL del_params (1, ianz, cpara, lpara, maxp) 
                  CALL ber_params (ianz, cpara, lpara, werte,     &
                  maxp)                                           
                  IF (ier_num.eq.0) THEN 
                     zmin = werte (1) 
                  IF (diffumax.ne.0) THEN 
                     ps_low = zmin / diffumax 
                  ELSE 
                     ps_low = 0.0 
                  ENDIF 
               ENDIF 
            ELSE 
               ier_num = - 11 
               ier_typ = ER_APPL 
            ENDIF 
         ELSE 
            ier_num = - 6 
            ier_typ = ER_COMM 
         ENDIF 
      ENDIF 
!                                                                       
!     Define output value 'value'                                       
!                                                                       
   ELSEIF (str_comp (befehl, 'value', 1, lbef, 5) ) THEN 
      CALL get_params (zeile, ianz, cpara, lpara, maxp, lp) 
      IF (ier_num.eq.0) THEN 
         CALL get_optional(ianz, MAXP, cpara, lpara, NOPTIONAL,  ncalc, &
                           oname, loname, opara, lopara, lpresent, owerte)
!------ ----Check if we want the average values <F> ?                   
         IF (cpara (1) (1:1) .eq.'<') THEN 
            ix = 2 
            laver = .true. 
         ELSE 
            ix = 1 
            laver = .false. 
         ENDIF 
!     ----Calculate intensity 'intensity'                               
         IF (cpara (1) (ix:ix + 1) .eq.'in') THEN 
            value = val_inten
!     ----Calculate amplitude 'amplitude'                               
         ELSEIF (cpara (1) (ix:ix) .eq.'a') THEN 
            value = val_ampli
!     ----Calculate phase 'phase'                                       
         ELSEIF (cpara (1) (ix:ix) .eq.'p') THEN 
            IF (ianz.eq.1) THEN 
               value = val_phase
            ELSEIF (ianz.eq.2.and.cpara (2) (1:1) .eq.'r')     &
            THEN                                               
               value = val_ranph
            ENDIF 
!     ----Calculate real part 'real'                                    
         ELSEIF (cpara (1) (ix:ix) .eq.'r') THEN 
            value = val_real
!     ----Calculate imaginary part 'imaginary'                          
         ELSEIF (cpara (1) (ix:ix + 1) .eq.'im') THEN 
            value = val_imag
!     ----Calculate I(Q)           'I(Q) =Inte/N'                       
         ELSEIF (cpara (1) (ix:ix + 3) .eq.'I(Q)') THEN 
            value = val_iq
!     ----Calculate S(Q)           'S(Q)     '                          
         ELSEIF (cpara (1) (ix:ix + 3) .eq.'S(Q)') THEN 
            value = val_sq
!     ----Calculate F(Q)=Q(S(Q)-1) 'F(Q)     '                          
         ELSEIF (cpara (1) (ix:ix + 3) .eq.'F(Q)') THEN 
            value = val_fq 
         ELSEIF (cpara (1) (ix:ix + 6) == 'f2averb') THEN
            value = val_f2averb
         ELSEIF (cpara (1) (ix:ix + 6) == 'faver2b') THEN
            value = val_faver2b
         ELSEIF (cpara (1) (ix:ix + 5) == 'f2aver') THEN
            value = val_f2aver
         ELSEIF (cpara (1) (ix:ix + 5) == 'faver2') THEN
            value = val_faver2
         ELSEIF (cpara (1) (ix:ix + 5) == 'faverb') THEN
            value = val_faverb
         ELSEIF (cpara (1) (ix:ix + 4) == 'faver') THEN
            value = val_faver
!     ----Calculate S(Q)           'N(Q) = S(Q) without thermal part    '                          
         ELSEIF (cpara (1) (ix:ix + 2) == 'PDF' ) THEN 
            value = val_pdf
         ELSEIF (cpara (1) (ix:ix + 4) == '3DPDF' ) THEN 
            CALL value_optional(lpresent, O_DSMAX, O_QMAX, O_HKLMAX, NOPTIONAL, opara, &
                                lopara, werte, l_val_limited, dsmax)
            value = val_3DPDF
            if(opara(O_SHARP)(1:6)=='normal') then
               pdf3d_mode = PDF3D_NORMAL
            elseif(opara(O_SHARP)(1:5)=='sharp') then
               pdf3d_mode = PDF3D_SHARP
            endif
         ELSEIF (cpara (1) (ix:ix + 5) == '3DBETA' ) THEN 
            CALL value_optional(lpresent, O_DSMAX, O_QMAX, O_HKLMAX, NOPTIONAL, opara, &
                                lopara, werte, l_val_limited, dsmax)
            value = val_3DBETA
            if(opara(O_SHARP)(1:6)=='normal') then
               pdf3d_mode = PDF3D_NORMAL
            elseif(opara(O_SHARP)(1:5)=='sharp') then
               pdf3d_mode = PDF3D_SHARP
            endif
         ELSE 
            ier_num = - 6 
            ier_typ = ER_COMM 
            value = 0 
         ENDIF 
!------ ----check lots and allowed output                               
!        IF (nlots.ne.1.and..NOT.(value==1 .OR. value==val_3dpdf).and..not.laver) THEN 
         IF (nlots.ne.1 .AND.  .NOT.laver .AND.   &
             (value==val_ampli .OR. value==val_phase .OR. value==val_ranph .OR. &
              value==val_real  .OR. value==val_imag)    ) THEN
            ier_num = - 60 
            ier_typ = ER_APPL 
            value = 0 
         ENDIF 
      ENDIF 
!                                                                       
!------  -waiting for user input                                        
!                                                                       
   ELSEIF (str_comp (befehl, 'wait', 3, lbef, 4) ) THEN 
      CALL do_input (zeile, lp) 
!                                                                       
!------ no valid subcommand found                                       
!                                                                       
   ELSE 
      ier_num = - 8 
      ier_typ = ER_COMM 
   ENDIF 
ENDIF 
ENDIF  main_if
!
IF (ier_num.ne.0) THEN 
   CALL errlist 
   IF (ier_sta.ne.ER_S_LIVE) THEN 
      IF (lmakro .OR. lmakro_error) THEN  ! Error within macro or termination errror
         IF(sprompt /= prompt ) THEN
            ier_num = -10
            ier_typ = ER_COMM
            ier_msg(1) = ' Error occured in output menu'
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
ENDIF 
IF(linteractive .or. lmakro) THEN
   GOTO 10
ELSE
   RETURN
ENDIF
 9999 CONTINUE 
prompt = orig_prompt
!                                                                       
 1015 FORMAT ( /1x,'Z-MIN = ',G20.6,/,1x,'Z-MAX = ',G20.6,//            &
     &                     1x,'Give new values zmin, zmax    : ')     
 3000 FORMAT( ' Output file                    : ',a) 
 3060 FORMAT(/' Bragg minimum                  : ',g13.6/                 &
     &        ' Bragg maximum                  : ',g13.6/                 &
     &        ' Diffuse minimum                : ',g13.6/                 &
     &        ' Diffuse maximum                : ',g13.6/                 &
     &        ' Diffuse average intensity      : ',g13.6/                 &
     &        ' Diffuse intensity sigma        : ',g13.6)                 
 3080 FORMAT(/' Maximum value for BITMAP'/                              &
     &        ' in % of highest diffuse value'/                         &
     &        ' and absolute                   : ',2x,f9.4,2x,g13.6)      
 3090 FORMAT( ' Minimum value for BITMAP'/                              &
     &        ' in % of highest diffuse value'/                         &
     &        ' and absolute                   : ',2x,f9.4,2x,g13.6)      
 3100 FORMAT( ' Graphicsformat                 : ',A) 
 3130 FORMAT( ' Graphicsformat                 : ',A,A) 
 3110 FORMAT( ' Output value                   : ',A) 
END SUBROUTINE do_niplps                      
!
!*****7*************************************************************************
!
SUBROUTINE output_range(zeile, length)
!-
!  Interpret range command, limit output range for 3D-maps
!+
!
USE output_mod
!
USE errlist_mod
USE get_params_mod
USE precision_mod
USE str_comp_mod
USE take_param_mod
!
IMPLICIT NONE
!
CHARACTER(LEN=*), INTENT(INOUT) :: zeile
INTEGER         , INTENT(INOUT) :: length
! 
INTEGER, PARAMETER :: MAXP = 3
!
CHARACTER(LEN=PREC_STRING), DIMENSION(MAXP) :: cpara
INTEGER                   , DIMENSION(MAXP) :: lpara
INTEGER                                     :: ianz
REAL(KIND=PREC_DP)        , DIMENSION(MAXP) :: werte
!
CHARACTER(LEN=PREC_STRING)                  :: ccpara
INTEGER                                     :: llpara
!
INTEGER, PARAMETER :: NOPTIONAL = 4
INTEGER, PARAMETER :: O_CENTER  = 1                  ! Center of map
INTEGER, PARAMETER :: O_PIXEL   = 2                  ! Pixels at plus minus
INTEGER, PARAMETER :: O_QUAD    = 3                  ! Select a quadrant
!INTEGER, PARAMETER :: O_HKLMAX  = 4                  ! Maximum hkl vector
CHARACTER(LEN=   6), DIMENSION(NOPTIONAL) :: oname   !Optional parameter names
CHARACTER(LEN=PREC_STRING), DIMENSION(NOPTIONAL) :: opara   !Optional parameter strings returned
INTEGER            , DIMENSION(NOPTIONAL) :: loname  !Lenght opt. para name
INTEGER            , DIMENSION(NOPTIONAL) :: lopara  !Lenght opt. para name returned
LOGICAL            , DIMENSION(NOPTIONAL) :: lpresent!opt. para is present
REAL(KIND=PREC_DP) , DIMENSION(NOPTIONAL) :: owerte   ! Calculated values
INTEGER, PARAMETER                        :: ncalc = 0 ! Number of values to calculate
!                                                                       
DATA oname  / 'center', 'pixel ', 'quad  ', 'hklmax'/
DATA loname /  6      ,  5      ,  4      ,  6      /
DATA opara  / 'middle', '0.00  ', 'rrr   ', '0.00  '/
DATA owerte /  0.0000 ,  0.000  ,  0.000  ,  0.000  /
!
CALL get_params (zeile, ianz, cpara, lpara, maxp, length) 
IF (ier_num /= 0) RETURN
!
IF(str_comp(cpara(1), 'full', 2, lpara(1), 4) ) THEN 
   out_lrange  = 0                                 ! Output range is not limited
!  out_lcenter = .FALSE.
!  out_lrange  = .FALSE.
   RETURN
ENDIF
!
CALL get_optional(ianz, MAXP, cpara, lpara, NOPTIONAL,  ncalc, &
                  oname, loname, opara, lopara, lpresent, owerte)
!
IF(lpresent(O_CENTER)) THEN                         ! User provided "center:[cx,cy,cz]" parameter
   IF(opara(O_CENTER)=='middle') THEN
      out_lrange  = 1                            ! Limit range
      out_lcenter = 0                            ! Center is at midpoint
   ELSE
      ccpara = opara(O_CENTER)
      llpara = lopara(O_CENTER)
      CALL get_optional_multi(MAXP, ccpara, llpara, werte, ianz)
      IF (ier_num /= 0) RETURN
      out_center(1) = NINT(werte(1))
      out_center(2) = NINT(werte(2))
      out_center(3) = NINT(werte(3))
      out_lrange  = 1                            ! Range is limited
      out_lcenter = 1                            ! Center is at user values
   ENDIF
ELSE
   out_lcenter = 0                               ! Center is at midpoint
ENDIF
!
IF(lpresent(O_PIXEL)) THEN                         ! User provided "pixel:[px,py,pz]" parameter
   ccpara = opara(O_PIXEL)
   llpara = lopara(O_PIXEL)
   CALL get_optional_multi(MAXP, ccpara, llpara, werte, ianz)
   IF (ier_num /= 0) RETURN
   out_pixel(1) = NINT(werte(1))
   out_pixel(2) = NINT(werte(2))
   out_pixel(3) = NINT(werte(3))
   out_lpixel  = 1                                 !User provided pixels
   out_lrange  = 1
ELSE
   out_pixel = 0                                   ! All pixels, may be at off center
ENDIF
!
IF(lpresent(O_QUAD)) THEN                          ! User provided "quad:rll" parameter
   out_quad = opara(O_QUAD)(1:3)
   out_lrange = 2
ELSE
   out_quad =  ' '
ENDIF
!
END SUBROUTINE output_range
!
!*****7*************************************************************************
!
SUBROUTINE form_optional(lpresent, ientry, MAXW, opara, &
                       lopara, werte)
!
USE ber_params_mod
USE precision_mod
USE str_comp_mod
!
IMPLICIT NONE
!
LOGICAL                            , INTENT(IN)    :: lpresent ! Optional parameter was present 
INTEGER                            , INTENT(IN)    :: ientry   ! Entry number
INTEGER                            , INTENT(IN)    :: MAXW     ! Dimension of werte
CHARACTER(LEN=*)  , DIMENSION(MAXW), INTENT(INOUT) :: opara    ! The string with optional values
INTEGER           , DIMENSION(MAXW), INTENT(INOUT) :: lopara   ! length of string
REAL(KIND=PREC_DP), DIMENSION(MAXW), INTENT(OUT)   :: werte    ! Numerical values
!
INTEGER :: ianz
!
ianz = ientry
IF(lpresent) THEN
   IF(str_comp(opara(ientry), 'auto', 4, lopara(ientry), 4) ) THEN
      werte(ientry) = -2.0D0
   elseif(str_comp(opara(ientry), 'bragg', 5, lopara(ientry), 5) ) THEN
      werte(ientry) = -2.0D0
   ELSEIF(str_comp(opara(ientry), 'data', 4, lopara(ientry), 4) ) THEN
      werte(ientry) = -1.0D0
   ELSE
      CALL ber_params (ianz, opara, lopara, werte, MAXW)
   ENDIF
ELSE
   werte(ientry) = -1.0D0
ENDIF
!
END SUBROUTINE form_optional
!
!*****7*************************************************************************
!
SUBROUTINE value_optional(lpresent, O_DSMAX, O_QMAX, O_HKLMAX, MAXW, opara, &
                          lopara, owerte, l_val_limited, dsmax)
!
!-
!  Sets limits for d-star to be applied prior to 3D-PDF output. all points in
!  reciprocal space outside a spher with dsmax radias are set to intensity zero
!+
!
USE metric_mod
!
USE errlist_mod
USE ber_params_mod
USE precision_mod
!USE str_comp_mod
USE take_param_mod
USE wink_mod
!
IMPLICIT NONE
!
INTEGER                            , INTENT(IN)    :: MAXW     ! Dimension of werte
LOGICAL           , DIMENSION(MAXW), INTENT(IN)    :: lpresent ! Optional parameter was present 
INTEGER                            , INTENT(IN)    :: O_DSMAX  ! Entry number
INTEGER                            , INTENT(IN)    :: O_QMAX   ! Entry number
INTEGER                            , INTENT(IN)    :: O_HKLMAX ! Entry number
CHARACTER(LEN=*)  , DIMENSION(MAXW), INTENT(INOUT) :: opara    ! The string with optional values
INTEGER           , DIMENSION(MAXW), INTENT(INOUT) :: lopara   ! length of string
REAL(KIND=PREC_DP), DIMENSION(MAXW), INTENT(OUT)   :: owerte   ! Numerical values
LOGICAL                            , INTENT(OUT)   :: l_val_limited ! Limitation is applied T/F
REAL(KIND=PREC_DP)                 , INTENT(OUT)   :: dsmax    ! Final maximum dstar
!
INTEGER, PARAMETER :: MAXP = 3
REAL(KIND=PREC_DP), DIMENSION(3) :: NULLV = (/0.00D0, 0.00D0, 0.00D0/)   ! Null vector
!
CHARACTER(LEN=1024)                 :: cpara
INTEGER                             :: lpara
REAL(KIND=PREC_DP), DIMENSION(MAXW) :: werte   ! Numerical values
REAL(KIND=PREC_DP), DIMENSION(MAXW) :: u       ! Numerical values
INTEGER :: ianz, i
!
werte         =  0.0D0
dsmax         = 0.0
l_val_limited = .FALSE.
!
IF(.NOT.(lpresent(O_DSMAX) .OR. lpresent(O_QMAX) .OR. lpresent(O_HKLMAX))) THEN
   l_val_limited = .FALSE.
   dsmax = 0.0
   RETURN
ELSEIF(lpresent(O_DSMAX) .AND. .NOT. (lpresent(O_QMAX) .OR. lpresent(O_HKLMAX))) THEN
   DO i=1,MAXW
      IF(i/=O_DSMAX) THEN
         opara(i)  = '0'
         lopara(i) = 1
      ENDIF
   ENDDO
   ianz = O_DSMAX
   CALL ber_params (ianz, opara, lopara, owerte, MAXW)
   IF(ier_num /= 0) RETURN
   dsmax         = owerte(O_DSMAX)
   l_val_limited = .TRUE.
ELSEIF(lpresent(O_QMAX) .AND. .NOT. (lpresent(O_DSMAX) .OR. lpresent(O_HKLMAX))) THEN
   DO i=1,MAXW
      IF(i/=O_QMAX) THEN
         opara(i)  = '0'
         lopara(i) = 1
      ENDIF
   ENDDO
   ianz = O_QMAX
   CALL ber_params (ianz, opara, lopara, owerte, MAXW)
   IF(ier_num /= 0) RETURN
   dsmax         = owerte(O_QMAX)/zpi
   l_val_limited = .TRUE.
ELSEIF(lpresent(O_HKLMAX) .AND. .NOT. (lpresent(O_DSMAX) .OR. lpresent(O_QMAX))) THEN
   cpara = opara (O_HKLMAX)
   lpara = lopara(O_HKLMAX)
   CALL get_optional_multi(MAXP, cpara, lpara, werte, ianz)
   u = werte
   dsmax = do_blen (.FALSE., u, NULLV)
   l_val_limited = .TRUE.
ELSE
   owerte        =  0.0D0
   dsmax         = 0.0
   l_val_limited = .FALSE.
   ier_num = -174
   ier_typ = ER_APPL
   ier_msg(1) = 'Only one of ''dsmax'', ''qmax'', ''hklmax'' allowed'
ENDIF
!
END SUBROUTINE value_optional
!
!*****7*****************************************************************
!
SUBROUTINE do_post (value, laver) 
!-                                                                      
!     Writes a POSTSCRIPT file                                          
!+                                                                      
USE discus_config_mod 
USE diffuse_mod 
USE output_mod 
USE qval_mod
USE envir_mod 
USE errlist_mod 
use precision_mod
USE support_mod
!
IMPLICIT none 
!
integer, intent(in) :: value
logical, intent(in) :: laver
!
ier_num = -6
ier_typ = ER_COMM
ier_msg(1) = "POSTSCRIPT output has been disabled. Please use KUPLOT"
ier_msg(2) = "to display the data and write the image as POSTSCRIPT format"
!                                                                       
!                                                                       
!INTEGER maxcol 
!PARAMETER (maxcol = 256) 
!!                                                                       
!CHARACTER(len=6) :: cout (maxqxy) 
!CHARACTER(len=6) :: cfarb (maxcol) 
!INTEGER :: i, ix, iy, iqqq, k, value 
!LOGICAL :: lread, laver 
!REAL(kind=PREC_DP) :: qqq 
!!                                                                       
!!                                                                       
!!     Check whether data are 2-dimensional                              
!!                                                                       
!      IF (.not. (out_inc (1) .gt.1.and.out_inc (2) .gt.1) ) THEN 
!         ier_num = - 50 
!         ier_typ = ER_APPL 
!         RETURN 
!      ENDIF 
!!                                                                       
!!-------Farbtabelle einlesen                                            
!!                                                                       
!      lread = .true. 
!      CALL oeffne (2, colorfile, 'old') 
!      IF (ier_num.ne.0) return 
!      DO i = 1, 255 
!      READ (2, 100, end = 20) cfarb (i) 
!  100 FORMAT      (1x,a6) 
!      ENDDO 
!   20 CONTINUE 
!      CLOSE (2) 
!      cfarb (256) = 'ffffff' 
!!                                                                       
!      lread = .false. 
!      CALL oeffne (2, outfile, 'unknown') 
!      IF (ier_num.ne.0) return 
!!                                                                       
!      WRITE (2, 1111) '%!PS-Adobe-2.0' 
!      WRITE (2, 1111) '%%Creator: DISCUS, Version 3.0' 
!      WRITE (2, 1111) '50  150 translate' 
!      WRITE (2, 1111) '288 288 scale' 
!      WRITE (2, 2000) nint (3.0 * out_inc (1) ), out_inc (1), out_inc ( &
!      2), i, 8, out_inc (1), 0, 0, out_inc (2), 0, 0                    
!!                                                                       
!      DO iy = 1, out_inc (2) 
!      DO ix = 1, out_inc (1) 
!      k = (ix - 1) * out_inc (2) + iy 
!      qqq = qval (ix,iy,1, value, ix, iy, laver) 
!      IF (qqq.lt.zmin) THEN 
!         qqq = zmin 
!      ELSEIF (qqq.gt.zmax) THEN 
!         qqq = zmax 
!      ENDIF 
!      iqqq = nint ( (maxcol - 2) * (qqq - zmin) / (zmax - zmin) )       &
!      + 1                                                               
!      WRITE (cout (ix), 1111) cfarb (iqqq) 
!      DO i = 1, 6 
!      IF (cout (ix) (i:i) .eq.' '.or.cout (ix) (i:i) .eq.' ') cout (ix) &
!      (i:i) = '0'                                                       
!      ENDDO 
!      ENDDO 
!      WRITE (2, 5000) (cout (ix), ix = 1, out_inc (1) ) 
!      ENDDO 
!      WRITE (2, 1111) 'showpage' 
!!                                                                       
!      CLOSE (2) 
!!                                                                       
! 1111 FORMAT (a) 
! 2000 FORMAT ('/DataString ',I4,' string def'/3(I3,1X),                 &
!     &        ' [ ',6(I3,1X),']'/'{'/                                   &
!     &        '  currentfile DataString readhexstring pop'/             &
!     &        ' }  false 3 colorimage')                                 
! 5000 FORMAT (10A6) 
!                                                                       
      END SUBROUTINE do_post                        
!
!*****7*****************************************************************
!
SUBROUTINE do_pgm (value, laver) 
!-                                                                      
!     Writes the data in a PGM format                                   
!-                                                                      
USE discus_config_mod 
USE diffuse_mod 
USE output_mod 
USE qval_mod
USE errlist_mod 
use precision_mod
USE support_mod
IMPLICIT none 
!                                                                       
!                                                                       
integer, intent(in) :: value
logical, intent(in) :: laver
!
ier_num = -6
ier_typ = ER_COMM
ier_msg(1) = "PGM output has been disabled. Please use KUPLOT"
ier_msg(2) = "to display the data and write the image as PGM format"
!                                                                       
!INTEGER maxcol 
!PARAMETER (maxcol = 255) 
!!                                                                       
!INTEGER :: iqqq (maxqxy), ncol, ix, iy, k, value 
!LOGICAL :: lread, laver 
!REAL(kind=PREC_DP) ::  qqq 
!!                                                                       
!!                                                                       
!!     Check whether data are 2-dimensional                              
!!                                                                       
!      IF (.not. (out_inc (1) .gt.1.and.out_inc (2) .gt.1) ) THEN 
!         ier_num = - 50 
!         ier_typ = ER_APPL 
!         RETURN 
!      ENDIF 
!!                                                                       
!      lread = .false. 
!      ncol = maxcol 
!!                                                                       
!      CALL oeffne (2, outfile, 'unknown') 
!      IF (ier_num.ne.0) return 
!!                                                                       
!      WRITE (2, 1111) 'P2' 
!      WRITE (2, 2000) out_inc (1), out_inc (2), ncol 
!!                                                                       
!      DO iy = out_inc (2), 1, - 1 
!      DO ix = 1, out_inc (1) 
!      k = (ix - 1) * out_inc (2) + iy 
!      qqq = qval (ix,iy,1, value, ix, iy, laver) 
!      IF (qqq.lt.zmin) THEN 
!         qqq = zmin 
!      ELSEIF (qqq.gt.zmax) THEN 
!         qqq = zmax 
!      ENDIF 
!      iqqq(ix) = nint(REAL(ncol - 1, kind=PREC_DP) * (qqq - zmin) / (zmax - zmin))
!      ENDDO 
!      WRITE (2, 5000) (iqqq (ix), ix = 1, out_inc (1) ) 
!      ENDDO 
!!                                                                       
!      CLOSE (2) 
!!                                                                       
! 1111 FORMAT     (a) 
! 2000 FORMAT    (2(1x,i4)/1x,i8) 
! 5000 FORMAT     (7i8) 
!!                                                                       
      END SUBROUTINE do_pgm                         
!
!*****7*****************************************************************
!
SUBROUTINE do_ppm (value, laver) 
!+                                                                      
!     Writes the data in a PGM format                                   
!-                                                                      
USE discus_config_mod 
USE diffuse_mod 
USE output_mod 
USE qval_mod
USE envir_mod 
USE errlist_mod 
USE support_mod
IMPLICIT none 
!                                                                       
integer, intent(in) :: value
logical, intent(in) :: laver
!
ier_num = -6
ier_typ = ER_COMM
ier_msg(1) = "PPM output has been disabled. Please use KUPLOT"
ier_msg(2) = "to display the data and write the image as PPM format"
!!                                                                       
!INTEGER maxcol 
!PARAMETER (maxcol = 255) 
!!                                                                       
!INTEGER :: iqqq (maxqxy), ncol, ix, iy, k, value 
!INTEGER :: icolor (maxcol, 3) 
!INTEGER :: i, j 
!LOGICAL :: lread, laver 
!REAL(kind=PREC_DP) :: qqq 
!!                                                                       
!      CHARACTER(6) cfarb (256) 
!!                                                                       
!!                                                                       
!!     Check whether data are 2-dimensional                              
!!                                                                       
!!      IF (.not. (out_inc (1) .gt.1.and.out_inc (2) .gt.1) ) THEN 
!         ier_num = - 50 
!         ier_typ = ER_APPL 
!!         RETURN 
!      ENDIF 
!!                                                                       
!!-------Farbtabelle einlesen                                            
!!                                                                       
!      CALL set_colorfeld (cfarb) 
!      lread = .true. 
!      CALL oeffne (2, colorfile, 'old') 
!      IF (ier_num.ne.0) return 
!      DO i = 1, 255 
!      READ (2, 100, end = 20) (icolor (i, j), j = 1, 3) 
!!  100 FORMAT      (1x,3z2) 
!      ENDDO 
!   20 CONTINUE 
!!      CLOSE (2) 
!!                                                                       
!      icolor (255, 1) = 255 
!!      icolor (255, 2) = 255 
!      icolor (255, 3) = 255 
!!                                                                       
!!      lread = .false. 
!      CALL oeffne (2, outfile, 'unknown') 
!!      IF (ier_num.ne.0) return 
!!                                                                       
!      ncol = maxcol 
!      WRITE (2, 1111) 'P3' 
!      WRITE (2, 2000) out_inc (1), out_inc (2), ncol 
!!                                                                       
!      DO iy = out_inc (2), 1, - 1 
!      DO ix = 1, out_inc (1) 
!      k = (ix - 1) * out_inc (2) + iy 
!      qqq = qval (ix,iy,1, value, ix, iy, laver) 
!      IF (qqq.lt.zmin) THEN 
!         qqq = zmin 
!      ELSEIF (qqq.gt.zmax) THEN 
!         qqq = zmax 
!      ENDIF 
!      iqqq (ix) = nint(REAL(ncol-1, kind=PREC_DP) * (qqq-zmin)/(zmax-zmin)) + 1
!      ENDDO 
!      WRITE(2, 5000) ((icolor(iqqq(ix), j), j=1, 3), ix=1, out_inc(1))
!      ENDDO 
!!                                                                       
!      CLOSE (2) 
!!                                                                       
! 1111 FORMAT     (a) 
! 2000 FORMAT    (1x,2i4/1x,i8) 
! 5000 FORMAT     (15i4) 
!                                                                       
END SUBROUTINE do_ppm                         
!
!*****7*****************************************************************
!
SUBROUTINE set_colorfeld (cfarb) 
!+                                                                      
!     This routine sets the pseudo color color map                      
!-                                                                      
use precision_mod
!
IMPLICIT NONE
!
INTEGER, PARAMETER ::maxcol = 256 
!                                                                       
CHARACTER (LEN=*),DIMENSION(maxcol) :: cfarb !(maxcol) 
CHARACTER(LEN=20),DIMENSION(maxcol) :: ccc   ! (256) 
INTEGER,          DIMENSION(3)      :: rgb (3) 
INTEGER                             :: i,ii,j,ifarb
REAL(kind=PREC_DP)                  :: rh, rp, rq, rt, rf
!                                                                       
      cfarb (1) = '000000' 
!                                                                       
      DO ifarb = 2, 256 
      rh = 0.1 + REAL(ifarb - 1) / 283.0 
      rh = 6.0 * rh 
      i = int (rh) 
      rf = rh - REAL(i) 
      rp = 0.0 
      rq = 1.0 - rf 
      rt = (1.0 - (1.0 - rf) ) 
!                                                                       
      IF (rt.gt.1.0) rt = 1.0 
      IF (rp.gt.1.0) rp = 1.0 
      IF (rq.gt.1.0) rq = 1.0 
!                                                                       
      IF (i.eq.0) THEN 
         rgb (1) = int (0.5 + 1. * 255.0) 
         rgb (2) = int (0.5 + rt * 255.0) 
         rgb (3) = int (0.5 + rp * 255.0) 
      ELSEIF (i.eq.1) THEN 
         rgb (1) = int (0.5 + rq * 255.0) 
         rgb (2) = int (0.5 + 1. * 255.0) 
         rgb (3) = int (0.5 + rp * 255.0) 
      ELSEIF (i.eq.2) THEN 
         rgb (1) = int (0.5 + rp * 255.0) 
         rgb (2) = int (0.5 + 1. * 255.0) 
         rgb (3) = int (0.5 + rt * 255.0) 
      ELSEIF (i.eq.3) THEN 
         rgb (1) = int (0.5 + rp * 255.0) 
         rgb (2) = int (0.5 + rq * 255.0) 
         rgb (3) = int (0.5 + 1. * 255.0) 
      ELSEIF (i.eq.4) THEN 
         rgb (1) = int (0.5 + rt * 255.0) 
         rgb (2) = int (0.5 + rp * 255.0) 
         rgb (3) = int (0.5 + 1. * 255.0) 
      ELSEIF (i.eq.5) THEN 
         rgb (1) = int (0.5 + 1. * 255.0) 
         rgb (2) = int (0.5 + rp * 255.0) 
         rgb (3) = int (0.5 + rq * 255.0) 
      ENDIF 
!                                                                       
      WRITE (ccc (ifarb), 1000) (rgb (j), j = 3, 1, - 1) 
      DO ii = 1, 6 
      IF (ccc (ifarb) (ii:ii) .eq.' ') ccc (ifarb) (ii:ii) = '0' 
      ENDDO 
!                                                                       
      WRITE (cfarb (ifarb), 1000) (rgb (j), j = 3, 1, - 1) 
      DO ii = 1, 6 
      IF (cfarb (ifarb) (ii:ii) .eq.' ') cfarb (ifarb) (ii:ii) = '0' 
      ENDDO 
      ENDDO 
      DO ifarb = 256, 2, - 1 
      WRITE (44, 2000) ccc (ifarb) 
      ENDDO 
      CLOSE (44) 
!                                                                       
 1000 FORMAT     (3(z2)) 
 2000 FORMAT    ('#',a6,'00') 
      END SUBROUTINE set_colorfeld                  
!
!*****7*****************************************************************
!
SUBROUTINE do_output (value, laver, valmax) 
!-                                                                      
!     Writes output in standard or GNUPLOT format                       
!+                                                                      
USE discus_config_mod 
USE crystal_mod 
USE diffuse_mod 
USE discus_nipl_header
USE discus_fft_mod
use discus_output_save_mod
USE fourier_sup
use metric_mod
USE output_mod 
USE qval_mod
!
USE envir_mod 
USE errlist_mod 
USE lib_length
use precision_mod
USE prompt_mod 
USE support_mod
use trig_degree_mod
!
IMPLICIT none 
!                                                                       
integer, intent(out) :: value       
logical, intent(in) :: laver
real(kind=PREC_DP), intent(in) :: valmax      
!real(kind=PREC_DP) :: valmax      
!                                                                       
INTEGER iff 
PARAMETER (iff = 2) 
real(kind=PREC_DP),  dimension(3), parameter :: NULLV =(/0.0_PREC_DP, 0.0_PREC_DP, 0.0_PREC_DP/)
!                                                                       
CHARACTER(LEN=2024) dummy_file
INTEGER HKLF4, LIST5, LIST9 , ASCII3D
PARAMETER (HKLF4 = 6, LIST5 = 7, LIST9 = 8, ASCII3D = 9) 
!                                                                       
INTEGER :: extr_ima, i, j, k, l, m
LOGICAL :: lread
REAL(kind=PREC_DP), dimension(3) ::  h (3) 
REAL(kind=PREC_DP) ::  sq, qq, out_fac 
!                                                                       
INTEGER shel_inc (3) 
INTEGER shel_value 
REAL   (kind=PREC_DP) ::  shel_eck (3, 4) 
REAL   (kind=PREC_DP) ::  shel_vi (3, 3) 
REAL   (kind=PREC_DP) ::  shel_000 
COMPLEX(kind=PREC_DP) ::  shel_csf
COMPLEX(kind=PREC_DP) ::  shel_acsf
REAL   (kind=PREC_DP) ::  shel_dsi
COMPLEX(kind=PREC_DP) ::  shel_tcsf
REAL   (kind=PREC_DP) ::  factor
!
INTEGER                            :: npkt1      ! Points in 1D files standard file format
INTEGER                            :: npkt2      ! Points in 2D files standard file format
INTEGER                            :: npkt3      ! Points in 3D files standard file format
INTEGER                            :: all_status ! Allocation status 
INTEGER                            :: is_dim     ! dimension of standard output file
INTEGER                            :: is_axis    ! Axis of standard output file
INTEGER, DIMENSION(1:3)            :: loop       ! Allows flexible loop index
INTEGER, DIMENSION(1:3)            :: out_index  ! Index that is written along axis 
REAL(kind=PREC_DP), DIMENSION(1:4)            :: ranges     ! xmin, xmax, ymin, ymax for NIPL files
REAL(kind=PREC_DP), DIMENSION(:), ALLOCATABLE :: xwrt  ! 'x' - values for standard 1D files
REAL(kind=PREC_DP), DIMENSION(:), ALLOCATABLE :: ywrt  ! 'x' - values for standard 1D files
REAL(kind=PREC_DP), DIMENSION(:,:), ALLOCATABLE :: zwrt  ! 'z' - values for standard 2D files
INTEGER :: nnew1, nnew2
!
logical :: four_log_user   ! Store user setting for Fourier log
CHARACTER (LEN=160), DIMENSION(:), ALLOCATABLE :: header_lines
INTEGER :: nheader
real(kind=PREC_DP) :: qqmax
real(kind=PREC_DP) :: ext_cor
!                                                                       
factor = 0.0
npkt3  = 1
m = nint(valmax)
out_fac = 1.0D0   ! Default to no scaling
qqmax = 0.0_PREC_DP
!                                                                       
!     If output type is shelx, calculate qval(000) for scaling          
!                                                                       
four_log_user = four_log
IF(m==-1 .and. (ityp.eq.HKLF4.or.ityp.eq.LIST5)) THEN    ! Scale with 000
!        DO i = 1, 3 
!        shel_inc (i) = inc (i) 
!        ENDDO 
!        DO i = 1, 3 
!        DO j = 1, 4 
!        shel_eck (i, j) = eck (i, j) 
!        ENDDO 
!        DO j = 1, 3 
!        shel_vi (i, j) = vi (i, j) 
!        ENDDO 
!        ENDDO 
   shel_inc = inc
   shel_eck = eck
   shel_vi  = vi 
   shel_tcsf = CMPLX(csf (1,1,1),KIND=KIND(0.0D0))
   shel_acsf = CMPLX(acsf(1,1,1),KIND=KIND(0.0D0))
   shel_dsi = dsi(1,1,1)
!        inc (1) = 1 
!        inc (2) = 1 
!        inc (3) = 1 
!        DO i = 1, 3 
!        DO j = 1, 4 
!        eck (i, j) = 0.0 
!        ENDDO 
!        DO j = 1, 3 
!        vi (i, j) = 0.0 
!        ENDDO 
!        ENDDO 
   inc = 1
   num = 1
   eck = 0.0D0
   vi  = 0.0D0
   IF (ityp.eq.HKLF4) THEN 
      value = 2 
   ELSEIF (ityp.eq.LIST5) THEN 
      value = 2 
   ENDIF 
   four_log = .false.    ! Turn Fourier log off
   CALL four_run 
   four_log = four_log_user  ! Turn Fourier log back to user setting
   shel_csf = CMPLX(csf (1,1,1), KIND=KIND(0.0D0))
   shel_000 = qval (1,1,1, value, 1, 1, laver) 
   qq = (qval(1,1,1, value, 1, 1, laver) / cr_icc(1)/cr_icc(2)/cr_icc(3))**2
   IF (ityp.eq.HKLF4) THEN 
      factor = max (int (log (qq) / log (10.0D0) ) - 3, 0) 
   ELSEIF (ityp.eq.LIST5) THEN 
      factor = max (int (log (qq) / log (10.0D0) ) - 2, 0) 
   ENDIF 
   out_fac = 10** ( - factor) 
   csf (1,1,1) = shel_tcsf
   acsf (1,1,1) = shel_acsf
   dsi (1,1,1) = shel_dsi
   inc = shel_inc
   num = shel_inc
   eck = shel_eck
   vi  = shel_vi 
!        DO i = 1, 3 
!        inc (i) = shel_inc (i) 
!        ENDDO 
!        DO i = 1, 3 
!        DO j = 1, 4 
!        eck (i, j) = shel_eck (i, j) 
!        ENDDO 
!        DO j = 1, 3 
!        vi (i, j) = shel_vi (i, j) 
!        ENDDO 
!        ENDDO 
elseif(m==-2 .and. (ityp.eq.HKLF4.or.ityp.eq.LIST5)) THEN    ! Scale with largest BRAG /= 000
   qqmax = 0.0D0
   if(ityp==HKLF4) then
!     shel_value = value 
      shel_value = 2 
   elseif(ityp==LIST5) then
      shel_value = 2 
   endif
   DO l = 1, out_inc (3) 
      DO j = 1, out_inc (2) 
         DO i = 1, out_inc (1) 
            DO k = 1, 3 
               h (k) = out_eck (k, 1) + out_vi (k, 1) * REAL(i - 1)   &
                                      + out_vi (k, 2) * REAL(j - 1)   &
                                      + out_vi (k, 3) * REAL(l - 1)
            ENDDO 
            IF( (INT(h(1)))**2 + (INT(h(2)))**2 + (INT(h(3)))**2 /= 0 ) THEN
!              k  = (i - 1) * out_inc (2) + j 
               k  = (i - 1) * out_inc (3) * out_inc (2) + (j-1) * out_inc (3) + l 
               qq = (qval(i,j,l, shel_value, i, j, laver) / cr_icc (1) / cr_icc (2) / cr_icc (3))**2
               qqmax = max(qqmax, qq)
            ENDIF
         ENDDO 
      ENDDO 
   ENDDO 
   if(ityp==HKLF4) then
      factor = max (int (log (qqmax) / log (10.0D0) ) - 3, 0) 
   elseif(ityp==LIST5) then
      factor =      int (log (qqmax) / log (10.0D0) ) - 3
   endif
   out_fac = 10** ( - factor) 
elseif(valmax>  0.0 .and. (ityp.eq.HKLF4.or.ityp.eq.LIST5)) THEN    ! Scale with User value
  out_fac = valmax
ENDIF 
!write(*,*) ' OUT_FAC ', out_fac, factor, qqmax, vect
!                                                                       
extr_ima = 6 - out_extr_abs - out_extr_ord 
!                                                                       
IF(ityp.eq.0) THEN      ! A standard file, allocate temporary arrays
!                               Write data to temporary data structure 
!                               This allows to copy to KUPLOT
   IF(    out_inc(2) == 1 .and. out_inc(3) == 1 ) THEN ! 1D file along axis 1
      is_axis = 1  ! Axis is 1
      is_dim  = 1  ! is a 1D file
   ELSEIF(out_inc(1) == 1 .and. out_inc(3) == 1 ) THEN ! 1D file along axis 2
      is_axis = 2  ! Axis is 2
      is_dim  = 1  ! is a 1D file
   ELSEIF(out_inc(1) == 1 .and. out_inc(2) == 1 ) THEN ! 1D file along axis 3
      is_axis = 3  ! Axis is 3
      is_dim  = 1  ! is a 1D file
   ELSEIF(out_inc(3) == 1 ) THEN                       ! 2d Normal to axis 3
      is_axis = 3  ! Axis is 3
      is_dim  = 2  ! is a 2D file
      npkt1   = out_inc(1)
      npkt2   = out_inc(2)
      ranges(1) = out_eck (out_extr_abs, 1)
      ranges(2) = out_eck (out_extr_abs, 2)
      ranges(3) = out_eck (out_extr_ord, 1)
      ranges(4) = out_eck (out_extr_ord, 3)
!  ELSEIF(out_inc(2) == 1 ) THEN                       ! 2d Normal to axis 2
!     is_axis = 2  ! Axis is 2
!     is_dim  = 2  ! is a 2D file
!     npkt1   = out_inc(1)
!     npkt2   = out_inc(3)
!     ranges(1) = out_eck (out_extr_abs, 1)
!     ranges(2) = out_eck (out_extr_abs, 2)
!     ranges(3) = out_eck (out_extr_ord, 1)
!     ranges(4) = out_eck (out_extr_ord, 4)
!  ELSEIF(out_inc(1) == 1 ) THEN                       ! 2d Normal to axis 1
!     is_axis = 1  ! Axis is 1
!     is_dim  = 2  ! is a 2D file
!     npkt1   = out_inc(2)
!     npkt2   = out_inc(3)
!     ranges(1) = out_eck (out_extr_abs, 1)
!     ranges(2) = out_eck (out_extr_abs, 3)
!     ranges(3) = out_eck (out_extr_ord, 1)
!     ranges(4) = out_eck (out_extr_ord, 4)
   ELSE
      is_dim  = 3
      npkt1   = out_inc(1)
      npkt2   = out_inc(2)
      npkt3   = out_inc(3)
      ranges(1) = out_eck (out_extr_abs, 1)
      ranges(2) = out_eck (out_extr_abs, 2)
      ranges(3) = out_eck (out_extr_ord, 1)
      ranges(4) = out_eck (out_extr_ord, 3)
   ENDIF
!
   IF(is_dim==1) THEN                           ! 1D output
      out_index(1) = out_extr_abs
      out_index(2) = out_extr_ord
      out_index(3) =     extr_ima
      npkt1 = out_inc(is_axis)
      ALLOCATE(xwrt(1:npkt1), STAT=all_status)  ! Allocate x-table
      ALLOCATE(ywrt(1:npkt1), STAT=all_status)  ! Allocate y-table
      loop = 1                                  ! Preset all loop indices to 1
      j    = 1
      DO i = 1, out_inc (is_axis)   ! loop along axis is_axis 
         loop(is_axis) = i
         DO k = 1, 3 
            h (k) = out_eck(k,1) + out_vi(k,1) * REAL(loop(1)-1, kind=PREC_DP)   &
                                 + out_vi(k,2) * REAL(loop(2)-1, kind=PREC_DP)   &
                                 + out_vi(k,3) * REAL(loop(3)-1, kind=PREC_DP)  
         ENDDO 
         xwrt(i) = h(out_extr_abs)
         ywrt(i) = qval (i, 1, 1, value, i, j, laver)
      ENDDO 
      CALL output_save_file_1d(outfile, npkt1, xwrt, ywrt, out_mode)
      DEALLOCATE(xwrt)
      DEALLOCATE(ywrt)
   ELSEIF(is_dim==2) THEN                       ! 2D output
      ALLOCATE(zwrt(1:npkt1,1:npkt2), STAT=all_status)  ! Allocate z-table
      l = 1
      DO j = 1, npkt2
         DO i=1,npkt1
!           zwrt(i,j) = (qval ( (i - 1) * out_inc(3)*out_inc (2) +        &
!                               (j - 1) * out_inc(3)             + l,     &
!                        value,  i, j, laver))
            zwrt(i,j) = (qval (i,j,1, & 
                         value,  i, j, laver))
         ENDDO 
      ENDDO 
      nnew1 = NPKT1
      nnew2 = NPKT2
!     IF(value==val_3Dpdf) THEN
!        nnew1 = 201
!        nnew2 = 201
!        pdf3d_eck(1,1) = -2.0
!        pdf3d_eck(2,1) = -2.0
!        pdf3d_eck(3,1) =  0.0
!        pdf3d_eck(1,2) =  2.0
!        pdf3d_eck(2,2) = -2.0
!        pdf3d_eck(3,2) =  0.0
!        pdf3d_eck(1,3) = -2.0
!        pdf3d_eck(2,3) =  2.0
!        pdf3d_eck(3,3) =  0.0
!        pdf3d_inc(1) = nnew1
!        pdf3d_inc(2) = nnew2
!        pdf3d_vi(1,1) = (pdf3d_eck(1,2)-pdf3d_eck(1,1))/REAL(nnew1-1)
!        pdf3d_vi(2,1) = (pdf3d_eck(2,2)-pdf3d_eck(2,1))/REAL(nnew1-1)
!        pdf3d_vi(3,1) = (pdf3d_eck(3,2)-pdf3d_eck(3,1))/REAL(nnew1-1)
!        pdf3d_vi(1,2) = (pdf3d_eck(1,3)-pdf3d_eck(1,1))/REAL(nnew2-1)
!        pdf3d_vi(2,2) = (pdf3d_eck(2,3)-pdf3d_eck(2,1))/REAL(nnew2-1)
!        pdf3d_vi(3,2) = (pdf3d_eck(3,3)-pdf3d_eck(3,1))/REAL(nnew2-1)
!        ALLOCATE(znew(nnew1, nnew2))
!        znew(:,:) = 0.0D0
!        CALL do_fft_2d_cos(npkt1, npkt2, zwrt, out_eck, out_vi, out_inc, &
!                           nnew1, nnew2, znew, pdf3d_eck, pdf3d_vi, pdf3d_inc)
!        DEALLOCATE(zwrt)
!        ALLOCATE(zwrt(nnew1, nnew2))
!        zwrt(:,:) = znew(:,:)
!        ranges(1) = pdf3d_eck (1, 1)
!        ranges(2) = pdf3d_eck (1, 2)
!        ranges(3) = pdf3d_eck (2, 1)
!        ranges(4) = pdf3d_eck (2, 3)
!     ENDIF
      CALL write_discus_nipl_header(header_lines, nheader, l)
      CALL output_save_file_2d(outfile, ranges, nnew1, nnew2, zwrt,       &
                               header_lines, nheader, out_mode)
      DEALLOCATE(header_lines)
      DEALLOCATE(zwrt)
   ELSEIF(is_dim==3) THEN                       ! 3D output into standard slices
      ALLOCATE(zwrt(1:npkt1,1:npkt2), STAT=all_status)  ! Allocate z-table
      DO l = 1, npkt3                           ! For all layers along 3rd axis
         CALL write_discus_nipl_header(header_lines, nheader, l)
         WRITE(dummy_file, 7777) outfile(1:len_str(outfile)),l  ! Modify file name
7777 FORMAT(a,'.PART_',i4.4)
         DO j = 1, npkt2                        ! Loop over points in 2D
            DO i=1,npkt1                        ! and copy into intensity file
!              zwrt(i,j) = (qval ( (i - 1) * out_inc(3)*out_inc (2) +        &
!                                  (j - 1) * out_inc(3)             + l,     &
!                           value,  i, j, laver))
               zwrt(i,j) = (qval (i,j,l, & 
                            value,  i, j, laver))
            ENDDO 
         ENDDO 
         CALL output_save_file_2d(dummy_file, ranges, npkt1, npkt2, zwrt,    &
                                  header_lines, nheader, out_mode)
         DEALLOCATE(header_lines)
      ENDDO 
      DEALLOCATE(zwrt)
   ENDIF
ELSE      ! Data types ityp==0 or ELSE ! Block for all but standard file formats
      lread = .false. 
      IF(.not.((out_inc(3) > 1 .and. ityp.eq.0) .OR.    &     ! NOT ( multiple layers in standard file typ OR
               ityp==13                             )) THEN   !       HDF5 file format                        )
         CALL oeffne (iff, outfile, 'unknown') 
      ENDIF
      IF (ier_num.eq.0) THEN 
         IF (out_inc (1) .gt.1.and.out_inc (2) .gt.1) THEN  ! 2D or 3D data
            IF (ityp.eq.0) THEN                             ! Standard file format
!               DO l=1, out_inc(3)
!                  IF(out_inc(3) > 1) THEN
!                     WRITE(dummy_file, 7777) outfile(1:len_str(outfile)),l
!7777 FORMAT(a,'.PART_',i4.4)
!                     CALL oeffne (iff, dummy_file, 'unknown') 
!                  ENDIF
!               WRITE (iff, * ) out_inc (1), out_inc(2)
!               WRITE (iff, * ) out_eck (out_extr_abs, 1), out_eck (out_extr_abs, 2), &
!                               out_eck (out_extr_ord, 1), out_eck (out_extr_ord, 3)
!               DO j = 1, out_inc (2) 
!               WRITE (iff, 4) (qval ( (i - 1) * out_inc(3)*out_inc (2) +        &
!                                      (j - 1) * out_inc(3)             + l,     &
!                                      value,  i, j, laver), i = 1, out_inc (1) )
!               WRITE (iff, 100) 
!               ENDDO 
!                  IF(out_inc(3) > 1) THEN
!                     CLOSE(iff)
!                  ENDIF
!               ENDDO 
            ELSEIF (ityp.eq.ASCII3D) THEN                   ! 3D "NIPL" file
               WRITE (iff, * ) out_inc (1), out_inc(2), out_inc(3)
               WRITE (iff, * ) out_eck (out_extr_abs, 1), out_eck (out_extr_abs, 2), &
                               out_eck (out_extr_ord, 1), out_eck (out_extr_ord, 3), &
                               out_eck (out_extr_top, 1), out_eck (out_extr_top, 4)
               DO l=1, out_inc(3)
                  DO j = 1, out_inc (2) 
!                 WRITE (iff, 4) (qval ( (i - 1) * out_inc(3)*out_inc (2) +        &
!                                        (j - 1) * out_inc(3)             + l,     &
!                                        value,  i, j, laver), i = 1, out_inc (1) )
                  WRITE (iff, 4) (qval (i,j,l, value,  i, j, laver), i = 1, out_inc (1) )
                  WRITE (iff, 100) 
                  ENDDO 
                  IF(out_inc(3) > 1) THEN
                     CLOSE(iff)
                  ENDIF
               ENDDO 
            ELSEIF (ityp.eq.HKLF4) THEN                     ! SHELXS HKL File
               DO l = 1, out_inc (3) 
               DO j = 1, out_inc (2) 
               DO i = 1, out_inc (1) 
               DO k = 1, 3 
               h (k) = out_eck (k, 1) + out_vi (k, 1) * REAL(i - 1)   &
                                      + out_vi (k, 2) * REAL(j - 1)   &
                                      + out_vi (k, 3) * REAL(l - 1)
               ENDDO 
            IF( (INT(h(1)))**2 + (INT(h(2)))**2 + (INT(h(3)))**2 /= 0 ) THEN
!
! Extinction correction is calculated using Amplitude, which has been scaled to one unit cell
! Intensity written as (Amplitude*extinction)**2 * scale_factor
               qq = qval(i,j,l, 2    , i, j, laver) / cr_icc(1) / cr_icc(2) / cr_icc(3)
               ext_cor = (1.0_PREC_DP/(1.0_PREC_DP + 0.001_PREC_DP * diff_exti * qq**2*rlambda**3/ &
                         sind(2.0*asind(0.5_PREC_DP*rlambda*do_blen(.false., h, NULLV))))**0.25_PREC_DP)!**2.0_PREC_DP
!if(nint(h(1))==0 .and. nint(h(2))==0 .and. nint(h(3)) ==-1) then
!write(*,*) h, qq, ext_cor, 2.0*asind(0.5_PREC_DP*rlambda*do_blen(.false., h, NULLV))
!write(*,*) cr_icc, ' | ', value, qval (i,j,l, value, i, j, laver), qval (i,j,l, 2, i, j, laver)/ cr_icc (1) / cr_icc (2) / cr_icc (3)
!write(*,*) 'WRITE ', (qq * ext_cor)**2 * out_fac, ' | ', qq, out_fac
!endif
               sq = max(qq * sqrt (ext_cor) * out_fac *0.10_PREC_DP, 0.01_PREC_DP)
               qq = (qq * ext_cor)**2 * out_fac
               WRITE (iff, '(3i4,2f8.2)') int(h(1)), int(h(2)), int(h(3)), qq, sq
            ENDIF
               ENDDO 
               ENDDO 
               ENDDO 
            ELSEIF (ityp.eq.LIST5) THEN                     ! SHELXS HKL Fobs Fcalc File
               DO l = 1, out_inc (3) 
               DO j = 1, out_inc (2) 
               DO i = 1, out_inc (1) 
               DO k = 1, 3 
               h (k) = out_eck (k, 1) + out_vi (k, 1) * REAL(i - 1)   &
                                      + out_vi (k, 2) * REAL(j - 1)   &
                                      + out_vi (k, 3) * REAL(l - 1)
               ENDDO 
            IF( (INT(h(1)))**2 + (INT(h(2)))**2 + (INT(h(3)))**2 /= 0 ) THEN
!
               shel_value = 2                         ! Calculate amplitude
               qq = qval(i,j,l, shel_value, i, j, laver) / cr_icc(1) / cr_icc(2) / cr_icc(3) ! Ampltude normalized to one unit cell
               ext_cor = (1.0_PREC_DP/(1.0_PREC_DP + 0.001_PREC_DP * diff_exti * qq**2*rlambda**3/ &
                         sind(2.0*asind(0.5_PREC_DP*rlambda*do_blen(.false., h, NULLV))))**0.25_PREC_DP)!**2.0_PREC_DP
               qq = (qq * ext_cor) * out_fac          ! Amplitude * Extinction * scale
               shel_value = 3                         ! Calculate phase angle
               sq = qval(i,j,l, shel_value, i, j, laver) 
               IF(sq < 0.0D0 ) sq = sq + 360.0D0
               WRITE (iff, 8) int (h (1) ), int (h (2) ), int (h (3) ), qq, qq, sq
            ENDIF
               ENDDO 
               ENDDO 
               ENDDO 
            ELSEIF (ityp.eq.LIST9) THEN                     ! SHELXS File
               shel_value = 3 
               DO l = 1, out_inc (3) 
               DO j = 1, out_inc (2) 
               DO i = 1, out_inc (1) 
               DO k = 1, 3 
               h (k) = out_eck (k, 1) + out_vi (k, 1) * REAL(i - 1)   &
                                      + out_vi (k, 2) * REAL(j - 1)   &
                                      + out_vi (k, 3) * REAL(l - 1)
               ENDDO 
            IF( (INT(h(1)))**2 + (INT(h(2)))**2 + (INT(h(3)))**2 /= 0 ) THEN
!              k  = (i - 1) * out_inc (3) * out_inc (2) + (j-1) * out_inc (3) + l 
               shel_value = 2 
!              qq = qval (k, value, i, j, laver) / cr_icc (1) / cr_icc (2) / cr_icc (3)
!              qq = qval (i,j,l, value, i, j, laver) / cr_icc (1) / cr_icc (2) / cr_icc (3)
               shel_value = 2                         ! Calculate amplitude
               qq = qval(i,j,l, shel_value, i, j, laver) / cr_icc(1) / cr_icc(2) / cr_icc(3) ! Ampltude normalized to one unit cell
               ext_cor = (1.0_PREC_DP/(1.0_PREC_DP + 0.001_PREC_DP * diff_exti * qq**2*rlambda**3/ &
                         sind(2.0*asind(0.5_PREC_DP*rlambda*do_blen(.false., h, NULLV))))**0.25_PREC_DP)!**2.0_PREC_DP
               qq = (qq * ext_cor) * out_fac          ! Amplitude * Extrinction * scale
               shel_value = 3 
               sq = qval (i,j,l, shel_value, i, j, laver) 
               IF(sq < 0.0 ) sq = sq + 360.0
               WRITE (iff, 9) h (1), h (2), h (3), qq, qq, sq 
            ENDIF
               ENDDO 
               ENDDO 
               ENDDO 
            ELSE             ! Should be GNU type == 3      ! Standard 2D File
               DO j = 1, out_inc (2) 
               DO i = 1, out_inc (1) 
               DO k = 1, 3 
               h (k) = out_eck (k, 1) + out_vi (k, 1) * REAL(i - 1)   &
               + out_vi (k, 2) * REAL(j - 1)                          
               ENDDO 
               k = (i - 1) * out_inc (2) + j 
               WRITE (iff, 5) h (out_extr_abs), h (out_extr_ord),       &
               qval (i,j,1, value, i, j, laver), h (extr_ima)               
               ENDDO 
               WRITE (iff, 100) 
               ENDDO 
            ENDIF 
         ELSEIF (out_inc (1) .eq.1) THEN                    ! 1D Files
            IF(ityp /= 0) THEN                              ! All BUT standard files
            i = 1 
            DO j = 1, out_inc (2) 
            DO k = 1, 3 
            h (k) = out_eck (k, 1) + out_vi (k, 1) * REAL(i - 1)      &
            + out_vi (k, 2) * REAL(j - 1)                             
            ENDDO 
!           k = (i - 1) * out_inc (2) + j 
            IF( (INT(h(1)))**2 + (INT(h(2)))**2 + (INT(h(3)))**2 /= 0 ) THEN
            IF (ityp.eq.HKLF4) THEN                         ! SHELXS HKL INTENSITY
!              qq = qval (k, value, i, j, laver) / cr_icc (1) / cr_icc (&
               qq = qval (1,j,1, value, i, j, laver) / cr_icc (1) / cr_icc (&
               2) / cr_icc (3) * out_fac                                
               sq = sqrt (qq) 
               WRITE (iff, 7) int (h (1) ), int (h (2) ), int (h (3) ), &
               qq, sq                                                   
            ELSEIF (ityp.eq.LIST5) THEN                     ! SHELXS Fobs Fcalc
               shel_value = 2 
!              qq = qval (k, value, i, j, laver) / cr_icc (1) / cr_icc (&
!              2) / cr_icc (3) * out_fac                                
               qq = qval (1,j,1, value, i, j, laver) / cr_icc (1) / cr_icc (2) / cr_icc (3) * out_fac                                
               shel_value = 3 
!              sq = qval (k, shel_value, i, j, laver) / cr_icc (1)      &
!              / cr_icc (2) / cr_icc (3) * out_fac                      
               sq = qval (1,j,1, shel_value, i, j, laver) 
               IF(sq < 0.0 ) sq = sq + 360.0
               WRITE (iff, 8) int (h (1) ), int (h (2) ), int (h (3) ), &
               qq, qq, sq                                               
            ELSEIF (ityp.eq.LIST9) THEN                     ! SHELXS
               shel_value = 2 
!              qq = qval (k, value, i, j, laver) / cr_icc (1) / cr_icc (&
!              2) / cr_icc (3)                                          
               qq = qval (1,j,1, value, i, j, laver) / cr_icc (1) / cr_icc (2) / cr_icc (3)
               shel_value = 3 
!              sq = qval (k, shel_value, i, j, laver) / cr_icc (1)      &
!              / cr_icc (2) / cr_icc (3)                                
               sq = qval (1,j,1, shel_value, i, j, laver) 
               IF(sq < 0.0 ) sq = sq + 360.0
               WRITE (iff, 9) h (1), h (2), h (3), qq, qq, sq 
            ENDIF 
            ENDIF 
               ENDDO 
            ELSE     ! Should be GNU type == 3              ! Standard File
               i = 1 
               DO j = 1, out_inc (2) 
                  DO k = 1, 3 
                     h(k) = out_eck(k,1) + out_vi(k,1) * REAL(i-1)    &
                                         + out_vi(k,2) * REAL(j-1)
                  ENDDO 
!                 k       = (i - 1) * out_inc (2) + j 
!                 WRITE(iff,6) h(out_extr_ord), qval(k,value,i,j,laver)
                  WRITE(iff,6) h(out_extr_ord), qval(1,j,1,value,i,j,laver)
                  xwrt(i) = h(out_extr_ord)
!                 ywrt(i) = qval (k, value, i, j, laver)
                  ywrt(i) = qval (1,j,1, value, i, j, laver)
               ENDDO 
            ENDIF 
         ELSEIF (out_inc (2) .eq.1) THEN 
            IF(ityp /= 0) THEN                              ! All BUT standard files
            j = 1 
            DO i = 1, out_inc (1) 
            DO k = 1, 3 
            h (k) = out_eck (k, 1) + out_vi (k, 1) * REAL(i - 1)      &
            + out_vi (k, 2) * REAL(j - 1)                             
            ENDDO 
!           k = (i - 1) * out_inc (2) + j 
            IF( (INT(h(1)))**2 + (INT(h(2)))**2 + (INT(h(3)))**2 /= 0 ) THEN
            IF (ityp.eq.HKLF4) THEN                         ! SHELXS Intensity
!              qq = qval (k, value, i, j, laver) / cr_icc (1) / cr_icc (&
!              2) / cr_icc (3) * out_fac                                
               qq = qval (i,1,1, value, i, j, laver) / cr_icc (1) / cr_icc ( 2) / cr_icc (3) * out_fac                                
               sq = sqrt (qq) 
               WRITE (iff, 7) int (h (1) ), int (h (2) ), int (h (3) ), qq, sq                                                   
            ELSEIF (ityp.eq.LIST5) THEN                     ! SHELS Fobs Fcalc
               shel_value = 2 
!              qq = qval (k, value, i, j, laver) / cr_icc (1) / cr_icc (&
!              2) / cr_icc (3) * out_fac                                
               qq = qval (i,1,1, value, i, j, laver) / cr_icc (1) / cr_icc (2) / cr_icc (3) * out_fac                                
               shel_value = 3                               ! Phase angle
               sq = qval (i,1,1, shel_value, i, j, laver) 
               IF(sq < 0.0 ) sq = sq + 360.0
               WRITE (iff, 8) int (h (1) ), int (h (2) ), int (h (3) ), &
               qq, qq, sq                                               
            ELSEIF (ityp.eq.LIST9) THEN                     ! SHELS File
               shel_value = 2 
!              qq = qval (k, value, i, j, laver) / cr_icc (1) / cr_icc (&
!              2) / cr_icc (3)                                          
               qq = qval (i,1,1, value, i, j, laver) / cr_icc (1) / cr_icc ( 2) / cr_icc (3)                                          
               shel_value = 3 
!              sq = qval (k, shel_value, i, j, laver) / cr_icc (1)      &
!              / cr_icc (2) / cr_icc (3)                                
!              sq = qval (k, shel_value, i, j, laver) 
               sq = qval (i,1,1, shel_value, i, j, laver) 
               IF(sq < 0.0 ) sq = sq + 360.0
               WRITE (iff, 9) h (1), h (2), h (3), qq, qq, sq 
            ENDIF 
            ENDIF 
            ENDDO 
            ELSE                                            ! Standard File
!              j = 1 
!              DO i = 1, out_inc (1) 
!                 DO k = 1, 3 
!                    h (k) = out_eck(k,1) + out_vi(k,1) * REAL(i-1)   &
!                                         + out_vi(k,2) * REAL(j-1)
!                 ENDDO 
!                 k = (i - 1) * out_inc (2) + j 
!                 WRITE(iff,6) h(out_extr_abs), qval(k,value,i,j,laver)
!              ENDDO 
            ENDIF 
         ENDIF 
      ENDIF 
!     if(ier_num.ne.0) THEN                                             
!       call errlist                                                    
!     endif                                                             
      IF (ityp.eq.HKLF4.or.ityp.eq.LIST5) THEN 
         WRITE (output_io, 1000) out_fac 
      ENDIF 
      CLOSE (iff) 
      ENDIF       ! DATA TYPES ityp == 0 or else 
      IF(ALLOCATED(xwrt)) DEALLOCATE(xwrt,STAT=all_status)
      IF(ALLOCATED(ywrt)) DEALLOCATE(ywrt,STAT=all_status)
!                                                                       
    4 FORMAT (5(1x,e12.5)) 
    5 FORMAT (4(1x,e12.5)) 
    6 FORMAT (2(1x,e12.5)) 
    7 FORMAT (3i4,2f8.2) 
    8 FORMAT (3i4,2f10.2,f7.2) 
    9 FORMAT (3(f10.6,1x),2(e12.5,1x),f7.2) 
  100 FORMAT () 
 1000 FORMAT    (' Data have been scaled by ',g18.8e3) 
!100      format(/)                                                     
!                                                                       
      END SUBROUTINE do_output                      
!
!*****7*****************************************************************
      SUBROUTINE set_output (linverse) 
!-                                                                      
!     Sets the proper output values for either Fourier or               
!     inverse Fourier and Patterson                                     
!+                                                                      
      USE discus_config_mod 
      USE diffuse_mod 
      USE output_mod 
      USE patters_mod 
      IMPLICIT none 
!                                                                       
!                                                                       
      LOGICAL linverse 
!                                                                       
      INTEGER i, j 
!                                                                       
      IF (linverse) THEN 
         out_extr_abs = rho_extr_abs 
         out_extr_ord = rho_extr_ord 
!                                                                       
         DO i = 1, 3 
            DO j = 1, 3 
               out_eck (i, j) = rho_eck (i, j) 
            ENDDO 
         ENDDO 
!                                                                       
         DO i = 1, 2 
            DO j = 1, 3 
               out_vi (j, i) = rho_vi (j, i) 
            ENDDO 
            out_inc (i) = rho_inc (i) 
         ENDDO 
         out_inc(3) = 1
      ELSE 
         out_extr_abs = extr_abs 
         out_extr_ord = extr_ord 
         out_extr_top = extr_top 
!                                                                       
         DO i = 1, 3 
         DO j = 1, 4 
         out_eck (i, j) = eck (i, j) 
         ENDDO 
         ENDDO 
!                                                                       
         DO i = 1, 3 
         DO j = 1, 3 
         out_vi (j, i) = vi (j, i) 
         ENDDO 
         out_inc (i) = inc (i) 
         ENDDO 
      ENDIF 
!
      END SUBROUTINE set_output                     
!
!*******************************************************************************
!
SUBROUTINE out_prep_3dpdf(laver, l_val_limited, dsmax)
!-
!  Calculate the 3DPDF value via FFT
!+
!
USE crystal_mod
USE diffuse_mod
USE metric_mod
USE output_mod
!
USE errlist_mod
USE param_mod
USE precision_mod
!
IMPLICIT NONE
!
LOGICAL           , INTENT(IN) :: laver
LOGICAL           , INTENT(IN) :: l_val_limited
REAL(KIND=PREC_DP), INTENT(IN) :: dsmax
!
real(kind=PREC_DP), parameter :: EPS = 1.0e-6
INTEGER :: isdim, i,j
INTEGER, DIMENSION(3) :: dsort
CHARACTER(LEN=PREC_STRING)   :: string
INTEGER               :: lcomm
REAL(KIND=PREC_DP), DIMENSION(3)    :: u
REAL(KIND=PREC_DP)    :: uu
!
IF(l_val_limited) THEN
   CALL output_limit(dsmax)
ENDIF
!
if(allocated(rpdf)) deallocate(rpdf)
allocate(rpdf(1:num(1), 1:num(2), 1:num(3)))
!
!  Determine sequence of array dimensions 
!
dsort(1)      = MAXLOC(num, 1)
num(dsort(1)) = -num(dsort(1))
dsort(2)      = MAXLOC(num, 1)
num(dsort(2)) = -num(dsort(2))
dsort(3)      = MAXLOC(num, 1)
num(dsort(3)) = -num(dsort(3))
num = -num
!
!  Determine dimensions that we need 
!
isdim = 3
IF(num(1)==1) isdim = isdim - 1
IF(num(2)==1) isdim = isdim - 1
IF(num(3)==1) isdim = isdim - 1
!read(*,*) i
IF(isdim==1) THEN
   CALL out_prep_3dpdf_1d(laver, dsort)
ELSEIF(isdim==2) THEN
   CALL out_prep_3dpdf_2d(laver, dsort)
ELSEIF(isdim==3) THEN
   CALL out_prep_3dpdf_3d(laver, dsort)
ENDIF
!
! Set increments and corners for out_vi, out_eck
!
!write(*,*) ' num ', num
!write(*,*) ' vi   ', vi(:,1)
!write(*,*) ' vi   ', vi(:,2)
!write(*,*) ' vi   ', vi(:,3)
DO i = 1, 3
   u(1) = INT((num(1)-1))*1.00D0*vi(1,i)
   u(2) = INT((num(2)-1))*1.00D0*vi(2,i)
   u(3) = INT((num(3)-1))*1.00D0*vi(3,i)
   uu = skalpro (u, u, cr_rten)
   IF( uu > 0.0) THEN
      WRITE(string,'(2(F16.9,'',''), F16.9)') u
      lcomm = LEN_TRIM(string)
      CALL d2r(string, lcomm, .FALSE.)
!write(*,*) 'ier_num ', ier_num, ier_typ
      out_vi(1,i) = res_para(4)
      out_vi(2,i) = res_para(5)
      out_vi(3,i) = res_para(6)
!write(*,*) ' res ', res_para(0:3)
!write(*,*) ' res ', res_para(4:6)
   ENDIF
ENDDO
out_vi      = 0.0D0
do i=1, 3
  do j= 1, 3
    if(abs(vi(j,i)) > EPS) out_vi(j,i) = 1.0D0/(vi(j,i)*(num(i)-1))
  enddo
enddo
!out_vi(1,1) = 0.125D0
!out_vi(2,2) = 0.125D0
!out_vi(3,3) = 0.125D0
DO i=1, 3
   out_eck(i,1) = (- out_vi(i,1)*INT((num(1)-1)/2) - out_vi(i,2)*INT((num(2)-1)/2) - out_vi(i,3)*INT((num(3)-1)/2))  ! lower left
   out_eck(i,2) = (+ out_vi(i,1)*INT((num(1)-1)/2) - out_vi(i,2)*INT((num(2)-1)/2) - out_vi(i,3)*INT((num(3)-1)/2))  ! lower left
   out_eck(i,3) = (- out_vi(i,1)*INT((num(1)-1)/2) + out_vi(i,2)*INT((num(2)-1)/2) - out_vi(i,3)*INT((num(3)-1)/2))  ! lower left
   out_eck(i,4) = (- out_vi(i,1)*INT((num(1)-1)/2) - out_vi(i,2)*INT((num(2)-1)/2) + out_vi(i,3)*INT((num(3)-1)/2))  ! lower left
ENDDO
!write(*,'(a,3i5   )') ' num ', num
!write(*,'(a,3f12.5)') '     abs  ',     vi(:,1)
!write(*,'(a,3f12.5)') '     ord  ',     vi(:,2)
!write(*,'(a,3f12.5)') '     top  ',     vi(:,3)
!write(*,'(a,3f12.5)') ' out_abs  ', out_vi(:,1)
!write(*,'(a,3f12.5)') ' out_ord  ', out_vi(:,2)
!write(*,'(a,3f12.5)') ' out_top  ', out_vi(:,3)
!write(*,'(a,3f12.5)') 'eck_ll  ', out_eck(:,1)
!write(*,'(a,3f12.5)') 'eck_lr  ', out_eck(:,2)
!write(*,'(a,3f12.5)') 'eck_ul  ', out_eck(:,3)
!write(*,'(a,3f12.5)') 'eck_tl  ', out_eck(:,4)
!write(*,'(a,3f12.5)') 'cal lru ', out_eck(1,1)+(num(1)-1)*out_vi(1,1)
!write(*,'(a,3f12.5)') 'cal ulv ', out_eck(2,1)+(num(2)-1)*out_vi(2,2)
!write(*,'(a,3f12.5)') 'cal tlw ', out_eck(3,1)+(num(3)-1)*out_vi(3,3)
!read(*,*) i
!
END SUBROUTINE out_prep_3dpdf
!
!*******************************************************************************
!
SUBROUTINE out_prep_3dpdf_1d(laver, dsort)
!-
!  Calculate the 3DPDF value via FFT
!+
!
USE diffuse_mod
!
USE errlist_mod
USE precision_mod
USE map_1dtofield
use lib_f90_fftw3
!
IMPLICIT NONE
!
LOGICAL              , INTENT(IN) :: laver
INTEGER, DIMENSION(3), INTENT(IN) :: dsort
!integer :: i
!
!COMPLEX(KIND=KIND(0.0D0)) , DIMENSION(:), ALLOCATABLE  :: pattern  ! the diffraction pattern
COMPLEX(KIND=C_DOUBLE_COMPLEX) , DIMENSION(:), ALLOCATABLE  :: in_pattern  ! the diffraction pattern
COMPLEX(KIND=C_DOUBLE_COMPLEX) , DIMENSION(:), ALLOCATABLE  ::out_pattern  ! the diffraction pattern
type(c_ptr) :: plan    ! FFWT3 plan
!
ALLOCATE(in_pattern(num(dsort(1))))
ALLOCATE(out_pattern(num(dsort(1))))
!
CALL maptofftfd(num, dsort, dsi(1:num(1),1,1), in_pattern)
!
!pattern = fft(pattern) / SQRT(REAL(num(1)))
!
plan = fftw_plan_dft_1d(num(1)        , in_pattern, out_pattern, FFTW_FORWARD, FFTW_ESTIMATE)
call   fftw_execute_dft(plan, in_pattern, out_pattern)
call   fftw_destroy_plan(plan)
!
out_pattern = out_pattern / sqrt(real(num(1), kind=PREC_DP))  ! Normalize PDF
!
CALL mapfftfdtoline(num, dsort, rpdf(1:num(1),1,1), out_pattern)
!
DEALLOCATE( in_pattern)
DEALLOCATE(out_pattern)
!
END SUBROUTINE out_prep_3dpdf_1d
!
!*******************************************************************************
!
SUBROUTINE out_prep_3dpdf_2d(laver, dsort)
!-
!  Calculate the 3DPDF value via FFT
!  Uses FFTW3 library
!+
!
USE diffuse_mod
!
USE errlist_mod
USE precision_mod
USE map_1dtofield
use lib_f90_profile
use lib_f90_fftw3
!
IMPLICIT NONE
!
LOGICAL              , INTENT(IN) :: laver
INTEGER, DIMENSION(3), INTENT(IN) :: dsort
!
COMPLEX(KIND=C_DOUBLE_COMPLEX)  , DIMENSION(:,:), ALLOCATABLE  :: in_pattern  ! the diffraction pattern!
COMPLEX(KIND=C_DOUBLE_COMPLEX)  , DIMENSION(:,:), ALLOCATABLE  :: out_pattern  ! the diffraction pattern!
type(c_ptr) :: plan    ! FFWT3 plan
!
ALLOCATE( in_pattern(num(dsort(1)) , num(dsort(2)) ))
ALLOCATE( out_pattern(num(dsort(1)), num(dsort(2)) ))
!
CALL maptofftfd(num, dsort, dsi(1:num(1),1:num(2),1), in_pattern)
!
plan = fftw_plan_dft_2d(num(dsort(2)), num(dsort(1)), in_pattern, out_pattern, FFTW_FORWARD, FFTW_ESTIMATE)
call   fftw_execute_dft(plan, in_pattern, out_pattern)
call   fftw_destroy_plan(plan)
!
out_pattern = out_pattern / sqrt(real(num(dsort(2))*num(dsort(1)), kind=PREC_DP))  ! Normalize PDF
!
CALL mapfftfdtoline(num, dsort, rpdf(1:num(1),1:num(2),1), out_pattern)
!
DEALLOCATE(in_pattern)
DEALLOCATE(out_pattern)
!
END SUBROUTINE out_prep_3dpdf_2d
!
!*******************************************************************************
!
SUBROUTINE out_prep_3dpdf_3d(laver, dsort)
!-
!  Calculate the 3DPDF value via FFT  3D version
!+
!
USE diffuse_mod
!
USE errlist_mod
USE precision_mod
USE map_1dtofield
use lib_f90_fftw3
!
IMPLICIT NONE
!
LOGICAL              , INTENT(IN) :: laver
INTEGER, DIMENSION(3), INTENT(IN) :: dsort
!
COMPLEX(KIND=KIND(0.0D0)) , DIMENSION(:,:,:), ALLOCATABLE  ::  in_pattern  ! the diffraction pattern
COMPLEX(KIND=KIND(0.0D0)) , DIMENSION(:,:,:), ALLOCATABLE  :: out_pattern  ! the diffraction pattern
type(c_ptr) :: plan    ! FFWT3 plan
!
ALLOCATE( in_pattern(num(dsort(1)), num(dsort(2)), num(dsort(3))))
ALLOCATE(out_pattern(num(dsort(1)), num(dsort(2)), num(dsort(3))))
!
CALL maptofftfd(num, dsort, dsi, in_pattern)
!
plan = fftw_plan_dft_3d(num(dsort(3)), num(dsort(2)), num(dsort(1)), in_pattern, out_pattern, FFTW_FORWARD, FFTW_ESTIMATE)
call   fftw_execute_dft(plan, in_pattern, out_pattern)
!
out_pattern = out_pattern / sqrt(real(num(1)*num(2)*num(3),kind=PREC_DP))
call   fftw_destroy_plan(plan)
!
CALL mapfftfdtoline(num, dsort, rpdf, out_pattern)
!
DEALLOCATE( in_pattern)
DEALLOCATE(out_pattern)
!
END SUBROUTINE out_prep_3dpdf_3d
!
!*******************************************************************************
!
SUBROUTINE output_reset
!
USE output_mod
!
IMPLICIT NONE
!
outfile      = 'fcalc.dat'
ityp         = 0
extr_abs     = 1
extr_ord     = 2
extr_top     = 3
rho_extr_abs = 1
rho_extr_ord = 2
out_extr_abs = 1
out_extr_ord = 2
out_extr_top = 3
out_inc(:)   = (/121, 121, 1/)
out_eck      = reshape((/ 0.0, 0.0,  0.0, &
                          5.0, 0.0,  0.0, &
                          0.0, 5.0,  0.0, &
                          0.0, 0.0,  0.0/),shape(out_eck))
out_vi       = reshape((/0.05, 0.00, 0.00, &
                         0.0 , 0.05, 0.00, &
                         0.00, 0.00, 0.00/),shape(out_vi))
cpow_form    = 'tth'
out_user_limits = .false.
out_user_values(:) = (/1.0, 10.0, 0.01/)
out_lrange  = 0
out_lcenter = 0
out_lpixel  = 0
out_center  = 0
out_pixel   = 0
!
END SUBROUTINE output_reset
!
!*******************************************************************************
!
SUBROUTINE output_limit(dsmax)
!-
!  Limits reciprocal space to a maximum dstar value
!+
USE crystal_mod
USE diffuse_mod
USE metric_mod
USE output_mod
!
IMPLICIT NONE
!
REAL(KIND=PREC_DP), INTENT(IN) :: dsmax
!
REAL(KIND=PREC_DP), DIMENSION(3), PARAMETER :: NULLV = (/0.0D0, 0.0D0, 0.0D0/)
!
!REAL(KIND=PREC_DP)               :: dstarmax
REAL(KIND=PREC_DP), DIMENSION(3) :: h
INTEGER :: i,j,k,l
INTEGER :: iii
!
DO l = 1, inc(3) 
   DO j = 1, inc(2) 
      DO i = 1, inc(1) 
         DO k = 1, 3 
            h (k) = eck(k, 1) + vi(k, 1) * REAL(i - 1)   &
                              + vi(k, 2) * REAL(j - 1)   &
                              + vi(k, 3) * REAL(l - 1)
!           IF(do_blen(.FALSE., h, NULLV)> dstarmax) THEN
            IF(do_blen(.FALSE., h, NULLV)> dsmax   ) THEN
               iii       =  (i - 1) * inc(3)*inc (2) +        &
                            (j - 1) * inc(3)             + l
               acsf(i,j,k) = CMPLX(0.0D0, 0.0D0)                     ! How does iii runs? How to rewrite this line for general dimensions?
                csf(i,j,k) = CMPLX(0.0D0, 0.0D0)                     ! How does iii runs? How to rewrite this line for general dimensions?
                dsi(i,j,k) = CMPLX(0.0D0, 0.0D0)                     ! How does iii runs? How to rewrite this line for general dimensions?
            ENDIF
         ENDDO
      ENDDO
   ENDDO
ENDDO
!
END SUBROUTINE output_limit
!
!*******************************************************************************
!
END MODULE output_menu
