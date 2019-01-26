MODULE demolec
!
USE demolec_mod
PRIVATE
PUBLIC :: demolecularize
!
CONTAINS
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
SUBROUTINE demolecularize
!
! Menu to remove a molecule status
!
USE calc_expr_mod
USE class_macro_internal
USE discus_allocate_appl_mod
!USE do_read_number_mod
USE do_eval_mod
USE do_wait_mod
USE doact_mod
USE errlist_mod
USE learn_mod
USE molecule_mod
USE prompt_mod
USE sup_mod
!
IMPLICIT NONE
!
CHARACTER (LEN=5)                       :: befehl! command on input line
CHARACTER(LEN=LEN(prompt))              :: orig_prompt  ! original prompt
CHARACTER (LEN=1024)                    :: line  ! input line
CHARACTER (LEN=1024)                    :: zeile ! remainder with parameters
INTEGER                                 :: indxg ! location of "="
INTEGER                                 :: lp    ! length of zeile
INTEGER                                 :: laenge
INTEGER                                 :: lbef
INTEGER                                 :: nmoletype ! temp variable for allocation
!
LOGICAL                                 :: sel   ! condition of select or deselect
LOGICAL                                 :: lend  ! condition of EOF
!
!INTEGER, PARAMETER :: NOPTIONAL = 4
!INTEGER, PARAMETER :: O_CX      = 1
!INTEGER, PARAMETER :: O_CY      = 2
!INTEGER, PARAMETER :: O_Cz      = 3
!INTEGER, PARAMETER :: O_KEEP    = 4
!CHARACTER(LEN=1024), DIMENSION(NOPTIONAL) :: oname   !Optional parameter names
!CHARACTER(LEN=1024), DIMENSION(NOPTIONAL) :: opara   !Optional parameter strings returned
!INTEGER            , DIMENSION(NOPTIONAL) :: loname  !Lenght opt. para name
!INTEGER            , DIMENSION(NOPTIONAL) :: lopara  !Lenght opt. para name returned
!LOGICAL            , DIMENSION(NOPTIONAL) :: lpresent!opt. para present
!REAL               , DIMENSION(NOPTIONAL) :: owerte   ! Calculated values
!INTEGER, PARAMETER                        :: ncalc = 3 ! Number of values to calculate 
!
!DATA oname  / 'cx', 'cy',  'cz',  'keep'   /
!DATA loname /  2,    2,     2  ,   4       /
!
!
INTEGER, EXTERNAL :: len_str
LOGICAL, EXTERNAL :: str_comp
!
orig_prompt = prompt
prompt = prompt (1:len_str (prompt) ) //'/demol'
!
nmoletype = MAX(DEM_MAX_MOLETYPE, mole_num_type)
IF(nmoletype>DEM_MAX_MOLETYPE) THEN
   CALL alloc_demol(nmoletype)
ENDIF
!
main_loop: DO
   CALL no_error
   CALL get_cmd (line, laenge, befehl, lbef, zeile, lp, prompt)
   no_err: IF(ier_num.eq.0) THEN
      no_com: IF(line /= ' '      .AND. line(1:1) /= '#' .AND.      &
                 line /= char(13) .AND. line(1:1) /= '!'        ) THEN
!                                                                       
!     ----search for "="                                                
!                                                                       
         indxg = index (line, '=') 
         is_math: IF(indxg.ne.0                                             &
                     .AND..NOT. (str_comp (befehl, 'echo', 2, lbef, 4) )    &
                     .AND..NOT. (str_comp (befehl, 'syst', 2, lbef, 4) )    &
                     .AND..NOT. (str_comp (befehl, 'help', 2, lbef, 4) .OR. &
                                 str_comp (befehl, '?   ', 2, lbef, 4) )    &
                     .AND. INDEX(line,'==') == 0                            ) THEN
!                                                                       
! ------evaluate an expression and assign the value to a variabble      
!                                                                       
               CALL do_math (line, indxg, laenge)
         ELSE                                    ! is_math, al other commands
!                                                                       
!------ ----execute a macro file                                        
!                                                                       
              IF (befehl (1:1) .eq.'@') THEN     ! macro, reset or all other commands
                  IF (laenge.ge.2) THEN 
                     CALL file_kdo (line (2:laenge), laenge-1) 
                  ELSE 
                     ier_num = - 13 
                     ier_typ = ER_MAC 
                  ENDIF 
!                                                                       
!     ----continues a macro 'continue'                                  
!                                                                       
              ELSEIF (str_comp (befehl, 'continue', 3, lbef, 8) ) THEN 
                  CALL macro_continue (zeile, lp) 
!                                                                       
!     ----Echo a string, just for interactive check in a macro 'echo'   
!                                                                       
              ELSEIF (str_comp (befehl, 'echo', 2, lbef, 4) ) THEN 
                  CALL echo (zeile, lp) 
!                                                                       
!      ---Evaluate an expression, just for interactive check 'eval'     
!                                                                       
              ELSEIF (str_comp (befehl, 'eval', 2, lbef, 4) ) THEN 
                  CALL do_eval (zeile, lp, .TRUE.) 
!                                                                       
!     ----exit 'exit'                                                   
!                                                                       
              ELSEIF (str_comp (befehl, 'exit', 2, lbef, 4) ) THEN 
                  lend = .true. 
                  EXIT main_loop
!                                                                       
!     ----help 'help' , '?'                                             
!                                                                       
              ELSEIF (str_comp (befehl, 'help', 1, lbef, 4) .OR.  &
                       str_comp (befehl, '?   ', 1, lbef, 4) ) THEN                                      
                  IF (str_comp (zeile, 'errors', 2, lp, 6) ) THEN 
                     lp = lp + 7 
                     CALL do_hel ('discus '//zeile, lp) 
                  ELSE 
                     lp = lp + 12 
                     CALL do_hel ('discus demo '//zeile, lp) 
                  ENDIF 
!                                                                       
!------- -Operating System Kommandos 'syst'                             
!                                                                       
              ELSEIF (str_comp (befehl, 'syst', 2, lbef, 4) ) THEN 
                  IF (zeile.ne.' '.and.zeile.ne.char (13) ) THEN 
                     CALL do_operating (zeile (1:lp), lp) 
                  ELSE 
                     ier_num = - 6 
                     ier_typ = ER_COMM 
                  ENDIF 
!                                                                       
!------  -----waiting for user input                                    
!                                                                       
              ELSEIF (str_comp (befehl, 'wait', 3, lbef, 4) ) THEN 
                  CALL do_input (zeile, lp) 
!
              ELSEIF (str_comp (befehl, 'reset', 3, lbef, 4)) THEN
!                 CALL demo_reset
              ELSE    ! macro, reset or all other commands
!
!opara  =  (/ '0.0000', '0.0000', '0.0000', 'inside' /)   ! Always provide fresh default values
!lopara =  (/  6,        6,        6      ,  6       /)
!owerte =  (/  0.0,      0.0,      0.0    ,  0.0     /)
!
!CALL get_optional(ianz, MAXW, cpara, lpara, NOPTIONAL,  ncalc, &
!                  oname, loname, opara, lopara, lpresent, owerte)
                 is_com: IF(str_comp (befehl, 'incl', 3, lbef, 4) ) THEN 
                    CALL demol_incl(zeile, lp)
                 ELSEIF(str_comp (befehl, 'run', 3, lbef, 3) ) THEN is_com
                    CALL demol_run
                 ELSEIF(str_comp (befehl, 'sel', 3, lbef, 3) .OR.  &
                        str_comp (befehl, 'des', 3, lbef, 3)     ) THEN is_com
!
!--Select
!
                    sel = str_comp (befehl, 'sel', 3, lbef, 3)
                    CALL demol_sel(zeile, lp, sel)
                 ELSEIF(str_comp (befehl, 'show', 3, lbef, 4) ) THEN is_com
!
!--show
!
                    CALL demol_show
!
              ELSE is_com
                 ier_num = -8
                 ier_typ = ER_COMM
!
              ENDIF is_com ! END IF BLOCK actual commands
              ENDIF
           ENDIF is_math   ! END IF BLOCK math equation or specific command
        ENDIF no_com       ! END IF BLOCK no comment
      ENDIF no_err         ! END IF BLOCK no error reading input
!
      IF (ier_num.ne.0) THEN 
         CALL errlist 
         IF (ier_sta.ne.ER_S_LIVE) THEN 
            IF (lmakro .OR. lmakro_error) THEN  ! Error within macro or termination errror
               IF(sprompt /= prompt ) THEN
                  ier_num = -10
                  ier_typ = ER_COMM
                  ier_msg(1) = ' Error occured in demolec menu'
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
!
!
   ENDDO main_loop     ! END DO main loop of menu 
!
   prompt = orig_prompt
!
END SUBROUTINE demolecularize
!
!*******************************************************************************
!
SUBROUTINE demol_incl(zeile, lp)
!
! Select molecules
!
USE get_params_mod
USE ber_params_mod
USE errlist_mod
!
IMPLICIT NONE
!
CHARACTER(LEN=*), INTENT(INOUT) :: zeile
INTEGER, INTENT(INOUT)          :: lp
!
INTEGER, PARAMETER                   :: MAXW = 2
CHARACTER(LEN=1024), DIMENSION(MAXW) :: cpara
INTEGER            , DIMENSION(MAXW) :: lpara
REAL               , DIMENSION(MAXW) :: werte
INTEGER                              :: ianz
!
LOGICAL, EXTERNAL :: str_comp
!
CALL get_params(zeile, ianz, cpara, lpara, MAXW, lp)
IF(ier_num==0) THEN
   IF(ianz==1) THEN
      IF(str_comp(cpara(1), 'all', 1, lpara(1), 3) ) THEN
         dem_incl(1) =  1
         dem_incl(2) = -1
      ELSE
         CALL ber_params(ianz, cpara, lpara, werte, MAXW)
         IF(ier_num==0) THEN
            dem_incl(1) = NINT(werte(1))
            dem_incl(2) = NINT(werte(1))
         ENDIF
      ENDIF
   ELSEIF(ianz==2) THEN
      CALL ber_params(ianz, cpara, lpara, werte, MAXW)
      IF(ier_num==0) THEN
         dem_incl(1) = NINT(werte(1))
         dem_incl(2) = NINT(werte(2))
      ENDIF
   ELSE
      ier_num = - 6
      ier_typ = ER_COMM
   ENDIF
ENDIF
!
END SUBROUTINE demol_incl
!
!*******************************************************************************
!
SUBROUTINE demol_sel(zeile, lp, sel)
!
! Select molecules
!
USE get_params_mod
USE take_param_mod
!
IMPLICIT NONE
!
CHARACTER(LEN=*), INTENT(INOUT) :: zeile
INTEGER, INTENT(INOUT)          :: lp
LOGICAL, INTENT(IN)             :: sel    ! =true for select
!
INTEGER, PARAMETER                   :: MAXW = 4
CHARACTER(LEN=1024), DIMENSION(MAXW) :: cpara
INTEGER            , DIMENSION(MAXW) :: lpara
REAL               , DIMENSION(MAXW) :: werte
INTEGER                              :: ianz
!
INTEGER, PARAMETER :: NOPTIONAL = 4
INTEGER, PARAMETER :: O_MOLTYPE = 1
INTEGER, PARAMETER :: O_CY      = 2
INTEGER, PARAMETER :: O_Cz      = 3
INTEGER, PARAMETER :: O_KEEP    = 4
CHARACTER(LEN=1024), DIMENSION(NOPTIONAL) :: oname   !Optional parameter names
CHARACTER(LEN=1024), DIMENSION(NOPTIONAL) :: opara   !Optional parameter strings returned
INTEGER            , DIMENSION(NOPTIONAL) :: loname  !Lenght opt. para name
INTEGER            , DIMENSION(NOPTIONAL) :: lopara  !Lenght opt. para name returned
LOGICAL            , DIMENSION(NOPTIONAL) :: lpresent!opt. para is present
REAL               , DIMENSION(NOPTIONAL) :: owerte   ! Calculated values
INTEGER, PARAMETER                        :: ncalc = 0 ! Number of values to calculate 
!
LOGICAL, EXTERNAL :: str_comp
!
DATA oname  / 'moletype', 'cy',  'cz',  'keep'   /
DATA loname /  8,    2,     2  ,   4       /
opara  =  (/ '0.0000', '0.0000', '0.0000', 'inside' /)   ! Always provide fresh default values
lopara =  (/  6,        6,        6      ,  6       /)
owerte =  (/  0.0,      0.0,      0.0    ,  0.0     /)
!
CALL get_params(zeile, ianz, cpara, lpara, MAXW, lp)
!
CALL get_optional(ianz, MAXW, cpara, lpara, NOPTIONAL,  ncalc, &
                  oname, loname, opara, lopara, lpresent, owerte)
IF(lpresent(O_MOLTYPE)) THEN                ! User provided molecule types
   IF(str_comp (opara(O_MOLTYPE), 'all', 3, lopara(O_MOLTYPE), 3) ) THEN
      dem_lmoletype(:) = sel
   ENDIF
ENDIF
!
END SUBROUTINE demol_sel
!
!*******************************************************************************
!
SUBROUTINE demol_show
!
USE molecule_mod
USE prompt_mod
!
IMPLICIT NONE
!
INTEGER, DIMENSION(1:DEM_MAX_MOLETYPE) :: moletype_list
INTEGER :: j,k
INTEGER :: i
!
j = 0
DO i=1, mole_num_type
  IF(dem_lmoletype(i)) THEN
     j = j+1
     moletype_list(j) = i
  ENDIF
  WRITE(output_io, 3300) (moletype_list(k),k=1,j)
ENDDO
IF(dem_incl(2)==-1) THEN
  WRITE(output_io, 3400)
ELSE
  WRITE(output_io, 3500) dem_incl(:)
ENDIF
!
3300 FORMAT( '   Selected molecule types   : ',2x,50(i4,1x))
3400 FORMAT( '   Range of incl. molecules  : ',                    &
     &          '  All molecules included')
3500 FORMAT( '   Range of incl. molecules  : ',i9,' to ',i9)
!
END SUBROUTINE demol_show
!
!*******************************************************************************
!
!
SUBROUTINE demol_run
!
! Performs the demolecularization
!
USE crystal_mod
USE molecule_mod
USE prop_para_mod
!
USE errlist_mod
!
IMPLICIT NONE
!
INTEGER :: i, j, k, mmol
INTEGER :: istart, iend
INTEGER :: nat                  ! >0 if there is a true atom in a molecule
LOGICAL :: isdel = .FALSE.
LOGICAL :: empty = .FALSE.
!
isdel = .FALSE.
istart = dem_incl(1)                      ! Cope include range into working variables
IF(dem_incl(2)==-1) THEN
   iend = mole_num_mole
ELSE
   iend   = dem_incl(2)
ENDIF
!
write(*,'(20i3)') mole_cont(:mole_num_atom)
IF(istart>=1    .AND. istart <= mole_num_mole .AND. &    ! Test include range
   istart<=iend .AND. iend   <= mole_num_mole       ) THEN
   DO i=istart, iend                    ! Loop over all included molecules
      IF(dem_lmoletype(mole_type(i))) THEN                 ! Moletype is selected
         DO j=1, mole_len(i)                               ! Loop over all atoms in molecule
            k = mole_cont(mole_off(i)+j)                   ! Absolute atom number
            cr_mole(k) = 0
            cr_prop (k ) = IBCLR (cr_prop (k ), PROP_MOLECULE)   ! UNFLAG THIS ATOM AS SURFACE ANCHOR
            mole_cont(mole_off(i)+j) = 0                   ! Clear molecule content
            isdel = .TRUE.
         ENDDO
      ENDIF
   ENDDO
ELSE
   ier_num = -63
   ier_typ = ER_APPL
   ier_msg(1) = 'Molecule include range is outside range'
ENDIF
!
write(*,*) 'MOLE_CONT  ', LBOUND(mole_cont), UBOUND(mole_cont), mole_num_mole, mole_num_type, mole_num_atom
write(*,'(20i3)') mole_cont(:mole_num_atom)
IF(isdel) THEN                     ! A molecule was deleted
   i = 1
   clean: DO 
      IF(i > mole_num_mole) EXIT clean 
      empty = .TRUE.
      isempty: DO j=1, mole_len(i)                       ! Loop over all atoms in molecule
         IF(mole_cont(mole_off(i)+j) > 0) THEN           ! Absolute atom number
            empty = .FALSE.
            EXIT isempty
         ENDIF
      ENDDO isempty
      IF(empty) THEN                                     ! Molecule is empty
write(*,*) 'molecule is empty ',i
         nat = mole_len(i)
         DO j=mole_off(i)+1, mole_num_atom - mole_len(i)                             ! Shift atoms down
            mole_cont(j) = mole_cont(j+mole_len(i))
            k = mole_cont(j)
            cr_mole(k) = cr_mole(k) - 1
         ENDDO
         DO mmol = i, mole_num_mole-1                    ! Shift molecule properties down
            mole_len (mmol) = mole_len (mmol + 1)
            mole_off (mmol) = mole_off (mmol + 1) - nat
            mole_type(mmol) = mole_type(mmol + 1)
            mole_file(mmol) = mole_file(mmol + 1)
            mole_char(mmol) = mole_char(mmol + 1)
         ENDDO
         mole_num_mole = mole_num_mole - 1
         mole_num_atom = mole_num_atom - nat
      ELSE
write(*,*) 'molecule is full  ',i
         i = i + 1
      ENDIF
   ENDDO clean
ENDIF
write(*,'(a,20i3)') 'CONTENT ', mole_cont(:mole_num_atom)
write(*,'(a,20i3)') 'OFFSET  ', mole_off (:mole_num_mole)
write(*,'(a,20i3)') 'LENGTH  ', mole_len (:mole_num_mole)
!
END SUBROUTINE demol_run
!
!*******************************************************************************
!
END MODULE demolec
