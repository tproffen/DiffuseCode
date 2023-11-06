MODULE domain_menu
!
CONTAINS
!
!*****7*****************************************************************
!                                                                       
SUBROUTINE do_domain
!-                                                                      
!     This subroutine reads a domain from a file and replaces the 
!     corresponding part of the structure.                              
!     This is the new version of the microdomain level.                  
!-                                                                      
USE discus_config_mod 
USE discus_allocate_appl_mod
USE crystal_mod 
USE domain_mod 
USE micro_mod 
USE domaindis_mod 
USE discus_show_menu
!
USE ber_params_mod
USE calc_expr_mod
USE doact_mod 
USe do_eval_mod
USE do_wait_mod
USE build_name_mod
USE errlist_mod 
USE get_params_mod
USE learn_mod 
USE lib_do_operating_mod
USE lib_echo
USE lib_errlist_func
USE lib_help
USE lib_length
USE lib_macro_func
USE class_macro_internal 
USE precision_mod
USE prompt_mod 
USE str_comp_mod
USE string_convert_mod
USE sup_mod
!                                                                       
IMPLICIT none 
!                                                                       
INTEGER, PARAMETER :: MAXW = 9 
!                                                                       
CHARACTER(len=PREC_STRING) :: line    ! Command string
integer                    :: lp      ! length of command string
!
CHARACTER(len=11) :: befehl            ! Command verb
CHARACTER(LEN=LEN(prompt)) :: orig_prompt  ! Prompt prior to this menu
CHARACTER(LEN=MAX(PREC_STRING,LEN(LINE))) :: zeile        ! Dummy string
CHARACTER(LEN=MAX(PREC_STRING,LEN(LINE))), dimension(MAXW) :: cpara    ! Dummy strings for parameters
INTEGER, dimension(MAXW) :: lpara                         ! length of parameter strings
INTEGER :: j, ianz, laenge, lbef 
INTEGER :: indxg 
INTEGER, SAVE    :: n_clu = 0      ! current number of clusters
LOGICAL, SAVE    :: linit = .true. ! do we need to initialize?
LOGICAL          :: lend           ! Finish menu yes/no
REAL(KIND=PREC_DP), dimension(MAXW) :: werte
!                                                                       
!                                                                       
lend = .false. 
!
!     Allocate the arrays in domain_mod.f90. As DISCUS cannot know the 
!     required number, of cluster types, we will just allocate a reasonable
!     amount and increase if necessary.
!                                                                       
!     BEGIN OF new domain code                                          
!                                                                       
orig_prompt = prompt
prompt = prompt (1:len_str (prompt) ) //'/domain' 
!                                                                       
loop_main: DO while (.not.lend)         ! Main "DOMAIN" loop
   CALL no_error 
   CALL get_cmd (line, laenge, befehl, lbef, zeile, lp, prompt) 
   if_cmd: IF (ier_num.eq.0) then 
      IF(line /= ' '      .and. line(1:1) /= '#' .and. &
         line /= char(13) .and. line(1:1) /= '!'        ) THEN
!                                                                       
!     ----search for "="                                                
!                                                                       
         indxg = index (line, '=') 
         IF (indxg.ne.0.AND..NOT. (str_comp (befehl, 'echo',   2, lbef, 4) ) &
                       .AND..NOT. (str_comp (befehl, 'system', 2, lbef, 6) )    &
                       .AND..NOT. (str_comp (befehl, 'help',   2, lbef, 4) .OR. &
                                   str_comp (befehl, '?   ',   2, lbef, 4) )    &
                       .AND. INDEX(line,'==') == 0                            ) THEN
!                                                                       
! ------evaluate an expression and assign the value to a variabble      
!                                                                       
            CALL do_math (line, indxg, laenge) 
         ELSE 
!                                                                       
!------ ----execute a macro file                                        
!                                                                       
            IF(befehl (1:1) .eq.'@') then 
               IF(laenge.ge.2) then 
                  line(1:laenge-1) = line(2:laenge)
                  line(laenge:laenge) = ' '
                  laenge = laenge - 1
                  CALL file_kdo(line, laenge)
               ELSE 
                  ier_num = - 13 
                  ier_typ = ER_MAC 
               ENDIF 
!                                                                       
!------ ----list asymmetric unit 'asym'                                 
!                                                                       
            ELSEIF (str_comp (befehl, 'asym', 2, lbef, 4) ) then 
               CALL show_asym 
!                                                                       
!------ ----list atoms present in the crystal 'chem'                    
!                                                                       
            ELSEIF (str_comp (befehl, 'chemistry', 2, lbef, 9) ) then 
               CALL show_chem 
!                                                                       
!     ----continues a macro 'continue'                                  
!                                                                       
            ELSEIF (str_comp (befehl, 'continue', 5, lbef, 8) ) then 
               CALL macro_continue (zeile, lp) 
!                                                                       
!     ----Echo a string, just for interactive check in a macro 'echo'   
!                                                                       
            ELSEIF (str_comp (befehl, 'echo', 2, lbef, 4) ) then 
               CALL echo (zeile, lp) 
!                                                                       
!      ---Evaluate an expression, just for interactive check 'eval'     
!                                                                       
            ELSEIF (str_comp (befehl, 'evaluate', 2, lbef, 8) ) then 
               CALL do_eval (zeile, lp, .TRUE.) 
!                                                                       
!     ----exit 'exit'                                                   
!                                                                       
            ELSEIF(str_comp (befehl, 'exit', 2, lbef, 4) ) then 
               lend = .true. 
!                                                                       
!     ----help 'help' , '?'                                             
!                                                                       
            ELSEIF(str_comp(befehl, 'help', 1, lbef, 4) .or.                    &
                   str_comp(befehl, '?   ', 1, lbef, 4)      ) then                                      
               IF (str_comp (zeile, 'errors', 2, lp, 6) ) then 
                  lp = lp + 7 
                  CALL do_hel ('discus '//zeile, lp) 
               ELSE 
                  lp = lp + 13 
                  CALL do_hel ('discus domain '//zeile, lp) 
               ENDIF 
!                                                                       
!------- -Operating System Kommandos 'syst'                             
!                                                                       
            ELSEIF (str_comp (befehl, 'system', 2, lbef, 6) ) then 
               IF (zeile.ne.' '.and.zeile.ne.char (13) ) then 
                  CALL do_operating (zeile (1:lp), lp) 
               ELSE 
                  ier_num = - 6 
                  ier_typ = ER_COMM 
               ENDIF 
!                                                                       
!------  -----waiting for user input                                    
!                                                                       
            ELSEIF (str_comp (befehl, 'wait', 3, lbef, 4) ) then 
               CALL do_input (zeile, lp) 
!
!     ----Reset the distribution of domains 'rese'                       
!
            ELSEIF (str_comp (befehl, 'reset', 2, lbef, 5) ) then 
               CALL domain_reset 

!                                                                       
!     ----Original domain commands                                      
!                                                                       
            ELSE
!
               IF(linit ) THEN
                  CALL alloc_domain ( clu_increment )
                  MK_MAX_SCAT = MAX(MK_MAX_SCAT, MAXSCAT)
                  MK_MAX_ATOM = MAX(MK_MAX_ATOM, NMAX)
                  CALL alloc_micro  ( MK_MAX_SCAT , MK_MAX_ATOM)
                  n_clu = CLU_MAX_TYPE
                  linit = .false.
               ENDIF
!                                                                       
!     assign properties to pseudoatoms 'assign'                         
!                                                                       
               IF (str_comp (befehl, 'assign', 1, lbef, 6) ) then 
                  call domain_assign(zeile, lp)
!                                                                       
!     ----set orientation for the current type 'orientation'            
!                                                                       
               ELSEIF (str_comp (befehl, 'orientation', 1, lbef, 11) ) then 
                  CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
                  IF (ier_num.eq.0) then 
                     CALL ber_params (ianz, cpara, lpara, werte, maxw) 
                     IF (ier_num.eq.0) then 
                        j = nint (werte (1) ) 
                        IF (0.lt.j.and.j.le.md_ori_n) then 
                           mv_orient = j 
                        ELSE 
                           ier_num = - 18 
                           ier_typ = ER_APPL 
                        ENDIF 
                     ELSE 
                        ier_num = - 6 
                        ier_typ = ER_COMM 
                     ENDIF 
                  ELSE 
                     ier_num = - 6 
                     ier_typ = ER_COMM 
                  ENDIF 
!                                                                       
!     define name of input file 'inputfile'                             
!                                                                       
               ELSEIF(str_comp(befehl, 'inputfile', 1, lbef, 9) ) then
                  CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
                  IF (ier_num.eq.0) then 
                     CALL do_build_name (ianz, cpara, lpara, werte,     &
                     maxw, 1)                                           
                     IF (ier_num.eq.0) then 
                        clu_infile = cpara (1) (1:lpara(1))
                        clu_infile_internal = clu_infile(1:8)=='internal'
                     ENDIF 
                  ENDIF 
!                                                                       
!     define interpretation of input file 'mode'                        
!                                                                       
               ELSEIF (str_comp (befehl, 'mode', 1, lbef, 4) ) then 
                  CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
                  IF (ier_num.eq.0) then 
                     IF(str_comp(cpara(1), 'domain', 1, lbef, 7)) then
                        clu_mode = CLU_IN_CLUSTER 
                     ELSEIF(str_comp(cpara(1), 'pseudo', 1, lbef, 6)) then
                        clu_mode = CLU_IN_PSEUDO 
                     ELSEIF(str_comp(cpara(1), 'irregular' , 1, lbef, 9)) then
                        clu_mode = CLU_IN_IRREG 
                     ENDIF 
                  ENDIF 
!                                                                       
!     ----Start the distribution of domains 'run'                       
!                                                                       
               ELSEIF (str_comp (befehl, 'run ', 2, lbef, 4) ) then 
                  IF(clu_mode == CLU_IN_CLUSTER) THEN
                     IF(clu_infile_internal) THEN
                        ier_num = -170
                        ier_typ = ER_APPL
                        ier_msg(1) = clu_infile(1:MIN(LEN_TRIM(clu_infile),LEN(ier_msg)))
                     ENDIF
                     clu_number  = 0
                     clu_index   = 0
                     clu_name    = ' '
                     clu_content = ' ' 
                  ENDIF
                  IF(ier_num ==0) THEN
                     call domain_test_expr      ! Test if "EXPR" were used
                     CALL micro_filereading 
                  ENDIF
!                                                                       
!     ----define various surface related settings 'set'                 
!                                                                       
               ELSEIF (str_comp (befehl, 'set', 2, lbef, 3) ) then 
                  CALL domain_do_set (zeile, lp) 
!                                                                       
!     ----show current settings 'show'                                  
!                                                                       
               ELSEIF (str_comp (befehl, 'show', 3, lbef, 4) ) then 
                  CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
                  IF (ier_num.eq.0) then 
                     IF (ianz.eq.0.or.ianz.eq.1) then 
                        CALL show_domain 
                     ELSE 
                        ier_num = -6 
                        ier_typ = ER_COMM 
                     ENDIF 
                  ELSE 
                     ier_num = -6 
                     ier_typ = ER_COMM 
                  ENDIF 
!                                                                       
!     ----unknown command                                               
!                                                                       
               ELSE 
                  ier_num = -8 
                  ier_typ = ER_COMM 
               ENDIF 
            ENDIF 
         ENDIF 
      ENDIF 
   ENDIF if_cmd
!                                                                       
!     --Goto_s and jumps are terrible, yet FORTRTAN does not            
!       have a break                                                    
!     --Jump here if an error occured                                   
!                                                                       
   IF (ier_num.ne.0) THEN 
      CALL errlist 
      IF (ier_sta.ne.ER_S_LIVE) THEN 
         IF (lmakro .OR. lmakro_error) THEN  ! Error within macro or termination errror
            IF(sprompt /= prompt ) THEN
               ier_num = -10
               ier_typ = ER_COMM
               ier_msg(1) = ' Error occured in domain menu'
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
            ier_num = -11 
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
ENDDO loop_main
!
      prompt = orig_prompt
!                                                                       
      END SUBROUTINE do_domain                      
!
!*****7*****************************************************************
!
subroutine domain_assign(zeile, lp)
!-
!  Perform the 'assign' task
!+
!
use domain_mod
use discus_allocate_appl_mod
!
use errlist_mod
use ber_params_mod
use get_params_mod
use build_name_mod
use precision_mod
use str_comp_mod
use string_convert_mod
!
character(len=*), intent(inout) :: zeile
integer         , intent(inout) :: lp
!
integer, parameter :: maxw = 15
!
character(LEN=PREC_STRING), dimension(4   ) :: dummy_expr
character(LEN=PREC_STRING), dimension(MAXW) :: cpara
integer                   , dimension(MAXW) :: lpara
real(kind=PREC_DP)        , dimension(MAXW) :: werte
integer :: i, j    ! Dummy loop indices
integer :: ianz    ! Number parameters
integer :: n_clu   ! Dummy for allocation
!
CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
IF(ier_num /= 0) return
!
IF (ianz.ge.3) then 
   CALL do_build_name(ianz, cpara, lpara, werte, maxw, 2)
   clu_index = 0 
   CALL do_cap (cpara (2) ) 
   DO i = 1, clu_number 
!                       IF (str_comp(cpara (2), clu_name (i), 2, lpara &
!                       (2), 4) ) then                                  
      IF(cpara(2)(1:4)==clu_name(i)(1:4)) THEN
         clu_index = i 
      ENDIF 
   ENDDO 
   IF (clu_index.eq.0) then 
!
!                          If necessary increase array sizes
!
      IF ( clu_number >= CLU_MAX_TYPE ) THEN
         i     = CLU_MAX_TYPE + clu_increment
         CALL alloc_domain ( i )
         n_clu = CLU_MAX_TYPE
      ENDIF
      clu_number = clu_number + 1 
      clu_index = clu_number 
      clu_name (clu_index) = ' ' 
      clu_name (clu_index) = cpara (2) (1:lpara(2))
   ENDIF 
   IF(str_comp(cpara(1), 'character', 2, lpara(1), 9) ) then                                 
      IF(str_comp(cpara(3), 'cube', 2, lpara(3), 4) ) then
         clu_character(clu_index) = CLU_CHAR_CUBE
      ELSEIF(str_comp(cpara(3), 'cylinder', 2, lpara(3), 8)) then                        
         clu_character(clu_index) = CLU_CHAR_CYLINDER
      ELSEIF(str_comp(cpara(3), 'fuzzy', 1, lpara(3), 5)) then
         clu_character(clu_index) = CLU_CHAR_FUZZY
      ELSEIF(str_comp(cpara(3), 'irregular', 1, lpara(3), 9)) then
         clu_character(clu_index) = CLU_CHAR_IRREG
      ELSEIF(str_comp(cpara(3), 'sphere', 2, lpara(3), 6)) then
         clu_character(clu_index) = CLU_CHAR_SPHERE
      ENDIF 
   ELSEIF(str_comp(cpara(1), 'contentfile', 2, lpara(1), 11)) then
      CALL del_params (2, ianz, cpara, lpara, maxw) 
      IF (ier_num.eq.0) then 
         CALL do_build_name(ianz, cpara, lpara, werte, maxw, 1)
         IF (ier_num.eq.0) then
            clu_content(clu_index) = cpara(1)(1:lpara(1))
         ENDIF 
      ENDIF 
   ELSEIF(str_comp(cpara(1), 'fuzzy', 1, lpara(1), 5) ) then
      CALL del_params (2, ianz, cpara, lpara, maxw) 
      IF (ier_num.eq.0) then 
         CALL ber_params (ianz, cpara, lpara, werte, maxw)
         IF (ier_num.eq.0) then 
            clu_fuzzy(clu_index) = werte (1) 
         ENDIF 
      ENDIF 
   ELSEIF(str_comp(cpara(1), 'orientation', 1, lpara(1), 11) ) then   ! Orientation of guest within host
      CALL del_params(2, ianz, cpara, lpara, maxw) 
      IF (ier_num.eq.0) then 
         werte(:) = 0.0
         do j = 1, 4                              ! Test if "EXPR" occurs in a parameter
            if(index(cpara(j+1),'EXPR')>0) then
               dummy_expr(j) = cpara(j+1)         ! Temporary copy
            else
               dummy_expr(j) = ' '                ! No occurence, nullify temporary
            endif
         enddo
         CALL ber_params(ianz, cpara, lpara, werte, maxw)
         IF (ier_num.eq.0) then 
            i = nint (werte (1) ) 
            DO j = 1, 4 
               clu_orient(clu_index, i, j) = werte(j + 1)        ! Copy result of berechne
               if(dummy_expr(j) == ' ') then          ! Dummy is empty use current value
                  write(clu_or_expr(clu_index, i,j), '(g16.8e3)') werte(j+1)
               else
                  clu_or_expr(clu_index, i,j) = dummy_expr(j)    ! Copy original EXPR string
               endif
!write(*,*) ' CLU_OR_EXPR ', clu_index, clu_or_expr(clu_index, i,j) &
!                           (1:len_trim(clu_or_expr(clu_index, i,j)))
            ENDDO 
         ENDIF 
      ENDIF 
   ELSEIF(str_comp(cpara(1), 'shape', 1, lpara(1) , 5) ) then                                 
      CALL del_params(2, ianz, cpara, lpara, maxw) 
      IF (ier_num.eq.0) then 
         werte(:) = 0.0
         CALL ber_params(ianz, cpara, lpara, werte, maxw)                              
         IF (ier_num.eq.0) then 
            i = nint (werte (1) ) 
            DO j = 1, 4 
               clu_shape(clu_index,i,j) = werte(j+1)                                   
            ENDDO 
            clu_sigma(clu_index,i) = ABS(werte(6))
         ENDIF 
      ENDIF 
   ENDIF 
ELSE 
   ier_num = - 6 
   ier_typ = ER_COMM 
endif 
!
end subroutine domain_assign
!
!*****7*****************************************************************
!
subroutine domain_test_expr
!-
!  Test if the orientations and/or the shapes are subject to "EXPR" calculation mode
!+
!
use domain_mod
!
integer :: i,j, k
!
clu_use_expr = .FALSE.                 ! Initialize all
!
loop_orient: do k=1, clu_number
   do i=1, 3
      do j=1, 4
         if(index(clu_or_expr(k, i,j),'EXPR') > 0) then
            clu_use_expr(k,1) = .TRUE.
            cycle loop_orient
         endif
      enddo
   enddo
enddo loop_orient
!do k=1, clu_number
!write(*,*) ' CLUSTER TYPE ', k, clu_use_expr(k, :)
!enddo
!
end subroutine domain_test_expr
!
!*****7*****************************************************************
!
      SUBROUTINE show_domain 
!-                                                                      
!     Displays the current parameters for domain distributions          
!+                                                                      
      USE discus_config_mod 
      USE discus_allocate_appl_mod 
      USE crystal_mod 
      USE domain_mod 
      USE surface_mod 
      USE errlist_mod 
USE lib_length
      USE prompt_mod 
      IMPLICIT none 
!                                                                       
       
!                                                                       
      CHARACTER(9) char_type ( -5:0) 
      CHARACTER(12) char_mode (0:2) 
      INTEGER i, ii, j 
      INTEGER              :: new_nscat ! DUMMY for allocation
!                                                                       
!                                                                       
DATA char_mode / 'domain     ', 'pseudoatoms', 'irregular' / 
DATA char_type / 'irreg    ', 'fuzzy    ', 'sphere   ', 'cylinder ',            &
                 'cube     ', 'undefined'                             /
!                                                                       
      IF (cr_nscat == MAXSCAT) then 
         new_nscat = MAX(MAXSCAT + 5, INT ( MAXSCAT * 1.025 ))
         CALL alloc_crystal_scat(new_nscat)
         CALL alloc_crystal_nmax(NMAX)
         CALL alloc_surf   (new_nscat)
      ENDIF 
      IF (MAXSCAT  > ubound(surf_in_dist,1)) THEN
         new_nscat = MAXSCAT
         CALL alloc_surf   (new_nscat)
      ENDIF 
      i = 0 
      i = len_str (clu_infile) 
!                                                                       
      WRITE (output_io, 3000) 
      WRITE (output_io, 3100) char_mode (clu_mode) 
      IF (i.gt.0) then 
         WRITE (output_io, 3200) clu_infile (1:i) 
      ELSE 
         WRITE (output_io, 3250) 
      ENDIF 
!                                                                       
      IF (clu_mode.eq.CLU_IN_PSEUDO) then 
         WRITE (output_io, 4000) 
         DO i = 1, clu_number 
         WRITE (output_io, 4100) clu_name (i) 
         WRITE (output_io, 4200) char_type (clu_character (i) ) 
         j = 0 
         j = len_str (clu_content (i) ) 
         IF (j.gt.0) then 
            WRITE (output_io, 4300) clu_content (i) (1:j) 
         ELSE 
            WRITE (output_io, 4350) 
         ENDIF 
         IF( clu_fuzzy(i) > 0 ) THEN
            WRITE (output_io, 4400) clu_fuzzy (i) 
         ELSE
            WRITE (output_io, 4410) 
         ENDIF
         DO ii = 1, 3 
         WRITE (output_io, 4500) (clu_shape (i, ii, j), j = 1, 4), &
              clu_sigma(i,ii)
         ENDDO 
         WRITE (output_io, * ) 
         DO ii = 1, 3 
         WRITE (output_io, 4600) (clu_orient (i, ii, j), j = 1, 4) 
         ENDDO 
         ENDDO 
      ENDIF 
      IF (SURF_MAXSCAT==0) THEN
        WRITE(output_io,*) ' No distances to surfaces have been defined yet'
        WRITE(output_io,*) ' Set distances within surface menu first'
        RETURN
      ENDIF
      WRITE (output_io, * ) 
      WRITE (output_io, 5100) 
      WRITE (output_io, * ) 
      DO i = 0, cr_nscat 
      WRITE (output_io, 5010) cr_at_lis (i), surf_in_dist (i) 
      ENDDO 
      IF (cr_nscat.lt.SURF_MAXSCAT) then 
         WRITE (output_io, 5010) 'new', surf_in_dist (SURF_MAXSCAT) 
      ENDIF 
!                                                                       
      WRITE (output_io, * ) 
!                                                                       
 3000 FORMAT(/'****************************************************'//  &
     &        ' Characterisation of domain distribution'/)              
 3100 FORMAT(' Inputfile contains     : ',a11) 
 3200 FORMAT(' Inputfile name         : ',a) 
 3250 FORMAT(' Inputfile name         : ','undefined') 
 4000 FORMAT(' Pseudoatom definitions : ',a) 
 4100 FORMAT(/,'     Pseudoatom name    : ',a4) 
 4200 FORMAT('     domain type        : ',a9) 
 4300 FORMAT('     content in file    : ',a) 
 4350 FORMAT('     content in file    : ','undefined') 
 4400 FORMAT('     fuzzy separation   : ',f15.7,' A') 
 4410 FORMAT('     fuzzy separation   : ','switched OFF old atoms are not removed')
 4500 FORMAT('     shape orientation  : ','(',3(f9.4,2x),2x,f9.4,') +- ',f9.4 ) 
 4600 FORMAT('     atom orientation   : ','(',3(f9.4,2x),2x,f9.4,')') 
 5100 FORMAT (' Distances between atom types and an internal surface') 
 5010 FORMAT (' Atom type   ',a4,' at ',f8.3,' Angstroem') 
      END SUBROUTINE show_domain                    
!*****7*****************************************************************
!
SUBROUTINE micro_filereading 
!
!  Main execution routine for the domains
!
USE discus_allocate_appl_mod
USE discus_config_mod 
USE discus_allocate_appl_mod
USE chem_mod
USE crystal_mod 
USE domain_mod 
USE domaindis_mod 
USE micro_mod 
USE molecule_mod
use prop_para_func
use prop_para_mod
USE read_internal_mod
USE discus_save_mod 
USE save_menu, ONLY: save_internal, save_store_setting, save_restore_setting, save_default_setting, save_struc, save_show
USE structur, ONLY: stru_readheader, test_file, read_to_internal
!
USE errlist_mod 
use matrix_mod
USE precision_mod
USE prompt_mod
USE support_mod
USE times_mod
!                                                                       
IMPLICIT none 
!                                                                       
INTEGER, PARAMETER :: IMD = 45 
integer, parameter :: MAXMASK = 4
!                                                                       
CHARACTER(LEN=PREC_STRING)           :: infile 
INTEGER, PARAMETER                   :: AT_MAXP = 16
INTEGER                              :: at_ianz
CHARACTER(LEN=8), DIMENSION(AT_MAXP) :: at_param
character(len=PREC_STRING)           :: host_file
character(len=PREC_STRING)           :: strucfile
character(len=PREC_STRING)           :: line
!
LOGICAL ::lend 
LOGICAL :: l_ok    ! the domain list file contains a correct input pseuodatom
LOGICAL :: lread 
LOGICAL :: lmetric 
REAL(PREC_DP), DIMENSION(4,4) :: mc_dimen
REAL(PREC_DP), DIMENSION(4,4) :: mc_idimen
REAL(PREC_DP), DIMENSION(4,4) :: mc_matrix
INTEGER :: natoms_old
INTEGER               :: i, j, k, l
integer               :: length
INTEGER, DIMENSION(5) :: maxdim = 0
INTEGER  :: natoms  ! Number of atoms in the input file
INTEGER  :: nscats  ! NUmber of atom types in the input file
INTEGER  :: n_gene  ! Number of generators for molecules in the input file
INTEGER  :: n_symm  ! Number of symm. ops. for molecules in the input file
INTEGER  :: n_mole  ! Number of molecules in the input file
INTEGER  :: n_type  ! Number of molecule types  in the input file
INTEGER  :: n_atom  ! Number of atoms within molecules in the input file
integer, dimension(3) :: n_cells  ! Number of unit cells in the input file
integer, dimension(:,:), allocatable :: clu_moles  ! No of molecule and types in each cluster
logical, dimension(0:MAXMASK) :: uni_mask
REAL(PREC_DP)               :: shortest
REAL(PREC_DP), DIMENSION(3) :: vv
REAL(PREC_DP) :: sgrand  ! grand time
!
nscats = 0
natoms = 0
n_mole = 0
n_type = 0
n_atom = 0
maxdim(:) = 0
clu_mole_num_mole  = 0   ! Number of Clusters in case of CLUSTER MODE
clu_mole_num_type  = 0   ! Number of Cluster  types in case of CLUSTER MODE
clu_mole_num_atom  = 0   ! Number of Cluster  atoms in case of CLUSTER MODE
uni_mask(0)   = .true.
uni_mask(1:3) = .true.
uni_mask(4)   = .false.
! Read all content files into internal storage
!
! --- Temporarily save current structure into internal
!
call save_store_setting             ! Backup user "save" setting
call save_default_setting           ! Default to full saving
!
line       = 'ignore, all'          ! Ignore all properties
length     = 11 
call property_select(line, length, sav_sel_prop)
line       = 'ignore, all'          ! Ignore all properties for global as well
length     = 11
call property_select(line, length,  cr_sel_prop)
host_file = 'internal_domain_host.stru'
call save_internal(host_file)
maxdim(1) = max(maxdim(1), cr_natoms, NMAX)
maxdim(2) = MAX(maxdim(2), nscats, MAXSCAT)
maxdim(3) = MAX(maxdim(3), n_mole, MOLE_MAX_MOLE)
maxdim(4) = MAX(maxdim(4), n_type, MOLE_MAX_TYPE)
maxdim(5) = MAX(maxdim(5), n_atom, MOLE_MAX_ATOM)
!
! -- Store list file internally
!
strucfile = clu_infile
!write(*,*) ' CLU_INFILE ', clu_infile(1:len_trim(clu_infile))
call read_to_internal(strucfile, 'internal_domain_list.')
if(ier_num/=0) return
!write(*,*) ' infile : ier_num, ier_typ ', ier_num, ier_typ
if(clu_mode==CLU_IN_CLUSTER) then
   CLU_MOLE_MAX_mole = mole_num_mole
   CLU_MOLE_MAX_type = mole_num_type
   CLU_MOLE_MAX_atom = mole_num_atom
   clu_mole_num_mole = mole_num_mole
   clu_mole_num_type = mole_num_type
   clu_mole_num_atom = mole_num_atom
endif

clu_infile = 'internal_domain_list.' // clu_infile(1:len_trim(clu_infile))
clu_infile_internal = .true.
maxdim(1) = max(maxdim(1), cr_natoms, NMAX)
maxdim(2) = MAX(maxdim(2), nscats, MAXSCAT)
maxdim(3) = MAX(maxdim(3), n_mole, MOLE_MAX_MOLE)
maxdim(4) = MAX(maxdim(4), n_type, MOLE_MAX_TYPE)
maxdim(5) = MAX(maxdim(5), n_atom, MOLE_MAX_ATOM)
MK_MAX_SCAT = MAX(MK_MAX_SCAT, MAXSCAT, nscats)
MK_MAX_ATOM = MAX(MK_MAX_ATOM, NMAX, natoms)
n_mole = MAX(n_mole, MOLE_MAX_MOLE)
n_type = MAX(n_type, MOLE_MAX_TYPE)
n_atom = MAX(n_atom, MOLE_MAX_ATOM)
!
if_clu:if(clu_mode==CLU_IN_CLUSTER) THEN          ! Cluster mode, read content files
   loop_clu_mole: do i=1, mole_num_mole
      loop_clu_srch: do j=1, clu_number
         if(clu_content(j)==mole_file(i)) cycle loop_clu_mole
      enddo loop_clu_srch
      clu_number = clu_number + 1
      clu_content(clu_number) = mole_file(i)
   enddo loop_clu_mole
endif if_clu
!write(*,*) ' READING GUEST FILES '
!
! --- Now read all content files from disk 
!
loop_read: do i=1, clu_number
   if(clu_content(i)(1:8)/='internal') then    ! not yet an internal file
      strucfile = clu_content(i)
      call read_to_internal(strucfile, 'internal_domain_guest.')
      if(ier_num/=0) then
         j = ier_num
         call readstru_internal(MAXMASK, host_file, uni_mask)
         call save_restore_setting           ! Restore user "save" settings
         ier_num = j
         write(ier_msg(1), '(a,i3)') ' Error in content file ', i
         ier_msg(2) = clu_content(i)(1:min(len(ier_msg), len_trim(clu_content(i))))
         return
      endif
      clu_content(i) = 'internal_domain_guest.' // clu_content(i)(1:len_trim(clu_content(i)))
   endif
enddo loop_read
!write(*,*) ' FINISHED  GUEST FILES '
!
call readstru_internal(MAXMASK, host_file, uni_mask)   ! Restore original host file
if(ier_num/=0) return
call save_restore_setting           ! Restore user "save" settings
if(ier_num/=0) return
!-------------------------------------------------------------------------------
!
!     TEST FILESIZE OF ALL internal CONTENT FILES
!
! Prepare list of molecule types with the clusters
!
if(allocated(clu_moles)) deallocate(clu_moles)
allocate(clu_moles(2,0:clu_number))
clu_moles  = 0
clu_moles(1,0) = mole_num_mole
clu_moles(2,0) = mole_num_type
DO i=1, clu_number
   infile = clu_content(i)
   IF(clu_content(i)(1:8)=='internal')THEN
         CALL testfile_internal(  infile , natoms, &
              nscats, n_mole, n_type, n_atom)
      clu_moles(1,i) = n_mole
      clu_moles(2,i) = n_type
      MK_MAX_SCAT = MAX(MK_MAX_SCAT, MAXSCAT, nscats)
      MK_MAX_ATOM = MAX(MK_MAX_ATOM, NMAX, natoms)
      CALL alloc_micro  ( MK_MAX_SCAT , MK_MAX_ATOM)
      IF(n_mole > MOLE_MAX_MOLE .or.  &
         n_type > MOLE_MAX_TYPE .or.  &
         n_atom > MOLE_MAX_ATOM       ) THEN
         n_gene = MAX( 1, MOLE_MAX_GENE)
         n_symm = MAX( 1, MOLE_MAX_SYMM)
         n_mole = MAX(n_mole, MOLE_MAX_MOLE)
         n_type = MAX(n_type, MOLE_MAX_TYPE)
         n_atom = MAX(n_atom, MOLE_MAX_ATOM)
         CALL alloc_molecule(n_gene, n_symm, n_mole, n_type, n_atom)
      ENDIF
         maxdim(1) = MAX(natoms, MK_MAX_ATOM, NMAX)
         maxdim(2) = MAX(nscats, MK_MAX_SCAT, MAXSCAT)
         maxdim(3) = MAX(n_mole, MOLE_MAX_MOLE)
         maxdim(4) = MAX(n_type, MOLE_MAX_TYPE)
         maxdim(5) = MAX(n_atom, MOLE_MAX_ATOM)
         CALL stru_readheader_internal (infile, MK_MAX_SCAT, mk_name,   &
         mk_spcgr, mk_spcgr_set, mk_set, mk_iset,                       &
         mk_at_lis, mk_nscat, mk_dw, mk_occ, mk_a0, mk_win,       &
         sav_ncell, sav_r_ncell, sav_ncatoms, spcgr_ianz, spcgr_para, &
         mk_GEN_ADD_MAX, mk_gen_add_n, mk_gen_add_power, mk_gen_add,  &
         mk_SYM_ADD_MAX, mk_sym_add_n, mk_sym_add_power, mk_sym_add )
         IF(ier_num /= 0) then
            ier_msg(1) = ' Guest file ' // infile(1:min(len_trim(infile),66))
            RETURN
         endif
   ELSE
      CALL test_file(infile ,natoms, &
              nscats, n_mole, n_type, n_atom, n_cells, -1 , .false., MAXMASK, uni_mask)
      clu_moles(1,i) = n_mole
      clu_moles(2,i) = n_type
      IF (ier_num.ne.0) THEN
         RETURN
      ENDIF
      MK_MAX_SCAT = MAX(MK_MAX_SCAT, MAXSCAT, nscats)
      MK_MAX_ATOM = MAX(MK_MAX_ATOM, NMAX, natoms)
      CALL alloc_micro  ( MK_MAX_SCAT , MK_MAX_ATOM)
      IF(n_mole > MOLE_MAX_MOLE .or.  &
         n_type > MOLE_MAX_TYPE .or.  &
         n_atom > MOLE_MAX_ATOM       ) THEN
         n_gene = MAX( 1, MOLE_MAX_GENE)
         n_symm = MAX( 1, MOLE_MAX_SYMM)
         n_mole = MAX(n_mole, MOLE_MAX_MOLE)
         n_type = MAX(n_type, MOLE_MAX_TYPE)
         n_atom = MAX(n_atom, MOLE_MAX_ATOM)
         CALL alloc_molecule(n_gene, n_symm, n_mole, n_type, n_atom)
      ENDIF
      maxdim(1) = MAX(natoms, MK_MAX_ATOM, NMAX)
      maxdim(2) = MAX(nscats, MK_MAX_SCAT, MAXSCAT)
      maxdim(3) = MAX(n_mole, MOLE_MAX_MOLE)
      maxdim(4) = MAX(n_type, MOLE_MAX_TYPE)
      maxdim(5) = MAX(n_atom, MOLE_MAX_ATOM)
      CALL oeffne (imd, infile, 'old') 
      IF (ier_num.ne.0) THEN
         CLOSE(imd)
         RETURN
      ENDIF
!                                                                       
!     Read the input file header                                        
!                                                                       
      CALL stru_readheader (imd, MK_MAX_SCAT, mk_name,     &
      mk_spcgr, mk_set, mk_at_lis, mk_nscat, mk_dw, mk_occ, mk_anis, mk_a0, mk_win, sav_ncell,   &
      sav_r_ncell, sav_ncatoms, mk_spcgr_ianz, mk_spcgr_para, &
      AT_MAXP, at_ianz, at_param)           
      CLOSE(imd)
      IF (ier_num.ne.0) THEN 
         ier_msg(1) = ' Guest file ' // infile(1:min(len_trim(infile),66))
         RETURN
      ENDIF
   ENDIF
!                                                                       
!------ Here we should include comparison of the unit cell              
!       space group etc.                                                
!                                                                       
   lmetric = abs (mk_a0 (1) - cr_a0 (1) ) .lt.0.0001 .AND. &
             abs (mk_a0 (2) - cr_a0 (2) ) .lt.0.0001 .AND. &
             abs (mk_a0 (3) - cr_a0 (3) ) .lt.0.0001 .AND. &
             abs (mk_win (1) - cr_win (1) ) .lt.0.001.AND. &
             abs (mk_win (2) - cr_win (2) ) .lt.0.001.AND. &
             abs (mk_win (3) - cr_win (3) ) .lt.0.001                                                    
   IF (.NOT.lmetric) THEN 
      ier_num = -87 
      ier_typ = ER_APPL 
      ier_msg (1) = 'The lattice constants in the content' 
      ier_msg (2) = 'file differ from those of the host' 
      RETURN 
   ENDIF 
   IF (mk_spcgr (1:1) .ne.cr_spcgr (1:1) ) then 
      ier_num = - 88 
      ier_typ = ER_APPL 
      ier_msg (1) = 'between microdomain content and host' 
      RETURN 
   ENDIF 
ENDDO
if(allocated(clu_mole_tab)) deallocate(clu_mole_tab)
allocate(clu_mole_tab(maxval(clu_moles(2,1:clu_number)),1:clu_number))
clu_mole_tab = 0
j = clu_moles(2,0)             ! molecule types prior
do i=1,clu_number
   do k=1,clu_moles(2,i)
     j = j + 1
     clu_mole_tab(k,i) = j
   enddo
enddo
natoms = max(natoms, maxdim(1))
nscats = max(nscats, maxdim(2))
n_mole = max(n_mole, maxdim(3))
n_type = max(n_type, maxdim(4), maxval(clu_mole_tab))   ! Previous and cluster molecule types
n_atom = max(n_atom, maxdim(5))
   IF(natoms > MAX(MK_MAX_ATOM, NMAX)    .or. &
      nscats > MAX(MK_MAX_SCAT, MAXSCAT) .or. &
      n_mole > MOLE_MAX_MOLE             .or. &
      n_type > MOLE_MAX_TYPE             .or. &
      n_atom > MOLE_MAX_ATOM                   ) THEN
      natoms = MAX(natoms, MK_MAX_ATOM, NMAX)
      nscats = MAX(nscats, MK_MAX_SCAT, MAXSCAT)
      n_mole = MAX(n_mole, MOLE_MAX_MOLE)
      n_type = MAX(n_type, MOLE_MAX_TYPE)
      n_atom = MAX(n_atom, MOLE_MAX_ATOM)
      CALL alloc_micro  ( nscats, natoms)
      MK_MAX_ATOM = natoms
      MK_MAX_SCAT = nscats
   ENDIF
!
CALL alloc_crystal_scat (nscats)
CALL alloc_crystal_nmax (natoms)
!
!
clu_remove_end = cr_natoms    ! Initially remove only atoms in original crystal
!                                                                       
lread = .true. 
IF ( clu_infile_internal ) THEN
!write(*,*) ' CLU_INFILE AAA ', clu_infile(1:len_trim(clu_infile))
   CALL stru_readheader_internal (clu_infile, MK_MAX_SCAT, mk_name,   &
   mk_spcgr, mk_spcgr_set, mk_set, mk_iset,                           &
   mk_at_lis, mk_nscat, mk_dw, mk_occ, mk_a0, mk_win,       &
   sav_ncell, sav_r_ncell, sav_ncatoms, spcgr_ianz, spcgr_para, &
   mk_GEN_ADD_MAX, mk_gen_add_n, mk_gen_add_power, mk_gen_add,  &
   mk_SYM_ADD_MAX, mk_sym_add_n, mk_sym_add_power, mk_sym_add )
   IF(ier_num /= 0) then
      ier_msg(1) = ' Input file ' // clu_infile(1:min(len_trim(clu_infile),66))
      return
   endif
   clu_iatom = 0
   at_param(1) = 'X'
   at_param(2) = 'Y'
   at_param(3) = 'Z'
   at_param(4) = 'BISO'
   at_param(5) = 'PROPERTY'
   at_param(6) = 'MOLENO'
   at_param(7) = 'MOLEAT'
   at_param(8) = 'OCC'
   at_ianz = 8
ELSE
   CALL test_file(clu_infile, natoms, nscats, n_mole, n_type, n_atom, n_cells, -1 , .false., MAXMASK, uni_mask)
   IF(natoms > MAX(MK_MAX_ATOM, NMAX)    .or. &
      nscats > MAX(MK_MAX_SCAT, MAXSCAT) .or. &
      n_mole > MOLE_MAX_MOLE             .or. &
      n_type > MOLE_MAX_TYPE             .or. &
      n_atom > MOLE_MAX_ATOM                   ) THEN   
      natoms = MAX(natoms, MK_MAX_ATOM, NMAX)
      nscats = MAX(nscats, MK_MAX_SCAT, MAXSCAT)
      n_mole = MAX(n_mole, MOLE_MAX_MOLE)
      n_type = MAX(n_type, MOLE_MAX_TYPE)
      n_atom = MAX(n_atom, MOLE_MAX_ATOM)
      CALL alloc_micro  ( nscats, natoms)
      MK_MAX_ATOM = natoms
      MK_MAX_SCAT = nscats
   ENDIF
   CALL oeffne (imd, clu_infile, 'old') 
   IF (ier_num.ne.0) THEN 
      CLOSE(imd)
      if(allocated(clu_moles)) deallocate(clu_moles)
      RETURN
   ENDIF
!                                                                       
!     Read the input file header                                        
!                                                                       
   CALL stru_readheader (imd, MK_MAX_SCAT, mk_name,     &
   mk_spcgr, mk_set, mk_at_lis, mk_nscat, mk_dw, mk_occ, mk_anis, mk_a0, mk_win, sav_ncell,   &
   sav_r_ncell, sav_ncatoms, mk_spcgr_ianz, mk_spcgr_para, &
   AT_MAXP, at_ianz, at_param)           
ENDIF
IF (ier_num.ne.0) THEN 
      ier_msg(1) = ' Input file ' // clu_infile(1:min(len_trim(clu_infile),66))
   CLOSE(imd)
      if(allocated(clu_moles)) deallocate(clu_moles)
   RETURN
ENDIF
!                                                                       
!------ Here we should include comparison of the unit cell              
!       space group etc.                                                
!                                                                       
lmetric = abs (mk_a0 (1) - cr_a0 (1) ) .lt.0.0001 .AND. &
          abs (mk_a0 (2) - cr_a0 (2) ) .lt.0.0001 .AND. &
          abs (mk_a0 (3) - cr_a0 (3) ) .lt.0.0001 .AND. &
          abs (mk_win (1) - cr_win (1) ) .lt.0.001.AND. &
          abs (mk_win (2) - cr_win (2) ) .lt.0.001.AND. &
          abs (mk_win (3) - cr_win (3) ) .lt.0.001                                                    
IF (.NOT.lmetric) THEN 
   ier_num = -87 
   ier_typ = ER_APPL 
   ier_msg (1) = 'The lattice constants in the domain' 
   ier_msg (2) = 'file differ from those of the host' 
   CLOSE (imd) 
      if(allocated(clu_moles)) deallocate(clu_moles)
   RETURN 
ENDIF 
IF (mk_spcgr (1:1) .ne.cr_spcgr (1:1) ) then 
   ier_num = - 88 
   ier_typ = ER_APPL 
   ier_msg (1) = 'between microdomain and host' 
   CLOSE (imd) 
      if(allocated(clu_moles)) deallocate(clu_moles)
   RETURN 
ENDIF 
!                                                                       
sgrand = seknds(0.0)
!satoms = 0.0
!shead  = 0.0
!srem   = 0.0
!sstep  = 0.0
IF (clu_mode.eq.CLU_IN_PSEUDO) then               ! Pseudo atom mode
   infile = ' '
   CALL micro_read_simple(imd, lend, l_ok, infile, mc_dimen, mc_idimen,&
   mc_matrix, MK_MAX_SCAT, mk_at_lis)                                                     
   IF (ier_num.ne.0) THEN 
      CLOSE(imd)
      if(allocated(clu_moles)) deallocate(clu_moles)
      RETURN
   ENDIF
!                                                                       
   read_domains: DO WHILE (.NOT.lend)
      pseudo_ok:  IF(l_ok) THEN
         IF(mc_type .lt.0) then 
            CALL micro_read_atoms(infile, mc_idimen, mc_matrix)
         ENDIF 
         IF (ier_num.ne.0) THEN 
            CLOSE(imd)
      if(allocated(clu_moles)) deallocate(clu_moles)
            RETURN
         ENDIF
!
         IF ( clu_infile_internal ) THEN
            CALL stru_readheader_internal (clu_infile, MK_MAX_SCAT, mk_name, &
            mk_spcgr, mk_spcgr_set, mk_set, mk_iset,                         &
            mk_at_lis, mk_nscat, mk_dw, mk_occ, mk_a0, mk_win,     &
            sav_ncell, sav_r_ncell, sav_ncatoms, spcgr_ianz, spcgr_para, &
            mk_GEN_ADD_MAX, mk_gen_add_n, mk_gen_add_power, mk_gen_add,  &
            mk_SYM_ADD_MAX, mk_sym_add_n, mk_sym_add_power, mk_sym_add )
         ENDIF 
!
      ENDIF  pseudo_ok
      CALL micro_read_simple (imd, lend, l_ok, infile, mc_dimen, mc_idimen,&
            mc_matrix, MK_MAX_SCAT, mk_at_lis)                                                     
      IF(ier_ctrlc) THEN
         CLOSE(imd)
      if(allocated(clu_moles)) deallocate(clu_moles)
         ier_num = -14
         ier_typ = ER_COMM
         RETURN
      ENDIF
      IF (ier_num.ne.0) THEN 
         CLOSE(imd)
      if(allocated(clu_moles)) deallocate(clu_moles)
         RETURN
      ENDIF
   ENDDO  read_domains 
!
!        if the user trusts that no domains overlap among each other, then
!        old atoms are removed only here at the end
!
   IF(clu_remove_mode == CLU_REMOVE_TRUST ) THEN
      mk_dim       = 0.
      natoms_old   = clu_remove_end
      mk_natoms    = cr_natoms - natoms_old 
      md_sep_fuz   = clu_remove_dist
      vv           = 0.0
      shortest     = 0.0
      mc_type      = MD_DOMAIN_FUZZY
      IF(md_sep_fuz > 0.0 ) THEN
!
         CALL micro_fuzzy_rem (mk_dim, natoms_old, mk_natoms, md_sep_fuz,  &
         vv, shortest, mc_type, MD_DOMAIN_FUZZY)                                                     
!
      ENDIF
   ENDIF
!                                                                       
elseif(clu_mode == CLU_IN_IRREG ) then               ! IRREGULAR domain mode
   close(imd)
   call domain_irregular !(imd, lend, infile, mc_dimen, mc_idimen, mc_matrix)
ELSE 
    CALL stru_readheader_internal (clu_infile, MK_MAX_SCAT, mk_name, &
            mk_spcgr, mk_spcgr_set, mk_set, mk_iset,                         &
            mk_at_lis, mk_nscat, mk_dw, mk_occ, mk_a0, mk_win,     &
            sav_ncell, sav_r_ncell, sav_ncatoms, spcgr_ianz, spcgr_para, &
            mk_GEN_ADD_MAX, mk_gen_add_n, mk_gen_add_power, mk_gen_add,  &
            mk_SYM_ADD_MAX, mk_sym_add_n, mk_sym_add_power, mk_sym_add )
   allocate(clu_at_lis  (  0:MK_MAX_SCAT      ))
   allocate(clu_mole_len(  0:CLU_MOLE_MAX_MOLE))
   allocate(clu_mole_off(  0:CLU_MOLE_MAX_MOLE))
   allocate(clu_mole_type( 0:CLU_MOLE_MAX_MOLE))
   allocate(clu_mole_char( 0:CLU_MOLE_MAX_MOLE))
   allocate(clu_mole_file( 0:CLU_MOLE_MAX_MOLE))
   allocate(clu_mole_dens( 0:CLU_MOLE_MAX_MOLE))
   allocate(clu_mole_biso( 0:CLU_MOLE_MAX_TYPE))
   allocate(clu_mole_clin( 0:CLU_MOLE_MAX_TYPE))
   allocate(clu_mole_cqua( 0:CLU_MOLE_MAX_TYPE))
   allocate(clu_mole_fuzzy(0:CLU_MOLE_MAX_MOLE))
   allocate(clu_mole_cont( 0:CLU_MOLE_MAX_ATOM))
   clu_at_lis = mk_at_lis
   call stru_internal_molecules(clu_infile, CLU_MOLE_MAX_MOLE,           &
              CLU_MOLE_MAX_TYPE, CLU_MOLE_MAX_ATOM, CLU_mole_num_mole, CLU_mole_num_type, &
              CLU_mole_num_atom, CLU_mole_len, CLU_mole_off, CLU_mole_type, CLU_mole_char,    &
              CLU_mole_file, CLU_mole_dens, CLU_mole_biso, CLU_mole_clin, CLU_mole_cqua,      &
              CLU_mole_fuzzy, CLU_mole_cont)
   loop_micro_atoms: do i=1, clu_mole_num_mole
      mc_matrix = 0.0
      mc_dimen  = 0.0
      mc_idimen = 0.0
      do j=1, 8
         k = clu_mole_cont(clu_mole_off(i)+j)
         call struc_read_one_atom_internal(clu_infile, k,  &
              clu_cr_pos, clu_cr_iscat, clu_cr_prop, clu_cr_surf, &
              clu_cr_magn, clu_cr_mole, clu_cr_moleatom )
         if(clu_at_lis(clu_cr_iscat)=='POSI') THEN
            mc_matrix(1:3,4) = clu_cr_pos(1:3)
         elseif(clu_at_lis(clu_cr_iscat)=='XAXI') THEN
            mc_matrix(1,1:3) = clu_cr_pos(1:3)
         elseif(clu_at_lis(clu_cr_iscat)=='YAXI') THEN
            mc_matrix(2,1:3) = clu_cr_pos(1:3)
         elseif(clu_at_lis(clu_cr_iscat)=='ZAXI') THEN
            mc_matrix(3,1:3) = clu_cr_pos(1:3)
         elseif(clu_at_lis(clu_cr_iscat)=='CENT') THEN
            mc_dimen (1:3,4) = clu_cr_pos(1:3)
         elseif(clu_at_lis(clu_cr_iscat)=='XDIM') THEN
            mc_dimen (1,1:3) = clu_cr_pos(1:3)
         elseif(clu_at_lis(clu_cr_iscat)=='YDIM') THEN
            mc_dimen (2,1:3) = clu_cr_pos(1:3)
         elseif(clu_at_lis(clu_cr_iscat)=='ZDIM') THEN
            mc_dimen (3,1:3) = clu_cr_pos(1:3)
         endif
      enddo
      mc_matrix(4,4) = 1.0
      mc_dimen (4,4) = 1.0
      DO k = 1, 3 
         DO l = 1, 3 
            mc_matrix(k, l) = mc_matrix(k, l) - mc_matrix(l, 4)
            mc_dimen (k, l) = mc_dimen (k, l) - mc_dimen (l, 4)
         ENDDO 
      ENDDO 
      mc_idimen = mc_dimen
!     CALL invmat4 (mc_idimen) 
      CALL matinv  (mc_dimen, mc_idimen) 
      infile = 'internal_domain_guest.' // clu_mole_file(i)(1:len_trim(clu_mole_file(i)))
      mc_type    = clu_mole_char(i)              ! Set character for this cluster
      md_sep_fuz = clu_mole_fuzzy(i)             ! Set fuzzy distance
      CALL micro_read_atoms(infile, mc_idimen, mc_matrix)                                                  
   enddo loop_micro_atoms
!
   deallocate(clu_at_lis)
   deallocate(clu_mole_len)
   deallocate(clu_mole_off)
   deallocate(clu_mole_type)
   deallocate(clu_mole_char)
   deallocate(clu_mole_file)
   deallocate(clu_mole_dens)
   deallocate(clu_mole_biso)
   deallocate(clu_mole_clin)
   deallocate(clu_mole_cqua)
   deallocate(clu_mole_fuzzy)
   deallocate(clu_mole_cont)
ENDIF 
IF ( .not. clu_infile_internal ) THEN
   CLOSE (imd) 
      if(allocated(clu_moles)) deallocate(clu_moles)
ENDIF 
!
!
chem_purge = .TRUE.     ! The crystal is most likely NOT periodic.
                        ! In the rare circumstances that is is the user
                        ! has to turn this on explicitly
chem_quick = .FALSE.
chem_period(:) = .FALSE.
if(allocated(clu_moles)) deallocate(clu_moles)
!
! Clean up internal files
!
call store_remove_single(host_file, ier_num)
if(clu_mode==CLU_IN_CLUSTER) THEN
   strucfile = clu_infile
   call store_remove_single(strucfile, ier_num)
endif
do i=1, clu_number
   if(clu_content(i)(1:22) == 'internal_domain_guest.') then
      strucfile = clu_content(i)
      call store_remove_single(strucfile, ier_num)
      if(ier_num==-113) ier_num = 0   ! File with equal name lreaned up prior
   endif
enddo
!                                                                       
END SUBROUTINE micro_filereading              
!
!*****7*****************************************************************
!
subroutine domain_irregular 
!-
!  Replace a contiguous irregularly shaped microdomain by the structure contained 
!  inside the cluster file
!+
!
!
use crystal_mod
use domain_mod 
use domaindis_mod 
use chem_mod
use prop_para_mod
use structur,  only:do_readstru
use spcgr_apply
USE read_internal_mod
!
use str_comp_mod
use string_convert_mod

!
implicit none
!
!
integer, parameter :: MAXMASK = 4
character(len=4) :: line   ! dummy line
integer :: i, ii         ! Dummy loop index
integer :: iatom         ! Dummy atom number
integer :: nprior        ! number of atoms prior to replacement
logical :: l_site
logical :: lout
logical, dimension(0:MAXMASK) :: uni_mask
!
lout = .false.
uni_mask(0)   = .true.
uni_mask(1:3) = .true.
uni_mask(4)   = .true.
!
!***  Read the list of pseudoatoms as structure 
!
l_site = .FALSE.
!
if(clu_infile(1:21) == 'internal_domain_list.') then
   call readstru_internal(MAXMASK, clu_infile, uni_mask)
else
   ier_num = -182
   ier_typ = ER_APPL
   ier_msg(1) = 'Domain files should have been internal'
   ier_msg(2) = 'Domain file error 01'
   return
!INTER 
!INTER call do_readstru(clu_infile, l_site)
!INTER call setup_lattice (cr_a0, cr_ar, cr_eps, cr_gten, cr_reps,             &
!INTER             cr_rten, cr_win, cr_wrez, cr_v, cr_vr, lout, cr_gmat,       &
!INTER             cr_fmat, cr_cartesian,                                      &
!INTER             cr_tran_g, cr_tran_gi, cr_tran_f, cr_tran_fi)
!INTER call get_symmetry_matrices
endif
!
do i=1,3
  if(cr_icc(i)==1) chem_period(i) = .FALSE.
end do
!
clu_remove_end = cr_natoms     ! Store atom number in original crystal
nprior = cr_natoms             ! Store atom number in original crystal
!
call domain_irreg_conn         ! Create connectivity with closest neighbors
do i=1, 1
   call domain_clean_conn         ! change atoms with 2 or less partners
   call domain_irreg_conn         ! Create connectivity with closest neighbors
enddo
!
loop_atoms: do iatom = 1, cr_natoms
   if(cr_iscat(iatom,1)==0) cycle loop_atoms ! NEEDS WORK 
!   if(btest(cr_prop(iatom), PROP_TEMP)) cycle loop_atoms   ! Temp property is set already, skip
!                              ! Determine cluster type for current pseudo atom
!                              ! Copy shape and orientation matrix
   line = cr_at_lis(cr_iscat(iatom,1))
   call do_cap(line(1:4))
   ii = 0 
   loop_types: DO i = 1, clu_number       ! Determine pseudo atom type
      if(str_comp(line(1:4), clu_name(i), 2, 4, 4) ) then
         ii = i 
         exit loop_types
      end if
   end do loop_types
   if(ii == 0) then       ! Pseudo atom type is unknown ERROR
      ier_num = -91 
      ier_typ = ER_APPL 
      ier_msg(1)  = ' Unknown pseudo atom name '//line(1:4)
      return 
   end if
!
!
   call domain_irreg_find(iatom, nprior, ii)  ! Replace the domain by guest
   if(ier_num/=0)  return                     ! Error like evaluating the "EXPR" 
end do loop_atoms
!
end subroutine domain_irregular
!
!*****7*****************************************************************
!
subroutine domain_irreg_conn
!-
!  Create a connectivity list for the pseudoatom structure
!+
!
use crystal_mod
use conn_mod
use conn_def_mod
!
use precision_mod
!
implicit none
!
character(len=PREC_STRING) :: line
integer :: i
integer :: lp
real(kind=PREC_DP) :: dist
!
call conn_reset        ! Reset all connectivities   (NEEDS WORK save old ?)
!
dist = cr_a0(1)*1.1    ! NEEDS WORK
do i=1, cr_nscat       ! Loop over all pseudo atoms
   write(line, '(a4, 4a,f8.4,a,i3.3)') '    ', &
            cr_at_lis(i)(1:len_trim(cr_at_lis(i))), ',', &
            cr_at_lis(i)(1:len_trim(cr_at_lis(i))), ', 0.1,', dist, &
            ', cd_',i
   lp = len_trim(line)
   call conn_do_set(code_add,line, lp)
enddo
line = ' '
call conn_create(line)
!
end subroutine domain_irreg_conn
!
!*****7*****************************************************************
!
subroutine domain_clean_conn
!-
!  Create a connectivity list for the pseudoatom structure
!+
!
use crystal_mod
use chem_mod
use atom_env_mod
use conn_sup_mod
use conn_def_mod
use do_find_mod
!
use precision_mod
!
implicit none
!
integer                    :: is1    ! central atom type
integer                    :: ino    ! number of connectivity list 
integer                    :: natoms ! number of atoms in connectivity list 
integer, dimension(:),   allocatable :: c_list ! Result of connectivity search
integer, dimension(:,:), allocatable :: c_offs ! Result of connectivity search
integer :: i, k
integer :: ianz      ! Number of general neighbors
integer :: iscat
logical                           :: fq
logical           , dimension(3)  :: fp
integer            , dimension(0:MAXSCAT) :: nscats
real(KIND=PREC_DP) , dimension(MAXSCAT)   :: werte
real(KIND=PREC_DP) , dimension(3)         :: pos
real(kind=PREC_DP) :: rmin
real(kind=PREC_DP) :: rmax
!
fp(1) = chem_period(1)
fp(2) = chem_period(2)
fp(3) = chem_period(3)
fq = chem_quick
!
do i=1, cr_natoms
  is1 = cr_iscat(i,1)
  ino = 1
  call get_connectivity_list(i, is1, ino, c_list, c_offs, natoms )
  if(natoms <=2) then
!write(*,*) ' AT ATOM ', i, is1, (' | ',c_list(k), cr_iscat(c_list(k,1)),k=1,natoms)
     nscats = 0
     ianz     = 1             ! Just one atom type
     werte(1) = -1            ! Find all atom types
     rmin = 0.01
     rmax = cr_a0(1) + 0.5    ! NEEDS WORK
     pos = cr_pos(:, i)
     call do_find_env(ianz, werte, MAXSCAT, pos, rmin, rmax, fq, fp)
     do k=1, atom_env(0)
        nscats(cr_iscat(atom_env(k),1)) = nscats(cr_iscat(atom_env(k),1)) + 1
!write(*,*) ' ENV ', i, k, atom_env(k), cr_iscat(atom_env(k),1)
     enddo
!write(*,*) ' replace ', i, cr_iscat(i,1), MAXLOC(nscats, 1)-1, '|', nscats
     iscat = MAXLOC(nscats, 1)-1
     cr_iscat(i,1) = iscat
!read(*,*) k
  endif
enddo
!
!
end subroutine domain_clean_conn
!
!*****7*****************************************************************
!
subroutine domain_irreg_find(iatom, nprior, iclu) !, infile, mc_dimen, mc_idimen, mc_matrix)
!-
!  Find all neighbors of starting atom
!  Replace this volume by guest
!+
use crystal_mod
use conn_sup_mod
use domain_mod
use domain_irreg_mod
use domaindis_mod
use prop_para_mod
use metric_mod
!
use errlist_mod
use precision_mod
!use lib_random_func
!
implicit none
!
integer,           intent(in) :: iatom            ! starting atom to work on
integer,           intent(in) :: nprior           ! Atom number in original pseudo atom list
integer,           intent(in) :: iclu             ! Cluster type
!
logical, parameter :: LSPACE=.TRUE.
!
character(len=PREC_STRING)    :: infile               ! Cluster content file
real(kind=PREC_DP)            :: mc_dimen (4, 4) 
real(kind=PREC_DP)            :: mc_idimen(4, 4) 
real(kind=PREC_DP)            :: mc_matrix(4, 4) 
!
integer                    :: i,j,k  ! Dummy counter
integer                    :: is1    ! central atom type
integer                    :: ino    ! number of connectivity list 
integer                    :: natoms ! number of atoms in connectivity list 
integer, dimension(:),   allocatable :: c_list ! Result of connectivity search
integer, dimension(:,:), allocatable :: c_offs ! Result of connectivity search
integer, dimension(:),   allocatable :: list   ! List of atoms in the current domain
integer, dimension(:,:), allocatable :: list_o ! List of atoms in the current domain offsets
logical                          :: lfound  ! Found new neighbors
real(kind=PREC_DP)               :: dist    ! dummy distance
real(kind=PREC_DP)               :: dmin    ! minimum distance
real(kind=PREC_DP), dimension(3) :: com     ! Center of mass for domain
real(kind=PREC_DP), dimension(3) :: uu      ! Dummy vector
!
if(cr_iscat(iatom,1)==0) return           ! Ignore VOIDS  NEEDS WORK ==> FLAG
!
!**** Get connectivity of this first atom in the new pseudo atom group
!
is1 = cr_iscat(iatom,1)
ino = 1
call get_connectivity_list(iatom, is1, ino, c_list, c_offs, natoms )
!
!if(natoms==0) return      !NEEDS WORK
!
allocate(list(1:nprior))             ! Long temporary list
allocate(list_o(3,1:nprior))         ! Long temporary list, including offsets
list = 0
list_o = 0
!
list(iatom) = 1                         ! Mark this atom as done
ngroup      = 1                         ! Currently one atom in the group
mgroup      = 1                         ! Currently one atom in the group
!write(*,*) ' FIRST ATOM ', natoms
do i=1, natoms
   if(c_list(i)>0) then
      list(c_list(i)) = -1 !Mark this as a new neighbor
!write(*,*) ' INITIAL NE ', c_list(i)
      list_o(:,c_list(i)) = c_offs(:, i)
   endif
enddo
!
!***  Search all atoms that are connected
!
lfound = .TRUE.
loop_neig: do
   if(.NOT. lfound) exit loop_neig
   lfound = .FALSE.
   loop_list: do j=1, nprior         ! Go through list to find new neighbors
      if(list(j)<0) then                ! This is a new neighbor
         is1 = cr_iscat(j,1)
         ino = 1
         call get_connectivity_list(j, is1, ino, c_list, c_offs, natoms )
         if(natoms==0) cycle loop_list  ! This one does not have a connectivity
         do i=1, natoms
           if(c_list(i)>0) then 
              if(list(c_list(i))==0) then 
                 list(c_list(i)) = -1 !Mark this as a new neighbor
!write(*,*) ' NEW NEIGHB ', c_list(i)
                 list_o(:,c_list(i)) = list_o(:,j) + c_offs(:,i)
               end if
            end if
         end do
         lfound = .TRUE.                        ! We found a new neighbor
         list(j) = 1
         ngroup = ngroup + 1
      endif
   end do loop_list
end do loop_neig
!
!write(*,*) ' NGROUP ', ngroup, nprior
!
allocate(short(1:ngroup))               ! make short group of members in this group
allocate(offs(3,ngroup))
allocate(coor(3,ngroup))
!
! Loop until all group members are replaced
!
k = 0                                  ! Modifier to handle no-replacement situation
loop_mgroup: do                        ! Loop over group until all atoms are replaced
   short = 0
   offs  = 0
   coor  = 0.0D0
!
   i     = 0
   com   = 0.0E0
!write(*,'(a,2x,a,4x,a,6x,a,4x,a)') ' member ', 'short', 'cr_pos', 'offs', 'effective'
   loop_build: do j=1, nprior
      if(list(j)==1) then
         if(btest(cr_prop(j), PROP_TEMP)) cycle loop_build
         i = i + 1
         short(i) = j
         offs(:,i) = list_o(:,j)
         com = com + cr_pos(:,j) + offs(:,i)
!write(*,'(a,i8,4(3f8.1,2x))') ' member ', short(i), cr_pos(:,short(i)), offs(:,(i))
      endif
   end do loop_build
   mgroup = i - k                      ! k is incremented if no atoms were replace
!write(*,*) 'MGROUP ', mgroup
   if(mgroup<=0) exit loop_mgroup
   com = com/real(mgroup, kind=PREC_DP)
   dmin = 1.0E9
   icent = 1
   loop_icent: do i=1, mgroup                    ! Find atom closest to COM
      if(btest(cr_prop(short(i)), PROP_TEMP)) cycle loop_icent
      uu = cr_pos(:,short(i)) + offs(:,i)
      dist = do_blen(lspace, uu, com)
!write(*,'(a,i5,3(3f6.1,2x))') 'dist',i, cr_pos(:,short(i)), offs(:,i), uu
      if(dist<dmin) then
         icent = i
         dmin = dist
      endif
   enddo loop_icent
!
!!write(*,*)
!write(*,*) ' group size ', ngroup, mgroup, ' ISCAT : ', cr_iscat(iatom,1)
!write(*,*) ' COM        ', com, icent
!
!                             ! Determine origin, coordinates of group elements
!
!call random_number(r1)
!icent = int(ngroup*r1) + 1
   icent_cr = short(icent)
!
   origin = cr_pos(:, icent_cr)
   orig_off = offs  (:, icent)
!write(*,'(a, (3f8.2,2x),3i5, 2i8)') ' Origin at ', origin, orig_off, icent, icent_cr
!read(*,*) i
   maxdim = 0.0D0
!write(*,'(a,2x,a,4x,a,6x,a,4x,a)') ' member ', 'short', 'cr_pos', 'offs', 'effective'
   do i=1, mgroup
      coor(:,i) = (cr_pos(:,short(i))+offs(:,(i))-origin(:)) - orig_off(:)
      maxdim(1,1) = min(maxdim(1,1), coor(1,i))
      maxdim(1,2) = max(maxdim(1,2), coor(1,i))
      maxdim(2,1) = min(maxdim(2,1), coor(2,i))
      maxdim(2,2) = max(maxdim(2,2), coor(2,i))
      maxdim(3,1) = min(maxdim(3,1), coor(3,i))
      maxdim(3,2) = max(maxdim(3,2), coor(3,i))
      cr_prop(short(i)) = ibclr(cr_prop(short(i)), PROP_TEMP)
!write(*,'(a,i8,4(3f8.1,2x),i5)') ' member ', short(i), cr_pos(:,short(i)), offs(:,(i)), &
! (cr_pos(:,short(i))-offs(:,(i))-origin(:)), coor(:,i), cr_prop(short(i))
   enddo
   maxdim(:,1) = maxdim(:,1) - 0.6
   maxdim(:,2) = maxdim(:,2) + 0.6
!
!***  Determine cluster type, shape and orientation matrix
!     The shape matrix is a unit matrix
!
   infile     = clu_content(iclu)        ! Guest structure
   mc_type    = clu_character(iclu)      ! Character of this domain
   md_sep_fuz = clu_fuzzy(iclu)          ! Distance to previous atoms, remove old at shorter
!
   call domain_set_matrix(iclu, mc_type, mc_dimen, mc_idimen, mc_matrix)
   if(ier_num/=0) then
      if(allocated(c_list)) deallocate(c_list)
      if(allocated(c_offs)) deallocate(c_offs)
      if(allocated(list))   deallocate(list)
      if(allocated(short))  deallocate(short)
      if(allocated(offs))   deallocate(offs)
      if(allocated(coor))   deallocate(coor)
      return                   ! Error evaluating the "EXPR"
   endif
   clu_current = iclu                       ! Transfer current cluster type
!
!if(iatom<5) then
!  write(*,'(a,4f9.6,1x)') 'MC ',mc_matrix(1,:)
!  write(*,'(a,4f9.6,1x)') 'MC ',mc_matrix(2,:)
!  write(*,'(a,4f9.6,1x)') 'MC ',mc_matrix(3,:)
!  write(*,'(a,4f9.6,1x)') 'MC ',mc_matrix(4,:)
!end if
!!!!!!!!!!!!!!!!!!!!!!!!!!
!write(*,*) ' PSEUDO ATOM ', iatom, cr_iscat(iatom,1), cr_at_lis(cr_iscat(iatom,1)), iclu, infile(1:len_trim(infile))
!                                         ! Now replace atoms
   mc_type = MD_DOMAIN_IRREG
!
!write(*,*) ' GROUP PRIOR REPLACE ', iatom, cr_natoms
i = cr_natoms
   call micro_read_atoms(infile, mc_idimen, mc_matrix) !
!write(*,*) ' GROUP AFTER REPLACE ', iatom, cr_natoms
!read(*,*) i
   if(i==cr_natoms) then                  ! No atoms were replaced, adapt modifier
      k = k + 1
      if(k>ngroup) then                    ! Failure could not replace any atoms
         if(allocated(c_list)) deallocate(c_list)
         if(allocated(c_offs)) deallocate(c_offs)
         if(allocated(list))   deallocate(list)
         if(allocated(short))  deallocate(short)
         if(allocated(offs))   deallocate(offs)
         if(allocated(coor))   deallocate(coor)
         ier_num = 1
         ier_typ = ER_FORT
         return
      endif
   else
      k = 0                              ! Reset modifier, as atoms were replaced
   endif
   do i=1,mgroup
      if(btest(cr_prop(short(i)), PROP_TEMP)) then
         cr_iscat(short(i),1) = 0
!  cr_prop(short(i)) = ibset(cr_prop(short(i)), PROP_TEMP)
!   write(*,*) ' PROP ', i, short(i), cr_iscat(short(i),1), cr_prop(short(i))
      endif
   enddo
!read(*,*) i
enddo loop_mgroup
!
!call domain_irreg_replace(infile)
!
if(allocated(c_list)) deallocate(c_list)
if(allocated(c_offs)) deallocate(c_offs)
if(allocated(list))   deallocate(list)
if(allocated(short))  deallocate(short)
if(allocated(offs))   deallocate(offs)
if(allocated(coor))   deallocate(coor)
!
end subroutine domain_irreg_find
!
!*****7*****************************************************************
!
subroutine domain_irreg_test_atom(vv, linside, ddd)
!-
!  Tests if an atom at coordinates vv is inside the current domain
!  Apply periodic boundary conditions to place near actual pseudo atom
!+
!
use domain_mod
use domain_irreg_mod
use crystal_mod
use metric_mod
use prop_para_mod
!
use precision_mod
!
implicit none
!
real(kind=PREC_DP), dimension(3), intent(inout) :: vv
logical                         , intent(out)   :: linside
real(kind=PREC_DP) :: ddd
!
logical, parameter :: lspace = .TRUE.
!
integer            :: i       ! dummy loop index
integer            :: ishort  ! Index of closest atom
real(kind=PREC_DP) :: dist    ! distance to pseudo atoms
real(kind=PREC_DP) :: d       ! distance to pseudo atoms
real(kind=PREC_DP), dimension(3) :: uu  ! Vector atom to pseudo
real(kind=PREC_DP), dimension(3) :: NULLV!  ! Null vector
 data NULLV / 0.0D0, 0.0D0, 0.0D0 /
!
linside = .false.                         ! Assume outside
!
!write(*,*) 'testing ', vv, origin
dist = 1.e9
if(maxdim(1,1)<vv(1) .and. vv(1)<maxdim(1,2)   .and.                            &
   maxdim(2,1)<vv(2) .and. vv(2)<maxdim(2,2)   .and.                            &
   maxdim(3,1)<vv(3) .and. vv(3)<maxdim(3,2)         ) then     ! Initial crude estimate
   ishort = 0
   loop_ishort: do i=1, mgroup                         ! Loop over all atoms to find closest
      if(btest(cr_prop(short(i)), PROP_TEMP)) cycle loop_ishort
      uu = vv - coor(:,i)
      d  = do_blen(lspace, uu, NULLV)
      if(d<dist) then
         dist = d
         ishort = i
      endif
      if(d<dist) then
         dist = d
         ishort = i
      endif
   enddo loop_ishort
!
   if(dist < clu_fuzzy(clu_current)) then
!if(cr_natoms < clu_remove_end+20) then
! write(*,'(a, 3f6.1,2i6, 6(3f6.1,1x),3i5,2f6.1)') ' vv ', vv, ishort, short(ishort), coor(:,ishort), &
! cr_pos(:,short(ishort)), offs(:,ishort), dist, clu_fuzzy(clu_current)
!!
!end if
      vv = vv + origin - offs(:,ishort) + orig_off
      linside = .true.
      cr_prop(short(ishort)) = ibset(cr_prop(short(ishort)), PROP_TEMP)
else
!if(cr_natoms < clu_remove_end+20) then
! write(*,'(a, 3f6.1,2i6, 6(3f6.1,1x),3i5,2f6.1)') ' xx ', vv, ishort, short(ishort), coor(:,ishort), &
! cr_pos(:,short(ishort)), offs(:,ishort), dist, clu_fuzzy(clu_current)
!endif
   endif
endif
ddd = dist
!
end subroutine domain_irreg_test_atom
!
!*****7*****************************************************************
!
subroutine domain_set_matrix(iclu, mc_type, mc_dimen, mc_idimen, mc_matrix)
!-
!  Set the current shape / orientation matrix 
!  either from clu_orient /clu_shape
!  or from "EXPR"
!+
!
use domain_mod
!
use ber_params_mod
use do_replace_expr_mod
use errlist_mod
use matrix_mod
use precision_mod
!
implicit none
!
integer, intent(in) :: iclu            ! Current cluster number
integer, intent(in) :: mc_type         ! Current cluster type
real(kind=PREC_DP)   , intent(out)           :: mc_dimen (4, 4) 
real(kind=PREC_DP)   , intent(out)           :: mc_idimen(4, 4) 
real(kind=PREC_DP)   , intent(out)           :: mc_matrix(4, 4) 
!
integer, parameter :: MAXW = 12
character(len=PREC_STRING), dimension(MAXW) :: cpara
integer                   , dimension(MAXW) :: lpara
real(kind=PREC_DP)        , dimension(MAXW) :: werte
real(kind=PREC_DP)        , dimension(4,4)  :: mat_dir
real(kind=PREC_DP)        , dimension(4,4)  :: mat_inv
!
integer :: i, j
integer :: ianz
!
call do_evaluate_expr                      ! Determine current EXPR value
if(ier_num/=0) return
if(clu_use_expr(iclu,1)) then              ! Orientation uses "EXPR", use current value
   ianz = 0
   do i = 1, 3
      do j = 1,4
         ianz = ianz + 1
         cpara(ianz) = clu_or_expr(iclu, i, j)  ! Transfer expressions into cpara for berechne
         lpara(ianz) = len_trim(cpara(ianz))
      enddo
   enddo
   call do_use_expr(ianz, cpara, lpara, MAXW)        ! Use evaluated EXPR
   call ber_params(ianz, cpara, lpara, werte, MAXW)  ! Evaluate, with current "EXPR"
   if(ier_num/=0)  then
      ier_msg(1) = 'Expression for orient erroneous '
      return
   endif
   ianz = 0
   do i = 1, 3
      do j = 1,4
         ianz = ianz + 1
         mc_matrix(i, j) = werte(ianz)         ! Transfer current results into orient matrix
      enddo
   enddo
else                                       ! Orientation uses clu_orient
   do i = 1, 3
      do j = 1,4
         mc_matrix(i, j) = clu_orient(iclu, i, j)
      end  do
      mc_matrix(i, 4) = clu_orient(iclu, i, 4) !+ cr_pos(i,icent_cr)
      mc_matrix(4, i) = 0.0
   enddo
endif
if(clu_use_expr(iclu,2) .and. mc_type/=CLU_IN_IRREG) then   ! Shape       uses "EXPR"; not an irregular character
   ianz = 0
   do i = 1, 3
      do j = 1,4
         ianz = ianz + 1
         cpara(ianz) = clu_or_expr(iclu, i, j)  ! Transfer expressions into cpara for berechne
         lpara(ianz) = len_trim(cpara(ianz))
      enddo
   enddo
   call do_use_expr(ianz, cpara, lpara, MAXW)        ! Use evaluated EXPR
   call ber_params(ianz, cpara, lpara, werte, MAXW)  ! Evaluate current "EXPR"
   if(ier_num/=0)  then
      ier_msg(1) = 'Expression for orient erroneous '
      return
   endif
   ianz = 0
   do i = 1, 3
      do j = 1,4
         ianz = ianz + 1
         mc_dimen (i, j) = werte(ianz)         ! Transfer current results into orient matrix
      enddo
   enddo
   mc_dimen(4,:)     = 0.0
   mc_dimen(4,4)     = 1.0
else                                       ! Shape       uses clu_shape
   if(mc_type==CLU_IN_IRREG) then          ! Irregular shape, use unit matrix
      mc_dimen  = 0.0
      mc_idimen = 0.0
      do i = 1, 4
         mc_dimen (i,i) = 1.0
         mc_idimen(i,i) = 1.0
      enddo
   else
      mc_dimen(1:3,1:4) = clu_shape(iclu, 1:3, 1:4)
      mc_dimen(4,:)     = 0.0
      mc_dimen(4,4)     = 1.0
   endif
endif
mat_dir = mc_dimen
call matinv4(mat_dir, mat_inv)
mc_idimen = mat_inv
!
!
end subroutine domain_set_matrix
!
!*****7*****************************************************************
!
SUBROUTINE micro_read_atoms (infile, mc_idimen, mc_matrix)
!-
! Read the guest atoms and insert into host structure
!+
!
USE discus_config_mod 
USE crystal_mod 
USE metric_mod
USE micro_mod 
USE read_internal_mod 
USE discus_save_mod 
USE structur, ONLY: stru_readheader
USE errlist_mod 
USE molecule_mod
USE support_mod
!                                                                       
IMPLICIT none 
       
!                                                                       
CHARACTER(len=*)     , intent(in) :: infile 
REAL(kind=PREC_DP) , dimension(4,4), intent(in) :: mc_idimen
REAL(kind=PREC_DP) , dimension(4,4), intent(in) :: mc_matrix
!                                                                       
INTEGER, PARAMETER  :: ist = 46
!                                                                       
INTEGER, PARAMETER                   :: AT_MAXP = 16
INTEGER                              :: at_ianz
CHARACTER(LEN=8), DIMENSION(AT_MAXP) :: at_param
LOGICAL  :: lread 
!                                                                       
lread = .true. 
IF(infile(1:8)=='internal') THEN
   mk_infile_internal = .true.
   CALL stru_readheader_internal (infile, MK_MAX_SCAT, mk_name,   &
         mk_spcgr, mk_spcgr_set, mk_set, mk_iset,                       &
         mk_at_lis, mk_nscat, mk_dw, mk_occ, mk_a0, mk_win,   &
         sav_ncell, sav_r_ncell, sav_ncatoms, spcgr_ianz, spcgr_para, &
         mk_GEN_ADD_MAX, mk_gen_add_n, mk_gen_add_power, mk_gen_add,  &
         mk_SYM_ADD_MAX, mk_sym_add_n, mk_sym_add_power, mk_sym_add )
   mk_iatom = 0
   at_param(1) = 'X'
   at_param(2) = 'Y'
   at_param(3) = 'Z'
   at_param(4) = 'BISO'
   at_param(5) = 'PROPERTY'
   at_param(6) = 'MOLENO'
   at_param(7) = 'MOLEAT'
   at_param(8) = 'OCC'
   at_ianz = 8
ELSE
         mk_infile_internal = .false.
         CALL oeffne (ist, infile, 'unknown') 
         IF (ier_num.ne.0) return 
!                                                                       
!     Read the input file header                                        
!                                                                       
         CALL stru_readheader (ist, MK_MAX_SCAT, mk_name,     &
         mk_spcgr, mk_set, mk_at_lis, mk_nscat, mk_dw, mk_occ, mk_anis, mk_a0, mk_win, sav_ncell,   &
         sav_r_ncell, sav_ncatoms, mk_spcgr_ianz, mk_spcgr_para, &
         AT_MAXP, at_ianz, at_param)           
         IF(ier_num /= 0) THEN
            ier_msg(1) = ' Guest file ' // infile
            CLOSE(IST)
            RETURN
         ENDIF
ENDIF
!                                                                       
CALL micro_read_atom(ist, infile, mc_idimen, mc_matrix, &
                     AT_MAXP, at_ianz, at_param) 
IF(.not.(infile(1:8)=='internal')) THEN
   CLOSE (ist) 
ENDIF
!                                                                       
END SUBROUTINE micro_read_atoms               
!
!*****7*****************************************************************
!
      SUBROUTINE micro_remove_old (mc_idimen) 
!                                                                       
!     This routine will become obsolete if the "fuzzy" property is      
!     always used.                                                      
!                                                                       
      USE discus_config_mod 
      USE crystal_mod 
      USE domaindis_mod 
      USE micro_mod 
      USE metric_mod
      USE prop_para_mod 
      USE errlist_mod 
use precision_mod
      IMPLICIT none 
!                                                                       
       
!                                                                       
      REAL(kind=PREC_DP) :: mc_idimen (4, 4) 
!                                                                       
!     INTEGER idim 
!     PARAMETER (idim = 4) 
!                                                                       
      INTEGER i, j 
      LOGICAL lspace 
      REAL(kind=PREC_DP) :: d 
      REAL(kind=PREC_DP) :: u (3), v (3) 
      REAL(kind=PREC_DP) :: u4 (4), v4 (4) 
      REAL(kind=PREC_DP) :: radius (3) 
      REAL(kind=PREC_DP), dimension(3), parameter :: NULLV = (/ 0.0D0, 0.0D0, 0.0D0 /)
!                                                                       
      IF (mc_type .eq.MD_DOMAIN_FUZZY) then 
         RETURN 
      ENDIF 
!                                                                       
      lspace = .true. 
      v (1) = 0.0 
      v (2) = 0.0 
      v (3) = 0.0 
      DO i = 1, 3 
      DO j = 1, 3 
      u (1) = 0.0 
      u (2) = 0.0 
      u (3) = 0.0 
      ENDDO 
      u (i) = 1.0 
      radius (i) = do_blen (lspace, u, v) 
      ENDDO 
!                                                                       
      DO i = 1, cr_natoms 
      DO j = 1, 3 
      u4 (j) = cr_pos (j, i) 
      ENDDO 
      u4 (4) = 1.0 
!     CALL trans (u4, mc_idimen, v4, idim) 
      v4 = matmul(mc_idimen, u4)
      DO j = 1, 3 
      v (j) = v4 (j) 
      ENDDO 
      IF (mc_type  .eq.MD_DOMAIN_SPHERE) then 
         d = do_blen (lspace, v, NULLV) 
         IF (d.lt.radius (2) ) then 
            cr_iscat (i,1) = 0 
            cr_prop (i) = IBCLR (cr_prop (i), PROP_NORMAL) 
            cr_prop (i) = IBSET (cr_prop (i), PROP_DOMAIN) 
         ENDIF 
      ELSEIF (mc_type  .eq.MD_DOMAIN_CUBE) then 
         IF (abs (v (1) ) .le.1.and.abs (v (2) ) .le.1.and.abs (v (3) ) &
         .le.1) then                                                    
            cr_iscat (i,1) = 0 
            cr_prop (i) = IBCLR (cr_prop (i), PROP_NORMAL) 
            cr_prop (i) = IBSET (cr_prop (i), PROP_DOMAIN) 
         ENDIF 
      ELSEIF (mc_type  .eq.MD_DOMAIN_CYLINDER) then 
         IF (abs (v (3) ) .le.1) then 
            v (3) = 0.0 
            d = do_blen (lspace, v, NULLV) 
            IF (d.lt.radius (2) ) then 
               cr_iscat (i,1) = 0 
               cr_prop (i) = IBCLR (cr_prop (i), PROP_NORMAL) 
               cr_prop (i) = IBSET (cr_prop (i), PROP_DOMAIN) 
            ENDIF 
         ENDIF 
      ENDIF 
      ENDDO 
!                                                                       
      END SUBROUTINE micro_remove_old               
!
!*****7*****************************************************************
!
SUBROUTINE micro_read_atom(ist, infile, mc_idimen, mc_matrix, &
                            AT_MAXP, at_ianz, at_param) 
!-
!  Reads the atoms of the current cluster type 
!+
!                                                                       
USE discus_config_mod 
USE discus_allocate_appl_mod 
USE crystal_mod 
USE domain_mod 
USE domaindis_mod 
use domain_irreg_mod
USE metric_mod
USE micro_mod 
USE molecule_mod 
USE prop_para_mod 
USE read_internal_mod
USE structur, ONLY: struc_mole_header !, read_atom_line
USE spcgr_apply, ONLY: mole_insert_current, mole_insert_explicit
USE surface_mod 
USE errlist_mod 
!
use blanks_mod
USE lib_errlist_func
USE lib_length
USE precision_mod
USE str_comp_mod
USE string_convert_mod
!                                                                       
IMPLICIT none 
!                                                                       
INTEGER            , INTENT(IN) :: ist 
CHARACTER (LEN = *), INTENT(IN) :: infile
REAL(kind=PREC_DP) , INTENT(IN) :: mc_idimen (4, 4) 
REAL(kind=PREC_DP) , INTENT(IN) :: mc_matrix (4, 4) 
INTEGER                                  , INTENT(IN)  :: AT_MAXP
INTEGER                                  , INTENT(IN)  :: at_ianz
CHARACTER(LEN=8), DIMENSION(AT_MAXP)     , INTENT(IN)  :: at_param
!
!                                                                       
!     INTEGER idim4 
!     PARAMETER (idim4 = 4) 
      INTEGER maxw 
      PARAMETER (maxw = 12) 
!                                                                       
      CHARACTER(10) befehl 
      CHARACTER(LEN=PREC_STRING) :: line, zeile 
      INTEGER i, j, k, ii, ibl, lbef 
      INTEGER natoms_old 
      INTEGER i_count 
      INTEGER              :: new_nscat ! DUMMY for allocation
      INTEGER              :: new_nmax  ! DUMMY for allocation
      LOGICAL              :: lcontent
      LOGICAL lspace 
      LOGICAL linside
      REAL(kind=PREC_DP) :: d 
      REAL(kind=PREC_DP) :: u (4), v (4), w (4) 
      REAL(kind=PREC_DP) :: vv (3) , ww(3)
      REAL(kind=PREC_DP) :: radius (3) 
      REAL(KIND=PREC_DP):: werte (maxw)
      REAL(kind=PREC_DP) :: dw1 
      REAL(kind=PREC_DP) ::  NULLV (3) 
      REAL(kind=PREC_DP) :: shortest 
      REAL(kind=PREC_DP), DIMENSION(3)  :: xyz ! Atom position
      INTEGER             :: dummy_iscat ! Atom type
      INTEGER             :: dummy_prop ! Atom property
      INTEGER             :: dummy_mole ! Atom is in this molecule
      INTEGER             :: dummy_moleatom ! Atom is in this molecul at number
      INTEGER, DIMENSION(0:3)  :: dummy_surf
      REAL(kind=PREC_DP)   , DIMENSION(0:3)  :: dummy_magn
      INTEGER             :: natoms ! Number of molecules
      INTEGER             :: nscat  ! Number of molecule types
      INTEGER             :: TEMP_MAX_MOLE ! Number of molecules
      INTEGER             :: TEMP_MAX_TYPE ! Number of molecule types
      INTEGER             :: TEMP_MAX_ATOM ! Number of atoms in molecules
      INTEGER             :: temp_num_mole ! Number of molecules
      INTEGER             :: temp_num_type ! Number of molecule types
      INTEGER             :: temp_num_atom ! Number of atoms in molecules
      INTEGER             :: istatus! Number of atoms in molecules
      INTEGER, DIMENSION(:), ALLOCATABLE :: temp_mole_len
      INTEGER, DIMENSION(:), ALLOCATABLE :: temp_mole_off
      INTEGER, DIMENSION(:), ALLOCATABLE :: temp_mole_type
      INTEGER, DIMENSION(:), ALLOCATABLE :: temp_mole_char
      CHARACTER (LEN=200), DIMENSION(:), ALLOCATABLE :: temp_mole_file
      REAL(kind=PREC_DP)   , DIMENSION(:), ALLOCATABLE :: temp_mole_dens
      REAL(kind=PREC_DP)   , DIMENSION(:), ALLOCATABLE :: temp_mole_biso
      REAL(kind=PREC_DP)   , DIMENSION(:), ALLOCATABLE :: temp_mole_clin
      REAL(kind=PREC_DP)   , DIMENSION(:), ALLOCATABLE :: temp_mole_cqua
      REAL(kind=PREC_DP)   , DIMENSION(:), ALLOCATABLE :: temp_mole_fuzzy
      INTEGER, DIMENSION(:), ALLOCATABLE :: temp_mole_cont
      LOGICAL, DIMENSION(:), ALLOCATABLE :: temp_present
      INTEGER, DIMENSION(:), ALLOCATABLE :: temp_in_crystal
      INTEGER, DIMENSION(:), ALLOCATABLE :: temp_in_mole
      INTEGER             :: n_mole
      INTEGER             :: n_type
      INTEGER             :: n_atom
      INTEGER             :: n_read_atoms
      INTEGER             :: n_mole_old = 0    ! previous number of molecules in structure
      INTEGER             :: iimole
      INTEGER             :: lline
      LOGICAL             :: linternal
      LOGICAL, SAVE       :: at_init = .TRUE.
logical :: lmole_type    ! The molecule has a 'molecule  type' line
integer :: is_mole_type  ! Molecule type provided on the 'molecule type line'
REAL(kind=PREC_DP) :: d100, d010, d001
!                                                                       
      DATA NULLV / 0.0D0, 0.0D0, 0.0D0 / 
!
!     Define distances to domain surfaces
!
      v(:) = 0.0
      v(1) = surf_in_dist(SURF_MAXSCAT)/cr_a0(1)
!     CALL trans (v, mc_idimen, w, idim4)
      w = matmul(mc_idimen, v)
      d100 = sqrt(w(1)*w(1) + w(2)*w(2) + w(3)*w(3) )
      v(:) = 0.0
      v(2) = surf_in_dist(SURF_MAXSCAT)/cr_a0(2)
!     CALL trans (v, mc_idimen, w, idim4)
      w = matmul(mc_idimen, v)
      d010 = sqrt(w(1)*w(1) + w(2)*w(2) + w(3)*w(3) )
      v(:) = 0.0
      v(3) = surf_in_dist(SURF_MAXSCAT)/cr_a0(3)
!     CALL trans (v, mc_idimen, w, idim4)
      w = matmul(mc_idimen, v)
      d001 = sqrt(w(1)*w(1) + w(2)*w(2) + w(3)*w(3) )
!
      n_read_atoms = 0
      lspace  = .true.
      linside = .false.
      v (1) = 0.0 
      v (2) = 0.0 
      v (3) = 0.0 
      DO i = 1, 3 
         u (1) = 0.0 
         u (2) = 0.0 
         u (3) = 0.0 
         u (i) = 1.0 
         radius (i) = do_blen (lspace, u(1:3), v(1:3)) 
      ENDDO 
      i_count = 0 
!                                                                       
      DO i = 1, 3 
         mk_dim (i, 1) = mc_matrix (i, 4) 
         mk_dim (i, 2) = mc_matrix (i, 4) 
      ENDDO 
!write(*,*) 'MC_MATRIX '
!do i=1,4
!  write(*,*)  mc_matrix(i,:)
!enddo
!write(*,*) 'MD_IDIM   '
!do i=1,4
!  write(*,*)  mc_idimen(i,:)
!enddo
!
!     Does the user want to remove existing atoms at insertion of each domain
!     which includes atoms inserted by other domains in this run or
!     only the atoms that existed previous to the 'run' command
!
      IF (     clu_remove_mode ==  CLU_REMOVE_STRICT  ) THEN
         natoms_old = cr_natoms 
      ELSE
         natoms_old = clu_remove_end
      ENDIF 
      at_init = .TRUE.
!                                                                       
      IF(mk_infile_internal) THEN
         mk_iatom = 0
         ALLOCATE ( temp_present   (0:MK_MAX_ATOM), STAT = istatus )
         ALLOCATE ( temp_in_crystal(0:MK_MAX_ATOM), STAT = istatus )
         ALLOCATE ( temp_in_mole   (0:MK_MAX_ATOM), STAT = istatus )
         temp_present    = .false.
         temp_in_crystal = 0
         temp_in_mole    = 0
      ENDIF
      IF(cr_natoms == 0) THEN
         n_mole_old = 0
      ELSE
         n_mole_old = mole_num_mole
      ENDIF
lmole_type = .false.                   ! Assume no 'molecule type' line
is_mole_type = 1
 1000 CONTINUE 
   ier_num = - 49 
   ier_typ = ER_APPL 
   line = ' ' 
      IF(mk_infile_internal) THEN
         mk_iatom =  mk_iatom + 1   ! Increment internal atom number
         CALL struc_read_one_atom_internal(infile,  mk_iatom,  &
              xyz, dummy_iscat, dummy_prop, dummy_surf, dummy_magn, dummy_mole, dummy_moleatom )
         IF( ier_num == -105 ) THEN  ! read "end_of_file" 
            ier_num = 0
            ier_typ = ER_NONE
             mk_iatom =  mk_iatom - 1   ! Increment internal atom number
            GOTO 2
         ELSEIF(ier_num/=0) then
            RETURN
         ENDIF
         dummy_mole = 0      ! Internal molecule atoms are set further down
         dummy_moleatom = 0  ! Internal molecule atoms are set further down
         WRITE(line, 3000) mk_at_lis(dummy_iscat), xyz, mk_dw(dummy_iscat), &
                           dummy_prop, dummy_mole, dummy_moleatom, 1.0 ! copy into line
!write(*,*) 'READING ATOMS ', mk_iatom,'>',line(1:len_trim(line))
         werte(1:3) = xyz
         werte(4)   = mk_dw(dummy_iscat)
      ELSE
         READ (ist, 2000, END = 2, ERR = 999) line 
      ENDIF
      lline = len_str (line)
      call tab2blank(line,lline)
blank1: IF (line.ne.' '.AND.line (1:1) .ne.'#'.AND.line(1:1)/='!'.AND.line.ne.char (13) ) THEN
         i_count = i_count + 1 
         ibl = index (line, ' ') + 1 
         lbef = 10 
         befehl = ' ' 
         ibl = index (line, ' ') 
         lbef = min (ibl - 1, lbef) 
         befehl = line (1:lbef) 
         lbef = len_str (befehl) 
         befehl = line (1:lbef) 
is_mole: IF (str_comp (befehl, 'molecule', 4, lbef, 8) .or. &
             str_comp (befehl, 'domain',   4, lbef, 6) .or. &
             str_comp (befehl, 'object',   4, lbef, 6)      ) then                                             
!                                                                       
!     ------Start/End of a molecule                                     
!                                                                       
            CALL no_error 
            zeile = line (ibl:200) 
            i = 200 - ibl 
            i = len_str (zeile) 
            CALL struc_mole_header (zeile, i, .false.,lcontent) 
            if(index(zeile,'type')>0) then
               lmole_type = .true.                    ! Input file provided 'molecule type'
               is_mole_type = mole_type(mole_num_mole)
            elseif(index(zeile,'end')>0) then          ! 'molecule end' line
               if(lmole_type) then                    ! Input file provided 'molecule type'
                  mole_type(mole_num_mole) = clu_mole_tab(is_mole_type, clu_current)
                  mole_num_type = max(mole_num_type, mole_type(mole_num_mole))
               else
                  mole_type(mole_num_mole) = clu_mole_tab(1           , clu_current)
                  mole_num_type = max(mole_num_type, mole_type(mole_num_mole))
                  lmole_type = .false.                ! Reset flag for next molecule
               endif
            endif
            IF (ier_num.ne.0) return 
         ELSE is_mole
!           READ (line (ibl:80), *, end = 999, err = 999) (werte (j), j = 1, 4)
!write(*,*) ' LINE ', line(1:len_trim(line)), ' INFILE ', infile(1:len_trim(infile))
!write(*,*) 'at_init', at_init, lline, ibl, at_ianz, AT_MAXP
!   write(*,*) 'werte, wwerte ', ubound(werte),ubound(at_param)
!           CALL read_atom_line (line, ibl, lline, as_natoms, maxw, werte, &
!                AT_MAXP, at_ianz, at_param, at_init)
            n_read_atoms = n_read_atoms + 1
            DO j = 1, 3 
               u (j) = werte (j) 
            ENDDO 
            u (4) = 1.0 
!           CALL trans (u, mc_matrix, v, idim4) 
!           CALL trans (v, mc_idimen, w, idim4) 
            v = matmul(mc_matrix, u)
            w = matmul(mc_idimen, v)
            DO j = 1, 3 
               vv (j) = w (j) 
               mk_dim (j, 1) = min (mk_dim (j, 1), v (j) ) 
               mk_dim (j, 2) = max (mk_dim (j, 2), v (j) ) 
            ENDDO 
            linternal = .FALSE.
            IF (mc_type  .eq.MD_DOMAIN_SPHERE) then 
               d = do_blen (lspace, vv, NULLV) 
               linside = (d.lt.radius (2) ) 
               IF(linside) THEN
                 linternal = (radius (2)-d <= surf_in_dist(SURF_MAXSCAT)/cr_a0(2))
               ENDIF
            ELSEIF (mc_type  .eq.MD_DOMAIN_CUBE) then 
               linside = ABS(vv(1)) <= 1. .AND. ABS(vv(2)) <= 1. .AND.   &
                         ABS(vv(3)) <= 1.                            
               IF(linside) THEN
                 linternal = (ABS(-1.-vv(1)) <= d100 ) .OR. &
                             (ABS( 1.-vv(1)) <= d100 ) .OR. &
                             (ABS(-1.-vv(2)) <= d010 ) .OR. &
                             (ABS( 1.-vv(2)) <= d010 ) .OR. &
                             (ABS(-1.-vv(3)) <= d001 ) .OR. &
                             (ABS( 1.-vv(3)) <= d001 )
               ENDIF
            ELSEIF (mc_type  .eq.MD_DOMAIN_CYLINDER) then 
               linside = .false. 
               d = 0.0
               IF (abs (vv (3) ) .le.1) then 
                  ww(:) = vv(:)
                  ww (3) = 0.0 
                  d = do_blen (lspace, ww, NULLV) 
                  linside = d.lt.radius (2) 
               ENDIF 
               IF(linside) THEN
                 linternal = (radius (2)-d <= surf_in_dist(SURF_MAXSCAT)/cr_a0(2)) .or. &
                             (ABS(-1.-vv(3)) <= d001 ) .or. &
                             (ABS( 1.-vv(3)) <= d001 )
               ENDIF
            ELSEIF (mc_type  .eq.MD_DOMAIN_FUZZY) then 
               linside = .true. 
            ELSEIF(mc_type ==  MD_DOMAIN_IRREG) then           ! Irregularly shaped domain
               call domain_irreg_test_atom(v(1:3), linside, d)        ! Test if atom is inside, including periodic boundaries
!write(*,*) ' TESTED ', v(1:3), linside
!if(cr_natoms < clu_remove_end+20 )then
!write(*,*) ' ATOM ', u(1:3), v(1:3)
!end if
!               linside = .true. 
            ENDIF 
inside:     IF (linside) then 
!write(*,*) ' INSIDE ', v
               IF (cr_natoms.eq.nmax) then 
                  new_nmax = MAX(nmax + 100, int( nmax*1.1 ))
                  CALL alloc_crystal_scat(MAXSCAT)
                  CALL alloc_crystal_nmax(new_nmax)
               ENDIF
               IF (cr_natoms >  nmax) then 
!                                                                       
!     --------Too many atoms in the structure file                      
!                                                                       
                  ier_num = - 10 
                  ier_typ = ER_APPL 
                  RETURN 
               ENDIF 
               cr_natoms = cr_natoms + 1 
               i = cr_natoms 
               DO j = 1, 3 
               cr_pos (j, i) = v (j)
               cr_dim (j, 1) = min (cr_dim (j, 1), cr_pos (j, i) ) 
               cr_dim (j, 2) = max (cr_dim (j, 2), cr_pos (j, i) ) 
               ENDDO 
               dw1 = werte (4) 
noblank:      IF (line (1:4) .ne.'    ') then 
                  ibl = ibl - 1 
                  CALL do_cap (line (1:ibl) ) 
                  DO j = 0, cr_nscat 
                  IF (line (1:ibl) .eq.cr_at_lis (j) .and.dw1.eq.cr_dw (&
                  j) ) then                                             
                     cr_iscat (i,1) = j 
                     cr_prop (i) = 0 
                     cr_prop (i) = ibset (cr_prop (i), PROP_NORMAL) 
                     cr_prop (i) = ibset (cr_prop (i), PROP_DOMAIN) 
                     IF(linternal) THEN
                        cr_prop (i) = ibset (cr_prop (i), PROP_SURFACE_INT) 
                     ENDIF
                     GOTO 11 
                  ENDIF 
                  ENDDO 
                  IF (cr_nscat == MAXSCAT) then 
                     new_nscat = MAX(MAXSCAT + 5, INT ( MAXSCAT * 1.025 ))
                     CALL alloc_crystal_scat(new_nscat)
                     CALL alloc_crystal_nmax(NMAX)
                     CALL alloc_surf   (new_nscat)
                  ENDIF 
                  IF (MAXSCAT  > ubound(surf_in_dist,1)) THEN
                     new_nscat = MAXSCAT
                     CALL alloc_surf   (new_nscat)
                  ENDIF 
                  IF (cr_nscat.eq.MAXSCAT) then 
!                                                                       
!     --------  Too many atom types in the structure file               
!                                                                       
                     ier_num = - 72 
                     ier_typ = ER_APPL 
                     RETURN 
                  ENDIF 
                  cr_nscat = cr_nscat + 1 
                  cr_iscat (i,1) = cr_nscat 
                  cr_prop (i) = 0 
                  cr_prop (i) = ibset (cr_prop (i), PROP_NORMAL) 
                  cr_prop (i) = ibset (cr_prop (i), PROP_DOMAIN) 
                  IF(linternal) THEN
                     cr_prop (i) = ibset (cr_prop (i), PROP_SURFACE_INT) 
                  ENDIF
                  cr_at_lis (cr_nscat) = line (1:ibl) 
                  cr_dw (cr_nscat) = dw1 
                  cr_occ(cr_nscat) = 1.00D0    !! WORK OCC
!                                                                       
                  IF (0.0.le.cr_pos (1, i) .and.cr_pos (1, i)           &
                  .lt.1.and.0.0.le.cr_pos (2, i) .and.cr_pos (2, i)     &
                  .lt.1.and.0.0.le.cr_pos (3, i) .and.cr_pos (3, i)     &
                  .lt.1) then                                           
                     as_natoms = as_natoms + 1 
                     as_at_lis (cr_nscat) = cr_at_lis (cr_nscat) 
                     as_iscat (as_natoms) = cr_iscat (i,1) 
                     as_prop (as_natoms) = cr_prop (i) 
                     as_dw (as_natoms) = cr_dw (cr_nscat) 
                     as_occ(as_natoms) = cr_occ(cr_nscat) 
                     DO j = 1, 3 
                     as_pos (j, as_natoms) = cr_pos (j, i) 
                     ENDDO 
                  ENDIF 
   11             CONTINUE 
                  IF(mk_infile_internal) THEN
                     temp_present   (mk_iatom) = .true.
                     temp_in_crystal(mk_iatom) = i
                  ENDIF
!                                                                       
!     --------If we are reading a molecule insert atom into current     
                  IF (mole_l_on .AND. .NOT. mk_infile_internal) then 
                     CALL mole_insert_current (cr_natoms, mole_num_curr) 
                     IF (ier_num.lt.0.and.ier_num.ne. - 49) then 
                        GOTO 999 
                     ENDIF 
                     cr_prop(cr_natoms) = IBSET(cr_prop(cr_natoms),PROP_MOLECULE)
                     cr_mole(cr_natoms) = mole_num_curr
                  ELSE             ! No molecule header, but explicit info on line
                     IF(NINT(werte(6))>0 .AND. NINT(werte(7))>0) THEN
                        iimole = NINT(werte(6)) + n_mole_old
                        CALL mole_insert_explicit(cr_natoms, iimole, NINT(werte(7)))
                        cr_prop(cr_natoms) = IBSET(cr_prop(cr_natoms),PROP_MOLECULE)
                        cr_mole(cr_natoms) = NINT(werte(6)) + n_mole_old
                     ENDIF
                  ENDIF 
               ENDIF  noblank
            ENDIF inside
         ENDIF is_mole
      ENDIF blank1
      GOTO 1000 
!                                                                       
    2 CONTINUE 
!
!   If the file was read from internal storage, we need to assign the molecules
!
mole_int: IF(mk_infile_internal) THEN
       CALL testfile_internal (infile, natoms, &        ! Get size of internal structure
              nscat, TEMP_MAX_MOLE, TEMP_MAX_TYPE, TEMP_MAX_ATOM)
       temp_num_mole = TEMP_MAX_MOLE
       temp_num_type = max(TEMP_MAX_TYPE, maxval(clu_mole_tab))
       temp_num_atom = TEMP_MAX_ATOM
       ALLOCATE ( temp_mole_len  (0:temp_num_mole), STAT = istatus) ! Allocate temporary
       ALLOCATE ( temp_mole_off  (0:temp_num_mole), STAT = istatus) ! molecule arrays
       ALLOCATE ( temp_mole_type (0:temp_num_mole), STAT = istatus)
       ALLOCATE ( temp_mole_char (0:temp_num_mole), STAT = istatus)
       ALLOCATE ( temp_mole_file (0:temp_num_mole), STAT = istatus)
       ALLOCATE ( temp_mole_dens (0:temp_num_mole), STAT = istatus)
       ALLOCATE ( temp_mole_biso (0:temp_num_type), STAT = istatus)
       ALLOCATE ( temp_mole_clin (0:temp_num_type), STAT = istatus)
       ALLOCATE ( temp_mole_cqua (0:temp_num_type), STAT = istatus)
       ALLOCATE ( temp_mole_fuzzy(0:temp_num_mole), STAT = istatus)
       ALLOCATE ( temp_mole_cont (0:temp_num_atom), STAT = istatus)
       temp_mole_len  (0:temp_num_mole)= 0        ! Allocate temporary
       temp_mole_off  (0:temp_num_mole)= 0        ! molecule arrays
       temp_mole_type (0:temp_num_mole)= 0       
       temp_mole_char (0:temp_num_mole)= 0       
       temp_mole_file (0:temp_num_mole)= ' '    
       temp_mole_dens (0:temp_num_mole)= 0.0
       temp_mole_biso (0:temp_num_type)= 0.0
       temp_mole_clin (0:temp_num_type)= 0.0
       temp_mole_cqua (0:temp_num_type)= 0.0
       temp_mole_fuzzy(0:temp_num_mole)= 0.0
       temp_mole_cont (0:temp_num_atom)= 0       
       CALL stru_internal_molecules(infile, TEMP_MAX_MOLE, TEMP_MAX_TYPE, & ! Read domain 
              TEMP_MAX_ATOM, temp_num_mole, temp_num_type, & ! molecules
              temp_num_atom, temp_mole_len, temp_mole_off, temp_mole_type, & ! into temp
              temp_mole_char,                                              &
              temp_mole_file, temp_mole_dens, temp_mole_biso,              &
              temp_mole_clin, temp_mole_cqua, temp_mole_fuzzy, temp_mole_cont)
       IF(MOLE_MAX_MOLE < mole_num_mole + temp_num_mole .or.  &    ! If necessary increase
          MOLE_MAX_TYPE < mole_num_type + temp_num_type .or.  &    ! size of crystal molecule
          MOLE_MAX_ATOM < mole_num_atom + temp_num_atom     ) THEN ! arrays
          n_mole = MAX(MOLE_MAX_MOLE , mole_num_mole + temp_num_mole)
          n_type = MAX(MOLE_MAX_TYPE , mole_num_type + temp_num_type)
          n_atom = MAX(MOLE_MAX_ATOM , mole_num_atom + temp_num_atom)
          call alloc_molecule( MOLE_MAX_GENE, MOLE_MAX_SYMM, n_mole, n_type, n_atom )
       ENDIF
       DO i=1, temp_num_mole ! Create a lookup table, atom is in molecule
          temp_mole_type(i) = clu_mole_tab(temp_mole_type(i), clu_current)  ! Replace by look up table
          DO j=1, temp_mole_len(i)
             k = temp_mole_cont(temp_mole_off(i)+j)
             temp_in_mole(k) = i
          ENDDO
       ENDDO
       DO i=1, temp_num_mole ! Create content of all molecules
          k  = mole_off(mole_num_mole) + mole_len(mole_num_mole) ! Offset of next molecule
          ii = 0
          DO j=1,mk_iatom    ! test all domain atoms
             IF(temp_present(j)) THEN  ! Domain atom is within the crystal
                IF(temp_in_mole(j)==i) THEN ! Domain atom is inside this molecule
                   ii = ii + 1         ! Increment no of atoms in this molecule
                   mole_len  (mole_num_mole+1) = mole_len  (mole_num_mole+1) + 1    ! Adjust length
                   mole_cont (k + mole_len  (mole_num_mole+1)) = temp_in_crystal(j) ! Set content
                   cr_prop (temp_in_crystal(j)) = IBSET (cr_prop (temp_in_crystal(j)), PROP_MOLECULE) 
                   cr_mole(temp_in_crystal(j)) = mole_num_mole+1
                ENDIF
             ENDIF
          ENDDO
          IF(ii> 0) THEN ! If any atom is inside current molecule, then set all parameters
             mole_num_mole             = mole_num_mole + 1
             mole_num_type             = MAX(mole_num_type,temp_mole_type(i))
             mole_off  (mole_num_mole) = k
             mole_type (mole_num_mole) = temp_mole_type (i)
             mole_char (mole_num_mole) = temp_mole_char (i)
             mole_file (mole_num_mole) = temp_mole_file (i)
             mole_dens (mole_num_mole) = temp_mole_dens (i)
             mole_biso (mole_type(mole_num_mole)) = temp_mole_biso (temp_mole_type(i))
             mole_clin (mole_type(mole_num_mole)) = temp_mole_clin (temp_mole_type(i))
             mole_cqua (mole_type(mole_num_mole)) = temp_mole_cqua (temp_mole_type(i))
             mole_fuzzy(mole_num_mole) = temp_mole_fuzzy(i)
          ENDIF
       ENDDO
       mole_num_atom = mole_off(mole_num_mole) + mole_len(mole_num_mole) ! set total no of atoms in the molecules
!
!      Deallocate temporary arrays
!
       DEALLOCATE ( temp_mole_len  , STAT = istatus)
       DEALLOCATE ( temp_mole_off  , STAT = istatus)
       DEALLOCATE ( temp_mole_type , STAT = istatus)
       DEALLOCATE ( temp_mole_char , STAT = istatus)
       DEALLOCATE ( temp_mole_file , STAT = istatus)
       DEALLOCATE ( temp_mole_dens , STAT = istatus)
       DEALLOCATE ( temp_mole_biso , STAT = istatus)
       DEALLOCATE ( temp_mole_clin , STAT = istatus)
       DEALLOCATE ( temp_mole_cqua , STAT = istatus)
       DEALLOCATE ( temp_mole_fuzzy, STAT = istatus)
       DEALLOCATE ( temp_mole_cont , STAT = istatus)
       DEALLOCATE ( temp_present   , STAT = istatus)
       DEALLOCATE ( temp_in_crystal, STAT = istatus)
       DEALLOCATE ( temp_in_mole   , STAT = istatus)
    ENDIF mole_int
!                                                                       
!     --Remove old atoms inside micromdomain if microdomain is fuzzy    
!DBG_RBN      This should becom standard for all domains!!!!!!!!!!!!    
!     -- Removal is only performed, if md_sep_fuz is > 0
!                                                                       
!DBG      if(mc_type.eq.MD_DOMAIN_FUZZY) then                   
      remove_strict: IF ( clu_remove_mode <  CLU_REMOVE_TRUST ) THEN
         IF ( md_sep_fuz > 0.00 ) THEN   ! Separation is > 0, remove atoms
            DO i = 1, 3 
               mk_dim (i, 1) = mk_dim (i, 1) - mc_matrix (i, 4) 
               mk_dim (i, 2) = mk_dim (i, 2) - mc_matrix (i, 4) 
            ENDDO 
!                                                                       
            mk_natoms = cr_natoms - natoms_old 
            vv (1) = mc_matrix (1, 4) 
            vv (2) = mc_matrix (2, 4) 
            vv (3) = mc_matrix (3, 4) 
            IF (MAXSCAT  > ubound(surf_in_dist,1)) THEN
               new_nscat = MAXSCAT
               CALL alloc_surf   (new_nscat)
            ENDIF 
            CALL micro_fuzzy_rem (mk_dim, natoms_old, mk_natoms, md_sep_fuz,  &
            vv, shortest, mc_type, MD_DOMAIN_FUZZY)                                                     
            IF (shortest.gt.1.5 * md_sep_fuz) then 
               ier_num = + 1 
               ier_typ = ER_APPL 
               WRITE (ier_msg (2), 3100) shortest 
               WRITE (ier_msg (3), 3200) md_sep_fuz 
               CALL errlist 
            ENDIF 
         ENDIF 
      ELSE  remove_strict
         clu_remove_dist = MAX(clu_remove_dist, md_sep_fuz)
      ENDIF remove_strict
!
      IF (ier_num.eq. - 49) then 
         CALL no_error 
      ENDIF 
!                                                                       
  999 CONTINUE 
      CLOSE (ist) 
      IF(n_read_atoms == 0) THEN
         ier_num = -127
         ier_typ = ER_APPL
         ier_msg(1) = 'There were no atoms in the content file'
         ier_msg(2)(1:44) = infile(1:44)
      ENDIF
!                                                                       
 2000 FORMAT    (a) 
 3000 FORMAT    (a4, 4(f16.8, ','), 3(i8, ','), f16.8)
 3100 FORMAT    ('Shortest Distance to host ',f12.4) 
 3200 FORMAT    ('Intended fuzzy distance   ',f12.4) 
!
END SUBROUTINE micro_read_atom                
!
!*****7*****************************************************************
!
SUBROUTINE micro_read_simple (imd, lend, l_ok,infile, mc_dimen,  &
      mc_idimen, mc_matrix, MK_MAX_SCAT,mk_at_lis)                                             
!
!     Reads a single line from the input file. This should contain a 
!     pseudo atom that is interpreted as cluster type
!                                                                       
USE discus_config_mod 
USE domain_mod 
USE domaindis_mod 
USE read_internal_mod 
USE errlist_mod 
USE lib_random_func
use matrix_mod
USE string_convert_mod
USE precision_mod
USE str_comp_mod
!
IMPLICIT none 
!                                                                       
!                                                                       
INTEGER                                    , INTENT(IN)  :: imd 
LOGICAL                                    , INTENT(OUT) :: lend 
LOGICAL                                    , INTENT(OUT) :: l_ok 
CHARACTER (LEN=*)                          , INTENT(OUT) :: infile 
REAL(kind=PREC_DP),DIMENSION(4,4)          , INTENT(OUT) :: mc_dimen
REAL(kind=PREC_DP),DIMENSION(4,4)          , INTENT(OUT) :: mc_idimen
REAL(kind=PREC_DP),DIMENSION(4,4)          , INTENT(OUT) :: mc_matrix
INTEGER                                    , INTENT(IN)  :: MK_MAX_SCAT
CHARACTER (LEN=*), DIMENSION(0:MK_MAX_SCAT), INTENT(IN)  :: mk_at_lis
!                                                                       
CHARACTER(LEN=PREC_STRING) :: line
!                                                                       
INTEGER  :: i, j, ii 
REAL(kind=PREC_DP), dimension(3)     :: xyz (3) 
REAL(kind=PREC_DP) :: size_sigma
INTEGER  :: dummy_iscat
INTEGER  :: dummy_prop
INTEGER  :: dummy_mole
INTEGER  :: dummy_moleatom
INTEGER, DIMENSION(0:3)  :: dummy_surf
REAL(kind=PREC_DP)   , DIMENSION(0:3)  :: dummy_magn
!                                                                       
!                                                                       
ii = 0
clu_current = 0
l_ok = .false.                 ! Assume pseudoatom is incorrect
IF(clu_infile_internal) THEN   ! Read pseudo atom from internal storage
   clu_iatom = clu_iatom + 1   ! Increment internal atom number
   CALL struc_read_one_atom_internal(clu_infile, clu_iatom,  &
        xyz, dummy_iscat, dummy_prop, dummy_surf, dummy_magn, dummy_mole, dummy_moleatom )
   IF( ier_num == -105 ) THEN  ! read "end_of_file" 
      lend    = .true.
      ier_num = 0
      ier_typ = ER_NONE
      RETURN
   ENDIF
   WRITE(line, 1000) mk_at_lis(dummy_iscat), xyz ! copy into line
!if(clu_infile=='internal/STRU/cubo.0002.0001.pt') goto 1234
ELSE
  READ (imd, '(a)', end = 999) line 
ENDIF
IF(line.ne.' '.and.line (1:1) .ne.'#'.and.line(1:1)/='!'.AND.line.ne.char(13)) THEN
   CALL do_cap(line(1:4))
   ii = 0 
   DO i = 1, clu_number 
      IF (str_comp (line (1:4), clu_name (i), 2, 4, 4) ) then 
         ii = i 
      ENDIF 
   ENDDO 
!read(*,*) i
!                                                                       
   IF (ii.eq.0) then 
      ier_num = - 91 
      ier_typ = ER_APPL 
      ier_msg(1)  = ' Unknown pseudo atom name '//line(1:4)
      RETURN 
   ENDIF 
!                                                                       
   infile = clu_content (ii) 
   mc_type    = clu_character (ii) 
   md_sep_fuz = clu_fuzzy (ii) 
   READ (line (5:52), * ) xyz 
!                                                                       
   DO i = 1, 3 
      size_sigma = MAX(1. + gasdev(DBLE(clu_sigma(ii,i))), 0.01D0)
      DO j = 1, 4 
         mc_dimen (i, j)  = MAX(clu_shape (ii, i, j) * size_sigma , 0.001D0 )
         mc_idimen (i, j) = MAX(clu_shape (ii, i, j) * size_sigma , 0.001D0 )
         mc_matrix (i, j) = clu_orient (ii, i, j) 
      ENDDO 
      mc_dimen (i, 4) = clu_shape (ii, i, 4) + xyz (i) 
      mc_idimen (i, 4) = clu_shape (ii, i, 4) + xyz (i) 
      mc_matrix (i, 4) = clu_orient (ii, i, 4) + xyz (i) 
      mc_dimen (4, i) = 0.0 
      mc_idimen (4, i) = 0.0 
      mc_matrix (4, i) = 0.0 
   ENDDO 
   mc_dimen (4, i) = 1.0 
   mc_idimen (4, i) = 1.0 
   mc_matrix (4, i) = 1.0 
!  CALL invmat4 (mc_idimen) 
   CALL matinv  (mc_dimen, mc_idimen) 
   lend = .false. 
   l_ok = .true.       ! Pseudoatom is fine
ENDIF 
!DBG                                                                    
!DBG      do i=1,4                                                      
!DBG      write(*,2222) 'MC_DIMEN  ',(mc_dimen (i,j),j=1,4)             
!DBG      ENDDO                                                         
!DBG      do i=1,4                                                      
!DBG      write(*,2222) 'MC_IDIMEN ',(mc_idimen(i,j),j=1,4)             
!DBG      ENDDO                                                         
!DBG      do i=1,4                                                      
!DBG      write(*,2222) 'MC_MATRIX ',(mc_matrix(i,j),j=1,4)             
!DBG      ENDDO                                                         
!DBG2222      format(a10,3(f8.3,2x),2x,f8.3)                            
clu_current = ii                       ! Transfer current cluster type
!write(*,*) ' CLU_IATOM', clu_iatom, '>>',infile(1:len_trim(infile)),'<<', dummy_iscat, ier_num, ier_typ
RETURN 
!                                                                       
999 CONTINUE 
lend = .true. 
!
1000  FORMAT(a4, 3(f16.8))
!                                                                       
END SUBROUTINE micro_read_simple              
!
!*****7*****************************************************************
!
      SUBROUTINE micro_fuzzy_rem (mk_dim, natoms_old, mk_natoms,        &
      md_sep_fuz, w, shortest, mc_type, MD_DOMAIN_FUZZY)                                          
!-                                                                      
!     Removes all atoms up to number natoms_old that are inside         
!     the microdomain of fuzzy boundary type                            
!+                                                                      
      USE discus_config_mod 
      USE crystal_mod 
      USE metric_mod
      USE molecule_mod 
      USE prop_para_mod 
      USE surface_mod 
!
use precision_mod
!
      IMPLICIT none 
!                                                                       
      REAL(kind=PREC_DP), DIMENSION(3,2), INTENT(IN) ::  mk_dim
      INTEGER             , INTENT(in) ::  natoms_old
      INTEGER             , INTENT(IN) ::  mk_natoms 
      REAL(kind=PREC_DP)  , INTENT(IN) ::  md_sep_fuz 
      REAL(kind=PREC_DP), DIMENSION(3)  , INTENT(IN) ::  w
      REAL(kind=PREC_DP)  , INTENT(OUT)::  shortest 
      INTEGER             , INTENT(IN) ::  mc_type
      INTEGER             , INTENT(IN) ::  MD_DOMAIN_FUZZY
       
!                                                                       
      INTEGER i, j, k 
      LOGICAL lspace 
REAL(kind=PREC_DP), dimension(3) :: u (3), v (3)
      REAL(kind=PREC_DP) :: a, b, c , d
      REAL(kind=PREC_DP) :: distance 
      REAL(kind=PREC_DP) :: separation
!                                                                       
      shortest = 1.0e03 
      separation = 1.0e03
      lspace = .true. 
!
      type_fuzzy: IF(mc_type.eq.MD_DOMAIN_FUZZY) THEN  !This is a fuzzy domain do fast loop
         fuzzy_main: DO i = 1, natoms_old 
            separation = 1.0e03 
            fuzzy_is_void: IF (cr_iscat (i,1) .ne.0) then 
               DO k = 1, 3 
                  u (k) = cr_pos (k, i) 
               ENDDO 
               fuzzy_inner: DO j = cr_natoms - mk_natoms + 1, cr_natoms 
            DO k = 1, 3 
            v (k) = cr_pos (k, j) 
            ENDDO 
            distance = do_blen (lspace, u, v) 
            IF (distance.le.md_sep_fuz) then 
!                                                                       
!     ----------Atom is too close to microdomain atom, remove.          
!                                                                       
               cr_iscat (i,1) = 0 
               cr_prop (i) = IBCLR (cr_prop (i), PROP_NORMAL) 
               cr_prop (i) = IBSET (cr_prop (i), PROP_DOMAIN) 
               shortest = distance 
            ELSE 
               shortest = min (shortest, distance) 
               separation = min (separation, distance) 
            ENDIF 
               ENDDO fuzzy_inner
            ENDIF fuzzy_is_void
         ENDDO fuzzy_main
         IF (separation.lt.surf_in_dist (cr_iscat (i,1) ) ) then 
            IF (BTEST (cr_prop (i), PROP_NORMAL) ) then 
               cr_prop (i) = IBSET (cr_prop (i), PROP_SURFACE_INT) 
            ENDIF 
         ENDIF 
      ELSE type_fuzzy ! Other domain types (Should probalbly changed as well...)
!                                                                       
!     Set up extra space around microdomain to allow for outside        
!     influence of separation                                           
!                                                                       
      d = MAX(md_sep_fuz, MAXVAL(surf_in_dist))
      a = 1.5 * d / cr_a0 (1) 
      b = 1.5 * d / cr_a0 (2) 
      c = 1.5 * d / cr_a0 (3) 
!                                                                       
!     Loop over all atoms previously inside the crystal                 
!                                                                       
      DO i = 1, natoms_old 
      separation = 1.0e03 
      IF (cr_iscat (i,1) .ne.0) then 
!                                                                       
!     ----Check if the atom is inside a domain or object,               
!         keep atom if inside a domain or object                        
         DO j = 1, mole_num_mole 
         IF (mole_char (j) .lt.MOLE_ATOM) then 
            DO k = 1, mole_len (j) 
            IF (mole_cont (mole_off (j) + k) .eq.i) then 
               GOTO 999 
            ENDIF 
            ENDDO 
         ENDIF 
         ENDDO 
         DO k = 1, 3 
         u (k) = cr_pos (k, i) 
         ENDDO 
         IF (mk_dim (1, 1) + w (1) - a.le.cr_pos (1, i)             .and. &
             cr_pos (1,  i)           .le.mk_dim (1, 2) + w (1) + a .and. &
             mk_dim (2, 1) + w (2) - b.le.cr_pos (2, i)             .and. &
             cr_pos (2, i)            .le.mk_dim (2, 2) + w (2) + b .and. &
             mk_dim (3, 1) + w (3) - c.le.cr_pos (3, i)             .and. &
             cr_pos (3, i)            .le.mk_dim (3, 2) + w (3) + c) then         
            DO j = cr_natoms - mk_natoms + 1, cr_natoms 
            DO k = 1, 3 
            v (k) = cr_pos (k, j) 
            ENDDO 
            distance = do_blen (lspace, u, v) 
            IF (distance.le.md_sep_fuz) then 
!                                                                       
!     ----------Atom is too close to microdomain atom, remove.          
!                                                                       
               cr_iscat (i,1) = 0 
               cr_prop (i) = IBCLR (cr_prop (i), PROP_NORMAL) 
               cr_prop (i) = IBSET (cr_prop (i), PROP_DOMAIN) 
               shortest = distance 
            ELSE 
               shortest = min (shortest, distance) 
               separation = min (separation, distance) 
            ENDIF 
            ENDDO 
         ENDIF 
      ENDIF 
  999 CONTINUE 
      IF (separation.lt.surf_in_dist (cr_iscat (i,1) ) ) then 
         IF (BTEST (cr_prop (i), PROP_NORMAL) ) then 
            cr_prop (i) = IBSET (cr_prop (i), PROP_SURFACE_INT) 
         ENDIF 
      ENDIF 
      ENDDO 
!
      ENDIF type_fuzzy
!                                                                       
!     DO j = cr_natoms - mk_natoms + 1, cr_natoms 
!     DO k = 1, 3 
!     v (k) = cr_pos (k, j) 
!     ENDDO 
!     separation = 1.0e03 
!     DO i = 1, natoms_old 
!     IF (cr_iscat (i,1) .ne.0) then 
!        IF (mk_dim (1, 1) + w (1) - a.le.cr_pos (1, i) .and.cr_pos (1, &
!        i) .le.mk_dim (1, 2) + w (1) + a.and.mk_dim (2, 1) + w (2)     &
!        - b.le.cr_pos (2, i) .and.cr_pos (2, i) .le.mk_dim (2, 2)      &
!        + w (2) + b.and.mk_dim (3, 1) + w (3) - c.le.cr_pos (3, i)     &
!        .and.cr_pos (3, i) .le.mk_dim (3, 2) + w (3) + c) then         
!           DO k = 1, 3 
!           u (k) = cr_pos (k, i) 
!           ENDDO 
!           distance = do_blen (lspace, u, v) 
!           separation = min (separation, distance) 
!        ENDIF 
!     ENDIF 
!     ENDDO 
!     IF (separation.lt.surf_in_dist (cr_iscat (j,1) ) ) then 
!        IF (BTEST (cr_prop (j), PROP_NORMAL) ) then 
!           cr_prop (j) = IBSET (cr_prop (j), PROP_SURFACE_INT) 
!        ENDIF 
!     ENDIF 
!     ENDDO 
                                                                        
!                                                                       
      END SUBROUTINE micro_fuzzy_rem                
!*****7*****************************************************************
      SUBROUTINE domain_do_set (zeile, length) 
!+                                                                      
!     This subroutine sets various parameters                           
!-                                                                      
      USE discus_config_mod 
      USE domain_mod
      USE modify_mod
      USE surface_func_mod
      USE errlist_mod 
      USE get_params_mod
USE precision_mod
USE str_comp_mod
      IMPLICIT none 
!                                                                       
!                                                                       
      INTEGER maxw 
      PARAMETER (maxw = 20) 
!                                                                       
      CHARACTER ( * ) zeile 
      INTEGER length 
!                                                                       
      CHARACTER(LEN=PREC_STRING) :: cpara (maxw) 
      INTEGER lpara (maxw) 
      INTEGER ianz 
      REAL(KIND=PREC_DP) :: werte (maxw) 
!                                                                       
!                                                                       
      CALL get_params (zeile, ianz, cpara, lpara, maxw, length) 
      IF (ier_num.ne.0) return 
      IF (ianz.le.0) return 
!                                                                       
      IF (str_comp (cpara (1) , 'distance', 2, lpara (1) , 8) ) then 
         CALL del_params (1, ianz, cpara, lpara, maxw) 
         CALL surf_set_fuzzy (ianz, cpara, lpara, werte, maxw, 1) 
      ELSEIF(str_comp (cpara (1) , 'remove', 2, lpara (1) , 2) ) then
         IF    (str_comp (cpara (2) , 'initial', 2, lpara (2) , 6) ) then
            clu_remove_mode = CLU_REMOVE_INITIAL
         ELSEIF(str_comp (cpara (2) , 'strict', 2, lpara (2) , 6) ) then
            clu_remove_mode = CLU_REMOVE_STRICT
         ELSEIF(str_comp (cpara (2) , 'trust', 2, lpara (2) , 5) ) then
            clu_remove_mode = CLU_REMOVE_TRUST
         ELSEIF(str_comp (cpara (2) , 'none', 2, lpara (2) , 5) ) then
            clu_remove_mode = CLU_REMOVE_NONE
         ELSE 
            ier_num = - 6 
            ier_typ = ER_COMM 
         ENDIF 
      ELSE 
         ier_num = - 6 
         ier_typ = ER_COMM 
      ENDIF 
!                                                                       
      END SUBROUTINE domain_do_set                  
!
!*******************************************************************************
!
SUBROUTINE domain_reset
!
USE discus_allocate_appl_mod
USE domain_mod
USE micro_mod
!
IMPLICIT NONE
!
CALL alloc_domain(1)
!
mk_name  = ' '
mk_spcgr = 'P1'
!
mk_a0(:)    = (/1.0, 1.0, 1.0/)
mk_win(:)   = (/90.0, 90.0, 90.0/)
mk_dim(:,:) = 0.0
mk_nscat    = 0
mk_natoms   = 0
!
mk_spcgr_ianz = 1
mk_spcgr_para = 0
mk_spcgr_no   = 1
mk_gen_add_n  = 0
mk_gen_add_power(:) = 1
mk_gen_add(:,:,:  ) = &
     RESHAPE((/(1.,(0.,0.,0.,0.,1.,ik=1,3),il=0,mk_GEN_ADD_MAX)/),(/4,4,mk_GEN_ADD_MAX+1/))
!
mk_sym_add_n                     = 0
mk_sym_add_power(:) = 1
mk_sym_add(:,:,:)   = &
     RESHAPE((/(1.,(0.,0.,0.,0.,1.,ik=1,3),il=0,mk_SYM_ADD_MAX)/),(/4,4,mk_SYM_ADD_MAX+1/))
!
mk_infile_internal = .FALSE.  ! File with domain atoms in stored internally
mk_iatom           = 1        ! Current atom number in internal file
!
IF(ALLOCATED(clu_content))    clu_content(:) = ' '   ! (CLU_MAX_TYPE)
IF(ALLOCATED(clu_name))       clu_name(:)    = ' '      ! (CLU_MAX_TYPE)
!
IF(ALLOCATED(clu_character))  clu_character(:) = -4  ! (CLU_MAX_TYPE)
IF(ALLOCATED(clu_fuzzy))      clu_fuzzy(:)     = 2.55  ! (CLU_MAX_TYPE)
IF(ALLOCATED(clu_orient))     clu_orient(:,:,:) = 0.00 ! (CLU_MAX_TYPE,3,4)
IF(ALLOCATED(clu_shape))      clu_shape(:,:,:) = 0.00  ! (CLU_MAX_TYPE,3,4)
IF(ALLOCATED(clu_sigma))      clu_sigma(:,:  ) = 0.00  ! (CLU_MAX_TYPE,3  )
!
clu_infile = ' '
clu_index = 0
clu_mode  = CLU_IN_PSEUDO
clu_number = 0
!
clu_surface = .FALSE.
clu_infile_internal = .false. ! Is infile an internal file ?
clu_iatom = 0
!
clu_remove_mode = CLU_REMOVE_STRICT   ! Remove initial atoms or include prev domains
clu_remove_end  = 0  ! Initial last atom to remove
clu_remove_dist = 0.0! Initial removal distance
!
END SUBROUTINE domain_reset
!
!*******************************************************************************
!
END MODULE domain_menu
