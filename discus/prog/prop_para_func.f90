MODULE prop_para_func
!
PRIVATE
!
PUBLIC property_menu
PUBLIC property_select
PUBLIC property_set_user
!
CONTAINS
!
!*****7*****************************************************************
!
SUBROUTINE property_menu 
!-                                                                      
!     Main menu for property related operations                         
!+                                                                      
USE discus_config_mod 
USE crystal_mod 
USE molecule_mod 
USE prop_para_mod 
USE surface_func_mod
!
      USE calc_expr_mod
      USE doact_mod 
      USE do_eval_mod
      USE do_wait_mod
      USE errlist_mod 
      USE learn_mod 
      USE class_macro_internal
      USE prompt_mod 
      USE sup_mod
      IMPLICIT none 
!                                                                       
       
!                                                                       
      INTEGER maxw 
!                                                                       
      CHARACTER(5) befehl 
      CHARACTER(LEN=LEN(prompt)) :: orig_prompt
      CHARACTER(1024) line, zeile
      INTEGER lp, length, lbef 
      INTEGER indxg
      LOGICAL lend
!                                                                       
      INTEGER len_str 
      LOGICAL str_comp 
!                                                                       
      maxw = MAXSCAT
      lend = .false. 
      CALL no_error 
      orig_prompt = prompt
      prompt = prompt (1:len_str (prompt) ) //'/prop' 
!                                                                       
      DO while (.not.lend) 
      CALL get_cmd (line, length, befehl, lbef, zeile, lp, prompt) 
      IF (ier_num.eq.0) then 
         IF (line /= ' '      .and. line(1:1) /= '#' .and. &
             line /= char(13) .and. line(1:1) /= '!'        ) THEN
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
!     ------evaluatean expression and assign the value to a variabble   
!                                                                       
               CALL do_math (line, indxg, length) 
            ELSE 
!                                                                       
!------ ----execute a macro file                                        
!                                                                       
               IF (befehl (1:1) .eq.'@') then 
                  IF (length.ge.2) then 
                     CALL file_kdo (line (2:length), length - 1) 
                  ELSE 
                     ier_num = - 13 
                     ier_typ = ER_MAC 
                  ENDIF 
!                                                                       
!     ----continues a macro 'continue'                                  
!                                                                       
               ELSEIF (str_comp (befehl, 'continue', 2, lbef, 8) ) then 
                  CALL macro_continue (zeile, lp) 
!                                                                       
!------ ----Echo a string, just for interactive check in a macro 'echo' 
!                                                                       
               ELSEIF (str_comp (befehl, 'echo', 2, lbef, 4) ) then 
                  CALL echo (zeile, lp) 
!                                                                       
!      ---Evaluate an expression, just for interactive check 'eval'     
!                                                                       
               ELSEIF (str_comp (befehl, 'eval', 2, lbef, 4) ) then 
                  CALL do_eval (zeile, lp, .TRUE.) 
!                                                                       
!     ----exit 'exit'                                                   
!                                                                       
               ELSEIF (str_comp (befehl, 'exit', 2, lbef, 4) ) then 
                  lend = .true. 
!                                                                       
!     ----help 'help','?'                                               
!                                                                       
      ELSEIF (str_comp (befehl, 'help', 2, lbef, 4) .or.str_comp (befehl&
     &, '?   ', 1, lbef, 4) ) then                                      
                  IF (str_comp (zeile, 'errors', 2, lp, 6) ) then 
                     lp = lp + 7 
                     CALL do_hel ('discus '//zeile, lp) 
                  ELSE 
                     lp = lp + 16 
                     CALL do_hel ('discus property '//zeile, lp) 
                  ENDIF 
!                                                                       
!------- -Operating System Kommandos 'syst'                             
!                                                                       
               ELSEIF (str_comp (befehl, 'syst', 2, lbef, 4) ) then 
                  IF (zeile.ne.' ') then 
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
!     ----Set/clear properties for types or (individual) atoms 'set'         
!
               ELSEIF (str_comp (befehl, 'clear', 2, lbef, 5) ) then 
                  CALL property_set_clr (zeile, lp, .false.) 
!
!     ----Set/clear properties for types or (individual) atoms 'set'         
!
               ELSEIF (str_comp (befehl, 'set', 2, lbef, 3) ) then 
                  CALL property_set_clr (zeile, lp, .true.) 
!                                                                       
!                                                                       
!     ----Define which properties have to be present 'property'         
!                                                                       
               ELSEIF (str_comp (befehl, 'property', 2, lbef, 8) ) then 
                  CALL property_select (zeile, lp, cr_sel_prop) 
!                                                                       
!     ----show current parameters 'show'                                
!                                                                       
               ELSEIF (str_comp (befehl, 'show', 2, lbef, 4) ) then 
                  CALL property_show 
               ELSE 
                  ier_num = - 8 
                  ier_typ = ER_COMM 
               ENDIF 
            ENDIF 
         ENDIF 
      ENDIF 
      IF (ier_num.ne.0) THEN 
         CALL errlist 
         IF (ier_sta.ne.ER_S_LIVE) THEN 
            IF (lmakro .OR. lmakro_error) THEN  ! Error within macro or termination errror
               IF(sprompt /= prompt ) THEN
                  ier_num = -10
                  ier_typ = ER_COMM
                  ier_msg(1) = ' Error occured in property menu'
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
      ENDDO 
!
      prompt = orig_prompt
!                                                                       
      END SUBROUTINE property_menu                  
!
!*****7*****************************************************************
!
      SUBROUTINE property_set_clr(zeile, lp, set_clr)
!+                                                                      
!     This subroutine sets or clears the property flag of an atom or atom type
!-                                                                      
      USE discus_config_mod 
      USE crystal_mod 
      USE get_iscat_mod
      USE prop_para_mod 
!
      USE ber_params_mod
      USE errlist_mod 
      USE get_params_mod
USE precision_mod
      IMPLICIT none 
!
      CHARACTER (LEN=*), INTENT(INOUT) :: zeile
      INTEGER          , INTENT(INOUT) :: lp
      LOGICAL          , INTENT(IN   ) :: set_clr

      INTEGER, PARAMETER                   :: TYPES = 0
      INTEGER, PARAMETER                   :: ATOMS = 1
      INTEGER, PARAMETER                   :: MAXW = 200 
      LOGICAL, PARAMETER                   :: LOLD = .true.
      CHARACTER(LEN=1024), DIMENSION(MAXW) :: cpara
      INTEGER            , DIMENSION(MAXW) :: lpara
      REAL(KIND=PREC_DP) , DIMENSION(MAXW) :: werte
      INTEGER                              :: ianz
      INTEGER                              :: i,is
      INTEGER                              :: ibit_nr = PROP_NORMAL
      INTEGER                              :: sel_mode
      INTEGER                              :: istart, iend
      LOGICAL, DIMENSION(:), ALLOCATABLE   :: latom
!
      LOGICAL str_comp
!
      CALL get_params (zeile, ianz, cpara, lpara, MAXW, lp) 
      IF (ier_num.ne.0) return 
      IF(str_comp(cpara(1)(1:lpara(1)),'domain', 3, lpara(1), 6)) THEN
         ibit_nr = PROP_DOMAIN
      ELSEIF(str_comp(cpara(1)(1:lpara(1)),'outside', 3, lpara(1), 7)) THEN
         ibit_nr = PROP_OUTSIDE
      ELSEIF(str_comp(cpara(1)(1:lpara(1)),'external', 3, lpara(1), 7)) THEN
         ibit_nr = PROP_SURFACE_EXT
      ELSEIF(str_comp(cpara(1)(1:lpara(1)),'internal', 3, lpara(1), 7)) THEN
         ibit_nr = PROP_SURFACE_INT
      ELSEIF(str_comp(cpara(1)(1:lpara(1)),'ligand', 3, lpara(1), 6)) THEN
         ibit_nr = PROP_LIGAND
      ELSE
         ier_num = -6
         ier_typ = ER_COMM
         ier_msg(1) = 'The 1st parameter must be either of '
         ier_msg(2) = '''domain'', ''outside'', ''external'', ''internal'' '
         ier_msg(3) = '''ligand'' '
         RETURN
      ENDIF
      IF(str_comp(cpara(2)(1:lpara(2)),'types', 3, lpara(2), 5)) THEN
         sel_mode = TYPES
      ELSEIF(str_comp(cpara(2)(1:lpara(2)),'atoms', 3, lpara(2), 5)) THEN
         sel_mode = ATOMS
      ELSE
         ier_num = -6
         ier_typ = ER_COMM
         RETURN
      ENDIF
      CALL del_params (2, ianz, cpara, lpara, MAXW) 
      IF (ier_num.ne.0) return 
      ALLOCATE(latom(0:MAXSCAT))
      latom(:) =.FALSE.
      typesel: IF(sel_mode == TYPES) THEN
         CALL get_iscat (ianz, cpara, lpara, werte, maxw, lold) 
         IF (ier_num.ne.0) return 
         IF(NINT(werte(1))==-1) THEN
            latom = .true.
         ELSE
            latom = .false.
            DO i = 1, ianz 
               is = nint (werte (i) ) 
               IF (is.ge.0.and.is.le.cr_nscat) then 
                  latom (is) = .true. 
               ELSE 
                  ier_num = - 27 
                  ier_typ = ER_APPL 
               ENDIF 
            ENDDO 
         ENDIF
         IF(set_clr) THEN
            DO i=1, cr_natoms
               IF(latom(cr_iscat(i))) THEN
                  cr_prop(i) = IBSET(cr_prop(i),ibit_nr)
               ENDIF
            ENDDO
         ELSE
            DO i=1, cr_natoms
               IF(latom(cr_iscat(i))) THEN
                  cr_prop(i) = IBCLR(cr_prop(i),ibit_nr)
               ENDIF
            ENDDO
         ENDIF
      ELSEIF(sel_mode == ATOMS) THEN typesel
         IF(str_comp(cpara(1)(1:lpara(1)),'all', 3, lpara(1), 3)) THEN
             istart=1
             iend  = cr_natoms
         ELSE
            CALL ber_params (ianz, cpara, lpara, werte, maxw) 
            IF (ier_num.ne.0) return 
            IF(ianz==1) THEN
               istart=nint(werte(1))
               iend  =nint(werte(1))
            ELSEIF(ianz==2) THEN
               istart=nint(werte(1))
               iend  =nint(werte(2))
            ELSE
               ier_num = -6
               ier_typ = ER_COMM
               RETURN
            ENDIF
         ENDIF
         IF(set_clr) THEN
            DO i=istart, iend
               cr_prop(i) = IBSET(cr_prop(i),ibit_nr)
            ENDDO
         ELSE
            DO i=istart, iend
               cr_prop(i) = IBCLR(cr_prop(i),ibit_nr)
            ENDDO
         ENDIF
      ENDIF typesel
!
      DEALLOCATE(latom)
!
      END SUBROUTINE property_set_clr
!
!*****7*****************************************************************
!
      SUBROUTINE property_show 
!+                                                                      
!     This subroutine shows the property settings                       
!-                                                                      
      USE discus_config_mod 
      USE crystal_mod 
      USE prop_char_mod 
      USE prop_para_mod 
USE surface_func_mod
!
      USE errlist_mod 
      USE prompt_mod 
      IMPLICIT none 
!                                                                       
       
!                                                                       
      CHARACTER(32) c_property 
      INTEGER length 
!                                                                       
      WRITE (output_io, 2000) 
      CALL char_prop_2 (c_property, cr_sel_prop (1), cr_sel_prop (0),   &
      length)                                                           
      WRITE (output_io, 2100) c_property (1:length) 
!
      CALL property_show_user
!                                                                       
      CALL surf_show 
!                                                                       
 2000 FORMAT     (/' Property related settings'/) 
                                                                        
 2100 FORMAT    (/' Atom properties              ',/,                   &
     &                   '   N = normal atom            ',/,            &
     &                   '   M = atom in a molecule     ',/,            &
     &                   '   D = atom in a domain       ',/,            &
     &                   '   O = atom outside boundaries',/,            &
     &                   '   E = atom near ext. surface ',/,            &
     &                   '   I = atom near int. surface ',/,            &
     &                   '   L = atom in ligand molecule',/,            &
     &                   '                              : ','NMDOEIL'/, &
     &                   '      absent=- ignored=.      : ',a)          
!                                                                       
      END SUBROUTINE property_show                  
!
!*****7*****************************************************************
!
      SUBROUTINE property_select (line, length, sel_mask) 
!-                                                                      
!     Sets the property bits for a select/deselect/mselecet/mdeselect   
!                                                                       
      USE prop_para_mod 
      USE errlist_mod 
      USE get_params_mod
USE precision_mod
      USE take_param_mod
      IMPLICIT none 
!                                                                       
!                                                                       
      CHARACTER ( * ) line 
      INTEGER length 
      INTEGER sel_mask (0:1) 
!                                                                       
      INTEGER maxw 
      PARAMETER (maxw = 25) 
!                                                                       
      CHARACTER(1024) cpara (maxw) 
      INTEGER lpara (maxw) 
      INTEGER ianz 
!                                                                       
      INTEGER, PARAMETER :: NOPTIONAL = 7
      CHARACTER(LEN=1024), DIMENSION(NOPTIONAL) :: oname   !Optional parameter names
      CHARACTER(LEN=1024), DIMENSION(NOPTIONAL) :: opara   !Optional parameter strings returned
      INTEGER            , DIMENSION(NOPTIONAL) :: loname  !Lenght opt. para name
      INTEGER            , DIMENSION(NOPTIONAL) :: lopara  !Lenght opt. para name returned
      LOGICAL            , DIMENSION(NOPTIONAL) :: lpresent!opt. para present
      REAL(KIND=PREC_DP) , DIMENSION(NOPTIONAL) :: owerte   ! Calculated values
      INTEGER, PARAMETER                        :: ncalc = 5 ! Number of values to calculate 
!
      LOGICAL str_comp 
!
      DATA oname  / 'nmin' , 'nmax'  , 'emin'  , 'emax'  , 'no  ', 'atom', 'conn'   /
      DATA loname /  4     ,  4      ,  4      ,  4      ,  2    ,  4    ,  4       /
      opara  =  (/ '0.0000', '192.00', '-1.000', '-1.000', '0.0000', 'none  ', 'none  '   /)   ! Always provide fresh default values
      lopara =  (/  6      ,  6      ,  6      ,  6      ,  6      ,  4      ,  4         /)
      owerte =  (/  0.0,      192.0  ,  -1.0   ,  -1.0   ,  0.00   , -2.00   ,  0.00      /)
!
      CALL get_params (line, ianz, cpara, lpara, maxw, length) 
      CALL get_optional(ianz, MAXW, cpara, lpara, NOPTIONAL,  ncalc, &
                        oname, loname, opara, lopara, lpresent, owerte)
      IF (ier_num.ne.0) RETURN 
      IF (ianz >  1) THEN
!                                                                       
      IF (str_comp (cpara (1) , 'ignore', 2, lpara (1) , 6) ) THEN 
         CALL del_params (1, ianz, cpara, lpara, maxw) 
         CALL property_set (ianz, cpara, lpara, maxw, sel_mask, .false.,&
         .false.)                                                       
         IF(cpara(1)=='all') THEN
            opara(5) = '-1'
            lopara(5) = 2
            owerte(5) = -1.
            CALL property_set_user( 'ignore', opara, lopara, owerte, NOPTIONAL)
         ENDIF
      ELSEIF (str_comp (cpara (1) , 'present', 2, lpara (1) , 7) ) THEN 
         CALL del_params (1, ianz, cpara, lpara, maxw) 
         CALL property_set (ianz, cpara, lpara, maxw, sel_mask, .true., &
         .true.)                                                        
      ELSEIF (str_comp (cpara (1) , 'absent', 2, lpara (1) , 6) ) THEN 
         CALL del_params (1, ianz, cpara, lpara, maxw) 
         CALL property_set (ianz, cpara, lpara, maxw, sel_mask, .true., &
         .false.)                                                       
      ELSE 
         ier_num = - 6 
         ier_typ = ER_COMM 
      ENDIF 
      ELSE 
         line = cpara(1)(1:lpara(1))
         CALL property_set_user( cpara(1), opara, lopara, owerte, NOPTIONAL)
      ENDIF 
!                                                                       
      END SUBROUTINE property_select                
!
!*****7*****************************************************************
!
      SUBROUTINE property_set (ianz, cpara, lpara, maxw, sel_field,     &
      lentry1, lentry2)                                                 
!-                                                                      
!     Sets the property bits for a select/deselect/mselecet/mdeselect   
!                                                                       
      USE prop_para_mod 
      USE modify_func_mod
      USE errlist_mod 
      IMPLICIT none 
!                                                                       
!                                                                       
      INTEGER maxw 
      INTEGER ianz 
      CHARACTER ( * ) cpara (maxw) 
      INTEGER lpara (maxw) 
      INTEGER sel_field (0:1) 
      LOGICAL lentry1 
      LOGICAL lentry2 
!                                                                       
      INTEGER i, j 
!                                                                       
      LOGICAL str_comp 
!     INTEGER bit_set 
!                                                                       
      DO i = 1, ianz 
      IF (str_comp (cpara (i) , 'all', 2, lpara (i) , 3) ) then 
         DO j = 0, MAXPROP - 1 
         sel_field (0) = bit_set (sel_field, 0, j, lentry1) 
         sel_field (1) = bit_set (sel_field, 1, j, lentry2) 
         ENDDO 
      ELSEIF (str_comp (cpara (i) , 'normal', 2, lpara (i) , 6) ) then 
         sel_field (0) = bit_set (sel_field, 0, PROP_NORMAL, lentry1) 
         sel_field (1) = bit_set (sel_field, 1, PROP_NORMAL, lentry2) 
      ELSEIF (str_comp (cpara (i) , 'molecule', 2, lpara (i) , 8) )     &
      then                                                              
         sel_field (0) = bit_set (sel_field, 0, PROP_MOLECULE, lentry1) 
         sel_field (1) = bit_set (sel_field, 1, PROP_MOLECULE, lentry2) 
      ELSEIF (str_comp (cpara (i) , 'domain', 2, lpara (i) , 6) ) then 
         sel_field (0) = bit_set (sel_field, 0, PROP_DOMAIN, lentry1) 
         sel_field (1) = bit_set (sel_field, 1, PROP_DOMAIN, lentry2) 
      ELSEIF (str_comp (cpara (i) , 'outside', 2, lpara (i) , 7) ) then 
         sel_field (0) = bit_set (sel_field, 0, PROP_OUTSIDE, lentry1) 
         sel_field (1) = bit_set (sel_field, 1, PROP_OUTSIDE, lentry2) 
      ELSEIF (str_comp (cpara (i) , 'external', 2, lpara (i) , 8) )     &
      then                                                              
         sel_field (0) = bit_set (sel_field, 0, PROP_SURFACE_EXT,       &
         lentry1)                                                       
         sel_field (1) = bit_set (sel_field, 1, PROP_SURFACE_EXT,       &
         lentry2)                                                       
      ELSEIF (str_comp (cpara (i) , 'internal', 2, lpara (i) , 8) )     &
      then                                                              
         sel_field (0) = bit_set (sel_field, 0, PROP_SURFACE_INT,       &
         lentry1)                                                       
         sel_field (1) = bit_set (sel_field, 1, PROP_SURFACE_INT,       &
         lentry2)                                                       
      ELSEIF (str_comp (cpara (i) , 'ligand', 2, lpara(i) , 6) ) THEN
         sel_field (0) = bit_set (sel_field, 0, PROP_LIGAND, lentry1)
         sel_field (1) = bit_set (sel_field, 1, PROP_LIGAND, lentry2)
      ELSE 
         ier_num = - 6 
         ier_typ = ER_COMM 
         ier_msg(1) = 'Selection must be: ignore, present, absent'
         RETURN 
      ENDIF 
      ENDDO 
!                                                                       
      END SUBROUTINE property_set                   
!
!*******************************************************************************
!
SUBROUTINE property_set_user( act, opara, lopara, owerte, NOPTIONAL)
!
USE discus_config_mod
USE conn_sup_mod
USE get_iscat_mod
USE prop_para_mod
!
USE ber_params_mod
USE errlist_mod
USE precision_mod
!
IMPLICIT NONE
!
INTEGER, INTENT(IN)                                :: NOPTIONAL
CHARACTER(LEN=*)                      , INTENT(IN) :: act
CHARACTER(LEN=*), DIMENSION(NOPTIONAL), INTENT(IN) :: opara
INTEGER         , DIMENSION(NOPTIONAL), INTENT(IN) :: lopara
REAL(KIND=PREC_DP),DIMENSION(NOPTIONAL), INTENT(IN) :: owerte
!
TYPE(prop_templ), DIMENSION(:), ALLOCATABLE :: p_temp
TYPE(prop_templ)                            :: pp
!
INTEGER, PARAMETER :: MAXWW= 10
INTEGER            :: MAXW = 10
CHARACTER(LEN=1024), DIMENSION(MAX(MAXWW, MAXSCAT)) :: cpara
INTEGER            , DIMENSION(MAX(MAXWW, MAXSCAT)) :: lpara
REAL(KIND=PREC_DP) , DIMENSION(MAX(MAXWW, MAXSCAT)) :: werte
INTEGER                              :: ianz
!
INTEGER, DIMENSION(:), ALLOCATABLE :: at_kind    ! Atom type(s) on "atom:" parameter
!
CHARACTER (LEN=256)  :: c_name   ! Connectivity name
INTEGER              :: c_name_l ! connectivity name length
INTEGER              :: ino      ! connectivity number
INTEGER              :: i,j      ! index 
INTEGER              :: id       ! Propert identifier
INTEGER              :: nuser=0  ! Number of user definitions
LOGICAL              :: lsearch  ! Need a search Yes / No
LOGICAL              :: lall     ! Work on all entries No
LOGICAL              :: lactive  ! There are properties with act /= 0
!
MAXW = MAX(MAXWW, MAXSCAT)
id = NINT(owerte(5))
IF(opara(6)=='none'.AND.id==0) RETURN          ! No atom selected nothing to do
!
IF(ALLOCATED(prop_user)) THEN                    ! Old defs exist
   IF((prop_user_no>= UBOUND(prop_user,1) .OR.     &
       id >UBOUND(prop_user,1)                 )   &
      .AND. act/='ignore'                             ) THEN   ! Need more space
      nuser = MAX(nuser, prop_user_no+10, UBOUND(prop_user,1), id)
      ALLOCATE(p_temp(nuser))                                  ! Make temporary storage
      p_temp(1:prop_user_no) = prop_user(1:prop_user_no)       ! Restore old data
      CALL MOVE_ALLOC(FROM=p_temp, TO=prop_user)
   ENDIF
ELSE
   ALLOCATE(prop_user(10))
ENDIF
!
IF(id == 0 ) THEN                      ! User did not specify number
   IF(act=='ignore') THEN
      lsearch = .TRUE.
      lall    = .FALSE.
   ELSE
      lsearch = .FALSE.
      lall    = .FALSE.
   ENDIF
ELSEIF(id==-1) THEN                     ! Usere requested all entries
   lsearch = .FALSE.
   lall    = .TRUE.                     ! User did specify number
ELSE
   lsearch = .FALSE.
   lall    = .FALSE.
ENDIF
!                                       ! Determine values and store in temporary pp
IF(act=='absent') THEN
   pp%act = -1
ELSEIF(act=='ignore') THEN
   pp%act =  0
ELSEIF(act=='present') THEN
   pp%act =  1
ELSE
   ier_num = -6
   ier_typ = ER_COMM
   ier_msg(1) = 'Selection must be: ignore, present, absent'
   RETURN
ENDIF
IF(opara(6)=='none') THEN
   IF(ALLOCATED(at_kind)) DEALLOCATE(at_kind)
   ALLOCATE(at_kind(0:1))
   at_kind(1:1) = -2
   at_kind(0) = 1
   pp%at_type = -2
ELSE
   cpara(1) = opara(6)
   lpara(1) = lopara(6)
   werte(:) = 0
   ianz     = 1
!  store requested atom types
   CALL get_iscat(ianz, cpara, lpara, werte, MAXW, .FALSE. )
   IF(ier_num/=0) THEN
      RETURN
   ENDIF
   IF(ALLOCATED(at_kind)) DEALLOCATE(at_kind)
   ALLOCATE(at_kind(0:ianz))
   at_kind(1:ianz) = NINT(werte(1:ianz))
   at_kind(0) = ianz
!
   pp%at_type = NINT(werte(1))
ENDIF
!
!   store requested connectivities
loop_type: DO j=1, at_kind(0)
   pp%at_type = at_kind(j)
!
   pp%conn_no   = 0
   pp%conn_name = ' '
   IF(opara(7)=='none') THEN
      pp%conn_no   = 0
      pp%conn_name = ' '
   ELSE
      cpara(1) = opara(7)
      lpara(1) = lopara(7)
      ianz     = 1
      CALL ber_params (ianz, cpara, lpara, werte, MAXW)
      IF(pp%at_type==-1) THEN      ! All atom type selected
         IF(ier_num==0) THEN
            pp%conn_no   = NINT(werte(1))
            pp%conn_name = ' '
         ELSE
            RETURN
         ENDIF
      ELSE
         IF(ier_num==0) THEN
            ino      = NINT(werte(1))
            c_name   = ' '
            c_name_l = 0
         ELSE
            ino      = 0
            c_name   = opara(6)
            c_name_l = lopara(6)
         ENDIF
         CALL get_connectivity_identity( pp%at_type, ino, c_name, c_name_l)
         IF(ier_num==0) THEN
            pp%conn_no = ino
            pp%conn_name = c_name
         ELSE
            pp%conn_no = ino
            pp%conn_name = c_name
            ier_num = 0
            ier_typ = 0
            ier_msg(:) = ' '
         ENDIF
      ENDIF
   ENDIF
   pp%n_min   = NINT(owerte(1))
   pp%n_max   = NINT(owerte(2))
   pp%e_min   = NINT(owerte(3))
   pp%e_max   = NINT(owerte(4))
   IF(pp%n_min > pp%n_max) THEN
      ier_msg(1) = 'Lower boundary of included range higher than upper'
      ier_num = -155
      ier_typ = ER_APPL
      RETURN
   ENDIF
   IF(pp%e_min > pp%e_max) THEN
      ier_msg(1) = 'Lower boundary of excluded range higher than upper'
      ier_num = -155
      ier_typ = ER_APPL
      RETURN
   ENDIF
!
! Values are determined, find entry and store
!
   IF((prop_user_no>= UBOUND(prop_user,1) .OR.     &
       id >UBOUND(prop_user,1)                 )   &
      .AND. act/='ignore'                             ) THEN   ! Need more space
      nuser = MAX(nuser, prop_user_no+10, UBOUND(prop_user,1), id)
      ALLOCATE(p_temp(nuser))                                  ! Make temporary storage
      p_temp(1:prop_user_no) = prop_user(1:prop_user_no)       ! Restore old data
      CALL MOVE_ALLOC(FROM=p_temp, TO=prop_user)
   ENDIF
   IF(lsearch) THEN     !User did not specify number, and act='ignore'
     id = -2
     search: DO i=1,prop_user_no
       IF(prop_user(i)%at_type   == pp%at_type    .AND. &
          prop_user(i)%conn_no   == pp%conn_no    .AND. &
          prop_user(i)%conn_name == pp%conn_name  .AND. &
          prop_user(i)%n_min     == pp%n_min      .AND. &
          prop_user(i)%n_max     == pp%n_max      .AND. &
          prop_user(i)%e_min     == pp%e_min      .AND. &
          prop_user(i)%e_max     == pp%e_max           ) THEN
!
           id = i
           EXIT search
        ENDIF
     ENDDO search
     IF(id==-2) THEN
        RETURN
     ENDIF
   ELSE
     IF(id==0) THEN     !User did not specify number, and act/='ignore' Find empty
        id = prop_user_no + 1                  ! Default to a new entry if no empty is found
        empty: DO i=1,prop_user_no
           IF(prop_user(i)%act==0) THEN
              id = i                           ! Empty entry was found, take this as id
              EXIT empty
           ENDIF
        ENDDO empty
      ENDIF
   ENDIF
   IF(pp%act==0) THEN                         ! Turn off the property
      IF(lall  ) THEN                         ! Turn off all
         DO id=1,prop_user_no
            prop_user(id)%act       = 0             ! Found entry back to default for all parameter
            prop_user(id)%at_type   = 0
            prop_user(id)%conn_no   = 0
            prop_user(id)%conn_name = ' '
            prop_user(id)%n_min     = 0
            prop_user(id)%n_max     = 192
            prop_user(id)%e_min     = -1
            prop_user(id)%e_max     = -1
         ENDDO
         prop_user_no = 0
      ELSE
         prop_user(id)%act       = 0             ! Found entry back to default for all parameter
         prop_user(id)%at_type   = 0
         prop_user(id)%conn_no   = 0
         prop_user(id)%conn_name = ' '
         prop_user(id)%n_min     = 0
         prop_user(id)%n_max     = 192
         prop_user(id)%e_min     = -1
         prop_user(id)%e_max     = -1
      ENDIF
   ELSE
      prop_user(id)%act       = pp%act
      prop_user(id)%at_type   = pp%at_type
      prop_user(id)%conn_no   = pp%conn_no
      prop_user(id)%conn_name = pp%conn_name
      prop_user(id)%n_min     = pp%n_min
      prop_user(id)%n_max     = pp%n_max
      prop_user(id)%e_min     = pp%e_min
      prop_user(id)%e_max     = pp%e_max
      prop_user_no = MAX(prop_user_no, id)
   ENDIF
   id = 0
!
ENDDO loop_type
!
lactive = .FALSE.
active: DO i=1,prop_user_no
   IF(prop_user(i)%act/=0) THEN
      lactive = .TRUE.
      EXIT active
   ENDIF
ENDDO active
IF(.NOT.lactive) prop_user_no = 0
!
END SUBROUTINE property_set_user
!
!*******************************************************************************
!
SUBROUTINE property_show_user
!
USE crystal_mod
USE atom_name
USE prop_para_mod
USE prompt_mod
!
IMPLICIT NONE
!
CHARACTER(LEN=7), DIMENSION(-1:1) :: char_l
CHARACTER(LEN=1024)               :: line
INTEGER :: i,ll
DATA char_l /'absent ', 'ignore ', 'present' /
!
IF(prop_user_no>0) THEN
   WRITE(output_io,'(a)')  ' USERS:  No., Status ,Atoms     , Neighbors  , Excluded   , Conn, Conn name'
   DO i=1, prop_user_no
      line = ' '
      ll = 0
      IF(prop_user(i)%act/=0) THEN
      WRITE(line(ll+1:ll+13),'(i4,'', '',a7)') i, char_l(prop_user(i)%act)
      ll = 13
      IF(prop_user(i)%at_type==-1) THEN
         WRITE(line(ll+1:ll+11), '('', '',a9)') 'all'
      ELSE
         WRITE(line(ll+1:ll+11), '('', '',a9)') at_name(prop_user(i)%at_type)
      ENDIF
      ll = 24
      WRITE(line(ll+1:ll+26), '('', ['',i4,'':'',i4,''], ['',i4,'':'',i4,'']'')') &
           prop_user(i)%n_min, prop_user(i)%n_max, &
           prop_user(i)%e_min, prop_user(i)%e_max
      ll = 50
      IF(prop_user(i)%conn_no==0) THEN
         WRITE(line(ll+1: ll+6), '(a)') ', none'
         ll = 56
      ELSE
         WRITE(line(ll+1:), '('', '',i4,'', '',a)') prop_user(i)%conn_no, &
         prop_user(i)%conn_name(1:LEN_TRIM(prop_user(i)%conn_name))
         ll = LEN_TRIM(line)
      ENDIF
         WRITE(output_io,'(a,a)') ' Users: ',line(1:ll)
      ENDIF
   ENDDO
!ELSE
!   WRITE(output_io,'(a)') 'No user properties defined'
ENDIF 
!
END SUBROUTINE property_show_user
END MODULE prop_para_func
