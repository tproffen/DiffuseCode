MODULE surface_func_mod
!
CONTAINS
!*****7*****************************************************************
!
SUBROUTINE surface_menu 
!-                                                                      
!     Main menu for surface related operations                          
!+                                                                      
USE discus_config_mod 
USE crystal_mod 
USE molecule_mod 
USE prop_para_mod 
!
USE calc_expr_mod
USE doact_mod 
USE do_eval_mod
USE do_wait_mod
USE errlist_mod 
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
USE sup_mod
!                                                                       
IMPLICIT none 
!                                                                       
CHARACTER(len=9) :: befehl 
CHARACTER(LEN=LEN(prompt)) :: orig_prompt
CHARACTER(LEN=PREC_STRING) :: line, zeile
INTEGER :: lp, length, lbef 
INTEGER :: indxg
LOGICAL :: lend
!                                                                       
!                                                                       
lend = .false. 
CALL no_error 
orig_prompt = prompt
prompt = prompt (1:len_str (prompt) ) //'/surf' 
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
              .AND..NOT. (str_comp (befehl, 'system', 2, lbef, 6) )    &
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
                     line(1:length-1) = line(2:length)
                     length = 1
                     CALL file_kdo(line, length)
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
               ELSEIF (str_comp (befehl, 'evaluate', 2, lbef, 4) ) then 
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
                     lp = lp + 15 
                     CALL do_hel ('discus surface '//zeile, lp) 
                  ENDIF 
!                                                                       
!------- -Operating System Kommandos 'syst'                             
!                                                                       
               ELSEIF (str_comp (befehl, 'system', 2, lbef, 6) ) then 
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
!     ----run a boundary 'boundary'                                     
!                                                                       
               ELSEIF (str_comp (befehl, 'boundary', 2, lbef, 8) ) then 
                  CALL boundary (zeile, lp) 
!                                                                       
!     ----Test the boundary character of an atom                        
!                                                                       
               ELSEIF (str_comp (befehl, 'character', 2, lbef, 9) ) then 
                  CALL do_surface_char(zeile, lp)
!                                                                       
!     ----define various surface related settings 'set'                 
!                                                                       
               ELSEIF (str_comp (befehl, 'set', 2, lbef, 3) ) then 
                  CALL surf_do_set (zeile, lp) 
!                                                                       
!     ----show current parameters 'show'                                
!                                                                       
               ELSEIF (str_comp (befehl, 'show', 2, lbef, 4) ) then 
                  CALL surf_show 
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
                  ier_msg(1) = ' Error occured in surface menu'
                  prompt_status = PROMPT_ON 
                  prompt = orig_prompt
                  RETURN
               ELSE
                  IF(lmacro_close) THEN
                     CALL macro_close(-1)
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
END SUBROUTINE surface_menu                   
!
!*****7*****************************************************************
!
SUBROUTINE surf_do_set (zeile, length) 
!+                                                                      
!     This subroutine sets various parameters                           
!-                                                                      
USE discus_config_mod 
USE errlist_mod 
USE get_params_mod
USE precision_mod
USE str_comp_mod
!                                                                       
IMPLICIT none 
!                                                                       
INTEGER, PARAMETER :: MAXW = 20
!                                                                       
CHARACTER(len=*), intent(inout) :: zeile 
INTEGER         , intent(inout) :: length 
!                                                                       
CHARACTER(LEN=PREC_STRING), dimension(MAXW) :: cpara ! (maxw) 
INTEGER, dimension(MAXW) :: lpara ! (maxw) 
INTEGER :: ianz 
REAL(KIND=PREC_DP), DIMENSION(MAXW) :: werte
!                                                                       
!                                                                       
CALL get_params (zeile, ianz, cpara, lpara, maxw, length) 
IF (ier_num.ne.0) return 
IF (ianz.le.0) return 
!                                                                       
IF (str_comp (cpara (1) , 'distance', 2, lpara (1) , 8) ) then 
   CALL del_params (1, ianz, cpara, lpara, maxw) 
   CALL surf_set_fuzzy (ianz, cpara, lpara, werte, maxw, 0) 
ELSEIF(str_comp(cpara(1), 'surface', 2, lpara (1) , 7) ) then 
   CALL del_params (1, ianz, cpara, lpara, maxw) 
   CALL surf_set_vector(MAXW, ianz, cpara, lpara, werte)
ELSE 
   ier_num = - 6 
   ier_typ = ER_COMM 
ENDIF 
!                                                                       
END SUBROUTINE surf_do_set                    
!
!*****7*****************************************************************
!
SUBROUTINE surf_set_fuzzy (ianz, cpara, lpara, werte, maxw, iflag) 
!+                                                                      
!     This subroutine sets the distances between atoms and a surface.   
!-                                                                      
USE discus_config_mod 
USE discus_allocate_appl_mod
USE crystal_mod 
USE get_iscat_mod
USE surface_mod 
USE berechne_mod
use do_replace_expr_mod
USE errlist_mod 
USE get_params_mod
USE precision_mod
USE str_comp_mod
!                                                                       
IMPLICIT none 
!                                                                       
INTEGER                            , intent(inout) :: ianz 
INTEGER                            , intent(in)    :: MAXW
CHARACTER(LEN=*)  , dimension(MAXW), intent(inout) :: cpara  !(maxw) 
INTEGER           , dimension(MAXW), intent(inout) :: lpara ! (maxw) 
REAL(KIND=PREC_DP), dimension(MAXW), intent(inout) :: werte ! (maxw) 
INTEGER                            , intent(in)    :: iflag 
!                                                                       
CHARACTER(LEN=MAX(PREC_STRING,LEN(cpara))) :: string 
INTEGER :: laenge 
INTEGER :: i 
LOGICAL :: lnew 
LOGICAL :: linternal = .true.
LOGICAL :: lexternal = .false.
REAL(kind=PREC_DP) ::  distance 
!                                                                       
!                                                                       
lnew = .false. 
!                                                                       
IF (ianz.eq.1) then 
   ier_num = - 6 
   ier_typ = ER_COMM 
   RETURN 
ENDIF 
!                                                                       
IF( MAXSCAT > SURF_MAXSCAT) THEN
   CALL alloc_surf ( MAXSCAT )
   IF ( ier_num < 0 ) THEN
      RETURN
   ENDIF
ENDIF
!
!     WRITE ( * , * ) 'IFLAG', iflag 
!     WRITE ( * , * ) 'IANZ ', ianz 
!     DO i = 1, ianz 
!     WRITE ( * , * ) '>', cpara (i) (1:lpara (i) ) , '<' 
!     ENDDO 
      IF (iflag.eq.SURF_SURFACE) THEN 
         IF (str_comp (cpara (1) , 'external', 2, lpara (1) , 8) ) then 
            lexternal = .true. 
            linternal = .false. 
         ELSEIF (str_comp (cpara (1) , 'internal', 2, lpara (1) , 8) )  &
         then                                                           
            lexternal = .false. 
            linternal = .true. 
         ELSEIF (str_comp (cpara (1) , 'all', 2, lpara (1) , 3) ) then 
            lexternal = .true. 
            linternal = .true. 
         ELSE 
            ier_num = - 6 
            ier_typ = ER_COMM 
            RETURN 
         ENDIF 
         CALL del_params (1, ianz, cpara, lpara, maxw) 
      ELSEIF (IFLAG.eq.SURF_DOMAIN) THEN 
         lexternal = .false. 
         linternal = .true. 
      ENDIF 
!                                                                       
      IF (ianz.gt.1) then 
         IF (str_comp (cpara (1) , 'new', 2, lpara (1) , 3) ) then 
            CALL del_params (1, ianz, cpara, lpara, maxw) 
            IF (ianz.eq.1) then 
               string = '('//cpara (ianz) (1:lpara (ianz) ) //')' 
               laenge = lpara (ianz) + 2 
               call do_replace_expr(string, laenge)
               distance = berechne (string, laenge) 
               DO i = cr_nscat + 1, SURF_MAXSCAT 
               IF (linternal) then 
                  surf_in_dist (i) = distance 
               ENDIF 
               IF (lexternal) then 
                  surf_ex_dist (i) = distance 
               ENDIF 
               ENDDO 
            ELSE 
               ier_num = - 6 
               ier_typ = ER_COMM 
               RETURN 
            ENDIF 
         ELSE 
            string = '('//cpara (ianz) (1:lpara (ianz) ) //')' 
            laenge = lpara (ianz) + 2 
            call do_replace_expr(string, laenge)
            distance = berechne (string, laenge) 
            IF (ier_num.ne.0) return 
            ianz = ianz - 1 
            CALL get_iscat (ianz, cpara, lpara, werte, maxw, lnew) 
            IF (ier_num.ne.0) return 
            IF(NINT(werte(1)) == -1) then 
               DO i = 0, cr_nscat 
               IF (linternal) then 
                  surf_in_dist (i) = distance 
               ENDIF 
               IF (lexternal) then 
                  surf_ex_dist (i) = distance 
               ENDIF 
               ENDDO 
            ELSE 
               DO i = 1, ianz 
               IF (0.le.nint (werte (i) ) .and.nint (werte (i) )        &
               .le.cr_nscat) then                                       
                  IF (linternal) then 
                     surf_in_dist (nint (werte (i) ) ) = distance 
                  ENDIF 
                  IF (lexternal) then 
                     surf_ex_dist (nint (werte (i) ) ) = distance 
                  ENDIF 
               ELSE 
                  ier_num = - 27 
                  ier_typ = ER_APPL 
                  RETURN 
               ENDIF 
               ENDDO 
            ENDIF 
         ENDIF 
      ELSE 
         IF (str_comp (cpara (1) , 'default', 2, lpara (1) , 7) ) then 
            DO i = 0, SURF_MAXSCAT 
            IF (lexternal) then 
               surf_ex_dist (i) = SURF_DIST_DEF 
            ENDIF 
            IF (linternal) then 
               surf_in_dist (i) = SURF_DIST_DEF 
            ENDIF 
            ENDDO 
         ELSEIF (str_comp (cpara (1) , 'off', 2, lpara (1) , 3) ) then 
            DO i = 0, SURF_MAXSCAT 
            IF (lexternal) then 
               surf_ex_dist (i) = 0.0 
            ENDIF 
            IF (linternal) then 
               surf_in_dist (i) = 0.0 
            ENDIF 
            ENDDO 
         ELSE 
            ier_num = - 6 
            ier_typ = ER_COMM 
         ENDIF 
      ENDIF 
!                                                                       
END SUBROUTINE surf_set_fuzzy                 
!
!*****7*****************************************************************
!
SUBROUTINE surf_set_vector(MAXW, ianz, cpara, lpara, werte)
!+
!  This subroutine allows the user to set a surface character for a 
!  given atom or atom type. 
!-
USE crystal_mod 
USE surface_mod 
!
USE get_params_mod
USE precision_mod
USE str_comp_mod
!
IMPLICIT NONE
!
INTEGER                             , INTENT(IN) :: MAXW
INTEGER                             , INTENT(INOUT) :: ianz
CHARACTER(LEN=*   ), DIMENSION(MAXW), INTENT(INOUT) :: cpara
INTEGER            , DIMENSION(MAXW), INTENT(INOUT) :: lpara
REAL(KIND=PREC_DP) , DIMENSION(MAXW), INTENT(INOUT) :: werte
!
INTEGER :: istart
INTEGER :: ifinish
!
IF(str_comp(cpara(1) , 'local', 3, lpara(1), 5)) THEN
!
   istart = 1
   ifinish = cr_natoms
   CALL surf_set_local(istart, ifinish)
ENDIF
!
END SUBROUTINE surf_set_vector                 
!
!*****7*****************************************************************
!
SUBROUTINE surf_set_local(istart, ifinish)                 
!+
!   Find a local surface
!-
USE crystal_mod
USE atom_env_mod
USE chem_mod
USE do_find_mod
USE metric_mod
USE prop_para_mod
USE surface_mod
!
USE errlist_mod
USE math_sup
USE param_mod
USE precision_mod
!
IMPLICIT NONE
!
INTEGER, INTENT(IN) :: istart
INTEGER, INTENT(IN) :: ifinish
!
INTEGER, PARAMETER    :: MAXW = 1
LOGICAL, PARAMETER    :: LSPACE = .TRUE.
REAL(KIND=PREC_DP), DIMENSION(1:3), PARAMETER :: VNULL=(/0.0D0, 0.0D0, 0.0D0/)
!
CHARACTER(LEN=PREC_STRING)   :: line
INTEGER               :: i, j
INTEGER               :: laenge
INTEGER               :: idiv
INTEGER               :: ianz
INTEGER, DIMENSION(2) :: n_neigh
INTEGER, DIMENSION(:,:), ALLOCATABLE :: neigh
LOGICAL, DIMENSION(3) :: fp
LOGICAL               :: fq
LOGICAL               :: lsurf
REAL(KIND=PREC_DP) , DIMENSION(1:MAXW) :: werte
REAL(KIND=PREC_DP) , DIMENSION(1:3   ) :: u
REAL(KIND=PREC_DP) , DIMENSION(1:3   ) :: x
REAL(KIND=PREC_DP) , DIMENSION(1:3,2 ) :: cofm     ! Center of mass for neighbors
REAL(KIND=PREC_DP) , DIMENSION(1:3,2 ) :: vect     ! Vectors to Center of mass for neighbors
REAL(kind=PREC_DP) :: rmin   = 0.01
REAL(kind=PREC_DP) :: rsmall = 2.50
REAL(kind=PREC_DP) :: rbig   = 5.00
REAL(kind=PREC_DP) :: alpha
REAL(kind=PREC_DP) :: alpha_max0 =   0.0
REAL(kind=PREC_DP) :: alpha_max1 =   0.0
REAL(kind=PREC_DP) :: alpha_max2 =   0.0
REAL(kind=PREC_DP) :: alpha_min0 = 200.0
REAL(kind=PREC_DP) :: alpha_min1 = 200.0
REAL(kind=PREC_DP) :: alpha_min2 = 200.0
!
!
fp(1) = .FALSE.
fp(2) = .FALSE.
fp(3) = .FALSE.
fq    = .FALSE.
n_neigh(:) = 0
alpha_max0 =   0.0
alpha_max1 =   0.0
alpha_max2 =   0.0
alpha_min0 = 200.0
alpha_min1 = 200.0
alpha_min2 = 200.0
ALLOCATE(neigh(1:1,2))
!
main_loop: DO i=istart, ifinish
!main_loop: DO i= 63, 63
!     If atom has surface property AND only new surface sites are to be found; CYCLE
   IF(IBITS(cr_prop(i),PROP_SURFACE_EXT,1).eq.1 .AND.   &  ! real Atom is near surface
      surf_local_new                                 ) CYCLE main_loop
!
   werte(1)  = -1                        ! Find all atom types
   x(:)      = cr_pos(:,i)
   cofm(:,:) = 0.0
!
   ianz = 1                              ! Find all neighbors in large shell
   CALL do_find_env(ianz, werte, MAXW, x, rmin, rbig, fq, fp)
   IF(atom_env(0) == 0) THEN
      CYCLE main_loop                    ! Atom has no neighbors, ignore
   ENDIF
   IF(atom_env(0)>UBOUND(neigh,1)) THEN  ! Adjust size of neighbor array
      IF(ALLOCATED(neigh)) DEALLOCATE(neigh)
      ALLOCATE(neigh(1:atom_env(0),2))
   ENDIF
!                                        ! Copy neighbors for later use
   neigh(1:atom_env(0),2) = atom_env(1:atom_env(0))
   n_neigh(2) = atom_env(0)
   DO j=1, atom_env(0)                   ! Calculate Center of mass 
      cofm(:,2) = cofm(:,2) + cr_pos(:,atom_env(j))
   ENDDO
   cofm(:,2) = cofm(:,2)/REAL(atom_env(0), kind=PREC_DP)
   vect(:,2) = x(:) - cofm(:,2)          ! Vector CofM ==> Atom
!
!  Same for a smaller shell
!
   CALL do_find_env(atom_env(0), werte, maxw, x, rmin, rsmall, fq, fp)
   IF(atom_env(0) == 0) THEN
      CYCLE main_loop                    ! Atom has no neighbors, ignore
   ENDIF
   neigh(1:atom_env(0),1) = atom_env(1:atom_env(0))
   n_neigh(1) = atom_env(0)
   DO j=1, atom_env(0)
      cofm(:,1) = cofm(:,1) + cr_pos(:,atom_env(j))
   ENDDO
   cofm(:,1) = cofm(:,1)/REAL(atom_env(0), kind=PREC_DP)
   vect(:,1) = x(:) - cofm(:,1)
!
   u(:) = vect(:,1)                      ! If |vect|==0, atom is surrounded, not at surface 
   IF(do_blen(LSPACE, u, VNULL)==0.0) CYCLE main_loop
   u(:) = vect(:,2)
   IF(do_blen(LSPACE, u, VNULL)==0.0) CYCLE main_loop
!
!  if the two vectors CofM==> atom are reasonable parallel we might be at a surface
   alpha = do_bang(LSPACE, vect(1:3,1), VNULL, vect(1:3,2))
   lsurf = alpha < 60.0
   IF(lsurf) alpha_max0 = MAX(alpha, alpha_max0)
   IF(lsurf) alpha_min0 = MIN(alpha, alpha_min0)
   IF(lsurf) THEN                            ! test for indented surface
!
!     At the surface all vectors atom==>neigh should make a "large"
!     angle to the vector CofM==>atom
      DO j=1, n_neigh(1)
         u(:)  = cr_pos(:,neigh(j,1)) - x(:)
         alpha = do_bang(LSPACE, u, VNULL, vect(1:3,1))
         lsurf = lsurf .AND. alpha > 60.0
         IF(.NOT.lsurf) CYCLE main_loop
         IF(lsurf) alpha_max1 = MAX(alpha, alpha_max1)
         IF(lsurf) alpha_min1 = MIN(alpha, alpha_min1)
      ENDDO
      DO j=1, n_neigh(2)
         u(:) = cr_pos(:,neigh(j,2)) - x(:)
         alpha = do_bang(LSPACE, u, VNULL, vect(1:3,2))
          lsurf = lsurf .AND. alpha > 60.0
         IF(.NOT.lsurf) CYCLE main_loop
         IF(lsurf) alpha_max2 = MAX(alpha, alpha_max2)
         IF(lsurf) alpha_min2 = MIN(alpha, alpha_min2)
      ENDDO
      IF(lsurf) THEN
!
!        Atom is at surface get a rough surface normal as 
!        CofM(large) ==> atom
!
!        May need to be refined by analysis of all surface neighbors only
!        This would require a separate loop
         cr_surf(0,i) = SURF_LOCAL
         WRITE(line,'(2(G16.8E3,a1),G16.8E3)') vect(1,2), ',', vect(2,2), ',', vect(3,2)
         laenge = 50
         CALL d2r(line, laenge, LSPACE)
         u(1:3) = INT(res_para(1:3))      ! Rough normal 
         CALL surface_normalize(u)
         cr_surf(1:3, i) = NINT(u(:))
!
         idiv = gcd(cr_surf(1, i), cr_surf(2, i), cr_surf(3, i))
         IF(idiv /= 0) THEN
             cr_surf(1:3,i) = cr_surf(1:3, i)/IABS(idiv)
         ENDIF
         cr_prop(i) = IBSET(cr_prop(i), PROP_SURFACE_EXT)
      ENDIF
   ENDIF
!
ENDDO main_loop
!
IF(ALLOCATED(neigh)) DEALLOCATE(neigh)
!
END SUBROUTINE surf_set_local                 
!
!*****7*****************************************************************
!
SUBROUTINE surface_normalize(u)
!
USE crystal_mod
USE metric_mod
!
USE precision_mod
!
IMPLICIT NONE
!
REAL(KIND=PREC_DP), DIMENSION(3), INTENT(INOUT) :: u
!
LOGICAL, PARAMETER    :: LSPACE = .FALSE.
REAL(KIND=PREC_DP), DIMENSION(1:3), PARAMETER :: VNULL=(/0.0D0, 0.0D0, 0.0D0/)
REAL(kind=PREC_DP) :: dstar
REAL(kind=PREC_DP) :: umax
!
dstar = do_blen (LSPACE, u, VNULL) 
IF(dstar > 0.0) THEN
   u(:) = u(:) / dstar
   umax=MAX(ABS(u(1)), ABS(u(2)), ABS(u(3)))
   IF(umax > 9.900) THEN
      u(:) = u(:) * 9.90/umax
   ENDIF
   u(:) = 10*u(:)
ENDIF
!
END SUBROUTINE surface_normalize
!
!*****7*****************************************************************
!
SUBROUTINE surf_show 
!+                                                                      
!     This subroutine shows the surface settings                        
!-                                                                      
USE discus_config_mod 
USE crystal_mod 
USE surface_mod 
USE discus_allocate_appl_mod
USE prompt_mod 
!
USE errlist_mod 
!
IMPLICIT none 
!                                                                       
INTEGER :: i 
!                                                                       
IF( MAXSCAT > SURF_MAXSCAT) THEN
   CALL alloc_surf ( MAXSCAT )
   IF ( ier_num < 0 ) THEN
      RETURN
   ENDIF
ENDIF
!                                                                       
IF (SURF_MAXSCAT==0) THEN
   WRITE(output_io,*) ' No distances to surfaces have been defined yet'
   WRITE(output_io,*) ' Set distances first'
   RETURN
ENDIF
WRITE (output_io, * ) 
WRITE (output_io, 3000) 
WRITE (output_io, * ) 
DO i = 0, cr_nscat 
   WRITE (output_io, 3010) cr_at_lis (i), surf_ex_dist (i) 
ENDDO 
IF (cr_nscat.lt.SURF_MAXSCAT) then 
   WRITE (output_io, 3010) 'new', surf_ex_dist (SURF_MAXSCAT) 
ENDIF 
WRITE (output_io, * ) 
WRITE (output_io, 3100) 
WRITE (output_io, * ) 
DO i = 0, cr_nscat 
   WRITE (output_io, 3010) cr_at_lis (i), surf_in_dist (i) 
ENDDO 
IF (cr_nscat.lt.SURF_MAXSCAT) then 
   WRITE (output_io, 3010) 'new', surf_in_dist (SURF_MAXSCAT) 
ENDIF 
!                                                                       
 3000 FORMAT (' Distances between atom types and an external surface') 
 3100 FORMAT (' Distances between atom types and an internal surface') 
 3010 FORMAT (' Atom type   ',a4,' at ',f8.3,' Angstroem') 
!
END SUBROUTINE surf_show                      
!
!*****7*****************************************************************
!
SUBROUTINE boundary (zeile, lp) 
!+                                                                      
!     This subroutine removes all atomes outside a given boundary.      
!-                                                                      
USE metric_mod
USE discus_config_mod 
USE crystal_mod 
use chem_aver_mod
USE discus_plot_mod
USE discus_plot_init_mod
USE point_grp
USE prop_para_mod 
USE symm_sup_mod
USE surface_mod
USE wyckoff_mod 
USE ber_params_mod
USE errlist_mod 
USE get_params_mod
USE lib_errlist_func
USE precision_mod
USE take_param_mod
USE str_comp_mod
!
IMPLICIT none 
!                                                                       
       
!                                                                       
INTEGER, PARAMETER :: maxw = 12
INTEGER, PARAMETER :: maxww = 3
INTEGER, PARAMETER :: NOPTIONAL = 11
!                                                                       
CHARACTER (LEN=* ), INTENT(INOUT) :: zeile 
INTEGER           , INTENT(INOUT) :: lp 
!                                                                       
real(kind=PREC_DP), parameter :: TOL = 1.0d-9
CHARACTER(LEN=MAX(PREC_STRING,LEN(ZEILE))) :: cpara (maxw) 
CHARACTER(LEN=MAX(PREC_STRING,LEN(ZEILE))) :: ccpara (maxw) 
CHARACTER(LEN=   6), DIMENSION(NOPTIONAL) :: oname   !Optional parameter names
CHARACTER(LEN=MAX(PREC_STRING,LEN(ZEILE))), DIMENSION(NOPTIONAL) :: opara   !Optional parameter strings returned
INTEGER            , DIMENSION(NOPTIONAL) :: loname  !Lenght opt. para name
INTEGER            , DIMENSION(NOPTIONAL) :: lopara  !Lenght opt. para name returned
LOGICAL            , DIMENSION(NOPTIONAL) :: lpresent ! Optional parameter is present
REAL(KIND=PREC_DP) , DIMENSION(NOPTIONAL) :: owerte   ! Calculated values
INTEGER, PARAMETER                        :: ncalc = 1 ! Number of values to calculate 
INTEGER, PARAMETER :: O_THICK   = 1
INTEGER, PARAMETER :: O_CENTX   = 2
INTEGER, PARAMETER :: O_CENTY   = 3
INTEGER, PARAMETER :: O_CENTZ   = 4
INTEGER, PARAMETER :: O_CENTER  = 5
INTEGER, PARAMETER :: O_KEEP    = 6
INTEGER, PARAMETER :: O_ACCUM   = 7
INTEGER, PARAMETER :: O_EXEC    = 8
INTEGER, PARAMETER :: O_THIRD   = 9
INTEGER, PARAMETER :: O_FIRST   = 10
INTEGER, PARAMETER :: O_AXIS    = 11
INTEGER lpara (maxw) 
integer, dimension(1:MAXWW) :: llpara
INTEGER :: i, j, k, ianz 
integer :: iianz            ! Number of parameters for cent 
INTEGER :: special_form
INTEGER :: special_n = 0
INTEGER :: nplanes       ! number of planes that atom is close to
INTEGER :: iplane        !  index of plane  that atom is closest to
LOGICAL lspace 
LOGICAL linside 
LOGICAL l_plane 
LOGICAL l_sphere 
LOGICAL l_form 
LOGICAL l_cyl 
LOGICAL l_ell 
LOGICAL l_local
LOGICAL l_special 
LOGICAL lwall           ! True if atom is close to cylinder wall
LOGICAL ltop            ! True if atom is close to cylinder top
LOGICAL :: lrem         ! True if cylinder need to remove an atom
!      LOGICAL l_new 
REAL(kind=PREC_DP), DIMENSION(3) :: wall  !local normal at cylinder wall
REAL(kind=PREC_DP), DIMENSION(3) :: top   !local normal at cylinder top
REAL(kind=PREC_DP) :: h (3), d, dstar, radius, height , dshort
REAL(kind=PREC_DP) :: hkl(4)!, hklw(4)
REAL(kind=PREC_DP), DIMENSION(3)    :: center       ! center of the shape
REAL(kind=PREC_DP), DIMENSION(3)    :: com   ! Center of Mass
REAL(kind=PREC_DP), DIMENSION(3, 2) :: special_hkl
REAL(kind=PREC_DP), DIMENSION(:,:), ALLOCATABLE :: point_hkl ! (3,48)
REAL(kind=PREC_DP), DIMENSION(:,:), ALLOCATABLE, SAVE :: accum_hkl ! (3,48)
REAL(kind=PREC_DP), DIMENSION(:,:), ALLOCATABLE ::  temp_hkl ! (3,48)
REAL(kind=PREC_DP), DIMENSION(:  ), ALLOCATABLE ::  dstars   ! d-stars of accumulated surfaces
REAL(KIND=PREC_DP), DIMENSION(4,4) :: m_comb    ! Rotation matrix for Cyllinder / Ellipsoid
REAL(KIND=PREC_DP), DIMENSION(4,4) :: m_combr   ! Rotation matrix for Cyllinder / Ellipsoid
REAL(KIND=PREC_DP), DIMENSION(4)   :: v4        ! Augmented atom vector
INTEGER               :: point_n
INTEGER, SAVE         :: accum_n
REAL(kind=PREC_DP), DIMENSION(3)    :: radius_ell
REAL(KIND=PREC_DP), dimension(3) :: v !(3) 
REAL(kind=PREC_DP), dimension(3) :: nullv !(3) 
REAL(kind=PREC_DP) :: thick
REAL(KIND=PREC_DP) :: werte (maxw) 
REAL(KIND=PREC_DP) :: wwerte (maxw) 
!                                                                       
!                                                                       
DATA nullv / 0.0D0, 0.0D0, 0.0D0 / 
DATA oname  / 'thick', 'centx', 'centy', 'centz',  'center', 'keep ',&
              'accum', 'exec ', 'third', 'first',  'axis' /
DATA loname /  5,       5,       5,       5     ,   6      ,  4     ,&
               5     ,  4     ,  5     ,  5     ,   4 /
!
opara  =  (/ '-2.550 ', '0.0000 ', '0.0000 ', '0.0000 ','0.0000 ','inside ',    &
             'init   ', 'run    ', '[0,0,1]', '[1,0,0]','[0,0,1]'  /)   ! Always provide fresh default values
lopara =  (/  6,         6,         6       ,  6       ,  6       ,  6,         &
              4       ,  3       ,  7       ,  7       ,  7        /)
owerte =  (/ -2.55,      0.0,       0.0     ,  0.00    ,  0.0     ,  0.0,        &
              0.0     ,  0.0     ,  0.0     ,  0.0     ,  0.0  /)
!
CALL symm_store                      ! Store current symmetry settings
!                                                                       
IF (cr_v.le.0.0) then 
   ier_num = - 35 
   ier_typ = ER_APPL 
   ier_msg (1) = 'A proper unit cell must be defined' 
   ier_msg (2) = 'for this command to operate       ' 
   RETURN 
ENDIF 
!                                                                       
l_plane  = .false. 
l_form   = .false. 
l_sphere = .false. 
l_cyl    = .false. 
l_ell    = .false. 
l_local  = .FALSE.
linside  = .true. 
radius   = 0.0
height   = 0.0
dstar    = 1.0
center(:) = 0.0     ! Default center to 0.0, 0.0, 0.0
CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
IF(ier_num /= 0) RETURN
!opara  =  (/ '0.0000', '0.0000', '0.0000', '-2.550','inside', 'init  ', 'run   ' /)   ! Always provide fresh default values
!lopara =  (/  6,        6,        6      ,  6      ,  6      ,  4      ,  3 /)
!owerte =  (/  0.0,      0.0,      0.0    , -2.55   ,  0.0    ,  0.0    ,  0.0 /)
CALL get_optional(ianz, MAXW, cpara, lpara, NOPTIONAL,  ncalc, &
                  oname, loname, opara, lopara, lpresent, owerte)
IF(ier_num /= 0) RETURN
!     Handle optional and global parameters
!     IF (ier_num.eq.0) then 
if(opara(O_CENTX)=='com' .or. opara(O_CENTY)=='com' .or. opara(O_CENTZ)=='com' .or. &
   opara(O_CENTER)=='com'                                                          ) then !'com' given
   call chem_com (com,.false.)
endif
   if(lpresent(O_CENTX)) then       ! User provided 'centx:'
      if(opara(O_CENTX)=='com') write(opara(O_CENTX),'(G20.12E3)') com(1)
   else                              ! User did not provide 'centx:'
      if(lpresent(O_CENTER)) then   ! User provided 'center:'
         if(opara(O_CENTER)=='com') then
            write(opara(O_CENTX),'(G20.12E3)') com(1)
         else
            opara(O_CENTX) = opara(O_CENTER)
         endif
      endif
   endif
!  Handle centy
   if(lpresent(O_CENTY)) then       ! User provided 'centy:'
      if(opara(O_CENTY)=='com') write(opara(O_CENTY),'(G20.12E3)') com(1)
   else                              ! User did not provide 'centy:'
      if(lpresent(O_CENTER)) then   ! User provided 'center:'
         if(opara(O_CENTER)=='com') then
            write(opara(O_CENTY),'(G20.12E3)') com(1)
         else
            opara(O_CENTY) = opara(O_CENTER)
         endif
      endif
   endif
!  Handle centz
   if(lpresent(O_CENTZ)) then       ! User provided 'centz:'
      if(opara(O_CENTZ)=='com') write(opara(O_CENTZ),'(G20.12E3)') com(1)
   else                              ! User did not provide 'centz:'
      if(lpresent(O_CENTER)) then   ! User provided 'center:'
         if(opara(O_CENTER)=='com') then
            write(opara(O_CENTZ),'(G20.12E3)') com(1)
         else
            opara(O_CENTZ) = opara(O_CENTER)
         endif
      endif
   endif
!  if(opara(O_CENTX)=='com' .or.(.not.lpresent(O_CENTX) .and. opara(O_CENTER)=='com')) write(opara(O_CENTX),'(G20.12E3)') com(1)
!  if(opara(O_CENTY)=='com' .or.(.not.lpresent(O_CENTY) .and. opara(O_CENTER)=='com')) write(opara(O_CENTY),'(G20.12E3)') com(2)
!  if(opara(O_CENTZ)=='com' .or.(.not.lpresent(O_CENTZ) .and. opara(O_CENTER)=='com')) write(opara(O_CENTZ),'(G20.12E3)') com(3)
!endif
ccpara(1) = opara(O_CENTX)
ccpara(2) = opara(O_CENTY)
ccpara(3) = opara(O_CENTZ)
llpara(1) = len_trim(ccpara(1))
llpara(2) = len_trim(ccpara(2))
llpara(3) = len_trim(ccpara(3))
iianz = 3
CALL ber_params (iianz, ccpara, llpara, wwerte, maxww) 
center(1:3) = wwerte(1:3)               ! As defaults are always provided, we can take it take blindly here
!
IF (str_comp (cpara (ianz) , 'outside', 3, lpara (ianz) , 7)) THEN
   linside = .false.
   ianz = ianz - 1
ELSEIF (str_comp (cpara (ianz) , 'inside', 3, lpara (ianz), 6) ) THEN
   linside = .true.
   ianz = ianz - 1
ELSE     ! Test optional parameter form
   IF (str_comp (opara(O_KEEP) , 'outside', 3, lopara(O_KEEP) , 7)) THEN
      linside = .false.
   ELSEIF (str_comp (opara(O_KEEP) , 'inside', 3, lopara(O_KEEP), 6) ) THEN
      linside = .true.
   ENDIF
ENDIF
IF (str_comp (cpara (1) , 'hkl',  3, lpara (1) , 3)  .OR.  &
    str_comp (cpara (1) , 'form', 3, lpara (1) , 4) ) THEN 
  l_form = str_comp (cpara (1) , 'form', 3, lpara (1) , 4)
!                                                                       
!     ----Crystal is limited by a plane or by a form of symmetrically equivalent planes
!                                                                       
   lspace = .false. 
   l_special = .FALSE.
   IF(str_comp(cpara(2),'cubeoct',7,lpara(2),7)) THEN
      IF(cr_syst==CR_CUBIC) THEN
         cpara(5) = cpara(3)
         lpara(5) = lpara(3)
         cpara(2) = '1.0'
         cpara(3) = '0.0'
         cpara(4) = '0.0'
         lpara(2) = 3
         lpara(3) = 3
         lpara(4) = 3
         ianz = 5
!                 l_special = .TRUE.
         special_form = 1
      ELSE
         ier_num = -142
         ier_typ = ER_APPL
         ier_msg(1) = 'Unpredictable forms would result in non-cubic'
         ier_msg(2) = 'systems. Use explicit forms and combinations'
         ier_msg(3) = 'in non-cubic systems.'
         RETURN
      ENDIF
   ELSE
      special_form = 0
   ENDIF
   IF (ianz.ge.4) then 
      CALL del_params (1, ianz, cpara, lpara, maxw) 
      IF (ier_num.ne.0) return 
      CALL ber_params (ianz, cpara, lpara, werte, maxw) 
      IF (ier_num.ne.0) return 
      DO i = 1, 3 
      h (i) = werte (i) 
      ENDDO 
      dstar = do_blen (lspace, h, nullv) 
      IF (dstar.le.0.0) then 
         ier_num = - 32 
         ier_typ = ER_APPL 
         RETURN 
      ENDIF 
      IF (ianz.eq.4) then 
         IF(abs(werte(4)) < TOL) THEN
            DO i = 1, 3 
               h (i) = h (i) / dstar *1.0E12
            ENDDO 
            dstar = do_blen (lspace, h, nullv) 
         ELSE
            DO i = 1, 3 
               h (i) = h (i) / dstar / werte (4) 
            ENDDO 
            dstar = do_blen (lspace, h, nullv) 
         ENDIF 
      ENDIF 
      IF (dstar.le.0.0) then 
         ier_num = - 32 
         ier_typ = ER_APPL 
         RETURN 
      ENDIF 
   ELSE 
      ier_num = - 6 
      ier_typ = ER_FORT 
      RETURN 
   ENDIF 
   l_plane  = .true. 
   l_sphere = .false. 
   l_cyl    = .false. 
   l_ell    = .false. 
   IF(l_form) THEN
      l_plane= .FALSE.
      l_form = .TRUE.
   ENDIF
   IF(special_form == 0) THEN
      special_hkl(:,1) = h(:)
      special_n = 1
   ELSEIF(special_form == 1) THEN            ! CUBEOCTAHEDRON h= 1,0,0
      special_hkl(:,1) = h(:)
      special_hkl(1,2) = h(1) / 2.!/sqrt(3.)
      special_hkl(2,2) = h(1) / 2.!/sqrt(3.)
      special_hkl(3,2) = h(1) / 2.!/sqrt(3.)
      special_n = 2
   ENDIF
!
!           Handle accumulation of flat surfaces
!
   IF(str_comp(opara(O_ACCUM), 'init', 4, lopara(O_ACCUM),4)) THEN
      IF(ALLOCATED(accum_hkl)) DEALLOCATE(accum_hkl)
      accum_n = 0
   ENDIF
   hkl(1:3) = special_hkl(:,1)
   hkl(4)   = 0
   CALL point_init(hkl, point_hkl, point_n)
   IF(l_form) THEN
      DO k = 1,special_n
         hkl(1:3) = special_hkl(:,k)
         hkl(4)   = 0
         CALL point_set(hkl, point_hkl, point_n)
      ENDDO
   ENDIF
   IF(.NOT.ALLOCATED(accum_hkl)) THEN
      ALLOCATE  (accum_hkl(1:4, 1:accum_n + 48))
      accum_hkl(:,:) = 0
      accum_n = 0
   ELSE
   IF(accum_n + point_n > UBOUND(accum_hkl,1)) THEN
      IF(ALLOCATED(accum_hkl)) THEN
         k = UBOUND(accum_hkl,2)
         ALLOCATE(temp_hkl(1:4,1:k))
         temp_hkl(:,:) = accum_hkl(:,:)
         DEALLOCATE(accum_hkl)
         ALLOCATE  (accum_hkl(1:4, 1:accum_n + 48))
         accum_hkl(1:4,1:UBOUND(temp_hkl,2)) = temp_hkl(:,:)
         DEALLOCATE(temp_hkl)
      ELSE
         ALLOCATE  (accum_hkl(1:4, 1:accum_n + 48))
         accum_hkl(:,:) = 0
      ENDIF
   ENDIF
   ENDIF
   accum_hkl(1:3,accum_n+1:accum_n+point_n) = point_hkl(1:3,1:point_n)
   accum_hkl(4  ,accum_n+1:accum_n+point_n) = owerte(O_THICK)
   accum_n = accum_n + point_n
ELSEIF (str_comp (cpara (1) , 'sphere', 3, lpara (1) , 6) ) THEN
!                                                                       
!     ----Crystal is limited by a sphere                                
!                                                                       
   IF (ianz.eq.2) then 
      CALL del_params (1, ianz, cpara, lpara, maxw) 
      IF (ier_num.ne.0) return 
      CALL ber_params (ianz, cpara, lpara, werte, maxw) 
      IF (ier_num.ne.0) return 
      radius = werte (1) 
      IF (radius.le.0.0) then 
         ier_num = - 28 
         ier_typ = ER_APPL 
         RETURN 
      ENDIF 
   ELSE 
      ier_num = - 6 
      ier_typ = ER_FORT 
      RETURN 
   ENDIF 
   l_plane  = .false. 
   l_form   = .false. 
   l_sphere = .true. 
   l_cyl    = .false. 
   l_ell    = .false. 
   l_local  = .FALSE.
ELSEIF(str_comp(cpara(1), 'cylinder', 3, lpara(1), 8)) THEN
!                                                                       
!     ----Crystal is limited by a cylinder                              
!                                                                       
   IF (ianz.eq.3) then 
      CALL del_params (1, ianz, cpara, lpara, maxw) 
      IF (ier_num.ne.0) return 
      CALL ber_params (ianz, cpara, lpara, werte, maxw) 
      IF (ier_num.ne.0) return 
      radius = werte (1) 
      height = werte (2) 
      IF (radius.le.0.0.or.height.le.0.0) then 
         ier_num = - 28 
         ier_typ = ER_APPL 
         RETURN 
      ENDIF 
   ELSE 
      ier_num = - 6 
      ier_typ = ER_FORT 
      RETURN 
   ENDIF 
   l_plane  = .false. 
   l_form   = .false. 
   l_sphere = .false. 
   l_cyl    = .true. 
   l_ell    = .false. 
   l_local  = .FALSE.
ELSEIF(str_comp(cpara(1), 'ellipsoid', 3, lpara(1), 9)) THEN
!                                                                       
!     ----Crystal is limited by a standardized ellipsoid
!                                                                       
   IF (ianz.eq.4) then 
      CALL del_params (1, ianz, cpara, lpara, maxw) 
      IF (ier_num.ne.0) return 
      CALL ber_params (ianz, cpara, lpara, werte, maxw) 
      IF (ier_num.ne.0) return 
      radius_ell(1:3) = werte(1:3)*0.5  ! User provides diameters
   ELSE 
      ier_num = - 6 
      ier_typ = ER_FORT 
      RETURN 
   ENDIF 
   l_plane  = .false. 
   l_form   = .false. 
   l_sphere = .false. 
   l_cyl    = .false. 
   l_ell    = .true. 
   l_local  = .FALSE.
ELSEIF(str_comp(cpara(1), 'local', 3, lpara(1), 5)) THEN
   l_plane  = .false. 
   l_form   = .false. 
   l_sphere = .false. 
   l_cyl    = .false. 
   l_ell    = .FALSE. 
   l_local  = .TRUE.
ELSE 
   ier_num = - 6 
   ier_typ = ER_COMM 
   RETURN
ENDIF 
!
IF(l_cyl .OR. l_ell) THEN
   CALL prep_rotation(lpresent(O_THIRD), lpresent(O_FIRST), lpresent(O_AXIS),   &
        opara(O_THIRD), opara(O_FIRST), opara(O_AXIS),                          &
        lopara(O_THIRD), lopara(O_FIRST), lopara(O_AXIS), center,               &
        m_comb, m_combr)
   IF(ier_num/=0) RETURN
ENDIF
!
IF (ier_num /= 0) RETURN
IF ((l_plane .OR. l_form) .AND. str_comp (opara(O_EXEC) , 'run', 3, lopara(O_EXEC) , 3)) THEN
   ALLOCATE(dstars(1:accum_n))
   DO j=1,accum_n
      dstars(j) = do_blen (lspace, accum_hkl(1:3,j), nullv)
   ENDDO
form_loop:     DO i = 1, cr_natoms 
      IF(BTEST(cr_prop(i), PROP_OUTSIDE)) cycle form_loop 
      IF(linside) THEN
         DO j=1,accum_n
!                        dstar = do_blen (lspace, accum_hkl(1:3,j), null)
             d = 1.0 - (cr_pos (1, i)-center(1)) * accum_hkl (1,j) &
                     - (cr_pos (2, i)-center(2)) * accum_hkl (2,j) &
                     - (cr_pos (3, i)-center(3)) * accum_hkl (3,j) 
             d = d / dstars(j) 
             h(1:3) = accum_hkl (1:3,j)
             CALL boundarize_atom (center, d, i, linside, SURF_PLANE, h, accum_hkl(4,j)) 
         ENDDO 
      ELSE
         dshort  = 1.E8
         thick   = -2.55
         nplanes = 0
         iplane  = 0 
         DO j=1,accum_n
!                        dstar = do_blen (lspace, accum_hkl(1:3,j), null)
             d = 1.0 - (cr_pos (1, i)-center(1)) * accum_hkl (1,j) &
                     - (cr_pos (2, i)-center(2)) * accum_hkl (2,j) &
                     - (cr_pos (3, i)-center(3)) * accum_hkl (3,j)
             d = d / dstars(j) 
             IF(accum_hkl(4,j) > 0.0) THEN
                IF(ABS(d)<ABS(accum_hkl(4,j)) ) THEN
                   nplanes = nplanes + 1
                   iplane  = j
                   thick   = MAX(thick, ABS(accum_hkl(4,j)))
                ENDIF
             ELSE
                IF(ABS(d)<surf_ex_dist(cr_iscat(1,i) ) ) THEN 
                   nplanes = nplanes + 1
                   iplane  = j
                   thick   = MAX(thick,surf_ex_dist(cr_iscat(1,i)))
                ENDIF
             ENDIF
             dshort = MIN(dshort, d)
         ENDDO 
         IF(nplanes>2) THEN            ! Atom is at a corner
            h(:) = NINT(100*cr_pos(:,i))
            CALL boundarize_atom (center, dshort, i, linside, SURF_CORNER, h, thick) 
!           cr_surf(0,i) = SURF_CORNER
         ELSEIF(nplanes==2) THEN            ! Atom is at an edge 
            h(:) = NINT(100*cr_pos(:,i))
            CALL boundarize_atom (center, dshort, i, linside, SURF_EDGE  , h, thick) 
!           cr_surf(0,i) = SURF_EDGE
         ELSEIF(nplanes==1) THEN            ! Atom is at a PLANE 
            h(1:3) = accum_hkl (1:3,iplane)        !WRONG NEEDS WORK
            CALL boundarize_atom (center, dshort, i, linside, SURF_PLANE, h, thick) 
!           cr_surf(0,i) = SURF_PLANE
         ELSE
            h(:) = NINT(100*cr_pos(:,i))
            CALL boundarize_atom (center, dshort, i, linside, SURF_NONE , h, thick) 
         ENDIF
      ENDIF
   ENDDO form_loop
   DEALLOCATE  (accum_hkl)      ! reset the accumulation list
   accum_n = 0                  ! reset the accumulation list
   DEALLOCATE(dstars)
ELSEIF (l_sphere) then 
sphere_loop:   DO i = 1, cr_natoms 
      IF(BTEST(cr_prop(i), PROP_OUTSIDE)) cycle sphere_loop 
      v (1) = cr_pos (1, i)-center(1)
      v (2) = cr_pos (2, i)-center(2)
      v (3) = cr_pos (3, i)-center(3)
      d = radius - sqrt (v(1) * v(1) * cr_gten(1, 1)        &
         +     v(2) * v(2) * cr_gten(2, 2)    &
         +     v(3) * v(3) * cr_gten(3, 3)    &
         + 2 * v(1) * v(2) * cr_gten(1, 2)    &
         + 2 * v(1) * v(3) * cr_gten(1, 3)    &
         + 2 * v(2) * v(3) * cr_gten(2, 3)    )                                       
      h(:) = v(:)
      CALL boundarize_atom (center, d, i, linside, SURF_SPHERE, h, -1.0D0)
   ENDDO sphere_loop
ELSEIF (l_cyl) THEN
cyl_loop:      DO i = 1, cr_natoms 
      IF(BTEST(cr_prop(i), PROP_OUTSIDE)) cycle cyl_loop 
      v4(1) = cr_pos(1, i)
      v4(2) = cr_pos(2, i)
      v4(3) = cr_pos(3, i)
      v4(4) = 1.0D0
      v4 = MATMUL(m_comb,v4)
      IF(linside) THEN
         v (1) = v4(1)        ! cr_pos (1, i)-center(1)
         v (2) = v4(2)        ! cr_pos (2, i)-center(2)
         v (3) = 0.0          !-center(3)
         d = radius - sqrt(v(1) * v(1) * cr_gten(1, 1)        &
            +     v(2) * v(2) * cr_gten(2, 2)    &
            +     v(3) * v(3) * cr_gten(3, 3)    &
            + 2 * v(1) * v(2) * cr_gten(1, 2)    &
            + 2 * v(1) * v(3) * cr_gten(1, 3)    &
            + 2 * v(2) * v(3) * cr_gten(2, 3)    )                                       
         h(:) = v(:)
         CALL boundarize_atom (center, d, i, linside, SURF_CYLINDER, h, -1.0D0) 
         IF(BTEST(cr_prop(i), PROP_OUTSIDE)) cycle cyl_loop 
         v (1) = 0.0          !-center(1)
         v (2) = 0.0          !-center(2)
         v (3) = v4(3)        ! cr_pos (3, i)-center(3)
         d = height - sqrt(v(1) * v(1) * cr_gten(1, 1)        &
            +     v(2) * v(2) * cr_gten(2, 2)    & 
            +     v(3) * v(3) * cr_gten(3, 3)    & 
            + 2 * v(1) * v(2) * cr_gten(1, 2)    &
            + 2 * v(1) * v(3) * cr_gten(1, 3)    & 
            + 2 * v(2) * v(3) * cr_gten(2, 3)    )                                       
         h(:) = v(:)
         CALL boundarize_atom (center, d, i, linside, SURF_PLANE, h, -1.0D0) 
      ELSE                             ! 
         dshort = 1.E8
         lwall =.FALSE.
         ltop  =.FALSE.
         v (1) = v4(1)         !cr_pos (1, i)-center(1)
         v (2) = v4(2)         !cr_pos (2, i)-center(2)
         v (3) = 0.0          ! -center(3)
         d = radius - sqrt(v(1) * v(1) * cr_gten(1, 1)        &
            +     v(2) * v(2) * cr_gten(2, 2)    &
            +     v(3) * v(3) * cr_gten(3, 3)    &
            + 2 * v(1) * v(2) * cr_gten(1, 2)    &
            + 2 * v(1) * v(3) * cr_gten(1, 3)    &
            + 2 * v(2) * v(3) * cr_gten(2, 3)    )
         lrem = d > 0.0        ! Cylinder wall indicates removal
         wall(:) = v(:)
         IF(ABS(d)<surf_ex_dist(cr_iscat(1,i) ) ) THEN 
            lwall = .TRUE.               ! Atom is close to wall
         ENDIF
         dshort = MIN(dshort, ABS(d))
         IF(BTEST(cr_prop(i), PROP_OUTSIDE)) cycle cyl_loop 
         v (1) = 0.0         ! -center(1)
         v (2) = 0.0         ! -center(2)
         v (3) = v4(3)       !cr_pos (3, i)-center(3)
         d = height - sqrt(v(1) * v(1) * cr_gten(1, 1)        &
            +     v(2) * v(2) * cr_gten(2, 2)    & 
            +     v(3) * v(3) * cr_gten(3, 3)    & 
            + 2 * v(1) * v(2) * cr_gten(1, 2)    &
            + 2 * v(1) * v(3) * cr_gten(1, 3)    & 
            + 2 * v(2) * v(3) * cr_gten(2, 3)    )
         lrem = lrem .AND. d > 0.0   ! Cylinder wall  and top indicates removal
         top(:) = v(:)
         IF(ABS(d)<surf_ex_dist(cr_iscat(1,i) ) ) THEN 
            ltop  = .TRUE.               ! Atom is close to wall
         ENDIF
         dshort = MIN(dshort, ABS(d))
         IF(.NOT.lrem) dshort = -ABS(dshort)
         IF(lwall .AND. .NOT.ltop) THEN   ! Atom is at wall only
            CALL boundarize_atom (center, dshort, i, linside, SURF_CYLINDER, wall, -1.0D0) 
         ELSEIF(.NOT.lwall .AND. ltop) THEN   ! Atom is at top  only
            CALL boundarize_atom (center, dshort, i, linside, SURF_PLANE, top, -1.0D0) 
         ELSEIF(lwall .AND. ltop) THEN   ! Atom is at edge
            h(:) = v4(1:3) !cr_pos(:,i)
            CALL boundarize_atom (center, dshort, i, linside, SURF_EDGE , h, -1.0D0) 
         ELSE
            h(:) = v4(1:3) !cr_pos(:,i)
            CALL boundarize_atom (center, dshort, i, linside, SURF_NONE , h, -1.0D0) 
         ENDIF
      ENDIF
      v4 = MATMUL(m_combr,v4)
      cr_pos(1:3,i) = v4(1:3)
      v4(1:3) = REAL(cr_surf(1:3, i), KIND=PREC_DP)
      v4(4)   = 0.0D0
      v4 = MATMUL(m_combr,v4)
      cr_surf(1:3,i) = NINT(v4(1:3))
   ENDDO cyl_loop
ELSEIF (l_ell) then 
   CALL plot_ini_trans (1.0D0,                        &
     pl_tran_g, pl_tran_gi, pl_tran_f, pl_tran_fi, &
     cr_gten, cr_rten, cr_eps)
   radius = (radius_ell(1)*radius_ell(2)*radius_ell(3))**(1./3.)
ell_loop:      DO i = 1, cr_natoms 
      IF(BTEST(cr_prop(i), PROP_OUTSIDE)) cycle ell_loop 
      v4(1) = cr_pos(1, i)
      v4(2) = cr_pos(2, i)
      v4(3) = cr_pos(3, i)
      v4(4) = 1.0D0
      v4 = MATMUL(m_comb,v4)
      v(:) = v4(1:3)  !  cr_pos(:, i)-center(:)
      v = MATMUL(pl_tran_f(1:3,1:3), v)
      d = (1. - sqrt((v(1)/radius_ell(1))**2   &
                    +(v(2)/radius_ell(2))**2   &
                    +(v(3)/radius_ell(3))**2   ))* radius
      h(:) = v(:)
      CALL boundarize_atom (center, d, i, linside, SURF_SPHERE, h, -1.0D0)
      v4 = MATMUL(m_combr,v4)
      cr_pos(1:3,i) = v4(1:3)
      v4(1:3) = REAL(cr_surf(1:3, i), KIND=PREC_DP)
      v4(4)   = 0.0D0
      v4 = MATMUL(m_combr,v4)
      cr_surf(1:3,i) = NINT(v4(1:3))
   ENDDO ell_loop
ELSEIF(l_local) THEN 
   CALL surf_set_local(1, cr_natoms)
ENDIF 
!        ENDIF 
!     ENDIF 
IF(ier_num==0) THEN
   IF(opara(O_EXEC)=='hold' ) THEN
      IF(lpresent(O_CENTX) .OR. lpresent(O_CENTY) .OR. lpresent(O_CENTZ) .OR. &
         lpresent(O_KEEP)) THEN
         ier_num = 7
         ier_typ = ER_APPL
         ier_msg(1) = 'Only one common center and keep status will' 
         ier_msg(2) = 'apply to an accumulated set of hkls. Best use'
         ier_msg(3) = 'cent and keep with exec:hold only.'
         CALL errlist
         CALL no_error
      ENDIF
   ENDIF
ENDIF
!
CALL symm_restore                      ! Store current symmetry settings
!                                                                       
END SUBROUTINE boundary                       
!
!*****7*****************************************************************
!
SUBROUTINE boundarize_atom(center, distance, iatom, linside, &
                           surface_type, normal, thick) 
!-                                                                      
!     This subroutine sets the boundary property of the atom            
!                                                                       
!     Version  : 1.0                                                    
!     Date     : 07 Sep 10                                              
!                                                                       
!     Author   : R.B. Neder (reinhard.neder@fau.de)      
!+                                                                      
USE discus_config_mod 
USE discus_allocate_appl_mod
USE crystal_mod 
USE metric_mod
USE prop_para_mod 
USE surface_mod 
!
USE errlist_mod 
USE math_sup
use precision_mod
!                                                                       
IMPLICIT none 
       
!                                                                       
REAL(kind=PREC_DP)   , DIMENSION(3), INTENT(IN) :: center 
REAL(kind=PREC_DP)                 , INTENT(IN) :: distance 
INTEGER              , INTENT(IN) :: iatom 
LOGICAL              , INTENT(IN) :: linside 
INTEGER              , INTENT(IN) :: surface_type
REAL(kind=PREC_DP)   , DIMENSION(3), INTENT(IN) :: normal
REAL(kind=PREC_DP)                 , INTENT(IN) :: thick
!
REAL(kind=PREC_DP), DIMENSION(3), PARAMETER :: VNULL = (/ 0.0D0, 0.0D0, 0.0D0 /)
REAL(kind=PREC_DP)              , PARAMETER :: TOLERANCE = 5.0   ! Accept a 5 degree tilt for same surface
LOGICAL           , PARAMETER :: LSPACE = .FALSE.
!
INTEGER               :: idiv      ! Largest common divisor for the normal
INTEGER, DIMENSION(3) :: hkl
REAL(kind=PREC_DP)  , DIMENSION(3) :: rhkl, u, l_normal
REAL(kind=PREC_DP)                 :: angle, r
!                                                                       
IF( cr_nscat > SURF_MAXSCAT) THEN
   CALL alloc_surf ( cr_nscat )
   IF ( ier_num < 0 ) THEN
      RETURN
   ENDIF
ENDIF
!
r = do_blen(.FALSE., VNULL, normal)
l_normal(1) = NINT(10.*normal(1)/r)
l_normal(2) = NINT(10.*normal(2)/r)
l_normal(3) = NINT(10.*normal(3)/r)
!
IF((     linside.AND.distance <  0) .OR.  &
   (.not.linside.AND.distance >  0)      ) THEN                            
   cr_iscat (1,iatom) = 0 
   cr_prop (iatom) = ibclr (cr_prop (iatom), PROP_NORMAL) 
   cr_prop (iatom) = ibset (cr_prop (iatom), PROP_OUTSIDE) 
   IF ((thick<0 .AND. abs (distance) .lt.surf_ex_dist (cr_iscat (1,iatom) )) .OR. &
       (thick>0 .AND. ABS(distance) <     thick                          )) then 
      cr_prop (iatom) = ibset (cr_prop (iatom), PROP_SURFACE_EXT) 
   ENDIF 
   cr_surf (:,iatom) = 0
ELSE 
   IF ((thick<0 .AND. abs (distance) .lt.surf_ex_dist (cr_iscat (1,iatom) )) .OR. &
       (thick>0 .AND. ABS(distance) < thick                              )) THEN 
      cr_prop (iatom) = ibset (cr_prop (iatom), PROP_SURFACE_EXT) 
      IF(cr_surf(0, iatom) == SURF_NONE) THEN  ! Atom was not yet at a surface
         cr_surf(0,   iatom) = surface_type 
               
         cr_surf(1:3, iatom) = nint(l_normal(:))
      ELSEIF(cr_surf(0, iatom) < SURF_EDGE) THEN    ! Atom was already at a plane, sphere, cylinder wall
         IF(surface_type > SURF_CYLINDER) THEN      ! New location is at least an edge
            u(:) = cr_pos(:,iatom) - center(:)
            CALL pos2hkl(u, cr_gten, hkl)
            cr_surf(0,   iatom) = SURF_CORNER
            cr_surf(1:3, iatom)  = hkl(1:3)
         ELSEIF(surface_type > SURF_NONE) THEN      ! Should be all other cases
            rhkl(:) = cr_surf(1:3, iatom)           ! get old normal
            angle = do_bang(LSPACE, rhkl, VNULL, normal)
            IF(angle>TOLERANCE) THEN
                u(:) = cr_pos(:,iatom) - center(:)
                CALL pos2hkl(u, cr_gten, hkl)
                cr_surf(0,   iatom) = SURF_EDGE
                cr_surf(1:3, iatom)  = hkl(1:3)
            ENDIF
         ENDIF
      ELSE                                          ! Atom was at edge or corner
         u(:) = cr_pos(:,iatom) - center(:)
         CALL pos2hkl(u, cr_gten, hkl)
         cr_surf(0,   iatom) = SURF_CORNER
         cr_surf(1:3, iatom)  = hkl(1:3)
      ENDIF
      IF(.NOT.linside) cr_surf(:, iatom) = -cr_surf(:, iatom)  !invert for inside pointing surface
   ENDIF 
ENDIF 
IF(cr_surf(0, iatom) /= 0) THEN
   u(1:3) = cr_surf(1:3,iatom)
   CALL surface_normalize(u)
   cr_surf(1:3,iatom) = NINT(u(1:3))
   idiv = gcd(cr_surf(1, iatom), cr_surf(2, iatom), cr_surf(3, iatom))
   IF(idiv /= 0) THEN
       cr_surf(1:3, iatom) = cr_surf(1:3, iatom)/IABS(idiv)
   ENDIF
ENDIF
!                                                                       
END SUBROUTINE boundarize_atom                
!
!*****7*****************************************************************
!
SUBROUTINE prep_rotation(l_third, l_first, l_axis, o_third, o_first, o_axis, &
                         lo_third, lo_first, lo_axis,                        &
                         center, m_comb, m_combr)
!-
!  Prepares the rotation of the crystal
!+
USE crystal_mod
USE metric_mod
USE symm_mod
use symm_menu
USE symm_sup_mod
use symm_menu
!
USE errlist_mod
USE ber_params_mod
USE get_params_mod
USE matrix_mod
USE precision_mod
!
IMPLICIT NONE
!
LOGICAL         , INTENT(IN)    :: l_third
LOGICAL         , INTENT(IN)    :: l_first
LOGICAL         , INTENT(IN)    :: l_axis
CHARACTER(LEN=*), INTENT(INOUT) :: o_third
CHARACTER(LEN=*), INTENT(INOUT) :: o_first
CHARACTER(LEN=*), INTENT(INOUT) :: o_axis
INTEGER         , INTENT(INOUT) :: lo_third
INTEGER         , INTENT(INOUT) :: lo_first
INTEGER         , INTENT(INOUT) :: lo_axis
REAL(KIND=PREC_DP), DIMENSION(3), INTENT(IN) :: center
REAL(KIND=PREC_DP), DIMENSION(4,4), INTENT(OUT) :: m_comb
REAL(KIND=PREC_DP), DIMENSION(4,4), INTENT(OUT) :: m_combr
!
!INTEGER, PARAMETER :: MAXW = 3
REAL(KIND=PREC_DP), DIMENSION(3), PARAMETER :: V_NULL = (/0.0D0, 0.0D0, 0.0D0/)
REAL(KIND=PREC_DP)              , PARAMETER :: EPS    = 1.0D-4
REAL(KIND=PREC_DP)              , PARAMETER :: TOL    = 5.0D0
INTEGER :: i
REAL(KIND=PREC_DP), DIMENSION(3)   :: v_third
REAL(KIND=PREC_DP), DIMENSION(3)   :: v_first
REAL(KIND=PREC_DP), DIMENSION(4,4) :: m_third
REAL(KIND=PREC_DP), DIMENSION(4,4) :: m_first
!
REAL(KIND=PREC_DP), DIMENSION(3)   :: u         ! Dummy vector
REAL(KIND=PREC_DP), DIMENSION(3)   :: v         ! Dummy vector
REAL(KIND=PREC_DP), DIMENSION(4)   :: v4        ! Dummy vector
REAL(KIND=PREC_DP), DIMENSION(3)   :: ww        ! Dummy vector
REAL(KIND=PREC_DP), DIMENSION(3)   :: www       ! Dummy vector
REAL(KIND=PREC_DP)                 :: alpha
!
m_comb  = 0.0D0                     ! setup a default matrix
m_combr = 0.0D0                     ! setup a default matrix
DO i=1,4
   m_comb(i,i)  = 1.0D0
   m_combr(i,i) = 1.0D0
ENDDO
!
IF(l_third .and. .not.l_axis) THEN                    ! third axis is present
   CALL prep_values(o_third, lo_third, v_third, 'third')
   IF(ier_num/=0) RETURN
elseif(.not.l_third .and. l_axis) THEN                    ! axis is present
   CALL prep_values(o_axis , lo_axis , v_third, 'axis')
   IF(ier_num/=0) RETURN
elseif(l_third .and. l_axis) THEN                    ! axis is present
   ier_num = -6
   ier_typ = ER_COMM
   ier_msg(1) = 'Only one of ''third'' and ''axis'' may be present'
   return
ELSE
   v_third(1) = 0.0D0
   v_third(2) = 0.0D0
   v_third(3) = 1.0D0
ENDIF
!
m_third = 0.0D0
DO i=1,4
   m_third(i,i) = 1.0D0
ENDDO
!
IF(l_first) THEN                    ! first axis is present
   CALL prep_values(o_first, lo_first, v_first, 'first')
   IF(ier_num/=0) RETURN
ELSE
   u(1) = 1.0D0
   u(2) = 0.0D0
   u(3) = 0.0D0
!  CALL trans(u, cr_rten, v_first, 3)   ! transform b* into direct space
!  v_first = matmul(real(cr_rten,KIND=PREC_DP), u)
   v_first = matmul(     cr_rten              , u)
ENDIF
!
m_first = 0.0D0
DO i=1,4
   m_first(i,i) = 1.0D0
ENDDO
!
IF(l_third .AND. l_first) THEN
!  alpha = do_bang (.TRUE., REAL(v_first,kind=PREC_DP), V_NULL, real(v_third,kind=PREC_DP))
   alpha = do_bang (.TRUE.,      v_first              , V_NULL,      v_third              )
   IF(ABS(alpha-90.0D0)>TOL) THEN                 ! Not at 90 degrees
      ier_num = -171
      ier_typ = ER_APPL
      WRITE(ier_msg(1),'(a,f10.5)') 'Angle is ', alpha
      RETURN
   ENDIF
ENDIF
!
u(1) = 0.0D0
u(2) = 0.0D0
u(3) = 1.0D0
!alpha = do_bang(.TRUE., real(u,kind=PREC_DP), V_NULL, real(v_third,kind=PREC_DP))     ! Angle third to [001]
alpha = do_bang(.TRUE.,      u              , V_NULL, v_third                   )     ! Angle third to [001]
!
IF(ABS(alpha)<EPS) THEN                        ! already parallel
   alpha = 0.0D0                               ! Matrix will remain unit
   m_third = 0.0D0
   DO i=1,4
      m_third(i,i) = 1.0D0
   ENDDO
   ww(1) = 0.0D0
   ww(1) = 0.0D0
   ww(3) = 1.0D0
ELSEIF(ABS(180.0D0-alpha)<EPS) THEN             ! at 180 degrees
   alpha = 180.0D0
   v(1) = 0.0D0
   v(2) = 1.0D0                             ! Dummy b* vector
   v(3) = 0.0D0
   ww = matmul(     cr_rten              , v)
ELSE                                        ! Need a rotation
   CALL vekprod(     v_third              ,      u              , &
                     www             , cr_eps, cr_rten)  ! Get rotation axis
   ww = www
ENDIF
!
!  Define rotation
!
   CALL symm_reset
   sym_angle      = alpha                      ! Rotation angle
   sym_uvw        = ww                         ! Rotation axis
   sym_orig(:)    = 0.0                        ! Rotate around crystal center
   sym_trans(:)   = 0.0                        ! No translation
   sym_power      =  1                         ! Need just a single rotation
   sym_type       = .TRUE.                     ! Proper rotation
!  CALL trans (sym_uvw, cr_gten, sym_hkl, 3)   ! Make reciprocal space axis
!  sym_hkl = matmul(real(cr_gten,kind=PREC_DP), sym_uvw)
   sym_hkl = matmul(     cr_gten              , sym_uvw)
   CALL symm_setup                             ! get all matrices
   m_third   = sym_mat                          ! copy rotation matrix
!
!  CALL matinv4(m_third, m_longr)               ! Invert matrix for backwards operation
!
v4(1:3) = v_first(1:3)
v4(4)   = 0.00D0
v4 = MATMUL(m_third, v4)                        ! Rotate first axis by op for third axis
v(1:3) = v4(1:3)                               ! Copy current first axis
u(1) = 1.0D0
u(2) = 0.0D0
u(3) = 0.0D0
!alpha = do_bang(.TRUE., u, V_NULL, v)             ! Angle first to [100]
!alpha = do_bang(.TRUE., real(u,kind=PREC_DP), V_NULL, real(v_first,kind=PREC_DP))     ! Angle first to [001]
alpha = do_bang(.TRUE., u, V_NULL, v_first)    ! Angle first to [100]
CALL symm_reset
sym_angle      = alpha                         ! Rotation angle
sym_uvw(1)     = 0.0                           ! Rotation axis [001]
sym_uvw(2)     = 0.0                           ! Rotation axis [001]
sym_uvw(3)     = 1.0                           ! Rotation axis [001]
sym_orig(:)    = 0.0                           ! Rotate around crystal center
sym_trans(:)   = 0.0                           ! No translation
sym_power      =  1                            ! Need just a single rotation
sym_type       = .TRUE.                        ! Proper rotation
!CALL trans (sym_uvw, cr_gten, sym_hkl, 3)      ! Make reciprocal space axis
!sym_hkl = matmul(real(cr_gten,KIND=PREC_DP), sym_uvw)
sym_hkl = matmul(     cr_gten              , sym_uvw)
CALL symm_setup                                ! get all matrices
m_first   = sym_mat                            ! copy rotation matrix
m_comb = MATMUL(m_first, m_third)               ! Combine rotations m_first X m_third
!
v4(1:3) = center(1:3)                          ! Finally run computation for center
v4(4) = 0.00D0
v4 = MATMUL(m_comb, v4)                        ! Rotate center of Cyl/Ellipse
m_comb(1:3,4) = -v4(1:3)                       ! Place negative result as shift
!
CALL matinv4(m_comb, m_combr)                  ! Invert matrix for backwards operation
!
END SUBROUTINE prep_rotation
!
!*****7*****************************************************************
!
SUBROUTINE prep_values(o_third, lo_third, v_third, cname )
!
! Interpret the values for the optional parameter third / first
!
USE crystal_mod
!
USE errlist_mod
USE ber_params_mod
USE get_params_mod
USE precision_mod
!
IMPLICIT NONE
!
CHARACTER(LEN=*)  , INTENT(INOUT) :: o_third
INTEGER           , INTENT(INOUT) :: lo_third
REAL(KIND=PREC_DP),DIMENSION(3), INTENT(OUT) ::  v_third
CHARACTER(LEN=*)  , INTENT(IN) :: cname 
!
INTEGER, PARAMETER :: MAXW = 4
!
CHARACTER(LEN=MAX(PREC_STRING,LEN(o_third))), DIMENSION(MAXW) :: cpara
INTEGER            , DIMENSION(MAXW) :: lpara
INTEGER                              :: ianz 
REAL(KIND=PREC_DP), DIMENSION(MAXW)  :: werte
INTEGER :: length
LOGICAL :: ldirect
!
ldirect = .TRUE.
length = LEN_TRIM(cname)
IF(o_third(1:1) == '[' .AND. o_third(lo_third:lo_third) == ']') THEN
   o_third(1:1) = ' '
   o_third(lo_third:lo_third) = ' '
   CALL get_params(o_third, ianz, cpara, lpara, MAXW, lo_third) 
   IF(ier_num==0) THEN
      IF(ianz==4) THEN          ! 'd', 'r' is present 
         IF(cpara(4) == 'd') THEN
            ldirect = .TRUE.
         ELSEIF(cpara(4) == 'r') THEN
            ldirect = .FALSE.
         ELSE
            ier_num = -6
            ier_typ = ER_FORT
            ier_msg(1) = 'Optional parameter '//cname(1:length)//' failed '
            ier_msg(4) = 'Fourth paramter must be ''d'' or ''r'''
            RETURN
         ENDIF
         ianz = 3
         cpara(4) = ' '
      ENDIF
         IF(ianz==3) THEN
         CALL ber_params (ianz, cpara, lpara, werte, MAXW) 
         IF(ier_num==0) THEN
            v_third(1) = werte(1)
            v_third(2) = werte(2)
            v_third(3) = werte(3)
            IF(.NOT.ldirect) THEN
               v_third = matmul(real(cr_rten, kind=PREC_DP), werte(1:3))
            ENDIF
         ELSE
            ier_msg(1) = 'Optional parameter '//cname(1:length)//' failed '
            RETURN
         ENDIF
      ELSE
         RETURN
      ENDIF
   ELSE
      ier_msg(1) = 'Optional parameter '//cname(1:length)//' failed get params'
      RETURN
   ENDIF
ELSE
   ier_num = -9
   ier_typ = ER_FORT
   ier_msg(1) = 'Multiple optional values must be '
   ier_msg(2) = 'enclosed by []'
   ier_msg(3) = 'Optional parameter '//cname(1:length)
   RETURN
ENDIF
!
END SUBROUTINE prep_values
!
!*****7*****************************************************************
!
SUBROUTINE pos2hkl(u, gten, hkl)
!
USE metric_mod
use precision_mod
!
IMPLICIT NONE
!
REAL(kind=PREC_DP)   , DIMENSION(3)  , INTENT(INOUT) :: u
REAL(kind=PREC_DP)   , DIMENSION(3,3), INTENT(IN )   :: gten
INTEGER, DIMENSION(3)  , INTENT(OUT)   :: hkl
!
REAL(kind=PREC_DP)   , DIMENSION(3) :: v
REAL(kind=PREC_DP) :: uu
!
uu = SQRT(skalpro (u, u, gten) )
!     uu = SQRT(u(1)**2 + u(2)**2 + u(3)**2)  ! Really rough length
u(:) = u(:)/uu
!     CALL trans (u, gten, v, 3)              ! Transform into reciprocal space
v = matmul(gten, u)
hkl(1:3) = NINT(10.*v(1:3))             ! Make an approximate integer vector
!
END SUBROUTINE pos2hkl
!
!*****7*****************************************************************
!
SUBROUTINE do_surface_char(zeile, lp)
!
USE crystal_mod 
USE prop_para_mod
USE param_mod
USE prompt_mod
USE ber_params_mod
USE errlist_mod
USE get_params_mod
USE precision_mod
USE str_comp_mod
IMPLICIT NONE
CHARACTER(LEN=*) , INTENT(INOUT) :: zeile
INTEGER          , INTENT(INOUT) :: lp
!
INTEGER, PARAMETER :: MAXW = 3
CHARACTER(LEN=MAX(PREC_STRING,LEN(zeile))), DIMENSION(1:MAXW) :: cpara
INTEGER            , DIMENSION(1:MAXW) :: lpara
REAL(KIND=PREC_DP) , DIMENSION(1:MAXW) :: werte
INTEGER               :: surf_char
LOGICAL               :: lshow, lequal
INTEGER                 :: i
INTEGER, DIMENSION(3,6) :: surf_normal
INTEGER, DIMENSION(3)   :: surf_kante
INTEGER, DIMENSION(6)   :: surf_weight
!
INTEGER :: ianz, iatom
!
!
cpara(:) = ' '
lpara(:) = 1
CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
IF(ier_num/=0) RETURN
!
! Check if output is desired
!
lshow = .FALSE.
IF (str_comp(cpara(ianz) , 'show', 2, lpara(ianz) , 4) ) then 
  lshow = .TRUE.
  cpara(ianz) = ' '
  lpara(ianz) = 1
  ianz = ianz - 1
ENDIF
!
! Check if Atoms are restricted to equal or any atom tpye
!
lequal = .TRUE.
IF (str_comp(cpara(2) , 'equal', 2, lpara(ianz) , 5) ) then 
  lequal = .TRUE.
  cpara(2) = ' '
  lpara(2) = 1
  ianz = ianz - 1
ELSEIF (str_comp(cpara(2) , 'any', 2, lpara(ianz) , 3) ) then 
  lequal = .FALSE.
  cpara(2) = ' '
  lpara(2) = 1
  ianz = ianz - 1
ENDIF
!
! calculate atom number
!
CALL ber_params (ianz, cpara, lpara, werte, maxw) 
IF(ier_num/=0) RETURN
!
iatom = NINT(werte(1))
CALL surface_character(iatom, 1, cr_natoms, surf_char, surf_normal, surf_kante, surf_weight, lequal)
!
i = 0
res_para(0) = 1
res_para(1) = surf_char
res_para(2:4) = surf_normal(1:3,1)
!IF(surf_char==1) THEN                 ! Planar surface
!   res_para(0) = 4
!   res_para(2:4) = surf_normal(1:3,1)
!ELSEIF(surf_char==2) THEN             ! An edge
!   res_para(0)   = 10
!   res_para(2:4) = surf_normal(1:3,1)
!   res_para(5:7) = surf_normal(1:3,2)
!   res_para(8:10)= surf_kante(1:3)
!ELSEIF(surf_char==3) THEN             ! A corner
!   i = 0
!   fill: DO
!      IF(surf_weight(i+1) == 0) EXIT fill
!      i = i+1
!      res_para((i-1)*3+2:(i-1)*3+4) = surf_normal(1:3,i)
!   ENDDO fill
!   res_para(0) = 3*i+1
!ENDIF
IF(lshow) THEN
   IF(surf_char==0) THEN
      IF(IBITS(cr_prop(iatom),PROP_SURFACE_EXT,1).eq.1 .and.        &  ! real Atom is near surface
         IBITS(cr_prop(iatom),PROP_OUTSIDE    ,1).eq.0       ) THEN    ! real Atom is near surface
         WRITE(output_io, 2000) 
      ELSE
         WRITE(output_io, 2050)
      ENDIF
   ELSEIF(surf_char==1) THEN
      WRITE(output_io, 2100) 'Planar      outside face ',surf_normal(:,1)
   ELSEIF(surf_char==2) THEN
      WRITE(output_io, 2100) 'Spherical   outside face ',surf_normal(:,1)
   ELSEIF(surf_char==3) THEN
      WRITE(output_io, 2100) 'Cylindrical outside face ',surf_normal(:,1)
   ELSEIF(surf_char==4) THEN
      WRITE(output_io, 2100) 'Outside edge,  normal is ',surf_normal(:,1)
   ELSEIF(surf_char==5) THEN
      WRITE(output_io, 2100) 'Outside corner,normal is ',surf_normal(:,1)
   ELSEIF(surf_char==-1) THEN
      WRITE(output_io, 2100) 'Planar       inside face ',surf_normal(:,1)
   ELSEIF(surf_char==-2) THEN
      WRITE(output_io, 2100) 'Spherical    inside face ',surf_normal(:,1)
   ELSEIF(surf_char==-3) THEN
      WRITE(output_io, 2100) 'Cylindrical  inside face ',surf_normal(:,1)
   ELSEIF(surf_char==-4) THEN
      WRITE(output_io, 2100) ' Inside edge,  normal is ',surf_normal(:,1)
   ELSEIF(surf_char==-5) THEN
      WRITE(output_io, 2100) ' Inside corner,normal is ',surf_normal(:,1)
   ENDIF
ENDIF
!
2000 FORMAT(' Surface character could not be determined')
2050 FORMAT(' Atom is not close to a surface')
2100 FORMAT(a,3(I4,2x))
!2200 FORMAT(' Edge   surface,  Edge, dominant normal: ',3(I4,2x), 2x, 3(I4,2x), ':', I4)
!2250 FORMAT('                       secondary normal: ',20x         , 3(I4,2x), ':', I4)
!2300 FORMAT(' Corner surface,  dominant normal: ',3(I4,2x), ':', I4)
!2350 FORMAT('                 secondary normal: ',3(I4,2x), ':', I4)
!3100 FORMAT(' Indented plane,  normal: ',3(I4,2x))
!3200 FORMAT(' Indented edge ,  Edge, dominant normal: ',3(I4,2x), 2x, 3(I4,2x), ':', I4)
!3300 FORMAT(' Indented corner, dominant normal: ',3(I4,2x), ':', I4)
END SUBROUTINE do_surface_char
!
!*****7*****************************************************************
!
SUBROUTINE surface_character(iatom, jstart, jend, surf_char, surf_normal,  &
                             surf_kante, surf_weight, equal)
!
use atom_name
USE atom_env_mod
USE chem_mod
USE crystal_mod
USE do_find_mod
USE metric_mod
USE prop_para_mod
!
USE param_mod
USE precision_mod
USE prompt_mod
USE math_sup
!
IMPLICIT NONE
INTEGER,                 INTENT(IN)  :: iatom
INTEGER,                 INTENT(IN)  :: jstart
INTEGER,                 INTENT(IN)  :: jend
INTEGER,                 INTENT(OUT) :: surf_char
INTEGER, DIMENSION(3,6), INTENT(OUT) :: surf_normal
INTEGER, DIMENSION(3)  , INTENT(OUT) :: surf_kante
INTEGER, DIMENSION(6)  , INTENT(OUT) :: surf_weight
LOGICAL,                 INTENT(IN)  :: equal
!
!INTEGER, PARAMETER :: SURF_IN_PLANE  = -1
INTEGER, PARAMETER :: SURF_NONE   =  0
INTEGER, PARAMETER :: SURF_PLANE  =  1
INTEGER, PARAMETER :: SURF_EDGE   =  2
INTEGER, PARAMETER :: SURF_CORNER =  3
!
INTEGER, PARAMETER                     :: MAXW = 1
!LOGICAL, PARAMETER :: LNEW = .false.
REAL(kind=PREC_DP)   , PARAMETER                     :: RADIUS_MIN  = 1.51
REAL(kind=PREC_DP)   , PARAMETER                     :: RADIUS_STEP = 1.50
REAL(kind=PREC_DP)   , PARAMETER , DIMENSION(3)      :: NULLV = 0.0D0
REAL(kind=PREC_DP)   , PARAMETER                     :: IS_OUTSIDE  = 80.0
REAL(kind=PREC_DP)   , PARAMETER                     :: IS_PARALLEL = 15.0
!
CHARACTER(LEN=PREC_STRING), DIMENSION(1:MAXW) :: cpara
INTEGER            , DIMENSION(1:MAXW) :: lpara
REAL(KIND=PREC_DP) , DIMENSION(1:MAXW) :: werte
!
CHARACTER(LEN=PREC_STRING)     :: line
INTEGER  , DIMENSION(:), ALLOCATABLE :: neigh
REAL(kind=PREC_DP)     , DIMENSION(:), ALLOCATABLE :: angles
REAL(kind=PREC_DP)     , DIMENSION(:), ALLOCATABLE :: sorted
INTEGER  , DIMENSION(:,:), ALLOCATABLE :: surfaces
INTEGER                 :: nsurface
INTEGER                 :: i,j,k,l, ianz
INTEGER                 :: counter
INTEGER                 :: laenge
INTEGER                 :: neigsurf
INTEGER                 :: divisor
LOGICAL  , DIMENSION(3) :: fp
LOGICAL                 :: fq
LOGICAL                 :: lspace
LOGICAL                 :: isfound
LOGICAL                 :: isfirst
REAL(kind=PREC_DP)      :: rmin, radius
REAL(kind=PREC_DP)      :: alpha, beta
REAL(kind=PREC_DP)      :: dstar
REAL(kind=PREC_DP)      :: dist, dmin
REAL(kind=PREC_DP)     , DIMENSION(3) :: x         ! Vector from center to atom
REAL(kind=PREC_DP)     , DIMENSION(3) :: center    ! Average position of neighboring atoms
INTEGER                , DIMENSION(3) :: rough     ! rough normal 
REAL(kind=PREC_DP)     , DIMENSION(3) :: u,v,w     ! Vectors from central atom to neighbors
INTEGER                , DIMENSION(3) :: tempsurf  ! Vectors from central atom to neighbors
REAL(kind=PREC_DP)     , DIMENSION(3) :: realsurf  ! Vectors from central atom to neighbors
!
INTEGER, DIMENSION(0:1) :: temp_sel_prop = 0
!
temp_sel_prop(:) = cr_sel_prop(:)   ! save user settings for property select
cr_sel_prop(:)   = 0                ! Ignore all property selection rules 
surf_char = SURF_NONE
surf_normal = 0
surf_kante  = 0
surf_weight = 0
fp (1) = chem_period (1)
fp (2) = chem_period (2)
fp (3) = chem_period (3)
fq = chem_quick
rmin = 0.5
nsurface = 0
isfirst = .TRUE.
!
IF(IBITS(cr_prop(iatom),PROP_SURFACE_EXT,1).eq.1 .and.        &  ! real Atom is near surface
   IBITS(cr_prop(iatom),PROP_OUTSIDE    ,1).eq.0       ) THEN    ! real Atom is near surface
   IF(cr_surf(0,iatom)/=0) THEN     ! Atom has surface character set
      surf_char = cr_surf(0, iatom)
      surf_normal(:, 1) = cr_surf(1:3, iatom)
      surf_kante(:)     = 0
      surf_weight(:)    = 0
      surf_weight(1)    = 1
      RETURN                        ! Success go back
   ELSE
      RETURN
   ENDIF
ELSE
   RETURN
ENDIF
!OLD CODE TO BE WORKED ON
!
!
IF(IBITS(cr_prop(iatom),PROP_SURFACE_EXT,1).eq.1 .and.        &  ! real Atom is near surface
   IBITS(cr_prop(iatom),PROP_OUTSIDE    ,1).eq.0       ) THEN    ! real Atom is near surface
   cpara(1) = 'all'
   lpara(1) = 3
   ianz     = 1
   x(1)     = cr_pos(1,iatom)              ! Vector from center to atom
   x(2)     = cr_pos(2,iatom)
   x(3)     = cr_pos(3,iatom)
   WRITE(line,'(2(G16.8E3,a1),G16.8E3)') x(1), ',', x(2), ',', x(3)
   lspace = .TRUE.
   laenge = 50
   CALL d2r(line, laenge, lspace)
   rough(1:3) = INT(res_para(1:3))      ! Rough normal 
!
   radius   = RADIUS_MIN
   counter = 0                          ! Prevent infinit loop
   grand: DO                            ! increase search radius if character is not yet found
      counter  = counter + 1
      radius   = radius + RADIUS_STEP
      ianz     = 1
      werte(1) = -1                     ! Find all atom types
!
!     Determine a local center of mass with a larger search radius
!
      center(:) = 0.0
      CALL do_find_env (ianz, werte, maxw, x, rmin, radius+RADIUS_STEP, fq, fp)
      DO i=1, atom_env(0)               ! Pick out surface atom types only
         IF(jstart<=atom_env(i) .AND. atom_env(i)<=jend)  &
               center(:) = center(:) + cr_pos(:,atom_env(i))
      ENDDO
      center(:) = center(:) / atom_env(0)
!
!     Find local environment to construct the surface normal
!
      CALL do_find_env (ianz, werte, maxw, x, rmin, radius, fq, fp)
      ALLOCATE(neigh(0:atom_env(0)))
      neigh(:) = 0
      neigsurf = 0
      DO i=1, atom_env(0)               ! Pick out surface atom types only
         IF(IBITS(cr_prop(atom_env(i)),PROP_SURFACE_EXT,1).eq.1 .and.        &  ! real Atom is near surface
            IBITS(cr_prop(atom_env(i)),PROP_OUTSIDE    ,1).eq.0       ) THEN    ! real Atom is near surface
            IF(.NOT.equal .OR. (equal .AND. cr_iscat(1,atom_env(i))==cr_iscat(1,iatom)) ) THEN
               neigsurf = neigsurf + 1
               neigh(neigsurf) = atom_env(i)
            ENDIF
         ENDIF
      ENDDO
!
      IF(neigsurf >= 3 ) THEN             ! Found three surface atoms as neighbors
!        Determine rough surface as vector from center of mass to actual atom
         WRITE(line,'(2(G16.8E3,a1),G16.8E3)') x(1)-center(1), ',', x(2)-center(2), ',', x(3)-center(3)
         lspace = .TRUE.
         laenge = 50
         CALL d2r(line, laenge, lspace)
         rough(1:3) = INT(10*res_para(1:3))      ! Rough local normal 
!                                         ! Make space for all possible surfaces
         ALLOCATE(surfaces(0:3,neigsurf*neigsurf/2))
         nsurface      = 0
         surfaces(:,:) = 0
         indented: DO i=1, neigsurf-1               ! Loop over all neigbor pairs
            u(:) = cr_pos(:,neigh(i))-x(:)
            inner: DO j=i+1, neigsurf
               v(:) = cr_pos(:,neigh(j))-x(:)
               WRITE(line,2000) u,v
2000 FORMAT(6(G16.8E3,','),'ddr')
               laenge = 105
               CALL vprod(line, laenge)
               tempsurf(:) = NINT(res_para(1:3))  ! Normal to atom triplet
                                                  ! If necessary invert direction
               lspace = .FALSE.
               IF(do_blen(lspace, NULLV, REAL(tempsurf, kind=PREC_DP    ))>0) THEN
               lspace = .FALSE.
               IF(do_bang(lspace, REAL(rough, kind=PREC_DP), NULLV, &
                                  REAL(tempsurf, kind=PREC_DP)) > 90.0) tempsurf(:) = -tempsurf(:)
               WRITE(line,2010) tempsurf(:)
               laenge = 50
               lspace = .FALSE.
               CALL d2r(line, laenge, lspace)     ! Calculate direction in real space
2010 FORMAT(2(G16.8E3,','),G16.8E3) 
               realsurf(:) = res_para(1:3)
!              Test to see if any neighbor is further away than central atom
!
               DO k=1, neigsurf
                  v(:) = cr_pos(:,neigh(k))-x(:)
                  lspace = .TRUE.
                  IF(do_bang(lspace, realsurf, NULLV, v)<IS_OUTSIDE) THEN
                     CYCLE inner ! proceed to next pair
                  ENDIF
               ENDDO
!              Not an indented surface, proceed to sort
               isfound = .FALSE.
               lspace = .FALSE.
               DO k=1, nsurface                  ! Loop over all previous surfaces
                  alpha = do_bang(lspace, REAL(tempsurf, kind=PREC_DP), NULLV, REAL(surfaces(1:3,k), kind=PREC_DP) )
                  IF(alpha>90.) alpha = 180.-alpha
!
                  IF(alpha<IS_PARALLEL) THEN
                     surfaces(0,k) = surfaces(0,k) + 1   ! increment count of triplets with identical normal
                     isfound = .TRUE.
                  ENDIF
               ENDDO
!
               IF(.NOT.isfound) THEN            ! New surface normal
                  nsurface = nsurface + 1
                  surfaces(0,nsurface) = 1
                  surfaces(1:3,nsurface) = tempsurf
               ENDIF
               ENDIF
            ENDDO inner
         ENDDO indented
!        If only one surface, test if neighbors are on all sides
         IF(nsurface==1) THEN
            ALLOCATE(angles(1:neigsurf+1))
            ALLOCATE(sorted(1:neigsurf+1))
            angles(:) = 0.0
            sorted(:) = 0.0
            WRITE(line,2020) cr_pos(:,neigh(1))-x(:), surfaces(1:3,1)
            laenge = 106
2020 FORMAT(6(G16.8E3,','),'drdd') 
            CALL do_proj(line, laenge)
            u(:) = res_para(4:6)
            angles(1) = 0
            angles(neigsurf+1) = 360.0
            lspace = .TRUE.
            DO k=2,neigsurf
               WRITE(line,2020) cr_pos(:,neigh(k))-x(:), surfaces(1:3,1)
               laenge = 106
               CALL do_proj(line, laenge)
               v(:) = res_para(4:6)
               alpha = do_bang(lspace, u, NULLV, v)
               WRITE(line,2000) u,v
               laenge = 105
               CALL vprod(line, laenge)
               w(:) = res_para(1:3)
               IF(do_blen(lspace, NULLV, w) > 0.0001) THEN
                  beta = do_bang(lspace, w, NULLV, realsurf)
                  IF(beta.gt.90) alpha = 360. - alpha
               ENDIF      
               angles(k) = alpha
            ENDDO
            sorted(1) =   0.0
            angles(1) = 400.
            sorted(neigsurf+1) = 360.0
            angles(neigsurf+1) = 400.
            DO k=2,neigsurf
               i=MINLOC(angles,1)
               sorted(k) = angles(i)
               angles(i) = 400.0
            ENDDO
            alpha = 0.0
            DO k=2, neigsurf+1
               beta = sorted(k) - sorted(k-1)
               IF(beta > alpha) alpha = beta
            ENDDO
            DEALLOCATE(angles)
            DEALLOCATE(sorted)
            IF(alpha > 190) THEN       ! Current atom is at the corner of plane!
               nsurface = 3
            ELSEIF(alpha > 170) THEN   ! Current atom is at the edge of plane!
               nsurface = 2
            ENDIF
         ENDIF
         IF(nsurface>2 .AND. isfirst) THEN
            isfirst = .FALSE.
            DEALLOCATE(neigh)
            DEALLOCATE(surfaces)
            CYCLE grand
         ENDIF
         IF(nsurface>0) EXIT grand
         DEALLOCATE(neigh)
         DEALLOCATE(surfaces)
      ELSE                             ! Found less than three neigbors
!         radius = radius + RADIUS_STEP ! Increment search radius
         DEALLOCATE(neigh)
         IF(ALLOCATED(surfaces)) DEALLOCATE(surfaces)
      ENDIF
      IF(counter == 3) THEN 
         EXIT grand   ! Too many trials, give up
      ENDIF
   ENDDO grand
ENDIF
! Normalize
lspace = .false.
DO i=1, MIN(nsurface, UBOUND(surf_weight,1))
   dstar=do_blen(lspace, NULLV, REAL(surfaces(1:3,i), kind=PREC_DP))
   IF(dstar > 0) THEN
      surfaces(1:3,i) = NINT(REAL(surfaces(1:3,i), kind=PREC_DP)*10./dstar)
      divisor = IABS(gcd(surfaces(1,i),surfaces(2,i),surfaces(3,i)))
      surfaces(1:3,i) = surfaces(1:3,i)/divisor
   ENDIF
   surf_weight(i) = surfaces(0,i)
ENDDO
!
! might be used to differentiate external from internal surfaces ???
!
!w(1:3)=surfaces(1:3,1)
!beta = do_bang(lspace, w, NULL, x)
!write(*,*) ' normal, coord, angle ', w,x, beta
!
ianz = -1
CALL do_find_env (ianz, werte, maxw, x, rmin, RADIUS+RADIUS_STEP, fq, fp)
dmin = 1
DO i=1, nsurface
   DO k=1, atom_env(0)
   
      dist  = 1.0 - (cr_pos (1, atom_env(k))-x(1)) * surfaces (1,i) &
                  - (cr_pos (2, atom_env(k))-x(2)) * surfaces (2,i)  &
                  - (cr_pos (3, atom_env(k))-x(3)) * surfaces (3,i)
      dmin = MIN(dmin, dist)
   ENDDO
ENDDO
IF(nsurface == 0) THEN
   surf_char = SURF_NONE
ELSEIF(nsurface == 1) THEN
   surf_char        = SURF_PLANE
   surf_normal(:,1) = surfaces(1:3,1)
   surf_kante(:)    = 0
ELSEIF(nsurface == 2) THEN
   surf_char = SURF_EDGE
   i = MAXLOC(surf_weight,1)
   surf_normal(:,1) = surfaces(1:3,i)
   i = 3-i                          !Choose the other normal
   surf_normal(:,2) = surfaces(1:3,i)
   WRITE(line,2000) surfaces(1:3,1),surfaces(1:3,2)
   laenge = 105
   CALL vprod(line, laenge)
   surf_kante(:) = NINT(res_para(1:3)*100)
   dstar=do_blen(lspace, NULLV, REAL(surf_kante(1:3), kind=PREC_DP))
   surf_kante(1:3) = NINT(REAL(surf_kante(1:3), kind=PREC_DP)*10./dstar)
ELSEIF(nsurface >  2) THEN
   surf_char = SURF_CORNER
   fill:DO l = 1,6
      i = 1
      k=surfaces(0,1)
      DO j=2,nsurface
         IF(surfaces(0,j) > k ) THEN
            i = j
            k = surfaces(0,j)
         ENDIF
      ENDDO
      IF(k==0) EXIT fill
      surf_normal(:,l) = surfaces(1:3,i)
      surfaces(0,i) = -1
   ENDDO fill
   surf_kante(:)  = 0
ELSE
   surf_char = SURF_NONE
ENDIF
IF(dmin < 0) surf_char = -surf_char   ! Indented surface, return negative character
!
! Divide by greatest common divisor
!
k = 0
DO i=1,3
   k = MAX(k, IABS(surf_kante(i)))
ENDDO
IF(k>0) THEN
   divisor = gcd(surf_kante(1), surf_kante(2),surf_kante(3))
   surf_kante(:) = surf_kante(:)/divisor
ENDIF
IF(ALLOCATED(neigh)) DEALLOCATE(neigh)
IF(ALLOCATED(surfaces)) DEALLOCATE(surfaces)
cr_sel_prop(:) = temp_sel_prop(:)   ! restore user settings for property select
!
END SUBROUTINE surface_character
!
!*****7*****************************************************************
!
!*****7*****************************************************************
!
END MODULE surface_func_mod
