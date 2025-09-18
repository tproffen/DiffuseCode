MODULE get_iscat_mod
!
PRIVATE
PUBLIC get_iscat
!
CONTAINS
!
!*****7*****************************************************************
!
SUBROUTINE get_iscat (ianz, cpara, lpara, werte, maxw, lnew) 
!-                                                                      
!     Determines the scattering type of the parameter                   
!  Parameter come as either:
!  Al, Be, Cs ...
!  or as a single parameter 
!  (Al, Be, Cs)
!+                                                                      
USE discus_config_mod 
use discus_allocate_appl_mod
USE crystal_mod 
!
USE charact_mod
USE berechne_mod
use do_replace_expr_mod
USE do_variable_mod
USE errlist_mod 
USE get_params_mod
USE precision_mod
USE string_convert_mod
!                                                                       
IMPLICIT none 
!                                                                       
INTEGER                            , INTENT(INOUT) :: ianz    ! Number of atom types found
INTEGER                            , INTENT(IN)    :: maxw    ! ARRAY Dimensions
CHARACTER (LEN=*) , DIMENSION(MAXW), INTENT(INOUT) :: cpara   ! parameters
INTEGER           , DIMENSION(MAXW), INTENT(INOUT) :: lpara   ! parameter length
REAL(KIND=PREC_DP), DIMENSION(MAXW), INTENT(OUT)   :: werte   ! Resulting atom type numbers
LOGICAL                            , INTENT(IN)    :: lnew    ! Make a new atom type YES/NO
!
CHARACTER(LEN=MAX(PREC_STRING,LEN(cpara))) :: zeile 
INTEGER i, j, l, jj, jp 
LOGICAL  , DIMENSION(MAX(PREC_STRING,LEN(cpara)),0:1) :: lmask
logical :: lfound   ! did find old atom type
INTEGER :: omask
!                                                                       
INTEGER :: length
INTEGER :: ianz1
integer :: nscat   ! Dummy number of scatterers
CHARACTER(LEN=PREC_STRING), DIMENSION(MAXW) :: cpara1(maxw) 
INTEGER                   , DIMENSION(MAXW) :: lpara1(maxw) 
!
!  If the first parameter is of type (Al, Be, Cs), split this one
!
IF(cpara(1)(1:1) == '(' .AND. cpara(1)(lpara(1):lpara(1))==')') THEN
   zeile = cpara(1)(2:lpara(1)-1)
   length = lpara(1) - 2
   CALL get_params(zeile, ianz1, cpara1, lpara1, maxw, length)
   DO i=1, ianz1
      cpara(i) = cpara1(i)
      lpara(i) = lpara1(i)
   ENDDO
   ianz = ianz1
ENDIF
!                                                                       
!     ----Select which atoms are included in the wave                   
!                                                                       
lmask = .TRUE.
!
!     Attempt to replace a string from a variable
!
DO i=1, ianz
   zeile = ' '
   zeile = cpara(i)(1:lpara(i))
   j     = lpara(i)
   CALL ersetz_variable (zeile, j, lmask, omask)
   jj = LEN_TRIM(zeile)
!  IF(ier_num == 0 .AND. jj>2 .AND. zeile(1:1)=='''' .AND.           &
!                                   zeile(jj:jj) == ''''    ) THEN
   if(ier_num==0 .and. zeile/=cpara(i)(1:lpara(i))) then
      cpara(i) = ' '
      cpara(i)(1:jj) = zeile(1:jj)
   ENDIF
ENDDO
!
ier_num = 0 
ier_typ = ER_NONE 
werte   = 0.0            ! werte(:)
jj = 1 
jp = 0 
j = 1 
DO WHILE(j <= ianz.AND.ier_num == 0) 
   i = ichar (cpara (j) (1:1) ) 
   IF(cpara(j) == 'all') THEN 
      werte (1) = - 1 
      RETURN 
   ELSEIF(((a <= i .AND. i <= z) .OR. (aa <= i .AND. i <= zz) ) .AND.         &
          (index (cpara (j) , '[')  == 0)                             ) THEN                            
      CALL do_cap(cpara(j)) 
      ier_num = -27 
      ier_typ = ER_APPL 
      lfound = .false.
      place: DO i = 0, cr_nscat 
         IF (cpara(j)(1:lpara(j))  == cr_at_lis(i) ) THEN 
            lfound = .true.
            DO l=1,jj-1
               IF(NINT(werte(l)) == i) THEN
                  ier_num = 0 
                  ier_typ = ER_NONE 
                  CYCLE place    ! avoid double placement
               ENDIF
            ENDDO
            werte (jj) = i 
            jj = jj + 1 
            jp = jp + 1 
            ier_num = 0 
            ier_typ = ER_NONE 
         ENDIF 
      ENDDO place
      IF(lnew .and. .not.lfound) THEN 
!        IF(j >  1) THEN 
            if(cr_nscat==ubound(cr_at_lis,1)) then   ! Allocate more atoms
               nscat = cr_nscat + 1
               call alloc_crystal_scat(nscat)
            endif
            IF(cr_nscat < ubound(cr_at_lis,1)) THEN 
               cr_nscat = cr_nscat + 1 
               cr_at_lis(cr_nscat) = cpara(j) 
               cr_dw  (cr_nscat) = cr_dw  (1)
               cr_occ (cr_nscat) = cr_occ (1)
               cr_magn(:,cr_nscat) = cr_magn(:,1)
               cr_valu(cr_nscat) = cr_valu(1)
               cr_scat_int(i)      = .true.
               cr_scat_equ(i)      = .false.
               werte(jj) = cr_nscat 
               ier_num = 0 
               ier_typ = ER_NONE 
               jj = jj + 1
            ELSE 
               ier_num = - 26 
               ier_typ = ER_APPL 
            ENDIF 
!        ENDIF 
      ENDIF 
   ELSE 
      zeile = ' ' 
      l = lpara(j) 
      zeile(1:1) = '(' 
      zeile(2:l + 1) = cpara(j)(1:lpara(j))
      zeile(l + 2:l + 2) = ')' 
      l = l + 2 
      call do_replace_expr(zeile, l)
      werte(jj) = berechne(zeile, l) 
      IF(ier_num == 0) THEN 
         IF(0 <= NINT(werte(jj)) .AND. NINT(werte(jj))  <= cr_nscat) THEN                                          
            jj = jj + 1 
            ier_num = 0 
            ier_typ = ER_NONE 
         ELSE 
            ier_num = - 27 
            ier_typ = ER_APPL 
         ENDIF 
      ENDIF 
   ENDIF 
   j = j + 1 
ENDDO 
ianz = max (ianz, jp) 
!                                                                       
END SUBROUTINE get_iscat                      
!
!*******************************************************************************
!
END MODULE get_iscat_mod
