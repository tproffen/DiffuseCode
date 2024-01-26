MODULE modify_func_mod
!
USE errlist_mod 
!
PRIVATE
PUBLIC atom_allowed
PUBLIC bit_set
PUBLIC check_select_status
!PUBLIC check_user_property
!
CONTAINS

!*****7*****************************************************************
LOGICAL FUNCTION atom_allowed (i, werte, ianz, maxw) 
!+                                                                      
!     checks if atom i is within the selected atom range in             
!     werte(ianz).                                                      
!-                                                                      
USE discus_config_mod 
USE crystal_mod 
USE precision_mod
IMPLICIT none 
!                                                                       
       
!                                                                       
INTEGER,                    INTENT(IN) :: i
INTEGER,                    INTENT(IN) :: ianz
INTEGER,                    INTENT(IN) :: maxw 
REAL(KIND=PREC_DP), DIMENSION(1:MAXW), INTENT(IN) :: werte
!                                                                       
INTEGER :: j 
LOGICAL :: ltype 
!                                                                       
IF (i.le.0) THEN 
   ltype = .FALSE. 
ELSE 
   IF(NINT(werte(1)) == -1) THEN 
      ltype = .TRUE. 
   ELSE 
      ltype = .FALSE. 
      DO j = 1, ianz 
         ltype = ltype .OR. cr_iscat(1,i)  == NINT(werte(j) ) 
      ENDDO 
   ENDIF 
ENDIF 
!
atom_allowed = ltype 
!
END FUNCTION atom_allowed                     
!
!*****7*****************************************************************
!
INTEGER FUNCTION bit_set (sel_field, col, bit, value) 
!                                                                       
IMPLICIT none 
!                                                                       
INTEGER, DIMENSION(0:1), INTENT(IN) :: sel_field (0:1) 
INTEGER                , INTENT(IN) :: col 
INTEGER                , INTENT(IN) :: bit 
LOGICAL                , INTENT(IN) :: value 
!                                                                       
IF (value) then 
   bit_set = IBSET(sel_field(col), bit) 
ELSE 
   bit_set = IBCLR(sel_field(col), bit) 
ENDIF 
!                                                                       
END FUNCTION bit_set                          
!*****7*****************************************************************
LOGICAL FUNCTION check_select_status (iatom, atom_status, property, sel_mask)
!-                                                                      
!     checks whether an atom was selected and whether the               
!     properties match to the selection mask                            
!                                                                       
!     for each property, the corresponding bit in "comp"                
!     is set, if the bits of the atom and row 1 of the selction         
!     mask are identical.                                               
!     "erg" is the logical "and" of "comp" and row 1 of the             
!     selection mask, i.e. the apply or ignore bit of the               
!     corresponding property.                                           
!     if "erg" is equal to row 1 of the selction mask, the atom         
!     fulfills the property requirements.                               
!+                                                                      
USE check_user_prop_mod
USE prop_para_mod 
!
IMPLICIT none 
!
INTEGER                             :: iatom
LOGICAL,                 INTENT(IN) :: atom_status 
INTEGER,                 INTENT(IN) :: property 
INTEGER, DIMENSION(0:1), INTENT(IN) :: sel_mask (0:1) 
!
INTEGER :: comp 
INTEGER :: erg 
INTEGER :: i
!LOGICAL check_user_property
!
IF (sel_mask (0) .eq.0) THEN 
   check_select_status = atom_status 
ELSE 
   comp = 0 
   DO i = 0, MAXPROP - 1 
      IF (IBITS(sel_mask(1), i, 1) == IBITS(property, i, 1)) THEN
         comp = IBSET (COMP, i) 
      ENDIF 
   ENDDO 
   erg = IAND (sel_mask (0), comp) 
!                                                                       
   check_select_status = atom_status.and.erg.eq.sel_mask (0) 
ENDIF 
!
IF(check_select_status) THEN
  check_select_status = check_select_status .AND. &
                        check_user_property(iatom)
ENDIF
!                                                                       
END FUNCTION check_select_status              
!
!*******************************************************************************
!
!LOGICAL FUNCTION check_user_property(iatom)
!
!USE crystal_mod
!USE conn_mod
!USE prop_para_mod
!
!IMPLICIT NONE
!
!INTEGER, INTENT(IN) :: iatom
!
!CHARACTER(LEN=256) :: c_name
!INTEGER            :: c_name_l
!INTEGER            :: ino
!INTEGER            :: i, ll
!INTEGER            :: natoms
!LOGICAL            :: test
!INTEGER, DIMENSION(:), ALLOCATABLE :: c_list
!INTEGER, DIMENSION(:,:), ALLOCATABLE :: c_offs
!!
!check_user_property = .TRUE.
!!
!test = .TRUE.
!check: DO i=1,prop_user_no
!   test = .TRUE.
!   IF(prop_user(i)%act/=0) THEN                               ! NOT ignore
!      test = test .AND. cr_iscat(iatom,1) == prop_user(i)%at_type  ! same atom type
!      ino  = prop_user(i)%conn_no
!      c_name   = prop_user(i)%conn_name
!      c_name_l = LEN_TRIM(c_name)
!      ll       = c_name_l
!      CALL get_connectivity_identity(cr_iscat(iatom,1), ino, c_name, c_name_l)
!      test = test .AND. prop_user(i)%conn_name(1:ll)==c_name(1:c_name_l)
!      CALL get_connectivity_list (iatom, cr_iscat(iatom,1), ino, c_list, c_offs, natoms )
!      test = test .AND. prop_user(i)%n_min<=natoms     &         ! Correct number of neighbors
!                  .AND.                     natoms<=prop_user(i)%n_max
!      test = test .AND. .NOT. (prop_user(i)%e_min<=natoms       &       ! Correct number of neighbors
!                               .AND.               natoms<=prop_user(i)%e_max)
!      IF(prop_user(i)%act==-1) test = .NOT.test               ! Absent invert the test
!   ENDIF
!   IF(.NOT.test) EXIT check
!ENDDO check
!check_user_property = test
!
!END FUNCTION check_user_property
!
!*******************************************************************************
!
END MODULE modify_func_mod
