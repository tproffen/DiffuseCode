MODULE modify_func_mod
CONTAINS

!*****7*****************************************************************
      LOGICAL FUNCTION atom_allowed (i, werte, ianz, maxw) 
!+                                                                      
!     checks if atom i is within the selected atom range in             
!     werte(ianz).                                                      
!-                                                                      
      USE config_mod 
      USE crystal_mod 
      IMPLICIT none 
!                                                                       
       
      include'errlist.inc' 
!                                                                       
      INTEGER i, ianz, maxw 
      REAL werte (maxw) 
!                                                                       
      INTEGER j 
      LOGICAL ltype 
!                                                                       
      IF (i.le.0) then 
         ltype = .false. 
      ELSE 
         IF (werte (1) .eq. - 1) then 
            ltype = .true. 
         ELSE 
            ltype = .false. 
            DO j = 1, ianz 
            ltype = ltype.or.cr_iscat (i) .eq.nint (werte (j) ) 
            ENDDO 
         ENDIF 
      ENDIF 
      atom_allowed = ltype 
      END FUNCTION atom_allowed                     
!*****7*****************************************************************
      INTEGER FUNCTION bit_set (sel_field, col, bit, value) 
!                                                                       
      IMPLICIT none 
!                                                                       
      INTEGER sel_field (0:1) 
      INTEGER col 
      INTEGER bit 
      LOGICAL value 
!                                                                       
      IF (value) then 
         bit_set = IBSET (sel_field (col), bit) 
      ELSE 
         bit_set = IBCLR (sel_field (col), bit) 
      ENDIF 
!                                                                       
      END FUNCTION bit_set                          
!*****7*****************************************************************
      LOGICAL FUNCTION check_select_status (atom_status, property,      &
      sel_mask)                                                         
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
      USE prop_para_mod 
      IMPLICIT none 
!                                                                       
!                                                                       
      LOGICAL atom_status 
      INTEGER property 
      INTEGER sel_mask (0:1) 
!                                                                       
      INTEGER comp 
      INTEGER erg 
      INTEGER i, j 
!                                                                       
      IF (sel_mask (0) .eq.0) then 
         check_select_status = atom_status 
      ELSE 
         comp = 0 
         DO i = 0, MAXPROP - 1 
         IF (IBITS (sel_mask (1), i, 1) .eq.IBITS (property, i, 1) )    &
         then                                                           
            comp = IBSET (COMP, i) 
         ENDIF 
         ENDDO 
         erg = IAND (sel_mask (0), comp) 
!                                                                       
         check_select_status = atom_status.and.erg.eq.sel_mask (0) 
      ENDIF 
!                                                                       
      END FUNCTION check_select_status              
END MODULE modify_func_mod
