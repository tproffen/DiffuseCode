MODULE celltoindex_mod
!
CONTAINS
!
!*****7*****************************************************************
!
SUBROUTINE celltoindex (icell, isite, iatom) 
!-                                                                      
!       calculates in which unit cell on which site the atom <ia> is    
!+
USE discus_config_mod 
USE crystal_mod 
USE errlist_mod 
!
IMPLICIT none 
!
INTEGER, DIMENSION(3), INTENT(IN)  :: icell
INTEGER              , INTENT(IN)  :: isite
INTEGER              , INTENT(OUT) :: iatom
!                                                                       
iatom = ( (icell(3) - 1) * cr_icc(1) * cr_icc(2) +       &
          (icell(2) - 1) * cr_icc(1) + (icell(1) - 1) )  &
        * cr_ncatoms + isite        
!
END SUBROUTINE celltoindex                    
!
!*****7*****************************************************************
!
SUBROUTINE indextocell (iatom, icell, isite) 
!-
!       calculates in which unit cell on which site the atom <ia> is    
!+
USE discus_config_mod 
USE crystal_mod 
USE errlist_mod 
!
IMPLICIT none 
!
INTEGER              , INTENT(IN)  :: iatom
INTEGER, DIMENSION(3), INTENT(OUT) :: icell
INTEGER              , INTENT(OUT) :: isite
!
INTEGER :: ia 
!
ia = iatom - 1 
!
icell(3) = INT(ia / cr_icc(1) / cr_icc(2) / cr_ncatoms) + 1
ia       = ia - (icell(3) - 1) * cr_icc (1) * cr_icc (2) * cr_ncatoms 
icell(2) = INT(ia / cr_icc(1) / cr_ncatoms) + 1 
ia       = ia - (icell(2) - 1) * cr_icc(1) * cr_ncatoms 
icell(1) = INT(ia / cr_ncatoms) + 1 
isite = ia - (icell(1) - 1) * cr_ncatoms + 1 
!
END SUBROUTINE indextocell                    
!
!*****7*****************************************************************
!
END MODULE celltoindex_mod
