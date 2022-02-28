MODULE trans_cart_mod
!
CONTAINS
!
!*****7*****************************************************************
!
SUBROUTINE trans_atoms_tocart (uvw_out, &
           NMAX, cr_natoms, cr_pos, pl_tran_f)
!-                                                                      
!     transforms atom coordinates into a cartesian space                
!     Warning, only the fractional coordinates are transformed,         
!     the unit cell and space group information is not touched.         
!+                                                                      
!     USE discus_config_mod 
!     USE crystal_mod 
!     USE discus_plot_mod 
!
USE trans_sup_mod
!
use precision_mod
!
IMPLICIT none 
!                                                                       
REAL(kind=PREC_DP) ,DIMENSION(1:3), INTENT(OUT) :: uvw_out !(3)
INTEGER                   , INTENT(IN)    :: NMAX
INTEGER                   , INTENT(IN)    :: cr_natoms
REAL(kind=PREC_DP), DIMENSION(3,NMAX), INTENT(INOUT) :: cr_pos
REAL(kind=PREC_DP), DIMENSION(4,4)   , INTENT(IN   ) :: pl_tran_f
!
INTEGER              ::  i
LOGICAL, PARAMETER   :: lscreen = .false. 
REAL(kind=PREC_DP), DIMENSION(1:4) :: uvw
REAL(kind=PREC_DP)             :: xmin
REAL(kind=PREC_DP)             :: xmax
REAL(kind=PREC_DP)             :: ymin
REAL(kind=PREC_DP)             :: ymax
REAL(kind=PREC_DP)             :: zmin
REAL(kind=PREC_DP)             :: zmax
!                                                                       
xmin = 0.0D0
xmax = 0.0D0
ymin = 0.0D0
ymax = 0.0D0
zmin = 0.0D0
zmax = 0.0D0
uvw(4) = 1.0D0
!         
DO i = 1, cr_natoms 
   uvw (1) = cr_pos (1, i) 
   uvw (2) = cr_pos (2, i) 
   uvw (3) = cr_pos (3, i) 
   CALL tran_ca (uvw, pl_tran_f, lscreen) 
   cr_pos (1, i) = uvw (1) 
   cr_pos (2, i) = uvw (2) 
   cr_pos (3, i) = uvw (3) 
   xmin = MIN(xmin,uvw(1))
   xmax = MAX(xmax,uvw(1))
   ymin = MIN(ymin,uvw(2))
   ymax = MAX(ymax,uvw(2))
   zmin = MIN(zmin,uvw(3))
   zmax = MAX(zmax,uvw(3))
ENDDO
uvw_out (1) = ABS(xmax-xmin)
uvw_out (2) = ABS(ymax-ymin)
uvw_out (3) = ABS(zmax-zmin) 
!                                                                       
END SUBROUTINE trans_atoms_tocart      
!
!*****7*****************************************************************
!
SUBROUTINE trans_atoms_fromcart (NMAX, cr_natoms, cr_pos, pl_tran_fi)
!-                                                                      
!     transforms atom coordinates from a cartesian space back           
!     to the original coordinates                                       
!     Warning, only the fractional coordinates are transformed,         
!     the unit cell and space group information is not touched.         
!+                                                                      
!     USE discus_config_mod 
!     USE crystal_mod 
!     USE discus_plot_mod 
!
USE trans_sup_mod
use precision_mod
!
IMPLICIT none 
!
INTEGER                   , INTENT(IN)    :: NMAX
INTEGER                   , INTENT(IN)    :: cr_natoms
REAL(kind=PREC_DP)   , DIMENSION(3,NMAX), INTENT(INOUT) :: cr_pos
REAL(kind=PREC_DP)   , DIMENSION(4,4)   , INTENT(IN   ) :: pl_tran_fi
!                                                                       
INTEGER              :: i 
LOGICAL, PARAMETER   :: lscreen = .false.
!                                                                       
REAL(kind=PREC_DP), DIMENSION(1:4) ::  uvw !(4) 
!                                                                       
!                                                                       
uvw(4) = 1.0
DO i = 1, cr_natoms 
   uvw (1) = cr_pos (1, i) 
   uvw (2) = cr_pos (2, i) 
   uvw (3) = cr_pos (3, i) 
   CALL tran_ca (uvw, pl_tran_fi, lscreen) 
   cr_pos (1, i) = uvw (1) 
   cr_pos (2, i) = uvw (2) 
   cr_pos (3, i) = uvw (3) 
ENDDO 
!                                                                       
END SUBROUTINE trans_atoms_fromcart    
!
!*****7*****************************************************************
!
END MODULE trans_cart_mod
