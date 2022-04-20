MODULE discus_plot_export_mod
!
CONTAINS
!*****7*****************************************************************
      SUBROUTINE plot_atoms (iff) 
!-                                                                      
!     Writes the selected atoms for plotting with ATOMS                 
!+                                                                      
      USE discus_config_mod 
      USE crystal_mod 
      USE metric_mod
      USE modify_func_mod
      USE discus_plot_mod 
      USE errlist_mod 
USE lib_length
use matrix_mod
use precision_mod
      USE wink_mod
      IMPLICIT none 
!                                                                       
!                                                                       
      REAL(kind=PREC_DP) :: d, dist, fac 
      REAL(kind=PREC_DP), dimension(3) ::  v (3) 
      INTEGER i, j, iff 
      LOGICAL lno_slice, latom 
real(kind=PREC_DP), dimension(3) :: pl_uvw_local
pl_uvw_local = real(pl_uvw, kind=PREC_DP)     ! Prioir to compete migratin to DP
!                                                                       
!                                                                       
      WRITE (iff, 500) cr_name (1:len_str (cr_name) ) 
      WRITE (iff, 510) (cr_a0 (i) * cr_icc (i), i = 1, 3),  &
                       (cr_win (i), i = 1, 3)
      WRITE (iff, 520) 'P1' 
      WRITE (iff, 530) 
!                                                                       
      latom = .false. 
      fac = 1.0 / (8.0 * REAL(pi)**2) 
!                                                                       
      IF (pl_hkl(1) .eq.0.and.pl_hkl(2) .eq.0.and.pl_hkl(3) .eq.0) then
         lno_slice = .true. 
         d = 1.0 
      ELSE 
         lno_slice = .false. 
         d = sqrt (skalpro (pl_uvw_local, pl_uvw_local, cr_gten) ) 
         DO j = 1, 3 
         pl_mat (j, 1) = pl_abs (j) 
         pl_mat (j, 2) = pl_ord (j) 
         pl_mat (j, 3) = pl_uvw (j) 
         ENDDO 
!        CALL invmat (pl_inv, pl_mat) 
         call matinv3(pl_mat, pl_inv)
         IF (ier_num.eq. - 1) then 
            RETURN 
         ENDIF 
      ENDIF 
!                                                                       
      DO i = 1, cr_natoms 
!                                                                       
!     --Select atom if:                                                 
!       type has been selected and                                      
!            all atoms are selected                            or       
!                                                                       
      IF (check_select_status (i, pl_latom (cr_iscat (i) ), cr_prop (i),   &
      pl_sel_prop) ) then                                               
!                                                                       
!     --Check dimensions of plotting space                              
!                                                                       
         IF (pl_dim (1, 1) .le.cr_pos (1, i) .and.     &
             cr_pos (1, i) .le.pl_dim (1, 2) .and.     &
             pl_dim (2, 1) .le.cr_pos (2, i) .and.     &
             cr_pos (2, i) .le.pl_dim (2, 2) .and.     &
             pl_dim (3, 1) .le.cr_pos (3, i) .and.     &
             cr_pos (3, i) .le.pl_dim (3, 2) ) then
!                                                                       
!     ------Determine distance to point in plot slice                   
!                                                                       
            DO j = 1, 3 
            v (j) = cr_pos (j, i) - pl_vec (j) 
            ENDDO 
            dist = abs (skalpro (pl_uvw_local, v, cr_gten) ) / d 
!                                                                       
!     ------Write atom if inside slice or if slice has been deselected  
!                                                                       
            IF (dist.le.pl_width.or.lno_slice) then 
!                                                                       
!     ------write atom position or projection onto abscissa/ordinate    
!                                                                       
               DO j = 1, 3 
               v (j) = cr_pos (j, i) 
               ENDDO 
!                                                                       
!     ------Write atom position in desired sequence of coordinates      
!                                                                       
               latom = .true. 
               WRITE (iff, 1000) cr_at_lis (cr_iscat (i) ),          &
                 cr_iscat (i), (v (j) / cr_icc (j), j = 1, 3),       &
                 fac * cr_dw (cr_iscat (i) ), 0.0, 0.0, 0.0, 0.0, 0.0                            
            ENDIF 
         ENDIF 
      ENDIF 
      ENDDO 
!                                                                       
      IF (.not.latom) then 
         ier_num = - 58 
         ier_typ = ER_APPL 
      ENDIF 
!                                                                       
  500 FORMAT    ('TITL ',a) 
  510 FORMAT    ('CELL ',6(f10.6,1x)) 
  520 FORMAT    ('SPGP ',a) 
  530 FORMAT    ('FIELDS LAB TYP COO TFU') 
 1000 FORMAT     (a4,1x,i4,3x,3(f11.6,1x),/,6(f8.6,1x)) 
      END SUBROUTINE plot_atoms                     
!*****7*****************************************************************
      SUBROUTINE plot_drawxtl (iff) 
!-                                                                      
!     Writes the selected atoms for plotting with DRAWxtl               
!+                                                                      
      USE discus_config_mod 
      USE crystal_mod 
      USE metric_mod
      USE modify_func_mod
      USE discus_plot_mod 
      USE errlist_mod 
USE lib_length
use matrix_mod
use precision_mod
      USE wink_mod
      IMPLICIT none 
!                                                                       
!                                                                       
      REAL(kind=PREC_DP) :: d, dist, fac 
      REAL(kind=PREC_DP) :: v (3) 
      INTEGER i, j, iff 
      LOGICAL lno_slice, latom 
!                                                                       
      CHARACTER(12) povcolor (15) 
!                                                                       
real(kind=PREC_DP), dimension(3) :: pl_uvw_local
!
      DATA povcolor / 'Scarlet', 'HuntersGreen', 'MediumBlue',          &
      'Magenta', 'Yellow', 'Black', 'IndianRed', 'DarkGreen',           &
      'MidnightBlue', 'Maroon', 'Gold', 'Gray20', 'Cyan', 'SkyBlue',    &
      'White' /                                                         
!                                                                       
pl_uvw_local = real(pl_uvw, kind=PREC_DP)     ! Prioir to compete migratin to DP
!
      WRITE (iff, 501) cr_name (1:len_str (cr_name) ) 
      WRITE (iff, 511) (cr_a0 (i) * cr_icc (i), i = 1, 3),              &
                       (cr_win (i), i = 1, 3)
      WRITE (iff, 521) 'P 1' 
      WRITE (iff, 531) 
!                                                                       
      latom = .false. 
      fac = 1.0 / (8.0 * REAL(pi)**2) 
!                                                                       
      IF (pl_hkl(1) .eq.0.and.pl_hkl(2) .eq.0.and.pl_hkl(3) .eq.0) then
         lno_slice = .true. 
         d = 1.0 
      ELSE 
         lno_slice = .false. 
         d = sqrt (skalpro (pl_uvw_local, pl_uvw_local, cr_gten) ) 
         DO j = 1, 3 
         pl_mat (j, 1) = pl_abs (j) 
         pl_mat (j, 2) = pl_ord (j) 
         pl_mat (j, 3) = pl_uvw (j) 
         ENDDO 
!        CALL invmat (pl_inv, pl_mat) 
         call matinv3(pl_mat, pl_inv)
         IF (ier_num.eq. - 1) then 
            RETURN 
         ENDIF 
      ENDIF 
!                                                                       
      DO i = 1, cr_natoms 
!                                                                       
!     --Select atom if:                                                 
!       type has been selected and                                      
!            all atoms are selected                            or       
!                                                                       
      IF (check_select_status (i, pl_latom (cr_iscat (i) ), cr_prop (i),   &
      pl_sel_prop) ) then                                               
!                                                                       
!     --Check dimensions of plotting space                              
!                                                                       
         IF (pl_dim (1, 1) .le.cr_pos (1, i) .and.      &
             cr_pos (1, i) .le.pl_dim (1, 2) .and.      &
             pl_dim (2, 1) .le.cr_pos (2, i) .and.      &
             cr_pos (2, i) .le.pl_dim (2, 2) .and.      &
             pl_dim (3, 1) .le.cr_pos (3, i) .and.      &
             cr_pos (3, i) .le.pl_dim (3, 2) ) then  
!                                                                       
!     ------Determine distance to point in plot slice                   
!                                                                       
            DO j = 1, 3 
            v (j) = cr_pos (j, i) - pl_vec (j) 
            ENDDO 
            dist = abs (skalpro (pl_uvw_local, v, cr_gten) ) / d 
!                                                                       
!     ------Write atom if inside slice or if slice has been deselected  
!                                                                       
            IF (dist.le.pl_width.or.lno_slice) then 
!                                                                       
!     ------write atom position or projection onto abscissa/ordinate    
!                                                                       
               DO j = 1, 3 
               v (j) = cr_pos (j, i) 
               ENDDO 
!                                                                       
!     ------Write atom position in desired sequence of coordinates      
!                                                                       
               latom = .true. 
               WRITE (iff, 1001) cr_at_lis (cr_iscat (i) ),   &
                                 cr_iscat(i),(v(j)/cr_icc (j), j = 1, 3)
               WRITE (iff, 1002) cr_at_lis (cr_iscat (i) ),      &
                     cr_iscat (i),fac * cr_dw (cr_iscat (i) ),   &
                     fac * cr_dw (cr_iscat (i) ),                &
                     fac * cr_dw (cr_iscat (i) ), 0.0, 0.0, 0.0, &
                     povcolor ( pl_color (i) )
!     &                          pl_rgb(1,i),pl_rgb(2,i),pl_rgb(3,i);   
            ENDIF 
         ENDIF 
      ENDIF 
      ENDDO 
!                                                                       
      IF (.not.latom) then 
         ier_num = - 58 
         ier_typ = ER_APPL 
      ENDIF 
!                                                                       
      WRITE (iff, 541) 
                                                                        
  501 FORMAT    ('title ',a) 
  511 FORMAT    ('cell ',6(f10.6,1x)) 
  521 FORMAT    ('spgp ',a) 
  531 FORMAT    ('ellipsoids 50') 
  541 FORMAT    ('end') 
 1001 FORMAT     ('atom ',a4,1x,i4,3x,3(f11.6,1x)) 
!1002      format ('uij ',a4,1x,i4,3x,6(f8.6,1x),3f7.3)                 
 1002 FORMAT     ('uij ',a4,1x,i4,3x,6(f8.6,1x),a) 
      END SUBROUTINE plot_drawxtl                   
!
!*****7*****************************************************************
!
SUBROUTINE plot_cif (iff, lplot, do_spcgr) 
!-                                                                      
!     Writes the selected atoms as cif file (one unit cell)             
!     If lpolt is true, a unit cell for plotting is written,
!     else the unit cell is as small as possible.
!+                                                                      
USE discus_config_mod 
USE crystal_mod 
USE metric_mod
USE modify_func_mod
USE discus_plot_mod 
USE errlist_mod 
use matrix_mod
use precision_mod
USE wink_mod
use string_convert_mod
!
IMPLICIT none 
!
INTEGER, INTENT(IN) :: iff
LOGICAL, INTENT(IN) :: lplot
CHARACTER(LEN=*), INTENT(IN) :: do_spcgr
!                                                                       
REAL(kind=PREC_DP)   , PARAMETER   :: eightpisq = 8.*3.14159265**2
!                                                                       
!                                                                       
CHARACTER(LEN=16)     :: do_spcgr_w = 'P1'
character(len=4)      :: atom_i
REAL(kind=PREC_DP)    :: d, dist , shift
REAL(kind=PREC_DP),    DIMENSION(3) :: v
INTEGER, DIMENSION(3) :: scalef
INTEGER               :: i, j
LOGICAL               :: lno_slice, latom 
!                                                                       
real(kind=PREC_DP), dimension(3) :: pl_uvw_local
!
pl_uvw_local = real(pl_uvw, kind=PREC_DP)     ! Prioir to compete migratin to DP
!                                                                       
!     scalef(1) = MAX(cr_icc(1), INT((cr_dim(1,2)-cr_dim(1,1)))+2)
!     scalef(2) = MAX(cr_icc(2), INT((cr_dim(2,2)-cr_dim(2,1)))+2)
!     scalef(3) = MAX(cr_icc(3), INT((cr_dim(3,2)-cr_dim(3,1)))+2)
IF(lplot) THEN
   scalef(1) =                INT((cr_dim(1,2)-cr_dim(1,1)))+2
   scalef(2) =                INT((cr_dim(2,2)-cr_dim(2,1)))+2
   scalef(3) =                INT((cr_dim(3,2)-cr_dim(3,1)))+2
   shift     = 0.01
ELSE
   DO i=1,3
      IF(NINT(cr_dim(i,2)-cr_dim(i,1))-(cr_dim(i,2)-cr_dim(i,1))== 0.000) THEN
         scalef(i) = MAX(1,NINT((cr_dim(i,2)-cr_dim(i,1)))+1)
      ELSE
         scalef(i) =  INT((cr_dim(i,2)-cr_dim(i,1))) + 1
      ENDIF
   ENDDO
   shift     = 0.00
ENDIF
IF(do_spcgr=='original') THEN
   do_spcgr_w = cr_spcgr   ! needs a copy as it is called with a fixed string from plot
ENDIF
!
WRITE (iff, 500) 
WRITE(iff,'(a,a,a)') '_data_chemical_name_common ''',cr_name(1:LEN_TRIM(cr_name)),''''
WRITE(iff,*)
!
WRITE (iff, 510) (cr_a0 (i) * scalef (i), i = 1, 3),  &
                 (cr_win (i), i = 1, 3), do_spcgr_w
!                                                                       
latom = .false. 
!                                                                       
IF (pl_hkl(1).eq.0.and.pl_hkl(2).eq.0.and.pl_hkl(3).eq.0) then
   lno_slice = .true. 
   d = 1.0 
ELSE 
   lno_slice = .false. 
   d = sqrt (skalpro (pl_uvw_local, pl_uvw_local, cr_gten) ) 
   DO j = 1, 3 
         pl_mat (j, 1) = pl_abs (j) 
         pl_mat (j, 2) = pl_ord (j) 
         pl_mat (j, 3) = pl_uvw (j) 
   ENDDO 
!  CALL invmat (pl_inv, pl_mat) 
   call matinv3(pl_mat, pl_inv) 
   IF (ier_num.eq. - 1) then 
      RETURN 
   ENDIF 
ENDIF 
!                                                                       
DO i = 1, cr_natoms 
!                                                                       
!     --Select atom if:                                                 
!       type has been selected and                                      
!            all atoms are selected                            or       
!                                                                       
   IF(check_select_status(i, pl_latom(cr_iscat(i)), cr_prop(i), pl_sel_prop) ) THEN
!                                                                       
!     --Check dimensions of plotting space                              
!                                                                       
      IF (pl_dim (1, 1) .le.cr_pos (1, i) .and.        &
          cr_pos (1, i) .le.pl_dim (1, 2) .and.        &
          pl_dim (2, 1) .le.cr_pos (2, i) .and.        &
          cr_pos (2, i) .le.pl_dim (2, 2) .and.        &
          pl_dim (3, 1) .le.cr_pos (3, i) .and.        &
          cr_pos (3, i) .le.pl_dim (3, 2) ) then
!                                                                       
!     ------Determine distance to point in plot slice                   
!                                                                       
         DO j = 1, 3 
            v (j) = cr_pos (j, i) - pl_vec (j) 
         ENDDO 
         dist = abs (skalpro (pl_uvw_local, v, cr_gten) ) / d 
!                                                                       
!     ------Write atom if inside slice or if slice has been deselected  
!                                                                       
         IF (dist.le.pl_width.or.lno_slice) then 
!                                                                       
!     ------write atom position or projection onto abscissa/ordinate    
!                                                                       
            DO j = 1, 3 
               v (j) = cr_pos (j, i) + shift
            ENDDO 
!                                                                       
!     ------Write atom position in desired sequence of coordinates      
!                                                                       
            latom = .true. 
            atom_i = cr_at_lis(cr_iscat(i))
            call do_str(atom_i)               ! Remove non-character 
            WRITE (iff, 1000) atom_i,                                &
                ( (v (j) - cr_dim (j, 1) ) / scalef (j), j = 1, 3),  &
                 cr_dw ( cr_iscat (i) )/eightpisq
         ENDIF 
      ENDIF 
   ENDIF 
ENDDO 
!                                                                       
IF (.not.latom) then 
   ier_num = - 58 
   ier_typ = ER_APPL 
ENDIF 
!                                                                       
  500 FORMAT ('data_00001',/,                                           &
     &        '_audit_creation_method   ''DISCUS''',/)                  
  510 FORMAT ('_cell_length_a     ',f15.4,/                             &
     &        '_cell_length_b     ',f15.4,/                             &
     &        '_cell_length_c     ',f15.4,/                             &
     &        '_cell_angle_alpha  ',f10.4,/                             &
     &        '_cell_angle_beta   ',f10.4,/                             &
     &        '_cell_angle_gamma  ',f10.4,//                            &
     &        '_symmetry_space_group_name_H-M   ''',a,'''',//           &
     &        'loop_',/                                                 &
     &        '_atom_site_label',/                                      &
     &        '_atom_site_fract_x',/                                    &
     &        '_atom_site_fract_y',/                                    &
     &        '_atom_site_fract_z',/                                    &
     &        '_atom_site_u_iso_or_equiv')                              
 1000 FORMAT (a4,3x,3(f11.6,1x),4x,f8.6) 
!
END SUBROUTINE plot_cif                       
!
!*****7*****************************************************************
      SUBROUTINE plot_kuplot (iff, lkupl) 
!-                                                                      
!     Writes the selected atoms in a format suitable for plotting       
!     using KUPLOT or GNUPLOT                                           
!+                                                                      
      USE discus_config_mod 
      USE crystal_mod 
      USE metric_mod
      USE modify_func_mod
      USE discus_plot_mod 
      USE errlist_mod 
use matrix_mod
use precision_mod
!
      IMPLICIT none 
!                                                                       
!     INTEGER idim 
!     PARAMETER (idim = 3) 
!                                                                       
       
!                                                                       
      REAL(kind=PREC_DP) :: d, dist, ps 
      REAL(kind=PREC_DP) ::  v (3), u (3) 
      INTEGER i, j, iff, pt, pc 
      LOGICAL lno_slice, latom, lkupl 
real(kind=PREC_DP), dimension(3) :: pl_uvw_local
real(kind=PREC_DP), dimension(3,3) :: pl_inv_local
!
pl_uvw_local = real(pl_uvw, kind=PREC_DP)     ! Prioir to compete migratin to DP
!                                                                       
!                                                                       
      latom = .false. 
      IF (pl_hkl(1).eq.0.and.pl_hkl(2).eq.0.and.pl_hkl(3).eq.0) then
         lno_slice = .true. 
         d = 1.0 
      ELSE 
         lno_slice = .false. 
         d = sqrt (skalpro (pl_uvw_local, pl_uvw_local, cr_gten) ) 
         DO j = 1, 3 
         pl_mat (j, 1) = pl_abs (j) 
         pl_mat (j, 2) = pl_ord (j) 
         pl_mat (j, 3) = pl_uvw (j) 
         ENDDO 
!        CALL invmat (pl_inv, pl_mat) 
         call matinv3(pl_mat, pl_inv)
         IF (ier_num.eq. - 1) then 
            RETURN 
         ENDIF 
      ENDIF 
pl_inv_local = real(pl_inv, kind=PREC_DP)     ! Prioir to compete migratin to DP
!                                                                       
      DO i = 1, cr_natoms 
!                                                                       
!     --Select atom if:                                                 
!       type has been selected and                                      
!            all atoms are selected                            or       
!                                                                       
      IF (check_select_status (i, pl_latom (cr_iscat (i) ), cr_prop (i),   &
      pl_sel_prop) ) then                                               
!                                                                       
!     --Check dimensions of plotting space                              
!                                                                       
         IF (pl_dim (1, 1) .le.cr_pos (1, i) .and.          &
             cr_pos (1, i) .le.pl_dim (1, 2) .and.          &
             pl_dim (2, 1) .le.cr_pos (2, i) .and.          &
             cr_pos (2, i) .le.pl_dim (2, 2) .and.          &
             pl_dim (3, 1) .le.cr_pos (3, i) .and.          &
             cr_pos (3, i) .le.pl_dim (3, 2) ) then
!                                                                       
!     ------Determine distance to point in plot slice                   
!                                                                       
            DO j = 1, 3 
            v (j) = cr_pos (j, i) - pl_vec (j) 
            ENDDO 
            dist = abs (skalpro (pl_uvw_local, v, cr_gten) ) / d 
!                                                                       
!     ------Write atom if inside slice or if slice has been deselected  
!                                                                       
            IF (dist.le.pl_width.or.lno_slice) then 
!                                                                       
!     ------write atom position or projection onto abscissa/ordinate    
!                                                                       
               IF (lno_slice) then 
                  DO j = 1, 3 
                  v (j) = cr_pos (j, i) 
                  ENDDO 
               ELSE 
                  DO j = 1, 3 
                  u (j) = cr_pos (j, i) 
                  ENDDO 
!                 CALL trans (u, pl_inv_local, v, idim) 
                  v = matmul( pl_inv_local, u)
               ENDIF 
!                                                                       
!     ------Write atom position in desired sequence of coordinates      
!                                                                       
               latom = .true. 
               pt = pl_typ (cr_iscat (i) ) 
               pc = pl_color (cr_iscat (i) ) 
               ps = pl_siz (cr_iscat (i) ) 
!                                                                       
               IF (pl_col.eq.'xyz') then 
                  CALL write_atom (lkupl, iff, i, v (1), v (2), v (3),  &
                  pt, pc, ps)                                           
               ELSEIF (pl_col.eq.'yzx') then 
                  CALL write_atom (lkupl, iff, i, v (2), v (3), v (1),  &
                  pt, pc, ps)                                           
               ELSEIF (pl_col.eq.'zxy') then 
                  CALL write_atom (lkupl, iff, i, v (3), v (1), v (2),  &
                  pt, pc, ps)                                           
               ELSEIF (pl_col.eq.'zyx') then 
                  CALL write_atom (lkupl, iff, i, v (3), v (2), v (1),  &
                  pt, pc, ps)                                           
               ELSEIF (pl_col.eq.'yxz') then 
                  CALL write_atom (lkupl, iff, i, v (2), v (1), v (3),  &
                  pt, pc, ps)                                           
               ELSEIF (pl_col.eq.'xzy') then 
                  CALL write_atom (lkupl, iff, i, v (1), v (3), v (2),  &
                  pt, pc, ps)                                           
               ENDIF 
            ENDIF 
         ENDIF 
      ENDIF 
      ENDDO 
!                                                                       
      IF (.not.latom) then 
         ier_num = - 58 
         ier_typ = ER_APPL 
      ENDIF 
!                                                                       
      END SUBROUTINE plot_kuplot                    
!*****7*****************************************************************
      SUBROUTINE plot_kuplot_mol (iff, lkupl) 
!-                                                                      
!     Writes the selected molecules in a format suitable                
!     for plotting using KUPLOT or GNUPLOT                              
!+                                                                      
      USE discus_config_mod 
      USE crystal_mod 
      USE metric_mod
      USE molecule_mod 
      USE discus_plot_mod 
      USE errlist_mod 
use matrix_mod
use precision_mod
      IMPLICIT none 
!                                                                       
!     INTEGER idim 
!     PARAMETER (idim = 3) 
!                                                                       
       
!                                                                       
      REAL(kind=PREC_DP) :: d, dist, ps 
      REAL(kind=PREC_DP) ::  v (3), u (3) 
      INTEGER i, j, k, i0, iff, pt, pc 
      INTEGER i_start, i_end, imol 
      LOGICAL lno_slice, latom, lkupl 
real(kind=PREC_DP), dimension(3) :: pl_uvw_local
real(kind=PREC_DP), dimension(3,3) :: pl_inv_local
!                                                                       
!                                                                       
pl_uvw_local = real(pl_uvw, kind=PREC_DP)     ! Prioir to compete migratin to DP
      latom = .false. 
      IF (pl_hkl(1).eq.0.and.pl_hkl(2).eq.0.and.pl_hkl(3).eq.0) then
         lno_slice = .true. 
         d = 1.0 
      ELSE 
         lno_slice = .false. 
         d = sqrt (skalpro (pl_uvw_local, pl_uvw_local, cr_gten) ) 
         DO j = 1, 3 
         pl_mat (j, 1) = pl_abs (j) 
         pl_mat (j, 2) = pl_ord (j) 
         pl_mat (j, 3) = pl_uvw (j) 
         ENDDO 
!        CALL invmat (pl_inv, pl_mat) 
         call matinv3(pl_mat, pl_inv)
         IF (ier_num.eq. - 1) then 
            RETURN 
         ENDIF 
      ENDIF 
!
pl_inv_local = real(pl_inv, kind=PREC_DP)     ! Prioir to compete migratin to DP
!                                                                       
      DO i = 1, mole_num_mole 
!                                                                       
!     --Select molecule if:                                             
!       type has been selected                                          
!                                                                       
      IF (pl_latom (mole_type (i) ) ) then 
!                                                                       
!     --Check dimensions of plotting space                              
!                                                                       
         i0 = mole_cont (mole_off (i) + 1) 
         IF (pl_dim (1, 1) .le.cr_pos (1, i0) .and.             &
             cr_pos (1, i0) .le.pl_dim (1, 2) .and.             &
             pl_dim (2, 1) .le.cr_pos (2, i0) .and.             &
             cr_pos (2, i0) .le.pl_dim (2, 2) .and.             &
             pl_dim (3, 1) .le.cr_pos (3, i0) .and.             &
             cr_pos (3, i0) .le.pl_dim (3, 2) ) then
!                                                                       
!     ------Determine distance to point in plot slice                   
!------ ------This is done only for the molecules origin                
!                                                                       
            DO j = 1, 3 
            v (j) = cr_pos (j, i0) - pl_vec (j) 
            ENDDO 
            dist = abs (skalpro (pl_uvw_local, v, cr_gten) ) / d 
!                                                                       
!     ------Write atom if inside slice or if slice has been deselected  
!                                                                       
            IF (dist.le.pl_width.or.lno_slice) then 
!                                                                       
!     ------Write atom position in desired sequence of coordinates      
!                                                                       
               latom = .true. 
               pt = pl_typ (mole_type (i) ) 
               pc = pl_color (mole_type (i) ) 
               ps = pl_siz (mole_type (i) ) 
!                                                                       
               IF (pl_mol_all) then 
                  i_start = 1 
                  i_end = mole_len (i) 
               ELSE 
                  i_start = 1 
                  i_end = 1 
               ENDIF 
!                                                                       
               DO imol = i_start, i_end 
               j = mole_cont (mole_off (i) + imol) 
!                                                                       
!     ------write atom position or projection onto abscissa/ordinate    
!                                                                       
               IF (lno_slice) then 
                  DO k = 1, 3 
                  v (k) = cr_pos (k, j) 
                  ENDDO 
               ELSE 
                  DO k = 1, 3 
                  u (k) = cr_pos (k, j) 
                  ENDDO 
!                 CALL trans (u, pl_inv_local, v, idim) 
                  v = matmul(pl_inv_local, u)
               ENDIF 
               IF (pl_col.eq.'xyz') then 
                  CALL write_atom (lkupl, iff, j, v (1), v (2), v (3),  &
                  pt, pc, ps)                                           
               ELSEIF (pl_col.eq.'yzx') then 
                  CALL write_atom (lkupl, iff, j, v (2), v (3), v (1),  &
                  pt, pc, ps)                                           
               ELSEIF (pl_col.eq.'zxy') then 
                  CALL write_atom (lkupl, iff, j, v (3), v (1), v (2),  &
                  pt, pc, ps)                                           
               ELSEIF (pl_col.eq.'zyx') then 
                  CALL write_atom (lkupl, iff, j, v (3), v (2), v (1),  &
                  pt, pc, ps)                                           
               ELSEIF (pl_col.eq.'yxz') then 
                  CALL write_atom (lkupl, iff, j, v (2), v (1), v (3),  &
                  pt, pc, ps)                                           
               ELSEIF (pl_col.eq.'xzy') then 
                  CALL write_atom (lkupl, iff, j, v (1), v (3), v (2),  &
                  pt, pc, ps)                                           
               ENDIF 
               ENDDO 
            ENDIF 
         ENDIF 
      ENDIF 
      ENDDO 
!                                                                       
      IF (.not.latom) then 
         ier_num = - 58 
         ier_typ = ER_APPL 
      ENDIF 
!                                                                       
      END SUBROUTINE plot_kuplot_mol                
!*****7*****************************************************************
SUBROUTINE write_atom (lkupl, iff, iatom, x, y, z, pt, pc, ps) 
!                                                                       
USE discus_config_mod 
USE chem_mod
USE crystal_mod 
USE celltoindex_mod
USE discus_plot_mod 
!                                                                       
IMPLICIT none 
!                                                                       
REAL(kind=PREC_DP)    :: x, y, z, ps 
INTEGER :: cr_end 
INTEGER :: i, pt, pc, iff, iatom, isite, icell (3) 
LOGICAL :: lkupl 
!                                                                       
cr_end = cr_ncatoms*cr_icc(1) * cr_icc(2)*cr_icc(3) + 1
!                                                                       
IF (pl_dens) then 
   IF(chem_quick) THEN      !Fast chem mode
      IF (iatom.lt.cr_end) then 
         CALL indextocell (iatom, icell, isite) 
      ELSE 
         DO i = 1, 3 
         icell (i) = int (cr_pos (i, iatom) - cr_dim0 (i, 1) ) + 1
         ENDDO 
      ENDIF 
      x = x - REAL(icell (1) - 1) - cr_dim0 (1, 1) 
      y = y - REAL(icell (2) - 1) - cr_dim0 (2, 1) 
      z = z - REAL(icell (3) - 1) - cr_dim0 (3, 1) 
   ELSE
      x = x - REAL(INT(x))
      y = y - REAL(INT(y))
      z = z - REAL(INT(z))
      x = x+1.25 - REAL(INT(x+1.25)) - 0.25
      y = y+1.25 - REAL(INT(y+1.25)) - 0.25
      z = z+1.25 - REAL(INT(z+1.25)) - 0.25
   ENDIF 
ENDIF 
!                                                                       
IF (lkupl) then 
   WRITE (iff, 1000) x, y, z, pt, pc, ps 
ELSE 
   WRITE (iff, 2000) x, y, z 
ENDIF 
!                                                                       
 1000 FORMAT  (3(2x,f12.6),2(2x,i2),2x,f6.2) 
 2000 FORMAT  (3(2x,f12.6)) 
!
END SUBROUTINE write_atom                     
!
!*****7*****************************************************************
!
      SUBROUTINE plot_xbs (iff) 
!-                                                                      
!     Writes the selected atoms in a format suitable for plotting       
!     using xbs                                                         
!+                                                                      
      USE discus_config_mod 
      USE crystal_mod 
      USE metric_mod
      USE modify_func_mod
      USE discus_plot_mod 
      USE discus_plot_init_mod
      USE trans_sup_mod
      USE errlist_mod 
use matrix_mod
      USE param_mod 
      USE precision_mod
      IMPLICIT none 
!                                                                       
integer, intent(in)  :: iff
!                                                                       
      CHARACTER(LEN=PREC_STRING) :: zeile 
      INTEGER :: laenge 
      INTEGER :: i, j, k! , iff 
      LOGICAL :: latom, lspace 
      LOGICAL :: lscreen 
REAL(kind=PREC_DP), dimension(4) :: uvw (4) 
REAL(kind=PREC_DP), dimension(3) :: absz (3) 
REAL(kind=PREC_DP), dimension(3) :: xmin (3), xmax (3)
REAL(kind=PREC_DP)               :: xx 
real(kind=PREC_DP), dimension(3) :: pl_uvw_local
real(kind=PREC_DP), dimension(3) :: pl_abs_local
real(kind=PREC_DP), dimension(3), PARAMETER :: NULLV = (/ 0.0D0, 0.0D0, 0.0D0 /)
!                                                                       
!                                                                       
      DATA lscreen / .false. / 
      DATA lspace / .true. / 
!                                                                       
      CALL plot_ini_trans (1.0D0,                                &
                 pl_tran_g, pl_tran_gi, pl_tran_f, pl_tran_fi, &
                 cr_gten, cr_rten, cr_eps)
!
pl_uvw_local = real(pl_uvw, kind=PREC_DP)     ! Prioir to compete migratin to DP
pl_abs_local = real(pl_abs, kind=PREC_DP)     ! Prioir to compete migratin to DP
!                                                                       
      latom = .false. 
      uvw (4) = 1.0 
      DO j = 1, 3 
      xmin (j) = 0.0 
      xmax (j) = 0.0 
      ENDDO 
      xx = 0.0 
!                                                                       
      DO i = 1, cr_natoms 
!                                                                       
!     --Select atom if:                                                 
!       type has been selected and                                      
!            all atoms are selected                            or       
!                                                                       
      IF (check_select_status (i, pl_latom (cr_iscat (i) ), cr_prop (i),   &
      pl_sel_prop) ) then                                               
!                                                                       
!     --Check dimensions of plotting space                              
!                                                                       
         IF (pl_dim (1, 1) .le.cr_pos (1, i) .and.      &
             cr_pos (1, i) .le.pl_dim (1, 2) .and.      &
             pl_dim (2, 1) .le.cr_pos (2, i) .and.      &
             cr_pos (2, i) .le.pl_dim (2, 2) .and.      &
             pl_dim (3, 1) .le.cr_pos (3, i) .and.      &
             cr_pos (3, i) .le.pl_dim (3, 2) ) then

!                                                                       
!     ------write atom position                                         
!                                                                       
            latom = .true. 
            DO j = 1, 3 
            uvw (j) = cr_pos (j, i) 
            ENDDO 
            CALL tran_ca (uvw, pl_tran_f, lscreen) 
            DO j = 1, 3 
            xmin (j) = min (xmin (j), uvw (j) ) 
            xmax (j) = max (xmax (j), uvw (j) ) 
            ENDDO 
            WRITE (iff, 2100) cr_at_lis(cr_iscat(i)),(uvw(j), j = 1, 3)
         ENDIF 
      ENDIF 
      ENDDO 
!                                                                       
!     determine scale                                                   
!                                                                       
      IF (pl_scale (0) .eq. - 1.0) then 
         xx = max (abs (xmax (1) - xmin (1) ), abs (xmax (2) - xmin (1))&
                 , abs (xmax (3) - xmin (1) ) )                                
         pl_scale (1) = 200. / xx 
      ENDIF 
!                                                                       
!     write Atom specifications                                         
!                                                                       
      WRITE (iff, * ) 
      DO i = 0, cr_nscat 
      IF (pl_latom (i) ) then 
         WRITE (iff, 2200) cr_at_lis (i), pl_siz (i),                   &
                          (pl_rgb (j, i), j = 1, 3)
      ENDIF 
      ENDDO 
!                                                                       
!     Check that normal and asbzissa are not parallel                   
!                                                                       
      xx = do_bang (lspace, pl_uvw_local, nullv, pl_abs_local) 
      IF (xx.eq.0.0) then 
         ier_num = - 80 
         ier_typ = ER_APPL 
         RETURN 
      ENDIF 
!                                                                       
      WRITE (zeile, 1000) pl_abs, pl_uvw 
      laenge = 82 
      CALL do_proj (zeile, laenge) 
      absz (1) = res_para (4) 
      absz (2) = res_para (5) 
      absz (3) = res_para (6) 
!                                                                       
      xx = do_blen (lspace, pl_uvw_local, nullv) 
      IF (xx.eq.0.0) then 
         ier_num = - 32 
         ier_typ = ER_APPL 
         RETURN 
      ENDIF 
      DO j = 1, 3 
      uvw (j) = pl_uvw (j) / xx 
      ENDDO 
      CALL tran_ca (uvw, pl_tran_f, lscreen) 
      DO j = 1, 3 
      pl_tran_gi (3, j) = uvw (j) 
      ENDDO 
!                                                                       
      xx = do_blen (lspace, absz, nullv) 
      IF (xx.eq.0.0) then 
         ier_num = - 32 
         ier_typ = ER_APPL 
         RETURN 
      ENDIF 
      DO j = 1, 3 
      uvw (j) = absz (j) / xx 
      ENDDO 
      CALL tran_ca (uvw, pl_tran_f, lscreen) 
      DO j = 1, 3 
      pl_tran_gi (1, j) = uvw (j) 
      ENDDO 
      pl_tran_gi (2, 1) = pl_tran_gi (3, 2) * pl_tran_gi (1, 3) -       &
      pl_tran_gi (3, 3) * pl_tran_gi (1, 2)                             
      pl_tran_gi (2, 2) = pl_tran_gi (3, 3) * pl_tran_gi (1, 1) -       &
      pl_tran_gi (3, 1) * pl_tran_gi (1, 3)                             
      pl_tran_gi (2, 3) = pl_tran_gi (3, 1) * pl_tran_gi (1, 2) -       &
      pl_tran_gi (3, 2) * pl_tran_gi (1, 1)                             
!                                                                       
      DO i = 1, 4 
      DO j = 1, 4 
      pl_tran_g (i, j) = pl_tran_gi (i, j) 
      ENDDO 
      ENDDO 
!                                                                       
!     CALL invmat4 (pl_tran_g) 
      call matinv(pl_tran_gi, pl_tran_g)
!                                                                       
!     write Bond types                                                  
!                                                                       
      WRITE (iff, * ) 
      DO i = 1, cr_nscat 
      DO j = 1, cr_nscat 
      IF (pl_bond (i, j) ) then 
         WRITE (iff, 2300) cr_at_lis (i), cr_at_lis (j),                &
               pl_bond_len (1, i, j), pl_bond_len (2, i, j),            &
               pl_bond_rad (i, j), (pl_bond_col (k, i, j), k = 1, 3)
      ENDIF 
      ENDDO 
      ENDDO 
!                                                                       
!     write trailer, for right now a standard                           
!                                                                       
      WRITE (iff, * ) 
      WRITE (iff, 2400) ( (pl_tran_g (i, j), i = 1, 3), j = 1, 3) 
      WRITE (iff, 2500) 12.0 
      WRITE (iff, 2600) 1.0 
      WRITE (iff, 2700) pl_scale (1) 
      WRITE (iff, 2800) 1.0 
      WRITE (iff, 2900) 1.0 
      WRITE (iff, 3000) 0.00, 0.00 
      WRITE (iff, 3100) 1, 0, 1, 0, 0, 1, 0, 0, 0 
!                                                                       
      IF (.not.latom) then 
         ier_num = - 58 
         ier_typ = ER_APPL 
      ENDIF 
!                                                                       
 1000 FORMAT    (6(f12.6,','),'dddd') 
 2100 FORMAT    ('atom      ',a4,3f11.3) 
 2200 FORMAT    ('spec      ',a4,f10.3,3f7.2) 
 2300 FORMAT    ('bonds     ',a4,4x,a4,6f9.3) 
 2400 FORMAT    ('tmat',9f7.3) 
 2500 FORMAT    ('dist   ',f7.3) 
 2600 FORMAT    ('inc    ',f7.3) 
 2700 FORMAT    ('scale  ',f7.3) 
 2800 FORMAT    ('rfac   ',f7.3) 
 2900 FORMAT    ('bfac   ',f7.3) 
 3000 FORMAT    ('pos    ',2f7.3) 
 3100 FORMAT    ('switches',9i2) 
!                                                                       
      END SUBROUTINE plot_xbs                       
!*****7*****************************************************************
      SUBROUTINE plot_frames (iff) 
!-                                                                      
!     Writes the selected atoms in a format suitable for plotting       
!     using xbs                                                         
!     this subroutine writes the frames for a movie                     
!+                                                                      
      USE discus_config_mod 
      USE crystal_mod 
!     USE metric_mod
      USE modify_func_mod
      USE discus_plot_mod 
      USE discus_plot_init_mod
      USE trans_sup_mod
      USE errlist_mod 
      USE param_mod 
use precision_mod
!
      IMPLICIT none 
!                                                                       
!                                                                       
      INTEGER i, j, iff 
      LOGICAL latom, lspace 
      LOGICAL lscreen 
REAL(kind=PREC_DP), dimension(4) :: uvw (4) 
REAL(kind=PREC_DP), dimension(3) :: xmin (3), xmax (3)!, nullv (3) 
REAL(kind=PREC_DP)               :: xx 
!                                                                       
      DATA lscreen / .false. / 
      DATA lspace / .true. / 
!     DATA nullv / 0.0, 0.0, 0.0 / 
!                                                                       
      CALL plot_ini_trans (1.0D0,                                &
                 pl_tran_g, pl_tran_gi, pl_tran_f, pl_tran_fi, &
                 cr_gten, cr_rten, cr_eps)
!                                                                       
      latom = .false. 
      uvw (4) = 1.0 
      DO j = 1, 3 
      xmin (j) = 0.0 
      xmax (j) = 0.0 
      ENDDO 
      xx = 0.0 
!                                                                       
      WRITE (iff, 1000) pl_title 
!                                                                       
      DO i = 1, cr_natoms 
!                                                                       
!     --Select atom if:                                                 
!       type has been selected and                                      
!            all atoms are selected                            or       
!                                                                       
      IF (check_select_status (i, pl_latom (cr_iscat (i) ), cr_prop (i),   &
      pl_sel_prop) ) then                                               
!                                                                       
!     --Check dimensions of plotting space                              
!                                                                       
         IF (pl_dim (1, 1) .le.cr_pos (1, i) .and.            &
             cr_pos (1, i) .le.pl_dim (1, 2) .and.            &
             pl_dim (2, 1) .le.cr_pos (2, i) .and.            &
             cr_pos (2, i) .le.pl_dim (2, 2) .and.            &
             pl_dim (3, 1) .le.cr_pos (3, i) .and.            &
             cr_pos (3, i) .le.pl_dim (3, 2) ) then 
!                                                                       
!     ------write atom position                                         
!                                                                       
            latom = .true. 
            DO j = 1, 3 
            uvw (j) = cr_pos (j, i) 
            ENDDO 
            CALL tran_ca (uvw, pl_tran_f, lscreen) 
            DO j = 1, 3 
            xmin (j) = min (xmin (j), uvw (j) ) 
            xmax (j) = max (xmax (j), uvw (j) ) 
            ENDDO 
            WRITE (iff, 2100) (uvw (j), j = 1, 3) 
         ENDIF 
      ENDIF 
      ENDDO 
      WRITE (iff, * ) 
!                                                                       
      IF (.not.latom) then 
         ier_num = - 58 
         ier_typ = ER_APPL 
      ENDIF 
!                                                                       
 1000 FORMAT    ('frame ',a) 
 2100 FORMAT    (3f11.3) 
!                                                                       
      END SUBROUTINE plot_frames                    
!*****7*****************************************************************
      SUBROUTINE plot_diamond (iff) 
!-                                                                      
!     Writes the selected atoms in a format suitable for plotting       
!     using diamond                                                     
!+                                                                      
      USE discus_config_mod 
      USE crystal_mod 
      USE modify_func_mod
      USE discus_plot_mod 
      USE discus_plot_init_mod
      USE trans_sup_mod
      USE blanks_mod
      USE errlist_mod 
      USE precision_mod
      IMPLICIT none 
!                                                                       
!                                                                       
      CHARACTER(LEN=PREC_STRING) :: zeile 
      INTEGER :: laenge 
      INTEGER :: i, j, iff 
      INTEGER :: natoms
      LOGICAL :: latom, lspace 
      LOGICAL :: lscreen 
REAL(kind=PREC_DP), dimension(4) ::  uvw (4) 
!                                                                       
!                                                                       
      DATA lscreen / .false. / 
      DATA lspace / .true. / 
!                                                                       
      CALL plot_ini_trans (1.0D0,                                &
                 pl_tran_g, pl_tran_gi, pl_tran_f, pl_tran_fi, &
                 cr_gten, cr_rten, cr_eps)
!                                                                       
      latom = .false. 
!
      natoms = 0
      DO i = 1, cr_natoms 
         IF (check_select_status(i, pl_latom(cr_iscat(i)), cr_prop (i),pl_sel_prop) ) THEN
            natoms = natoms + 1
         ENDIF
      ENDDO
!                                                                       
!     Write number of atoms                                             
!                                                                       
      WRITE (zeile, 1000) natoms 
      laenge = 20 
      CALL rem_bl (zeile, laenge) 
      WRITE (iff, 1100) zeile (1:laenge) 
!                                                                       
!     Write title                                                       
!                                                                       
      WRITE (iff, 1100) cr_name 
!                                                                       
      DO i = 1, cr_natoms 
!                                                                       
!     --Select atom if:                                                 
!       type has been selected and                                      
!            all atoms are selected                            or       
!                                                                       
      IF (check_select_status (i, pl_latom (cr_iscat (i) ), cr_prop (i),   &
      pl_sel_prop) ) then                                               
!                                                                       
!     --Check dimensions of plotting space                              
!                                                                       
         IF (pl_dim (1, 1) .le.cr_pos (1, i) .and.       &
             cr_pos (1, i) .le.pl_dim (1, 2) .and.       &
             pl_dim (2, 1) .le.cr_pos (2, i) .and.       &
             cr_pos (2, i) .le.pl_dim (2, 2) .and.       &
             pl_dim (3, 1) .le.cr_pos (3, i) .and.       &
             cr_pos (3, i) .le.pl_dim (3, 2) ) then
!                                                                       
!     ------write atom position                                         
!                                                                       
            latom = .true. 
            DO j = 1, 3 
            uvw (j) = cr_pos (j, i) 
            ENDDO 
            CALL tran_ca (uvw, pl_tran_f, lscreen) 
            WRITE (iff, 2100) cr_at_lis (cr_iscat (i) ),  &
                              (uvw (j), j = 1, 3)
         ENDIF 
      ENDIF 
      ENDDO 
!                                                                       
      IF (.not.latom) then 
         ier_num = - 58 
         ier_typ = ER_APPL 
      ENDIF 
!                                                                       
 1000 FORMAT    (i20) 
 1100 FORMAT    (a) 
 2100 FORMAT    (a4,3f11.3) 
!                                                                       
      END SUBROUTINE plot_diamond                   
!
END MODULE discus_plot_export_mod
