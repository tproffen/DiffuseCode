MODULE chem_multi_mod
!
! Routines to determine multiple distortions for MMC
!
CONTAINS
!
!*****7*****************************************************************
!
SUBROUTINE chem_disp_multi (ianz, cpara, lpara, werte, maxw, MAX_ATOM_ENV_L, lout) 
!+                                                                      
!     Calculates distortions within the crystal                         
!-                                                                      
      USE discus_allocate_appl_mod
      USE discus_config_mod 
      USE crystal_mod 
      USE atom_name 
      USE chem_mod 
use chem_neig_multi_mod
      USE get_iscat_mod
      USE metric_mod
      USE mc_mod 
      USE mmc_mod 
      USE modify_mod
      USE modify_func_mod
      USE errlist_mod 
USE precision_mod
      USE prompt_mod 
USE support_mod
      IMPLICIT none 
!                                                                       
       
!                                                                       
!     INTEGER maxatom 
!                                                                       
!     PARAMETER (maxatom = chem_max_neig) 
!                                                                       
      INTEGER ianz, maxw 
      CHARACTER ( * ) cpara (maxw) 
      INTEGER lpara (maxw) 
      REAL(KIND=PREC_DP) ::  werte (maxw) 
INTEGER, INTENT(IN) :: MAX_ATOM_ENV_L
      LOGICAL lout 
!                                                                       
      CHARACTER(9) at_name_i, at_name_j 
      INTEGER atom (0:MAX_ATOM_ENV_L, MMC_MAX_CENT) 
      LOGICAL  :: tatom (0:MAX_ATOM_ENV_L, MMC_MAX_CENT) 
      INTEGER natom (MMC_MAX_CENT) 
      INTEGER iianz, i, j, k, is, js, ic 
      INTEGER ja, je 
      INTEGER icent, ncent 
      INTEGER bl_anz (0:maxscat, 0:maxscat) 
      REAL(kind=PREC_DP) :: bl_sum (0:maxscat, 0:maxscat) 
      REAL(kind=PREC_DP) :: bl_s2 (0:maxscat, 0:maxscat) 
      REAL(kind=PREC_DP) :: patom (3, 0:MAX_ATOM_ENV_L, MMC_MAX_CENT) 
      REAL(kind=PREC_DP) ::  u (3), v (3), d (3), di 
      LOGICAL lfile 
!                                                                       
!     allocate displacement arrays
!
      CALL alloc_chem_disp(CHEM_MAX_COR, MAXSCAT)
!                                                                       
      iianz = 2 
      CALL get_iscat (iianz, cpara, lpara, werte, maxw, .false.) 
      IF (ier_num.ne.0) return 
!                                                                       
!------ write output                                                    
!                                                                       
      IF (lout) then 
         WRITE (output_io, 1000) cpara (1) (1:lpara (1) ), cpara (2)    &
         (1:lpara (2) )                                                 
         IF (ianz.gt.2) then 
            CALL oeffne (37, cpara (3) , 'unknown')
         ENDIF 
      ENDIF 
!                                                                       
      lfile = lout.and. (ianz.gt.2) 
!                                                                       
!------ loop over all defined correlations                              
!                                                                       
      DO ic = 1, chem_ncor 
!                                                                       
!     --reset the average and sigma values. Necessary when the          
!         correlations or                                               
!       periodicity has been changed                                    
!                                                                       
      DO i = 0, cr_nscat 
      DO j = i, cr_nscat 
      chem_disp_ave (ic, i, j) = 0.0 
      chem_disp_sig (ic, i, j) = 0.0 
      ENDDO 
      ENDDO 
      IF (chem_ctyp (ic) .eq.CHEM_VEC.or.chem_ctyp (ic)                 &
      .eq.CHEM_ENVIR.or.chem_ctyp (ic) .eq.CHEM_RANGE.or.chem_ctyp (ic) &
      .eq.CHEM_DIST) then                                               
!                                                                       
!------ --- reset counters, vectors, ..                                 
!                                                                       
         DO i = 0, cr_nscat 
         DO j = 0, cr_nscat 
         bl_sum (i, j) = 0.0 
         bl_s2 (i, j) = 0.0 
         bl_anz (i, j) = 0 
         ENDDO 
         ENDDO 
!                                                                       
!------ --- calculate distortions                                       
!                                                                       
         DO i = 1, cr_natoms 
         IF (atom_allowed (i, werte, ianz, maxw) ) then 
            CALL chem_neighbour_multi (i, ic, atom, patom, tatom, natom, ncent,&
            MAXW, MAX_ATOM_ENV_L)                                                    
            DO icent = 1, ncent 
            IF (natom (icent) .gt.0) then 
               IF (i.eq.atom (0, icent) ) then 
!                                                                       
!     --------- The selected atom is the central atom, check all atoms  
!                                                                       
                  ja = 1 
                  je = natom (icent) 
                  DO k = 1, 3 
                  u (k) = patom (k, 0, icent) 
                  ENDDO 
                  is = cr_iscat (1,i) 
               ELSE 
!                                                                       
!     --------- The selected atom is a neighbour, check central         
!                 atom only                                             
!                                                                       
                  ja = 0 
                  je = 0 
                  DO k = 1, 3 
                  u (k) = cr_pos (k, i) 
                  ENDDO 
                  is = cr_iscat (1,i) 
               ENDIF 
               DO j = ja, je 
               IF (atom_allowed (atom (j, icent), werte, ianz, maxw) )  &
               then                                                     
                  DO k = 1, 3 
                  v (k) = patom (k, j, icent) 
                  d (k) = v (k) - u (k) 
                  ENDDO 
                  di = do_blen (.true., u, v) 
                  js = cr_iscat (1,atom (j, icent) ) 
                  IF (lfile) WRITE (37, 3000) d, is, js 
                  bl_sum (is, js) = bl_sum (is, js) + di 
                  bl_s2 (is, js) = bl_s2 (is, js) + di**2 
                  bl_anz (is, js) = bl_anz (is, js) + 1 
               ENDIF 
               ENDDO 
            ENDIF 
            ENDDO 
         ENDIF 
         ENDDO 
!                                                                       
!------ - write results and save to res_para block                      
!                                                                       
         DO i = 0, cr_nscat 
         DO j = i, cr_nscat 
         IF (bl_anz (i, j) .ne.0.or.bl_anz (j, i) .ne.0) then 
            chem_disp_ave (ic, i, j) = (bl_sum (i, j) + bl_sum (j, i) ) &
            / (bl_anz (i, j) + bl_anz (j, i) )                          
            chem_disp_sig (ic, i, j) = (bl_s2 (i, j) + bl_s2 (j, i) )   &
            / (bl_anz (i, j) + bl_anz (j, i) )                          
            chem_disp_sig (ic, i, j) = (chem_disp_sig (ic, i, j)        &
            - (chem_disp_ave (ic, i, j) **2) )                          
            IF (chem_disp_sig (ic, i, j) .gt.0) then 
               chem_disp_sig (ic, i, j) = sqrt (chem_disp_sig (ic, i, j)&
               )                                                        
            ELSE 
               chem_disp_sig (ic, i, j) = 0.0 
            ENDIF 
            IF (lout) then 
               at_name_i = at_name (i) 
               at_name_j = at_name (j) 
               WRITE (output_io, 2000) ic, at_name_i, at_name_j,        &
               chem_disp_ave (ic, i, j), chem_disp_sig (ic, i, j),      &
               bl_anz (i, j) + bl_anz (j, i)                            
            ENDIF 
         ENDIF 
         ENDDO 
         ENDDO 
         IF (ier_num.ne.0) return 
      ENDIF 
      ENDDO 
!                                                                       
      IF (lfile) close (37) 
!                                                                       
 1000 FORMAT (  ' Calculating distortions ',/,                          &
     &          '    Atom types : A = ',A4,' and B = ',A4,' ',//,       &
     &          '    Neig.  Atom A      Atom B       distance',         &
     &          '   sigma     # pairs',/,4x,60('-'))                    
 2000 FORMAT (4x,i3,3x,a9,3x,a9,5x,f7.3,3x,f7.3,3x,i8) 
 3000 FORMAT (3(f12.5,1x),3x,2(i3,1x)) 
!                                                                       
      END SUBROUTINE chem_disp_multi                
!
!*****7*****************************************************************
!
      SUBROUTINE chem_angle_multi (ianz, cpara, lpara, werte, uerte,    &
      verte, maxw, MAX_ATOM_ENV_L, lout)                                                
!+                                                                      
!     Calculates angular distortions within the crystal                 
!-                                                                      
      USE discus_allocate_appl_mod
      USE discus_config_mod 
      USE crystal_mod 
      USE atom_name 
      USE chem_mod 
use chem_neig_multi_mod
      USE get_iscat_mod
      USE metric_mod
      USE mc_mod 
      USE mmc_mod 
      USE modify_mod   
      USE modify_func_mod   
      USE errlist_mod 
      USE get_params_mod
USE precision_mod
      USE prompt_mod 
USE support_mod
      IMPLICIT none 
!                                                                       
       
!                                                                       
!     INTEGER maxatom 
!                                                                       
!     PARAMETER (maxatom = chem_max_neig) 
!                                                                       
      INTEGER ianz, maxw 
      CHARACTER ( * ) cpara (maxw) 
      INTEGER lpara (maxw) 
      REAL(KIND=PREC_DP) :: werte (maxw) 
      REAL(KIND=PREC_DP) :: uerte (maxw) 
      REAL(KIND=PREC_DP) :: verte (maxw) 
INTEGER, INTENT(IN) :: MAX_ATOM_ENV_L
      LOGICAL lout 
!                                                                       
      CHARACTER(9) at_name_i, at_name_j 
      CHARACTER(9) name_1, name_2, name_3 
      INTEGER atom (0:MAX_ATOM_ENV_L, MMC_MAX_CENT) 
LOGICAL :: tatom (0:MAX_ATOM_ENV_L, MMC_MAX_CENT) 
      INTEGER natom (MMC_MAX_CENT) 
      INTEGER iianz, i, j, k, is, js, ic 
      INTEGER jj 
      INTEGER jjanz, kkanz 
      INTEGER icent, ncent 
      INTEGER lname_1, lname_2, lname_3 
      INTEGER ba_anz (0:maxscat, 0:maxscat) 
      REAL(kind=PREC_DP) :: ba_sum (0:maxscat, 0:maxscat) 
      REAL(kind=PREC_DP) :: ba_s2 (0:maxscat, 0:maxscat) 
      REAL(kind=PREC_DP) :: patom (3, 0:MAX_ATOM_ENV_L, MMC_MAX_CENT) 
      REAL(kind=PREC_DP) :: u (3), v (3), w (3), wi, wis 
      LOGICAL lfile 
!                                                                       
!     allocate displacement arrays
!
      CALL alloc_chem_disp(CHEM_MAX_COR, MAXSCAT)
!                                                                       
      iianz = 1 
      CALL get_iscat (iianz, cpara, lpara, werte, maxw, .false.) 
      name_1 = cpara (1) 
      lname_1 = lpara (1) 
      CALL del_params (1, ianz, cpara, lpara, maxw) 
      jjanz = 1 
      CALL get_iscat (jjanz, cpara, lpara, uerte, maxw, .false.) 
      name_2 = cpara (1) 
      lname_2 = lpara (1) 
      CALL del_params (1, ianz, cpara, lpara, maxw) 
      kkanz = 1 
      CALL get_iscat (kkanz, cpara, lpara, verte, maxw, .false.) 
      name_3 = cpara (1) 
      lname_3 = lpara (1) 
      IF (ier_num.ne.0) return 
!                                                                       
!------ write output                                                    
!                                                                       
      IF (lout) then 
         WRITE (output_io, 1000) name_1 (1:lname_1), name_2 (1:lname_2),&
         name_3 (1:lname_3)                                             
         IF (ianz.gt.1) then 
            CALL oeffne (37, cpara (2) , 'unknown') 
         ENDIF 
      ENDIF 
!                                                                       
      lfile = lout.and. (ianz.gt.1) 
!                                                                       
!------ loop over all defined correlations                              
!                                                                       
      DO ic = 1, chem_ncor 
      IF (chem_ctyp (ic) .eq.CHEM_ANG.or.chem_ctyp (ic) .eq.CHEM_ENVIR) &
      then                                                              
!                                                                       
!------ - reset counters, vectors, ..                                   
!                                                                       
         DO i = 0, cr_nscat 
         DO j = 0, cr_nscat 
         ba_sum (i, j) = 0.0 
         ba_s2 (i, j) = 0.0 
         ba_anz (i, j) = 0 
         ENDDO 
         ENDDO 
!                                                                       
!------ - calculate angles                                              
!                                                                       
         DO i = 1, cr_natoms 
         IF (atom_allowed (i, werte, iianz, maxw) ) then 
            CALL chem_neighbour_multi (i, ic, atom, patom, tatom, natom, ncent,&
            MAXW, MAX_ATOM_ENV_L)                                                    
            DO icent = 1, ncent 
            IF (natom (icent) .gt.1) then 
               DO k = 1, 3 
               u (k) = patom (k, 0, icent) 
               ENDDO 
               DO j = 1, natom (icent) - 1 
               DO jj = j + 1, natom (icent) 
               IF (atom_allowed (atom (j, icent), uerte, jjanz, maxw)   &
               .and.atom_allowed (atom (jj, icent), verte, kkanz, maxw) &
               .or.atom_allowed (atom (jj, icent), uerte, jjanz, maxw)  &
               .and.atom_allowed (atom (j, icent), verte, kkanz, maxw) )&
               then                                                     
                  DO k = 1, 3 
                  v (k) = patom (k, j, icent) 
                  ENDDO 
                  DO k = 1, 3 
                  w (k) = patom (k, jj, icent) 
                  ENDDO 
                  wi = do_bang (.true., v, u, w) 
                  is = cr_iscat (1,atom (j, icent)) 
                  js = cr_iscat (1,atom (jj, icent)) 
                  wis = mmc_target_corr (ic, MC_ANGLE, is, js) 
!DBG                                                                    
!DBG      if(cr_iscat(i,1).eq.11 .or.                                     
!DBG     &              (is.eq.12 .and. js.eq.28) ) then                
!DBG        write(*,*) '==============================================' 
!DBG        write(*,*) 'Selected   , Type',i,cr_iscat(i,1)                
!DBG        write(*,*) 'Position zentral ',u                            
!DBG        write(*,*) 'Zentral Atom     ',atom(0 ,icent)               
!DBG        write(*,*) 'Correlation      ',ic                           
!DBG        write(*,*) 'Nachbarn         ',natom(icent)                 
!DBG        write(*,*) 'Erster           ',atom(j ,icent)               
!DBG        write(*,*) 'Position         ',v                            
!DBG        write(*,*) 'Zweiter          ',atom(jj,icent)               
!DBG        write(*,*) 'Position         ',w                            
!DBG        write(*,*) 'Winkel           ',wi                           
!DBG        write(*,*) 'Sollwinkel       ',wis                          
!DBG        write(*,*) 'is,js            ',is,js                        
!DBG      endif                                                         
                  IF (wis.le.90.) then 
                     IF (wi.gt.1.5 * wis) then 
                        wi = mod (wi + wis / 2., wis) - wis / 2. 
                     ENDIF 
                  ENDIF 
                  IF (lfile) WRITE (37, 3000) u, is, js 
                  ba_sum (is, js) = ba_sum (is, js) + wi 
                  ba_s2 (is, js) = ba_s2 (is, js) + wi**2 
                  ba_anz (is, js) = ba_anz (is, js) + 1 
!DBG      if(cr_iscat(i,1).eq.11 .or.                                     
!DBG     &              (is.eq.12 .and. js.eq.28) ) then                
!DBG        write(*,*) 'i,cr_iscat(i,1),is,js,ba_anz(is,js)   ',          
!DBG      i,cr_iscat(i,1),is,js,ba_anz(is,js)                             
!DBG        write(*,*) 'ba_sum(is,js),wi   ',ba_sum(is,js),wi           
!DBG      endif                                                         
               ENDIF 
               ENDDO 
               ENDDO 
            ENDIF 
            ENDDO 
         ENDIF 
         ENDDO 
!                                                                       
!------ - write results and save to res_para block                      
!                                                                       
         DO i = 0, cr_nscat 
         DO j = i, cr_nscat 
         IF (ba_anz (i, j) .ne.0.or.ba_anz (j, i) .ne.0) then 
            chem_disp_ave (ic, i, j) = (ba_sum (i, j) + ba_sum (j, i) ) &
            / (ba_anz (i, j) + ba_anz (j, i) )                          
            chem_disp_sig (ic, i, j) = (ba_s2 (i, j) + ba_s2 (j, i) )   &
            / (ba_anz (i, j) + ba_anz (j, i) )                          
            chem_disp_sig (ic, i, j) = (chem_disp_sig (ic, i, j)        &
            - (chem_disp_ave (ic, i, j) **2) )                          
            IF (chem_disp_sig (ic, i, j) .gt.0) then 
               chem_disp_sig (ic, i, j) = sqrt (chem_disp_sig (ic, i, j)&
               )                                                        
            ELSE 
               chem_disp_sig (ic, i, j) = 0.0 
            ENDIF 
            IF (lout) then 
               at_name_i = at_name (i) 
               at_name_j = at_name (j) 
               WRITE (output_io, 2000) ic, at_name_i, at_name_j,        &
               chem_disp_ave (ic, i, j), chem_disp_sig (ic, i, j),      &
               ba_anz (i, j) + ba_anz (j, i)                            
            ENDIF 
         ENDIF 
         ENDDO 
         ENDDO 
         IF (ier_num.ne.0) return 
      ENDIF 
      ENDDO 
!                                                                       
      IF (lfile) close (37) 
!                                                                       
 1000 FORMAT ( ' Calculating angles ',/,                                &
     &         '    Zentral atom   = ',A4,/,                            &
     &         '    Atom types : A = ',A4,' and B = ',A4,' ',//,        &
     &         '    Neig.  Atom A      Atom B       angle   ',          &
     &         '   sigma     # pairs',/,4x,60('-'))                     
 2000 FORMAT (4x,i3,3x,a9,3x,a9,5x,f7.3,3x,f7.3,3x,i8) 
 3000 FORMAT (3(f12.5,1x),3x,2(i3,1x)) 
!                                                                       
      END SUBROUTINE chem_angle_multi               
!
!*****7*****************************************************************
!
END MODULE chem_multi_mod
