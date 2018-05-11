MODULE thermal_mod
!
CONTAINS
!*****7*****************************************************************
!                                                                       
      SUBROUTINE ther_displ (line, ll) 
!-                                                                      
!           Displaces atoms in a completely random fashion. The mean    
!           square displacement is given by the temperature factor.     
!     Molecules are displaced according to the B value of the           
!       atom at the origin.                                             
!+                                                                      
      USE discus_config_mod 
      USE crystal_mod 
      USE atom_name
      USE molecule_mod 
      USE update_cr_dim_mod
      USE trafo_mod
!                                                                       
      USE errlist_mod 
      USE get_params_mod
      USE prompt_mod 
!
      IMPLICIT none 
       
!                                                                       
      INTEGER, PARAMETER :: maxw = 2 
!                                                                       
      CHARACTER (LEN=*), INTENT(INOUT) :: line 
      INTEGER          , INTENT(INOUT) :: ll
!
      CHARACTER(LEN=1024) :: cpara (maxw) 
      CHARACTER(LEN=9)    :: at_name_i 
      INTEGER lpara (maxw) 
      INTEGER i, j, k, is, ii, ianz 
      INTEGER uc_n (0:maxscat) 
      REAL uc (3), up (3) 
      REAL uc_max (3, 0:maxscat) 
      REAL uc_su2 (3, 0:maxscat) 
      REAL pi2, bfac, a 
      LOGICAL flag_all, flag_mol 
!
      flag_all=.true.
      flag_mol=.false.
!                                                                       
!------ get parameters                                                  
!                                                                       
      CALL get_params (line, ianz, cpara, lpara, maxw, ll) 
      IF (ier_num /= 0) RETURN 
!                                                                       
      IF (ianz.eq.0) THEN 
         flag_all = .true. 
         flag_mol = .false. 
      ELSEIF (ianz == 1) THEN 
         flag_mol = (cpara (1) (1:1) .eq.'m') 
         flag_all = .true. 
         IF (flag_mol) THEN 
            IF (ianz == 2) THEN 
               flag_all = (cpara (2) (1:1) .ne.'2') 
            ENDIF 
         ELSE 
            flag_all = (cpara (1) (1:1) .ne.'2') 
         ENDIF 
      ENDIF 
!                                                                       
!------ reset counters                                                  
!                                                                       
!     DO i = 0, maxscat 
!        uc_n (i) = 0 
!        DO j = 1, 3 
!           uc_max (j, i) = 0.0 
!           uc_su2 (j, i) = 0.0 
!        ENDDO 
!     ENDDO 
      uc_n(:)   = 0
      uc_max(:, :) = 0.0 
      uc_su2(:, :) = 0.0 
!                                                                       
      pi2  = (4.0 * atan (1.0) ) **2 
      bfac = 1.0 / (8.0 * pi2) 
!                                                                       
!------ Displace molecules                                              
!                                                                       
      IF (flag_mol) THEN 
         DO i = 1, mole_num_mole 
            is = mole_cont (mole_off (i) + 1) 
!           a = sqrt (bfac * cr_dw (cr_iscat (is) ) ) 
            a = sqrt (bfac * mole_biso (mole_type (i) ) ) 
            CALL ther_vec(flag_all, a, uc, up)
            is = mole_type (i) 
!        DO ii = 1, 3 
!        IF (flag_all.or.cr_icc (ii) .ne.1) THEN 
!           up (ii) = gasdev (a) 
!        ELSE 
!           up (ii) = 0.0 
!        ENDIF 
!        ENDDO 
            CALL trans (up, cr_gmat, uc, 3) 
            DO j = 1, 3 
               DO k = 1, mole_len (i) 
                  ii = mole_cont (mole_off (i) + k) 
                  cr_pos (j, ii) = cr_pos (j, ii) + uc (j) 
               ENDDO 
               uc_su2 (j, is) = uc_su2 (j, is) + up (j) **2 
               uc_max (j, is) = max (uc_max (j, is), abs (up (j) ) ) 
            ENDDO 
            uc_n (is) = uc_n (is) + 1 
         ENDDO 
!                                                                       
!------ displace atoms                                                  
!                                                                       
      ELSE 
         DO i = 1, cr_natoms 
            is = cr_iscat (i) 
            a  = sqrt (bfac * cr_dw (is) ) 
            CALL ther_vec(flag_all, a, uc, up)
!        DO ii = 1, 3 
!        IF (flag_all.or.cr_icc (ii) .ne.1) THEN 
!           up (ii) = gasdev (a) 
!        ELSE 
!           up (ii) = 0.0 
!        ENDIF 
!        ENDDO 
!        CALL trans (up, cr_gmat, uc, 3) 
            DO j = 1, 3 
               cr_pos (j, i)  = cr_pos (j, i) + uc (j) 
               uc_max (j, is) = max (uc_max (j, is), abs (up (j) ) ) 
               uc_su2 (j, is) = uc_su2 (j, is) + up (j) **2 
            ENDDO 
            uc_n (is) = uc_n (is) + 1 
         ENDDO 
      ENDIF 
!                                                                       
!     Update crystal dimensions                                         
!                                                                       
      CALL update_cr_dim 
!                                                                       
!------ print statistics of 'therm' command                             
!                                                                       
      IF (flag_mol) THEN 
         WRITE (output_io, 1000) 'Molecule' 
         DO i = 1, mole_num_type 
            IF (uc_n (i) .ne.0) THEN 
            WRITE (output_io, 1050) i, cr_dw (i), (bfac * cr_dw (i) ),  &
            ( (uc_su2 (j, i) / uc_n (i) ), j = 1, 3), (uc_max (j, i),   &
            j = 1, 3)                                                   
            ENDIF 
         ENDDO 
      ELSE 
         WRITE (output_io, 1000) 'Atom    ' 
         DO i = 0, cr_nscat 
            IF (uc_n (i) .ne.0) THEN 
               at_name_i = at_name (i) 
            WRITE (output_io, 1150) at_name_i, cr_dw (i), (bfac * cr_dw &
            (i) ), ( (uc_su2 (j, i) / uc_n (i) ), j = 1, 3), (uc_max (j,&
            i), j = 1, 3)                                               
            ENDIF 
         ENDDO 
      ENDIF 
!                                                                       
 1000 FORMAT (' Thermal displacements summary :',//,                    &
     &        17x,'Input',15x,'Achieved',13x,'Maximum displacement',/   &
     &        3x,a8,4x,'B',3x,'<u**2>',5x,'<ux**2> <uy**2> <uz**2>',    &
     &        8x,'x',6x,'y',6x,'z',/,3x,75('-'))                        
 1050 FORMAT (3x,i5,4x,f5.2,2x,f6.4,6x,3(f6.4,2x),3x,3(f6.4,1x)) 
 1150 FORMAT (3x,a9   ,f5.2,2x,f6.4,6x,3(f6.4,2x),3x,3(f6.4,1x)) 
!
      END SUBROUTINE ther_displ                     
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
SUBROUTINE ther_vec(flag_all, a, uc, up)
!
USE crystal_mod
USE metric_mod
USE trafo_mod
USE random_mod
!
IMPLICIT NONE
!
LOGICAL           , INTENT(IN ) :: flag_all
REAL              , INTENT(IN ) :: a
REAL, DIMENSION(3), INTENT(OUT) :: uc
REAL, DIMENSION(3), INTENT(OUT) :: up
!
INTEGER            :: i
REAL               :: disp
REAL               :: length
REAL               :: up_length
!
REAL :: ran1
REAL :: gasdev
!
search: DO
   DO i=1,3
      IF (flag_all.or.cr_icc (i) .ne.1) THEN 
         up(i) = 2.0*ran1(idum) - 1.0
      ELSE
         up(i) = 0
      ENDIF
   ENDDO
   up_length = up(1)**2 + up(2)**2 + up(3)**2
   IF(up_length <= 1.0) EXIT search
ENDDO search
up_length = SQRT(up_length)
CALL trans (up, cr_gmat, uc, 3) 
length = SQRT(skalpro(uc,uc,cr_gten))
disp   = gasdev(a)
uc(:)  = uc(:) /length * disp
up(:)  = up(:) /up_length * disp
!
END SUBROUTINE ther_vec
END MODULE thermal_mod
