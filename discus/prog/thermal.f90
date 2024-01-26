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
!                                                                       
      USE errlist_mod 
      USE get_params_mod
USE lib_random_func
      USE prompt_mod 
      USE random_mod
USE precision_mod
use wink_mod
!
      IMPLICIT none 
       
!                                                                       
      INTEGER, PARAMETER :: maxw = 2 
!                                                                       
      CHARACTER (LEN=*), INTENT(INOUT) :: line 
      INTEGER          , INTENT(INOUT) :: ll
!
      CHARACTER(LEN=PREC_STRING) :: cpara (maxw) 
      CHARACTER(LEN=9)    :: at_name_i 
      INTEGER lpara (maxw) 
      INTEGER i, j, k, is, ii, ianz 
      INTEGER uc_n (0:maxscat) 
      REAL(kind=PREC_DP) :: uc (3), up (3) , vv(3)
      REAL(kind=PREC_DP) :: uc_max (3, 0:maxscat) 
      REAL(kind=PREC_DP) :: uc_su2 (3, 0:maxscat) 
      REAL(kind=PREC_DP) :: bfac, a 
      LOGICAL flag_all, flag_mol 
real(kind=prec_dp) :: sigma
!
!prin = 0.0D0
!prin(1,1) = 1.0D0     ! Simple a-axis
!prin(2,2) = 1.0D0     ! Simple b-axis
!prin(3,3) = 1.0D0     ! Simple c-axis
!prin(4,1) =     (1.30/8.0/PI**2)
!prin(4,2) =     (0.90/8.0/PI**2)
!prin(4,3) =     (0.90/8.0/PI**2)
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
      ELSEIF (ianz >= 1) THEN 
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
uc_n   = 0
uc_max = 0.0D0 
uc_su2 = 0.0D0 
!                                                                       
!     pi2  = (4.0D0 * atan (1.0D0) ) **2 
bfac = 1.0D0 / (8.0D0 * pi**2) 
!                                                                       
!------ Displace molecules                                              
!                                                                       
cond_mol:      IF (flag_mol) THEN 
         DO i = 1, mole_num_mole 
            is = mole_cont (mole_off (i) + 1) 
!           a = sqrt (bfac * cr_dw (cr_iscat (1,is) ) ) 
            a = sqrt (bfac * mole_biso (mole_type (i) ) ) 
!           CALL ther_vec(flag_all, a, uc, up)
            is = mole_type (i) 
         DO ii = 1, 3 
         IF (flag_all.or.cr_icc (ii) .ne.1) THEN 
            up (ii) = gasdev (DBLE(a)) 
         ELSE 
            up (ii) = 0.0 
         ENDIF 
         ENDDO 
!           CALL trans (up, cr_gmat, uc, 3) 
            uc = matmul(cr_gmat, up)
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
ELSE  cond_mol
!do i=1, cr_ncatoms
!  write(*,*)
!  write(*,'(a,3f9.5,i8, i4)' ) ' Atom ', cr_pos(:,i), i, cr_iscat(1,i)
!  write(*,'(a,3f9.5,a, f9.5)') ' Eigen',cr_prin(1:3,1,i), ' val ', cr_prin(4,1,i) 
!  write(*,'(a,3f9.5,a, f9.5)') ' Eigen',cr_prin(1:3,2,i), ' val ', cr_prin(4,2,i) 
!  write(*,'(a,3f9.5,a, f9.5)') ' Eigen',cr_prin(1:3,3,i), ' val ', cr_prin(4,3,i) 
!enddo
!write(*,*) ' THERMAL ', cr_nanis
!do i=1, cr_nanis
!write(*,'(a,3f9.5,a,f9.5,a,3f9.5)') 'E(1,0,0) => ', cr_prin(1:3,1,i),  ' VAL ', cr_prin(4,1,i), ' HEX ',matmul(cr_emat,cr_prin(1:3,1,i))
!write(*,'(a,3f9.5,a,f9.5,a,3f9.5)') 'E(0,1,0) => ', cr_prin(1:3,2,i),  ' VAL ', cr_prin(4,2,i), ' HEX ',matmul(cr_emat,cr_prin(1:3,2,i))
!write(*,'(a,3f9.5,a,f9.5,a,3f9.5)') 'E(0,0,1) => ', cr_prin(1:3,3,i),  ' VAL ', cr_prin(4,3,i), ' HEX ',matmul(cr_emat,cr_prin(1:3,3,i))
!write(*,*)
!uc = matmul(cr_gmat, up)
!enddo
!do i=1, cr_ncatoms
!write(*,*) ' ATOM ', i, cr_iscat(:,i)
!enddo
   loop_atom_disp:DO i = 1, cr_natoms 
      is = cr_iscat (1,i) 
      a  = sqrt (bfac * cr_dw (is) ) 
!           CALL ther_vec(flag_all, a, uc, up)
      cond_anis:if(is>=0) then    ! FOR GE 
         ii = cr_ianis(i)
         ii = cr_iscat(3,i)
         call ther_anis(cr_prin(:,:,ii),uc, up, i)
!if(i<=cr_ncatoms) then 
!   write(*,'(a,i4, i4, i4    )') 'THRMAL ANISO', i, is, ii
!   write(*,'(a,4f8.4,2x, f8.4)') 'PRIN 1    ', cr_prin(:,1,ii), sqrt(cr_prin(4,1,ii))
!   write(*,'(a,4f8.4,2x, f8.4)') 'PRIN 2    ', cr_prin(:,2,ii), sqrt(cr_prin(4,2,ii))
!   write(*,'(a,4f8.4,2x, f8.4)') 'PRIN 3    ', cr_prin(:,3,ii), sqrt(cr_prin(4,3,ii))
!   write(*,'(a,3f8.4,2x, f8.4)') 'VECTOR uc ', uc
!   write(*,'(a,3f8.4,2x, f8.4)') 'VECTOR up ', up
!vv(1) = 1.0D0
!vv(2) = 0.0D0
!vv(3) = 0.0D0
!   write(*,'(a,2(3f8.4,2x, f8.4))') '100       ', vv, matmul(cr_emat, vv)
!vv(1) = 0.0D0
!vv(2) = 1.0D0
!vv(3) = 0.0D0
!   write(*,'(a,2(3f8.4,2x, f8.4))') '010       ', vv, matmul(cr_emat, vv)
!vv(1) = 0.0D0
!vv(2) = 0.0D0
!vv(3) = 1.0D0
!   write(*,'(a,2(3f8.4,2x, f8.4))') '001       ', vv, matmul(cr_emat, vv)
!endif
      else cond_anis
         DO ii = 1, 3 
sigma = a
!if(cr_iscat(1,i)==1) then
!  if(ii==1) sigma = a*1.15
!  if(ii==2) sigma = a*1.00
!  if(ii==3) sigma = a*1.00
!endif
         IF (flag_all.or.cr_icc (ii) .ne.1) THEN 
!           up (ii) = gasdev (DBLE(a)) 
            up (ii) = gasdev (DBLE(sigma)) 
         ELSE 
            up (ii) = 0.0 
         ENDIF 
         ENDDO 
!        CALL trans (up, cr_gmat, uc, 3) 
         uc = matmul(cr_emat, up)
!if(i==5) then
!write(*,'(a,1f8.4,2x, f8.4, i8)') 'THRMAL ANISO', cr_dw(is), a, is
!write(*,'(a,3f8.4,2x, f8.4)') 'VECTOR up ', up
!write(*,'(a,3f8.4,2x, f8.4)') 'VECTOR uc ', uc
!endif
      endif cond_anis
!write(*,*) ' ATOM ', i, is, cr_iscat(1,i), uc, is<=9
      DO j = 1, 3 
         cr_pos (j, i)  = cr_pos (j, i) + uc (j) 
         uc_max (j, is) = max (uc_max (j, is), abs (up (j) ) ) 
         uc_su2 (j, is) = uc_su2 (j, is) + up (j) **2 
      ENDDO 
      uc_n (is) = uc_n (is) + 1 
   ENDDO loop_atom_disp 
ENDIF  cond_mol
!                                                                       
!     Update crystal dimensions                                         
!                                                                       
      CALL update_cr_dim 
!                                                                       
!------ print statistics of 'therm' command                             
!                                                                       
cond_pr_mol: IF (flag_mol) THEN 
   WRITE (output_io, 1000) 'Molecule' 
   DO i = 1, mole_num_type 
      IF (uc_n (i) .ne.0) THEN 
         WRITE(output_io, 1050) i, cr_dw(i), (bfac * cr_dw (i) ),  &
            ((uc_su2(j, i) / uc_n(i)), j = 1, 3), (uc_max(j, i), j = 1, 3)                                                   
      ENDIF 
   ENDDO 
ELSE cond_pr_mol
!  write(*,*) ' uc_n   ', uc_n(1:2)
!  write(*,*) ' uc_max 1 ', uc_max(1:3, 1), uc_su2(1:3,1)
!  write(*,*) ' uc_max 2 ', uc_max(1:3, 2), uc_su2(1:3,2)
   WRITE (output_io, 1000) 'Atom    ' 
   DO i = 0, cr_nscat 
      IF (uc_n (i) .ne.0) THEN 
         at_name_i = at_name (i) 
         WRITE (output_io, 1150) at_name_i, cr_dw (i), (bfac * cr_dw &
            (i) ), ( (uc_su2 (j, i) / uc_n (i) ), j = 1, 3), (uc_max (j,&
            i), j = 1, 3)                                               
      ENDIF 
   ENDDO 
ENDIF cond_pr_mol  
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
USE lib_random_func
USE random_mod
USE precision_mod
!
IMPLICIT NONE
!
LOGICAL           , INTENT(IN ) :: flag_all
REAL(kind=PREC_DP), INTENT(IN ) :: a
REAL(kind=PREC_DP), DIMENSION(3), INTENT(OUT) :: uc
REAL(kind=PREC_DP), DIMENSION(3), INTENT(OUT) :: up
!
INTEGER            :: i
REAL(kind=PREC_DP) :: disp
REAL(kind=PREC_DP) :: length
REAL(kind=PREC_DP) :: up_length
!
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
!CALL trans (up, cr_gmat, uc, 3) 
uc = matmul(cr_gmat, up)
length = SQRT(skalpro(uc,uc,cr_gten))
disp   = gasdev(DBLE(a))
uc(:)  = uc(:) /length * disp
up(:)  = up(:) /up_length * disp
!
END SUBROUTINE ther_vec
!
!*******************************************************************************
!
subroutine ther_anis(prin, uc, up, jj)
!-
!  Calculate an anisotropic displacement
!+
!
use crystal_mod
use metric_mod
use precision_mod
!use random_mod
use lib_random_func
use wink_mod
!
real(kind=PREC_DP), dimension(4, 3), intent(in)  :: prin  ! prin(i,1:3) vectors, prin(i,4) <u>**2
real(kind=PREC_DP), dimension(3)   , intent(out) :: uc    ! Crystal space
real(kind=PREC_DP), dimension(3)   , intent(out) :: up    ! "plot" space == cartesian
integer :: jj
!
integer                            :: i     ! Loop index
real(kind=PREC_DP), dimension(3,3) :: ud    ! Displacement vector along prinipal axis
real(kind=prec_dp) :: sigma                 ! Sigma for Gaussian
real(kind=prec_dp) :: disp                  ! actual displacement
!
ud = 0.0D0
up = 0.0D0
!if(jj<=2) then
!  write(*,*) ' ATOM ', jj
!  write(*,'(a,3f9.5,a, f9.5)') ' Eigen',prin(1,1:3), ' val ', prin(1,4) 
!  write(*,'(a,3f9.5,a, f9.5)') ' Eigen',prin(2,1:3), ' val ', prin(2,4) 
!  write(*,'(a,3f9.5,a, f9.5)') ' Eigen',prin(3,1:3), ' val ', prin(3,4) 
!endif
!
do i=1, 3
   ud(:,i) = prin(1:3,i)    ! Copy principal vector
!  length = sqrt(dot_product(u, matmul(gten, u)))
!  ud(i,:) = ud(i,:) / length         ! Normalize to unit length
!
   sigma = sqrt(prin(4,i))
   disp  = gasdev(sigma)
   ud(:,i) = ud(:,i) * disp
!if(jj<=2) then
!write(*,'(a,4f8.4,2x, f8.4)') ' DISP ', prin(4,i), sigma, disp
!write(*,'(a,4f8.4,2x, f8.4)') ' vect ', ud(:,i)
!endif
   up = up + ud(:,i)
enddo
!
uc = matmul(cr_emat, up)
!if(jj<=2) then
!write(*,'(a,4f8.4,2x, f8.4)') ' UP   ', up(:)
!write(*,'(a,4f8.4,2x, f8.4)') ' UC   ', uc(:)
!endif
!
end subroutine ther_anis
!
!*******************************************************************************
!
END MODULE thermal_mod
