MODULE thermal_mod
!
CONTAINS
!*****7*****************************************************************
!                                                                       
SUBROUTINE ther_displ(line, ll) 
!-                                                                      
!  Displaces atoms in a completely random fashion. The mean    
!  square displacement is given by the temperature factor.     
!  Molecules are displaced according to the B value of the           
!  atom at the origin.                                             
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
INTEGER, PARAMETER :: MAXW = 2 
!                                                                       
CHARACTER(LEN=*), INTENT(INOUT) :: line 
INTEGER         , INTENT(INOUT) :: ll
!
CHARACTER(LEN=PREC_STRING), dimension(MAXW) :: cpara
CHARACTER(LEN=9)                            :: at_name_i 
INTEGER, dimension(MAXW)                    :: lpara
INTEGER :: i, j, k, is, ii, ianz 
INTEGER, dimension(0:MAXSCAT) ::  uc_n !(0:maxscat) 
REAL(kind=PREC_DP), dimension(3) :: uc   ! Displacement vector in crystal space
REAL(kind=PREC_DP), dimension(3) :: up   ! Displacement vector in "plot" == cartesian space
REAL(kind=PREC_DP), dimension(3, 0:MAXSCAT) :: uc_max ! Maximum displacements for statistics
REAL(kind=PREC_DP), dimension(3, 0:MAXSCAT) :: uc_su2 ! displacement sigmas for statistics
REAL(kind=PREC_DP) :: bfac       ! Conversio B == u
real(kind=PREC_DP) :: a          ! Sigma for Gaussian distribution
LOGICAL            :: flag_all   ! 'therm 3d' if true
LOGICAL            :: flag_mol   ! Work on molecules if true
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
!------ reset counters                                                  
!                                                                       
uc_n   = 0
uc_max = 0.0D0 
uc_su2 = 0.0D0 
!                                                                       
bfac = 1.0D0 / (8.0D0 * pi**2) 
!                                                                       
!------ Displace molecules                                              
!                                                                       
cond_mol:      IF (flag_mol) THEN    ! Molecuse still isotropic
   DO i = 1, mole_num_mole 
      is = mole_cont(mole_off(i) + 1) 
!
      a = sqrt(bfac * mole_biso(mole_type(i) ) ) 
!     CALL ther_vec(flag_all, a, uc, up)
      is = mole_type (i) 
      DO ii = 1, 3 
         IF(flag_all.or.cr_icc (ii) .ne.1) THEN 
            up(ii) = gasdev(DBLE(a)) 
         ELSE 
            up(ii) = 0.0 
         ENDIF 
      ENDDO 
!
      uc = matmul(cr_gmat, up)
      DO j = 1, 3 
         DO k = 1, mole_len(i) 
            ii = mole_cont(mole_off(i) + k) 
            cr_pos(j, ii) = cr_pos(j, ii) + uc(j) 
         ENDDO 
         uc_su2(j, is) = uc_su2 (j, is) + up(j) **2 
         uc_max(j, is) = max(uc_max(j, is), abs(up (j) ) ) 
      ENDDO 
      uc_n (is) = uc_n (is) + 1 
   ENDDO 
!                                                                       
!------ displace atoms                                                  
!                                                                       
ELSE  cond_mol
   loop_atom_disp:DO i = 1, cr_natoms 
      is = cr_iscat(1,i) 
      ii = cr_iscat(3,i)
      call ther_anis(cr_prin(:,:,ii),uc, up, cr_emat, cr_icc, flag_all)
!OLD      a  = sqrt (bfac * cr_dw (is) ) 
!OLD!           CALL ther_vec(flag_all, a, uc, up)
!OLD         DO ii = 1, 3 
!OLDsigma = a
!OLD         IF (flag_all.or.cr_icc (ii) .ne.1) THEN 
!OLD!           up (ii) = gasdev (DBLE(a)) 
!OLD            up (ii) = gasdev (DBLE(sigma)) 
!OLD         ELSE 
!OLD            up (ii) = 0.0 
!OLD         ENDIF 
!OLD         ENDDO 
!OLD         uc = matmul(cr_emat, up)
!OLD      endif cond_anis
!write(*,*) ' ATOM ', i, is, cr_iscat(1,i), uc, is<=9
      DO j = 1, 3 
         cr_pos (j, i)  = cr_pos (j, i) + uc (j) 
         uc_max (j, is) = max (uc_max (j, is), abs (up (j) ) ) 
         uc_su2 (j, is) = uc_su2 (j, is) + up (j) **2 
      ENDDO 
      uc_n(is) = uc_n(is) + 1 
   ENDDO loop_atom_disp 
ENDIF  cond_mol
!                                                                       
!     Update crystal dimensions                                         
!                                                                       
CALL update_cr_dim 
!                                                                       
!------ print statistics of 'therm' command                             
!                                                                       
cond_pr_mol: IF(flag_mol) THEN 
   WRITE (output_io, 1000) 'Molecule' 
   DO i = 1, mole_num_type 
      IF (uc_n(i) .ne.0) THEN 
         WRITE(output_io, 1050) i, cr_dw(i), (bfac * cr_dw (i) ),  &
            ((uc_su2(j, i) / uc_n(i)), j = 1, 3), (uc_max(j, i), j = 1, 3)                                                   
      ENDIF 
   ENDDO 
ELSE cond_pr_mol
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
subroutine ther_anis(prin, uc, up, emat, icc, flag_all)
!-
!  Calculate an anisotropic displacement
!+
!
!use crystal_mod
use metric_mod
use precision_mod
!use random_mod
use lib_random_func
use wink_mod
!
real(kind=PREC_DP), dimension(4, 3), intent(in)  :: prin  ! prin(i,1:3) vectors, prin(i,4) <u>**2
real(kind=PREC_DP), dimension(3)   , intent(out) :: uc    ! Crystal space
real(kind=PREC_DP), dimension(3)   , intent(out) :: up    ! "plot" space == cartesian
real(kind=PREC_DP), dimension(3, 3), intent(in)  :: emat  ! Conversion matrix crystal to cartesian bases
integer           , dimension(3)   , intent(in)  :: icc   !  number of unit cells (for therm 2D)
logical                            , intent(in)  :: flag_all ! 'therm 3D' ==true
!integer :: jj
!
integer                            :: i     ! Loop index
real(kind=PREC_DP), dimension(3,3) :: ud    ! Displacement vector along prinipal axis
real(kind=prec_dp) :: sigma                 ! Sigma for Gaussian
real(kind=prec_dp) :: disp                  ! actual displacement
!
ud = 0.0D0
up = 0.0D0
!
do i=1, 3
   ud(:,i) = prin(1:3,i)    ! Copy principal vector
!
   sigma = sqrt(prin(4,i))
   disp  = gasdev(sigma)
   ud(:,i) = ud(:,i) * disp   ! Accumulate displacement along this vector
   up = up + ud(:,i)   ! Accumulate displacement along this vector
enddo
!
!  Restrict dimension if 'therm 2d'
!
do i=1,3
   if(flag_all.or.icc(i) /= 1) then 
      continue
   else
      up(i) = 0.0_PREC_DP
   endif
enddo
!
uc = matmul(emat, up)
!
end subroutine ther_anis
!
!*******************************************************************************
!
END MODULE thermal_mod
