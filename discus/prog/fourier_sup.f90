MODULE fourier_sup
!
USE errlist_mod 
USE precision_mod
!
COMPLEX(KIND=PREC_DP), DIMENSION(:,:), ALLOCATABLE :: fft_field
COMPLEX(KIND=PREC_DP), DIMENSION(:,:), ALLOCATABLE :: fft_sum  
INTEGER, DIMENSION(2) :: fft_dim      ! Size of fft_field
INTEGER               :: fft_grid     ! No of grid points in a unit cell
!
CONTAINS
!**********************************************************************
!**********************************************************************
!
SUBROUTINE four_run 
!+                                                                      
!     Main routine to calculate the Fourier transform for the           
!     given plane in reciprocal space. Based on the program             
!     'diffuse' written by B.D. Butler.                                 
!-                                                                      
USE discus_config_mod 
USE crystal_mod 
USE diffuse_mod 
USE fourier_conv_mod
USE four_strucf_mod
USE fourier_lmn_mod
!                                                                       
use param_mod
USE prompt_mod 
USE precision_mod 
!USE times_mod
USE support_mod
!
IMPLICIT none 
!                                                                       
REAL(KIND=PREC_DP) :: ss!, seknds
REAL(KIND=PREC_DP) :: dnorm
INTEGER :: lbeg (3), csize (3) 
INTEGER :: iscat, nlot, ncell, i! ,j,k
INTEGER :: ii          ! Dummy variable
logical :: four_is_new  ! The reciprocal space dimensions have changed
!                                                                       
!ier_num = 0    ! STRANGE BUG on MacAir with M1 chip ??? 2022-May-11
!                                                                       
!------ preset some values                                              
!                                                                       
ss = seknds (0.0) 
!
!                                                                       
!------ preset some values                                              
!                                                                       
call four_layer (four_is_new)
!                                                                       
!------ zero some arrays
!                                                
!
csf (1:num(1),1:num(2),1:num(3)) = cmplx (0.0D0, 0.0D0, KIND=KIND(0.0D0))
acsf(1:num(1),1:num(2),1:num(3)) = cmplx (0.0D0, 0.0D0, KIND=KIND(0.0D0))
dsi (1:num(1),1:num(2),1:num(3)) = 0.0d0
!                                                                       
!------ preset some tables, calculate average structure                 
!                                                                       
call fourier_lmn(eck,vi,inc,lmn,off_shift)
call four_cexpt 
call four_stltab 
!write(*,*) ' eck ll ', eck(:,1)
!write(*,*) ' eck lr ', eck(:,2)
!write(*,*) ' eck ul ', eck(:,3)
!write(*,*) ' eck tl ', eck(:,4)
!write(*,*) ' vi abs ', vi(:,1), uin
!write(*,*) ' vi ord ', vi(:,2), vin
!write(*,*) ' vi top ', vi(:,3), win
!write(*,*) ' OFF  x ', off_shift(:,1)
!write(*,*) ' OFF  y ', off_shift(:,2)
!write(*,*) ' OFF  z ', off_shift(:,3)
!write(*,*) ' inc num', inc, num
!write(*,*) ' lmn    ', lmn
IF (ier_num.ne.0) return 
!
call four_formtab
!
!----- Test for Friedels law, and reduce calculation time
!
call four_test_friedel
call four_csize (cr_icc, csize, lperiod, ls_xyz) 
call four_aver (ilots, fave, csize) 
!call four_fft_prep
!                                                                       
!------ loop over crystal regions (lots)                                
!                                                                       
loop_lots: DO nlot = 1, nlots 
!
   csf(1:num(1),1:num(2),1:num(3)) = cmplx (0.0D0, 0.0D0, KIND=KIND(0.0D0))
!                                                                       
   IF(lot_all) THEN
      ii      = nlot - 1 
      lbeg(3) = (ii  )/(csize(1)*csize(2)) + 1
      ii      = (ii  ) - (lbeg(3)-1)*(csize(1)*csize(2))
      lbeg(2) = (ii  )/(csize(1)) + 1
      ii      = (ii  ) - (lbeg(2)-1)*(csize(1))
      lbeg(1) = MOD(ii,csize(1)) + 1
   ELSE
      call four_ranloc (csize, lbeg) 
   ENDIF
   IF (four_log) then 
      IF (ilots.ne.LOT_OFF) then 
         WRITE (output_io, 2000) nlot, nlots, (lbeg (i), i = 1, 3) 
      ELSE 
         WRITE (output_io, 2010) nlot, nlots 
      ENDIF 
   ENDIF 
!                                                                       
!------ - loop over all different atom types                            
!                                                                       
   loop_atoms: DO iscat = 1, cr_nscat 
      call four_getatm (iscat, ilots, lbeg, ncell) 
      call four_strucf (iscat, .true.) 
!                                                                       
!------ --- Add this part of the structur factor to the total           
!                 
      csf (1:num(1),1:num(2),1:num(3)) = csf (1:num(1),1:num(2),1:num(3)) + tcsf (1:num(1),1:num(2),1:num(3))
!
      IF (four_log) then 
         WRITE (output_io, 3000) cr_at_lis (iscat), nxat 
      ENDIF 
      IF(ier_ctrlc) THEN
         ier_num = -14
         ier_typ = ER_COMM
         RETURN
      ENDIF
      IF(ier_num/=0) RETURN      ! An error occured or CTRL-C
   ENDDO loop_atoms
!                                                                       
!------ - subtract average structure factor, add intensity              
!
   csf(1:num(1),1:num(2),1:num(3)) = csf(1:num(1),1:num(2),1:num(3)) - &
                                    acsf(1:num(1),1:num(2),1:num(3)) 
   dsi(1:num(1),1:num(2),1:num(3)) = dsi(1:num(1),1:num(2),1:num(3)) + &
                          DBLE(      csf(1:num(1),1:num(2),1:num(3)) * &
                               conjg(csf(1:num(1),1:num(2),1:num(3)) ) ) 
!
   call four_apply_friedel
ENDDO loop_lots
!
!----- Reset for Friedels law, and reduce calculation time
!
call four_rese_friedel
!                                                                       
!------ if we had lots, normalise the intensities                       
!****************************************************************************************************************************************************************************                                                                       
IF (nlots.ne.1) then 
   dnorm = DBLE(cr_icc(1) * cr_icc(2) * cr_icc(3) ) / &
           DBLE(ncell * nlots)                                                 
   dsi(1:num(1),1:num(2),1:num(3))  = dnorm * dsi(1:num(1),1:num(2),1:num(3))
ENDIF 
!
call four_conv           ! Convolute diffraction pattern
!
call four_accumulate
call four_weight               ! Correct the relative weight of Bragg and diffuse
!                                                                       
call four_qinfo 
ss = seknds (ss) 
IF (four_log) then 
   WRITE (output_io, 4000) ss 
   res_para(1) = ss
   res_para(0) = 1
ENDIF 
!
!
!write(*,*) ' eck ll ', eck(:,1)
!write(*,*) ' eck lr ', eck(:,2)
!write(*,*) ' eck ul ', eck(:,3)
!write(*,*) ' eck tl ', eck(:,4)
!write(*,*) ' vi abs ', vi(:,1), uin
!write(*,*) ' vi ord ', vi(:,2), vin
!write(*,*) ' vi top ', vi(:,3), win
!write(*,*) ' OFF  x ', off_shift(:,1)
!write(*,*) ' OFF  y ', off_shift(:,2)
!write(*,*) ' OFF  z ', off_shift(:,3)
!write(*,*) ' inc num', inc, num
!write(*,*) ' lmn    ', lmn
!
!                                                                       
 2000 FORMAT     (/,' Lot # ',I4,'/',I4,' : starting at unit cell (',   &
     &                       3(I3,1X),')')                              
 2010 FORMAT     (/,' Lot # ',I4,'/',I4,' : complete crystal') 
 3000 FORMAT (  '                   Atom typ = ',A4,7X,'(# ',I9,' )') 
 4000 FORMAT     (/,' Elapsed time    : ',G13.6,' sec') 
!
END SUBROUTINE four_run
!
!**********************************************************************
!
subroutine four_run_nufft
!-
!  Interface to Non-uniform FFT 
!+
!
use crystal_mod
use diffuse_mod
use four_finufft 
USE fourier_conv_mod
use fourier_lmn_mod
use symmetrize_mod
!
use do_lanczos_mod
use errlist_mod
use param_mod
use precision_mod
use prompt_mod
use support_mod
use lib_write_mod
!
implicit none
!
real(kind=PREC_DP), parameter :: EPS=1.0D-8
real(kind=PREC_DP), parameter :: MAXSCALE=3.0D0      ! Maximum scale allowed by FINUFFT
integer               :: i, j, k  ! Dummy index
integer               :: ii       ! Dummy index
integer               :: iscat    ! Dummy scattering type
integer, dimension(:), allocatable :: nat      ! Atom number
integer               :: is_dim   ! Crystal is of this dimension
integer, dimension(3) :: ncells   ! Number unit cells in relevant region
integer, dimension(3) :: idims    ! Number of data points in FFT
integer, dimension(3) :: NQXYZ    ! Actual data points in FFT array
logical, dimension(3) :: ll_dim   ! This dimension is not flat
logical               :: four_is_new  ! The reciprocal space dimensions have changed
!
integer              , dimension(3)              :: iscales  ! Scaling if DELTA (hkl) /= 1/cr_icc
real(kind=PREC_DP)   , dimension(3)              :: scales  ! Scaling if DELTA (hkl) /= 1/cr_icc
real(kind=PREC_DP)   , dimension(3)              :: shift   ! Shift atom positions towards Center of mass
real(KIND=PREC_DP) :: ss                                    ! Timing variable
real(kind=PREC_DP)   , dimension(:,:), allocatable :: xpos    ! atom coordinates, fractional
real(kind=PREC_DP)   , dimension(:,:), allocatable :: ypos    ! atom coordinates, fractional
real(kind=PREC_DP)   , dimension(:,:), allocatable :: zpos    ! atom coordinates, fractional
complex(kind=PREC_DP), dimension(:,:,:), allocatable :: fcsf    ! complex structure factor from FINUFFT
!
real(kind=PREC_DP)   , dimension(:    ), allocatable :: infield_1d    ! input into LANCZOS
real(kind=PREC_DP)   , dimension(:,:  ), allocatable :: infield_2d    ! input into LANCZOS
real(kind=PREC_DP)   , dimension(:,:,:), allocatable :: infield_3d    ! input into LANCZOS
real(kind=PREC_DP)   , dimension(:    ), allocatable ::outfield_1d    ! input into LANCZOS
real(kind=PREC_DP)   , dimension(:,:  ), allocatable ::outfield_2d    ! input into LANCZOS
real(kind=PREC_DP)   , dimension(:,:,:), allocatable ::outfield_3d    ! input into LANCZOS
!
character(len=1), dimension(3) :: c_hkl
character(len=8), dimension(3) :: c_axes
data c_axes /'abscissa', 'ordinate', 'top axis'/
data c_hkl  /'h', 'k', 'l'/
!
ss = seknds (0.0) 
!
call four_cexpt 
!
is_dim = 0           ! Assume zero dimensional crystal
ll_dim = .false.     ! Assume all dimensions to be flat
do i=1,3
   if(cr_dim(i,2)-cr_dim(i,1)>eps  ) then
      is_dim = is_dim + 1                    ! Found a non-flat dimension
      ll_dim(i) = .true.                     ! dimension i is non-flat
   endif
enddo
!                                                         ! Get point ratio for reciprocal space
idims = 1                                                 ! Default to 1 data point along each axis in reciprocal space
iscales = 1                                               ! Default to scale 1
scales(1) = abs(vi(1,1)*real(cr_icc(1), kind=PREC_DP))    ! Currently parallel a*!
scales(2) = abs(vi(2,2)*real(cr_icc(2), kind=PREC_DP))    ! Currently parallel b*!
scales(3) = abs(vi(3,3)*real(cr_icc(3), kind=PREC_DP))    ! Currently parallel c*!
do j=1,3
   if(scales(j)>MAXSCALE) then                            ! Scales must be < MAXSCALE
      iscales(j) = int((scales(j)+EPS)/MAXSCALE) + 1      ! Adapt iscales to smallest possible integer
      scales(j)  = scales(j) / real(iscales(j),kind=PREC_DP) ! Corrrect scales> iscales*scales=original_scales
      idims(j)   = iscales(j)*(inc(j)-1)+1                ! Adapt data points for FINUFFT
   else
      iscales(j) = 1                                      ! Scales is small enough
      idims(j)   = iscales(j)*(inc(j)-1)+1                ! == idims = inc
   endif
enddo
!
!  Error checks: Scale = cr_icc*vi must be integer
!                1/vi              must be integer
if(fave/=0.0D0) then
   do i=1, 3
      if(    abs(vi(i,i)*real(cr_icc(i), kind=PREC_DP))-                        &
         nint(abs(vi(i,i)*real(cr_icc(i), kind=PREC_DP)))>0.0D0) then
        if(scales(i)-nint(scales(i))>EPS) then
           ier_num = -186
           ier_typ = ER_APPL
           ier_msg(1) = 'Increment vector along ' // c_axes(i)
           return
        endif
         if(1.D0/vi(i,i)-nint(1.D0/(vi(i,i)))>EPS ) then
            ier_num = -187
            ier_typ = ER_APPL
            ier_msg(1) = 'Increment vector along ' // c_axes(i)
            write(ier_msg(2),'(a7, a1, a4,f10.6)') '1/step(',c_hkl(i),') = ',1.D0/vi(i,i)
            return
         endif
         if(mod(inc(i),2)==1) then
            j = 1
         else
            j = 0
         endif
         if(abs(eck(i,1) + (inc(i)-j)/2*vi(i,i))>EPS) then  
            ier_num = -188
            ier_typ = ER_APPL
            ier_msg(1) = 'Increment vector along ' // c_axes(i)
            write(ier_msg(2), '(a17,a1,a3,f10.6)') 'Shift corners by ', c_hkl(i),  &
                              ' = ',eck(i,1) + (inc(i)-j)/2*vi(i,i)
            return
         endif
      endif
   enddo
endif
!
NQXYZ   = inc
!
! Calculate average structure factor
!
if(fave > 0.0D0) then
   call four_aver_finufft(is_dim, ll_dim)
   if(ier_num/=0) return
else
!
   acsf(1:num(1),1:num(2),1:num(3)) = cmplx (0.0D0, 0.0D0, KIND=KIND(0.0D0))
endif
!                                                                                                                                              
!------ preset some values                                              
!                                                                       
call four_layer(four_is_new)   ! copy eck, vi
call fourier_lmn(eck,vi,inc,lmn,off_shift)
call four_stltab               ! set up sin(theta)/lambda table
call four_formtab              ! Calculate atomic form factors
!
!  Clear Fourier arrays 
!
csf(1:num(1),1:num(2),1:num(3))  = cmplx (0.0D0, 0.0D0, KIND=KIND(0.0D0))
dsi(1:num(1),1:num(2),1:num(3))  = 0.0d0
!
!------------------------------------------------------------------------------------------------------------------------------------------------------
!q: What is the difference between num(1) and inc(1)?
!a: num(1) is the number of data points in the reciprocal space    !!!!!! ASK NEDER !!!!!!
!   inc(1) is the number of data points in the real space
!   num(1) = inc(1) + 1
!
shift(1) = -real(int((cr_dim(1,2)+cr_dim(1,1))*0.5D0), kind=PREC_DP)
shift(2) = -real(int((cr_dim(2,2)+cr_dim(2,1))*0.5D0), kind=PREC_DP)
shift(3) = -real(int((cr_dim(3,2)+cr_dim(3,1))*0.5D0), kind=PREC_DP)
cond_dim: if(is_dim==3) then            ! 3-D crystal 3333333333333333333333333333
!
   allocate(fcsf(1:inc(1),1:inc(2),1:inc(3)))
!
   allocate(nat (             0:cr_nscat))
   allocate(xpos(0:cr_natoms, 0:cr_nscat))
   allocate(ypos(0:cr_natoms, 0:cr_nscat))
   allocate(zpos(0:cr_natoms, 0:cr_nscat))
   ncells = cr_icc
   xpos = 0.0D0
   ypos = 0.0D0
   zpos = 0.0D0
   nat = 0
   loop_atom3: do i=1, cr_natoms               ! Loop over all atoms
      k = cr_iscat(i,1)
      nat(k) = nat(k) + 1
      xpos(nat(k),k) = cr_pos(1,i) + shift(1)
      ypos(nat(k),k) = cr_pos(2,i) + shift(2)
      zpos(nat(k),k) = cr_pos(3,i) + shift(3)
   enddo loop_atom3
   loop_scat3:do iscat=1,cr_nscat                 ! Loop over all atom types
!
      call four_strucf_3d(cr_natoms, nat(iscat), ncells, idims, NQXYZ, &
                         xpos(:,iscat), ypos(:,iscat),     &
                         zpos(:,iscat), cr_occ(iscat),iscales, scales, fcsf)
      call tcsf_form(iscat, .true., NQXYZ, fcsf) ! Multiply with atomic form factor
      csf(1:num(1),1:num(2),1:num(3)) = csf(1:num(1),1:num(2),1:num(3)) + fcsf(1:num(1),1:num(2),1:num(3))                           ! Add to complex structure factor
   enddo loop_scat3
   deallocate(nat)
   deallocate(xpos)
   deallocate(ypos)
   deallocate(zpos)
elseif(is_dim==2) then  cond_dim                 ! 2-D crystal 222222222222222222222222
   allocate(fcsf(1:inc(1),1:inc(2),1:inc(3)))
  ! 
   j = findloc(inc   , 1      , 1)      ! Find     flat dimension in reciprocal space
   if(j==3) then             ! 3 is flat, set dimension 3 to 1
      iscales(3) = 1
      scales(3)  = 1
      idims(3)   = 1
   elseif(j==2) then         ! 2 is flat, 1 remains, 3 is copied down to 2
      iscales(2) = iscales(3)
      iscales(3) = 1
      scales(2)  = scales(3)
      scales(3)  = 1
      idims(2)   = idims(3)
      idims(3)   = 1
   elseif(j==1) then         ! 1 is flat, 2;3 each copied down one element
      iscales(1) = iscales(2)
      iscales(2) = iscales(3)
      iscales(3) = 1
      scales(1)  = scales(2)
      scales(2)  = scales(3)
      scales(3)  = 1
      idims(1)   = idims(2)
      idims(2)   = idims(3)
      idims(3)   = 1
   endif
   idims(3) = 1
   j = findloc(ll_dim, .false., 1)      ! Find     flat dimension in direct space
   allocate(nat (             0:cr_nscat))
   allocate(xpos(0:cr_natoms, 0:cr_nscat))
   allocate(ypos(0:cr_natoms, 0:cr_nscat))
   allocate(zpos(0:cr_natoms, 0:cr_nscat))
   xpos = 0.0D0
   ypos = 0.0D0
   zpos = 0.0D0
   nat = 0
   loop_atom2: do i=1, cr_natoms               ! Loop over all atoms
      k = cr_iscat(i,1)
      nat(k) = nat(k) + 1
      xpos(nat(k),k) = cr_pos(1,i) + shift(1)
      ypos(nat(k),k) = cr_pos(2,i) + shift(2)
      zpos(nat(k),k) = cr_pos(3,i) + shift(3)
   enddo loop_atom2
   loop_scat2:do iscat=1,cr_nscat                 ! Loop over all atom types
!
      fcsf = cmplx (0.0D0, 0.0D0, KIND=KIND(0.0D0))
!
      if(j==3) then                               ! x-y crystal
         ncells(1) = cr_icc(1)
         ncells(2) = cr_icc(2)
         ncells(3) = 1
         call four_strucf_2d(cr_natoms, nat(iscat), ncells(1:2), idims(1:2), NQXYZ,    &
                             xpos(:,iscat), ypos(:,iscat), cr_occ(iscat),iscales(1:2), scales(1:2), fcsf)
      elseif(j==2) then                           ! x-z crystal
         ncells(1) = cr_icc(1)
         ncells(2) = cr_icc(3)
         ncells(3) = 1
         call four_strucf_2d(cr_natoms, nat(iscat), ncells(1:2), idims(1:2), NQXYZ,    &
                             xpos(:,iscat), zpos(:,iscat), cr_occ(iscat),iscales(1:2), scales(1:2), fcsf)
      elseif(j==3) then                           ! y-z crystal
         ncells(1) = cr_icc(2)
         ncells(2) = cr_icc(3)
         ncells(3) = 1
         call four_strucf_2d(cr_natoms, nat(iscat), ncells(1:2), idims(1:2), NQXYZ,    &
                             ypos(:,iscat), zpos(:,iscat), cr_occ(iscat),iscales(1:2), scales(1:2), fcsf)
      endif
!call tofile((/num(1),num(2)/),'nufft_2d.dat',real(fcsf(:,:,1),kind=PREC_DP),(/-4.0D0, -4.0D0/),(/-0.1D0, -0.1D0/))
      call tcsf_form(iscat, .true., NQXYZ, fcsf) ! Multiply with atomic form factor
!                                                ! Add to complex structure factor
      csf(1:num(1),1:num(2),1:num(3)) = csf(1:num(1),1:num(2),1:num(3)) + fcsf(1:num(1),1:num(2),1:num(3))
   enddo loop_scat2
!
   deallocate(nat)
   deallocate(xpos)
   deallocate(ypos)
   deallocate(zpos)
!
elseif(is_dim==1) then cond_dim                   ! 1-D crystal 11111111111111111111111
   !
   allocate(fcsf(1:inc(1),1:inc(2),1:inc(3))) 
   j = findloc(ll_dim, .true., 1)                 ! Find non-flat dimension
   allocate(nat(0:cr_nscat))
   allocate(xpos(0:cr_natoms, 0:cr_nscat))
!
      nat = 0
      xpos = 0.0D0
      loop_atom1: do i=1, cr_natoms               ! Loop over all atoms
         k = cr_iscat(i,1)
         nat(k) = nat(k) + 1                         ! Increment atom number
         xpos(nat(k),k) = cr_pos(j,i) + shift(j)    ! Copy atoms
      enddo loop_atom1
   loop_scat1:do iscat=1,cr_nscat                 ! Loop over all atom types
      !
      fcsf = cmplx (0.0D0, 0.0D0, KIND=KIND(0.0D0))
      !
      call four_strucf_1d(cr_natoms, nat(iscat), cr_icc(j), idims(j), NQXYZ, xpos(:,iscat),     &
                          cr_occ(iscat),iscales(j), scales(j), fcsf)
      call tcsf_form(iscat, .true., NQXYZ, fcsf)  ! Multiply with atomic form factor
!                                                 ! Add to complex structure factor
      !
      csf(1:num(1),1:num(2),1:num(3)) = csf(1:num(1),1:num(2),1:num(3)) + fcsf(1:num(1),1:num(2),1:num(3))   ! And yet, here we go again to num(i).
      !
   enddo loop_scat1
   deallocate(nat)
   deallocate(xpos)
!
endif cond_dim
!
deallocate(fcsf)          
!
!
! Calculate average structure factor
!
!AV!if(fave > 0.0D0) then
!AV!   call four_aver_finufft(is_dim, ll_dim)
!AV!write(*,*) ' DONE  FINUFFT 2D AVERAGE', ier_num
!AV!   if(ier_num/=0) return
!AV!else
!AV!   acsf    = cmplx(0.0D0, 0.0D0, kind=kind(1.0D0))
!AV!endif
csf(1:num(1),1:num(2),1:num(3)) = csf(1:num(1),1:num(2),1:num(3)) - acsf(1:num(1),1:num(2),1:num(3))
do i=1, num(1)
   do j=1, num(2)
      do k=1, num(3)
         dsi(i,j,k) = dble(csf(i,j,k) * conjg(csf(i,j,k)))
      enddo
   enddo
enddo
!
if(four_filter==FOUR_FILTER_LANCZOS) then
   cond_dim_b: if(is_dim==3) then            ! 3-D crystal 3333333333333333333333333333
      allocate( infield_3d(num(1), num(2), num(3)))
      allocate(outfield_3d(num(1), num(2), num(3)))
      ii = 0
      do i=1, num(1)
         do j=1, num(2)
            do k=1, num(3)
               ii = ii + 1
               infield_3d(i,j,k) = dsi(i,j,k)
            enddo
         enddo
      enddo
      call do_lanczos(four_rscale, four_damp, four_width, num, infield_3d,   &
                   num, outfield_3d, .true.)
      ii = 0
      do i=1, num(1)
         do j=1, num(2)
            do k=1, num(2)
               ii = ii + 1
               dsi(i,j,k) = outfield_3d(i,j,k) 
            enddo
         enddo
      enddo
      deallocate( infield_3d)
      deallocate(outfield_3d)
   elseif(is_dim==2) then  cond_dim_b
      allocate( infield_2d(num(1), num(2)))
      allocate(outfield_2d(num(1), num(2)))
      ii = 0
      do i=1, num(1)
         do j=1, num(2)
            ii = ii + 1
            infield_2d(i,j) = dsi(i,j,1)
         enddo
      enddo
      call do_lanczos(four_rscale, four_damp, four_width, num(1:2), infield_2d,   &
                   num(1:2), outfield_2d, .true.)
      ii = 0
      do i=1, num(1)
         do j=1, num(2)
            ii = ii + 1
            dsi(i,j,1) = outfield_2d(i,j) 
         enddo
      enddo
      deallocate( infield_2d)
      deallocate(outfield_2d)
   elseif(is_dim==1) then  cond_dim_b
      allocate( infield_1d(num(1)))
      allocate(outfield_1d(num(1)))
      ii = 0
      do i=1, num(1)
         ii = ii + 1
         infield_1d(i) = dsi(i,1,1)
      enddo
      call do_lanczos(four_rscale, four_damp, four_width, num(1), infield_1d,   &
                   num(1), outfield_1d, .true.)
      ii = 0
      do i=1, num(1)
         ii = ii + 1
         dsi(i,1,1) = outfield_1d(i) 
      enddo
      deallocate( infield_1d)
      deallocate(outfield_1d)
   endif cond_dim_b
endif
!
call four_conv           ! Convolute diffraction pattern
call four_accumulate
call four_weight               ! Correct the relative weight of Bragg and diffuse
!
ss = seknds (ss) 
if(four_log) then 
   write(output_io, 4000) ss 
   res_para(1) = ss
   res_para(0) = 1
endif 
!
 4000 FORMAT     (/,' Elapsed time    : ',G13.6,' sec') 
!
end subroutine four_run_nufft
!
!**********************************************************************
!**********************************************************************
!
subroutine four_aver_finufft(is_dim, ll_dim)
!-
!  Calculate average structure factor and place into acsf
!  Required steps
!  Determine required integer range of HKL
!  calculate structure factors for these hkl
!  Determine indices of these hkl and place into acsf
!
!  For right now vectors vi must be parallel to a*, b*, c*
!+
use crystal_mod
use diffuse_mod
!
use errlist_mod
use matrix_mod
!
implicit none
!
integer                         , intent(in) :: is_dim  ! Reciprocal dimension is 1, 2, 3
logical           , dimension(3), intent(in) :: ll_dim  ! This dimension is not flat
!
real(kind=PREC_DP), parameter :: EPS = 1.0D-9
integer :: i,j,   ii, jj, kk
integer :: ih,ik,il
integer :: h,k,l                              ! Bragg indices
!real(kind=PREC_DP)                 :: rjj     ! Real valued version of jj
real(kind=PREC_DP), dimension(3,3) :: vi_tmp  ! temporary vi, augmented to give a determinant /= 0
real(kind=PREC_DP), dimension(3,3) :: vi_inv  ! inverse matrix to vi
real(kind=PREC_DP), dimension(3)   :: hkl     ! Indices    hkl
real(kind=PREC_DP), dimension(3)   :: vi_ijk  ! Indices of hkl
complex(kind=PREC_DP), dimension(:,:,:), allocatable :: acsf_hkl   ! complex structure factor for Bragg
!
!write(*,*) ' IS_DIM, LL_DIM', is_dim, ll_dim
!write(*,*) ' inc            ', inc
!write(*,*) ' ecḱ ', eck(:,1)
!write(*,*) ' ecḱ ', eck(:,2)
!write(*,*) ' ecḱ ', eck(:,3)
!write(*,*) ' ecḱ ', eck(:,4)
if(is_dim==3) then             ! 3D-crystal  3333333333333333
   vi_tmp = vi                 ! No need for further action
   if(inc(1)==1) vi_tmp(1,1) = 1.0D0
   if(inc(2)==1) vi_tmp(2,2) = 1.0D0
   if(inc(3)==1) vi_tmp(3,3) = 1.0D0
elseif(is_dim==2) then         ! 2D-crystal  2222222222222222
   vi_tmp = vi                 ! Set most elements
   if(inc(1)==1) vi_tmp(1,1) = 1.0D0
   if(inc(2)==1) vi_tmp(2,2) = 1.0D0
   vi_tmp(3,3) = 1.0D0
elseif(is_dim==1) then         ! 1D-crystal  1111111111111111
   vi_tmp = vi                 ! Set most elements
   vi_tmp(2,2) = 1.0D0
   vi_tmp(3,3) = 1.0D0
endif
!write(*,*) ' VI ', vi(:,1)
!write(*,*) ' VI ', vi(:,2)
!write(*,*) ' VI ', vi(:,3)
call matinv(vi_tmp, vi_inv)    ! Inverse matrix to place integer hkl indices into full acsf
if(ier_num/=0) then
   ier_msg(1) = 'Increment vectors cannot be inverted'
   ier_msg(2) = 'Are increment vectors parallel ?'
   write(ier_msg(3),'(a,i1,a)') 'The crystal is ',is_dim,'D'
   return
endif
!
do j = 1, 4
   do i =1, 3
      diff_eck_hkl(i,j) = real(int(eck(i,j)),kind=PREC_DP)  ! HKL corners are integer positions
   enddo
enddo
diff_vi_hkl       =  0.0D0     ! can always be a unit matrix
diff_vi_hkl(1,1)  =  1.0D0
diff_vi_hkl(2,2)  =  1.0D0
diff_vi_hkl(3,3)  =  1.0D0
diff_inc_hkl(1)   = nint(diff_eck_hkl(1,2)-diff_eck_hkl(1,1)) + 1
diff_inc_hkl(2)   = nint(diff_eck_hkl(2,3)-diff_eck_hkl(2,1)) + 1
diff_inc_hkl(3)   = nint(diff_eck_hkl(3,4)-diff_eck_hkl(3,1)) + 1
!
call four_aver_hkl
!
! Zero averace complex structure factor arrays, copy Bragg into  local acsf_hkl
allocate(acsf_hkl(1:diff_inc_hkl(1),1:diff_inc_hkl(2),1:diff_inc_hkl(3)))    ! Always a 3D array 
acsf_hkl(1:diff_inc_hkl(1),1:diff_inc_hkl(2),1:diff_inc_hkl(3)) &            ! Copy temporary acsf with Bragg into local acsf_hkl
  = acsf(1:diff_inc_hkl(1),1:diff_inc_hkl(2),1:diff_inc_hkl(3))              ! sub-array size is for Bragg only
acsf(1:num(1),1:num(2),1:num(3)) = cmplx (0.0D0, 0.0D0, KIND=KIND(0.0D0))    ! Zero the full acsf array
!
!
!  The triple loop runs over all Bragg reflections
!  hkl contains the real valued vector from the left lower bottom corner to Bragg reflection HKL
!  ih,ik,il are the indices of the Bragg reflection h,k,l in "acsf_hkl"
!  vi_ij    are the coordinated of Bragg reflection h,k,l in "acsf"
!  Only if these real valued indices are close enough to integer, the Bragg reflection is stored
do h = nint(diff_eck_hkl(1,1)), nint(diff_eck_hkl(1,2))        ! Loop over H component all Bragg reflection
   hkl(1) = h - eck(1,1)
   ih     = h - nint(diff_eck_hkl(1,1)) + 1                    ! Index of H in acsf_hkl
   do k = nint(diff_eck_hkl(2,1)), nint(diff_eck_hkl(2,3))     ! Loop over K component all Bragg reflection
      hkl(2) = k - eck(2,1)
      ik     = k - nint(diff_eck_hkl(2,1)) + 1                 ! Index of K in acsf_hkl
      do l = nint(diff_eck_hkl(3,1)), nint(diff_eck_hkl(3,4))  ! Loop over L component all Bragg reflection
         hkl(3) = l - eck(3,1)
         il     = l - nint(diff_eck_hkl(3,1)) + 1              ! Index of L in acsf_hkl
!
         vi_ijk = matmul(vi_inv, hkl) + 1                      ! Calculate coordinates of HKL in full acsf
         ii = nint(vi_ijk(1))
         jj = nint(vi_ijk(2))
         kk = nint(vi_ijk(3))
         if(abs(vi_ijk(1)-ii)<EPS .and. abs(vi_ijk(2)-jj)<EPS .and. abs(vi_ijk(3)-kk)<EPS) then
            acsf(ii, jj, kk) = acsf_hkl(ih,ik,il)
         endif
      enddo
   enddo
enddo
!
deallocate(acsf_hkl)
!
end subroutine four_aver_finufft
!
!****************************************************************************************************
!****************************************************************************************************
!
!****************************************************************************************************
! PENDING: Ask Neder about K, dimensions of fcsf, cfact, istl, iscatt, and the dimensions of the loops.
!****************************************************************************************************
!
subroutine tcsf_form(iscat, lform, NQXY, fcsf)
!-
!------ Now we multiply with formfactor                                 
!+
!
use diffuse_mod
!
implicit none
!
integer,               intent(in) :: iscat  ! Current scattering type
logical,               intent(in) :: lform  ! Use form factor
integer, dimension(3), intent(in) :: NQXY   ! Dimensions in reciprocal space
complex(kind=PREC_DP), dimension(NQXY(1), NQXY(2), NQXY(3)), intent(inout) :: fcsf  ! complex structure factor from FINUFFT
!
!
integer :: i,j,k   ! dummy loop index
!
if(lform) then
   !
!  ii = 0
   do i = 1, num(1)
      do j = 1, num(2)
         do k = 1, num(3)
!           ii = ii + 1
            fcsf(i,j,k) = fcsf(i,j,k) * cfact(istl(i,j,k), iscat) 
         enddo
      enddo
   enddo    
endif
!
end subroutine tcsf_form
!
!****************************************************************************************************
!****************************************************************************************************
!
!OBSOLETE subroutine tcsf_form_3D(iscat, lform, idims, ffcsf)
!OBSOLETE !-
!OBSOLETE !------ Now we multiply with formfactor                                 
!OBSOLETE !+
!OBSOLETE !
!OBSOLETE use diffuse_mod
!OBSOLETE !
!OBSOLETE implicit none
!OBSOLETE !
!OBSOLETE integer, intent(in) :: iscat        ! Current scattering type
!OBSOLETE logical, intent(in) :: lform        ! Use form factor
!OBSOLETE integer, dimension(3), intent(in) :: idims            ! Current scattering type
!OBSOLETE complex(kind=PREC_DP), dimension(idims(1), idims(2), idims(3)), intent(inout) :: ffcsf  ! complex structure factor from FINUFFT
!OBSOLETE !
!OBSOLETE integer :: i   ! dummy loop index
!OBSOLETE integer :: h,k,l ! dummy loop index
!OBSOLETE !
!OBSOLETE if(lform) then 
!OBSOLETE    i = 0
!OBSOLETE    do h=1, inc(1)
!OBSOLETE    do k=1, inc(2)
!OBSOLETE    do l=1, inc(3)
!OBSOLETE    i = i + 1
!OBSOLETE       ffcsf(h,k,l) = ffcsf(h,k,l) * cfact(istl(i), iscat)                                ! Something similar to what I do is done here, maybe it works.
!OBSOLETE    enddo
!OBSOLETE    enddo
!OBSOLETE    enddo
!OBSOLETE endif
!OBSOLETE !
!OBSOLETE end subroutine tcsf_form_3D
!
!**********************************************************************
!**********************************************************************
!
SUBROUTINE four_aver (lots, ave, csize) 
!+                                                                      
!     This routine calculates the average structure factor <F>.         
!-                                                                      
USE discus_config_mod 
USE crystal_mod 
USE celltoindex_mod
USE diffuse_mod 
USE four_strucf_mod
!     USE modify_mod
USE random_mod
!                                                                       
USE lib_random_func
USE prompt_mod 
USE precision_mod
!
IMPLICIT none 
       
!
INTEGER              , INTENT(IN) :: lots
REAL(kind=PREC_DP)   , INTENT(IN) :: ave 
INTEGER, DIMENSION(3), INTENT(IN) :: csize
!                                                                       
REAL(KIND=PREC_DP) :: norm
INTEGER :: isite, iatom, iscat, icell (3) 
INTEGER :: scell, ncell, i, j, k, ii, jj, kk 
LOGICAL, DIMENSION(:,:,:), ALLOCATABLE :: sel_cell
!                                                                       
      ave_is_zero: IF (ave.eq.0.0) then 
         RETURN 
      ELSE ave_is_zero  
         IF (four_log) then 
            WRITE (output_io, 1000) 100.0 * ave 
         ENDIF 
         scell = csize (1) * csize (2) * csize (3) 
         ncell = 0 
         ALLOCATE(sel_cell(csize(1), csize(2), csize(3)))
         sel_cell(:,:,:) = .FALSE.
         DO ii = 1, csize (1) 
            DO jj = 1, csize (2) 
               DO kk = 1, csize (3) 
                  sel_cell(ii,jj,kk) = (ran1 (idum) .le.ave) 
               ENDDO
            ENDDO
         ENDDO
!                                                                       
!------ ----- Loop over all atom types                                  
!                                                                       
            loop_iscat: DO iscat = 1, cr_nscat 
               nxat = 0 
!                                                                       
!------ - Loop over all unit cells                                      
!                                                                       
         cell_x: DO ii = 1, csize (1) 
         cell_y: DO jj = 1, csize (2) 
         cell_z: DO kk = 1, csize (3) 
         icell (1) = ii 
         icell (2) = jj 
         icell (3) = kk 
!                                                                       
!------ --- get only 'ave'% of those unit cells and compute average     
!------ --- unit cell                                                   
!                                                                       
!        sel = (ran1 (idum) .le.ave) 
!        IF (sel) then 
         IF (sel_cell(ii,jj,kk)) then 
            ncell = ncell + 1 
!RBN!                                                                       
!RBN!------ ----- Loop over all atom types                                  
!RBN!                                                                       
!RBN            loop_iscat: DO iscat = 1, cr_nscat 
!RBN               nxat = 0 
!                                                                       
!------ ------- Loop over all sites within unit cell                    
!                                                                       
               DO isite = 1, cr_ncatoms 
                  call celltoindex (icell, isite, iatom) 
                  IF (cr_iscat (iatom,1) .eq.iscat) then 
                     nxat = nxat + 1 
                     DO j = 1, 3 
                     xat (nxat, j) = cr_pos (j, iatom) - REAL(icell (j)     &
                     - 1) - cr_dim0 (j, 1)                                    
                     ENDDO 
                  ENDIF 
               ENDDO 
!RBN               call four_strucf_aver (iscat, .true.) 
!RBN               DO j = 1, num (1) * num (2) *num (3)
!RBN                  acsf (j) = acsf (j) + tcsf (j) 
!RBN               ENDDO 
!RBN            ENDDO  loop_iscat
         ENDIF 
         ENDDO cell_z
         ENDDO cell_y
         ENDDO cell_x
               call four_strucf (iscat, .true.)
               !
               ! Neders original declaration
               !
               !DO j = 1, num (1) * num (2) *num (3)
               !   acsf (j) = acsf (j) + tcsf (j) 
               !ENDDO
               !
               ! My declaration
               ! 
               DO i = 1, num (1)
                  DO j = 1, num (2)
                     DO k = 1, num (3)
                        acsf (i,j,k) = acsf (i,j,k) + tcsf (i,j,k)                                
                     ENDDO
                  ENDDO
               ENDDO
               !    
            ENDDO  loop_iscat
         ncell = ncell/cr_nscat
         DEALLOCATE(sel_cell)
!                                                                       
!------ - now compute the interference function of the lot shape        
!                                                                       
         call four_getav (lots) 
!        call four_strucf_aver (0, .false.) 
         call four_strucf      (0, .false.) 
         IF(ncell >0) THEN
            norm = DBLE(1.0D0 / ncell)
         ELSE
            norm = 0.0D0
            ier_num = +1
            ier_typ = ER_FOUR
            ier_msg(1) = 'Does the crystal consist of just 1 unit cell?'
            ier_msg(2) = 'Increase the percentage for >set aver<'
         ENDIF
         !
         ! Neders original declaration
         !
         !DO j = 1, num (1) * num (2) * num (3)
         !   acsf (j) = acsf (j) * tcsf (j) * cmplx ( norm, 0.0D0, KIND=KIND(0.0D0))
         !ENDDO 
!        !
!        ! My declaration
         !
         DO i = 1, num (1)
            DO j = 1, num (2)
               DO k = 1, num (3)
                  acsf (i,j,k) = acsf (i,j,k) * tcsf (i,j,k) * cmplx ( norm, 0.0D0, KIND=KIND(0.0D0))                              
               ENDDO
            ENDDO
         ENDDO     
         !                                                         
!------ - write how much of the crystal we actually used                
!                                                                       
         IF (four_log) then 
            WRITE (output_io, 2000) (REAL(ncell) / REAL(scell) )    &
            * 100.0                                                     
         ENDIF 
!                                                                       
      ENDIF ave_is_zero
!                                                                       
 1000 FORMAT     (' Calculating <F> with ',F5.1,' % of the crystal ...') 
 2000 FORMAT     (' Used ',F5.1,' % of the crystal to calculate <F> ...') 
      END SUBROUTINE four_aver                      
!
!**********************************************************************
!**********************************************************************
!
subroutine four_aver_hkl
!-
!  Calculate average structure factor; just at HKL integer for finufft
!+
!
USE chem_mod
use chem_aver_mod
use crystal_mod
use diffuse_mod
USE fourier_lmn_mod
use four_strucf_mod
!
implicit none
!
integer :: iscat    ! Loop index scattering types
integer :: iatom    ! Loop index scattering types
integer           , dimension(1:3)      ::  tmp_inc      ! Temporary storage to preserve user settings
real(kind=PREC_DP), dimension(1:3, 1:4) ::  tmp_eck
real(kind=PREC_DP), dimension(1:3, 1:3) ::  tmp_vi
logical :: four_is_new  ! The reciprocal space dimensions have changed
!OBSOLETE logical :: do_aver
!
!OBSOLETE do_aver = .false.
!OBSOLETE !do_aver = .true.
!
tmp_inc = inc
tmp_eck = eck
tmp_vi  = vi
inc     = diff_inc_hkl
eck     = diff_eck_hkl
vi      = diff_vi_hkl
!                                                                       
!------ preset some values                                              
!                                                                       
call four_layer(four_is_new)   ! copy eck, vi
call fourier_lmn(eck,vi,inc,lmn,off_shift)
call four_stltab               ! set up sin(theta)/lambda table
call four_formtab              ! Calculate atomic form factors
!
!  Clear Fourier arrays 
!
acsf(1:num(1),1:num(2),1:num(3)) = cmplx (0.0D0, 0.0D0, KIND=KIND(1.0D0))
!
! Using the averaged structure does not seem to subtract the Bragg reflections properly
! the cond_aver construction is thus obsolete
!OBSOLETE !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!OBSOLETE cond_aver: if(do_aver) then
!OBSOLETE !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!OBSOLETE call chem_aver(.true., .false.)
!OBSOLETE    loop_atoms2: do iatom = 1, cr_ncatoms
!OBSOLETE       do k = 1, chem_ave_n(iatom)
!OBSOLETE       nxat = 1
!OBSOLETE       xat(nxat,:) = chem_ave_posit(:,iatom,k)
!OBSOLETE       iscat = chem_ave_iscat(iatom, k)
!OBSOLETE       call four_strucf(iscat, .true.)
!OBSOLETE       acsf = acsf + tcsf * chem_ave_bese(iatom, k)
!OBSOLETE       enddo
!OBSOLETE    enddo loop_atoms2
!OBSOLETE !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!OBSOLETE else cond_aver
!OBSOLETE !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!
!------ ----- Loop over all atom types                                  
!
loop_iscat: do iscat = 1, cr_nscat       ! Loop over all scattering types
   nxat = 0 
   loop_atoms: do iatom = 1, cr_natoms   ! Loop over all atoms in crystal
     if(cr_iscat(iatom,1) == iscat) then   ! If correct type, copy into fourier position
         nxat = nxat + 1
         xat(nxat,:) = cr_pos(:, iatom) !-floor(cr_pos(:, iatom)) ! - real(icell(:) - 1,kind=PREC_DP) - cr_dim0(:, 1)
     endif
   enddo loop_atoms
   call four_strucf(iscat, .true.) 
   acsf(1:num(1), 1:num(2), 1:num(3)) = acsf(1:num(1), 1:num(2), 1:num(3)) + &
                                        tcsf(1:num(1), 1:num(2), 1:num(3))
enddo loop_iscat
!OBSOLETE !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!OBSOLETE endif cond_aver
!OBSOLETE !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!
inc = tmp_inc
eck = tmp_eck
vi  = tmp_vi
!
end subroutine four_aver_hkl                  
!
!**********************************************************************
!**********************************************************************
!
subroutine four_weight
!-
!  Correct the weight of Bragg versus diffuse if fave /= 100 and
!  the grid size or lot size does not correspond to 1/cr_icc.
!
!  This is needed only if the full 3D-PDF is calculated.  
!
!  If the volume of a voxel in reciprocal space does not 
!  correspond to the volume of the reciprocal unit cell divided by
!  the number of unit cell dimensions**3 we need to scale
!
!  2D 
!  If the area of the 2D-Pixel does not correspond to the area of
!  the area of the reciprocal unit cell section divided by the 
!  number of unit cell dimensions**2 we need to scale
!
!  1D
!  If the length os a step in reciprocal space does not
!  correspond to the length of the reciprocal unit cell divided by
!  the number of unit cells we need to scale
!
!+
!
use crystal_mod
use diffuse_mod
use metric_mod
!
use lib_errlist_func
use precision_mod
!
implicit none
!
integer               :: isdim      ! Fourier was calculated in 1,2,3 dimensions
integer               :: ncells
integer               :: i, j, k, ii
real(kind=PREC_DP)    :: scalef     ! scale factor to apply to the Bragg reflections
real(kind=PREC_DP), dimension(3  ) :: hkl     ! base vectors to calculate a volume
real(kind=PREC_DP), dimension(3,3) :: bases   ! base vectors to calculate a volume
real(kind=PREC_DP)                 :: voxel
!
scalef = 1.0D0
if_fave:if(fave == 0.000) then
!
   isdim = 3
   IF(num(1)==1) isdim = isdim - 1
   IF(num(2)==1) isdim = isdim - 1
   IF(num(3)==1) isdim = isdim - 1
!write(*,*) ' IN four_weight ISDIM ', num, isdim, dsort
!write(*,*) ' eck ll ', eck(:,1)
!write(*,*) ' eck lr ', eck(:,2)
!write(*,*) ' eck ul ', eck(:,3)
!write(*,*) ' eck tl ', eck(:,4)
!write(*,*) ' vi abs ', vi(:,1), uin
!write(*,*) ' vi ord ', vi(:,2), vin
!write(*,*) ' vi top ', vi(:,3), win
!write(*,*) ' OFF  x ', off_shift(:,1)
!write(*,*) ' OFF  y ', off_shift(:,2)
!write(*,*) ' OFF  z ', off_shift(:,3)
!write(*,*) ' inc num', inc, num
!
   ncells = 1
   if(isdim==1) then              ! 1D correction
      if(num(1)/=1) then
         bases(:,1) = uin
         bases(:,2) = 0.0D0
         bases(:,3) = 0.0D0
         if(bases(1,1)/=0D0) then
             bases(2,2)=1.0D0
             bases(3,3)=1.0D0
         elseif(bases(2,1)/=0D0) then
             bases(3,2)=1.0D0
             bases(1,3)=1.0D0
         elseif(bases(3,1)/=0D0) then
             bases(1,2)=1.0D0
             bases(2,3)=1.0D0
         endif
         if(nlots==1) then
            ncells = cr_icc(1)
         else
            ncells = ls_xyz(1)
         endif
      elseif(num(2)==1) then
         bases(:,1) = 0.0D0
         bases(:,2) = vin
         bases(:,3) = 0.0D0
         if(bases(1,2)/=0D0) then
             bases(2,1)=1.0D0
             bases(3,3)=1.0D0
         elseif(bases(2,2)/=0D0) then
             bases(3,1)=1.0D0
             bases(1,3)=1.0D0
         elseif(bases(3,2)/=0D0) then
             bases(1,1)=1.0D0
             bases(2,3)=1.0D0
         endif
         if(nlots==1) then
            ncells = cr_icc(2)
         else
            ncells = ls_xyz(2)
         endif
      elseif(num(3)==1) then
         bases(:,1) = 0.0D0
         bases(:,2) = 0.0D0
         bases(:,3) = win
         if(bases(1,3)/=0D0) then
             bases(2,1)=1.0D0
             bases(3,2)=1.0D0
         elseif(bases(2,3)/=0D0) then
             bases(3,1)=1.0D0
             bases(1,2)=1.0D0
         elseif(bases(3,3)/=0D0) then
             bases(1,1)=1.0D0
             bases(2,2)=1.0D0
         endif
         if(nlots==1) then
            ncells = cr_icc(3)
         else
            ncells = ls_xyz(3)
         endif
      endif
   elseif(isdim==2) then          ! 2D correction
      if(num(3)==1) then
         bases(:,1) = uin
         bases(:,2) = vin
         bases(:,3) = 1.0D0
         if(bases(1,1)/=0D0 .or. bases(1,2)/=0.0D0) bases(1,3)=0.0D0
         if(bases(2,1)/=0D0 .or. bases(2,2)/=0.0D0) bases(2,3)=0.0D0
         if(bases(3,1)/=0D0 .or. bases(3,2)/=0.0D0) bases(3,3)=0.0D0
         if(nlots==1) then
            ncells = cr_icc(1)*cr_icc(2)
         else
            ncells = ls_xyz(1)*ls_xyz(2)
         endif
      elseif(num(2)==1) then
         bases(:,1) = uin
         bases(:,2) = 1.0D0
         bases(:,3) = win
         if(bases(1,1)/=0D0 .or. bases(1,3)/=0.0D0) bases(1,2)=0.0D0
         if(bases(2,1)/=0D0 .or. bases(2,3)/=0.0D0) bases(2,2)=0.0D0
         if(bases(3,1)/=0D0 .or. bases(3,3)/=0.0D0) bases(3,2)=0.0D0
         if(nlots==1) then
            ncells = cr_icc(1)*cr_icc(3)
         else
            ncells = ls_xyz(1)*ls_xyz(3)
         endif
      elseif(num(3)==1) then
         bases(:,1) = 1.0D0
         bases(:,2) = vin
         bases(:,3) = win
         if(bases(1,2)/=0D0 .or. bases(1,3)/=0.0D0) bases(1,1)=0.0D0
         if(bases(2,2)/=0D0 .or. bases(2,3)/=0.0D0) bases(2,1)=0.0D0
         if(bases(3,2)/=0D0 .or. bases(3,3)/=0.0D0) bases(3,1)=0.0D0
         if(nlots==1) then
            ncells = cr_icc(2)*cr_icc(3)
         else
            ncells = ls_xyz(2)*ls_xyz(3)
         endif
      endif
!write(*,*) ' Base 1 ', bases(:,1)
!write(*,*) ' Base 2 ', bases(:,2)
!write(*,*) ' Base 3 ', bases(:,3)
!   voxel = do_volume(.FALSE., bases(:,1), bases(:,2), bases(:,3))
!   scalef = real(ncells)/cr_vr*voxel
!   write(*,*)  ' Voxel ', voxel
!   write(*,*)  ' Vr/Voxel ', cr_vr/voxel
!   write(*,*)  ' Ratio    ', real(ncells)/cr_vr*voxel, scalef
!
   elseif(isdim==3) then          ! 3D correction
      bases(:,1) = uin
      bases(:,2) = vin
      bases(:,3) = win
      if(nlots==1) then
         ncells = cr_icc(1)*cr_icc(2)*cr_icc(3)
      else
         ncells = ls_xyz(1)*ls_xyz(2)*ls_xyz(3)
   endif
   endif
!
   voxel = do_volume(.FALSE., bases(:,1), bases(:,2), bases(:,3))
   if(voxel>0.0) then
      scalef = real(ncells)/cr_vr*voxel
   else
      call no_error
   endif
endif if_fave
!
if(allocated(dsi3d)) deallocate(dsi3d)
allocate(dsi3d(1:ubound(dsi,1),1:ubound(dsi,2),1:ubound(dsi,3)))
if(fave==0.0 .and. abs(scalef-1.0)>1.0E-5) then
!write(*,*) 'base1  ', bases(:,1)
!write(*,*) 'base2  ', bases(:,2)
!write(*,*) 'base3  ', bases(:,3)
!write(*,*) 'Scalef ', scalef, ncells, cr_vr, voxel
!read(*,*) i
   dsi3d = dsi
   ii= 0
   do k=0,num(3)-1
      do j=0,num(2)-1
         do i=0,num(1)-1
            ii = ii + 1
            hkl(1) = eck(1,1) + vi(1,1)*i + vi(1,2)*j + vi(1,3)*k
            hkl(2) = eck(2,1) + vi(2,1)*i + vi(2,2)*j + vi(2,3)*k
            hkl(3) = eck(3,1) + vi(3,1)*i + vi(3,2)*j + vi(3,3)*k
            if(abs(hkl(1)-nint(hkl(1)))<1.0e-4 .and.                               &
               abs(hkl(2)-nint(hkl(2)))<1.0e-4 .and.                               &
               abs(hkl(3)-nint(hkl(3)))<1.0e-4       ) then
               dsi3d(i+1,j+1,k+1) = dsi(i+1,j+1,k+1)/scalef
            endif
         enddo
      enddo
   enddo
else
   dsi3d = dsi
endif
!
end subroutine four_weight
!
!**********************************************************************
!**********************************************************************
      SUBROUTINE four_getatm (iscat, lots, lbeg, ncell) 
!+                                                                      
!     This routine creates an atom list of atoms of type 'iscat'        
!     which are within the current lot.                                 
!-                                                                      
      USE discus_config_mod 
      USE crystal_mod 
      USE celltoindex_mod
      USE diffuse_mod 
!     USE modify_mod
!
use precision_mod
!
      IMPLICIT none 
!                                                                       
      INTEGER, INTENT(IN)  :: iscat
      INTEGER, INTENT(IN)  :: lots
      INTEGER, DIMENSION(3), INTENT(IN)  :: lbeg
      INTEGER, INTENT(OUT) :: ncell
!                                                                       
      REAL(kind=PREC_DP) :: offset (3) 
      REAL(kind=PREC_DP) :: x0 (3), xtest (3) 
      INTEGER cr_end 
      INTEGER icell (3), jcell (3) 
      INTEGER i, j, ir, ii, jj, kk, is, ia 
!                                                                       
      nxat = 0 
      ncell = 0 
      cr_end = cr_ncatoms * cr_icc (1) * cr_icc (2) * cr_icc (3) + 1
!                                                                       
!------ No lots required, return all atoms of type 'iscat'              
!                                                                       
      IF (lots.eq.LOT_OFF) then 
         ncell = cr_icc (1) * cr_icc (2) * cr_icc (3) 
         DO i = 1, cr_natoms 
         IF (cr_iscat (i,1) .eq.iscat) then 
            nxat = nxat + 1 
            DO j = 1, 3 
            xat (nxat, j) = cr_pos (j, i) 
            ENDDO 
         ENDIF 
         ENDDO 
!                                                                       
!------ Box shaped lot                                                  
!                                                                       
      ELSEIF (lots.eq.LOT_BOX) then 
         DO kk = 0, ls_xyz (3) - 1 
         icell (3) = kk + lbeg (3) 
         offset (3) = cr_dim0 (3, 1) + lbeg (3) - 1 
         IF (icell (3) .gt.cr_icc (3) ) then 
            icell (3) = icell (3) - cr_icc (3) 
            offset (3) = offset (3) - REAL(cr_icc (3) ) 
         ENDIF 
         DO jj = 0, ls_xyz (2) - 1 
         icell (2) = jj + lbeg (2) 
         offset (2) = cr_dim0 (2, 1) + lbeg (2) - 1 
         IF (icell (2) .gt.cr_icc (2) ) then 
            icell (2) = icell (2) - cr_icc (2) 
            offset (2) = offset (2) - REAL(cr_icc (2) ) 
         ENDIF 
         DO ii = 0, ls_xyz (1) - 1 
         icell (1) = ii + lbeg (1) 
         offset (1) = cr_dim0 (1, 1) + lbeg (1) - 1 
         IF (icell (1) .gt.cr_icc (1) ) then 
            icell (1) = icell (1) - cr_icc (1) 
            offset (1) = offset (1) - REAL(cr_icc (1) ) 
         ENDIF 
!                                                                       
         ncell = ncell + 1 
         DO is = 1, cr_ncatoms 
         call celltoindex (icell, is, ia) 
         IF (cr_iscat (ia,1) .eq.iscat) then 
            nxat = nxat + 1 
            DO j = 1, 3 
            xat (nxat, j) = cr_pos (j, ia) - offset (j) 
            ENDDO 
         ENDIF 
         ENDDO 
!                                                                       
         DO ir = cr_end, cr_natoms 
         DO j = 1, 3 
         jcell (j) = icell (j) + nint (cr_dim (j, 1) ) - 1 
         ENDDO 
         IF (int (cr_pos (1, ir) ) .eq.jcell (1) .and.int (cr_pos (2,   &
         ir) ) .eq.jcell (2) .and.int (cr_pos (3, ir) ) .eq.jcell (3)   &
         .and.cr_iscat (ir,1) .eq.iscat) then                             
            nxat = nxat + 1 
            DO j = 1, 3 
            xat (nxat, j) = cr_pos (j, ir) - offset (j) 
            ENDDO 
         ENDIF 
         ENDDO 
!                                                                       
         ENDDO 
         ENDDO 
         ENDDO 
!                                                                       
!------ Ellipsoid shaped lot                                            
!                                                                       
      ELSEIF (lots.eq.LOT_ELI) then 
         DO ii = 1, 3 
         x0 (ii) = REAL(ls_xyz (ii) ) / 2.0 
         ENDDO 
!                                                                       
         DO kk = 0, ls_xyz (3) - 1 
         icell (3) = kk + lbeg (3) 
         offset (3) = cr_dim0 (3, 1) + lbeg (3) - 1 
         xtest (3) = (REAL(kk) - x0 (3) + 0.5) **2 / x0 (3) **2 
         IF (icell (3) .gt.cr_icc (3) ) then 
            icell (3) = icell (3) - cr_icc (3) 
            offset (3) = offset (3) - REAL(cr_icc (3) ) 
         ENDIF 
         DO jj = 0, ls_xyz (2) - 1 
         icell (2) = jj + lbeg (2) 
         offset (2) = cr_dim0 (2, 1) + lbeg (2) - 1 
         xtest (2) = (REAL(jj) - x0 (2) + 0.5) **2 / x0 (2) **2 
         IF (icell (2) .gt.cr_icc (2) ) then 
            icell (2) = icell (2) - cr_icc (2) 
            offset (2) = offset (2) - REAL(cr_icc (2) ) 
         ENDIF 
         DO ii = 0, ls_xyz (1) - 1 
         icell (1) = ii + lbeg (1) 
         offset (1) = cr_dim0 (1, 1) + lbeg (1) - 1 
         xtest (1) = (REAL(ii) - x0 (1) + 0.5) **2 / x0 (1) **2 
         IF (icell (1) .gt.cr_icc (1) ) then 
            icell (1) = icell (1) - cr_icc (1) 
            offset (1) = offset (1) - REAL(cr_icc (1) ) 
         ENDIF 
!                                                                       
         IF ( (xtest (1) + xtest (2) + xtest (3) ) .le.1.0) then 
            ncell = ncell + 1 
            DO is = 1, cr_ncatoms 
            call celltoindex (icell, is, ia) 
            IF (cr_iscat (ia,1) .eq.iscat) then 
               nxat = nxat + 1 
               DO j = 1, 3 
               xat (nxat, j) = cr_pos (j, ia) - offset (j) 
               ENDDO 
            ENDIF 
            ENDDO 
!                                                                       
            DO ir = cr_end, cr_natoms 
            DO j = 1, 3 
            jcell (j) = icell (j) + nint (cr_dim (j, 1) ) - 1 
            ENDDO 
            IF (int (cr_pos (1, ir) ) .eq.jcell (1) .and.int (cr_pos (2,&
            ir) ) .eq.jcell (2) .and.int (cr_pos (3, ir) ) .eq.jcell (3)&
            .and.cr_iscat (ir,1) .eq.iscat) then                          
               nxat = nxat + 1 
               DO j = 1, 3 
               xat (nxat, j) = cr_pos (j, ir) - offset (j) 
               ENDDO 
            ENDIF 
            ENDDO 
         ENDIF 
!                                                                       
         ENDDO 
         ENDDO 
         ENDDO 
      ENDIF 
!                                                                       
      END SUBROUTINE four_getatm                    
!**********************************************************************
!**********************************************************************
      SUBROUTINE four_getav (lots) 
!+                                                                      
!     This routine computes the lattice inside the lot                  
!-                                                                      
      USE discus_config_mod 
      USE crystal_mod 
      USE diffuse_mod 
use precision_mod
!
      IMPLICIT none 
!                                                                       
      INTEGER, INTENT(IN) :: lots
!                                                                       
      REAL(kind=PREC_DP) :: xtest (3), x0 (3) 
      INTEGER ii, jj, kk 
!                                                                       
      nxat = 0 
!                                                                       
!------ No lots required, return complete lattice                       
!                                                                       
      IF (lots.eq.LOT_OFF) then 
         DO ii = 0, cr_icc (1) - 1 
         DO jj = 0, cr_icc (2) - 1 
         DO kk = 0, cr_icc (3) - 1 
         nxat = nxat + 1 
         xat (nxat, 1) = REAL(ii) + cr_dim0 (1, 1) 
         xat (nxat, 2) = REAL(jj) + cr_dim0 (2, 1) 
         xat (nxat, 3) = REAL(kk) + cr_dim0 (3, 1) 
         ENDDO 
         ENDDO 
         ENDDO 
!                                                                       
!------ Box shaped lot                                                  
!                                                                       
      ELSEIF (lots.eq.LOT_BOX) then 
         DO kk = 0, ls_xyz (3) - 1 
         DO jj = 0, ls_xyz (2) - 1 
         DO ii = 0, ls_xyz (1) - 1 
         nxat = nxat + 1 
         xat (nxat, 1) = REAL(ii) 
         xat (nxat, 2) = REAL(jj) 
         xat (nxat, 3) = REAL(kk) 
         ENDDO 
         ENDDO 
         ENDDO 
!                                                                       
!------ Ellipsoid shaped lot                                            
!                                                                       
      ELSEIF (lots.eq.LOT_ELI) then 
         DO ii = 1, 3 
         x0 (ii) = REAL(ls_xyz (ii) ) / 2.0 
         ENDDO 
!                                                                       
         DO kk = 0, ls_xyz (3) - 1 
         xtest (3) = (REAL(kk) - x0 (3) + 0.5) **2 / x0 (3) **2 
         DO jj = 0, ls_xyz (2) - 1 
         xtest (2) = (REAL(jj) - x0 (2) + 0.5) **2 / x0 (2) **2 
         DO ii = 0, ls_xyz (1) - 1 
         xtest (1) = (REAL(ii) - x0 (1) + 0.5) **2 / x0 (1) **2 
         IF ( (xtest (1) + xtest (2) + xtest (3) ) .le.1.0) then 
            nxat = nxat + 1 
            xat (nxat, 1) = REAL(ii) 
            xat (nxat, 2) = REAL(jj) 
            xat (nxat, 3) = REAL(kk) 
         ENDIF 
         ENDDO 
         ENDDO 
         ENDDO 
      ENDIF 
!                                                                       
      END SUBROUTINE four_getav                     
!**********************************************************************
!**********************************************************************
      SUBROUTINE four_ranloc (csize, lbeg) 
!+                                                                      
!     Returns a random cell from within the simulated crystal           
!     which is limited by 'csize'.                                      
!-                                                                      
USE lib_random_func
      USE random_mod
      IMPLICIT none 
!                                                                       
!                                                                       
      INTEGER, DIMENSION(3), INTENT(IN)  :: csize
      INTEGER, DIMENSION(3), INTENT(OUT) :: lbeg
!                                                                       
      lbeg (1) = int (ran1 (idum) * csize (1) ) + 1 
      lbeg (2) = int (ran1 (idum) * csize (2) ) + 1 
      lbeg (3) = int (ran1 (idum) * csize (3) ) + 1 
!                                                                       
      END SUBROUTINE four_ranloc                    
!**********************************************************************
!**********************************************************************
      SUBROUTINE four_csize (cr_icc, csize, lperiod, ls_xyz) 
!+                                                                      
!     Limits crystal size in case of no periodic boundary               
!     conditions.                                                       
!-                                                                      
      IMPLICIT none 
!                                                                       
      INTEGER, DIMENSION(3), INTENT(IN)  :: cr_icc
      INTEGER, DIMENSION(3), INTENT(OUT) :: csize
      LOGICAL,               INTENT(IN)  :: lperiod 
      INTEGER, DIMENSION(3), INTENT(IN)  :: ls_xyz
!                                                                       
      IF (lperiod) then 
         csize (1) = cr_icc (1) 
         csize (2) = cr_icc (2) 
         csize (3) = cr_icc (3) 
      ELSE 
         csize (1) = MAX(1, cr_icc (1) - ls_xyz (1) + 1 )
         csize (2) = MAX(1, cr_icc (2) - ls_xyz (2) + 1 )
         csize (3) = MAX(1, cr_icc (3) - ls_xyz (3) + 1 )
      ENDIF 
!                                                                       
      END SUBROUTINE four_csize                     
!**********************************************************************
!**********************************************************************
!
SUBROUTINE four_layer (four_is_new)
!+                                                                      
!     This routine writes the corners of the plane to be calculated     
!     in the 'diffuse_mod.f90' variables.                                   
!-                                                                      
USE discus_config_mod 
USE diffuse_mod 
!
IMPLICIT none 
!                                                                       
logical, intent(out) :: four_is_new
!
real(kind=PREC_DP), parameter :: EPS = 1.0D-6
real(kind=PREC_DP), dimension(1:3, 1:4), save :: prev_eck = 0.0D0
real(kind=PREC_DP), dimension(1:3, 1:3), save :: prev_vi  = 0.0D0
integer           , dimension(1:3)     , save :: prev_inc = 1
!                                                                       
INTEGER :: i 
!
if(any(abs(prev_eck-eck)>EPS) .or. any(abs(prev_vi-vi)>EPS) .or. any(prev_inc/=inc)) then
   four_is_new = .true.
else
   four_is_new = .false.
endif
!                                                                       
DO i = 1, 3 
   xm  (i) = eck (i, 1) 
   uin (i) = vi  (i, 1) 
   vin (i) = vi  (i, 2) 
   win (i) = vi  (i, 3) 
ENDDO 
!                                                                       
DO i = 1, 3 
   num (i) = inc (i) 
ENDDO 
!                                                                       
END SUBROUTINE four_layer                     
!
!**********************************************************************
!**********************************************************************
      SUBROUTINE four_cexpt 
!+                                                                      
!     This routine initialises the complex exponent table and           
!     is called only at the first Fourier run.                          
!-                                                                      
      USE discus_config_mod 
      USE diffuse_mod 
      USE prompt_mod 
      USE precision_mod 
      USE wink_mod
      IMPLICIT none 
!                                                                       
!                                                                       
      REAL(PREC_DP) xmult, xarg 
      INTEGER i 
!                                                                       
      IF (.not.ffour) then 
         WRITE (output_io, 1000) 
!                                                                       
!        zpi = 8.0d0 * datan (1.0d0) 
!                                                                       
         DO i = 0, MASK 
            xmult   = (dble (i) * 1.0d0) / dble (I2PI) 
            xarg    = zpi * xmult 
            cex (i) = CMPLX (DBLE( COS (xarg)), DBLE( SIN (xarg)), KIND=KIND(0.0D0) ) 
         ENDDO 
         ffour = .true. 
      ENDIF 
!                                                                       
 1000 FORMAT     (' Computing complex exponent table ...') 
      END SUBROUTINE four_cexpt                     
!
!**********************************************************************
!**********************************************************************
!
SUBROUTINE four_stltab 
!+                                                                      
!     Sets up an integer array containing the corresponding integers    
!     to the formfactor table for each sin(theta)/lambda.               
!-                                                                      
USE discus_config_mod 
USE crystal_mod 
USE diffuse_mod
!USE quad_mod 
use metric_mod ,only: skalpro
USE prompt_mod 
use precision_mod
!                                                                       
IMPLICIT none 
!                                                                       
REAL(kind=PREC_DP) :: q2, h (3) 
INTEGER :: i, j, k, l 
!                                                                       
IF (four_log) then 
   WRITE (output_io, 1000) 
ENDIF 
!                                                                       
DO l = 1, num (3) 
   DO j = 1, num (2) 
      DO i = 1, num (1) 
      DO k = 1, 3 
      h (k) = REAL(xm (k) + uin (k) * REAL(i - 1, KIND=KIND(0.0D0)) &
                          + vin (k) * REAL(j - 1, KIND=KIND(0.0D0)) &
                          + win (k) * REAL(l - 1, KIND=KIND(0.0D0)))
      ENDDO 
!     q2 = quad (h, h, cr_rten) / 4.0 
      q2 = skalpro (h, h, cr_rten) / 4.0D0
!     k  = (i - 1) * num (2) + j 
!     k  = (i - 1) * num (3)* num (2) + (j - 1) * num (3) + l 
!RBN_3D
!write(*,2000) i,j,l, k
!2000 format( 'ijl:k ', 3i4,' : ', i6)
      istl(i,j,l) = nint (sqrt (q2) * (1.0D0 / CFINC) ) 
!                                                                       
      IF (istl (i,j,l) .gt.CFPKT) then 
         ier_num = - 3 
         ier_typ = ER_FOUR 
         RETURN 
      ENDIF 
!                                                                       
      ENDDO 
   ENDDO 
ENDDO 
!                                                                       
 1000 FORMAT     (' Computing sin(theta)/lambda table ...') 
END SUBROUTINE four_stltab                    
!
!**********************************************************************
!**********************************************************************
SUBROUTINE four_formtab
!+                                                                      
!     This routine sets up the complex formfactor lookup table          
!     for all atom types. The range in sin(theta)/lambda is             
!     0 -> 2 in steps of 0.001. These values can be changed             
!     in the 'diffuse_mod.f90' file.                                        
!-                                                                      
USE discus_config_mod 
USE crystal_mod 
USE diffuse_mod 
!                                                                       
USE prompt_mod 
USE precision_mod
!
IMPLICIT none 
!
!                                                                       
REAL(kind=PREC_DP)  :: q2
REAL(KIND=PREC_DP) :: sb, sf, sfp, sfpp 
INTEGER iq, iscat 
!                                                                       
IF (four_log) then 
   WRITE (output_io, 1000) 
ENDIF 
!                                                                       
DO iscat = 1, cr_nscat 
   DO iq = 0, CFPKT 
      q2 =      (REAL(iq, KIND=KIND(0.0D0)) * CFINC) **2
      sf = DBLE(form (iscat, cr_scat, lxray, q2, diff_power) )
!                                                                       
      IF (ano) then 
         sfp  = DBLE(cr_delfr ( (iscat) ) )
         sfpp = DBLE(cr_delfi ( (iscat) ) )
      ELSE 
         sfp  = 0.0D0 
         sfpp = 0.0D0
      ENDIF 
!                                                                       
      IF (ldbw) then 
         sb = exp ( - DBLE(cr_dw ( (iscat) ) * q2)) * DBLE(cr_occ(iscat))
      ELSE 
         sb = 1.0D0 * DBLE(cr_occ(iscat))
      ENDIF 
!                                                                       
      cfact     (iq, iscat) = cmplx (sb * (sf + sfp), sb * sfpp, KIND=KIND(0.0D0)) 
      cfact_pure(iq, iscat) = cmplx (     (sf + sfp),      sfpp, KIND=KIND(0.0D0)) 
   ENDDO 
ENDDO 
!                                                                       
 1000 FORMAT     (' Computing formfactor lookup table ...') 
END SUBROUTINE four_formtab                   
!
!**********************************************************************
!
subroutine four_dbtab
!-
!  Set up Debye-Waller table for all elements and all points in reciprocal space
!+
!
use crystal_mod
use diffuse_mod 
!
use precision_mod
use prompt_mod
!
implicit none
!
integer :: i,j,k,l  ! Dummy indices
integer :: iscat
real(kind=PREC_DP), dimension(3) :: h     ! Reciprocal space vector
real(kind=PREC_DP)               :: arg   ! Dummy argument
!
if(four_log) then 
   write(output_io, '(a)') ' Computing DebeyWaller lookup table ...' 
endif 
!
do k=1, num(3)
   do j=1, num(2)
      do i=1, num(1)
         h(l) = eck(l, 1) + vi(l, 1) * REAL(i - 1,kind=PREC_DP) + &
                            vi(l, 2) * REAL(j - 1,kind=PREC_DP) + &
                            vi(l, 3) * REAL(k - 1,kind=PREC_DP)
         arg = (h(1)*cr_anis_full(1,iscat)) * (h(1)*cr_anis_full(1,iscat)) *   &
               (h(2)*cr_anis_full(2,iscat)) * (h(2)*cr_anis_full(2,iscat)) *   &
               (h(3)*cr_anis_full(3,iscat)) * (h(3)*cr_anis_full(3,iscat)) *   &
               (h(2)*cr_anis_full(4,iscat)) * (h(3)*cr_anis_full(4,iscat)) * 2.0D0 * &
               (h(1)*cr_anis_full(5,iscat)) * (h(3)*cr_anis_full(5,iscat)) * 2.0D0 * &
               (h(1)*cr_anis_full(6,iscat)) * (h(2)*cr_anis_full(6,iscat)) * 2.0D0
!        four_dbw(i,j,k,iscat) = exp(-arg)
      enddo
   enddo
enddo
!
end subroutine four_dbtab
!**********************************************************************
!
SUBROUTINE four_qinfo 
!+                                                                      
!     Gives information about max/min values for diffuse and            
!     Bragg scattering.                                                 
!-                                                                      
      USE discus_config_mod 
      USE diffuse_mod 
!                                                                     
      USE precision_mod
      USE param_mod
      USE prompt_mod 
!
      IMPLICIT none 
!                                                                       
      INTEGER i, j, k, l, nd 
      REAL(kind=PREC_DP) :: h (3), dsum, dsum2 
      LOGICAL lbragg 
!                                                                       
      braggmax = - 9999.0 
      braggmin = 1E19 
      diffumax = - 9999.0 
      diffumin = 1E19 
!                                                                       
      nd = 0 
      dsum = 0.0 
      dsum2 = 0.0 
!                                                                       
      DO l = 1, num (3) 
      DO j = 1, num (2) 
      DO i = 1, num (1) 
      lbragg = .true. 
      DO k = 1, 3 
      h (k) = eck (k, 1) + vi (k, 1) * REAL(i - 1) + &
                           vi (k, 2) * REAL(j - 1) + &
                           vi (k, 3) * REAL(l - 1)
      lbragg = lbragg.and.mod(h(k), 1.0D0) .eq.0.0 
      ENDDO 
!     k = (i - 1) * num (2) + j 
      k = (i - 1) * num (3) * num (2) + (j -1) * num(3) + l 
      IF (lbragg) then 
         braggmax = max(braggmax, dsi(i,j,l) ) 
         braggmin = min(braggmin, dsi(i,j,l) ) 
      ELSE 
         diffumax = max(diffumax, dsi(i,j,l) ) 
         diffumin = min(diffumin, dsi(i,j,l) ) 
!                                                                       
         dsum  = dsum  + dsi(i,j,l)
         dsum2 = dsum2 + dsi(i,j,l)**2 
         nd = nd+1 
      ENDIF 
      ENDDO 
      ENDDO 
      ENDDO 
!                                                                       
      IF (nd.ne.0) then 
         diffuave = dsum / REAL(nd, kind=PREC_DP) 
         diffusig = sqrt(dsum2 / REAL(nd, kind=PREC_DP) - diffuave**2) 
      ELSE 
         diffuave = 0.0 
         diffusig = 0.0 
      ENDIF 
!                                                                       
      IF (four_log) then 
         WRITE (output_io, 1000) 
         IF (braggmax.ge.0.0) then 
            WRITE (output_io, 1010) braggmin, braggmax 
         ENDIF 
         WRITE (output_io, 1020) diffumin, diffumax 
         WRITE (output_io, 1030) diffuave, diffusig 
      ENDIF 
!
      res_para(0) = 6
      res_para(1) = diffumin
      res_para(2) = diffumax
      res_para(3) = diffuave
      res_para(4) = diffusig
      res_para(5) = braggmin
      res_para(6) = braggmax
!                                                                       
 1000 FORMAT     (/,' ') 
 1010 FORMAT     (  ' Bragg scat.     : ',G13.6,'  -> ',G13.6) 
 1020 FORMAT     (  ' Diffuse scat.   : ',G13.6,'  -> ',G13.6) 
 1030 FORMAT     (  '      Average    : ',G13.6,'  +- ',G13.6) 
END SUBROUTINE four_qinfo                     
!
!**********************************************************************
!**********************************************************************
!
REAL(kind=PREC_DP) FUNCTION form (ll, scat, lxray, h2, power) 
!+                                                                      
!       calculates the form factor                                      
!-                                                                      
USE discus_config_mod
USE element_data_mod
USE precision_mod
IMPLICIT none 
!                                                                       
INTEGER, INTENT(IN) :: ll
LOGICAL, INTENT(IN) :: lxray 
REAL(kind=PREC_DP)   , DIMENSION(11,0:MAXSCAT), INTENT(INOUT) :: scat ! (11, 0:maxscat) 
REAL(kind=PREC_DP)              , INTENT(IN)    :: h2
INTEGER, INTENT(IN) :: power
!
INTEGER   :: i 
!                                                                       
form = scat (1, ll) 
IF (lxray) then 
   DO i = 1, power 
      form = form + scat (2 * i, ll) * exp ( - scat (2 * i + 1, ll)   * h2)
   ENDDO 
ENDIF 
!
END FUNCTION form                             
!
!**********************************************************************
!**********************************************************************
!
SUBROUTINE dlink (ano, lambda, rlambda, renergy, l_energy, &
                  diff_radiation, diff_table, diff_power) 
!-                                                                      
!     This routine reads wavelength symbols, wavelength values 
!     and atomic form factors from module "element_data_mod"
!     It replaces the old dlink and LAZY_PULVERIX
!+                                                                      
USE discus_config_mod 
USE charact_mod
USE crystal_mod 
USE element_data_mod
USE chem_aver_mod
USE param_mod
use precision_mod
!
IMPLICIT none 
!                                                                       
!                                                                       
LOGICAL             , INTENT(IN)   :: ano             ! Anomalous scattering TRUE/FALSE
CHARACTER (LEN = * ), INTENT(IN)   :: lambda          ! Wavelength symbol
REAL(kind=PREC_DP)  , INTENT(OUT)  :: rlambda         ! Wave length value
REAL(kind=PREC_DP)  , INTENT(INOUT):: renergy         ! Wave length energy
LOGICAL             , INTENT(IN)   :: l_energy        ! Wave length specified as energy TRUE/FALSE
INTEGER             , INTENT(IN)   :: diff_radiation  ! xray/neutron/electron
INTEGER             , INTENT(IN)   :: diff_table      ! International / Waasmaier
INTEGER             , INTENT(OUT)  :: diff_power      ! 4 or 5 parameters
!
INTEGER, PARAMETER :: RAD_XRAY = 1
INTEGER, PARAMETER :: RAD_NEUT = 2
INTEGER, PARAMETER :: RAD_ELEC = 3
INTEGER, PARAMETER :: RAD_INTER = 0
LOGICAL, PARAMETER :: LOUT = .FALSE.
!
CHARACTER (LEN = 4 ) :: element 
INTEGER    :: i
INTEGER    :: j
REAL(kind=PREC_DP)   , DIMENSION(1:11)  :: temp_scat  ! a1,b1,---a4,b4,c
REAL(kind=PREC_DP)   , DIMENSION(1:2)   :: temp_delf  ! delfr, delfi
REAL(kind=PREC_DP)                      :: temp_bcoh  ! b_choherent
!
!                                                                       
ier_num = -77 
ier_typ = ER_APPL 
!
call get_wave ( lambda, rlambda, renergy, l_energy,diff_radiation, &
                ier_num, ier_typ )
!
call chem_elem(lout)
!                                                                       
IF (ier_num.ne.0) RETURN 
!                                                                       
any_element: IF (cr_nscat.gt.0) THEN 
   DO i = 1, cr_nscat 
      IF (cr_at_lis (i) /= 'XAXI' .AND. cr_at_lis (i) /= 'YAXI' .and. & 
          cr_at_lis (i) /= 'ZAXI') THEN                  
         IF(res_para(i+1)> 0.0) THEN
            IF (cr_scat_int (i) ) then 
               SELECTCASE(diff_radiation)
               CASE(RAD_NEUT)        !  neutron scattering
                  IF (cr_scat_equ (i) ) then 
                     element =         cr_at_equ (i) 
                  ELSE 
                     element =         cr_at_lis (i) 
                  ENDIF 
!                                                                       
                  call symbf ( element, j)
                  IF ( j /= 0 ) THEN 
                     call get_scat_neut ( j, temp_bcoh )
                     cr_scat(:,i) = 0.0
                     cr_scat(1,i) = temp_bcoh
                     cr_delfr (i) = 0.0 
                     cr_delfi (i) = 0.0 
                     diff_power   = PER_RAD_POWER(diff_radiation, 0)
                     ier_num = 0
                     ier_typ = ER_NONE 
                  ELSE
                     ier_num = -20 
                     ier_typ = ER_APPL 
                  ENDIF
               CASE(RAD_XRAY)        !  Xray diffraction
                  IF (cr_scat_equ (i) ) then 
                     element =         cr_at_equ (i) 
                  ELSE 
                     element =         cr_at_lis (i) 
                  ENDIF 
!
                  ier_num = -20 
                  ier_typ = ER_APPL 
!                                                                       
                  call symbf ( element, j)
                  IF ( j /= 0 ) THEN 
                     if(diff_table==RAD_INTER) then
                        call get_scat_xray ( j, temp_scat )
                     else
                        call get_scat_xray_waas(j, temp_scat)
                     endif
                     cr_scat(:,i) = 0.0
                     cr_scat(:,i) = temp_scat(:)   ! copy temp array into 1st column
!                    cr_delfr (i) = 0.0   ! DEVELOPMENT !
!                    cr_delfi (i) = 0.0   ! DEVELOPMENT !
                     IF(ano) THEN
                        IF(cr_delf_int (i) ) then 
                           call get_scat_ano ( j, lambda, temp_delf )
                           cr_delfr (i) = temp_delf(1)
                           cr_delfi (i) = temp_delf(2)
                        ELSE
                           cr_delfr(i) =   cr_delfr_u(i)
                           cr_delfi(i) =   cr_delfi_u(i)
                        ENDIF 
                     ELSE
                        cr_delfr(i) = 0.0   ! DEVELOPMENT !
                        cr_delfi(i) = 0.0   ! DEVELOPMENT !
                     ENDIF 
                     diff_power   = PER_RAD_POWER(diff_radiation, diff_table)
                     ier_num = 0
                  ELSE
                     ier_typ = ER_NONE 
                     ier_typ = ER_APPL 
                  ENDIF
!
               CASE(RAD_ELEC)        !  Electron diffraction
                  IF (cr_scat_equ (i) ) then 
                     element =         cr_at_equ (i) 
                  ELSE 
                     element =         cr_at_lis (i) 
                  ENDIF 
!
                  ier_num = -20 
                  ier_typ = ER_APPL 
!                                                                       
                  call symbf ( element, j)
                  IF ( j /= 0 ) THEN 
                     call get_scat_elec ( j, temp_scat)
                     cr_scat(:,i) = 0.0
                     cr_scat(:,i) = temp_scat(:)   ! copy temp array into 1st column
                     cr_delfr (i) = 0.0
                     cr_delfi (i) = 0.0
                     diff_power   = PER_RAD_POWER(diff_radiation, 0)
                     ier_num = 0
                  ELSE
                     ier_typ = ER_NONE 
                     ier_typ = ER_APPL 
                  ENDIF
                  j = LEN_TRIM(element)
                  IF(MAXVAL(cr_scat(:,i))==0.0 .AND. &
                     (element(j:j)=='+' .OR. element(j:j)=='-')) THEN
                     ier_num = -173
                     ier_typ = ER_APPL
                     WRITE(ier_msg(1),'(3a,i4,a)') 'Element ',element(1:j), '(',i,')'
                     WRITE(ier_msg(2),'(a)') 'Replace by neutral atom or use in fourier'
                     WRITE(ier_msg(3),'(a)') 'scat <ion_name>, <neutron_name>'
                     RETURN
                  ENDIF
            END SELECT
!                                                                       
            ENDIF 
         ELSE
            cr_scat(:,i) = 0.00
            cr_delfr (i) = 0.0
            cr_delfi (i) = 0.0
         ENDIF
      ENDIF 
   ENDDO 
ELSE any_element
   ier_num = - 21 
   ier_typ = ER_APPL 
ENDIF any_element
!                                                                       
END SUBROUTINE dlink                          
!
!************************************************************************************
! Neder's original code
!************************************************************************************
!
! SUBROUTINE calc_000 (rhkl) 
! !-                                                                      
! !     Calculates the value of F(rh,rk,rl)                               
! !+                                                                      
!       USE discus_config_mod 
!       USE crystal_mod 
!       USE diffuse_mod 
! !
!       USE param_mod 
!       USE precision_mod
! !
!       IMPLICIT none 
! !                                                                       
       
! !                                                                       
!       REAL(kind=PREC_DP), DIMENSION(3), INTENT(IN) :: rhkl
! !
!       INTEGER i, j 
!       INTEGER shel_inc (3) 
!       REAL(kind=PREC_DP) :: shel_eck (3, 4) 
!       REAL(kind=PREC_DP) :: shel_vi (3, 3) 
!       COMPLEX (KIND=PREC_DP) :: shel_acsf                ! Not clear on the meaning of this lines (shel_acsf, shel_dsi, shel_tcsf)
!       REAL    (KIND=PREC_DP) :: shel_dsi 
!       COMPLEX (KIND=PREC_DP) :: shel_tcsf 
! !                                                                       
!       DO i = 1, 3 
!          shel_inc (i) = inc (i) 
!       ENDDO 
! !
!       DO i = 1, 3 
!          DO j = 1, 4 
!             shel_eck (i, j) = eck (i, j) 
!          ENDDO 
!          DO j = 1, 3 
!             shel_vi (i, j) = vi (i, j) 
!          ENDDO 
!       ENDDO 
! !
!       shel_tcsf = csf (1)                              ! Relation between shel_tcsf and shel_acsf to csf and acsf is not clear.
!       shel_acsf = acsf (1) 
!       shel_dsi  = dsi (1) 
!       inc (1) = 1 
!       inc (2) = 1 
!       inc (3) = 1 
!       DO i = 1, 3 
!          DO j = 1, 4 
!             eck (i, j) = rhkl (i) 
!          ENDDO 
!          DO j = 1, 3 
!             vi (i, j) = 0.0 
!          ENDDO 
!       ENDDO 
! !
!       four_log = .false. 
!       call four_run 
!       res_para (1) =      REAL(csf (1) , KIND=KIND(0.0D0) ) 
!       res_para (2) = REAL(AIMAG(csf (1)), KIND=KIND(0.0D0) ) 
!       res_para (3) = res_para (1) / cr_icc (1) / cr_icc (2) / cr_icc (3) 
!       res_para (4) = res_para (2) / cr_icc (1) / cr_icc (2) / cr_icc (3) 
!       res_para (0) = 4 
!       csf (1) = shel_tcsf                             ! Relation between shel_tcsf and shel_acsf to csf and acsf is not clear.
!       acsf(1) = shel_acsf 
!       dsi (1) = shel_dsi 
! !
!       DO i = 1, 3 
!          inc (i) = shel_inc (i) 
!       ENDDO 
!       DO i = 1, 3 
!          DO j = 1, 4 
!             eck (i, j) = shel_eck (i, j) 
!          ENDDO 
!          DO j = 1, 3 
!             vi (i, j) = shel_vi (i, j) 
!          ENDDO 
!       ENDDO 
! !                                                                       
! END SUBROUTINE calc_000                       
!
!************************************************************************************
! My code ( Copilot Generated )
!************************************************************************************

SUBROUTINE calc_000 (rhkl) 
!-                                                                      
!     Calculates the value of F(rh,rk,rl)                               
!+                                                                      
USE discus_config_mod 
USE crystal_mod 
USE diffuse_mod 
!
USE param_mod 
USE precision_mod
!
IMPLICIT none 
!                                                                       
       
!                                                                       
REAL(kind=PREC_DP), DIMENSION(3), INTENT(IN) :: rhkl
!
INTEGER i, j 
INTEGER shel_inc (3) 
REAL(kind=PREC_DP) :: shel_eck (3, 4) 
REAL(kind=PREC_DP) :: shel_vi (3, 3) 
COMPLEX (KIND=PREC_DP) :: shel_acsf                ! treat acsf as 3D array
REAL    (KIND=PREC_DP) :: shel_dsi                  ! treat dsi as 3D array
COMPLEX (KIND=PREC_DP) :: shel_tcsf                 ! treat tcsf as 3D array
!                                                                       
DO i = 1, 3 
   shel_inc (i) = inc (i) 
ENDDO 
!
DO i = 1, 3 
   DO j = 1, 4 
      shel_eck (i, j) = eck (i, j) 
   ENDDO 
   DO j = 1, 3 
      shel_vi (i, j) = vi (i, j) 
   ENDDO 
ENDDO 
!
shel_tcsf = csf (1,1,1)
shel_acsf = acsf(1,1,1) 
shel_dsi  = dsi (1,1,1) 
inc (1) = 1 
inc (2) = 1 
inc (3) = 1 
DO i = 1, 3 
   DO j = 1, 4 
      eck (i, j) = rhkl (i) 
   ENDDO 
   DO j = 1, 3 
      vi (i, j) = 0.0 
   ENDDO 
ENDDO 
!
four_log = .false. 
call four_run 
res_para (1) =      REAL(csf (1,1,1) , KIND=KIND(0.0D0) ) 
res_para (2) = REAL(AIMAG(csf (1,1,1)), KIND=KIND(0.0D0) ) 
res_para (3) = res_para (1) / cr_icc (1) / cr_icc (2) / cr_icc (3) 
res_para (4) = res_para (2) / cr_icc (1) / cr_icc (2) / cr_icc (3) 
res_para (0) = 4 
csf (1,1,1) = shel_tcsf
acsf(1,1,1) = shel_acsf
dsi (1,1,1) = shel_dsi 
!
DO i = 1, 3 
   inc (i) = shel_inc (i) 
ENDDO 
DO i = 1, 3 
   DO j = 1, 4 
      eck (i, j) = shel_eck (i, j) 
   ENDDO 
   DO j = 1, 3 
      vi (i, j) = shel_vi (i, j) 
   ENDDO 
ENDDO 
!                                                                       
END SUBROUTINE calc_000    
!
!************************************************************************************
!************************************************************************************
! Neders original code
! 
! SUBROUTINE calc_hkl(infile,infile_l, calcfile, calcfile_l,scale,style)
! !
!       USE crystal_mod 
!       USE diffuse_mod 
!       USE discus_allocate_appl_mod
!       USE get_params_mod
!       USE param_mod
! USE prompt_mod
! USE precision_mod
! USE support_mod
! !
!       IMPLICIT NONE
! !
!       CHARACTER(LEN=*), INTENT(IN) :: infile
!       INTEGER         , INTENT(IN) :: infile_l
!       CHARACTER(LEN=*), INTENT(IN) :: calcfile
!       INTEGER         , INTENT(IN) :: calcfile_l
!       REAL(KIND=PREC_DP),INTENT(IN) :: scale
!       INTEGER         , INTENT(IN) :: style
! !
!       INTEGER, PARAMETER :: ird = 54
!       INTEGER, PARAMETER :: iwr = 55
!       INTEGER, PARAMETER :: HKLF4 = 4
!       INTEGER, PARAMETER :: CIF   = 1
! !
!       CHARACTER(LEN=PREC_STRING) :: line
!       CHARACTER(LEN=PREC_STRING), DIMENSION(:), ALLOCATABLE   :: ccpara
!       INTEGER            , DIMENSION(:), ALLOCATABLE   :: llpara
!       INTEGER            :: n_qxy, n_natoms,n_nscat
!       INTEGER            :: iostatus
!       INTEGER            :: ih,ik,il
!       INTEGER            :: ih_min,ik_min,il_min
!       INTEGER            :: ih_max,ik_max,il_max
!       INTEGER            :: n_refl
!       INTEGER            :: indx
!       INTEGER            :: startline
!       INTEGER            :: j
!       INTEGER            :: nentries
!       INTEGER            :: length, ianz
!       INTEGER            ::   j_h      = 0
!       INTEGER            ::   j_k      = 0
!       INTEGER            ::   j_l      = 0
!       INTEGER            ::   j_iobs   = 0
!       INTEGER            ::   j_icalc  = 0
!       INTEGER            ::   j_sigi   = 0
!       INTEGER            ::   j_fobs   = 0
!       INTEGER            ::   j_fcalc  = 0
!       INTEGER            ::   j_sigf   = 0
!       INTEGER            ::   j_flag   = 0
!       REAL(kind=PREC_DP)               :: rint, sint, wert
!       REAL(kind=PREC_DP), DIMENSION(7) :: values
!       REAL(kind=PREC_DP), DIMENSION(3) :: rhkl
! !
!       n_qxy    = 1
!       n_natoms = 1
!       n_nscat  = 1
! !
!       call oeffne(ird, infile(1:infile_l),   'old') 
!       call oeffne(iwr, calcfile(1:calcfile_l), 'unknown') 
!       inc(:) = 1
! !
!       ih_min  = 0
!       ih_max  = 0
!       ik_min  = 0
!       ik_max  = 0
!       il_min  = 0
!       il_max  = 0
!       n_refl  = 0
!       startline = 0
!       IF(style==HKLF4) THEN
!          startline = 1
! check:   DO     ! First read, get size and extrema
!             READ(ird,'(a)', IOSTAT=iostatus) line
!             IF(IS_IOSTAT_END(iostatus)) EXIT check
!             IF(line(1:1)=='#' .OR. line(1:1)=='!' .OR. line(1:1)=='d') THEN
!                ier_num = -9999
!                ier_typ = ER_APPL
!                CLOSE(ird) 
!                RETURN
!             ENDIF
!             READ(line,1000, IOSTAT=iostatus) ih,ik,il, rint, sint
!             ih_min = MIN(ih_min,ih)
!             ih_max = MAX(ih_max,ih)
!             ik_min = MIN(ik_min,ik)
!             ik_max = MAX(ik_max,ik)
!             il_min = MIN(il_min,il)
!             il_max = MAX(il_max,il)
!             n_refl = n_refl + 1
!          ENDDO check
!       ELSEIF(style==CIF) THEN
! find:    DO
!             READ(ird,'(a)', IOSTAT=iostatus) line
!             IF(IS_IOSTAT_END(iostatus)) THEN
!                ier_num = -9999
!                ier_typ = ER_APPL
!                CLOSE(ird) 
!                RETURN
!             ENDIF
!             startline = startline + 1
!             IF(line(1:14) == '_refln_index_h') EXIT find
!          ENDDO find
!          nentries = 0
!          j_h      = 0
!          j_k      = 0
!          j_l      = 0
!          j_iobs   = 0
!          j_icalc  = 0
!          j_sigi   = 0
!          j_fobs   = 0
!          j_fcalc  = 0
!          j_sigf   = 0
!          j_flag   = 0
!          j = 1
! entries: DO
!             IF(line(1:1)=='#' .or. line == ' ') THEN
!             startline = startline + 1
!                CYCLE entries
!             ENDIF
!             IF(line(1:14) == '_refln_index_h')         j_h = j
!             IF(line(1:14) == '_refln_index_k')         j_k = j
!             IF(line(1:14) == '_refln_index_l')         j_l = j
!             IF(line(1:21) == '_refln_F_squared_calc')  j_icalc = j
!             IF(line(1:21) == '_refln_F_squared_meas')  j_iobs  = j
!             IF(line(1:22) == '_refln_F_squared_sigma') j_sigi  = j
!             IF(line(1:22) == '_refln_observed_status') j_flag  = j
!             IF(line(1:7) /= '_refln_') EXIT entries
!             nentries = nentries + 1
!             READ(ird,'(a)', IOSTAT=iostatus) line
!             IF(IS_IOSTAT_END(iostatus)) THEN
!                ier_num = -9999
!                ier_typ = ER_APPL
!                CLOSE(ird) 
!                RETURN
!             ENDIF
!             j = j + 1
!             startline = startline + 1
!          ENDDO entries
!          ALLOCATE(ccpara(1:nentries))
!          ccpara(:) = ' '
!          ALLOCATE(llpara(1:nentries))
!          length = LEN_TRIM(line)
! check2:  DO
!             call get_params_blank(line,ianz, ccpara,llpara, nentries, length)
!             READ(ccpara(j_h)(1:llpara(j_h)),*) ih 
!             READ(ccpara(j_k)(1:llpara(j_k)),*) ik 
!             READ(ccpara(j_l)(1:llpara(j_l)),*) il 
!             ih_min = MIN(ih_min,ih)
!             ih_max = MAX(ih_max,ih)
!             ik_min = MIN(ik_min,ik)
!             ik_max = MAX(ik_max,ik)
!             il_min = MIN(il_min,il)
!             il_max = MAX(il_max,il)
!             n_refl = n_refl + 1
!             READ(ird,'(a)', IOSTAT=iostatus) line
! !        READ(ird,*   , IOSTAT=iostatus) ih,ik,il, rint, sint
!             IF(IS_IOSTAT_END(iostatus)) EXIT check2
!             length = LEN_TRIM(line)
!          ENDDO check2
!       ELSE
!          WRITE(output_io,*) ' WRONG STYLE ', style
!          RETURN
!       ENDIF
!       vi(:,:) = 0
!       vi(1,1) = 1.0
!       vi(2,2) = 1.0
!       vi(3,3) = 1.0
!       inc(1)   = ih_max - ih_min             + 1
!       inc(2)   = ik_max - ik_min             + 1
!       inc(3)   = il_max - il_min             + 1
!       eck(1,1) = ih_min             ! minimum H
!       eck(2,1) = ik_min             ! minimum K
!       eck(3,1) = il_min             ! minimum L
!       eck(1,2) = ih_max             ! maximum H
!       eck(2,2) = ik_min             ! minimum K
!       eck(3,2) = il_min             ! minimum L
!       eck(1,3) = ih_min             ! minimum H
!       eck(2,3) = ik_max             ! maximum K
!       eck(3,3) = il_min             ! minimum L
!       eck(1,4) = ih_min             ! minimum H
!       eck(2,4) = ik_min             ! minimum K
!       eck(3,4) = il_max             ! maximum L
! !
!       IF (ier_num == 0) then 
!          IF (inc(1) * inc(2) * inc(3) .gt. MAXQXY  .OR.   &
!              cr_natoms > DIF_MAXAT                 .OR.   &
!              cr_nscat>DIF_MAXSCAT              ) THEN
!             n_qxy    = MAX(n_qxy,inc(1) * inc(2)*inc(3),MAXQXY)
!             n_natoms = MAX(n_natoms,cr_natoms,DIF_MAXAT)
!             n_nscat  = MAX(n_nscat,cr_nscat,DIF_MAXSCAT)
!             call alloc_diffuse (n_qxy, n_nscat, n_natoms)
!             IF (ier_num /= 0) THEN
!                RETURN
!             ENDIF
!          ENDIF
!          call dlink (ano, lambda, rlambda, renergy, l_energy, &
!                      diff_radiation, diff_table, diff_power) 
!          call four_run
!          rhkl (1) =  0. 
!          rhkl (2) =  0. 
!          rhkl (3) =  0. 
!          REWIND(ird)
!          DO j=1,startline -1 
!             READ(ird,'(a)') line
!          ENDDO
! main:    DO
!             IF(style==HKLF4) THEN
!                READ(ird,1000, IOSTAT=iostatus) ih,ik,il, rint, sint
!             ELSEIF(style==CIF) THEN
!                READ(ird,'(a)',IOSTAT=iostatus) line
!                length = LEN_TRIM(line)
!                call get_params_blank(line,ianz, ccpara,llpara, nentries, length)
!                DO j=1,nentries-1
!                   READ(ccpara(j)(1:llpara(j)),*) values(j)
!                ENDDO
!                ih = NINT(values(1))
!                ik = NINT(values(2))
!                il = NINT(values(3))
! !              READ(ird,*   , IOSTAT=iostatus) ih,ik,il, rint, sint
!             ENDIF
!             IF(IS_IOSTAT_END(iostatus)) EXIT main
!             indx = (ih-ih_min)*inc(3)*inc(2) + (ik-ik_min)*inc(3) + (il-il_min)  + 1       ! This entire portion looks complicated. The dimensions here are going to be tricky
!             wert = REAL(csf(indx)*CONJG(csf(indx)),KIND=KIND(0.0D0))                       ! Discuss with Neder!!!!
!             sint = SQRT(ABS(wert))
!             IF(ih==0 .AND. ik==0 .AND. IL==0) THEN
!                WRITE(iwr,1000) ih,ik,il, 0.00, 0.00
!             ELSE
!                IF(style==HKLF4) THEN
!                   WRITE(iwr,1000) ih,ik,il, scale*wert, sqrt(scale)*sint 
!                ELSEIF(style==CIF) THEN
!                   IF(j_icalc/=0) THEN
!                      values(j_icalc) = scale*wert
!                   ENDIF
!                   WRITE(iwr, 2000) ih,ik,il, (values(j),j=4, nentries-1)
!                ENDIF
!             ENDIF
!          END DO main
!       ENDIF 
!       IF(ALLOCATED(ccpara))   DEALLOCATE(ccpara)
!       IF(ALLOCATED(llpara))   DEALLOCATE(llpara)
! 1000  FORMAT(3I4,F8.2,F8.2)
! 2000  FORMAT(3I4,F12.2,F12.2, F12.2)
!       CLOSE(ird)
!       CLOSE(iwr)
! END SUBROUTINE calc_hkl
!
!************************************************************************************
! My code ( Copilot Generated )
!************************************************************************************
SUBROUTINE calc_hkl(infile,infile_l, calcfile, calcfile_l,scale,style)
!
USE crystal_mod 
USE diffuse_mod 
USE discus_allocate_appl_mod
USE get_params_mod
USE param_mod
USE prompt_mod
USE precision_mod
USE support_mod
!
IMPLICIT NONE
!
CHARACTER(LEN=*), INTENT(IN) :: infile
INTEGER         , INTENT(IN) :: infile_l
CHARACTER(LEN=*), INTENT(IN) :: calcfile
INTEGER         , INTENT(IN) :: calcfile_l
REAL(KIND=PREC_DP),INTENT(IN) :: scale
INTEGER         , INTENT(IN) :: style
!
INTEGER, PARAMETER :: ird = 54
INTEGER, PARAMETER :: iwr = 55
INTEGER, PARAMETER :: HKLF4 = 4
INTEGER, PARAMETER :: CIF   = 1
!
CHARACTER(LEN=PREC_STRING) :: line
CHARACTER(LEN=PREC_STRING), DIMENSION(:), ALLOCATABLE   :: ccpara
INTEGER            , DIMENSION(:), ALLOCATABLE   :: llpara
INTEGER, dimension(3) :: n_qxy
INTEGER            :: n_natoms,n_nscat
INTEGER            :: iostatus
INTEGER            :: ih,ik,il
INTEGER            :: ih_min,ik_min,il_min
INTEGER            :: ih_max,ik_max,il_max
INTEGER            :: n_refl
INTEGER            :: indx
INTEGER            :: startline
INTEGER            :: j
INTEGER            :: nentries
INTEGER            :: length, ianz
INTEGER            ::   j_h      = 0
INTEGER            ::   j_k      = 0
INTEGER            ::   j_l      = 0
INTEGER            ::   j_iobs   = 0
INTEGER            ::   j_icalc  = 0
INTEGER            ::   j_sigi   = 0
INTEGER            ::   j_fobs   = 0
INTEGER            ::   j_fcalc  = 0
INTEGER            ::   j_sigf   = 0
INTEGER            ::   j_flag   = 0
REAL(kind=PREC_DP)               :: rint, sint, wert
REAL(kind=PREC_DP), DIMENSION(7) :: values
REAL(kind=PREC_DP), DIMENSION(3) :: rhkl
!
n_qxy    = 1
n_natoms = 1
n_nscat  = 1
!
call oeffne(ird, infile(1:infile_l),   'old') 
call oeffne(iwr, calcfile(1:calcfile_l), 'unknown') 
inc(:) = 1
!
ih_min  = 0
ih_max  = 0
ik_min  = 0
ik_max  = 0
il_min  = 0
il_max  = 0
n_refl  = 0
startline = 0
IF(style==HKLF4) THEN
   startline = 1
check:   DO     ! First read, get size and extrema
      READ(ird,'(a)', IOSTAT=iostatus) line
      IF(IS_IOSTAT_END(iostatus)) EXIT check
      IF(line(1:1)=='#' .OR. line(1:1)=='!' .OR. line(1:1)=='d') THEN
         ier_num = -9999
         ier_typ = ER_APPL
         CLOSE(ird) 
         RETURN
      ENDIF
      READ(line,1000, IOSTAT=iostatus) ih,ik,il, rint, sint
      ih_min = MIN(ih_min,ih)
      ih_max = MAX(ih_max,ih)
      ik_min = MIN(ik_min,ik)
      ik_max = MAX(ik_max,ik)
      il_min = MIN(il_min,il)
      il_max = MAX(il_max,il)
      n_refl = n_refl + 1
   ENDDO check
ELSEIF(style==CIF) THEN
find:    DO
      READ(ird,'(a)', IOSTAT=iostatus) line
      IF(IS_IOSTAT_END(iostatus)) THEN
         ier_num = -9999
         ier_typ = ER_APPL
         CLOSE(ird) 
         RETURN
      ENDIF
      startline = startline + 1
      IF(line(1:14) == '_refln_index_h') EXIT find
   ENDDO find
   nentries = 0
   j_h      = 0
   j_k      = 0
   j_l      = 0
   j_iobs   = 0
   j_icalc  = 0
   j_sigi   = 0
      j_fobs   = 0
      j_fcalc  = 0
      j_sigf   = 0
      j_flag   = 0
      j = 1
   entries: DO
      IF(line(1:1)=='#' .or. line == ' ') THEN
            startline = startline + 1
         CYCLE entries
      ENDIF
      IF(line(1:14) == '_refln_index_h')         j_h = j
      IF(line(1:14) == '_refln_index_k')         j_k = j
      IF(line(1:14) == '_refln_index_l')         j_l = j
      IF(line(1:21) == '_refln_F_squared_calc')  j_icalc = j
      IF(line(1:21) == '_refln_F_squared_meas')  j_iobs  = j
      IF(line(1:22) == '_refln_F_squared_sigma') j_sigi  = j
      IF(line(1:22) == '_refln_observed_status') j_flag  = j
      IF(line(1:7) /= '_refln_') EXIT entries
      nentries = nentries + 1
      READ(ird,'(a)', IOSTAT=iostatus) line
      IF(IS_IOSTAT_END(iostatus)) THEN
         ier_num = -9999
         ier_typ = ER_APPL
         CLOSE(ird) 
         RETURN
      ENDIF
      j = j + 1
      startline = startline + 1
   ENDDO entries
   ALLOCATE(ccpara(1:nentries))
   ccpara(:) = ' '
   ALLOCATE(llpara(1:nentries))
   length = LEN_TRIM(line)
   check2:  DO
      call get_params_blank(line,ianz, ccpara,llpara, nentries, length)
      READ(ccpara(j_h)(1:llpara(j_h)),*) ih 
      READ(ccpara(j_k)(1:llpara(j_k)),*) ik 
      READ(ccpara(j_l)(1:llpara(j_l)),*) il 
      ih_min = MIN(ih_min,ih)
      ih_max = MAX(ih_max,ih)
      ik_min = MIN(ik_min,ik)
      ik_max = MAX(ik_max,ik)
      il_min = MIN(il_min,il)
      il_max = MAX(il_max,il)
      n_refl = n_refl + 1
      READ(ird,'(a)', IOSTAT=iostatus) line
         !        READ(ird,*   , IOSTAT=iostatus) ih,ik,il, rint, sint
      IF(IS_IOSTAT_END(iostatus)) EXIT check2
      length = LEN_TRIM(line)
   ENDDO check2
ELSE
   WRITE(output_io,*) ' WRONG STYLE ', style
   RETURN
ENDIF
vi(:,:) = 0
vi(1,1) = 1.0
vi(2,2) = 1.0
vi(3,3) = 1.0
inc(1)   = ih_max - ih_min             + 1
inc(2)   = ik_max - ik_min             + 1
inc(3)   = il_max - il_min             + 1
eck(1,1) = ih_min             ! minimum H
eck(2,1) = ik_min             ! minimum K
eck(3,1) = il_min             ! minimum L
eck(1,2) = ih_max             ! maximum H
eck(2,2) = ik_min             ! minimum K
eck(3,2) = il_min             ! minimum L
eck(1,3) = ih_min             ! minimum H
eck(2,3) = ik_max             ! maximum K
eck(3,3) = il_min             ! minimum L
eck(1,4) = ih_min             ! minimum H
eck(2,4) = ik_min             ! minimum K
eck(3,4) = il_max             ! maximum L
         !
IF (ier_num == 0) then 
!                 IF (inc(1) * inc(2) * inc(3) .gt. MAXQXY  .OR.   &
   if(any(num/=ubound(csf))) then
      n_qxy = num
      call alloc_diffuse_four (n_qxy )
      if(ier_num/=0) return
   endif
   if(cr_nscat/=ubound(cfact,2)) then
      call alloc_diffuse_scat(cr_nscat)
      if(ier_num/=0) return
   endif
   if(cr_natoms/=ubound(xat,1)) then
      call alloc_diffuse_atom(cr_natoms)
      if(ier_num/=0) return
   endif
!   if(inc(1)>MAXQXY(1) .or. inc(2)>MAXQXY(2) .or. inc(3)>MAXQXY(3) .or. &
!      cr_natoms > DIF_MAXAT                 .OR.   &
!      cr_nscat>DIF_MAXSCAT              ) THEN
!!     n_qxy    = MAX(n_qxy,inc,MAXQXY)
!      n_natoms = MAX(n_natoms,cr_natoms,DIF_MAXAT)
!     n_nscat  = MAX(n_nscat,cr_nscat,DIF_MAXSCAT)
!     call alloc_diffuse_four (n_qxy)
!     call alloc_diffuse_scat (n_nscat)
!      call alloc_diffuse_atom (n_natoms)
!      IF (ier_num /= 0) THEN
!         RETURN
!      ENDIF
!   ENDIF
   call dlink (ano, lambda, rlambda, renergy, l_energy, &
               diff_radiation, diff_table, diff_power) 
   call four_run
   rhkl (1) =  0. 
   rhkl (2) =  0. 
   rhkl (3) =  0. 
   REWIND(ird)
   DO j=1,startline -1 
      READ(ird,'(a)') line
   ENDDO
   main:    DO
      IF(style==HKLF4) THEN
         READ(ird,1000, IOSTAT=iostatus) ih,ik,il, rint, sint
      ELSEIF(style==CIF) THEN
         READ(ird,'(a)',IOSTAT=iostatus) line
         length = LEN_TRIM(line)
         call get_params_blank(line,ianz, ccpara,llpara, nentries, length)
         DO j=1,nentries-1
            READ(ccpara(j)(1:llpara(j)),*) values(j)
         ENDDO
         ih = NINT(values(1))
         ik = NINT(values(2))
         il = NINT(values(3))
         !              READ(ird,*   , IOSTAT=iostatus) ih,ik,il, rint, sint
      ENDIF
      IF(IS_IOSTAT_END(iostatus)) EXIT main
      indx = (ih-ih_min)*inc(3)*inc(2) + (ik-ik_min)*inc(3) + (il-il_min)  + 1       
      wert = REAL(csf(ih,ik,il)*CONJG(csf(ih,ik,il)),KIND=KIND(0.0D0))                      
      sint = SQRT(ABS(wert))
      IF(ih==0 .AND. ik==0 .AND. IL==0) THEN
         WRITE(iwr,1000) ih,ik,il, 0.00, 0.00
      ELSE
         IF(style==HKLF4) THEN
           WRITE(iwr,1000) ih,ik,il, scale*wert, sqrt(scale)*sint 
         ELSEIF(style==CIF) THEN
           IF(j_icalc/=0) THEN
             values(j_icalc) = scale*wert
           ENDIF
           WRITE(iwr, 2000) ih,ik,il, (values(j),j=4, nentries-1)
         ENDIF
      ENDIF
   END DO main
ENDIF 
IF(ALLOCATED(ccpara))   DEALLOCATE(ccpara)
IF(ALLOCATED(llpara))   DEALLOCATE(llpara)
1000  FORMAT(3I4,F8.2,F8.2)
2000  FORMAT(3I4,F12.2,F12.2, F12.2)
CLOSE(ird)
CLOSE(iwr)
!
END SUBROUTINE calc_hkl
!
!**********************************************************************
!**********************************************************************
!
SUBROUTINE four_strucf_aver (iscat, lform) 
!+                                                                      
!     Here the complex structure factor of 'nxat' identical atoms       
!     from array 'xat' is computed.                                     
!
!     The phase "iarg0" is calculated via integer math as offset from 
!     phase = 0 at hkl=0.
!-                                                                      
      USE discus_config_mod 
      USE diffuse_mod 
      USE precision_mod
      IMPLICIT none 
!                                                                       
      INTEGER, INTENT(IN) :: iscat 
      LOGICAL, INTENT(IN) :: lform 
!                                                                       
      REAL(PREC_DP)        ::        xincu, xincv , xincw
      REAL(PREC_DP)        ::        oincu, oincv , oincw
      INTEGER (KIND=PREC_INT_LARGE)   :: i, ii, j, k, l, iarg, iarg0, iincu, iincv, iincw
      INTEGER (KIND=PREC_INT_LARGE)   ::                              jincu, jincv, jincw
      INTEGER (KIND=PREC_INT_LARGE), PARAMETER :: shift = -6
!
      INTEGER IAND, ISHFT 
!
!------ zero fourier array                                              
!     
      ! Neder's original code
      !                                                                  
      !tcsf = cmplx (0.0D0, 0.0D0, KIND=KIND(0.0D0)) 
      !
      ! My declaration
      !
      DO i = 1, num (1)
         DO j = 1, num (2)
            DO k = 1, num (3)
               tcsf(i,j,k) = cmplx (0.0D0, 0.0D0, KIND=KIND(0.0D0))                                                                    ! I am using num(i) here.
            ENDDO
         ENDDO
      ENDDO
      !
!                                                                       
!------ Loop over all atoms in 'xat'                                    
!                                                                       
      DO k = 1, nxat 
!        xarg0 = xm (1)        * xat(k, 1) + xm (2)         * xat(k, 2) + xm (3)         * xat(k, 3)
         xincu = uin(1)        * xat(k, 1) + uin(2)         * xat(k, 2) + uin(3)         * xat(k, 3)
         xincv = vin(1)        * xat(k, 1) + vin(2)         * xat(k, 2) + vin(3)         * xat(k, 3)
         xincw = win(1)        * xat(k, 1) + win(2)         * xat(k, 2) + win(3)         * xat(k, 3)
         oincu = off_shift(1,1)* xat(k, 1) + off_shift(2,1) * xat(k, 2) + off_shift(3,1) * xat(k, 3)
         oincv = off_shift(1,2)* xat(k, 1) + off_shift(2,2) * xat(k, 2) + off_shift(3,2) * xat(k, 3)
         oincw = off_shift(1,3)* xat(k, 1) + off_shift(2,3) * xat(k, 2) + off_shift(3,3) * xat(k, 3)
!                                                                       
!        iarg0 = nint (64 * I2PI * (xarg0 - int (xarg0) + 0.0d0) ) 
         iincu = nint (64 * I2PI * (xincu - int (xincu) + 0.0d0) ) 
         iincv = nint (64 * I2PI * (xincv - int (xincv) + 0.0d0) ) 
         iincw = nint (64 * I2PI * (xincw - int (xincw) + 0.0d0) ) 
         jincu = nint (64 * I2PI * (oincu - int (oincu) + 0.0d0) ) 
         jincv = nint (64 * I2PI * (oincv - int (oincv) + 0.0d0) ) 
         jincw = nint (64 * I2PI * (oincw - int (oincw) + 0.0d0) ) 
         iarg0 =  lmn(1)*iincu + lmn(2)*iincv + lmn(3)*iincw + &
                  lmn(4)*jincu + lmn(5)*jincv + lmn(6)*jincw
         iarg = iarg0 
!                                                                       
!------ - Loop over all points in Q. 'iadd' is the address of the       
!------ - complex exponent table. 'IADD' divides out the 64 and         
!------ - ISHFT acts as MOD so that the argument stays in the table     
!------ - boundaries.                                                   
!                 iadd      = ISHFT (iarg, - 6) 
!                 iadd      = IAND  (iadd, MASK) 
!                 tcsf (ii) = tcsf (ii) + cex (iadd, MASK) )
!                 iarg      = iarg + iincw
!
!*********************************************************************************
! Neders original code
!**********************************************************************************                                                                       
!          ii = 0 
! !                                                                       
!           DO j = 0, num (1) - 1
!              DO i = 0, num (2) - 1
!                 iarg = iarg0 + iincu*j + iincv*i 
!                 DO h = 1, num (3) 
!                    ii       = ii + 1 
!                    tcsf(ii) = tcsf (ii) + cex (IAND  (ISHFT(iarg, shift), MASK) )                           ! Treat cearfully, dimensions here are complicated, and run differently for each array. 
!                    iarg     = iarg + iincw
!                 ENDDO 
!              ENDDO 
!           ENDDO 
!       ENDDO
!**********************************************************************************
! My code ( Copilot Generated )
!**********************************************************************************

         ii = 0 
   !                                                                       
         DO l = 1, num(3)
            DO j = 1, num(2)
               DO i = 1, num(1)
                  iarg = iarg0 + iincu*(j-1) + iincv*(i-1) 
                  tcsf(i,j,l) = tcsf(i,j,l) + cex(IAND(ISHFT(iarg, shift), MASK))
                  iarg = iarg + iincw
               ENDDO 
            ENDDO 
         ENDDO 
      ENDDO 
!**********************************************************************************
!
!------ Now we multiply with formfactor                                 
!                                                                       
      IF (lform) then 
         !DO  i = 1, num (1) * num (2) * num(3)
!        FORALL( i = 1: num (1) * num (2) * num(3))   !!! DO Loops seem to be faster!
            !tcsf (i) = tcsf (i) * cfact (istl (i), iscat)                                                   ! Neders original code
         DO i = 1, num (1)
            DO j = 1, num (2)
               DO l = 1, num (3)
                  tcsf (i,j,l) = tcsf (i,j,l) * cfact (istl (i,j,l), iscat)          ! My declaration
               ENDDO
            ENDDO
         ENDDO
!        END FORALL
         !END DO
      ENDIF 
!                                                                       
END SUBROUTINE four_strucf_aver
!
!
!*******************************************************************************
!
SUBROUTINE four_test_friedel
!
USE diffuse_mod
!
USE precision_mod
!
IMPLICIT NONE
!
REAL(PREC_DP), DIMENSION(3) :: rut   ! right upper top corner
REAL(PREC_DP), DIMENSION(3) :: dia   ! sum of all vi's
REAL(PREC_DP), DIMENSION(3) :: point ! sum of (left lower bottom) and (right upper top)
REAL(PREC_DP)               :: EPS = 1.E-6
!
diff_l_friedel = .FALSE.
IF(.NOT.ano) THEN        ! anomalous scattering is off, test space
   rut(1) = eck(1,1) + vi(1,1)*(inc(1)-1) + vi(1,2)*(inc(2)-1) + vi(1,3)*(inc(3)-1)
   rut(2) = eck(2,1) + vi(2,1)*(inc(1)-1) + vi(2,2)*(inc(2)-1) + vi(2,3)*(inc(3)-1)
   rut(3) = eck(3,1) + vi(3,1)*(inc(1)-1) + vi(3,2)*(inc(2)-1) + vi(3,3)*(inc(3)-1)
!
   diff_l_even(1) = MOD(inc(1),2)==0
   diff_l_even(2) = MOD(inc(2),2)==0
   diff_l_even(3) = MOD(inc(3),2)==0
   point(1) = eck(1,1) + rut(1)
   point(2) = eck(2,1) + rut(2)
   point(3) = eck(3,1) + rut(3)
   dia(1) = vi(1,1) + vi(1,2) + vi(1,3)
   dia(2) = vi(2,1) + vi(2,2) + vi(3,3)
   dia(3) = vi(3,1) + vi(3,2) + vi(3,3)
   IF(.NOT.diff_l_even(1) .AND. .NOT.diff_l_even(2) .AND. .NOT.diff_l_even(3)) THEN   ! All ODD increments
      diff_l_even(0) = .FALSE.    ! All odd
      IF(ABS(point(1))<EPS .AND. ABS(point(3))<EPS .AND.                           &
         ABS(point(2))<EPS                               ) THEN
!                         !left lower bottom and right upper top add up to zero
         diff_eck_u = eck         ! Back up user settings
         diff_vi_u  = vi          ! Back up user settings
         diff_inc_u = inc         ! Back up user settings
         IF(inc(1) > 1) THEN      ! We have an abscissa vector, cut this one
            inc(1) = (inc(1) + 1)/2
            diff_l_friedel = .TRUE.    ! A suitable solution
            diff_idim      = 1
         ELSEIF(inc(1)==1 .AND. inc(2)>1) THEN
            inc(2) = (inc(2) + 1)/2
            diff_l_friedel = .TRUE.    ! A suitable solution
            diff_idim      = 2
         ELSEIF(inc(1)==1 .AND. inc(2)==1 .AND. inc(3)>1) THEN
            inc(3) = (inc(3) + 1)/2
            diff_l_friedel = .TRUE.    ! A suitable solution
            diff_idim      = 3
         ELSE
            diff_l_friedel = .FALSE.   ! Not a suitable solution
            diff_idim      = 0
         ENDIF
      ENDIF
!  ELSEIF(diff_l_even(1) .AND. diff_l_even(2) .AND. diff_l_even(3)) THEN              ! All Even increments
!     diff_l_even(0) = .TRUE.    ! All even
!     IF(ABS(point(1)-dia(1))<EPS .AND. ABS(point(3)-dia(2))<EPS .AND.                           &
!        ABS(point(2)-dia(3))<EPS                                     ) THEN
!                         !left lower bottom and right upper top add up to sum of vi
!        diff_eck_u = eck         ! Back up user settings
!        diff_vi_u  = vi          ! Back up user settings
!        diff_inc_u = inc         ! Back up user settings
!        IF(inc(1) > 1) THEN      ! We have an abscissa vector, cut this one
!           inc(1) = inc(1)/2 + 1
!        ELSEIF(inc(1)==1 .AND. inc(2) > 1) THEN      ! We have an abscissa vector, cut this one
!           inc(2) = inc(2)/2 + 1
!        ENDIF
!     ENDIF
   ELSE
           diff_l_friedel = .FALSE.   ! Not a suitable solution
   ENDIF
!write(*,*) ' llb ', eck(:,1)
!write(*,*) ' rut ', rut(:)
!write(*,*) ' point',point(:)
!write(*,*) ' RES ', diff_l_friedel
!write(*,*) ' INC ', inc(:)
ENDIF 
END SUBROUTINE four_test_friedel
!
!*******************************************************************************
! Neders original code
!*******************************************************************************

! SUBROUTINE four_apply_friedel
! !-
! !  Applies Friedels law to the Fourier calculation
! !+
! USE diffuse_mod
! !
! IMPLICIT NONE
! !
! INTEGER :: i, j
! !
! IF(diff_l_friedel) THEN
!    j = diff_inc_u(1)*diff_inc_u(2)*diff_inc_u(3) + 1                             ! This entire portion looks complicated. The dimensions here are going to be tricky. Who is j ??
!    DO i = 1, num(1) * num(2) * num(3)
!      acsf(j-i) = conjg(acsf(i))
!       csf(j-i) = conjg( csf(i))
!       dsi(j-i) =        dsi(i)
!    ENDDO 
! ENDIF
! !
! END SUBROUTINE four_apply_friedel
!*******************************************************************************
! My code ( Copilot Generated )
!*******************************************************************************
SUBROUTINE four_apply_friedel
!-
!  Applies Friedels law to the Fourier calculation
!+
USE diffuse_mod
!
IMPLICIT NONE
!
INTEGER :: i, j, k !, jj
!
IF(diff_l_friedel) THEN
!  jj = diff_inc_u(1)*diff_inc_u(2)*diff_inc_u(3) + 1                             ! This entire portion looks complicated. The dimensions here are going to be tricky. Who is j ??
   DO i = 1, num(1)
      DO j = 1, num(2)
         DO k = 1, num(3)
            acsf(num(1)+1-i,num(2)+1-j,num(3)+1-k) = conjg(acsf(i,j,k))
             csf(num(1)+1-i,num(2)+1-j,num(3)+1-k) = conjg( csf(i,j,k))
!           csf(j-i,k-j,i-k) = conjg(csf(i,j,k))
!           dsi(j-i,k-j,i-k) = dsi(i,j,k)
            dsi(num(1)+1-i,num(2)+1-j,num(3)+1-k) = dsi(i,j,k    )
         ENDDO
      ENDDO
   ENDDO
ENDIF
!
END SUBROUTINE four_apply_friedel
!*******************************************************************************
SUBROUTINE four_rese_friedel
!
! Set the corners, increments back to user values
!
USE diffuse_mod
!
IMPLICIT NONE
!
IF(diff_l_friedel) THEN
   eck = diff_eck_u
   vi  = diff_vi_u
   inc = diff_inc_u
   num = diff_inc_u
ENDIF
!
END SUBROUTINE four_rese_friedel
!
!*******************************************************************************
!
!SUBROUTINE four_fft_prep
!
! Obsolete with NUFFT
!
!USE crystal_mod
!USE diffuse_mod
!!
!USE precision_mod
!!
!IMPLICIT NONE
!!
!fft_grid = 100
!IF(ilots==LOT_OFF) THEN
!   fft_dim(1) = cr_icc(1)*fft_grid
!   fft_dim(2) = cr_icc(2)*fft_grid
!ELSE
!   fft_dim(1) = ls_xyz(1)*fft_grid
!   fft_dim(2) = ls_xyz(2)*fft_grid
!ENDIF
!ALLOCATE(fft_field(1:fft_dim(1),1:fft_dim(2)))
!ALLOCATE(fft_sum  (1:fft_dim(1),1:fft_dim(2)))
!fft_field = CMPLX(0.0_PREC_DP, 0.0_PREC_DP)
!fft_sum   = CMPLX(0.0_PREC_DP, 0.0_PREC_DP)
!!
!!END SUBROUTINE four_fft_prep
!!
!!*******************************************************************************
!!
!SUBROUTINE four_strucf_fft(iscat, lform, nlot)
!!
!!
!!
!USE crystal_mod
!USE diffuse_mod
!!
!USE singleton
!!
!IMPLICIT NONE
!!
!INTEGER, INTENT(IN) :: iscat
!LOGICAL, INTENT(IN) :: lform
!!
!INTEGER :: nlot
!INTEGER :: loop
!INTEGER :: i,j
!!
!fft_field = CMPLX(0.0_PREC_DP, 0.0_PREC_DP)
!!
!DO loop=1, nxat
!   i = NINT(fft_grid*xat(loop,1))
!   j = NINT(fft_grid*xat(loop,2))
!   fft_field(i,j) = CMPLX(1.0_PREC_DP, 0.0_PREC_DP)
!ENDDO
!!
!fft_field = fft(fft_field)/REAL(fft_dim(1)*fft_dim(2),KIND=PREC_DP)
!fft_sum = fft_sum + fft_field*CONJG(fft_field)
!!
!END SUBROUTINE four_strucf_fft
!!
!!*******************************************************************************
!!
!SUBROUTINE four_fft_finalize
!!
!USE crystal_mod
!!
!IMPLICIT NONE
!!
!IF(ALLOCATED(fft_field)) DEALLOCATE(fft_field)
!IF(ALLOCATED(fft_sum  )) DEALLOCATE(fft_sum  )
!!
!END SUBROUTINE four_fft_finalize
!
!*******************************************************************************
!Pending !!!!!!! ( discuss csf_sum)
!*******************************************************************************
!
SUBROUTINE four_accumulate
!-
!  Accumulates the structure factor / intensity 
!+ 
USE diffuse_mod
!
USE errlist_mod
!
IMPLICIT NONE
!
IF(four_accum==FOUR_ACCUM_SINGLE) THEN
   IF(four_symm) THEN
      IF(ALLOCATED(csf_sum)) DEALLOCATE(csf_sum)
      IF(ALLOCATED(dsi_sum)) DEALLOCATE(dsi_sum)
      ALLOCATE(csf_sum(LBOUND(csf,1):UBOUND(csf,1), LBOUND(csf,2):UBOUND(csf,2), LBOUND(csf,3):UBOUND(csf,3)))
      csf_sum = CMPLX(0.0D0, 0.0D0)
      ALLOCATE(dsi_sum(LBOUND(dsi,1):UBOUND(dsi,1), LBOUND(dsi,2):UBOUND(dsi,2), LBOUND(dsi,3):UBOUND(dsi,3)))                                                     ! Discuss with Neder this _sum arrays.
      dsi_sum = 0.0D0
      csf_sum = csf
      dsi_sum = dsi
      call four_fill_csf                         ! Apply symm, if Set aver was used, fill Bragg
      call four_fill_dsi                         ! Apply symm, if Set aver was used, fill Bragg
      csf = csf_sum
      dsi = dsi_sum
   ELSE
      RETURN                                      ! Single diffracation pattern mode no symmetry average
   ENDIF
ELSEIF(four_accum==FOUR_ACCUM_INIT) THEN          ! Initialization, clear *_sum arays
   IF(ALLOCATED(csf_sum)) DEALLOCATE(csf_sum)
   IF(ALLOCATED(dsi_sum)) DEALLOCATE(dsi_sum)
ELSEIF(four_accum==FOUR_ACCUM_ACCUM) THEN         ! Accumulate into current sum
   IF(.NOT.ALLOCATED(csf_sum)) THEN
      ALLOCATE(csf_sum(LBOUND(csf,1):UBOUND(csf,1), LBOUND(csf,2):UBOUND(csf,2), LBOUND(csf,3):UBOUND(csf,3)))
      csf_sum = CMPLX(0.0D0, 0.0D0)
   ENDIF
   IF(.NOT.ALLOCATED(dsi_sum)) THEN
!     ALLOCATE(dsi_sum(LBOUND(dsi,1):UBOUND(dsi,1)))
      ALLOCATE(dsi_sum(LBOUND(dsi,1):UBOUND(dsi,1), LBOUND(dsi,2):UBOUND(dsi,2), LBOUND(dsi,3):UBOUND(dsi,3)))
      dsi_sum = 0.0D0
   ENDIF
   csf_sum = csf_sum + csf
   dsi_sum = dsi_sum + dsi
ELSEIF(four_accum==FOUR_ACCUM_FINISHED) THEN      ! Finished, copy everything into the original
   IF(.NOT.ALLOCATED(csf_sum) .OR. .NOT.ALLOCATED(dsi_sum)) THEN
      ier_num = -175
      ier_typ = ER_APPL
      RETURN
   ENDIF
   call four_fill_csf                         ! Apply sym, if Set aver was used, fill Bragg
   call four_fill_dsi                         ! Apply sym, if Set aver was used, fill Bragg
   csf = csf_sum
   dsi = dsi_sum
   DEALLOCATE(csf_sum)
   DEALLOCATE(dsi_sum)
ENDIF
!
END SUBROUTINE four_accumulate
!
!*******************************************************************************
! This one is already 3D
!*******************************************************************************
!
SUBROUTINE four_fill_csf
!-
!  If set aver was used with faver > 0, we need to fill the Bragg intensities
!+
USE diffuse_mod
use precision_mod
!
IMPLICIT NONE
!
COMPLEX(KIND=KIND(0.0D0)), DIMENSION(:,:,:), ALLOCATABLE :: csf_3d
!
INTEGER :: h,k,l
INTEGER :: ii
INTEGER :: np
REAL(KIND=PREC_DP    ), DIMENSION(3) :: hkl
REAL(KIND=PREC_DP    ), DIMENSION(3) :: step
COMPLEX(KIND=PREC_DP    )            :: aaa
!
IF(fave==0.0 .AND. .NOT. four_symm) RETURN          ! nothing to do
!
step(1) = MAX(vi(1,1),vi(1,2), vi(1,3))*0.005      ! Get a step size that
step(2) = MAX(vi(2,1),vi(2,2), vi(2,3))*0.005      ! is used to check if hkl is integer
step(3) = MAX(vi(3,1),vi(3,2), vi(3,3))*0.005
IF(num(1)==1) step(1) = 0.0001
IF(num(2)==1) step(2) = 0.0001
IF(num(3)==1) step(3) = 0.0001
!
ALLOCATE(csf_3d(num(1), num(2), num(3)))
!
ii = 0
      DO h = 1, num (1) 
   DO k = 1, num (2) 
DO l = 1, num (3) 
!        ii       = ii + 1 
         csf_3d(h,k,l) = csf_sum(h,k,l)
      ENDDO
   ENDDO
ENDDO
!
IF(four_symm) call four_symm_csf(num, csf_3d, eck, vi)
!
!
IF(fave < 0.0) THEN
         DO h = 1, num (1) 
      DO k = 1, num (2) 
   DO l = 1, num (3) 
            hkl(1) = eck(1,1) + (h-1)*vi(1,1) + (k-1)*vi(1,2) + (l-1)*vi(1,3)
            hkl(2) = eck(2,1) + (h-1)*vi(2,1) + (k-1)*vi(2,2) + (l-1)*vi(2,3)
            hkl(3) = eck(3,1) + (h-1)*vi(3,1) + (k-1)*vi(3,2) + (l-1)*vi(3,3)
            IF(ABS(hkl(1)-NINT(hkl(1)))<step(1) .AND.  &
               ABS(hkl(2)-NINT(hkl(2)))<step(2) .AND.  &
               ABS(hkl(3)-NINT(hkl(3)))<step(3)       ) THEN       ! Got integer hkl
               aaa = 0.0
               np  = 0
               IF(h>1) THEN                        ! Add point at h-1
                  aaa = aaa + csf_3d(h-1,k,l)
                  np  = np + 1
               ENDIF
               IF(h<num(1)) THEN                   ! Add point at h+1
                  aaa = aaa + csf_3d(h+1,k,l)
                  np  = np + 1
               ENDIF
               IF(k>1) THEN                        ! Add point at h-1
                  aaa = aaa + csf_3d(h,k-1,l)
                  np  = np + 1
               ENDIF
               IF(k<num(2)) THEN                   ! Add point at h+1
                  aaa = aaa + csf_3d(h,k+1,l)
                  np  = np + 1
               ENDIF
               IF(l>1) THEN                        ! Add point at h-1
                  aaa = aaa + csf_3d(h,k,l-1)
                  np  = np + 1
               ENDIF
               IF(l<num(3)) THEN                   ! Add point at h+1
                  aaa = aaa + csf_3d(h,k,l+1)
                  np  = np + 1
               ENDIF
               IF(np>0) THEN
                  csf_3d(h,k,l) = aaa/REAL(np)
               ENDIF
            ENDIF
         ENDDO
      ENDDO
   ENDDO
ENDIF
!
ii = 0
      DO h = 1, num (1) 
   DO k = 1, num (2) 
DO l = 1, num (3) 
!        ii       = ii + 1 
         csf_sum(h,k,l) = csf_3d(h,k,l)
      ENDDO
   ENDDO
ENDDO
!
DEALLOCATE(csf_3d)
!
END SUBROUTINE four_fill_csf
!
!*******************************************************************************
!
SUBROUTINE four_fill_dsi
!-
!  If set aver was used with faver > 0, we need to fill the Bragg intensities
!+
USE diffuse_mod
!
IMPLICIT NONE
!
REAL(KIND=KIND(0.0D0))   , DIMENSION(:,:,:), ALLOCATABLE :: dsi_3d
!
INTEGER :: h,k,l
INTEGER :: ii
INTEGER :: np
REAL(kind=PREC_DP), DIMENSION(3) :: hkl
REAL(kind=PREC_DP), DIMENSION(3) :: step
REAL(KIND=PREC_DP)               :: aaa
!
IF(fave==0.0 .AND. .NOT. four_symm) RETURN          ! nothing to do
!
step(1) = MAX(vi(1,1),vi(1,2), vi(1,3))*0.005
step(2) = MAX(vi(2,1),vi(2,2), vi(2,3))*0.005
step(3) = MAX(vi(3,1),vi(3,2), vi(3,3))*0.005
IF(num(1)==1) step(1) = 0.0001
IF(num(2)==1) step(2) = 0.0001
IF(num(3)==1) step(3) = 0.0001
!
ALLOCATE(dsi_3d(num(1), num(2), num(3)))
!
ii = 0
      DO h = 1, num (1) 
   DO k = 1, num (2) 
DO l = 1, num (3) 
         ii       = ii + 1 
         dsi_3d(h,k,l) = dsi_sum(h,k,l)
      ENDDO
   ENDDO
ENDDO
!
IF(four_symm) call four_symm_dsi(num, dsi_3d, eck, vi)
!
IF(fave < 0.0) THEN
         DO h = 1, num (1) 
      DO k = 1, num (2) 
   DO l = 1, num (3) 
            hkl(1) = eck(1,1) + (h-1)*vi(1,1) + (k-1)*vi(1,2) + (l-1)*vi(1,3)
            hkl(2) = eck(2,1) + (h-1)*vi(2,1) + (k-1)*vi(2,2) + (l-1)*vi(2,3)
            hkl(3) = eck(3,1) + (h-1)*vi(3,1) + (k-1)*vi(3,2) + (l-1)*vi(3,3)
            IF(ABS(hkl(1)-NINT(hkl(1)))<step(1) .AND.  &
               ABS(hkl(2)-NINT(hkl(2)))<step(2) .AND.  &
               ABS(hkl(3)-NINT(hkl(3)))<step(3)       ) THEN       ! Got integer hkl
               aaa = 0.0
               np  = 0
               IF(h>1) THEN                        ! Add point at h-1
                  aaa = aaa + dsi_3d(h-1,k,l)
                  np  = np + 1
               ENDIF
               IF(h<num(1)) THEN                   ! Add point at h+1
                  aaa = aaa + dsi_3d(h+1,k,l)
                  np  = np + 1
               ENDIF
               IF(k>1) THEN                        ! Add point at h-1
                  aaa = aaa + dsi_3d(h,k-1,l)
                     np  = np + 1
               ENDIF
               IF(k<num(2)) THEN                   ! Add point at h+1
                  aaa = aaa + dsi_3d(h,k+1,l)
                  np  = np + 1
               ENDIF
               IF(l>1) THEN                        ! Add point at h-1
                  aaa = aaa + dsi_3d(h,k,l-1)
                  np  = np + 1
               ENDIF
               IF(l<num(3)) THEN                   ! Add point at h+1
                  aaa = aaa + dsi_3d(h,k,l+1)
                  np  = np + 1
               ENDIF
               IF(np>0) THEN
                  dsi_3d(h,k,l) = aaa/REAL(np)
               ENDIF
            ENDIF
         ENDDO
      ENDDO
   ENDDO
ENDIF
!
ii = 0
      DO h = 1, num (1) 
   DO k = 1, num (2) 
DO l = 1, num (3) 
         ii       = ii + 1 
         dsi_sum(h,k,l) = dsi_3d(h,k,l)
      ENDDO
   ENDDO
ENDDO
!
DEALLOCATE(dsi_3d)
!
END SUBROUTINE four_fill_dsi
!
!*******************************************************************************
! Also already 3D
!*******************************************************************************
!
SUBROUTINE four_symm_csf(num, csf_3d, eck, vi)
!-
!   Apply reciprocal space symmetry of point group
!+
USE crystal_mod
use wyckoff_mod
!
USE matrix_mod
use precision_mod
!
IMPLICIT NONE
!
INTEGER                  , DIMENSION(3)                     , INTENT(IN)    :: num
COMPLEX(KIND=KIND(0.0D0)), DIMENSION(num(1), num(2), num(3)), INTENT(INOUT) :: csf_3d
REAL(kind=PREC_DP)       , DIMENSION(1:3, 1:4)           ::  eck
REAL(kind=PREC_DP)       , DIMENSION(1:3, 1:3)           ::  vi
!
INTEGER :: n_center     ! Only use this fraction of symmetry operations
INTEGER :: i
INTEGER :: h,k,l
INTEGER :: h1,k1,l1
INTEGER, DIMENSION(3) :: izero                                    ! Indices of point 0,0,0
REAL(KIND=PREC_DP), DIMENSION(3) :: hkl
REAL(KIND=PREC_DP), DIMENSION(3) :: uvw
COMPLEX(KIND=KIND(0.0D0)), DIMENSION(:,:,:), ALLOCATABLE :: csf_sym
INTEGER                  , DIMENSION(:,:,:), ALLOCATABLE :: weight
REAL(KIND=PREC_DP)       , DIMENSION(3, 3)               :: sym_mat
REAL(KIND=PREC_DP)       , DIMENSION(3, 3)               :: tmp_mat
!
!sym_mat(:,1) = vi(:,1)              ! Find pixels at which hkl = (0,0,0)
!sym_mat(:,2) = vi(:,2)
!sym_mat(:,3) = vi(:,3)
!sym_mat(:,3) = (/ 0.0, 0.0, 1.0/)
!
if(num(1)>1) then
   sym_mat(:,1) = vi(:,1)              ! Find pixels at which hkl = (0,0,0)
else
   sym_mat(:,1) = (/ 1.0, 0.0, 0.0/)
endif
!
if(num(2)>1) then
   sym_mat(:,2) = vi(:,2)              ! Find pixels at which hkl = (0,0,0)
else
   sym_mat(:,2) = (/ 0.0, 1.0, 0.0/)
endif
!
if(num(3)>1) then
   sym_mat(:,3) = vi(:,3)              ! Find pixels at which hkl = (0,0,0)
else
   sym_mat(:,3) = (/ 0.0, 0.0, 1.0/)
endif
!
call matinv3(sym_mat, tmp_mat)
!
hkl(:) = -eck(:,1)
uvw = MATMUL(tmp_mat, hkl)
do i=1,3
   if(num(i)>1) then
      izero(i) = NINT(1+uvw(i))            ! Pixel at hkl = (0,0,0)
   else
      izero(i) = 1
   endif
enddo
!
!                                     ! Determine number of primitive operations
!
n_center = 1
IF (cr_spcgr (1:1) .eq.'P') THEN
   n_center = 1
ELSEIF (cr_spcgr (1:1) .eq.'A') THEN  ! Normal space group can be used
   n_center = 2
ELSEIF (cr_spcgr (1:1) .eq.'B') THEN  ! as n_center is identical for
   n_center = 2
ELSEIF (cr_spcgr (1:1) .eq.'C') THEN  ! A B and C, orthorhombic alternative setting
   n_center = 2
ELSEIF (cr_spcgr (1:1) .eq.'I') THEN
   n_center = 2
ELSEIF (cr_spcgr (1:1) .eq.'F') THEN
   n_center = 4
ELSEIF (cr_spcgr (1:1) .eq.'R'.and.cr_syst.eq.6) THEN
   n_center = 3
ENDIF
!
ALLOCATE(csf_sym(num(1), num(2), num(3)))
ALLOCATE(weight (num(1), num(2), num(3)))      ! Weight counts the number that symm adds intensity to a pixel
weight = 1
!
csf_sym = csf_3d                     ! Effectively operation 1
!
!
DO i = 2, spc_n/n_center             ! only apply numbers 2 for P centered part
   sym_mat=spc_mat(1:3,1:3,i)        ! Copy symmetry matrix to local copy
   DO l = 1, num (3) 
      hkl(3) = l - izero(3)
      DO k = 1, num (2) 
         hkl(2) = k - izero(2)
         DO h = 1, num (1) 
            hkl(1) = h - izero(1)
            uvw    = MATMUL(hkl, sym_mat)
            h1 = izero(1) + NINT(uvw(1))
            k1 = izero(2) + NINT(uvw(2))
            l1 = izero(3) + NINT(uvw(3))
            IF(h1>=1 .AND. h1<=num(1) .AND.          &
               k1>=1 .AND. k1<=num(2) .AND.          &
               l1>=1 .AND. l1<=num(3)      ) THEN
               csf_sym(h1,k1,l1) = csf_sym(h1, k1, l1) + csf_3d(h,k,l)
               weight (h1,k1,l1) = weight (h1, k1, l1) + 1
            ENDIF
         ENDDO
      ENDDO
   ENDDO
ENDDO
!
! Apply symmetry weight and divide by number of Symmetry operations
!
csf_3d = csf_sym/weight/REAL(spc_n/n_center)! Copy back into original
!
DEALLOCATE(csf_sym)
DEALLOCATE(weight )
!
END SUBROUTINE four_symm_csf
!
!*******************************************************************************
!
SUBROUTINE four_symm_dsi(num, dsi_3d, eck, vi)
!-
!   Apply reciprocal space symmetry of point group
!+
USE crystal_mod
use wyckoff_mod
!
use precision_mod
USE matrix_mod
!
IMPLICIT NONE
!
INTEGER           , DIMENSION(3)                     , INTENT(IN)    :: num
REAL(KIND=PREC_DP), DIMENSION(num(1), num(2), num(3)), INTENT(INOUT) :: dsi_3d
REAL(kind=PREC_DP), DIMENSION(1:3, 1:4)           ::  eck
REAL(kind=PREC_DP), DIMENSION(1:3, 1:3)           ::  vi
!
INTEGER :: n_center     ! Only use this fraction of symmetry operations
INTEGER :: i
INTEGER :: h,k,l
INTEGER :: h1,k1,l1
INTEGER, DIMENSION(3) :: izero                                    ! Indices of point 0,0,0
REAL(KIND=PREC_DP), DIMENSION(3) :: hkl
REAL(KIND=PREC_DP), DIMENSION(3) :: uvw
REAL(KIND=KIND(0.0D0)), DIMENSION(:,:,:), ALLOCATABLE :: dsi_sym
INTEGER               , DIMENSION(:,:,:), ALLOCATABLE :: weight
REAL(KIND=PREC_DP), DIMENSION(3, 3)                   :: sym_mat
REAL(KIND=PREC_DP), DIMENSION(3, 3)                   :: tmp_mat
!
if(num(1)>1) then
   sym_mat(:,1) = vi(:,1)              ! Find pixels at which hkl = (0,0,0)
else
   sym_mat(:,1) = (/ 1000.0, 0.0, 0.0/)
endif
!
if(num(2)>1) then
   sym_mat(:,2) = vi(:,2)              ! Find pixels at which hkl = (0,0,0)
else
   sym_mat(:,2) = (/ 0.0, 1000.0, 0.0/)
endif
!
if(num(3)>1) then
   sym_mat(:,3) = vi(:,3)              ! Find pixels at which hkl = (0,0,0)
else
   sym_mat(:,3) = (/ 0.0, 0.0, 1000.0/)
endif
!
call matinv3(sym_mat, tmp_mat)
!
hkl(:) = -eck(:,1)
uvw = MATMUL(tmp_mat, hkl)
do i=1,3
   if(num(i)>1) then
      izero(i) = NINT(1+uvw(i))            ! Pixel at hkl = (0,0,0)
   else
      izero(i) = 1
   endif
enddo
!
!                                     ! Determine number of primitive operations
!write(*,*) ' IZERO ', izero
!
n_center = 1
IF (cr_spcgr (1:1) .eq.'P') THEN
   n_center = 1
ELSEIF (cr_spcgr (1:1) .eq.'A') THEN  ! Normal space group can be used
   n_center = 2
ELSEIF (cr_spcgr (1:1) .eq.'B') THEN  ! as n_center is identical for
   n_center = 2
ELSEIF (cr_spcgr (1:1) .eq.'C') THEN  ! A B and C, orthorhombic alternative setting
   n_center = 2
ELSEIF (cr_spcgr (1:1) .eq.'I') THEN
   n_center = 2
ELSEIF (cr_spcgr (1:1) .eq.'F') THEN
   n_center = 4
ELSEIF (cr_spcgr (1:1) .eq.'R'.and.cr_syst.eq.6) THEN
   n_center = 3
ENDIF
!
ALLOCATE(dsi_sym(num(1), num(2), num(3)))
ALLOCATE(weight (num(1), num(2), num(3)))      ! Weight counts the number that symm adds intensity to a pixel
weight = 1
!
dsi_sym = dsi_3d                     ! Effectively operation 1
!!!
DO i = 2, spc_n/n_center             ! only apply numbers 2 for P centered part
   sym_mat=spc_mat(1:3,1:3,i)        ! Copy symmetry matrix to local copy
!write(*,'(a, 3F6.1)') ' SYMM ', sym_mat(1,:)
!write(*,'(a, 3F6.1)') ' SYMM ', sym_mat(2,:)
!write(*,'(a, 3F6.1)') ' SYMM ', sym_mat(3,:)
   DO l = 1, num (3) 
      hkl(3) = l - izero(3)
      DO k = 1, num (2) 
         hkl(2) = k - izero(2)
         DO h = 1, num (1) 
            hkl(1) = h - izero(1)
            uvw    = MATMUL(hkl, sym_mat)
            h1 = izero(1) + NINT(uvw(1))
            k1 = izero(2) + NINT(uvw(2))
            l1 = izero(3) + NINT(uvw(3))
            IF(h1>=1 .AND. h1<=num(1) .AND.          &
               k1>=1 .AND. k1<=num(2) .AND.          &
               l1>=1 .AND. l1<=num(3)      ) THEN
               dsi_sym(h1,k1,l1) = dsi_sym(h1, k1, l1) + dsi_3d(h,k,l)
               weight (h1,k1,l1) = weight (h1, k1, l1) + 1
            ENDIF
         ENDDO
      ENDDO
   ENDDO
ENDDO
!
! Apply symmetry weight and divide by number of Symmetry operations
!
dsi_3d = dsi_sym/weight/REAL(spc_n/n_center,kind=PREC_DP)! Copy back into original
!
DEALLOCATE(dsi_sym)
DEALLOCATE(weight )
!
END SUBROUTINE four_symm_dsi
!
!*******************************************************************************
!
END MODULE fourier_sup
