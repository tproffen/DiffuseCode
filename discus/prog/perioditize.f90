module perioditize_mod
!-
!  Attempt to rearrange the structure as periodic structure in the standard
!  DISCUS sequence
!
use precision_mod
!
private
public perioditize_menu
!
integer :: pdt_ncells     ! Number of cells in periodic crystal volume
integer :: pdt_nsite = 0  ! Number of sites in an averaged unit cell
integer :: pdt_usite = 0  ! Number of sites in an averaged unit cell, defined by user with 'set site'
integer           , dimension(:,:,:), allocatable :: pdt_cell     ! Atoms per cell
integer           , dimension(3)                  :: pdt_ilow     ! Unit cell dimensions in periodic
integer           , dimension(3)                  :: pdt_ihig     ! low and high inidces
logical                                           :: pdt_usr_nsite  ! True if user fixed number of sites
logical                                           :: pdt_usr_atom   ! True if user fixed atoms at sites
real(kind=PREC_SP), dimension(3, 2)               :: pdt_dims ! Crystal dimensions
real(kind=PREC_SP), dimension(:,:)  , allocatable :: pdt_pos  ! Atom positions in average cell
!
contains
!
!*******************************************************************************
!
subroutine perioditize_menu
!-
!  Main menu routine for periodtizing
!+
!
use discus_kdo_common_mod
!
use calc_expr_mod
use class_macro_internal
use doact_mod
use errlist_mod
use lib_length
use lib_errlist_func
use lib_macro_func
use macro_mod
!use macro_internal_mod
use precision_mod
use prompt_mod
use str_comp_mod
use sup_mod
!
implicit none
!
character(len=PREC_STRING) :: zeile        ! Input line
character(len=PREC_STRING) :: line         ! Input line
character(len=6          ) :: befehl       ! The actual command
character(len=len(prompt)) :: orig_prompt
!
integer :: indxg
integer :: length            ! Input line length
integer :: lbef              ! Command length
integer :: lp                ! Line length
logical :: lend              ! if true terminate the menu
logical :: success
!
call no_error
!
orig_prompt = prompt
prompt = prompt (1:len_str (prompt) ) //'/period'
lend   = .FALSE.
!
loop_menu: do
   if(ier_num == 0) then
      CALL get_cmd (line, length, befehl, lbef, zeile, lp, prompt)  ! get new command
   endif
   if(lend) exit loop_menu
!
   IF(ier_num /= 0) THEN   ! Error in previous command or in get_cmd
      CALL errlist
      IF (ier_sta.ne.ER_S_LIVE) THEN
         IF (lmakro .OR. lmakro_error) THEN  ! Error within macro or termination errror
            IF(sprompt /= prompt ) THEN
               ier_num = -10
               ier_typ = ER_COMM
               ier_msg(1) = ' Error occured in symmetry menu'
               prompt_status = PROMPT_ON 
               prompt = orig_prompt
               RETURN
            ELSE
               IF(lmacro_close) THEN
                  CALL macro_close 
                  prompt_status = PROMPT_ON 
               ENDIF 
            ENDIF 
         ENDIF 
         IF (lblock) THEN 
            ier_num = - 11 
            ier_typ = ER_COMM 
            prompt_status = PROMPT_ON 
            prompt = orig_prompt
            RETURN 
         ENDIF 
         CALL no_error 
         lmakro_error = .FALSE.
         sprompt = ' '
      ENDIF 
      cycle loop_menu
   ENDIF 
!
!  Normal menu commands
!
   IF (line == ' '      .or. line(1:1) == '#' .or. &
       line == char(13) .or. line(1:1) == '!'        ) cycle loop_menu
!                                                                       
!     ----search for "="                                                
!                                                                       
   indxg = index (line, '=') 
   IF (indxg.ne.0.AND..NOT. (str_comp (befehl, 'echo', 2, lbef, 4) ) &
                 .AND..NOT. (str_comp (befehl, 'syst', 2, lbef, 4) )    &
                 .AND..NOT. (str_comp (befehl, 'help', 2, lbef, 4) .OR. &
                             str_comp (befehl, '?   ', 2, lbef, 4) )    &
                 .AND. INDEX(line,'==') == 0                            ) THEN
!                                                                       
!     ---evaluate an expression and assign the value to a variabble   
!                                                                       
      CALL do_math (line, indxg, length) 
      cycle loop_menu
   endif
!
! - Test common menu entries
!
   CALL discus_kdo_common(befehl, lbef, line, length, zeile, lp, 'symm' , &
                          lend, success)
   if(lend) exit loop_menu
   if(success) cycle loop_menu
!                                                                       
!     ----run symmetry 'rese'                                            
!                                                                       
   if_rese: IF(str_comp (befehl, 'rese', 2, lbef, 4) ) THEN 
      CALL perioditize_reset
      cycle loop_menu
   endif if_rese
!-------------------------------------------------------------------------------
!-------- Perioditize specific commands if no common command was found
!-------------------------------------------------------------------------------
!
!  IF( cr_nscat > SYM_MAXSCAT .or. mole_num_type > SYM_MAXSCAT) THEN
!     nscat = max ( cr_nscat, mole_num_type)
!     nsite = max ( nsite, cr_ncatoms, SYM_MAXSITE)
!     CALL alloc_symmetry ( nscat, nsite )
!     IF ( ier_num < 0 ) THEN
!        RETURN
!     ENDIF
!  ENDIF
!                                                                       
!     ----exit menu        'exit'                                     
!                                                                       
   if(str_comp (befehl, 'run', 2, lbef, 3) ) then
      call perioditize_run
   elseif(str_comp (befehl, 'set', 2, lbef, 3) ) then
      call perioditize_set(zeile, lp)
   else
      ier_num = -8
      ier_typ = ER_COMM
   endif
   if(lend) exit loop_menu
enddo loop_menu
!
prompt = orig_prompt
!
end subroutine perioditize_menu
!
!*******************************************************************************
!
subroutine perioditize_set(zeile, lp)
!-
!  Set user defined parameters
!+
use crystal_mod
use errlist_mod
use get_params_mod
use precision_mod
use str_comp_mod
!
implicit none
!
character(len=*), intent(inout) :: zeile      ! input line
integer         , intent(inout) :: lp         ! input length
!
integer            :: maxw                        ! Actual parameter number
integer, parameter :: MIN_PARA =  2               ! Minimum parameter numbers
character(len=PREC_STRING), dimension(MAX(MIN_PARA,MAXSCAT+1)) :: cpara
integer                   , dimension(MAX(MIN_PARA,MAXSCAT+1)) :: lpara
real(kind=PREC_DP)        , dimension(MAX(MIN_PARA,MAXSCAT+1)) :: werte
integer :: ianz          ! Number of parameters
!
!
!
MAXW = MAX(MIN_PARA,MAXSCAT+1)
!
call get_params (zeile, ianz, cpara, lpara, maxw, lp)
if(ier_num/=0) return
!
!call get_optional(ianz, MAXW, cpara, lpara, NOPTIONAL,  &
!                  ncalc, oname, loname, opara, lopara, lpresent, owerte)
!
if(str_comp (cpara(1), 'natoms', 2, lpara(1), 6) ) then
   call perioditize_set_natoms(ianz, cpara, lpara, werte, maxw)
elseif(str_comp (cpara(1), 'site', 2, lpara(2), 4) ) then
   call perioditize_set_site  (ianz, cpara, lpara, werte, maxw)
endif
!
end subroutine perioditize_set
!
!*******************************************************************************
!
subroutine perioditize_set_natoms(ianz, cpara, lpara, werte, maxw)
!-
!  Set the number of atoms per unit cell
!+
!
use ber_params_mod
use errlist_mod
use precision_mod
!
implicit none
!
integer :: MAXW          ! max number of parameters
integer :: ianz          ! Number of parameters
character(len=PREC_STRING), dimension(MAXW), intent(inout)     :: cpara
integer                   , dimension(MAXW), intent(inout)     :: lpara
real(kind=PREC_DP)        , dimension(MAXW), intent(inout)     :: werte
!
cpara(1) = '0'
lpara(1) = 1
!
if(cpara(2) == 'auto') then
   pdt_nsite = 1
   pdt_usr_nsite = .FALSE.                 ! Determine sites automatically
else
   cpara(1) = '0'
   lpara(1) = 1
   call ber_params(ianz, cpara, lpara, werte, maxw)
   if(ier_num/=0) return
   pdt_nsite = nint(werte(2))
   pdt_usr_nsite = .TRUE.                  ! User fixed number of sites
endif
end subroutine perioditize_set_natoms
!
!*******************************************************************************
!
subroutine perioditize_set_site(ianz, cpara, lpara, werte, maxw)
!-
!  Set explicit positions expected in the average unit cell
!+
!
use crystal_mod
!
use ber_params_mod
use errlist_mod
use get_iscat_mod
use get_params_mod
use precision_mod
use take_param_mod
!
implicit none
!
integer :: MAXW          ! max number of parameters
integer :: ianz          ! Number of parameters
character(len=PREC_STRING), dimension(MAXW), intent(inout)     :: cpara
integer                   , dimension(MAXW), intent(inout)     :: lpara
real(kind=PREC_DP)        , dimension(MAXW), intent(inout)     :: werte
!
!
integer, parameter :: MAXU= 3                   ! Max params for pos:[}
character(len=PREC_STRING) :: line
character(len=PREC_STRING), dimension(MAXSCAT) :: ccpara
integer                   , dimension(MAXSCAT) :: llpara
real(kind=PREC_DP)        , dimension(MAXSCAT) :: wwerte
real(kind=PREC_DP)        , dimension(MAXU)    :: uwerte
integer :: iianz          ! Number of parameters
integer :: jjanz          ! Number of parameters
integer :: maxww          ! Number of parameters
integer :: length
integer :: osite
integer :: nsite
!
!
integer, parameter :: NOPTIONAL = 2
integer, parameter :: O_ATOM    = 1
integer, parameter :: O_POS     = 2
character(LEN=   6), dimension(NOPTIONAL) :: oname   !Optional parameter names
character(LEN=PREC_STRING), dimension(NOPTIONAL) :: opara   !Optional parameter strings returned
integer            , dimension(NOPTIONAL) :: loname  !Lenght opt. para name
integer            , dimension(NOPTIONAL) :: lopara  !Lenght opt. para name returned
logical            , dimension(NOPTIONAL) :: lpresent!opt. para is present
real(kind=PREC_DP) , dimension(NOPTIONAL) :: owerte   ! Calculated values
integer, parameter                        :: ncalc = 1 ! Number of values to calculate 
!
data oname  / 'atoms ', 'pos'   /
data loname /  5,    3          /
opara  =  (/ 'all    ', '[0,0,0]' /)   ! Always provide fresh default values
lopara =  (/  3,        76        /)
owerte =  (/  0.0,      0.00      /)
!
MAXWW = MAXSCAT
!
call get_optional(ianz, MAXW, cpara, lpara, NOPTIONAL,  ncalc, &
                  oname, loname, opara, lopara, lpresent, owerte)
cpara(1) = '0'
lpara(1) = 1
!
if(cpara(2) == 'auto') then
   pdt_nsite = 1
   pdt_usr_atom  = .FALSE.                 ! Determine sites automatically
   pdt_usite     = 0                       ! No sites fixed by user
   return
else
   if(lpresent(O_ATOM)) then               ! atom:[] is present
      if(opara(O_ATOM)=='all' ) then       ! Allow all atom ptypes at this site
         pdt_usr_atom = .TRUE.             ! User did flag atom positions
      else
         if(opara(O_ATOM)(1:1)=='[' .and. opara(O_ATOM)(lopara(O_ATOM):lopara(O_ATOM))==']') then
            line = opara(O_ATOM)(2:lopara(O_ATOM)-1)
            length = lopara(O_ATOM)-2
         else
            line = opara(O_ATOM)
            length = lopara(O_ATOM)
         endif
         call get_params(line, iianz, ccpara, llpara, maxwW, length)
         call get_iscat (iianz, ccpara, llpara, wwerte, MAXWW, .FALSE. )   ! No new atom types
      endif
   else
      ier_num = -6
      ier_typ = ER_COMM
      ier_msg(1) = '''atoms:'' parameter is missing'
      return
   endif
   if(lpresent(O_POS )) then               ! atom:[] is present
      line = opara(O_POS)
      length = lopara(O_POS)
      jjanz = 3
      call get_optional_multi(MAXU, line, length, uwerte, jjanz)
   else
      ier_num = -6
      ier_typ = ER_COMM
      ier_msg(1) = '''pos:'' parameter is missing'
      return
   endif
endif
!
if(pdt_usite==ubound(pdt_pos,2)) then 
   osite = ubound(pdt_pos, 2)
   nsite = max(pdt_nsite, pdt_usite) + 5
   call perioditize_realloc(pdt_pos, osite, nsite)
endif
   
!
write(*,*) 'SITE: ', wwerte(1:iianz), uwerte(1:jjanz)
!
end subroutine perioditize_set_site
!
!*******************************************************************************
!
subroutine perioditize_reset
!
!
implicit none
!
pdt_ncells = 0
pdt_nsite  = 0
pdt_usite  = 0
pdt_ilow   = 0
pdt_ihig   = 0
pdt_dims   = 0.0D0
pdt_usr_nsite = .FALSE.
pdt_usr_atom  = .FALSE.
!
end subroutine perioditize_reset
!
!*******************************************************************************
!
subroutine perioditize_run
!-
! Perform the actual perioditization
!+
use precision_mod
!
implicit none
!
real(kind=PREC_SP) :: aver
real(kind=PREC_SP) :: sigma
!
if(.not. (pdt_usr_nsite) ) then            ! User did not define number of atoms
   call estimate_ncatom(aver, sigma)
endif
if(.not. (pdt_usr_atom) ) then             ! User did not specify unit cell content
   call find_average(aver)
endif
call map_to_aver
if(allocated(pdt_pos)) deallocate(pdt_pos)
!
end subroutine perioditize_run
!
!*******************************************************************************
!
subroutine estimate_ncatom(aver, sigma)
!-
!  Try to estimate the number of atoms per unit cell
!+
use crystal_mod
!
use precision_mod
!
implicit none
!
real(kind=PREC_SP), intent(out) :: aver
real(kind=PREC_SP), intent(out) :: sigma
!
integer :: i,j,k, l   ! Dummy counter
integer :: natoms     ! Total atoms in search window
integer :: maxn       ! Max number per unit cell
integer           , dimension(3)                  :: ixyz         ! Atom is in this unit cell
integer           , dimension(0:100)              :: counter
!
pdt_dims(:,1) = cr_pos(:,1)           ! Initialize the min/ max dimensions
pdt_dims(:,2) = cr_pos(:,1)
!
do i=1, cr_natoms
   pdt_dims(1,1) = min(pdt_dims(1,1), cr_pos(1,i))
   pdt_dims(1,2) = max(pdt_dims(1,2), cr_pos(1,i))
   pdt_dims(2,1) = min(pdt_dims(1,1), cr_pos(2,i))
   pdt_dims(2,2) = max(pdt_dims(2,2), cr_pos(2,i))
   pdt_dims(3,1) = min(pdt_dims(3,1), cr_pos(3,i))
   pdt_dims(3,2) = max(pdt_dims(3,2), cr_pos(3,i))
enddo
!
write(*,'(a, f9.5, 2x, f9.5)') ' dims x ', pdt_dims(1,:)
write(*,'(a, f9.5, 2x, f9.5)') ' dims y ', pdt_dims(2,:)
write(*,'(a, f9.5, 2x, f9.5)') ' dims z ', pdt_dims(3,:)
!
pdt_ilow = 0
pdt_ihig = 0
do i=1, 3
   if(nint(pdt_dims(i,2)-pdt_dims(i,1))>2) then
   pdt_ilow(i) = nint(pdt_dims(i,1))! + 1
   pdt_ihig(i) = nint(pdt_dims(i,2))! - 1
   else
      pdt_ilow(i) = nint(pdt_dims(i,1))
      pdt_ihig(i) = nint(pdt_dims(i,2))
   endif
write(*,'( i3,  i3, i3 )' ) i, pdt_ilow(i), pdt_ihig(i)
enddo
allocate(pdt_cell(pdt_ilow(1):pdt_ihig(1), pdt_ilow(2):pdt_ihig(2),pdt_ilow(3):pdt_ihig(3)))
pdt_cell = 0
do i=1, cr_natoms
   ixyz(1) = int(cr_pos(1,i) + 0.5 - pdt_ilow(1)) + pdt_ilow(1)             ! add 0.5 to avoid loosing atoms at the low edge
   ixyz(2) = int(cr_pos(2,i) + 0.5 - pdt_ilow(2)) + pdt_ilow(2)
   ixyz(3) = int(cr_pos(3,i) + 0.5 - pdt_ilow(3)) + pdt_ilow(3)
   if(ixyz(1)>=pdt_ilow(1) .and. ixyz(1)<=pdt_ihig(1)   .and.          &
      ixyz(2)>=pdt_ilow(2) .and. ixyz(2)<=pdt_ihig(2)   .and.          &
      ixyz(3)>=pdt_ilow(3) .and. ixyz(3)<=pdt_ihig(3)         ) then
      pdt_cell(ixyz(1), ixyz(2), ixyz(3)) = pdt_cell(ixyz(1), ixyz(2), ixyz(3)) + 1
!write(*,'(i4, 3(2x,f8.3), 3i4)') i, cr_pos(:,i), ixyz
   endif
enddo
maxn = maxval(pdt_cell)
natoms = sum(pdt_cell)
pdt_ncells =  (pdt_ihig(1)-pdt_ilow(1)+1)*(pdt_ihig(2)-pdt_ilow(2)+1)*(pdt_ihig(3)-pdt_ilow(3)+1)
write(*,*) ' Nunit ', minval(pdt_cell), maxval(pdt_cell), sum(pdt_cell),  &
  (pdt_ihig(1)-pdt_ilow(1)+1)*(pdt_ihig(2)-pdt_ilow(2)+1)*(pdt_ihig(3)-pdt_ilow(3)+1)
counter = 0
do k=pdt_ilow(3), pdt_ihig(3)
   do j=pdt_ilow(2), pdt_ihig(2)
      do i=pdt_ilow(1), pdt_ihig(1)
         l = pdt_cell(i,j,k)
         counter(l) = counter(l) + 1
      enddo
   enddo
enddo
aver  = 0.0
sigma = 0.0
write(*,*) counter(0:maxn)
do i=0,maxn
  aver = aver + real(i*counter(i))
enddo
write(*,*) ' Aver, natoms ', aver, natoms, pdt_ncells
aver = aver/real(pdt_ncells)
do i=0,maxn
   sigma = sigma + counter(i)*(aver-i)**2
enddo
sigma = sqrt(sigma)/real(pdt_ncells*(pdt_ncells-1))
write(*,*) ' Av sig ', aver, sigma, nint(aver)*pdt_ncells, cr_natoms


!
deallocate(pdt_cell)
!
end subroutine estimate_ncatom
!
!*******************************************************************************
!
subroutine find_average(aver)
!-
!  Find the clusters that form the sites in the average unit cell
!+
!
use crystal_mod
use metric_mod
!
implicit none
!
real(kind=PREC_SP), intent(in) :: aver          ! Average number of atoms per unit cell
!
logical, parameter :: lspace =.TRUE.
integer :: i, j, ic                    ! Dummy counters
real(kind=PREC_SP) :: eps              ! Distance to say we are in the same cluster
real(kind=PREC_SP), dimension(3,2*nint(aver)) :: pdt_temp      ! Positions of the average cell
real(kind=PREC_SP), dimension(3)              :: uvw           ! fractional coordinates of atom
real(kind=PREC_SP), dimension(3)              :: u             ! fractional coordinates of atom
real(kind=PREC_SP) :: dist
real(kind=PREC_SP) :: dmin
integer           , dimension(  2*nint(aver)) :: pdt_temp_n    ! Number of atoms per site
!
eps = 0.9                              ! Start with 0.5 Angstroem
pdt_nsite = 1                              ! No sites yet
pdt_temp(:, 1) =    (cr_pos(:,1) + 0.0 + 2*nint(abs(pdt_dims(:,1)))) - &
                nint(cr_pos(:,1) + 0.0 + 2*nint(abs(pdt_dims(:,1))))            ! Arbitrarily put atom one onto site 1
pdt_temp_n(1)  = 1                                                              ! One atom at this site
write(*,*) 'SITE 1 ', pdt_temp(:,1)
write(*,*) 'Atom 1 ', cr_pos(:,1)
write(*,*) 'Atom 1 ', cr_pos(:,1) + 2*nint(abs(pdt_dims(:,1)))
write(*,*) 'CELL 1 ', nint(cr_pos(:,1) + 2*nint(abs(pdt_dims(:,1))))
!
do i=2, cr_natoms
!  do j=1,3
   uvw(:) =  (cr_pos(:,i) + 0.0 + 2*nint(abs(pdt_dims(:,1)))) - &
         nint(cr_pos(:,i) + 0.0 + 2*nint(abs(pdt_dims(:,1))))
!  enddo
   dmin = 1e12
   ic   = 0
   do j=1, pdt_nsite
      u = pdt_temp(:,j)
      dist = do_blen(lspace, u, uvw)
      if(dist < dmin) then                 ! Found new shortest distance
         dmin = dist
         ic   = j                          ! Archive site
      endif
   enddo
      if(dmin < eps) then                  ! Atom belongs to cluster
         pdt_temp(:, ic) = (pdt_temp(:,ic)*pdt_temp_n(ic) + uvw) / (pdt_temp_n(ic) + 1)
         pdt_temp_n( ic) = pdt_temp_n( ic) + 1
      else                                 ! To far , new cluster
         pdt_nsite = pdt_nsite + 1
         pdt_temp(:, pdt_nsite) = uvw
         pdt_temp_n( pdt_nsite) = 1
      endif
enddo
!
allocate(pdt_pos(3,pdt_nsite))                 ! Atom positions in the average unit cell
do j=1, pdt_nsite
  pdt_pos(:,j) = pdt_temp(:,j)
  write(*,*) ' SITE ', j, pdt_pos(:,j), pdt_temp_n(j)
enddo
write(*,*) ' NSITE ', pdt_nsite*pdt_ncells, pdt_nsite*pdt_ncells==cr_natoms

end subroutine find_average
!
!*******************************************************************************
!
subroutine map_to_aver
!-
!  Map the current structure to an averaged and periodic crystal
!+
!
use crystal_mod
use chem_mod
use metric_mod
use molecule_mod
!
implicit none
!
logical, parameter :: lspace =.TRUE.
real(kind=PREC_SP), parameter :: EPS = 2.5
integer :: i,j,k,l,m,n           ! dummy counters
integer :: nprior                ! Number of atoms in crystal
integer :: natoms                ! Number of atoms in new crystal
real(kind=PREC_SP) :: dmin
real(kind=PREC_SP) :: dist
real(kind=PREC_SP), dimension(3) :: u             ! fractional coordinates of atom
real(kind=PREC_SP), dimension(3) :: v             ! fractional coordinates of atom
!
integer, dimension(  :), allocatable ::  tmp_iscat  ! (  1:NMAX)  !Atom type 0 to cr_nscat
integer, dimension(  :), allocatable ::  tmp_prop   ! (  1:NMAX)  !Property flag
integer, dimension(  :), allocatable ::  tmp_mole   ! (  1:NMAX)  !Atom is in this molecule
integer, dimension(:,:), allocatable ::  tmp_surf   ! (  1:NMAX)  !Atom is on this surface 
real   , dimension(:,:), allocatable ::  tmp_magn   ! (  1:NMAX)  !Magnetic moment 
real   , dimension(:,:), allocatable ::  tmp_pos    ! (3,1:NMAX)  !Atom coordinates
!
nprior = cr_natoms
natoms = pdt_nsite*pdt_ncells
!
allocate(tmp_iscat(     1:natoms))
allocate(tmp_prop (     1:natoms))
allocate(tmp_mole (     1:natoms))
allocate(tmp_surf (0:3, 1:natoms))
allocate(tmp_magn (0:3, 1:natoms))
allocate(tmp_pos  (3  , 1:natoms))
!
tmp_iscat  = -1                        ! Initialize with an error flag
!
!  Build the periodic crystal
!
m = 0
do k=pdt_ilow(3), pdt_ihig(3)
   do j=pdt_ilow(2), pdt_ihig(2)
      do i=pdt_ilow(1), pdt_ihig(1)
         do l=1, pdt_nsite
            m = m + 1
            tmp_pos(1,m) = real(i) + pdt_pos(1,l)
            tmp_pos(2,m) = real(j) + pdt_pos(2,l)
            tmp_pos(3,m) = real(k) + pdt_pos(3,l)
         enddo
      enddo
   enddo
enddo
do m=1,pdt_nsite
write(*,*) ' LL ', tmp_pos(:,m)
enddo
do m=natoms+1-pdt_nsite, natoms
write(*,*) ' TR ', tmp_pos(:,m)
enddo
!
! Map the old crystal onto the new one
!
! use : k index of closest atom
n = 0
loop_new: do i=1, natoms
   if(tmp_iscat(i)>=0) cycle loop_new          ! Already occupied
   k = 0
   dmin = 1.0E12
   dist = 1.0E12
   v = tmp_pos(:,i)
   loop_old:do m=1, nprior
      if( cr_iscat(m) <0) cycle loop_old          ! Already transfered
      u = cr_pos(:, m)
      dist = do_blen(lspace, u, v)                ! Calculate distance
      if(dist < dmin) then                 ! Found new shortest distance
         dmin = dist
         k    = m                          ! Archive site
      endif
   enddo loop_old
   if(dmin < EPS) then                     ! Found close atom
      tmp_iscat(   i) = cr_iscat(   k)
      tmp_prop (   i) = cr_prop (   k)
      tmp_mole (   i) = cr_mole (   k)
      tmp_surf (:, i) = cr_surf (:, k)
      tmp_magn (:, i) = cr_magn (:, k)
      tmp_pos  (:, i) = cr_pos  (:, k)
      cr_iscat(k) = -1                     ! Flag as done
      n = n + 1
      do j=mole_off(cr_mole(k)), mole_off(cr_mole(k))+mole_len(cr_mole(k))
         if(mole_cont(j) == k) then
            mole_cont(j) = i           ! Change atom number in Molecule
         endif
      enddo
   endif
enddo loop_new
write(*,*) 'Map ', n
write(*,*) 'OLD ', maxval(cr_iscat(1:cr_natoms))
write(*,*) 'NEW ', minval(tmp_iscat(1:natoms))
!
cr_iscat = 0
cr_prop  = 0
cr_mole  = 0
cr_surf  = 0
cr_magn  = 0.0
cr_pos   = 0.0
cr_iscat(  1:natoms) = tmp_iscat(  1:natoms)
cr_prop (  1:natoms) = tmp_prop (  1:natoms)
cr_mole (  1:natoms) = tmp_mole (  1:natoms)
cr_surf (:,1:natoms) = tmp_surf (:,1:natoms)
cr_magn (:,1:natoms) = tmp_magn (:,1:natoms)
cr_pos  (:,1:natoms) = tmp_pos  (:,1:natoms)
cr_natoms = natoms
!
cr_icc (1) = pdt_ihig(1) - pdt_ilow(1) + 1
cr_icc (2) = pdt_ihig(2) - pdt_ilow(2) + 1
cr_icc (3) = pdt_ihig(3) - pdt_ilow(3) + 1
cr_ncatoms = pdt_nsite
chem_purge = .FALSE.     ! Crystal dimension should allow periodic boundary
if(cr_icc(1)>1) chem_period(1) = .TRUE.
if(cr_icc(2)>1) chem_period(2) = .TRUE.
if(cr_icc(3)>1) chem_period(3) = .TRUE.
chem_quick = .TRUE.

!
deallocate(tmp_iscat)
deallocate(tmp_prop )
deallocate(tmp_mole )
deallocate(tmp_surf )
deallocate(tmp_magn )
deallocate(tmp_pos  )
!
end subroutine map_to_aver
!
!*******************************************************************************
!
subroutine perioditize_realloc(old, osize, nsize)
!-
! Reallocate an arrac
!+
implicit none
real(kind=PREC_SP), dimension(:,:), allocatable :: old
integer, intent(in) :: osize
integer, intent(in) :: nsize
!
real(kind=PREC_SP), dimension(:,:)  , allocatable :: tmp      ! Atom positions in average cell
!
allocate(tmp(3, nsize))
tmp(:,1:osize) = old(:,1:osize)
call move_alloc(tmp, old)
!
end subroutine perioditize_realloc
end module perioditize_mod
