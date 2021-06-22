MODULE do_find_top
!
CONTAINS
!
!*****7*****************************************************************
!
SUBROUTINE do_find (line, laenge) 
!-                                                                      
!     Finds the environment around an atom                              
!+                                                                      
USE discus_config_mod 
USE charact_mod 
USE charact_mod 
USE crystal_mod 
USE get_iscat_mod
USE chem_mod 
!
use do_find_mod
USE get_params_mod
USE ber_params_mod
USE errlist_mod 
USE get_params_mod
USE precision_mod
USE str_comp_mod
use take_param_mod
!                                                                       
IMPLICIT none 
!                                                                       
INTEGER , PARAMETER :: maxw = 200
INTEGER , PARAMETER :: mmaxw = 5 
!                                                                       
CHARACTER(len= * ) :: line 
CHARACTER(LEN=MAX(PREC_STRING,LEN(line))) cpara (maxw) 
CHARACTER(LEN=MAX(PREC_STRING,LEN(line))) ccpara (mmaxw) 
INTEGER lpara (maxw) 
INTEGER llpara (maxw) 
INTEGER i, ii, ianz, iianz, laenge 
integer :: iatom
integer, dimension(5) :: vect
LOGICAL lnew, fq, fp (3) 
REAL rmin 
REAL radius 
REAL(KIND=PREC_DP) :: werte (maxw) 
REAL(KIND=PREC_DP) :: wwerte (maxw) 
REAL x (3) 
!                                                                       
PARAMETER (lnew = .false.) 
!
integer, parameter :: NOPTIONAL = 1
integer, parameter :: O_MODE    = 1
character(LEN=   4), dimension(NOPTIONAL) :: oname   !Optional parameter names
character(LEN=PREC_STRING), dimension(NOPTIONAL) :: opara   !Optional parameter strings returned
integer            , dimension(NOPTIONAL) :: loname  !Lenght opt. para name
integer            , dimension(NOPTIONAL) :: lopara  !Lenght opt. para name returned
logical            , dimension(NOPTIONAL) :: lpresent!opt. para is present
real(kind=PREC_DP) , dimension(NOPTIONAL) :: owerte   ! Calculated values
integer, parameter                        :: ncalc = 0 ! Number of values to calculate 
!
data oname  / 'mode' /
data loname /  4     /
opara  =  (/ 'single'/)   ! Always provide fresh default values
lopara =  (/  6      /)
owerte =  (/  0.0    /)
!
!                                                                       
!                                                                       
fp (1) = chem_period (1) 
fp (2) = chem_period (2) 
fp (3) = chem_period (3) 
fq = chem_quick 
!                                                                       
CALL get_params (line, ianz, cpara, lpara, maxw, laenge) 
IF (ier_num.eq.0) then 
!
   IF(str_comp(cpara(1), 'env', 1, lpara(1), 3) ) then 
!                                                                       
!     ----Find environment                                              
!                                                                       
      IF (ianz.ge.7) then 
!                                                                       
!     ------copy last five parameters for evaluation                    
!                                                                       
         DO i = ianz - 4, ianz 
            ii = i - (ianz - 5) 
            ccpara (ii) = cpara (i) 
            llpara (ii) = lpara (i) 
         ENDDO 
         iianz = 5 
         CALL ber_params (iianz, ccpara, llpara, wwerte, mmaxw) 
         x (1) = wwerte (1) 
         x (2) = wwerte (2) 
         x (3) = wwerte (3) 
         rmin = wwerte (4) 
         radius = wwerte (5) 
         IF (ier_num.eq.0) then 
!                                                                       
!     -------- shift remaining parameters one left                      
!                                                                       
            DO i = 2, ianz - 5 
               cpara (i - 1) = cpara (i) 
               lpara (i - 1) = lpara (i) 
            ENDDO 
            ianz = ianz - 6 
!                                                                       
!     --------Get scattering curves                                     
!                                                                       
            CALL get_iscat(ianz, cpara, lpara, werte, maxw, lnew) 
            IF (ier_num.eq.0) then 
               CALL do_find_env(ianz, werte, maxw, x, rmin, radius, fq, fp)
            ENDIF 
         ENDIF 
      ELSE 
         ier_num = - 6 
         ier_typ = ER_COMM 
      ENDIF 
   ELSEIF(str_comp(cpara(1) , 'menv', 1, lpara(1), 4) ) then 
!                                                                       
!     ----Find molecular environment                                    
!                                                                       
      IF (ianz.ge.7) then 
!                                                                       
!     ------copy last three parameters for evaluation                   
!                                                                       
         DO i = ianz - 4, ianz 
            ii = i - (ianz - 5) 
            ccpara (ii) = cpara (i) 
            llpara (ii) = lpara (i) 
         ENDDO 
         iianz = 5 
         CALL ber_params (iianz, ccpara, llpara, wwerte, mmaxw) 
         x (1) = wwerte (1) 
         x (2) = wwerte (2) 
         x (3) = wwerte (3) 
         rmin = wwerte (4) 
         radius = wwerte (5) 
         IF (ier_num.eq.0) then 
!                                                                       
!     -------- shift remaining parameters one left                      
!                                                                       
            DO i = 2, ianz - 5 
               cpara (i - 1) = cpara (i) 
               lpara (i - 1) = lpara (i) 
            ENDDO 
            ianz = ianz - 6 
!                                                                       
!     --------Get allowed molecule types                                
!                                                                       
            IF(str_comp(cpara(1), 'all', 1, lpara(1), 3) ) then
               ianz = - 1 
            ELSE 
               CALL ber_params (ianz, cpara, lpara, werte, maxw) 
            ENDIF 
            IF (ier_num.eq.0) then 
               CALL do_find_mol (ianz, werte, maxw, x, rmin, radius)
            ENDIF
         ENDIF 
      ELSE 
         ier_num = - 6 
         ier_typ = ER_COMM 
      ENDIF 
   elseif(str_comp(cpara(1), 'conn', 1, lpara(1), 4) ) then 
      call get_optional(ianz, MAXW, cpara, lpara, NOPTIONAL,  ncalc, &
                        oname, loname, opara, lopara, lpresent, owerte)
      if(ianz>=2) then
         call del_params(1, ianz, cpara, lpara, maxw)
         if(ier_num==0) then
            if(ianz>=3) then
               if(str_comp(cpara(ianz), 'all', 3, lpara(ianz), 3)) then    ! Flag that all connectivities are to be used
                  cpara(ianz) = '-1'
                  lpara(ianz) =  2
               endif
            endif
            call ber_params (ianz, cpara, lpara, werte, maxw) 
            if(ier_num==0) then
               if(opara(O_MODE) == 'single' .and. ianz==2) then
                  call do_find_neig_conn(ianz, MAXW, werte)   ! Get conn's for atom werte(1)
               elseif(opara(O_MODE) == 'rec' .and. ianz==1) then
                  call do_find_neig_conn_all(ianz, MAXW, werte)   ! Get recursive conn's for atom werte(1)
               else 
                  ier_num = - 6 
                  ier_typ = ER_COMM 
               endif 
            endif 
         endif 
      ELSE 
         ier_num = - 6 
         ier_typ = ER_COMM 
      ENDIF 
   elseif(str_comp(cpara(1), 'neig', 1, lpara(1), 4) ) then 
      if(ianz==3) then
         call del_params(1, ianz, cpara, lpara, maxw)
         if(ier_num==0) then
            call ber_params (ianz, cpara, lpara, werte, maxw) 
            if(ier_num==0) then
               call do_find_neig(nint(werte(1)), nint(werte(2)))
            endif 
         endif 
      ELSE 
         ier_num = - 6 
         ier_typ = ER_COMM 
      ENDIF 
   elseif(str_comp(cpara(1), 'vect', 1, lpara(1), 4) ) then 
      if(ianz==7) then
         call del_params(1, ianz, cpara, lpara, maxw)
         if(ier_num==0) then
            call ber_params (ianz, cpara, lpara, werte, maxw) 
            if(ier_num==0) then
               iatom = nint(werte(1))
               vect  = nint(werte(2:6))
               call do_find_nei_vect(iatom, vect)
            endif 
         endif 
      ELSE 
         ier_num = - 6 
         ier_typ = ER_COMM 
      ENDIF 
   ELSE 
      ier_num = - 6 
      ier_typ = ER_COMM 
   ENDIF 
ELSE 
   ier_num = - 10 
   ier_typ = ER_APPL 
ENDIF 
!                                                                       
END SUBROUTINE do_find                        
!
!*****7*****************************************************************
!
subroutine do_find_neig(jatom, ic)
!-
!  Find neighbors for a neighborhood as defined in mmc
!+
!
use atom_env_mod
use chem_neig_multi_mod
!
use precision_mod
!
implicit none
!
integer , intent(in) ::  jatom
integer , intent(in) ::  ic
!
integer, parameter :: MAX_CENT = 2
integer           , dimension(:, :)   , allocatable :: iatom ! indices of neighbs
real(kind=PREC_SP), dimension(:, :, :), allocatable :: patom ! Coordinates
logical           , dimension(:, :)   , allocatable :: tatom ! indices of neighbs
integer           , dimension(:)      , allocatable :: natom ! no of neigh
integer                                             :: ncent ! no of central atoms
integer :: i      ! Dummy index
integer :: icent  ! Fixe central atom to 1
!
allocate(iatom(  0:MAX_ATOM_ENV, MAX_CENT))
allocate(patom(3,0:MAX_ATOM_ENV, MAX_CENT))
allocate(tatom(  0:MAX_ATOM_ENV, MAX_CENT))
allocate(natom(                  MAX_CENT)) 
!
call chem_neighbour_multi(jatom, ic, iatom, patom, tatom, natom, ncent, MAX_ATOM_ENV, MAX_CENT)
icent = 1
do i=1, natom(icent)
   atom_env(i)   = iatom(  i,icent)
   atom_pos(:,i) = patom(:,i,icent)
enddo
!
deallocate(iatom)
deallocate(patom)
deallocate(tatom)
deallocate(natom)
!
end subroutine do_find_neig
!
!*****7*****************************************************************
!
subroutine do_find_nei_vect(iatom, vect)
!-
!  Find the neighbor under a given vector
!+
!
use crystal_mod
use atom_env_mod
use check_bound_mod
use chem_mod
use celltoindex_mod
!
use precision_mod
!
implicit none
!
integer              , intent(in) :: iatom    ! The central==start atom
integer, dimension(5), intent(in) :: vect     ! (isite, jsite, cx, cy, cz)
!
integer, dimension(3) :: icell                ! Unit cell number for start atom
integer, dimension(3) :: jcell                ! Unit cell number for end   atom
integer               :: isite                ! Site number for start atom
integer               :: jsite                ! Site number for end   atom
integer               :: jatom                ! RTarget atom number
logical               :: lok                  ! Flag if boundary is OK
real(kind=PREC_SP), dimension(3) :: offset    ! in case of periodic boundary conditions, the offset
!
call indextocell(iatom, icell, isite)         ! Determine unit cell location
!
atom_env(0) = 0
if(isite==vect(1)) then                       ! Start atom is at correct site
   jcell = icell + vect(3:5)                  ! Add unit cell shift to atom
   lok   = .false.
   call check_bound(jcell, offset, chem_period, lok)
   if(lok) then
      jsite = vect(2)                         ! Target site
      call celltoindex(jcell, jsite, jatom)   ! Determine unit cell location
      atom_env(0) = 1
      atom_env(1) = jatom
      atom_pos(:,1) = cr_pos(:,jatom) + offset
   endif
endif
!
end subroutine do_find_nei_vect
!
!*****7*****************************************************************
!
end MODULE do_find_top
