MODULE mole_surf_mod
!
USE crystal_mod
USE deco_mod
USE prop_para_mod
!
PRIVATE
PUBLIC do_place_molecule
!
CONTAINS
   SUBROUTINE do_place_molecule
!
!   Main menu and procedures to place molecules onto a surface
!
   USE modify_mod
   USE modify_func_mod
!
   IMPLICIT none
!
   include'doact.inc'
   include'errlist.inc'
   include'learn.inc'
   include'macro.inc'
   include'prompt.inc'
!
   INTEGER, PARAMETER   :: MAXW = 9
!
!  CHARACTER (LEN=*)   :: line
   CHARACTER (LEN=5)                       :: befehl! command on input line
   CHARACTER (LEN=50)                      :: prom  ! Menu prompt
   CHARACTER (LEN=1024)                    :: line  ! input line
   CHARACTER (LEN=1024)                    :: zeile ! remainder with parameters
   CHARACTER (LEN=1024), DIMENSION(1:MAXW) :: cpara ! parameter array
   INTEGER             , DIMENSION(1:MAXW) :: lpara ! length of parameters
   INTEGER                                 :: indxg ! location of "="
   INTEGER                                 :: lp    ! lengtz of zeile
   INTEGER i, j, ii, ianz, m, l, k, laenge, lbef
   LOGICAL                                 :: lend  ! condition of EOF
   LOGICAL                                 :: lold =.true. ! condition of EOF
   REAL                , DIMENSION(1:MAXW) :: werte ! parameter values
!
   INTEGER, EXTERNAL :: len_str
   LOGICAL, EXTERNAL :: str_comp
!
   lend = .false.
!
!  IF INITIALIZATION is needed, set property flags to external surface
!
   IF (dc_init) THEN
      dc_sel_prop(0) = bit_set ( dc_sel_prop, 0, PROP_SURFACE_EXT, .true.)
      dc_sel_prop(1) = bit_set ( dc_sel_prop, 1, PROP_SURFACE_EXT, .true.)
      dc_use_conn    = 1 ! while developing!!!
      dc_n_molecules = 1 ! while developing!!!
      dc_init = .false.
   ENDIF
!
   main_loop: do
     CALL no_error
      prom = prompt (1:len_str (prompt) ) //'/deco'
      CALL get_cmd (line, laenge, befehl, lbef, zeile, lp, prom)
      no_err: IF (ier_num.eq.0) THEN
         no_com: IF (line /= ' '      .and. line(1:1) /= '#' .and.      &
             line /= char(13) .and. line(1:1) /= '!'        ) THEN
!                                                                       
!     ----search for "="                                                
!                                                                       
            indxg = index (line, '=')
            is_math: IF (indxg.ne.0.and.                            &
               .not. (str_comp (befehl, 'echo', 2, lbef, 4) ) .and. &
               .not. (str_comp (befehl, 'syst', 2, lbef, 4) ) .and. &
               .not. (str_comp (befehl, 'help', 2, lbef, 4)   .or.  &
                      str_comp (befehl, '?   ', 2, lbef, 4) ) ) THEN
!                                                                       
! ------evaluate an expression and assign the value to a variabble      
!                                                                       
               CALL do_math (line, indxg, laenge)
            ELSE
!                                                                       
!------ ----execute a macro file                                        
!                                                                       
              is_com: IF (befehl (1:1) .eq.'@') THEN 
                  IF (laenge.ge.2) THEN 
                     CALL file_kdo (line (2:laenge), laenge-1) 
                  ELSE 
                     ier_num = - 13 
                     ier_typ = ER_MAC 
                  ENDIF 
!                                                                       
!     ----continues a macro 'continue'                                  
!                                                                       
              ELSEIF (str_comp (befehl, 'continue', 5, lbef, 8) ) THEN 
                  CALL macro_continue (zeile, lp) 
!                                                                       
!     ----Echo a string, just for interactive check in a macro 'echo'   
!                                                                       
              ELSEIF (str_comp (befehl, 'echo', 2, lbef, 4) ) THEN 
                  CALL echo (zeile, lp) 
!                                                                       
!      ---Evaluate an expression, just for interactive check 'eval'     
!                                                                       
              ELSEIF (str_comp (befehl, 'eval', 2, lbef, 4) ) THEN 
                  CALL do_eval (zeile, lp) 
!                                                                       
!     ----exit 'exit'                                                   
!                                                                       
              ELSEIF (str_comp (befehl, 'exit', 2, lbef, 4) ) THEN 
                  lend = .true. 
                  EXIT main_loop
!                                                                       
!     ----help 'help' , '?'                                             
!                                                                       
              ELSEIF (str_comp (befehl, 'help', 1, lbef, 4) .or.  &
                       str_comp (befehl, '?   ', 1, lbef, 4) ) THEN                                      
                  IF (str_comp (zeile, 'errors', 2, lp, 6) ) THEN 
                     lp = lp + 7 
                     CALL do_hel ('discus '//zeile, lp) 
                  ELSE 
                     lp = lp + 12 
                     CALL do_hel ('discus deco '//zeile, lp) 
                  ENDIF 
!                                                                       
!------- -Operating System Kommandos 'syst'                             
!                                                                       
              ELSEIF (str_comp (befehl, 'syst', 2, lbef, 4) ) THEN 
                  IF (zeile.ne.' '.and.zeile.ne.char (13) ) THEN 
                     CALL do_operating (zeile (1:lp), lp) 
                  ELSE 
                     ier_num = - 6 
                     ier_typ = ER_COMM 
                  ENDIF 
!                                                                       
!------  -----waiting for user input                                    
!                                                                       
              ELSEIF (str_comp (befehl, 'wait', 3, lbef, 4) ) THEN 
                  CALL do_input (zeile, lp) 
!
!     ----Original decorate commands                                      
!
!                                                                       
!------ --Handle property settings 'property'                           
!                                                                       
              ELSEIF (str_comp (befehl, 'property', 4, lbef, 8) ) then
!                                                                       
                  CALL property_select (zeile, lp, dc_sel_prop)
!
              ELSEIF (str_comp (befehl, 'run', 3, lbef, 3)) THEN
                  CALL deco_run
!
              ELSEIF (str_comp (befehl, 'sel', 3, lbef, 3) .or.  &
                      str_comp (befehl, 'del', 3, lbef, 3)     ) THEN
                  CALL atom_select (zeile, lp, 0,  DC_MAXSCAT,  dc_latom, &
                  dc_sel_atom , lold  ,   &
                  str_comp (befehl, 'sel', 2, lbef, 3) )
!
!
              ELSEIF (str_comp (befehl, 'show', 3, lbef, 4)) THEN
                  CALL deco_show
!
                  CALL find_surface_character
!
              ELSE is_com
                 ier_num = -8
                 ier_typ = ER_COMM
!
              ENDIF is_com ! END IF BLOCK actual commands
           ENDIF is_math   ! END IF BLOCK math equation or specific command
        ENDIF no_com       ! END IF BLOCK no comment
      ENDIF no_err         ! END IF BLOCK no error reading input
!
   ENDDO main_loop     ! END DO main loop of menu
!
   END SUBROUTINE do_place_molecule
!
!######################################################################
!
   SUBROUTINE find_surface_character
!
   USE crystal_mod
   USE atom_env_mod
   USE chem_mod
   USE modify_mod
!
   IMPLICIT NONE
!
   INTEGER, PARAMETER :: MAXW = 1
   LOGICAL, PARAMETER :: LNEW = .false.
!
   CHARACTER  (LEN=4), DIMENSION(1:MAXW) :: cpara
   INTEGER  ,          DIMENSION(1:MAXW) :: lpara
   REAL     ,          DIMENSION(1:MAXW) :: werte
   INTEGER   :: ia ! dummy index
   INTEGER   :: ib ! dummy index
   INTEGER   :: istart,iend ! dummy index
   INTEGER                 :: ianz  ! number of parameters in cpara
   INTEGER                 :: istat
   LOGICAL  , DIMENSION(3) :: fp
   LOGICAL                 :: fq
   REAL     , DIMENSION(3) :: x
   REAL      :: rmin
   REAL      :: radius
   REAL                      :: sigma
   REAL     , DIMENSION(3,3) :: direct
   REAL     , DIMENSION(3,3) :: recipr
   REAL     , DIMENSION(3  ) :: vector
   REAL     , DIMENSION(3  ) :: finalv
   REAL     , DIMENSION(:), ALLOCATABLE :: deviat  ! distances to plane
!
!
      fp (1) = chem_period (1)
      fp (2) = chem_period (2)
      fp (3) = chem_period (3)
      fq = chem_quick
      ALLOCATE(deviat(0:MAX_ATOM_ENV), STAT = istat)
!
ia = 184
istart = 1
iend   = cr_natoms
!istart = 9211
!iend   = 9211
   atomloop: DO ia=istart,iend
!
     IF(IBITS(cr_prop(ia),PROP_SURFACE_EXT,1).eq.1 .and.        &  ! real Atom is near surface
        IBITS(cr_prop(ia),PROP_OUTSIDE    ,1).eq.0       ) THEN    ! real Atom is near surface
write(*,*)
write(*,*) ' Central : ',ia, cr_iscat(ia),cr_pos(:,ia), cr_prop(ia)
        cpara(1) = 'all'
        lpara(1) = 3
        ianz     = 1
        x(1)     = cr_pos(1,ia)
        x(2)     = cr_pos(2,ia)
        x(3)     = cr_pos(3,ia)
        rmin     = 0.0
        radius   = 6.5
        CALL get_iscat (ianz, cpara, lpara, werte, maxw, lnew)
        ianz = 1
        werte(1) = cr_iscat(ia)
        call do_find_env (ianz, werte, maxw, x, rmin,&
                          radius, fq, fp)
        IF(atom_env(0) > 3 ) THEN
        direct = 0.0
        vector = 0.0
        DO ib = 1, atom_env(0) 
           IF(IBITS(cr_prop(atom_env(ib)),PROP_SURFACE_EXT,1).eq.1 ) THEN          ! real Atom is near surface
           WRITE(*,*) '   neig:',atom_env(ib),cr_iscat(atom_env(ib)),' at ',atom_pos(:,ib), cr_prop(atom_env(ib))
           direct(1,1) = direct(1,1) + atom_pos(1,ib)**2             ! Accumulate x^2 etc
           direct(2,2) = direct(2,2) + atom_pos(2,ib)**2
           direct(3,3) = direct(3,3) + atom_pos(3,ib)**2
           direct(1,2) = direct(1,2) + atom_pos(1,ib)*atom_pos(2,ib)
           direct(1,3) = direct(1,3) + atom_pos(1,ib)*atom_pos(3,ib)
           direct(2,3) = direct(2,3) + atom_pos(2,ib)*atom_pos(3,ib)
           vector(1  ) = vector(1  ) + atom_pos(1,ib)
           vector(2  ) = vector(2  ) + atom_pos(2,ib)
           vector(3  ) = vector(3  ) + atom_pos(3,ib)
           ENDIF
        ENDDO
        direct(2,1) = direct(1,2) ! copy into transposed elements
        direct(3,1) = direct(1,3)
        direct(3,2) = direct(2,3)
!
        write(*,*) ' row 1 ',direct(1,:), ' V ',vector(1)
        write(*,*) ' row 2 ',direct(2,:), ' V ',vector(2)
        write(*,*) ' row 3 ',direct(3,:), ' V ',vector(3)
        call invmat( recipr, direct )  ! calculate inverse matrix
        write(*,*) ' row 1 ',recipr(1,:)
        write(*,*) ' row 2 ',recipr(2,:)
        write(*,*) ' row 3 ',recipr(3,:)
        finalv      = MATMUL(recipr, vector) ! finalv contains the hkl of the plane
        sigma = 0.0
        DO ib = 1, atom_env(0) 
     IF(IBITS(cr_prop(atom_env(ib)),PROP_SURFACE_EXT,1).eq.1 ) THEN          ! real Atom is near surface
              deviat(ib) = (atom_pos(1,ib)*finalv(1)+      &  ! accumulate the distances
                            atom_pos(2,ib)*finalv(2)+      &
                            atom_pos(3,ib)*finalv(3)-1)
              sigma      = sigma + deviat(ib)**2              ! calculate the sigma of the distance distribution
write(*,*) ' neigbor ',ib,' deviates by ',deviat(ib)
     ENDIF
        ENDDO
        sigma = SQRT(sigma)
        write(*,*) ' Normal ',finalv,' sigma ',sigma
        IF(sigma <    0.05) THEN
           cr_iscat(ia) = cr_iscat(ia) + 4 !cr_nscat
           CYCLE atomloop
        ENDIF
        ENDIF
!
!       Atom is not in plane, try a row
!
     ENDIF
!
   ENDDO atomloop
!
   DEALLOCATE(deviat, STAT = istat)
!
   END SUBROUTINE find_surface_character
!
!######################################################################
!
   SUBROUTINE deco_run
!
!  Performs the actual decoration
!
   USE conn_mod
   USE modify_func_mod
!
   IMPLICIT none
!
   INTEGER, PARAMETER         :: MAXW = 2000  ! not ideal, should be dynamic ....
!
   INTEGER   :: istatus          ! status
   INTEGER   :: i  ! dummy index
   INTEGER   :: ia ! dummy index
   INTEGER   :: ib ! dummy index
   INTEGER   :: istart,iend ! dummy index
   INTEGER   :: is ! scattering number of surface atom
   INTEGER   :: idef ! connectivity definition number
   INTEGER, DIMENSION(1:MAXW) :: c_list
   INTEGER   :: natoms                 ! number of atoms in connectivity
!
integer ::itype
real, dimension(3) :: posit
integer iprop
!
!  Load the molecules into temporary structures to reduce disk I/O
!
write(*,*) ' LOADING THE MOLECULES'
   CALL deco_get_molecules
write(*,*) ' Content of molecules'
DO i=1, dc_n_molecules
   DO ia=1, dc_molecules(i)%get_natoms()
   CALL dc_molecules(i)%get_cryst_atom(ia, itype, posit, iprop)
   write(*,*) itype, posit,iprop
   ENDDO
ENDDO

!
   istart = 1
   iend   = cr_natoms
write(*,*) ' STARTING MAIN LOOP', istart, iend
   main_loop: DO ia=istart,iend
     is_sel: IF(check_select_status (dc_latom (cr_iscat (ia) ), cr_prop (ia),   &
                                     dc_sel_prop)              ) THEN
!write(*,*)
!write(*,*) ' Central : ',ia, cr_iscat(ia),cr_pos(:,ia), cr_prop(ia)
        is   = cr_iscat(ia)                  ! get scattering type central
        idef = dc_use_conn(is)               ! use this connectivity list for surface
        CALL get_connectivity_list (ia, is, idef, maxw, c_list, natoms )
     ENDIF is_sel ! END IF BLOCK is selected
   ENDDO main_loop   ! END DO main loop over all atoms
write(*,*) 'CLEANING UP'
!
!  Clean up temporary arrays
!
   DO i=1, dc_n_molecules
      CALL dc_molecules(i)%finalize_atoms()
   ENDDO
   DEALLOCATE(dc_molecules)
   DEALLOCATE(m_ntypes, STAT=istatus)
   DEALLOCATE(m_length, STAT=istatus)
!
   END SUBROUTINE deco_run
!
!######################################################################
   SUBROUTINE deco_get_molecules
!
   USE structur, ONLY: test_file
   IMPLICIT none
!
   CHARACTER (LEN=1024)                :: strufile
!
   INTEGER  :: i                ! dummy index
   INTEGER  :: istatus          ! status
   INTEGER  :: natoms           ! number of atoms in file
   INTEGER  :: ntypes           ! number of atom types in file
   INTEGER  :: n_mole           ! number of molecules
   INTEGER  :: n_type           ! number of molecule types
   INTEGER  :: n_atom           ! number of atoms in molecules
   INTEGER  :: init   = -1      ! Initialize atom names/types
   LOGICAL  :: lcell  = .true.  ! Treat atoms with equal name and B as one type
!
   ALLOCATE(m_ntypes(dc_n_molecules), STAT=istatus)
   ALLOCATE(m_length(dc_n_molecules), STAT=istatus)
!
   ALLOCATE(dc_molecules(dc_n_molecules), STAT = istatus)
!
   DO i=1,dc_n_molecules        ! load all molecules
      strufile = dc_input(i)
write(*,*) ' STRUCFILE ',strufile
      CALL test_file(strufile, natoms, ntypes, n_mole, n_type, &
                     n_atom, init, lcell)
      CALL dc_molecules(i)%alloc_arrays(natoms, ntypes, n_mole, &
                     n_type, n_atom)
!      CALL read_crystal ( dc_molecules(i), strufile )
      m_length(i) = natoms
      m_ntypes(i) = ntypes
   ENDDO
!
   END SUBROUTINE deco_get_molecules
!
!######################################################################
!
   SUBROUTINE deco_show
!
   IMPLICIT none
!
   END SUBROUTINE deco_show
!
!******************************************************************************
!   SUBROUTINE read_crystal ( this, infile)
!!
!!  Read a crystal structure from file
!!  This procedure interfaces to the old "reastru" in "structur.f90"
!
!   USE inter_readstru
!   USE structur, ONLY: readstru
!!
!   IMPLICIT none
!!
!   TYPE (cl_cryst)                  :: this   ! The current crystal
!   CHARACTER (LEN=*   ), INTENT(IN)  :: infile ! Disk file
!   INTEGER                           :: inum   ! dummy index
!   REAL   , DIMENSION(3)             :: posit  ! dummy position vector
!   INTEGER                           :: istat  ! status variable
!!
!   rd_strucfile = infile
!   rd_NMAX      = this%cr_natoms
!   rd_MAXSCAT   = this%cr_nscat
!!
!!  Allocate temporary arrays
!!
!write(*,*) ALLOCATED(rd_cr_dw    )
!write(*,*) ALLOCATED(rd_cr_at_lis)
!   ALLOCATE ( rd_cr_dw    (0:rd_MAXSCAT)    , STAT = istat )
!   ALLOCATE ( rd_cr_at_lis(0:rd_MAXSCAT)    , STAT = istat )
!   ALLOCATE ( rd_cr_pos   (1:3,1:rd_NMAX)   , STAT = istat )
!   ALLOCATE ( rd_cr_iscat (1:rd_NMAX)       , STAT = istat )
!   ALLOCATE ( rd_cr_prop  (1:rd_NMAX)       , STAT = istat )
!   ALLOCATE ( rd_as_at_lis(0:rd_MAXSCAT)    , STAT = istat )
!   ALLOCATE ( rd_as_dw    (0:rd_MAXSCAT)    , STAT = istat )
!   ALLOCATE ( rd_as_pos   (1:3,1:rd_MAXSCAT), STAT = istat )
!   ALLOCATE ( rd_as_iscat (1:rd_MAXSCAT)    , STAT = istat )
!   ALLOCATE ( rd_as_prop  (1:rd_MAXSCAT)    , STAT = istat )
!write(*,*) 'STATUS ',istat
!!
!   rd_cr_dw   = 0.0
!   rd_cr_at_lis = ' '
!   rd_cr_pos    = 0.0
!   rd_cr_iscat  = 0
!   rd_cr_prop   = 0
!   rd_as_at_lis = ' '
!   rd_as_dw     = 0.0
!   rd_as_pos    = 0.0
!   rd_as_iscat  = 0
!   rd_as_prop   = 0
!!
!!WRITE(*,*) ' 0   ', rd_cr_iscat(1)
!!WRITE(*,*) ' 0   ', rd_cr_iscat(1), rd_cr_pos(1,1)
!!WRITE(*,*) ' 0   ', rd_cr_iscat(1), rd_cr_pos(1,1),rd_cr_pos(2,1)
!!WRITE(*,*) ' 0   ', rd_cr_iscat(1), rd_cr_pos(1,1),rd_cr_pos(2,1), rd_cr_pos(3,1)
!!WRITE(*,*) ' 0   ', rd_cr_iscat(1), rd_cr_pos(1,1),rd_cr_pos(2,1), rd_cr_pos(3,1), rd_cr_prop(1)
!!write(*,*) ' IN read_crystal ', rd_NMAX, rd_MAXSCAT
!   CALL readstru (rd_NMAX, rd_MAXSCAT, rd_strucfile, rd_cr_name,        &
!               rd_cr_spcgr, rd_cr_a0, rd_cr_win, rd_cr_natoms, rd_cr_nscat, rd_cr_dw,     &
!               rd_cr_at_lis, rd_cr_pos, rd_cr_iscat, rd_cr_prop, rd_cr_dim, rd_as_natoms, &
!               rd_as_at_lis, rd_as_dw, rd_as_pos, rd_as_iscat, rd_as_prop, rd_sav_ncell,  &
!               rd_sav_r_ncell, rd_sav_ncatoms, rd_spcgr_ianz, rd_spcgr_para)
!write(*,*) ' BACK '
!WRITE(*,*) ' nmax', rd_NMAX, rd_MAXSCAT
!write(*,*) ' rd_cr_name  ', rd_cr_name
!write(*,*) ' rd_cr_spcgr ', rd_cr_spcgr
!write(*,*) ' rd_cr_a0    ', rd_cr_a0   
!write(*,*) ' rd_cr_win   ', rd_cr_win  
!
!write(*,*) ' spcgr_ianz   ',rd_spcgr_ianz
!write(*,*) ' spcgr_para   ',rd_spcgr_para
!
!write(*,*) ' cr_at_lis ', rd_cr_at_lis
!WRITE(*,*) ' 1   ', rd_cr_iscat(1)
!WRITE(*,*) ' 1   ', rd_cr_pos(1,1)
!WRITE(*,*) ' 1   ', rd_cr_pos(2,1)
!WRITE(*,*) ' 1   ', rd_cr_pos(3,1)
!WRITE(*,*) ' 1   ', rd_cr_prop(1)
!write(*,*)
!   DO inum=1, rd_NMAX
!     posit = rd_cr_pos(:,inum)
!write(*,*) rd_cr_iscat(inum), posit, rd_cr_prop(inum)
!     CALL this%atoms(inum)%set_atom ( rd_cr_iscat(inum), posit, rd_cr_prop(inum) )
!   ENDDO
!!
!   DEALLOCATE ( rd_cr_dw    , STAT = istat )
!   DEALLOCATE ( rd_cr_at_lis, STAT = istat )
!   DEALLOCATE ( rd_cr_pos   , STAT = istat )
!   DEALLOCATE ( rd_cr_iscat , STAT = istat )
!   DEALLOCATE ( rd_cr_prop  , STAT = istat )
!   DEALLOCATE ( rd_as_at_lis, STAT = istat )
!   DEALLOCATE ( rd_as_dw    , STAT = istat )
!   DEALLOCATE ( rd_as_pos   , STAT = istat )
!   DEALLOCATE ( rd_as_iscat , STAT = istat )
!   DEALLOCATE ( rd_as_prop  , STAT = istat )
!!
!   END SUBROUTINE read_crystal
END MODULE mole_surf_mod
