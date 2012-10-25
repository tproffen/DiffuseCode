MODULE mole_surf_mod
!
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
   CALL find_surface_character
!
   END SUBROUTINE do_place_molecule
!
!
   SUBROUTINE find_surface_character
!
   USE crystal_mod
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
   INTEGER                 :: ianz  ! number of parameters in cpara
   LOGICAL  , DIMENSION(3) :: fp
   LOGICAL                 :: fq
   REAL     , DIMENSION(3) :: x
   REAL      :: rmin
   REAL      :: radius
!
      fp (1) = chem_period (1)
      fp (2) = chem_period (2)
      fp (3) = chem_period (3)
      fq = chem_quick
!
   DO ia=1,cr_natoms
!
     IF(IBITS(cr_prop(ia),PROP_SURFACE_EXT,1).eq.1 .and.        &  ! real Atom is near surface
        IBITS(cr_prop(ia),PROP_OUTSIDE    ,1).eq.0       ) THEN    ! real Atom is near surface
        write(*,*) cr_pos(:,ia), cr_prop(ia)
        cpara(1) = 'all'
        lpara(1) = 3
        ianz     = 1
        x(1)     = cr_pos(1,ia)
        x(2)     = cr_pos(2,ia)
        x(3)     = cr_pos(3,ia)
        rmin     = 0.5
        radius   = 5.5
        CALL get_iscat (ianz, cpara, lpara, werte, maxw, lnew)
        call do_find_env (ianz, werte, maxw, x, rmin,&
                          radius, fq, fp)
     ENDIF
!
   ENDDO
!
   END SUBROUTINE find_surface_character
!
END MODULE mole_surf_mod
