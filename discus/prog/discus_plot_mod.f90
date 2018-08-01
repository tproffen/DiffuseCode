MODULE discus_plot_mod
!+
!
!     variables needed for the plotting of structures
!-
USE discus_config_mod
!
SAVE
!
INTEGER, PRIVATE    :: ik
INTEGER                                   :: PL_MAXSCAT = 1
!
CHARACTER(LEN=1024)  :: pl_jmol  = ' '
CHARACTER(LEN= 200)  :: pl_out   = 'plot.cif'
CHARACTER(LEN=  80)  :: pl_title = 'crystal structure'
CHARACTER(LEN=  12)  :: pl_prog  = 'cif'
CHARACTER(LEN=   8)  :: pl_vrml  = 'u'
CHARACTER(LEN=   4)  :: pl_col   = 'xyz'
LOGICAL                                   :: pl_init      ! Plot was initialized
REAL                                      :: pl_width     = 1e12
REAL   , DIMENSION(3,2)                   :: pl_dim       = reshape((/(-1e14,ik=1,3),(1e12,ik=1,3)/),shape(pl_dim)) ! (3,2)
REAL   , DIMENSION(3)                     :: pl_hkl       = (/0.0,0.0,1.0/)
REAL   , DIMENSION(3)                     :: pl_uvw       = (/0.0,0.0,1.0/)
REAL   , DIMENSION(3)                     :: pl_vec       = 0.0
REAL   , DIMENSION(3)                     :: pl_abs       = (/1.0,0.0,0.0/)
REAL   , DIMENSION(3)                     :: pl_ord       = (/0.0,1.0,0.0/)
REAL   , DIMENSION(3,3)                   :: pl_mat       = &
         RESHAPE((/1.,(0.,0.,0.,1.,ik=1,2)/),SHAPE(pl_mat ))
REAL   , DIMENSION(3,3)                   :: pl_inv       = &
         RESHAPE((/1.,(0.,0.,0.,1.,ik=1,2)/),SHAPE(pl_inv ))
REAL   , DIMENSION(:),     ALLOCATABLE    :: pl_siz       ! (0:MAXSCAT)
REAL   , DIMENSION(:,:),   ALLOCATABLE    :: pl_rgb       ! (3,0:MAXSCAT)
INTEGER, DIMENSION(3)                     :: pl_back = &  ! plot background
         (/ 240, 240, 240 /)
REAL   , DIMENSION(:,:,:), ALLOCATABLE    :: pl_bond_len  ! (2,0:MAXSCAT,0:MAXSCAT)
REAL   , DIMENSION(  :,:), ALLOCATABLE    :: pl_bond_rad  ! (  0:MAXSCAT,0:MAXSCAT)
REAL   , DIMENSION(:,:,:), ALLOCATABLE    :: pl_bond_col  ! (3,0:MAXSCAT,0:MAXSCAT)
REAL                                      :: pl_vrml_scaling = 0.05
REAL   , DIMENSION(4,4)                   :: pl_tran_g    = &
         RESHAPE((/1.,(0.,0.,0.,0.,1.,ik=1,3)/),SHAPE(pl_tran_g ))
REAL   , DIMENSION(4,4)                   :: pl_tran_gi   = &
         RESHAPE((/1.,(0.,0.,0.,0.,1.,ik=1,3)/),SHAPE(pl_tran_gi))
REAL   , DIMENSION(4,4)                   :: pl_tran_f    = &
         RESHAPE((/1.,(0.,0.,0.,0.,1.,ik=1,3)/),SHAPE(pl_tran_f ))
REAL   , DIMENSION(4,4)                   :: pl_tran_fi   = &
         RESHAPE((/1.,(0.,0.,0.,0.,1.,ik=1,3)/),SHAPE(pl_tran_fi))
REAL   , DIMENSION(0:1)                   :: pl_scale     = (/-1.0, 1.0/)
INTEGER, DIMENSION(0:1)                   :: pl_sel_prop  = (/0,0/)
INTEGER, DIMENSION(:), ALLOCATABLE        :: pl_typ       ! (0:MAXSCAT)
INTEGER, DIMENSION(:), ALLOCATABLE        :: pl_color     ! (0:MAXSCAT)
LOGICAL, DIMENSION(:), ALLOCATABLE        :: pl_latom     ! (0:MAXSCAT)
LOGICAL, DIMENSION(:), ALLOCATABLE        :: pl_batom_a   ! (0:MAXSCAT)
LOGICAL, DIMENSION(:), ALLOCATABLE        :: pl_batom_e   ! (0:MAXSCAT)
INTEGER                                   :: pl_poly_n    ! Number of polyhedra definitions
LOGICAL, DIMENSION(:), ALLOCATABLE        :: pl_poly_c    ! (0:MAXSCAT)
LOGICAL, DIMENSION(:), ALLOCATABLE        :: pl_poly_o    ! (0:MAXSCAT)
REAL                                      :: pl_poly_dmin ! Minimum distance to neighbor for polhedra
REAL                                      :: pl_poly_dmax ! Maximum distance to neighbor for polhedra
INTEGER                                   :: pl_poly_nmin ! Minimum neighbors for polhedra
INTEGER                                   :: pl_poly_nmax ! Maximum neighbors for polhedra
LOGICAL                                   :: pl_poly_face ! Face style flat/collapsed
LOGICAL                                   :: pl_poly_hue  ! Face style solid / transparent
CHARACTER(LEN=128)                        :: pl_poly_col  ! Face color
LOGICAL                                   :: pl_dens      = .false.
LOGICAL                                   :: pl_sel_atom  = .true.
LOGICAL                                   :: pl_mol_all   = .true.
LOGICAL, DIMENSION(:,:), ALLOCATABLE      :: pl_bond      ! (0:MAXSCAT,0:MAXSCAT)
LOGICAL                                   :: pl_append    = .false.
LOGICAL, DIMENSION(3)                     :: pl_ext_all   = .true.
INTEGER                                   :: pl_size_of   = 0      ! Bytes allocates for plot
!
END MODULE discus_plot_mod
