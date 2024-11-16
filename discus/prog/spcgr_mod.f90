MODULE spcgr_mod
!+
!     Variables needed for the spacegroups
!-
PUBLIC
PRIVATE i
SAVE
!
INTEGER, PARAMETER  ::  SPCGR_MAX  =  315
!
character(len=12), dimension(9), parameter  ::  space_group_crystal_system = &
   (/'triclinic   ', 'monoclinic  ', 'monoclinic  ', 'orthorhombic', 'tetragonal  ',  &
     'trigonal    ', 'trigonal    ', 'hexagonal   ', 'cubic       ' /)
CHARACTER(LEN=16), DIMENSION(SPCGR_MAX)     ::  spcgr_name      ! Normal space group name
CHARACTER(LEN=16), DIMENSION(16:74    ,1:6) ::  spcgr_name_set  ! alternative setting in orthorhombic only
CHARACTER(LEN=16), DIMENSION(SPCGR_MAX)     ::  spcgr_point     ! Point group name
CHARACTER(LEN=16), DIMENSION(SPCGR_MAX)     ::  spcgr_laue      ! Laue  group name
INTEGER          , DIMENSION(SPCGR_MAX,2)   ::  spcgr_num       ! Space group number
INTEGER          , DIMENSION(SPCGR_MAX)     ::  spcgr_syst      ! Crystal system to which spce group belongs to
integer          , dimension(SPCGR_MAX)     ::  spcgr_num_sym   ! Number of symmetry element each spacgroup
integer          , dimension(SPCGR_MAX)     ::  spcgr_origin    ! Origin choice                   spacgroup
INTEGER                                     :: i
!
!                                   abc           bac           cab           cba           bca           acb
DATA (spcgr_name_set(16,i),i=1,6) /'P 2 2 2   ', 'P 2 2 2   ', 'P 2 2 2   ', 'P 2 2 2   ', 'P 2 2 2   ', 'P 2 2 2   '/
DATA (spcgr_name_set(17,i),i=1,6) /'P 2 2 21  ', 'P 2 2 21  ', 'P 21 2 2  ', 'P 21 2 2  ', 'P 2 21 2  ', 'P 2 21 2  '/
DATA (spcgr_name_set(18,i),i=1,6) /'P 21 21 2 ', 'P 21 21 2 ', 'P 2 21 21 ', 'P 2 21 21 ', 'P 21 2 21 ', 'P 21 2 21 '/
DATA (spcgr_name_set(19,i),i=1,6) /'P 21 21 21', 'P 21 21 21', 'P 21 21 21', 'P 21 21 21', 'P 21 21 21', 'P 21 21 21'/
DATA (spcgr_name_set(20,i),i=1,6) /'C 2 2 21  ', 'C 2 2 21  ', 'A 21 2 2  ', 'A 21 2 2  ', 'B 2 21 2  ', 'B 2 21 2  '/
DATA (spcgr_name_set(21,i),i=1,6) /'C 2 2 2   ', 'C 2 2 2   ', 'A 2 2 2   ', 'A 2 2 2   ', 'B 2 2 2   ', 'B 2 2 2   '/
DATA (spcgr_name_set(22,i),i=1,6) /'F 2 2 2   ', 'F 2 2 2   ', 'F 2 2 2   ', 'F 2 2 2   ', 'F 2 2 2   ', 'F 2 2 2   '/
DATA (spcgr_name_set(23,i),i=1,6) /'I 2 2 2   ', 'I 2 2 2   ', 'I 2 2 2   ', 'I 2 2 2   ', 'I 2 2 2   ', 'I 2 2 2   '/
DATA (spcgr_name_set(24,i),i=1,6) /'I 21 21 21', 'I 21 21 21', 'I 21 21 21', 'I 21 21 21', 'I 21 21 21', 'I 21 21 21'/
!
!                                   a b c         b a c         c a b         c b a         b c a         a c b 
DATA (spcgr_name_set(25,i),i=1,6) /'P m m 2   ', 'P m m 2   ', 'P 2m m    ', 'P 2m m    ', 'P m 2m    ', 'P m 2m    '/
DATA (spcgr_name_set(26,i),i=1,6) /'P m c 21  ', 'P c m 21  ', 'P 21 m a  ', 'P 21 a m  ', 'P b 21 m  ', 'P m 21 b  '/
DATA (spcgr_name_set(27,i),i=1,6) /'P c c 2   ', 'P c c 2   ', 'P 2a a    ', 'P 2a a    ', 'P b 2b    ', 'P b 2b    '/
DATA (spcgr_name_set(28,i),i=1,6) /'P m a 2   ', 'P b m 2   ', 'P 2m b    ', 'P 2c m    ', 'P c 2m    ', 'P m 2a    '/
DATA (spcgr_name_set(29,i),i=1,6) /'P c a 21  ', 'P b c 21  ', 'P 21 a b  ', 'P 21 c a  ', 'P c 21 b  ', 'P b 21 a  '/
DATA (spcgr_name_set(30,i),i=1,6) /'P n c 2   ', 'P c n 2   ', 'P 2n a    ', 'P 2a n    ', 'P b 2n    ', 'P n 2b    '/
DATA (spcgr_name_set(31,i),i=1,6) /'P m n 21  ', 'P n m 21  ', 'P 21 m n  ', 'P 21 n m  ', 'P n 21 m  ', 'P m 21 n  '/
DATA (spcgr_name_set(32,i),i=1,6) /'P b a 2   ', 'P b a 2   ', 'P 2c b    ', 'P 2c b    ', 'P c 2a    ', 'P c 2a    '/
DATA (spcgr_name_set(33,i),i=1,6) /'P n a 21  ', 'P b n 21  ', 'P 21 n b  ', 'P 21 c n  ', 'P c 21 n  ', 'P n 21 a  '/
DATA (spcgr_name_set(34,i),i=1,6) /'P n n 2   ', 'P n n 2   ', 'P 2n n    ', 'P 2n n    ', 'P n 2n    ', 'P n 2n    '/
DATA (spcgr_name_set(35,i),i=1,6) /'C m m 2   ', 'C m m 2   ', 'A 2m m    ', 'A 2m m    ', 'B m 2m    ', 'B m 2m    '/
DATA (spcgr_name_set(36,i),i=1,6) /'C m c 21  ', 'C c m 21  ', 'A 21 m a  ', 'A 21 a m  ', 'B b 21 m  ', 'B m 21 b  '/
DATA (spcgr_name_set(37,i),i=1,6) /'C c c 2   ', 'C c c 2   ', 'A 2a a    ', 'A 2a a    ', 'B b 2b    ', 'B b 2b    '/
DATA (spcgr_name_set(38,i),i=1,6) /'A m m 2   ', 'B m m 2   ', 'B 2m m    ', 'C 2m m    ', 'C m 2m    ', 'A m 2m    '/
DATA (spcgr_name_set(39,i),i=1,6) /'A e m 2   ', 'B m e 2   ', 'B 2e m    ', 'C 2m e    ', 'C m 2e    ', 'A e 2m    '/
DATA (spcgr_name_set(40,i),i=1,6) /'A m a 2   ', 'B b m 2   ', 'B 2m b    ', 'C 2c m    ', 'C c 2m    ', 'A m 2a    '/
DATA (spcgr_name_set(41,i),i=1,6) /'A e a 2   ', 'B b e 2   ', 'B 2e b    ', 'C 2c e    ', 'C c 2e    ', 'A e 2a    '/
DATA (spcgr_name_set(42,i),i=1,6) /'F m m 2   ', 'F m m 2   ', 'F 2m m    ', 'F 2m m    ', 'F m 2m    ', 'F m 2m    '/
DATA (spcgr_name_set(43,i),i=1,6) /'F d d 2   ', 'F d d 2   ', 'F 2d d    ', 'F 2d d    ', 'F d 2d    ', 'F d 2d    '/
DATA (spcgr_name_set(44,i),i=1,6) /'I m m 2   ', 'I m m 2   ', 'I 2m m    ', 'I 2m m    ', 'I m 2m    ', 'I m 2m    '/
DATA (spcgr_name_set(45,i),i=1,6) /'I b a 2   ', 'I b a 2   ', 'I 2c b    ', 'I 2c b    ', 'I c 2a    ', 'I c 2a    '/
DATA (spcgr_name_set(46,i),i=1,6) /'I m a 2   ', 'I b m 2   ', 'I 2m b    ', 'I 2c m    ', 'I c 2m    ', 'I m 2a    '/
!
!                                   a b c         b a c         c a b         c b a         b c a         a c b 
DATA (spcgr_name_set(47,i),i=1,6) /'P m m m   ', 'P m m m   ', 'P m m m   ', 'P m m m   ', 'P m m m   ', 'P m m m   '/
DATA (spcgr_name_set(48,i),i=1,6) /'P n n n   ', 'P n n n   ', 'P n n n   ', 'P n n n   ', 'P n n n   ', 'P n n n   '/
DATA (spcgr_name_set(49,i),i=1,6) /'P c c m   ', 'P c c m   ', 'P m a a   ', 'P m a a   ', 'P b m b   ', 'P b m b   '/
DATA (spcgr_name_set(50,i),i=1,6) /'P b a n   ', 'P b a n   ', 'P n c b   ', 'P n c b   ', 'P c n a   ', 'P c n a   '/
DATA (spcgr_name_set(51,i),i=1,6) /'P m m a   ', 'P m m b   ', 'P b m m   ', 'P c m m   ', 'P m c m   ', 'P m a m   '/
DATA (spcgr_name_set(52,i),i=1,6) /'P n n a   ', 'P n n b   ', 'P b n n   ', 'P c n n   ', 'P n c n   ', 'P n a n   '/
DATA (spcgr_name_set(53,i),i=1,6) /'P m n a   ', 'P n m b   ', 'P b m n   ', 'P c n m   ', 'P n c m   ', 'P m a n   '/
DATA (spcgr_name_set(54,i),i=1,6) /'P c c a   ', 'P c c b   ', 'P b a a   ', 'P c a a   ', 'P b c b   ', 'P b a b   '/
DATA (spcgr_name_set(55,i),i=1,6) /'P b a m   ', 'P b a m   ', 'P m c b   ', 'P m c b   ', 'P c m a   ', 'P c m a   '/
DATA (spcgr_name_set(56,i),i=1,6) /'P c c n   ', 'P c c n   ', 'P n a a   ', 'P n a a   ', 'P b n b   ', 'P b n b   '/
DATA (spcgr_name_set(57,i),i=1,6) /'P b c m   ', 'P c a m   ', 'P m c a   ', 'P m a b   ', 'P b m a   ', 'P c m b   '/
DATA (spcgr_name_set(58,i),i=1,6) /'P n n m   ', 'P n n m   ', 'P m n n   ', 'P m n n   ', 'P n m n   ', 'P n m n   '/
DATA (spcgr_name_set(59,i),i=1,6) /'P m m n   ', 'P m m n   ', 'P n m m   ', 'P n m m   ', 'P m n m   ', 'P m n m   '/
DATA (spcgr_name_set(60,i),i=1,6) /'P b c n   ', 'P c a n   ', 'P n c a   ', 'P n a b   ', 'P b n a   ', 'P c n b   '/
DATA (spcgr_name_set(61,i),i=1,6) /'P b c a   ', 'P c a b   ', 'P b c a   ', 'P c a b   ', 'P b c a   ', 'P c a b   '/
DATA (spcgr_name_set(62,i),i=1,6) /'P n m a   ', 'P m n b   ', 'P b n m   ', 'P c m n   ', 'P m c n   ', 'P n a m   '/
DATA (spcgr_name_set(63,i),i=1,6) /'C m c m   ', 'C c m m   ', 'A m m a   ', 'A m a m   ', 'B b m m   ', 'B m m b   '/
DATA (spcgr_name_set(64,i),i=1,6) /'C m c e   ', 'C c m e   ', 'A e m a   ', 'A e a m   ', 'B b e m   ', 'B m e b   '/
DATA (spcgr_name_set(65,i),i=1,6) /'C m m m   ', 'C m m m   ', 'A m m m   ', 'A m m m   ', 'B m m m   ', 'B m m m   '/
DATA (spcgr_name_set(66,i),i=1,6) /'C c c m   ', 'C c c m   ', 'A m a a   ', 'A m a a   ', 'B b m b   ', 'B b m b   '/
DATA (spcgr_name_set(67,i),i=1,6) /'C m m e   ', 'C m m b   ', 'A e m m   ', 'A c m m   ', 'B m e m   ', 'B m a m   '/
DATA (spcgr_name_set(68,i),i=1,6) /'C c c e   ', 'C c c b   ', 'A e a a   ', 'A c a a   ', 'B b e b   ', 'B b a b   '/
DATA (spcgr_name_set(69,i),i=1,6) /'F m m m   ', 'F m m m   ', 'F m m m   ', 'F m m m   ', 'F m m m   ', 'F m m m   '/
DATA (spcgr_name_set(70,i),i=1,6) /'F d d d   ', 'F d d d   ', 'F d d d   ', 'F d d d   ', 'F d d d   ', 'F d d d   '/
DATA (spcgr_name_set(71,i),i=1,6) /'I m m m   ', 'I m m m   ', 'I m m m   ', 'I m m m   ', 'I m m m   ', 'I m m m   '/
DATA (spcgr_name_set(72,i),i=1,6) /'I b a m   ', 'I b a m   ', 'I m c b   ', 'I m c b   ', 'I c m a   ', 'I c m a   '/
DATA (spcgr_name_set(73,i),i=1,6) /'I b c a   ', 'I c a b   ', 'I b c a   ', 'I c a b   ', 'I b c a   ', 'I c a b   '/
DATA (spcgr_name_set(74,i),i=1,6) /'I m m a   ', 'I m m b   ', 'I b m m   ', 'I c m m   ', 'I m c m   ', 'I m a m   '/
!
data spcgr_num_sym /                          &
  1,  2,  2,  2,  4,  2,  2,  4,  4,  4,   4,  8,  4,  4,  8,  4,  4,  4,  4,  8, & !   1   20
  8, 16,  8,  8,  4,  4,  4,  4,  4,  4,   4,  4,  4,  4,  8,  8,  8,  8,  8,  8, & !  21   40
  8, 16, 16,  8,  8,  8,  8,  8,  8,  8,   8,  8,  8,  8,  8,  8,  8,  8,  8,  8, & !  41   60
  8,  8, 16, 16, 16, 16, 16, 16, 32, 32,  16, 16, 16, 16,  4,  4,  4,  4,  8,  8, & !  61   80
  4,  8,  8,  8,  8,  8, 16, 16,  8,  8,   8,  8,  8,  8,  8,  8, 16, 16,  8,  8, & !  81  100
  8,  8,  8,  8,  8,  8, 16, 16, 16, 16,   8,  8,  8,  8,  8,  8,  8,  8, 16, 16, & ! 101  120
 16, 16, 16, 16, 16, 16, 16, 16, 16, 16,  16, 16, 16, 16, 16, 16, 16, 16, 32, 32, & ! 121  140
 32, 32,  3,  3,  3,  9,  6, 18,  6,  6,   6,  6,  6,  6, 18,  6,  6,  6,  6, 18, & ! 141  160
 18, 12, 12, 12, 12, 36, 36,  6,  6,  6,   6,  6,  6,  6, 12, 12, 12, 12, 12, 12, & ! 161  180
 12, 12, 12, 12, 12, 12, 12, 12, 12, 12,  24, 24, 24, 24, 12, 48, 24, 12, 24, 24, & ! 181  200
 24, 96, 96, 48, 24, 48, 24, 24, 96, 96,  48, 24, 24, 48, 24, 96, 48, 24, 96, 48, & ! 201  220
 48, 48, 48, 48,192,192,192,192, 96, 96,   2,  2,  4,  4,  4,  4,  4,  2,  2,  2, & ! 221  240
  2,  2,  2,  4,  4,  4,  4,  4,  4,  4,   4,  4,  4,  4,  4,  8,  8,  8,  8,  8, & ! 241  260
  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,   8,  8,  8,  8,  8,  8,  8,  8, 16, 32, & ! 261  280
  8,  8, 16, 16, 16, 16, 16, 16, 16, 16,  16, 32, 32,  3,  6,  6,  6,  6, 12, 12, & ! 281  300
 24, 96, 48, 48,192,192,  8,  8, 16,  8,   8, 16, 16, 16, 16                      & ! 301  315
     /
!
data spcgr_origin /                          &
  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,   1,  1,  1,  1,  1,  1,  1,  1,  1,  1, & !   1   20
  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,   1,  1,  1,  1,  1,  1,  1,  1,  1,  1, & !  21   40
  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,   1,  1,  1,  1,  1,  1,  1,  1,  1,  1, & !  41   60
  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,   1,  1,  1,  1,  1,  1,  1,  1,  1,  1, & !  61   80
  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,   1,  1,  1,  1,  1,  1,  1,  1,  1,  1, & !  81  100
  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,   1,  1,  1,  1,  1,  1,  1,  1,  1,  1, & ! 101  120
  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,   1,  1,  1,  1,  1,  1,  1,  1,  1,  1, & ! 121  140
  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,   1,  1,  1,  1,  1,  1,  1,  1,  1,  1, & ! 141  160
  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,   1,  1,  1,  1,  1,  1,  1,  1,  1,  1, & ! 161  180
  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,   1,  1,  1,  1,  1,  1,  1,  1,  1,  1, & ! 181  200
  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,   1,  1,  1,  1,  1,  1,  1,  1,  1,  1, & ! 201  220
  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,   1,  1,  1,  1,  1,  1,  1,  1,  1,  1, & ! 221  240
  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,   1,  1,  1,  1,  1,  1,  1,  1,  1,  1, & ! 241  260
  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,   1,  1,  1,  1,  1,  2,  2,  2,  2,  2, & ! 261  280
  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,   2,  2,  2,  2,  2,  2,  2,  2,  2,  2, & ! 281  300
  2,  2,  2,  2,  2,  2,  1,  1,  1,  1,   1,  1,  1,  1,  2                      & ! 301  315
     /
!
END MODULE spcgr_mod
