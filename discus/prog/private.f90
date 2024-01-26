module private_mod
!-
!  A module in which private calculations can be performed
!+
contains
!
!*******************************************************************************
!
subroutine do_private(line)
!-
!  Build lines of equal atom types
!+
use crystal_mod
use atom_env_mod
use chem_mod
use do_find_mod
!
use errlist_mod
use get_params_mod
use precision_mod
use param_mod
use lib_random_func
use take_param_mod
!
implicit none
!
character(len=*), intent(inout) :: line
!
integer, parameter :: MIN_PARA = 2
integer            :: MAXW
!
integer                             :: iat
integer                             :: ity
integer                             :: length
integer                             :: slow   ! Loop index
integer                             :: upper  ! Loop boundary
real(kind=PREC_DP), dimension(3, 6) :: direcs
real(kind=PREC_DP), dimension(3, 6) :: second
real(kind=PREC_DP), dimension(3)    :: vector
real(kind=PREC_DP), dimension(3)    :: other
real(kind=PREC_DP), dimension(3)    :: pos
real(kind=PREC_DP) :: r1
real(kind=PREC_DP) :: rmin
real(kind=PREC_DP) :: rmax
real(kind=PREC_DP) :: P_len
real(kind=PREC_DP) :: P_slen
logical                           :: fq
logical           , dimension(3)  :: fp
integer :: i, j, k
integer :: ioff
integer :: lp
integer :: ianz
!
character(LEN=PREC_STRING), DIMENSION(MAX(MIN_PARA,MAXSCAT+1)) :: cpara
real(KIND=PREC_DP) , DIMENSION(MAX(MIN_PARA,MAXSCAT+1)) :: werte
INTEGER            , DIMENSION(MAX(MIN_PARA,MAXSCAT+1)) :: lpara
logical, dimension(:), allocatable :: lchoose
!
INTEGER, PARAMETER :: NOPTIONAL = 2
INTEGER, PARAMETER :: O_PLEN    = 1
INTEGER, PARAMETER :: O_PSIG    = 2
CHARACTER(LEN=   6)       , DIMENSION(NOPTIONAL) :: oname    !Optional parameter names
CHARACTER(LEN=PREC_STRING), DIMENSION(NOPTIONAL) :: opara    !Optional parameter strings returned
INTEGER                   , DIMENSION(NOPTIONAL) :: loname   !Lenght opt. para name
INTEGER                   , DIMENSION(NOPTIONAL) :: lopara   !Lenght opt. para name returned
LOGICAL                   , DIMENSION(NOPTIONAL) :: lpresent !opt. para is present
REAL(KIND=PREC_DP)        , DIMENSION(NOPTIONAL) :: owerte   ! Calculated values
INTEGER, PARAMETER                               :: ncalc = 2 ! Number of values to calculate 
!
DATA oname  / 'plen', 'pslen'/
DATA loname /  4     ,  5      /
!
opara  =  (/ '0.0000', '0.0000' /)   ! Always provide fresh default values
lopara =  (/  6,        6       /)
owerte =  (/  0.0,      0.0     /)
!
MAXW = MAX(MIN_PARA,MAXSCAT+1)
lp = len_trim(line)
call get_params(line, ianz, cpara, lpara, maxw, lp)
if(ier_num/=0) return
call get_optional(ianz, MAXW, cpara, lpara, NOPTIONAL,  ncalc, &
                  oname, loname, opara, lopara, lpresent, owerte)
if(ier_num/=0) return
!
direcs(1, 1) =  0.0
direcs(2, 1) =  0.5
direcs(3, 1) =  0.5
!
direcs(1, 2) =  0.0
direcs(2, 2) = -0.5
direcs(3, 2) =  0.5
!
direcs(1, 3) =  0.5
direcs(2, 3) =  0.0
direcs(3, 3) =  0.5
!
direcs(1, 4) = -0.5
direcs(2, 4) =  0.0
direcs(3, 4) =  0.5
!
direcs(1, 5) =  0.5
direcs(2, 5) =  0.5
direcs(3, 5) =  0.0
!
direcs(1, 6) = -0.5
direcs(2, 6) =  0.5
direcs(3, 6) =  0.0
!
second(1, 1) =  0.5
second(2, 1) =  0.5
second(3, 1) =  0.0
!
second(1, 2) =  0.5
second(2, 2) =  0.0
second(3, 2) =  0.5
!
second(1, 3) =  0.0
second(2, 3) =  0.5
second(3, 3) =  0.5
!
second(1, 4) =  0.5
second(2, 4) =  0.5
second(3, 4) =  0.0
!
second(1, 5) =  0.5
second(2, 5) =  0.0
second(3, 5) =  0.5
!
second(1, 6) =  0.0
second(2, 6) =  0.5
second(3, 6) =  0.5
!
fp(1) = chem_period(1)
fp(2) = chem_period(2)
fp(3) = chem_period(3)
fq = chem_quick
!
P_len  = owerte(O_PLEN) !cr_icc(1)*2 ! * 2./3.
P_slen = owerte(O_PSIG ) ! 3.0
!
upper = 20*cr_natoms
!
allocate(lchoose(1:cr_natoms))
!
do k = 1,2
   lchoose = .true.
loop_main: do slow=1, cr_natoms
!loop_main: do slow=1, upper
!   r1  = ran1(0)
!   iat = int(cr_natoms*r1) + 1
   r1   = ran1(0)
   ioff = int(cr_natoms*r1)
   iat  = 1
   loop_find: do j=1, cr_natoms
      iat = mod(slow + ioff + j, cr_natoms) + 1
      if(lchoose(j)) then
         exit loop_find
      endif
   end do loop_find
   r1  = ran1(0)
   ity = int(cr_nscat*r1) + 1
   length = nint(P_len + gasdev(P_slen))
!
   vector = direcs(:, ity)
   other  = second(:, ity)
   loop_inner: do i=1, length
      cr_iscat(1,iat) = ity
      pos = cr_pos(:,iat) + vector
      do j=1,3
         if(pos(j) > cr_dim(j,2)) pos(j) = pos(j) - cr_icc(j)
      enddo
      ianz     = 1             ! Just one atom type
      werte(1) = -1            ! Find all atom types
      rmin = 0.0
      rmax = 0.5
      call do_find_env(ianz, werte, MAXSCAT, pos, rmin, rmax, fq, fp)
      iat = atom_env(1)
      pos = cr_pos(:,iat) + other
      do j=1,3
         if(pos(j) > cr_dim(j,2)) pos(j) = pos(j) - cr_icc(j)
      enddo
      ianz     = 1             ! Just one atom type
      werte(1) = -1            ! Find all atom types
      rmin = 0.0
      rmax = 0.5
      call do_find_env(ianz, werte, MAXSCAT, pos, rmin, rmax, fq, fp)
      cr_iscat(1,atom_env(1)) = ity
  end do loop_inner
end do loop_main
end do 
!
deallocate(lchoose)
!
end subroutine do_private
!
!*******************************************************************************
!
end module private_mod
