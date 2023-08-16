module lib_use_coor_mod
!
!  Determine indices 1,2,3 for Abs, Ord, Top
!
private
public lib_get_use_coor
!
contains
!
!*******************************************************************************
!
subroutine lib_get_use_coor(vectors, calc_coor, use_coor)
!-
!  The routine determines which indices vary along the  Abscissa, Ordinate
!  and top axis
!+
!
use precision_mod
!
implicit none
!
real(kind=PREC_DP), dimension(3,3), intent(in)  :: vectors   ! Vectors abs, ord, top
logical                           , intent(out) :: calc_coor  ! Need to calculate coordinates
integer           , dimension(3)  , intent(out) :: use_coor   ! Indices for A, O, T
!
real(kind=PREC_DP), dimension(3) :: cosa  ! Rough cosine to [100], [010], [001]
!
!  Get direction cosine to [100], [010], [001]
cosa(1) = vectors(1,1) / sqrt(vectors(1,1)**2 + vectors(2,1)**2 + vectors(3,1)**2)
cosa(2) = vectors(2,2) / sqrt(vectors(1,2)**2 + vectors(2,2)**2 + vectors(3,2)**2)
cosa(3) = vectors(3,3) / sqrt(vectors(1,3)**2 + vectors(2,3)**2 + vectors(3,3)**2)
!
use_coor(1) = 1      ! Default for regular coordinates
use_coor(2) = 2
use_coor(3) = 3
!
cond_need: if(abs(abs(cosa(1))-1.00D0)>0.05D0 .or. &
              abs(abs(cosa(2))-1.00d0)>0.05D0 .or. &
              abs(abs(cosa(3))-1.00D0)>0.05D0     ) then    ! Need to signal coordinate calculation
  calc_coor = .true.
  cosa(1) = vectors(1,1) / sqrt(vectors(1,1)**2 + vectors(2,1)**2 + vectors(3,1)**2)
  cosa(2) = vectors(2,1) / sqrt(vectors(1,1)**2 + vectors(2,1)**2 + vectors(3,1)**2)
  cosa(3) = vectors(3,1) / sqrt(vectors(1,1)**2 + vectors(2,1)**2 + vectors(3,1)**2)
!
  use_coor(1) = maxloc(abs(cosa),1)     ! This is axis closest to [1, 0, 0]
  cosa(1) = vectors(1,2) / sqrt(vectors(1,2)**2 + vectors(2,2)**2 + vectors(3,2)**2)
  cosa(2) = vectors(2,2) / sqrt(vectors(1,2)**2 + vectors(2,2)**2 + vectors(3,2)**2)
  cosa(3) = vectors(3,2) / sqrt(vectors(1,2)**2 + vectors(2,2)**2 + vectors(3,2)**2)
!
  if(maxloc(abs(cosa),1) /= use_coor(1)) then
     use_coor(2) = maxloc(abs(cosa),1)   ! We can use axis closest to [0, 1, 0]
  else
     cosa(use_coor(1)) =  0.0D0   ! Invalidate the "x" axis
     use_coor(2) = maxloc(abs(cosa),1)   ! We have to use another axis
  endif
  cosa(1) = vectors(1,3) / sqrt(vectors(1,3)**2 + vectors(2,3)**2 + vectors(3,3)**2)
  cosa(2) = vectors(2,3) / sqrt(vectors(1,3)**2 + vectors(2,3)**2 + vectors(3,3)**2)
  cosa(3) = vectors(3,3) / sqrt(vectors(1,3)**2 + vectors(2,3)**2 + vectors(3,3)**2)
!
  if(maxloc(abs(cosa),1) /= use_coor(1) .and. maxloc(abs(cosa),1) /= use_coor(2)) then
     use_coor(3) = maxloc(abs(cosa),1)   ! We can use axis closest to [0, 0, 1]
  else
     if(use_coor(1) == 1) then    ! X-index is taken
        if(use_coor(2) == 2) then   ! Y-index is taken
           use_coor(3) = 3
        else
           use_coor(3) = 2
        endif
     elseif(use_coor(1) == 2) then    ! Y-index is taken
        if(use_coor(2) == 1) then   ! X-index is taken
           use_coor(3) = 3
        else
           use_coor(3) = 1
        endif
     elseif(use_coor(1) == 3) then    ! Z-index is taken
        if(use_coor(2) == 1) then   ! X-index is taken
           use_coor(3) = 2
        else
           use_coor(3) = 1
        endif
     endif
  endif
else cond_need                    ! Reguar sequence 1, 2, 3 
  calc_coor = .false.
  use_coor(1) = 1
  use_coor(2) = 2
  use_coor(3) = 3
endif cond_need
!
end subroutine lib_get_use_coor
!
!*******************************************************************************
!
end module lib_use_coor_mod
