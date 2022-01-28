module discus_3dpdf_mod
!-
!  Routines to interpret a 3D-PDF
!+
!
use precision_mod
!
implicit none
!
character(LEN=PREC_STRING)                        :: three_data         ! Input file name
real(KIND=PREC_DP), dimension(:,:,:), allocatable :: three_map          ! Original map
real(KIND=PREC_DP), dimension(:,:,:), allocatable :: three_map_temp     ! Copy of map
integer           , dimension(3)                  :: three_dims         ! Dimensions in pixel
real(KIND=PREC_DP), dimension(3)                  :: three_llims        ! Lower left bottom corner
real(KIND=PREC_DP), dimension(3,3)                :: three_steps        ! Steps in direct space
!                                                                       ! 1.st dim: [hkl], 2nd:[abs, ord, top]
!
integer           , dimension(3)                  :: three_width        ! Window in which value must be extreme
real(KIND=PREC_DP), dimension(:,:,:), allocatable :: three_peaks        ! (x:z,1:npeaks, neg:pos)
integer                                           :: three_max_peaks    ! Maximum number of peaks
integer           , dimension(2)                  :: three_npeaks       ! Number of peaks (neg:pos)
!
integer                                           :: three_npatt = 0    ! Number of crystal Patterson vectors
real(kind=PREC_DP), dimension(:,:), allocatable   :: three_patt         ! Crystal Patterson vectors
integer           , dimension(:,:), allocatable   :: three_patt_isite   ! Which sites start/end of patterson vector
logical           , dimension(:,:), allocatable   :: three_matched_xtal ! Matched crystal Patterson vector
logical           , dimension(:,:), allocatable   :: three_matched_3dpdf! Matched experimental Patterson vector
logical           , dimension(:,:), allocatable   :: three_matched_site ! This site is matched by 3d-PDF peak (-,+)
!
integer :: three_node_num
!
type :: three_match
   integer, dimension(2)      :: isite                          ! sites at start/end of Patterson vector
   real(kind=PREC_DP)         :: angle                          ! Angle to crystal vector
   real(kind=PREC_DP)         :: dist                           ! Deviation in distance
   type(three_match), pointer :: next
end type three_match
!
type :: three_list
   type(three_match), pointer :: ptr
end type three_list
!
type(three_list), dimension(:,:), allocatable :: three_res      ! Results for experimental vectors
!
!type(three_match) :: head
!type(three_match) :: temp
!
private
public  three_main
!
contains
!
!*******************************************************************************
!
subroutine three_main !(infile)
!
!-
!  Main menu for 3D-PDF interpretations
!
use errlist_mod
!
use build_name_mod
use calc_expr_mod
use class_macro_internal
use doact_mod
use get_params_mod
use lib_errlist_func
use lib_help
use lib_macro_func
use prompt_mod
use str_comp_mod
use sup_mod
!
implicit none
!
integer, parameter :: MIN_PARA = 20
integer            :: maxw               
!
character(len=PREC_STRING), dimension(MIN_PARA) :: cpara   ! Array for parameter strings
integer                   , dimension(MIN_PARA) :: lpara   ! Array for parameter lengths
real(kind=PREC_DP)        , dimension(MIN_PARA) :: werte   ! Array for parameter lengths
!
character(len=5)           :: befehl           ! The main command verb
character(len=len(prompt)) :: orig_prompt      ! Prompt in previous menu
character(len=PREC_STRING) :: line             ! Dummy character string
character(len=PREC_STRING) :: zeile            ! Dummy character string
!
integer                    :: length           ! length of a character string
integer                    :: lbef             ! length of main command verb 'bef'
integer                    :: lp               ! length of a character string
integer                    :: indxg            ! Location of '=' if anystring
integer                    :: ianz             ! Number of parameters
logical                    :: lend             ! end of operations, go back to upper menu
logical, save              :: linit = .true.
!
maxw = MIN_PARA
!
three_max_peaks = 50
three_npeaks = 0
!
orig_prompt = prompt
prompt = prompt (1:len_trim(prompt) ) //'/3dpdf'
lend   = .false.
!
loop_main: do while(.not.lend)
!
   if_error: IF (ier_num.ne.0) then
      call errlist
      IF (ier_sta.ne.ER_S_LIVE) then
         IF (lmakro .OR. lmakro_error) then  ! Error within macro or termination errror
            IF(sprompt /= prompt ) then
               ier_num = -10
               ier_typ = ER_COMM
               ier_msg(1) = ' Error occured in stack menu'
               prompt_status = PROMPT_ON
               prompt = orig_prompt
               RETURN
            ELSE
               IF(lmacro_close) then
                  CALL macro_close
                  prompt_status = PROMPT_ON
               ENDIF
            ENDIF
         ENDIF
         IF (lblock) then
            ier_num = -11
            ier_typ = ER_COMM
            prompt_status = PROMPT_ON
            prompt = orig_prompt
            RETURN
         ENDIF
         call no_error
         lmakro_error = .FALSE.
         sprompt = ' '
      ENDIF
   ENDIF if_error
!
   call get_cmd(line, length, befehl, lbef, zeile, lp, prompt)
   if(ier_num /= 0) cycle loop_main
!
   if(line == ' '      .or. line(1:1) == '#' .or. &
      line == char(13) .or. line(1:1) == '!'        ) cycle loop_main
!                                                                       
!     ----search for "="                                                
!                                                                       
   indxg = index (line, '=')
   IF (indxg.ne.0                                              &
        .and..not. (str_comp (befehl, 'echo', 2, lbef, 4) )    &
        .and..not. (str_comp (befehl, 'syst', 2, lbef, 4) )    &
        .and..not. (str_comp (befehl, 'help', 2, lbef, 4) .or. &
                    str_comp (befehl, '?   ', 2, lbef, 4) )    &
        .and. index(line,'==') == 0                            ) then
!                                                                       
!     ------evaluate an expression and assign the value to a variable   
!                                                                       
            call do_math(line, indxg, length)
      cycle loop_main
   endif
!
   if(linit) then           ! Initialization do basic allocation
      
      linit = .false.       ! Initialization is done until a 'reset'
   endif
!
!---  exit 'exit'
!
   if(str_comp(befehl, 'exit', 3, lbef, 4)) then
      lend = .true.
!
!     ----help 'help','?'
!
   elseif(str_comp(befehl, 'help', 2, lbef, 4) .or.            &
          str_comp(befehl, '?   ', 1, lbef, 4)        ) then
      if(str_comp(zeile, 'errors', 2, lp, 6) ) then
         lp = lp + 7
         call do_hel('discus '//zeile, lp)
      else
         lp = lp + 12
         call do_hel('discus three '//zeile, lp)
      endif
!
!--- Reset everything 'rese'
!
   elseif(str_comp(befehl, 'rese', 3, lbef, 4)) then
      call three_reset
      linit = .true.
!
!---  Define data file 'data'
!
   elseif(str_comp(befehl, 'data', 3, lbef, 4)) then
      call get_params(zeile, ianz, cpara, lpara, maxw, lp)
      if(ier_num/=0) cycle loop_main
      call do_build_name(ianz, cpara, lpara, werte, maxw, 1)
      if(ier_num/=0) cycle loop_main
!
      three_data = cpara(1)
!
!---  run the interpretation 'run'
!
   elseif(str_comp(befehl, 'run', 3, lbef, 3)) then
      call three_run
   endif
!
enddo loop_main
!
end subroutine three_main
!
!*******************************************************************************
!
subroutine three_run
!-
!  Execute the Patterson peak search
!+
!
if(allocated(three_peaks)) deallocate(three_peaks)
allocate(three_peaks(5,three_max_peaks,2))    ! (x_z_value_length,1:npeaks, neg:pos)
three_peaks = 0.0D0
!
call three_load(three_data)
!call three_glatt
call three_find
call three_show
call three_inter
!
end subroutine three_run
!
!*******************************************************************************
!
subroutine three_load(infile)
!-
!  Loads the input file into the H5 via KUPLOT and create a working copy
!
use kuplot_load_h5
!
use errlist_mod
use precision_mod
!
implicit none
!
character(len=PREC_STRING), intent(in) :: infile
character(len=PREC_STRING) :: string
integer                    :: length
integer                    :: node_number
integer, dimension(3)      :: h5_dims
!integer, dimension(3)      ::    dims
logical                    :: lecho
!
!infile = 'PDF/example2.3DPDF'
!infile = 'PDF/M4_replace_cent_100_srp.3DPDF'
!infile = 'PDF/M4_replace_carb_100_srp.3DPDF'
!infile = 'PDF/M4_replace_o1_100_srp.3DPDF'
!infile = 'PDF/M3_rotate_carb_100_nrm.3DPDF'
!read(*,'(a)') infile
string = 'h5, ' // infile(1:len_trim(infile))
length = len_trim(string)
lecho  = .true.

call do_load (string, length, lecho)
node_number = 1
call hdf5_find_node(node_number, ier_num, ier_typ)
call hdf5_get_dims(node_number, h5_dims)
three_dims(1) = h5_dims(3)
three_dims(2) = h5_dims(2)
three_dims(3) = h5_dims(1)
call hdf5_get_llims(node_number, three_llims)
call hdf5_get_steps(node_number, three_steps)
!call hdf5_get_steps(node_number, three_steps)
!three_steps      = 0.0D0
!three_steps(1,1) = steps(1)
!three_steps(2,2) = steps(2)
!three_steps(3,3) = steps(3)
!
if(allocated(three_map)) deallocate(three_map)
!
allocate(three_map(three_dims(1),three_dims(2), three_dims(3)))
call hdf5_get_tmap(three_dims, three_map)
!
if(allocated(three_map_temp)) deallocate(three_map_temp)
allocate(three_map_temp(three_dims(1),three_dims(2), three_dims(3)))
three_map_temp = three_map
!
end subroutine three_load
!
!*******************************************************************************
!
subroutine three_stat
!-
!  Get intensity statistics on 3D-DATA set
!+
!
implicit none
!
end subroutine three_stat
!
!*******************************************************************************
!
subroutine three_glatt
!+
! Simple smoothing to be replaced by a Filter
!-
implicit none
!
integer, dimension(3,2) :: limits
integer :: i,j,k
integer :: ii,jj,kk
integer, dimension(3) :: width
logical :: lsmooth
real(kind=PREC_DP) :: summe
real(kind=PREC_DP) :: weight
!
!
width = 1
limits(1,1) = min(1+width(1), three_dims(1)        )
limits(1,2) = max(1         , three_dims(1)-width(1)-1)
limits(2,1) = min(1+width(2), three_dims(2)        )
limits(2,2) = max(1         , three_dims(2)-width(2)-1)
limits(3,1) = min(1+width(3), three_dims(3)        )
limits(3,2) = max(1         , three_dims(3)-width(3)-1)
if(limits(1,2)==1) width(1) = 0
if(limits(2,2)==1) width(2) = 0
if(limits(3,2)==1) width(3) = 0
weight = (2*width(1)+1)*(2*width(2)+1)* (2*width(3)+1)
!
!
lsmooth = .false.
if(lsmooth) then
do k=limits(3,1),limits(3,2)
   do j=limits(2,1),limits(2,2)
      do i=limits(1,1),limits(1,2)
         summe = 0.0D0
         do kk = k-width(3), k+width(3)
            do jj = j-width(2), j+width(2)
               do ii = i-width(1), i+width(1)
                  summe = summe + three_map(ii,jj,kk)
               enddo
            enddo
         enddo
         three_map_temp(i,j,k) = summe/weight
      enddo
   enddo
enddo
endif
!
end subroutine three_glatt
!
!*******************************************************************************
!
subroutine three_find
!-
! Find peaks
!+
!
use metric_mod
!
implicit none
!
real(kind=PREC_DP), dimension(3), parameter :: NULLV = (/ 0.0D0, 0.0D0, 0.0D0/)
!
integer, dimension(3,2) :: limits      ! Limits for coordinates [0,1] plus safety
integer, dimension(3)   :: ipos        ! Window in which value must be extreme
integer :: i,j,k
!integer :: ii,jj,kk
logical :: lpeak
real(kind=PREC_DP) :: aver
real(kind=PREC_DP) :: summe
real(kind=PREC_DP) :: sigma
real(kind=PREC_DP) :: weight
real(kind=PREC_DP) :: v_length
real(kind=PREC_DP), dimension(3) :: pos
real(kind=PREC_DP), dimension(:,:,:), allocatable ::  field
real(kind=PREC_DP), dimension(3) :: vect
!
call three_get_ipos(three_dims, three_steps, three_llims, ipos)
!
!ipos(1) = nint((0.0D0-three_llims(1))/three_steps(1,1)) + 1
!ipos(2) = nint((0.0D0-three_llims(2))/three_steps(2,2)) + 1
!ipos(3) = 1
!write(*,*) ' IPOS ', ipos, three_map_temp(ipos(1), ipos(2), ipos(3))
!
limits = 1
do i=1,3
   if(three_dims(i)>1) then
      limits(i,1) = max(1            ,nint((-0.9-three_llims(i))/three_steps(i,i)) + 1)
      limits(i,2) = min(three_dims(i),nint(( 0.9-three_llims(i))/three_steps(i,i)) + 1)
   endif
enddo
!
! Determine three_width from FWHM(0,0,0)
!
three_width = 0 !1
if(limits(1,2)==1) then 
   three_width(1) = 0
else
   loop_fw1:do i=ipos(1), limits(1,2)
   if(three_map_temp(i, ipos(2), ipos(3))<0.5*three_map_temp(ipos(1), ipos(2), ipos(3))) exit loop_fw1
      three_width(1) = three_width(1) + 1
   enddo loop_fw1
endif
if(limits(2,2)==1) then
   three_width(2) = 0
else
   loop_fw2:do j=ipos(2), limits(2,2)
   if(three_map_temp(ipos(2), j, ipos(3))<0.5*three_map_temp(ipos(1), ipos(2), ipos(3))) exit loop_fw2
      three_width(2) = three_width(2) + 1
   enddo loop_fw2
endif
if(limits(3,2)==1) then
   three_width(3) = 0
else
   loop_fw3:do k=ipos(3), limits(3,2)
   if(three_map_temp(ipos(2), ipos(3), k)<0.5*three_map_temp(ipos(1), ipos(2), ipos(3))) exit loop_fw3
      three_width(3) = three_width(3) + 1
   enddo loop_fw3
endif
weight = (2*three_width(1)+1)*(2*three_width(2)+1)* (2*three_width(3)+1)
!
allocate(field(-three_width(1):three_width(1),                                  &
               -three_width(2):three_width(2),                                  &
               -three_width(3):three_width(3)))
! 
write(*,*) 'ZERO ',ipos,three_map_temp(ipos(1), ipos(2), ipos(3))
!
write(*,*) 'Peak search from to ', limits(1,:), three_width(1), weight
write(*,*) 'Peak search from to ', limits(2,:), three_width(2)
write(*,*) 'Peak search from to ', limits(3,:), three_width(3)
!
do k=limits(3,1),limits(3,2)
   do j=limits(2,1),limits(2,2)
      loop_i: do i=limits(1,1),limits(1,2)
         summe = 0.0D0
         lpeak = .true.
         field = three_map_temp(i-three_width(1): i+three_width(1), &
                                j-three_width(2): j+three_width(2), &
                                k-three_width(3): k+three_width(3))
         aver  = sum(field)/weight
         sigma = sqrt(sum((field-aver )**2)/(weight-1.))
!
         if(abs(three_map_temp(i,j,k)) > abs(aver)          .and.  &
            abs(three_map_temp(i,j,k)/three_map_temp(ipos(1), ipos(2), ipos(3)))>0.05 .and. &
            (  (three_map_temp(i,j,k))==maxval(field) .or.         &
               (three_map_temp(i,j,k))==minval(field)     ) .and.  &
            abs(three_map_temp(i,j,k)-aver)/sigma> 1.0D0) then
            pos(1) = three_llims(1)+(i-1)*three_steps(1,1)
            pos(2) = three_llims(2)+(j-1)*three_steps(2,2)
            pos(3) = three_llims(3)+(k-1)*three_steps(3,3)
!write(*,'(a,3i4,3f6.2,4g12.3e3)') ' Peak at ', i,j,k, pos, three_map_temp(i,j,k), &
!aver, sigma,abs(three_map_temp(i,j,k)-aver)/sigma
            vect = pos
            v_length = do_blen(.true., vect, NULLV)
            call three_insert(three_max_peaks, three_npeaks,                    &
                 three_map_temp(i, j, k), pos, v_length, three_peaks)
         endif
      enddo loop_i
   enddo
enddo
end subroutine three_find
!
!*******************************************************************************
!
subroutine three_get_ipos(dims, steps, llims, ipos)
!-
!  Determine indices of (0,0,0)
!+
use matrix_mod
use precision_mod
!
implicit none
!
integer           , dimension(3), intent(in)  :: dims
real(kind=PREC_DP), dimension(3,3), intent(in)  :: steps
real(kind=PREC_DP), dimension(3), intent(in)  :: llims
integer           , dimension(3), intent(out) :: ipos
!
integer :: ndims
real(kind=PREC_DP), dimension(3,3) :: imat_vi
real(kind=PREC_DP), dimension(2  ) :: llims2
real(kind=PREC_DP), dimension(2,2) :: mat2
real(kind=PREC_DP), dimension(2,2) :: imat2
!
ndims = 0
if(dims(1)>1) ndims = ndims +1
if(dims(2)>1) ndims = ndims +1
if(dims(3)>1) ndims = ndims +1
!
if(ndims==3) then                                ! Full 3D
   call matinv3(steps, imat_vi)
   ipos = -nint(matmul(imat_vi, llims)) + 1
elseif(ndims==2) then                            ! 2D data set
   if(dims(3) == 1) then                         ! x-y field
      mat2 = steps(1:2, 1:2)
      call matinv2(mat2, imat2)
      ipos(1:2) = -nint(matmul(imat2, llims(1:2))) + 1
      ipos(3) = 1
   elseif(dims(2) == 1) then                     ! x-z field
      mat2 = steps(1:3:2, 1:3:2)
      call matinv2(mat2, imat2)
      ipos(1:3:2) = -nint(matmul(imat2, llims(1:3:2))) + 1
      ipos(2) = 1
   elseif(dims(1) == 1) then                     ! y-z field
      llims2 = llims(2:3)
      mat2 = steps(2:3, 2:3)
      call matinv2(mat2, imat2)
      ipos(2:3) = -nint(matmul(imat2, llims(2:3))) + 1
      ipos(1) = 1
   endif
endif
!
!
end subroutine three_get_ipos
!
!*******************************************************************************
!
subroutine three_show
!-
! List the obtained peaks, general show
!+
!use metric_mod
use prompt_mod
!
implicit none
!
!real(kind=PREC_SP), dimension(3), parameter :: NULLV = (/ 0.0D0, 0.0D0, 0.0D0/)
integer :: i      ! Dummy loop index
!real(kind=PREC_SP), dimension(3) :: vect
!
write(output_io,'(a)')
write(output_io,'(a)') ' Patterson Search results '
write(output_io,'(a)')
write(output_io,'(27x,a)') 'Pix     u          v          w'
write(output_io,'(a,i6,3(2x,f9.4))') ' Search width abscissa :', three_width(1), three_steps(:,1)
write(output_io,'(a,i6,3(2x,f9.4))') ' Search width ordinate :', three_width(2), three_steps(:,2)
write(output_io,'(a,i6,3(2x,f9.4))') ' Search width top      :', three_width(3), three_steps(:,3)
write(output_io,'(a)')

write(output_io,'(a,i5,5x,a)')  ' Number of peaks > 0   : ', three_npeaks(2), &
              'u          v          w           Length    Height      rel'
loop_pos: do i=1, three_npeaks(2)
!  vect=three_peaks(1:3,i,2)
      write(output_io,'(a,i6, 3(2x,f9.4),3x,F9.4,2x,g12.3e3,2x,i3,a)') ' Peak at               :', i,     &
   three_peaks(1:5,i,2), &
                      nint(100.*abs(three_peaks(5,i,2)/three_peaks(5,1,2))),'%'
!  do_blen(.true.,vect, NULLV), &
!  three_peaks(4,i,2),nint(100.*abs(three_peaks(4,i,2)/three_peaks(4,1,2))),'%'
enddo loop_pos
!
write(output_io,'(a)')
write(output_io,'(a,i5,5x,a)')  ' Number of peaks < 0   : ', three_npeaks(1), &
              'u          v          w           Length    Height      rel'
loop_neg: do i=1, three_npeaks(1)
!  vect=three_peaks(1:3,i,1)
   write(output_io,'(a,i6, 3(2x,f9.4),3x,F9.4,2x,g12.3e3,2x,i3,a)') ' Peak at               :', i,     &
      three_peaks(1:5,i,1), &
                      nint(100.*abs(three_peaks(5,i,1)/three_peaks(5,1,2))),'%'
!  do_blen(.true.,vect, NULLV), &
!  three_peaks(4,i,1),nint(100.*abs(three_peaks(4,i,1)/three_peaks(4,1,2))),'%'
enddo loop_neg
!
write(output_io,'(a)')
!
end subroutine three_show
!
!*******************************************************************************
!
subroutine three_insert(NMAX, npeaks, dpdf, pos, v_length, peaks)
!-
!  Sorts the peaks by height into the array
!+
!
implicit none
!
integer                                  , intent(in)    :: NMAX      ! Maximum peak number
integer           , dimension(2)         , intent(inout) :: npeaks    ! Actual peak number
real(kind=PREC_DP)                       , intent(in)    :: dpdf      ! 3D-PDF value
real(kind=PREC_DP), dimension(3)         , intent(in)    :: pos       ! Current position
real(kind=PREC_DP)                       , intent(in)    :: v_length  ! 3D-PDF value
real(kind=PREC_DP), dimension(5, NMAX, 2), intent(inout) :: peaks     ! Stored positions
integer            :: i,j,k     ! Dummy loop indices
integer            :: nn        ! Dummy loop indices
integer            :: isig      ! 1 for negative ; 2 for positive
real(kind=PREC_DP) :: rsig      ! -1 for neg, 1 for pos
!
if(dpdf<0.0D0) then             ! Negative 3D-PDF value
   isig =  1
   rsig = -1.0D0
else                            ! Positive 3D-PDF value
   isig =  2
   rsig =  1.0D0
endif
!
nn = npeaks(isig)
k  = 0
loop_outer:do i=1, nn
   if(rsig*peaks(5,i,isig) < rsig*dpdf) then        ! This peak has smaller absolute value
      loop_inner: do j=nn, i, -1
         peaks(:, j+1, isig) = peaks(:, j, isig)
      enddo loop_inner
      exit loop_outer
   endif
   k = k + 1
enddo loop_outer
!
peaks(1:3,k+1, isig) = pos
peaks( 4 ,k+1, isig) = v_length
peaks( 5 ,k+1, isig) = dpdf
npeaks(isig) = npeaks(isig) + 1
!
end subroutine three_insert
!
!*******************************************************************************
!
subroutine three_inter
!-
!  Attempt an (atomatic/guided) interpretation of the 3DPDF
!+
!
use crystal_mod
use chem_mod
use chem_aver_mod
use metric_mod
!
character(len=48), dimension(2), parameter :: c_types = (/ &
         ' --------------- Negative peaks ----------------',  &
         ' +++++++++++++++ Positive peaks ++++++++++++++++' /)
!
logical, parameter :: lout  = .true.         ! Display average structure yes/no
logical, parameter :: lsite = .true.         ! Treat all atoms on one site as one
real(kind=PREC_DP), dimension(3), parameter :: NULLV = (/ 0.0D0, 0.0D0, 0.0D0/)
!
character(len=PREC_STRING)                  :: line       ! Dummy line
real(kind=PREC_DP), dimension(3)            :: vect       ! Dummy vector
real(kind=PREC_DP)                          :: sigma_length  ! Uncertainty in match Patterson vector length
integer :: i,j,k   ! Dummy loop indices
integer :: np      ! Dummy loop index
!
type(three_match), pointer :: head
type(three_match), pointer :: temp
type(three_match), pointer :: tail
!
sigma_length = 0.0D0                         ! Determine uncertainty in Patterson lengths
do i=1,3
   vect = three_steps(:,i)
   sigma_length = max(sigma_length, do_blen(.true., vect, NULLV))
enddo
sigma_length = 0.75*sigma_length
!
call three_clear_nodes                       ! Celan up previous nodes
nullify(head)
nullify(temp)
nullify(tail)
!
call chem_aver (lout, lsite)                 ! Get the average structure
!
three_npatt = cr_ncatoms*(cr_ncatoms-1)/2    ! Number of Patterson vectors 
write(*,*)  three_npatt, cr_ncatoms
if(allocated(three_patt)) deallocate(three_patt)
if(allocated(three_patt_isite)) deallocate(three_patt_isite)
if(allocated(three_matched_xtal )) deallocate(three_matched_xtal)
if(allocated(three_matched_3dpdf)) deallocate(three_matched_3dpdf)
if(allocated(three_matched_site )) deallocate(three_matched_site )
allocate(three_patt(4, three_npatt))
allocate(three_patt_isite(2, three_npatt))
allocate(three_matched_xtal (three_npatt,2))
allocate(three_matched_3dpdf(maxval(three_npeaks(:),1),2))     ! WORK 
allocate(three_matched_site(1:cr_ncatoms,2))
three_patt       = 0.0D0
three_patt_isite = 0
three_matched_xtal  = .false.                      ! No vectors matched yet
three_matched_3dpdf = .false.                      ! No vectors matched yet
three_matched_site  = .false.                      ! No vectors matched yet
!
np = 0
do i=1, cr_ncatoms-1
   do j=i+1, cr_ncatoms
      np = np + 1
      vect = chem_ave_pos(:,j) - chem_ave_pos(:,i)
      three_patt(1:3, np) = chem_ave_pos(:,j) - chem_ave_pos(:,i)
      three_patt(4  , np) = do_blen(.true., vect, NULLV)
      three_patt_isite(1,np) = i
      three_patt_isite(2,np) = j
   enddo
enddo
!
if(allocated(three_res)) deallocate(three_res)
allocate(three_res(1:maxval(three_npeaks,1),2))
do i=1, maxval(three_npeaks,1)
  nullify(three_res(i,1)%ptr)
  nullify(three_res(i,2)%ptr)
enddo
write(*,*) ' THREE_RES ', ubound(three_res)
!
!  Attempt immediate match of an experimental vector to any Patterson vector
!  in the unit cell
!
call three_inter_match(sigma_length)
!
!
! Test is a site is matched by any 3D-PDF Peak
!
loop_isite_match:do i=1, cr_ncatoms
   do j=1, three_npatt
     if((three_patt_isite(1,j)==i .or. three_patt_isite(2,j)==i) .and.          &
        (three_matched_xtal(j,1)                               )                & 
       ) then
        three_matched_site(i,1) = .true.
     endif
     if((three_patt_isite(1,j)==i .or. three_patt_isite(2,j)==i) .and.          &
        (                              three_matched_xtal(j,2) )                & 
       ) then
        three_matched_site(i,2) = .true.
     endif
        if(all(three_matched_site(i,:))) cycle loop_isite_match
   enddo
enddo loop_isite_match
!
write(*,*) ' SIGMA ', sigma_length
write(*,*)
write(*,'(a)') ' List of Patterson vectors in structure'
write(*,*)
do i=1, three_npatt
write(*,'(a, 3(2x,f9.4),3x,f9.4,2x, l1,2x,l1,2(2x,i3))') ' Patt ', three_patt(:,i), three_matched_xtal (i,:), &
   three_patt_isite(:,i)
enddo
line = ' '
j = 0
do i=1, cr_ncatoms
   if(.not. three_matched_site(i,2)) then
      write(line(5*j+1:5*j+4),'(i4)') i
      j = j + 1
   endif
enddo
write(*,*)
write(*,'(a,a)') ' Positive Sites not matched by experimental vectors ', line(1:len_trim(line))
line = ' '
j = 0
do i=1, cr_ncatoms
   if(.not. three_matched_site(i,1)) then
      write(line(5*j+1:5*j+4),'(i4)') i
      j = j + 1
   endif
enddo
write(*,*)
write(*,'(a,a)') ' Negative Sites not matched by experimental vectors ', line(1:len_trim(line))
!
write(*,*) 
write(*,*)  ' List of experimental 3D-PDF vectors '
do k=2, 1, -1
write(*,*)
write(*,'(a)') c_types(k)
write(*,*)
do i=1, three_npeaks(k) ! List recognized experimental peaks
   if(three_peaks(4,i,k)> 0.0001D0) then
      write(*,'(a, i3, 3(2x,f9.4),4x,f9.4,2x,g12.3e3)') ' 3D-PDF vector      ',i, three_peaks(1:5,i,k)
      temp => three_res(i,k)%ptr
      loop_nodes: do
         if(associated(temp)) then
            write(*,'(a,4x,2i4)')  ' Sites in unit cell ', temp%isite
            write(*,'(a,3x,2x,f9.4,26x,f9.4)') ' Angle, deviation   ', temp%angle, temp%dist
            write(*,*) 
            temp => temp%next
         else
            exit loop_nodes
         endif
      enddo loop_nodes
   endif
enddo
enddo
!
!  Attempt to match vectors by rotation / displacement of rigid units
!
call three_inter_rigid(sigma_length)
!
end subroutine three_inter
!
!*******************************************************************************
!
subroutine three_inter_match(sigma_length)
!-
!  Attempt  to match negative peaks by shift / rotation of ridig units
!+
!
use metric_mod
!
use precision_mod
!
implicit none
!
real(kind=PREC_DP), intent(in) :: sigma_length
!
real(kind=PREC_DP), dimension(3), parameter :: NULLV = (/ 0.0D0, 0.0D0, 0.0D0/)
!
integer :: i,j,k   ! dummy loop indices
!
real(kind=PREC_DP), dimension(3)            :: vect
real(kind=PREC_DP), dimension(3)            :: wect
!
type(three_match), pointer :: head
type(three_match), pointer :: temp
type(three_match), pointer :: tail
!
nullify(head)
nullify(temp)
nullify(tail)
!
do k=2, 1, -1
do i=1, three_npeaks(k)                         ! Try to match positive/negative peaks
!  head = three_res(i,2)%ptr                    ! Point to first node for element i Positive maxima
!  temp = three_res(i,2)%ptr
   nullify(head)
   nullify(temp)
   do j= 1, three_npatt                         ! Test against Patterson vectors
      if(abs(three_peaks(4,i,k)- three_patt(4,j))< sigma_length) then
         three_matched_xtal (j,k) = .true.
         three_matched_3dpdf(i,k) = .true.
         allocate(temp)                          ! Make a new node
         temp%isite(:) = three_patt_isite(:,j)
         vect = three_peaks(1:3,i,k)
         wect = three_patt(1:3,j)
         temp%dist  = do_blen(.true., vect, NULLV) - do_blen(.true., wect, NULLV)
         temp%angle = do_bang(.true., vect, NULLV, wect)
         nullify(temp%next)
         if(.not. associated(head)) then
            head => temp
            tail => temp
            three_res(i,k)%ptr => temp
         else
            tail%next => temp
         endif
      endif
   enddo
enddo
enddo
!
end subroutine three_inter_match
!
!*******************************************************************************
!
subroutine three_inter_rigid(sigma_length)
!-
!  Attempt immediate match of an experimental vector to any Patterson vector
!  in the unit cell
!+
!
!use chem_aver_mod
use crystal_mod
use chem_mod
use metric_mod
use molecule_mod
use symm_mod
!
use precision_mod
!
implicit none
!
real(kind=PREC_DP), intent(in) :: sigma_length
!
!real(kind=PREC_SP), dimension(3), parameter :: NULLV = (/ 0.0D0, 0.0D0, 0.0D0/)
!
!integer, dimension(4) :: ibest
integer :: i,j,l ! dummy loop indices
integer :: pmc, pmg  ! Loop indices plus minus crystal and ghost
integer, dimension(2) :: ibest
integer :: jmole
integer :: iatom
!
real(kind=PREC_DP)                          :: quality
real(kind=PREC_DP)                          :: optimum
real(kind=PREC_DP), dimension(3)            :: vect
real(kind=PREC_DP), dimension(3)            :: wect
real(kind=PREC_DP), dimension(:,:), allocatable :: extra
!real(kind=PREC_SP)                          :: dist
!real(kind=PREC_SP)                          :: angle
!real(kind=PREC_SP)                          :: sigma 
!real(kind=PREC_SP)                          :: short 
!
!type(three_match), pointer :: head
type(three_match), pointer :: temp
!type(three_match), pointer :: tail
!
!nullify(head)
nullify(temp)
!nullify(tail)
!sigma = 0.002
!short = 1e9
!ibest = 0
!
optimum = 1.E9
quality = 1.E9
ibest = 0
pmg = 1
pmc = 1
!write(*,*) ' START RIGID '
!
!  Loop over all positive peaks if they correspond to a rotated Patterson vector
!  Align molecule accordingly
!
loop_npeaks: do i=2, three_npeaks(2)                         ! Loop over positive peaks
   temp => three_res(i,2)%ptr
   if(.not.associated(temp)) cycle loop_npeaks                        ! This peak doesn't corresponds to a Patterson vector
   if(temp%angle < 15.0 .or. temp%angle > 165.0) cycle loop_npeaks    ! If rotation is not in significant range 
!write(*,'(a, i3,3(2x,f9.4),2(4x,f9.4),2i4)') ' Positive peak', i, three_peaks(1:4,i,2), temp%angle, temp%isite
   jmole = cr_mole(temp%isite(1))
   if(jmole==0) cycle loop_npeaks
!  loop_pmc: do pmc=1,1 !-1,1,2                                   ! Loop over minus plus Patterson
   vect = pmc*(cr_pos(:,temp%isite(2)) - cr_pos(:, temp%isite(1))) ! take crystal Patterson vector at +- 
!     loop_pmg: do pmg=1,1 !! -1,1,2                                   ! Loop over minus plus Patterson
   wect = pmg*three_peaks(1:3,i,2)                       ! Take extra Patterson vector
!write(*,'(a,2(2x,(3(2x,f9.4))))') ' Patt, extra ', vect, wect
   loop_iatom: do j=1,mole_len(jmole)                    ! Loop over all atoms in current molecule
      iatom = mole_cont(mole_off(jmole)+j)
      if(three_matched_site(iatom,2)) cycle loop_iatom
      call three_setup_rot(3, vect, wect, iatom, jmole, extra, .false.)
      call three_test_neg(iatom, ubound(extra,2),extra, quality)   ! Test if the extra match patterson vectors at negative 3DPDF poits
!write(*,*) ' unmatched atom ', iatom, ' peak ', i, ' qual ', quality, optimum
      if(quality<optimum) then
         optimum = quality
         ibest(1) = i
         ibest(2) = j
      endif
   enddo loop_iatom
!     enddo loop_pmg
!  enddo loop_pmc
enddo loop_npeaks
!
!write(*,*) ' Optimum solution ', ibest, optimum
temp => three_res(ibest(1),2)%ptr
jmole = cr_mole(temp%isite(1))
vect = pmc*(cr_pos(:,temp%isite(2)) - cr_pos(:, temp%isite(1))) ! take crystal Patterson vector at +
wect = pmg*three_peaks(1:3,ibest(1),2)                       ! Take extra Patterson vector
iatom = mole_cont(mole_off(jmole)+ibest(2))
!
call three_setup_rot(3, vect, wect, iatom, jmole, extra, .false.)
!
write(*,*) 
write(*,'(a)')     ' Good solution to negative peaks with'
write(*,'(a,2i5)') ' Molecule rotated: number, fixed atom', jmole, iatom
write(*,'(a,3(2x,f9.4))') ' Rotation axis  direct space         ', sym_uvw
write(*,'(a,3(2x,f9.4))') ' Rotation axis  reciprocal space     ', sym_hkl
write(*,'(a,3(2x,f9.4))') ' Rotation angle                      ', sym_angle
l = 0
loop_inmole: do i=1, ubound(extra,2) + 1
   if(mole_cont(mole_off(jmole)+i)==iatom) cycle loop_inmole
   j = mole_cont(mole_off(jmole)+i)
   l = l + 1
   write(*,'(a,i6, 3(2x,f9.4))') ' Atom number; new position     ', j, extra(1:3,l)
enddo loop_inmole
!write(*,*) ' RIGID DONE '
!
end subroutine three_inter_rigid
!
!*******************************************************************************
!
subroutine three_setup_rot(NDIM, vect, wect, iatom, jmole, extra, lout)
!-
!  Set up the rotation matrix to rotate vect into wekt
!+
!
use crystal_mod
use chem_mod
use metric_mod
use molecule_mod
use symm_menu ,only:symm_show
use symm_mod
use symm_sup_mod
!
use param_mod
use param_mod
use precision_mod
!
implicit none
!
integer, intent(in) :: NDIM              ! Vector dimensions
real(kind=PREC_DP), dimension(3), intent(in) :: vect
real(kind=PREC_DP), dimension(3), intent(in) :: wect
integer                         , intent(in) :: iatom      ! Atom at origin
integer                         , intent(in) :: jmole      ! Atom at origin
real(kind=PREC_DP), dimension(:,:), allocatable :: extra
logical, intent(in) :: lout
!
real(kind=PREC_DP), dimension(3), parameter :: NULLV = (/ 0.0D0, 0.0D0, 0.0D0/)
!
integer :: k,l    ! Dummy loop index
real(kind=PREC_DP), dimension(3) :: w
real(kind=PREC_DP), dimension(3) :: pos
!
!
call vekprod (vect, wect, w, cr_eps, cr_rten)
!
sym_uvw    = w
sym_angle  = do_bang(.true., vect, NULLV, wect)
sym_orig   = chem_ave_pos(:, iatom)                    ! CAREFULL WITH AVERAGE VERSUS UNIT CELL!
sym_orig_type  = 0
sym_trans  = 0.0d0
sym_type   = .true.             ! Proper rotation
sym_power  = 1
sym_sel_atom   = .true.
sym_start  = 1
sym_end    = 1
sym_latom  = .true.
!
!
call symm_setup
if(lout) call symm_show
!
if(allocated(extra)) deallocate(extra)
allocate(extra(3,mole_len(jmole)-1))
l = 0
loop_inmole: do k=1,mole_len(jmole)               ! Loop over all atoms in molecule to generate new sites
   if(mole_cont(mole_off(jmole)+k)==iatom) cycle loop_inmole
   pos = chem_ave_pos(:,mole_cont(mole_off(jmole)+k))
   call symm_ca_single(pos, .true., .false.  )
   l = l + 1
   extra(:,l) = res_para(1:3)
!write(*,'(a,3(2x,f9.4))') ' Ghost positions ', extra(:,l)
enddo loop_inmole
!
end subroutine three_setup_rot
!
!*******************************************************************************
!
subroutine three_test_neg(iatom, NDIM, extra, quality)                ! Test if the extra match patterson vectors at negative 3DPDF poits
!-
!   Generate Patterson vectors between real atoms and ghost atoms
!   and see if these match the negative positions
!+
!
use crystal_mod
use chem_mod
use metric_mod
!
use precision_mod
use trig_degree_mod
!
implicit none
!
integer                              , intent(in) :: iatom      ! Atom at origin
integer                              , intent(in)  :: NDIM
real(kind=PREC_DP), dimension(3,NDIM), intent(in)  :: extra
real(kind=PREC_DP)                   , intent(out) :: quality
!
real(kind=PREC_DP), dimension(3), parameter :: NULLV = (/ 0.0D0, 0.0D0, 0.0D0/)
!
integer :: i,j,k  ! Dummy loop indices
integer :: iopt
real(kind=PREC_DP), dimension(3) :: vect
real(kind=PREC_DP), dimension(3) :: wect
real(kind=PREC_DP)               :: dist
real(kind=PREC_DP)               :: dmin
real(kind=PREC_DP)               :: angle
!
quality = 0.0
         dmin = 1.e9
do i=1, NDIM
   loop_ncatom: do j=1, cr_ncatoms
      if(j==iatom) cycle loop_ncatom
      vect = extra(:,i) - chem_ave_pos(:,j)
!write(*,'(a,2i4,3(2x,f9.4))') ' Virtual vector ', i, j, vect
      iopt = 0
      angle = 360.
      dmin = 1.e9
      do k=1, three_npeaks(1) 
         dist = abs(do_blen(.true., vect, NULLV)-three_peaks(4,k,1))
!write(*,'(a,11x,3(2x, f9.4),2x, 3(2x,f9.4))') ' teste vers> ', three_peaks(1:3,k,1), dmin, dist, 0.00000
         if(dist<dmin) then       ! Found better match
            wect = three_peaks(1:3,k,1)
            angle = do_bang(.true., vect, NULLV, wect)
            iopt = k
            dmin = dist
         endif
      enddo
!
!write(*,'(a,11x,3(2x, f9.4),2x, 3(2x,f9.4),i4)') ' BEST        ', &
!   three_peaks(1:3,iopt,1), dmin, angle, sind(angle), iopt
!write(*,*) 'Length ', three_peaks(4,iopt,1)
!
      quality = quality + dmin**2 + (sind(angle)/three_peaks(4,iopt,1))**2
   enddo loop_ncatom
enddo
!write(*,*) ' Quality ', quality
!
end subroutine three_test_neg
!
!*******************************************************************************
!
subroutine three_inter_rigid_b(sigma_length)
!-
!  Attempt immediate match of an experimental vector to any Patterson vector
!  in the unit cell
!+
!
!use chem_aver_mod
use crystal_mod
use chem_mod
use metric_mod
use molecule_mod
!
use precision_mod
!
implicit none
!
real(kind=PREC_DP), intent(in) :: sigma_length
!
real(kind=PREC_DP), dimension(3), parameter :: NULLV = (/ 0.0D0, 0.0D0, 0.0D0/)
!
integer, dimension(4) :: ibest
integer :: i,j,k,l ! dummy loop indices
integer :: jmole
integer :: iatom
!
real(kind=PREC_DP), dimension(3)            :: vect
real(kind=PREC_DP), dimension(3)            :: wect
real(kind=PREC_DP)                          :: dist
real(kind=PREC_DP)                          :: angle
real(kind=PREC_DP)                          :: sigma 
real(kind=PREC_DP)                          :: short 
!
type(three_match), pointer :: head
type(three_match), pointer :: temp
type(three_match), pointer :: tail
!
nullify(head)
nullify(temp)
nullify(tail)
sigma = 0.002
short = 1e9
ibest = 0
!
!write(*,*) ' START RIGID '
!
!  Add a negative vector as +- to all sites in unit cell
!
loop_main1: do i=1, three_npeaks(1)                         ! Try to match negative peaks
   loop_pm: do k= -1, 1, 2                                  ! add with sign -1; +1
      loop_sites: do j=1, cr_ncatoms
         if(.not.three_matched_site(j,2)) cycle loop_sites
         vect = k*three_peaks(1:3, i,1) + chem_ave_pos(1:3,j)
         do l=1,three_npatt                                 ! Compare length to all Patterson vectors
            dist = abs(do_blen(.true., NULLV, vect)-three_patt(4,l))
         if(dist<sigma_length) then                         ! Possible match
            wect = three_patt(1:3,j)
            angle = do_bang(.true., vect, NULLV, wect)
            if(dist<short) then
               ibest(1) = i
               ibest(2) = k
               ibest(3) = j
               ibest(4) = l
               short = dist
            endif
write(*,'(a, i3, 2(a,3(2x,f9.4)))') ' added ', k,' * ',three_peaks(1:3, i,1), ' + ', chem_ave_pos(1:3,j)
write(*,'(a,4i4)') '   Indices   ', i,k,j,l
write(*,'(a, 6x,3(2x,f9.4),4x,f9.4 )')       '   Vect', vect, do_blen(.true., NULLV, vect)
write(*,'(a, 2(2x,f9.4))') '   Diff angle', dist, angle
         endif
         enddo
      enddo loop_sites
   enddo loop_pm
enddo loop_main1
!
write(*,*) ' best solution ', ibest
!
! Is the site within a molecule ?
!
write(*,*) ' Best site ', ibest(3), cr_prop(ibest(3)), cr_mole(ibest(3))
jmole = cr_mole(ibest(3))
if_inmole: if(jmole>0) then             ! The site is involved in a molecule
   do i=1,mole_len(jmole)
      iatom = mole_cont(mole_off(jmole)+i)
write(*,*) ' IATOM ', iatom, three_matched_site(iatom,2)
   enddo
endif if_inmole
!
write(*,*) ' DONE  RIGID '
!
end subroutine three_inter_rigid_b
!
!
!*******************************************************************************
!
subroutine three_clear_nodes
!-
! Clear and deallocate the list of nodes
!+
integer :: i, j                             ! Dummy loop indices
type(three_match), pointer :: head
type(three_match), pointer :: temp
!
if(allocated(three_res)) then
   do j=1,2
   loop_entries:do i=1,ubound(three_res,1)
      head => three_res(i,j)%ptr              ! Point to first node of element i
      temp => three_res(i,j)%ptr
      do
         if(associated(temp)) then          ! Further nodes exist
            temp => temp%next               ! Point to next node
            deallocate(head)                ! Deallocate current node
            head => temp                    ! Follow with current node
         else                               ! No more node
            if(associated(head)) deallocate(head) ! Clear last node
            cycle loop_entries
         endif
      enddo
   enddo loop_entries
   enddo
endif
nullify(head)
nullify(temp)
!
end subroutine three_clear_nodes
!
!*******************************************************************************
!
subroutine three_reset
!-
!  Reset all variables to system start
!+
if(allocated(three_map)     ) deallocate(three_map)
if(allocated(three_map_temp)) deallocate(three_map_temp)
if(allocated(three_peaks)   ) deallocate(three_peaks)
three_max_peaks =  1
three_npeaks    =  1
three_dims  = 0
three_llims = 0.0D0
three_steps = 0.0D0
!
end subroutine three_reset
!
!*******************************************************************************
!
end module discus_3dpdf_mod
