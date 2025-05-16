module trans_to_short_mod
!
! Save old structure transform current structure such that 
! the shortest step in powder complete is parallel to c*
!
use precision_mod
!
private
public trans_to_short
public trans_to_short_reset
!
character(len=PREC_STRING)      , save :: trn_structure     ! Original structure
!character(len=PREC_STRING)      , save :: trn_stacklist     ! Original stacking fault list
!character(len=PREC_STRING)      , save :: trn_stacksimple   ! Original stacking fault list
character(len=PREC_STRING), dimension(:)    , allocatable, save :: trn_stacklayer ! Original stacking fault layers
real(kind=PREC_DP)        , dimension(3)                 , save :: trn_pow_hkl_del   ! Original hkl steps
real(kind=PREC_DP)        , dimension(:,:)  , allocatable, save :: trn_stack_origins ! Original stacking fault origins
real(kind=PREC_DP)        , dimension(:,:,:), allocatable, save :: trn_stack_trans   ! Original stacking fault translations
logical                         , save :: trn_changed       ! Stucture was changed
!
contains
!
!*******************************************************************************
!
subroutine trans_to_short(pow_hkl_del, pow_four_mode_is_stack, pow_four_mode_is_four, crystal_is_stack)
!-
! Main transformation routine
!+
!
use crystal_mod
use discus_allocate_appl_mod, only: alloc_transfrm
use modify_mod, ONLY: atom_select
use read_internal_mod
use save_menu, ONLY: save_internal_backup
use stack_mod
use transform_menu
use transfrm_mod
!
use precision_mod
!
implicit none
!
real(kind=PREC_DP), dimension(3), intent(inout) :: pow_hkl_del
logical                         , intent(in)    :: pow_four_mode_is_stack
logical                         , intent(in)    :: pow_four_mode_is_four 
logical                         , intent(in)    :: crystal_is_stack
!
logical, parameter            :: lold = .FALSE.
integer, parameter            :: MAXMASK=4
character(len=PREC_STRING)    :: line   ! Dummy line
character(len=PREC_STRING)    :: temp_file
logical, dimension(0:MAXMASK) :: uni_mask
integer                       :: length ! dummy length
integer                       :: i, j   ! Dummy index
real(kind=PREC_DP), dimension(4) ::  usym(4)  ! Dummy vectors
real(kind=PREC_DP), dimension(4) ::  ures(4)  ! Dummy vectors
!
!if(pow_four_mode_is_stack .or. crystal_is_stack) then       ! Stacking fault mode
!write(*,*) st_nlayer, st_nchem, st_ntypes
!write(*,*) ' LAYER '
!do length=1, st_ntypes
!  write(*,'(i3, 2x, a)') length, st_layer(length)(1:len_trim(st_layer(length)))
!enddo
!write(*,*) ' LAYER '
!do length=1, st_nchem
!  write(*,'(i3, 2x,a)') length, st_layer_c(length)(1:len_trim(st_layer_c(length)))
!enddo
!read(*,*) length
!endif
!
trn_changed = .false.
!
if(    pow_hkl_del(3)<=pow_hkl_del(2) .and. pow_hkl_del(2)<=pow_hkl_del(1)) then ! 3<2<1 nothing to do
!  write(*,*) 'Sequence 3 < 2 < 1'
  return 
elseif(pow_hkl_del(3)<=pow_hkl_del(1) .and. pow_hkl_del(1)< pow_hkl_del(2)) then ! 3<1<2 nothing to do
!  write(*,*) 'Sequence 3 < 1 < 2'
  return 
endif
!
!  Backup original crystal structure
!
trn_structure   = 'internal.backup.main' ! internal_backup 
call save_internal_backup(trn_structure) ! Save as full backup
!
trn_changed = .true.
!
trn_pow_hkl_del = pow_hkl_del        ! Save users original steps
!
tran_g = 0.0_PREC_DP
!
if(    pow_hkl_del(2)< pow_hkl_del(3) .and. pow_hkl_del(3)< pow_hkl_del(1)) then ! 2 < 3 < 1
  tran_g(1,1) = -1.0_PREC_DP
  tran_g(2,3) =  1.0_PREC_DP
  tran_g(3,2) =  1.0_PREC_DP
  pow_hkl_del(1) = trn_pow_hkl_del(1)
  pow_hkl_del(2) = trn_pow_hkl_del(3)
  pow_hkl_del(3) = trn_pow_hkl_del(2)
!  write(*,*) 'Sequence 2 < 3 < 1', pow_hkl_del
elseif(pow_hkl_del(2)< pow_hkl_del(3) .and. pow_hkl_del(3)==pow_hkl_del(1)) then ! 2 < 3 == 1
  tran_g(1,3) =  1.0_PREC_DP
  tran_g(2,1) =  1.0_PREC_DP
  tran_g(3,2) =  1.0_PREC_DP
  pow_hkl_del(1) = trn_pow_hkl_del(3)
  pow_hkl_del(2) = trn_pow_hkl_del(1)
  pow_hkl_del(3) = trn_pow_hkl_del(2)
!  write(*,*) 'Sequence 2 < 3 == 1', pow_hkl_del
elseif(pow_hkl_del(2)< pow_hkl_del(1) .and. pow_hkl_del(1)< pow_hkl_del(3)) then ! 2 < 1 < 3
  tran_g(1,3) =  1.0_PREC_DP
  tran_g(2,1) =  1.0_PREC_DP
  tran_g(3,2) =  1.0_PREC_DP
  pow_hkl_del(1) = trn_pow_hkl_del(3)
  pow_hkl_del(2) = trn_pow_hkl_del(1)
  pow_hkl_del(3) = trn_pow_hkl_del(2)
!  write(*,*) 'Sequence 2 < 3 < 1', pow_hkl_del
elseif(pow_hkl_del(1)< pow_hkl_del(2) .and. pow_hkl_del(2)< pow_hkl_del(3)) then ! 1 < 2 < 3
  tran_g(1,3) =  1.0_PREC_DP
  tran_g(2,2) = -1.0_PREC_DP
  tran_g(3,1) =  1.0_PREC_DP
  pow_hkl_del(1) = trn_pow_hkl_del(3)
  pow_hkl_del(2) = trn_pow_hkl_del(2)
  pow_hkl_del(3) = trn_pow_hkl_del(1)
!  write(*,*) 'Sequence 1 < 2 < 3', pow_hkl_del
elseif(pow_hkl_del(1)< pow_hkl_del(2) .and. pow_hkl_del(2)==pow_hkl_del(3)) then ! 1 < 2 == 3
  tran_g(1,2) =  1.0_PREC_DP
  tran_g(2,3) =  1.0_PREC_DP
  tran_g(3,1) =  1.0_PREC_DP
  pow_hkl_del(1) = trn_pow_hkl_del(2)
  pow_hkl_del(2) = trn_pow_hkl_del(3)
  pow_hkl_del(3) = trn_pow_hkl_del(1)
!  write(*,*) 'Sequence 1 < 2 == 3', pow_hkl_del
elseif(pow_hkl_del(1)< pow_hkl_del(3) .and. pow_hkl_del(3)< pow_hkl_del(2)) then ! 1 < 3 < 2
  tran_g(1,2) =  1.0_PREC_DP
  tran_g(2,3) =  1.0_PREC_DP
  tran_g(3,1) =  1.0_PREC_DP
  pow_hkl_del(1) = trn_pow_hkl_del(2)
  pow_hkl_del(2) = trn_pow_hkl_del(3)
  pow_hkl_del(3) = trn_pow_hkl_del(1)
!  write(*,*) 'Sequence 1 < 3 < 2', pow_hkl_del
endif
!
tran_inp = TRAN_INP_G
call tran_setup
tran_start = 1  ! Include all atoms
tran_end = - 1
!
if(pow_four_mode_is_stack) then ! .or. crystal_is_stack) then       ! Stacking fault mode
!write(*,*) st_nlayer, st_nchem, st_ntypes
!read(*,*) length
!
!
   allocate(trn_stacklayer(1:st_nchem))
   do i=1, st_nchem
!
      trn_stacklayer(i) = st_layer_c(i)                ! actual layers 
      call readstru_internal(MAXMASK, trn_stacklayer(i), uni_mask)
      temp_file = trn_stacklayer(i)(1:len_trim(trn_stacklayer(i))) // '.backup'
      call save_internal_backup(temp_file) ! Save as full backup
!write(*,*) ' TRANS  SOURCE ' , trn_stacklayer(i)(1:len_trim(trn_stacklayer(i)))
!write(*,*) ' TRANS  TARGET ' , temp_file(1:len_trim(temp_file))
!
      call alloc_transfrm ( cr_nscat, cr_ncatoms )
      line ='all'     ! Select all atom types
      length = 3
      call atom_select(line, length, 0, TRAN_MAXSCAT, tran_latom, tran_lsite, 0, &
                       TRAN_MAXSITE, tran_sel_atom, lold,.TRUE.)
!do j=1, min(cr_natoms,10)
!  write(*,'(a,i3, 3f11.6)') cr_at_lis(cr_iscat(1,j)), cr_iscat(1,j), cr_pos(:,j)
!enddo
!   call transfrm_show
!read(*,*) length
      call tran_op(.FALSE.)  ! Transform silently
!do j=1, min(cr_natoms,10)
!  write(*,'(a,i3, 3f11.6)') cr_at_lis(cr_iscat(1,j)), cr_iscat(1,j), cr_pos(:,j)
!enddo
!write(*,*) ' SAVE AS ', trn_stacklayer(i)(1:len_trim(trn_stacklayer(i)))
!read(*,*) length
      call save_internal_backup(trn_stacklayer(i)) ! Save as full backup
!
   enddo
!
   allocate(trn_stack_origins(3, 1:st_nlayer))
   trn_stack_origins = st_origin
   usym(4) = 1.0_PREC_DP
   do i=1, st_nlayer   ! Transform all stacking fault layers
!
      usym(1:3) = st_origin(:,i)
      ures      = matmul(tran_f, usym)
      st_origin(:,i) = ures(1:3)
!if(i<10 .or. i>st_nlayer-10) then
!write(*,'(a,3f12.6, a, 3f12.6)'), ' Origin ', usym(1:3), ' ==> ', st_origin(:,i)
!endif
!
   enddo
endif
!
if(pow_four_mode_is_stack .or. cr_is_stack) then   ! Stacking Falut mode or final crystal
!write(*,*) ' TRANSFORMING TRANSITION '
   allocate(trn_stack_trans(st_ntypes,st_ntypes,3))
   trn_stack_trans(1:st_ntypes,1:st_ntypes,:) = st_trans(1:st_ntypes,1:st_ntypes,:)
   do i=1, st_ntypes
      do j=1, st_ntypes
         usym(1:3) = st_trans(i,j,:)
         ures      = matmul(tran_f, usym)
         st_trans(i,j,:) = ures(1:3)
!write(*,'(a,2i3,3f12.6, a, 3f12.6)') ' TRANS ', i,j, usym(1:3), ' ==> ', st_trans(i,j,:)
      enddo
   enddo
!read(*,*) length
endif
!
if(pow_four_mode_is_four .or. cr_is_stack) then  ! Fourier mode or final crystal build from stacks
!
!write(*,*) ' TRANSFORMING STRUCTURE  '
   call alloc_transfrm ( cr_nscat, cr_ncatoms )
   line ='all'     ! Select all atom types
   length = 3
   call atom_select(line, length, 0, TRAN_MAXSCAT, tran_latom, tran_lsite, 0, &
                    TRAN_MAXSITE, tran_sel_atom, lold,.TRUE.)
!   call transfrm_show
    call tran_op(.FALSE.)  ! Transform silently
endif
!
end subroutine trans_to_short
!
!*******************************************************************************
!
subroutine trans_to_short_reset(pow_hkl_del, pow_four_mode_is_stack, pow_four_mode_is_four, crystal_is_stack)
!-
! Restore transformation routine
!+
!
use crystal_mod
use read_internal_mod
use save_menu, ONLY: save_internal_backup
use stack_mod
!
use precision_mod
!
implicit none
!
real(kind=PREC_DP), dimension(3), intent(inout) :: pow_hkl_del
logical                         , intent(in)    :: pow_four_mode_is_stack
logical                         , intent(in)    :: pow_four_mode_is_four 
logical                         , intent(in)    :: crystal_is_stack
!
integer, parameter :: MAXMASK=4
!
character(len=PREC_STRING) :: backup_file
integer                       :: i, j
logical, dimension(0:MAXMASK) :: uni_mask
!
cond_changed:if(trn_changed) then        ! Structure was changed, restore backup
!
!write(*,*) ' RESTORING STRUCTURE ', pow_four_mode_is_stack, pow_four_mode_is_four, crystal_is_stack
   cond_stack: if(pow_four_mode_is_stack) then ! .or. cr_is_stack) then       ! Stacking fault mode
!write(*,*) ' RESTORING STACK '
      st_origin = trn_stack_origins    ! Restore origins
!do i=1, 10
!write(*,'(a,3f12.6, a, 3f12.6)'), ' Origin ', st_origin(:,i)
!enddo
!read(*,*) j
!
      uni_mask(0)   = .false.
      uni_mask(1:3) = .true.
      uni_mask(4)   = .false.
!
      do i=1, st_nchem                 ! Restore all layers 
         backup_file = trn_stacklayer(i)(1:len_trim(trn_stacklayer(i))) // '.backup'
!write(*,*) ' BACKUP SOURCE ' , backup_file(1:len_trim(backup_file))
!write(*,*) ' BACKUP TARGET ' , trn_stacklayer(i)(1:len_trim(trn_stacklayer(i)))
         call readstru_internal(MAXMASK, backup_file, uni_mask)   ! Read  backup file
         call save_internal_backup(trn_stacklayer(i)) ! Save as full backup
!        call store_remove_single(trn_stacklayer(i), ier_num)
!do j=1, min(cr_natoms,10)
!  write(*,'(a,i3, 3f11.6)') cr_at_lis(cr_iscat(1,j)), cr_iscat(1,j), cr_pos(:,j)
!enddo
!read(*,*) j
      enddo
!
      deallocate(trn_stacklayer)
      deallocate(trn_stack_origins)
!
   endif cond_stack
!
   cond_trans: if(pow_four_mode_is_stack .or. cr_is_stack) then   ! Stacking Falut mode or final crystal
!write(*,*) ' RESTORE TRANSITION ', allocated(trn_stack_trans), allocated(st_trans)
      st_trans  = trn_stack_trans      ! Restore transition vectors
      deallocate(trn_stack_trans)
!do i=1, st_ntypes
!do j=1, st_ntypes
!   write(*,'(a,2i3, 3f10.6)') ' Restored trans ', i,j, st_trans(i,j,:)
!enddo
!enddo
   endif cond_trans
!
   backup_file = trn_structure(1:len_trim(trn_structure)) ! 
   call readstru_internal(MAXMASK, backup_file, uni_mask)   ! Read  backup file
   pow_hkl_del = trn_pow_hkl_del
!  call store_remove_single(trn_structure, ier_num)
!write(*,*) ' RESTORING IS DONE'
endif cond_changed
!
end subroutine trans_to_short_reset
!
!*******************************************************************************
!
end module trans_to_short_mod
