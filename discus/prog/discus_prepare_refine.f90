module prep_refine_mod
!
use precision_mod
!
private
public write_refine_single_part1
public write_refine_single_part2
public write_refine_single_part3
public write_refine_single_part4
public write_refine_powder_part1
public write_refine_powder_part4
public write_refine_pdf_part4
public write_diffev_generic
public write_diffev_single_part1
public write_diffev_powder_part1
public write_diffev_single_part2
public write_diffev_single_part3
public write_diffev_single_part4
public write_diffev_powder_part4
public write_diffev_pdf_part1
public write_diffev_pdf_part4
!
integer, parameter :: IDI_MA  = 36    ! DISCUS_MAIN MACRO PART
integer, parameter :: IRE_MA  = 37    ! REFINE_MAIN PART
integer, parameter :: IDF_MA  = 38    ! DIFFEV_MAIN PART
integer, parameter :: IDI_PO  = 39    ! DISCUS_POWDER
integer, parameter :: IDF_NP  = 40    ! DIFFEV_NEWPARA
integer, parameter :: IRE_NP  = 41    ! DIFFEV_NEWPARA
integer, parameter :: IDF_SI  = 42    ! DIFFEV_SINGLE
integer, parameter :: IDF_PO  = 43    ! DIFFEV_POWDER
integer, parameter :: IDF_BS  = 44    ! diffev_best main
integer, parameter :: IDF_BS2 = 45    ! diffev_best part 2
!
character(len=PREC_STRING) :: discus_file
character(len=PREC_STRING) ::    hkl_file
character(len=PREC_STRING) ::    fcf_file
integer :: diffev_dimension
logical :: lexist_hkl_file
logical :: lexist_fcf_file
!
logical :: is_opened
!
contains
!
!*****7**************************************************************** 
!
subroutine write_refine_single_part1(substance, ilist, fv_1)
!
!  Write main discus macro; single crystal PART 1
!
!
use precision_mod
!
implicit none
!
character(len=*)   , intent(in) :: substance
integer            , intent(in) :: ilist
real(kind=PREC_DP) , intent(in) :: fv_1
!
character(len=PREC_STRING) :: refine_file
integer :: i
!
lexist_hkl_file = .false.
lexist_fcf_file = .false.
i= len_trim(substance)
discus_file = substance(1:i)//'_main.mac'
if(ilist == 0) then
   hkl_file = substance(1:i)//'.hkl'
   fcf_file = ' '
   inquire(file=hkl_file, exist=lexist_hkl_file)
else
   hkl_file = substance(1:i)//'.hkl'
   fcf_file = substance(1:i)//'.fcf'
   inquire(file=hkl_file, exist=lexist_hkl_file)
   inquire(file=fcf_file, exist=lexist_fcf_file)
   if(lexist_fcf_file) then
      hkl_file = substance(1:i)//'_merged.hkl'    ! hkl file will be merged
      lexist_hkl_file = .false.          ! FCF exists, do not use original HKL file
   endif
endif
open(unit=IDI_MA, file=discus_file, status='unknown')
open(unit=IRE_NP, file='newpara.mac', status='unknown')
!
write(IDI_MA,'(a )') 'branch discus'
write(IDI_MA,'(a )') '#'
write(IDI_MA,'(a )') 'read'
write(IDI_MA,'(3a)') '  stru ', substance(1:i), '.cell'
!
if(lexist_hkl_file) then       ! An hkl file exists
   refine_file = 'refine_'//substance(1:i)//'_hkl.mac'
   open(unit=IRE_MA, file=refine_file, status='unknown')
   write(IRE_MA,'(a )') 'refine'
   write(IRE_MA,'(a )') 'rese'
   write(IRE_MA,'(2a)') 'data hklf4, ',hkl_file(1:len_trim(hkl_file))
elseif(lexist_fcf_file) then    ! An fcf file exists, we will have merged hkl file
   refine_file = 'refine_'//substance(1:i)//'_merged.mac'
   open(unit=IRE_MA, file=refine_file, status='unknown')
   write(IRE_MA,'(a )') 'refine'
   write(IRE_MA,'(a )') 'rese'
   write(IRE_MA,'(2a)') 'data hklf4, ',hkl_file(1:len_trim(hkl_file))
else                           !NO hkl fikel assume an h5 file will be used
   refine_file = 'refine_'//substance(1:i)//'_dat.mac'
   open(unit=IRE_MA, file=refine_file, status='unknown')
   write(IRE_MA,'(a )') 'refine'
   write(IRE_MA,'(a )') 'rese'
   write(IRE_MA,'(3a)') 'data h5, DATA/',hkl_file(1:len_trim(hkl_file)-4), '.h5'
endif
write(IRE_MA,'(a        )') '#'
write(IRE_MA,'(a        )') '@newpara.mac'
write(IRE_MA,'(a        )') '#'
!
write(IRE_NP,'(a,f9.4,a )') 'newpara P_scale, value:', fv_1 , ', points:3, shift:0.001, status:free'
!
end subroutine write_refine_single_part1
!
!*****7**************************************************************** 
!
subroutine write_refine_single_part2(j, l, jj, iscat, lcontent, natoms, c_atom, uij_l, c_flag)
!
!  Write main discus macro; single crystal PART 2
!  Writes Uij or Ueqv
!  can be used for SINGLE and POWDER
!
!
use precision_mod
!
implicit none
!
integer            , intent(in) :: j
integer            , intent(in) :: l
integer            , intent(in) :: jj
integer            , intent(in) :: iscat
integer            , intent(in) :: lcontent
integer            , intent(in) :: natoms
character(len=4), dimension(lcontent), intent(in)  :: c_atom    ! Atom types 
real(kind=PREC_DP), dimension(6,natoms, natoms), intent(in) :: uij_l
character(len=*)   , intent(in) :: c_flag
!
integer :: m, k
!
if(l==1) then
   write(IDI_MA,'(a,i3,4a)') 'anis type:', jj,', values:[', 'U_',c_atom(iscat)(1:len_trim(c_atom(iscat))),'_1]'
   write(IRE_NP,'(4a,g15.8e3,2a)') 'newpara U_',c_atom(iscat)(1:len_trim(c_atom(iscat))),'_1', &
                 ', value:', uij_l(1,iscat, j), ', points:3, shift:0.001, status:', c_flag 
elseif(l==6) then
   write(IDI_MA,'(a,i3,a,5(3a,i1.1,a2),3a)') 'anis type:', jj,', values:[', &
        ('U_',c_atom(iscat)(1:len_trim(c_atom(iscat))),'_',k,', ',k=1,5), &
         'U_',c_atom(iscat)(1:len_trim(c_atom(iscat))),'_6]'
   do m=1,6
      write(IRE_NP,'(3a,i1.1,a,g15.8e3,2a)') 'newpara U_',c_atom(iscat)(1:len_trim(c_atom(iscat))),'_',m, &
                ', value:', uij_l(m,iscat, j), ', points:3, shift:0.001, status:', c_flag 
   enddo
endif
!
end subroutine write_refine_single_part2
!
!*****7**************************************************************** 
!
subroutine write_refine_single_part3(line       , iatom, iscat, lcontent, natoms, c_atom, posit, c_flag)
!
!  Write main discus macro; single crystal PART 3
!  can be used for SINGLE and POWDER
!
use spcgr_apply, only: get_wyckoff
use wyckoff_mod
!
use precision_mod
!
implicit none
!
character(len=*)   , intent(in) :: line
integer            , intent(in) :: iatom
integer            , intent(in) :: iscat
integer            , intent(in) :: lcontent
integer            , intent(in) :: natoms
character(len=4), dimension(lcontent), intent(in)  :: c_atom    ! Atom types 
real(kind=PREC_DP), dimension(3,natoms), intent(in) :: posit
character(len=*)   , intent(in) :: c_flag
!
character(len=1), dimension(3), parameter :: xyz_name = (/ 'x', 'y', 'z'/)
character(len=2), dimension(3), parameter :: xyz_appe = (/ '_x', '_y', '_z'/)
integer, parameter :: mode = 0
logical, parameter :: loutput = .false.
!
character(len=PREC_STRING) :: string
integer :: i, j, k, l
real(kind=PREC_DP), dimension(3) :: coordinates
!
coordinates = posit(:, iatom)
call get_wyckoff(coordinates, loutput, mode)
!
!     line = content(jc)(1:4)
i = len_trim(line)
loop_xyz: do j=1, 3         ! Loop over x, y, z
   if(wyc_fix(j)) then               ! Fixed or related
      if(all(wyc_fix_mat(j,1:3)==0)) cycle loop_xyz    ! Numerically fixed coordinate
!                                                      ! related coordinate
      write(string,'(2a1,i3,a4)') xyz_name(j),'[',iatom, '] = ' 
      l = 10   ! next free field in string
      do k =1, j-1
         if(wyc_fix_mat(j,k)/=0) then
            write(string(l:l+9), '(i2,3a1,i3,a1)') wyc_fix_mat(j,k), '*',xyz_name(k),'[',iatom, ']'
            l = l + 10
         endif
      enddo
      if(wyc_fix_mat(j,4)>0) then
         write(string(l:l+14), '(a2,f12.8)') ' +',    wyc_fix_mat(j,4)  / 24.0_prec_dp
      elseif(wyc_fix_mat(j,4)<0) then
         write(string(l:l+14), '(a2,f12.8)') ' -',abs(wyc_fix_mat(j,4)) / 24.0_prec_dp
      endif
      write(IDI_MA, '(a)') string(1:len_trim(string))
   else                                                ! Free coordinate
      write(IDI_MA,'(2a,i3,3a)') xyz_name(j),'[',iatom, '] = P_',line(1:i), xyz_appe(j)
      write(IRE_NP,'(4a,f12.8,2a)') 'newpara P_',c_atom(iscat)(1:len_trim(c_atom(iscat))), &
         xyz_appe(j),', value:',                                                        &
         posit(j, iatom), ', points:3, shift:0.001, status:', c_flag
   endif
!
enddo loop_xyz
!
end subroutine write_refine_single_part3
!
!*****7**************************************************************** 
!
subroutine write_refine_single_part4( &
        substance,  ilist, c_form, lambda, slambda, P_exti, hkl_max, P_ncycle)
!
!  Write main discus macro; single crystal PART 4
!  REFINE version
!
use blanks_mod
use precision_mod
!
implicit none
!
character(len=*)                , intent(in) :: substance
integer                         , intent(in) :: ilist
character(len=*)                , intent(in) :: c_form
character(len=*)                , intent(in) ::  lambda
character(len=*)                , intent(in) :: slambda
real(kind=PREC_DP)              , intent(in) :: P_exti
integer           , dimension(3), intent(inout) :: hkl_max
integer                         , intent(in) :: P_ncycle
!
character(len=PREC_STRING) :: aspher_file
character(len=PREC_STRING) :: line
character(len=PREC_STRING) :: line1
character(len=PREC_STRING) :: line2
integer :: i,j, k
integer :: length
integer :: ios
logical :: lexist
real(kind=PREC_DP) :: rlambda
!
read(lambda, *) rlambda
!
j = 0
i = len_trim(substance)
call single_main_part4(IDI_MA, .true., substance, c_form, rlambda, slambda, P_exti)
!
write(IDI_MA,'(a )') 'branch kuplot'
write(IDI_MA,'(a )') 'reset'
if(lexist_hkl_file) then                      ! Single HKL refinement
   write(IDI_MA,'(a )') 'load hklf4, calc.hkl'
elseif(lexist_fcf_file) then                      ! Single HKL refinement
   write(IDI_MA,'(a )') 'load hklf4, calc.hkl'
else
   write(IDI_MA,'(a )') 'load h5, calc.h5'
endif
write(IDI_MA,'(a )') 'exit  ! Back to DISCUS'
write(IDI_MA,'(a )') 'exit  ! Back to REFINE'
write(IDI_MA,'(a )') 'finished'
close(IDI_MA)
!
call single_k_inter_macro(IDI_MA, .true., substance)
!
if(P_exti> 0.00001) then
   write(IRE_NP,'(a,f8.5,a)') 'newpara P_exti, value:', P_exti, ', points:3, shift:0.001, status:free'
else
   write(IRE_NP,'(a,f8.5,a)') 'newpara P_exti, value:', P_exti, ', points:3, shift:0.001, status:fixed'
endif
write(IRE_MA,'(a,i5 )') 'set cycle, ', P_ncycle
write(IRE_MA,'(a )') 'set conver, status:on, dchi:0.050, chisq:1.10, pshift:2.0, conf:1.0, lambda:65000.'
write(IRE_MA,'(2a)') '@', discus_file(1:len_trim(discus_file))
write(IRE_MA,'( a)') 'branch kuplot'
if(lexist_hkl_file) then                      ! Single HKL refinement
   write(IRE_MA,'( a)') '  @k_fobs_fcalc.mac'
   write(IRE_MA,'(3a)') 'run ', discus_file(1:len_trim(discus_file)), ', plot:k_fobs_fcalc.mac'
elseif(lexist_fcf_file) then                      ! Single HKL refinement
   write(IRE_MA,'( a)') '  @k_fobs_fcalc.mac'
   write(IRE_MA,'(3a)') 'run ', discus_file(1:len_trim(discus_file)), ', plot:k_fobs_fcalc.mac'
else
   write(IRE_MA,'( a)') '  @k_inter.mac'
   write(IRE_MA,'(3a)') 'run ', discus_file(1:len_trim(discus_file)), ', plot:k_inter.mac'
endif
write(IRE_MA,'(a )') '#'
write(IRE_MA,'(a )') '@final_cell_cif.mac'
write(IRE_MA,'(a )') '#'
write(IRE_MA,'(a )') 'exit  ! Back to SUITE'
close(IRE_MA)
!
open(IRE_MA, file='final_cell_cif.mac', status='unknown')
i= len_trim(substance)
write(IRE_MA,'( a)') 'branch discus'
write(IRE_MA,'( a)') 'read'
write(IRE_MA,'(3a)') '  stru internal.' , substance(1:i), '.cell'
write(IRE_MA,'( a)') 'save'
write(IRE_MA,'(3a)') '  outf ', substance(1:i), '.cell'
write(IRE_MA,'( a)') '  write all'
write(IRE_MA,'( a)') '  run'
write(IRE_MA,'( a)') 'exit'
write(IRE_MA,'(3a)') 'export cif, ', substance(1:i), '.cif, spcgr:original'
write(IRE_MA,'( a)') 'exit  ! back to REFINE'
close(IRE_MA)
close(IRE_NP)
!
call write_prepare_hkl_full(substance, ilist, hkl_max, rlambda, slambda)
!
!  Interpret "aspher.json"
!
if(c_form=='discamb') then
   aspher_file = '../BUILD_TSC/aspher.json'  ! Use the special fixed name for aspher.json
   inquire(file=aspher_file, exist = lexist)
   cond_exist: if(lexist) then
      open(unit=IRE_MA, file=aspher_file, status='old')
      read(IRE_MA, '(a)', iostat=ios) line
      if(is_iostat_end(ios)) exit cond_exist
      if(ios/= 0           ) exit cond_exist
      if(line=='}'         ) exit cond_exist
      length=len_trim(line)
      call rem_leading_bl(line, length)
      loop_aspher: do
         i =index(line, 'structure')
         if(i > 0) then
            i = index(line, ':')
            j = index(line(i+1:len_trim(line)), '"')
            k = index(line, '"', .true.)
            line1 = line(i+j+1:k-1)
            line2 = 'cp ../ESSENTIAL_INPUT/' // line1(1:len_trim(line1)) // ' ../BUILD_TSC/'
            call execute_command_line(line2, wait=.false.)
            exit cond_exist
         endif
         read(IRE_MA, '(a)', iostat=ios) line
         if(is_iostat_end(ios)) exit cond_exist
         if(ios/= 0           ) exit cond_exist
         if(line=='}'         ) exit cond_exist
         length=len_trim(line)
         call rem_leading_bl(line, length)
      enddo loop_aspher
   endif cond_exist
   close(IRE_MA)
endif
!
write(line,'(3a)') 'cp ', substance(1:len_trim(substance)),'.cell CELL'
call execute_command_line(line, wait=.false.)
!
end subroutine write_refine_single_part4
!
!
!*****7**************************************************************** 
!
subroutine single_main_part4(IWR, l_refine, substance, c_form, rlambda, slambda, P_exti)
!-
!  Write the main discus macro Single crystal refinement 
!  DIFFEV and REFINE
!+
implicit none
!
integer           , intent(in) :: IWR
logical           , intent(in) :: l_refine     ! True == REFINE false==DIFFEV
character(len=*)  , intent(in) :: substance
character(len=*)  , intent(in) :: c_form
real(kind=PREC_DP), intent(in) :: rlambda
character(len=*)  , intent(in) :: slambda
real(kind=PREC_DP), intent(in) :: P_exti
!
character(len=PREC_STRING) :: line
integer :: i, j
!
j = 0
i = len_trim(substance)
write(IWR,'(a )') 'save'
write(IWR,'(3a)') '  outfile internal.',substance(1:i),'.cell'
write(IWR,'(a )') '  write all'
write(IWR,'(a )') '  run'
write(IWR,'(a )') 'exit'
write(IWR,'(a )') 'read'
write(IWR,'(3a)') '  cell internal.',substance(1:i),'.cell'
write(IWR,'(a )') 'variable integer, points_abs'
write(IWR,'(a )') 'variable integer, points_ord'
write(IWR,'(a )') 'variable integer, points_top'
write(IWR,'(a )') 'if(abs(F_XSTP)>0.0) then'
write(IWR,'(a )') '  points_abs = (F_XMAX-F_XMIN) / F_XSTP + 1'
write(IWR,'(a )') 'else'
write(IWR,'(a )') '  points_abs = 1'
write(IWR,'(a )') 'endif'
write(IWR,'(a )') 'if(abs(F_YSTP)>0.0) then'
write(IWR,'(a )') '  points_ord = (F_YMAX-F_YMIN) / F_YSTP + 1'
write(IWR,'(a )') 'else'
write(IWR,'(a )') '  points_ord = 1'
write(IWR,'(a )') 'endif'
write(IWR,'(a )') 'if(abs(F_ZSTP)>0.0) then'
write(IWR,'(a )') '  points_top = (F_ZMAX-F_ZMIN) / F_ZSTP + 1'
write(IWR,'(a )') 'else'
write(IWR,'(a )') '  points_top = 1'
write(IWR,'(a )') 'endif'
write(IWR,'(a )') 'fourier'
if(lexist_hkl_file .or. lexist_fcf_file) then     ! Shelx HKL file exists
   j= len_trim(hkl_file) - 4
   if(c_form == 'discamb') then
      write(IWR,'(3a)') '  xray table:discamb, file:', hkl_file(1:j), '.tsc'
      if(lexist_fcf_file) then   
         write(line, '(5a)') 'cp ', substance(1:i), '.ins', substance(1:i), '_merged.ins'
         call execute_command_line(line, wait=.false.)
      endif
   else
      write(IWR,'(2a)') '  xray table:', c_form(1:len_trim(c_form))
   endif
else
   write(IWR,'(2a)') '  xray table:', c_form(1:len_trim(c_form))
endif
write(IWR,'(a )') '  temp use'
write(IWR,'(a )') '  disp off'
if(slambda == ' ') then
   write(IWR,'(a, f7.5 )') '  wvle ', rlambda
else
   write(IWR,'(2a      )') '  wvle ', slambda(1:len_trim(slambda))
endif
write(IWR,'(a )') '  set aver, 0'
if(P_exti> 0.00001) then
   write(IWR,'(a )') '  set exti:P_exti'
endif
if(lexist_hkl_file .or. lexist_fcf_file) then                      ! Single HKL refinement
   write(IWR,'(a )') '  set technique:turbo'
   if(l_refine) then                          ! REFINE
      write(IWR,'(3a)') '  hkl in:',hkl_file(1:len_trim(hkl_file)), ', out:calc.hkl, scale:P_scale, style:hklf4'
   else
      write(IWR,'( a)') '  hkl in:"%c/DATA/%c", DATADIR, DATAFILE, out:"%c/INDI/indi.%4D", INDIDIR, REF_KID, scale:P_scale, style:hklf4'
   endif
   write(IWR,'(a )') 'exit'
else
   write(IWR,'(a )') '  set technique:nufft'
   write(IWR,'(a )') '  set scale:P_scale'
   write(IWR,'(a )') '  ll F_XMIN, F_YMIN, F_ZMIN'
   write(IWR,'(a )') '  lr F_XMAX, F_YMIN, F_ZMIN'
   write(IWR,'(a )') '  ul F_XMIN, F_YMAX, F_ZMIN'
   write(IWR,'(a )') '  tl F_XMIN, F_YMIN, F_ZMAX'
   write(IWR,'(a )') '  na points_abs'
   write(IWR,'(a )') '  no points_ord'
   write(IWR,'(a )') '  nt points_top'
   write(IWR,'(a )') '  !set filter:lanczos, damp:0.5, width:7, scale:1.0'
   write(IWR,'(a )') '  !set symmetry:apply'
   write(IWR,'(a )') '  run ! sigabs:[0.001, 1.0,0.0,0.0, sigord:[0.001, 0.0,1.0,0.0], sigtop:[0.001, 0.0,0.0,1.0]'
   write(IWR,'(a )') 'exit'
   write(IWR,'(a )') 'output'
   if(l_refine) then
      write(IWR,'(a )') '  outf calc.h5'
   else
      write(IWR,'(a )') '  outf "%c/INDI/indi.%4D", INDIDIR, REF_KID'
   endif
   write(IWR,'(a )') '  form hdf5'
   write(IWR,'(a )') '  value inte'
   write(IWR,'(a )') '  run'
   write(IWR,'(a )') 'exit'
endif
!
end subroutine single_main_part4
!
!*******************************************************************************
!
subroutine single_k_inter_macro(IWR, l_refine, substance)
!-
!  Write the k_fobs_fcalc.mac' or 'k_inter.mac'
!  REFINE version
!+
!
implicit none
!
integer         , intent(in) :: IWR
logical         , intent(in) :: l_refine    ! REFINE==true; DIFFEV==false
character(len=*), intent(in) :: substance   ! Compound name
!
integer :: i
!
i= len_trim(substance)
if(lexist_hkl_file .or. lexist_fcf_file) then                      ! Single HKL refinement
   if(l_refine) then
      open(IWR, file='k_fobs_fcalc.mac', status='unknown')
      write(IWR,'(a )') 'reset'
      if(lexist_hkl_file) then
         write(IWR,'(3a)') 'load csv, ',substance(1:len_trim(substance)), '.hkl, colx:4, coly:5, separator:[4,4,4,8,8], skip:0'
      elseif(lexist_fcf_file) then                      ! Single HKL refinement
         write(IWR,'(3a)') 'load csv, ',substance(1:len_trim(substance)), '_merged.hkl, colx:4, coly:5, separator:[4,4,4,8,8], skip:0'
      endif
      write(IWR,'(a )') 'load csv, calc.hkl, colx:5, coly:4, separator:[4,4,4,8,8], skip:0'
   else
      open(IWR, file='ksingle.mac', status='unknown')
      write(IWR,'(a )') 'reset'
      write(IWR,'(3a)') 'load csv, DATA/',substance(1:i), '_merged.hkl, colx:4, coly:5, separator:[4,4,4,8,8], skip:0'
      write(IWR,'(a )') 'load csv, "FINAL/final.%4D", $1, colx:5, coly:4, separator:[4,4,4,8,8], skip:0'
   endif
   write(IWR,'(a )') 'ccal mul, wy, 1, 0.0'
   write(IWR,'(a )') 'kcal add, 1, 2'
   write(IWR,'(a )') 'kfra 1, 3'
   write(IWR,'(a )') 'scale 0, max(xmax[3], ymax[3]), 0.0, max(xmax[3], ymax[3])'
   write(IWR,'(a )') 'aver 1'
   write(IWR,'(a )') 'mark'
   write(IWR,'(a )') 'ltyp 3, 0'
   write(IWR,'(a )') 'mtyp 3, 3'
   write(IWR,'(a )') 'fnam off'
   write(IWR,'(a )') 'fset 2'
   write(IWR,'(a )') 'grid on'
   write(IWR,'(a )') 'achx Iobs'
   write(IWR,'(a )') 'achy Icalc'
   write(IWR,'(a )') 'mcol 3, black'
   if(l_refine) then
      if(lexist_hkl_file) then
         write(IWR,'(3a)') 'load csv, ',substance(1:len_trim(substance)), '.hkl, colx:5, coly:4, separator:[4,4,4,8,8], skip:0'
      elseif(lexist_fcf_file) then                      ! Single HKL refinement
         write(IWR,'(3a)') 'load csv, ',substance(1:len_trim(substance)), '_merged.hkl, colx:5, coly:4, separator:[4,4,4,8,8], skip:0'
      endif
      write(IWR,'(a )') 'load csv, calc.hkl, colx:5, coly:4, separator:[4,4,4,8,8], skip:0'
   else
      write(IWR,'(3a)') 'load csv, DATA/',substance(1:i), '_merged.hkl, colx:5, coly:4, separator:[4,4,4,8,8], skip:0'
      write(IWR,'(a )') 'load csv, "FINAL/final.%4D", $1, colx:5, coly:4, separator:[4,4,4,8,8], skip:0'
   endif
   write(IWR,'(a )') 'do LOOP = 1, np[4]'
   write(IWR,'(a )') '  dy[4,LOOP] = x[4,LOOP]'
   write(IWR,'(a )') 'enddo'
   write(IWR,'(a )') 'rval 4,5, dat'
   write(IWR,'(2a)') 'tit1 Substance: ', substance(1:i)
   write(IWR,'(a )') 'tit2 "wR-value %8.5f",res[2]'
   write(IWR,'(a )') 'plot'
   if(l_refine) then
      write(IWR,'(a )') 'exit '
   endif
else
   if(l_refine) then
      open(IWR, file='k_inter.mac', status='unknown')
      write(IWR,'(a )') 'reset'
      write(IWR,'(3a)') 'load h5, DATA/',substance(1:len_trim(substance)),'.h5, layer:middle'
      write(IWR,'(a )') 'load h5, calc.h5, layer:middle'
   else
      open(IWR, file='ksingle.mac', status='unknown')
      write(IWR,'(a )') 'reset'
      write(IWR,'(3a)') 'load h5, DATA/',substance(1:len_trim(substance)),'.h5, layer:middle'
      write(IWR,'(a )') 'load h5, "FINAL/final.%4D", $1, layer:middle'
   endif
   write(IWR,'(a )') 'rval 1, 2, one'
   write(IWR,'(a )') 'r[12] = res[1]'
   write(IWR,'(a )') 'kcal sub, 1, 2'
   write(IWR,'(a )') 'nfra 3'
   write(IWR,'(a )') 'kfra 1, 1'
   write(IWR,'(a )') 'kfra 2, 2'
   write(IWR,'(a )') 'kfra 3, 3'
   write(IWR,'(a )') 'sfra 1, 0.00, 0.45, 0.55, 1.00'
   write(IWR,'(a )') 'sfra 2, 0.45, 0.45, 1.00, 1.00'
   write(IWR,'(a )') 'sfra 3, 0.00, 0.00, 0.55, 0.55'
   write(IWR,'(a )') '#'
   write(IWR,'(a )') 'afra 1'
   write(IWR,'(a )') 'aver data:1'
   write(IWR,'(a )') 'angle data:1'
   write(IWR,'(a )') 'scale'
   write(IWR,'(a )') 'mark'
   write(IWR,'(a )') 'fnam off'
   write(IWR,'(2a)') 'tit1 Substance: ', substance(1:len_trim(substance))
   write(IWR,'(a )') 'tit2 Observed intensity'
   write(IWR,'(a )') 'achx [H,0,0]'
   write(IWR,'(a )') 'achy [0,K,0]'
   write(IWR,'(a )') 'achz Intensity'
   write(IWR,'(a )') 'hart 1, 2'
   write(IWR,'(a )') 'cmap fire'
   write(IWR,'(a )') 'hlin 1, 0.1, 0.1, 100, %'
   write(IWR,'(a )') 'fset 3'
   write(IWR,'(a )') '#'
   write(IWR,'(a )') 'afra 2'
   write(IWR,'(a )') 'aver data:1'
   write(IWR,'(a )') 'angle data:1'
   write(IWR,'(a )') 'scale'
   write(IWR,'(a )') 'mark'
   write(IWR,'(a )') 'fnam off'
   write(IWR,'(a )') 'tit1 Calculated intensity'
   write(IWR,'(2a)') 'tit2 "R-value: %10.5f", r[12]'
   write(IWR,'(a )') 'achx [H,0,0]'
   write(IWR,'(a )') 'achy [0,K,0]'
   write(IWR,'(a )') 'achz Intensity'
   write(IWR,'(a )') 'hart 2, 2'
   write(IWR,'(a )') 'cmap fire'
   write(IWR,'(a )') 'hlin 1, 0.1, 0.1, 100, %'
   write(IWR,'(a )') 'fset -3'
   write(IWR,'(a )') '#'
   write(IWR,'(a )') 'afra 3'
   write(IWR,'(a )') 'aver data:1'
   write(IWR,'(a )') 'angle data:1'
   write(IWR,'(a )') 'scale'
   write(IWR,'(a )') 'mark'
   write(IWR,'(a )') 'fnam off'
   write(IWR,'(2a)') 'tit1  '
   write(IWR,'(a )') 'tit2 Difference Obs - calc'
   write(IWR,'(a )') 'achz Intensity'
   write(IWR,'(a )') 'achx [H,0,0]'
   write(IWR,'(a )') 'achy [0,K,0]'
   write(IWR,'(a )') 'hart 3, 2'
   write(IWR,'(a )') 'cmap pdf'
   write(IWR,'(a )') 'if(zmin[3]<0.0 .and.zmax[3]>0) then'
   write(IWR,'(a )') '  r[10] = -max( 0.50*abs(zmin[3]), 0.50*zmax[3])'
   write(IWR,'(a )') 'elseif(zmin[3]>=0.0 .and.zmax[3]>0) then'
   write(IWR,'(a )') '  r[10] = -   (                    0.50*zmax[3])'
   write(IWR,'(a )') 'elseif(zmin[3]< 0.0 .and.zmax[3]<=0) then'
   write(IWR,'(a )') '  r[10] = -   ( 0.50*abs(zmin[3])              )'
   write(IWR,'(a )') 'endif'
   write(IWR,'(a )') 'r[11] = -0.02*r[10]'
   write(IWR,'(a )') 'hlin 1, r[10], r[11], 101'
   write(IWR,'(a )') 'fset 3'
   write(IWR,'(a )') '#'
   write(IWR,'(a )') 'plot'
   if(l_refine) then
      write(IWR,'(a )') 'exit'
   endif
endif
close(IWR)
!
end subroutine single_k_inter_macro
!
!*****7**************************************************************** 
!
subroutine write_refine_powder_part1(substance, ilist, P_scale, &
           spcgr_syst, lattice_para, c_style)
!
!  Write main discus macro; powder diffraction PART 1
!
!
use precision_mod
!
implicit none
!
character(len=*)   , intent(in) :: substance
integer            , intent(in) :: ilist
real(kind=PREC_DP) , intent(in) :: P_scale
integer            , intent(in) :: spcgr_syst
real(kind=PREC_DP) , dimension(6), intent(in) :: lattice_para
character(len=PREC_STRING), intent(in) ::    c_style
!
character(len=PREC_STRING) :: refine_file
character(len=PREC_STRING) :: line
integer :: i
!
line = 'mkdir -p DATA'
call execute_command_line(line, wait=.false.)
!
i= len_trim(substance)
!
if(c_style=='powder') then
   line = 'mkdir -p POWDER'
   call execute_command_line(line, wait=.false.)
   discus_file = substance(1:i)//'_main_powder.mac'
   refine_file = 'refine_'//substance(1:i)//'_powder.mac'
!
   if(ilist == 0) then
      hkl_file = substance(1:len_trim(substance))//'.tth'
   else
      hkl_file = substance(1:len_trim(substance))//'.Q'
   endif
elseif(c_style=='pdf') then
   line = 'mkdir -p GRCALC'
   call execute_command_line(line, wait=.false.)
   discus_file = substance(1:i)//'_main_pdf.mac'
   refine_file = 'refine_'//substance(1:i)//'_pdf.mac'
!
   hkl_file = substance(1:len_trim(substance))//'.grobs'
endif
open(unit=IDI_MA, file=discus_file, status='unknown')
open(unit=IRE_MA, file=refine_file, status='unknown')
open(unit=IRE_NP, file='newpara.mac', status='unknown')
!
i= len_trim(substance)
write(IDI_MA,'(a )') 'branch discus'
write(IDI_MA,'(a )') '#'
write(IDI_MA,'(a )') 'read'
write(IDI_MA,'(3a)') '  stru CELL/', substance(1:i), '.cell'
!
write(IRE_MA,'(a )') 'refine'
write(IRE_MA,'(a )') 'rese'
write(IRE_MA,'(2a)') 'data xy, DATA/',hkl_file(1:len_trim(hkl_file))
write(IRE_MA,'(a        )') '#'
write(IRE_MA,'(a )') '@newpara.mac'
write(IRE_MA,'(a        )') '#'
!
write(IRE_NP,'(a,f9.4,a )') 'newpara P_scale, value:', P_scale , ', points:3, shift:0.001, status:free'
write(IRE_NP,'(a,f9.4,a )') 'newpara P_zero , value:', 0.00000 , ', points:3, shift:0.001, status:free'
write(IRE_NP,'(a )') '#'
write(IRE_NP,'(a,f9.4,a )') 'newpara P_eta  , value:', 0.50000 , ', points:3, shift:0.001, status:fixed'
write(IRE_NP,'(a,f9.4,a )') 'newpara P_eta_l, value:', 0.00000 , ', points:3, shift:0.001, status:fixed'
write(IRE_NP,'(a,f9.4,a )') 'newpara P_eta_q, value:', 0.00000 , ', points:3, shift:0.001, status:fixed'
write(IRE_NP,'(a,f9.4,a )') 'newpara P_u    , value:', 0.00000 , ', points:3, shift:0.001, status:fixed'
write(IRE_NP,'(a,f9.4,a )') 'newpara P_v    , value:', 0.00000 , ', points:3, shift:0.001, status:fixed'
write(IRE_NP,'(a,f9.4,a )') 'newpara P_w    , value:', 0.00700 , ', points:3, shift:0.001, status:fixed'
write(IRE_NP,'(a )') '#'
write(IRE_NP,'(a,f9.4,a )') 'newpara P_as1_c, value:', 0.00000 , ', points:3, shift:0.001, status:fixed'
write(IRE_NP,'(a,f9.4,a )') 'newpara P_as2_c, value:', 0.00000 , ', points:3, shift:0.001, status:fixed'
write(IRE_NP,'(a,f9.4,a )') 'newpara P_as1_i, value:', 0.00000 , ', points:3, shift:0.001, status:fixed'
write(IRE_NP,'(a,f9.4,a )') 'newpara P_as2_i, value:', 0.00000 , ', points:3, shift:0.001, status:fixed'
write(IRE_NP,'(a,f9.4,a )') 'newpara P_as1_l, value:', 0.00000 , ', points:3, shift:0.001, status:fixed'
write(IRE_NP,'(a,f9.4,a )') 'newpara P_as2_l, value:', 0.00000 , ', points:3, shift:0.001, status:fixed'
write(IRE_NP,'(a,f9.4,a )') 'newpara P_as1_q, value:', 0.00000 , ', points:3, shift:0.001, status:fixed'
write(IRE_NP,'(a,f9.4,a )') 'newpara P_as2_q, value:', 0.00000 , ', points:3, shift:0.001, status:fixed'
write(IRE_NP,'(a )') '#'
!
if(spcgr_syst==1) then         ! Triclinic
   write(IRE_NP,'(a,f9.4,a )') 'newpara P_alat , value:', lattice_para(1) , ', points:3, shift:0.001, status:free'
   write(IRE_NP,'(a,f9.4,a )') 'newpara P_blat , value:', lattice_para(2) , ', points:3, shift:0.001, status:free'
   write(IRE_NP,'(a,f9.4,a )') 'newpara P_clat , value:', lattice_para(3) , ', points:3, shift:0.001, status:free'
   write(IRE_NP,'(a,f9.4,a )') 'newpara P_alpha, value:', lattice_para(4) , ', points:3, shift:0.001, status:free'
   write(IRE_NP,'(a,f9.4,a )') 'newpara P_beta , value:', lattice_para(5) , ', points:3, shift:0.001, status:free'
   write(IRE_NP,'(a,f9.4,a )') 'newpara P_gamma, value:', lattice_para(6) , ', points:3, shift:0.001, status:free'
elseif(spcgr_syst==2) then         ! Monoclinic-B
   write(IRE_NP,'(a,f9.4,a )') 'newpara P_alat , value:', lattice_para(1) , ', points:3, shift:0.001, status:free'
   write(IRE_NP,'(a,f9.4,a )') 'newpara P_blat , value:', lattice_para(2) , ', points:3, shift:0.001, status:free'
   write(IRE_NP,'(a,f9.4,a )') 'newpara P_clat , value:', lattice_para(3) , ', points:3, shift:0.001, status:free'
   write(IRE_NP,'(a,f9.4,a )') 'newpara P_beta , value:', lattice_para(5) , ', points:3, shift:0.001, status:free'
elseif(spcgr_syst==3) then         ! Monoclinic-C
   write(IRE_NP,'(a,f9.4,a )') 'newpara P_alat , value:', lattice_para(1) , ', points:3, shift:0.001, status:free'
   write(IRE_NP,'(a,f9.4,a )') 'newpara P_blat , value:', lattice_para(2) , ', points:3, shift:0.001, status:free'
   write(IRE_NP,'(a,f9.4,a )') 'newpara P_clat , value:', lattice_para(3) , ', points:3, shift:0.001, status:free'
   write(IRE_NP,'(a,f9.4,a )') 'newpara P_gamma, value:', lattice_para(5) , ', points:3, shift:0.001, status:free'
elseif(spcgr_syst==4) then         ! Orthorhombic
   write(IRE_NP,'(a,f9.4,a )') 'newpara P_alat , value:', lattice_para(1) , ', points:3, shift:0.001, status:free'
   write(IRE_NP,'(a,f9.4,a )') 'newpara P_blat , value:', lattice_para(2) , ', points:3, shift:0.001, status:free'
   write(IRE_NP,'(a,f9.4,a )') 'newpara P_clat , value:', lattice_para(3) , ', points:3, shift:0.001, status:free'
elseif(spcgr_syst==5) then         ! tetragonal  
   write(IRE_NP,'(a,f9.4,a )') 'newpara P_alat , value:', lattice_para(1) , ', points:3, shift:0.001, status:free'
   write(IRE_NP,'(a,f9.4,a )') 'newpara P_clat , value:', lattice_para(3) , ', points:3, shift:0.001, status:free'
elseif(spcgr_syst==6) then         ! trigonal  
   write(IRE_NP,'(a,f9.4,a )') 'newpara P_alat , value:', lattice_para(1) , ', points:3, shift:0.001, status:free'
   write(IRE_NP,'(a,f9.4,a )') 'newpara P_clat , value:', lattice_para(3) , ', points:3, shift:0.001, status:free'
elseif(spcgr_syst==7) then         ! trigonal  
   write(IRE_NP,'(a,f9.4,a )') 'newpara P_alat , value:', lattice_para(1) , ', points:3, shift:0.001, status:free'
   write(IRE_NP,'(a,f9.4,a )') 'newpara P_alpha, value:', lattice_para(4) , ', points:3, shift:0.001, status:free'
elseif(spcgr_syst==8) then         ! hexagonal  
   write(IRE_NP,'(a,f9.4,a )') 'newpara P_alat , value:', lattice_para(1) , ', points:3, shift:0.001, status:free'
   write(IRE_NP,'(a,f9.4,a )') 'newpara P_clat , value:', lattice_para(3) , ', points:3, shift:0.001, status:free'
elseif(spcgr_syst==9) then         ! cubic      
   write(IRE_NP,'(a,f9.4,a )') 'newpara P_alat , value:', lattice_para(1) , ', points:3, shift:0.001, status:free'
endif
!
end subroutine write_refine_powder_part1
!
!*****7**************************************************************** 
!
subroutine write_refine_powder_part4(substance,                      &
        ilist, c_form,  lambda, slambda, P_exti, spcgr_syst, lattice_para, hkl_max)
!
!  Write main discus macro; powder diffraction  PART 4
!
use blanks_mod
use lib_conv_shelx_mod
use precision_mod
!
implicit none
!
character(len=*)                , intent(in) :: substance
integer                         , intent(in) :: ilist
character(len=*)                , intent(in) :: c_form
character(len=*)                , intent(in) ::  lambda 
character(len=*)                , intent(in) :: slambda 
real(kind=PREC_DP)              , intent(in) :: P_exti
integer            , intent(in) :: spcgr_syst
real(kind=PREC_DP) , dimension(6), intent(in) :: lattice_para
integer           , dimension(3), intent(inout) :: hkl_max
!
character(len=PREC_STRING) :: aspher_file
character(len=PREC_STRING) :: line
character(len=PREC_STRING) :: line1
character(len=PREC_STRING) :: line2
integer :: i,j, k
integer :: length
integer :: ios
logical :: lexist
real(kind=PREC_DP)  :: rlambda
!
read(lambda, *) rlambda
!
i= len_trim(substance)
j= len_trim(hkl_file) - 4
!
call write_powder_lattice(IDI_MA, spcgr_syst)
call write_macro_powder(c_form, hkl_file, rlambda, slambda)
!
write(IDI_MA,'(a )') 'save'
write(IDI_MA,'(3a)') '  outfile internal.',substance(1:i), '.cell'
write(IDI_MA,'(a )') '  write all'
write(IDI_MA,'(a )') '  run'
write(IDI_MA,'(a )') 'exit'
write(IDI_MA,'(a )') 'read'
write(IDI_MA,'(3a)') '  cell internal.',substance(1:i), '.cell'
write(IDI_MA,'(a )') 'variable real, rlambda'
write(IDI_MA,'(a )') 'rlambda = 1.5409'
write(IDI_MA,'(a )') '@powder.mac'
!
write(IDI_MA,'(a )') 'output'
write(IDI_MA,'(3a)') '  outf POWDER/', substance(1:i), '.raw'
write(IDI_MA,'(a )') '  value inte'
write(IDI_MA,'(a )') '  form powder, tth, F_XMIN, F_XMAX, F_XSTP'
write(IDI_MA,'(a )') '  run'
write(IDI_MA,'(a )') 'exit'
write(IDI_MA,'(a )') 'branch kuplot'
write(IDI_MA,'(a )') 'reset'
write(IDI_MA,'(3a)') 'load xy, POWDER/', substance(1:i), '.raw'
write(IDI_MA,'(3a)') 'load xy, DATA/', substance(1:i), '.tth'
write(IDI_MA,'(a )') 'fit 2'
write(IDI_MA,'(a )') '  func back, 1, 5'
write(IDI_MA,'(a )') '  para 1, 0, 1.0'
write(IDI_MA,'(a )') '  para 2, 1, ymin[2]'
write(IDI_MA,'(a )') '  para 3, 1, 0.0'
write(IDI_MA,'(a )') '  para 4, 1, 0.0'
write(IDI_MA,'(a )') '  para 5, 1, 0.0'
write(IDI_MA,'(a )') '  cycle 10'
write(IDI_MA,'(a )') '  run'
write(IDI_MA,'(a )') 'exit'
write(IDI_MA,'(a )') 'ksav n[1]-1'
write(IDI_MA,'(3a)') '  outf POWDER/', substance(1:i),'.tth'
write(IDI_MA,'(a )') '  form xy'
write(IDI_MA,'(a )') '  run'
write(IDI_MA,'(a )') 'n[1] = n[1]-1'
write(IDI_MA,'(a )') 'exit  ! Back to DISCUS'
write(IDI_MA,'(a )') 'exit  ! Back to REFINE'
write(IDI_MA,'(a )') 'finished'
close(IDI_MA)
!
i= len_trim(substance)
open(IDI_PO, file='k_inter.mac', status='unknown')
write(IDI_PO,'(a )') 'reset'
write(IDI_PO,'(3a)') 'load xy, DATA/', substance(1:i), '.tth'
write(IDI_PO,'(3a)') 'load xy, POWDER/', substance(1:i), '.tth'
write(IDI_PO,'(a )') 'kcal sub, 1, 2'
write(IDI_PO,'(a )') 'r[0] = min(ymin[1], ymin[2]) - ymax[3]'
write(IDI_PO,'(a )') 'ccal add, wy, 3, r[0]'
write(IDI_PO,'(a )') 'aver'
write(IDI_PO,'(a )') 'scale'
write(IDI_PO,'(a )') 'mark'
write(IDI_PO,'(a )') 'ltyp 1, 1'
write(IDI_PO,'(a )') 'ltyp 2, 1'
write(IDI_PO,'(a )') 'ltyp 3, 1'
write(IDI_PO,'(a )') 'lcol 1, blue'
write(IDI_PO,'(a )') 'lcol 2, red'
write(IDI_PO,'(a )') 'lcol 3, black'
write(IDI_PO,'(a )') 'mtyp 1, 0'
write(IDI_PO,'(a )') 'mtyp 2, 0'
write(IDI_PO,'(a )') 'mtyp 3, 0'
write(IDI_PO,'(a )') 'fnam off'
write(IDI_PO,'(a )') 'fset 2'
write(IDI_PO,'(a )') 'grid on'
write(IDI_PO,'(a,f9.7 )') 'achx 2Theta @ ', rlambda
write(IDI_PO,'(a )') 'achy Intensity'
write(IDI_PO,'(a )') 'rval 1,2, dat'
write(IDI_PO,'(2a)') 'tit1 Substance: ', substance(1:i)
write(IDI_PO,'(a )') 'tit2 "wR-value %8.5f",res[2]'
write(IDI_PO,'(a )') 'plot'
write(IDI_PO,'(a )') 'exit '
close(IDI_PO)
!
write(IRE_MA,'(a )') 'set cycle,   5'
write(IRE_MA,'(a )') 'set conver, status:on, dchi:0.050, chisq:1.10, pshift:2.0, conf:1.0, lambda:65000.'
write(IRE_MA,'(a )') 'set relax, start:0.25'
write(IRE_MA,'(2a)') '@', discus_file(1:len_trim(discus_file))
write(IRE_MA,'( a)') 'branch kuplot'
write(IRE_MA,'( a)') '  @k_inter.mac'
write(IRE_MA,'(3a)') 'run ', discus_file(1:len_trim(discus_file)), ', plot:k_inter.mac'
write(IRE_MA,'(a )') '#'
write(IRE_MA,'(a )') '@final_cell_cif.mac'
write(IRE_MA,'(a )') '#'
write(IRE_MA,'(a )') 'exit  ! Back to SUITE'
close(IRE_MA)
!
open(IRE_MA, file='final_cell_cif.mac', status='unknown')
write(IRE_MA,'( a)') 'branch discus'
write(IRE_MA,'( a)') 'read'
write(IRE_MA,'(3a)') '  stru internal.' , substance(1:i), '.cell'
write(IRE_MA,'( a)') 'save'
write(IRE_MA,'(3a)') '  outf ', substance(1:i), '.cell'
write(IRE_MA,'( a)') '  write all'
write(IRE_MA,'( a)') '  run'
write(IRE_MA,'( a)') 'exit'
write(IRE_MA,'(3a)') 'export cif, ', substance(1:i), '.cif, spcgr:original'
write(IRE_MA,'( a)') 'exit  ! back to REFINE'
close(IRE_MA)
close(IRE_NP)
!
!  Interpret "aspher.json"
!
if(c_form=='discamb') then
   aspher_file = '../BUILD_TSC/aspher.json'  ! Use the special fixed name for aspher.json
   inquire(file=aspher_file, exist = lexist)
   cond_exist: if(lexist) then
      open(unit=IRE_MA, file=aspher_file, status='old')
      read(IRE_MA, '(a)', iostat=ios) line
      if(is_iostat_end(ios)) exit cond_exist
      if(ios/= 0           ) exit cond_exist
      if(line=='}'         ) exit cond_exist
      length=len_trim(line)
      call rem_leading_bl(line, length)
      loop_aspher: do
         i =index(line, 'structure')
         if(i > 0) then
            i = index(line, ':')
            j = index(line(i+1:len_trim(line)), '"')
            k = index(line, '"', .true.)
            line1 = line(i+j+1:k-1)
            line2 = 'cp ../ESSENTIAL_INPUT/' // line1(1:len_trim(line1)) // ' ../BUILD_TSC/'
            call execute_command_line(line2, wait=.false.)
            exit cond_exist
         endif
         read(IRE_MA, '(a)', iostat=ios) line
         if(is_iostat_end(ios)) exit cond_exist
         if(ios/= 0           ) exit cond_exist
         if(line=='}'         ) exit cond_exist
         length=len_trim(line)
         call rem_leading_bl(line, length)
      enddo loop_aspher
   endif cond_exist
   close(IRE_MA)
endif
!
write(line,'(3a)') 'cp ', substance(1:len_trim(substance)),'.cell CELL'
call execute_command_line(line, wait=.false.)
!
end subroutine write_refine_powder_part4
!
!*****7**************************************************************** 
!
subroutine write_refine_pdf_part4(substance,                      &
        ilist, c_form,  lambda, slambda, P_exti, spcgr_syst, lattice_para, hkl_max)
!
!  Write main discus macro; PDF PART 4
!
use blanks_mod
use lib_conv_shelx_mod
use precision_mod
!
implicit none
!
character(len=*)                , intent(in) :: substance
integer                         , intent(in) :: ilist
character(len=*)                , intent(in) :: c_form
character(len=*)                , intent(in) ::  lambda 
character(len=*)                , intent(in) :: slambda 
real(kind=PREC_DP)              , intent(in) :: P_exti
integer            , intent(in) :: spcgr_syst
real(kind=PREC_DP) , dimension(6), intent(in) :: lattice_para
integer           , dimension(3), intent(inout) :: hkl_max
!
character(len=PREC_STRING) :: line
integer :: i,j
real(kind=PREC_DP)  :: rlambda
!
read(lambda, *) rlambda
!
line = 'mkdir -p GRCALC'
call execute_command_line(line, wait=.false.)
!
i= len_trim(substance)
j= len_trim(hkl_file) - 4
!
call write_powder_lattice(IDI_MA, spcgr_syst)
call write_macro_debye(c_form, hkl_file, rlambda, slambda)
!
write(IDI_MA,'(a )') 'save'
write(IDI_MA,'(3a)') '  outfile internal.',substance(1:i), '.cell'
write(IDI_MA,'(a )') '  write all'
write(IDI_MA,'(a )') '  run'
write(IDI_MA,'(a )') 'exit'
write(IDI_MA,'(a )') 'variable integer, ncell_a'
write(IDI_MA,'(a )') 'variable integer, ncell_b'
write(IDI_MA,'(a )') 'variable integer, ncell_c'
write(IDI_MA,'(a )') 'ncell_a = int(P_diam_a/lat[1]) + 4'
write(IDI_MA,'(a )') 'ncell_b = int(P_diam_b/lat[2]) + 4'
write(IDI_MA,'(a )') 'ncell_c = int(P_diam_c/lat[3]) + 4'
write(IDI_MA,'(a )') 'read'
write(IDI_MA,'(3a)') '  cell internal.',substance(1:i), '.cell, ncell_a, ncell_b, ncell_c'
write(IDI_MA,'(a )') 'surface'
write(IDI_MA,'(a )') '  boundary ellipsoid, P_diam_a, P_diam_b, P_diam_c'
write(IDI_MA,'(a )') 'exit'
write(IDI_MA,'(a )') 'purge type_yes'
write(IDI_MA,'(a )') 'variable real, rlambda'
write(IDI_MA,'(a,f11.6 )') 'rlambda = ', rlambda
write(IDI_MA,'(a )') '@debye.mac'
!
write(IDI_MA,'(a )') 'output'
write(IDI_MA,'(3a)') '  outf GRCALC/', substance(1:i), '.grcalc'
write(IDI_MA,'(a )') '  value PDF'
write(IDI_MA,'(a )') '  form pdf, r, F_XMIN, F_XMAX, F_XSTP'
write(IDI_MA,'(a )') '  run'
write(IDI_MA,'(a )') 'exit'
write(IDI_MA,'(a )') 'branch kuplot'
write(IDI_MA,'(a )') 'reset'
write(IDI_MA,'(3a)') 'load xy, GRCALC/', substance(1:i), '.grcalc'
write(IDI_MA,'(a )') 'exit  ! Back to DISCUS'
write(IDI_MA,'(a )') 'exit  ! Back to REFINE'
write(IDI_MA,'(a )') 'finished'
close(IDI_MA)
!
i= len_trim(substance)
open(IDI_PO, file='k_inter.mac', status='unknown')
write(IDI_PO,'(a )') 'reset'
write(IDI_PO,'(3a)') 'load xy, DATA/', substance(1:i), '.grobs'
write(IDI_PO,'(3a)') 'load xy, GRCALC/', substance(1:i), '.grcalc'
write(IDI_PO,'(a )') 'kcal sub, 1, 2'
write(IDI_PO,'(a )') 'r[0] = min(ymin[1], ymin[2]) - ymax[3]'
write(IDI_PO,'(a )') 'ccal add, wy, 3, r[0]'
write(IDI_PO,'(a )') 'aver'
write(IDI_PO,'(a )') 'scale'
write(IDI_PO,'(a )') 'mark'
write(IDI_PO,'(a )') 'ltyp 1, 1'
write(IDI_PO,'(a )') 'ltyp 2, 1'
write(IDI_PO,'(a )') 'ltyp 3, 1'
write(IDI_PO,'(a )') 'lcol 1, blue'
write(IDI_PO,'(a )') 'lcol 2, red'
write(IDI_PO,'(a )') 'lcol 3, black'
write(IDI_PO,'(a )') 'mtyp 1, 0'
write(IDI_PO,'(a )') 'mtyp 2, 0'
write(IDI_PO,'(a )') 'mtyp 3, 0'
write(IDI_PO,'(a )') 'fnam off'
write(IDI_PO,'(a )') 'fset 2'
write(IDI_PO,'(a )') 'grid on'
write(IDI_PO,'(a )') 'achx Distance [\A]'
write(IDI_PO,'(a )') 'achy PDF'
write(IDI_PO,'(a )') 'rval 1,2, dat'
write(IDI_PO,'(2a)') 'tit1 Substance: ', substance(1:i)
write(IDI_PO,'(a )') 'tit2 "wR-value %8.5f",res[2]'
write(IDI_PO,'(a )') 'plot'
write(IDI_PO,'(a )') 'exit '
close(IDI_PO)
!
write(IRE_NP,'(a,f9.4,a )') 'newpara P_diam_a , value:30.0000, points:3, shift:0.050, status:free'
write(IRE_NP,'(a,f9.4,a )') 'newpara P_diam_b , value:40.0000, points:3, shift:0.050, status:free'
write(IRE_NP,'(a,f9.4,a )') 'newpara P_diam_c , value:50.0000, points:3, shift:0.050, status:free'
write(IRE_MA,'(a )') 'set cycle,   5'
write(IRE_MA,'(a )') 'set conver, status:on, dchi:0.050, chisq:1.10, pshift:2.0, conf:1.0, lambda:65000.'
write(IRE_MA,'(a )') 'set relax, start:0.25'
write(IRE_MA,'(2a)') '@', discus_file(1:len_trim(discus_file))
write(IRE_MA,'( a)') 'branch kuplot'
write(IRE_MA,'( a)') '  @k_inter.mac'
write(IRE_MA,'(3a)') 'run ', discus_file(1:len_trim(discus_file)), ', plot:k_inter.mac'
write(IRE_MA,'(a )') '#'
write(IRE_MA,'(a )') '@final_cell_cif.mac'
write(IRE_MA,'(a )') '#'
write(IRE_MA,'(a )') 'exit  ! Back to SUITE'
close(IRE_MA)
close(IRE_NP)
!
open(IRE_MA, file='final_cell_cif.mac', status='unknown')
write(IRE_MA,'( a)') 'branch discus'
write(IRE_MA,'( a)') 'read'
write(IRE_MA,'(3a)') '  stru internal.' , substance(1:i), '.cell'
write(IRE_MA,'( a)') 'save'
write(IRE_MA,'(3a)') '  outf ', substance(1:i), '.cell'
write(IRE_MA,'( a)') '  write all'
write(IRE_MA,'( a)') '  run'
write(IRE_MA,'( a)') 'exit'
write(IRE_MA,'(3a)') 'export cif, ', substance(1:i), '.cif, spcgr:original'
write(IRE_MA,'( a)') 'exit  ! back to REFINE'
close(IRE_MA)
!
write(line,'(3a)') 'cp ', substance(1:len_trim(substance)),'.cell CELL'
call execute_command_line(line, wait=.false.)
!
end subroutine write_refine_pdf_part4
!
!*****7**************************************************************** 
!
subroutine write_powder_lattice(IDI, spcgr_syst)
!-
! Write the lattice parameters into the discus_main file 
! used by refine and diffev
!+
!
implicit none
!
integer                         , intent(in) :: IDI
integer                         , intent(in) :: spcgr_syst
!
if(spcgr_syst==1) then         ! Triclinic
   write(IDI,'(a )') 'lat[1] = P_alat'
   write(IDI,'(a )') 'lat[2] = P_blat'
   write(IDI,'(a )') 'lat[3] = P_clat'
   write(IDI,'(a )') 'lat[4] = P_alpha'
   write(IDI,'(a )') 'lat[5] = P_beta'
   write(IDI,'(a )') 'lat[6] = P_gamma'
elseif(spcgr_syst==2) then         ! monoclinic B
   write(IDI,'(a )') 'lat[1] = P_alat'
   write(IDI,'(a )') 'lat[2] = P_blat'
   write(IDI,'(a )') 'lat[3] = P_clat'
   write(IDI,'(a )') 'lat[4] = 90.0000'
   write(IDI,'(a )') 'lat[5] = P_beta'
   write(IDI,'(a )') 'lat[6] = 90.0000'
elseif(spcgr_syst==3) then         ! monoclinic C
   write(IDI,'(a )') 'lat[1] = P_alat'
   write(IDI,'(a )') 'lat[2] = P_blat'
   write(IDI,'(a )') 'lat[3] = P_clat'
   write(IDI,'(a )') 'lat[4] = 90.0000'
   write(IDI,'(a )') 'lat[5] = 90.0000'
   write(IDI,'(a )') 'lat[6] = P_gamma'
elseif(spcgr_syst==4) then         ! orthorhombic
   write(IDI,'(a )') 'lat[1] = P_alat'
   write(IDI,'(a )') 'lat[2] = P_blat'
   write(IDI,'(a )') 'lat[3] = P_clat'
   write(IDI,'(a )') 'lat[4] = 90.0000'
   write(IDI,'(a )') 'lat[5] = 90.0000'
   write(IDI,'(a )') 'lat[6] = 90.0000'
elseif(spcgr_syst==5) then         ! tetragonal
   write(IDI,'(a )') 'lat[1] = P_alat'
   write(IDI,'(a )') 'lat[2] = P_alat'
   write(IDI,'(a )') 'lat[3] = P_clat'
   write(IDI,'(a )') 'lat[4] = 90.0000'
   write(IDI,'(a )') 'lat[5] = 90.0000'
   write(IDI,'(a )') 'lat[6] = 90.0000'
elseif(spcgr_syst==6 .or. spcgr_syst==8) then         ! trigonal, hexagonal
   write(IDI,'(a )') 'lat[1] = P_alat'
   write(IDI,'(a )') 'lat[2] = P_alat'
   write(IDI,'(a )') 'lat[3] = P_clat'
   write(IDI,'(a )') 'lat[4] = 90.0000'
   write(IDI,'(a )') 'lat[5] = 90.0000'
   write(IDI,'(a )') 'lat[6] =120.0000'
elseif(spcgr_syst==7) then         ! Rhombohedral
   write(IDI,'(a )') 'lat[1] = P_alat'
   write(IDI,'(a )') 'lat[2] = P_alat'
   write(IDI,'(a )') 'lat[3] = P_alat'
   write(IDI,'(a )') 'lat[4] = P_alpha'
   write(IDI,'(a )') 'lat[5] = P_alpha'
   write(IDI,'(a )') 'lat[6] = P_alpha'
elseif(spcgr_syst==9) then         ! Cubic
   write(IDI,'(a )') 'lat[1] = P_alat'
   write(IDI,'(a )') 'lat[2] = P_alat'
   write(IDI,'(a )') 'lat[3] = P_alat'
   write(IDI,'(a )') 'lat[4] = 90.0000'
   write(IDI,'(a )') 'lat[5] = 90.0000'
   write(IDI,'(a )') 'lat[6] = 90.0000'
endif
end subroutine write_powder_lattice
!
!*****7**************************************************************** 
!
subroutine write_diffev_generic(substance, compute, c_style, ilist )
!-
!  Writes the general part for a DIFFEV refinement
!+
!
use precision_mod
!
implicit none
!
character(len=*), intent(in) :: substance      ! The substance name
character(len=*), intent(in) :: compute        ! 'serial' or 'parallel'
character(len=*), intent(in) :: c_style        ! 'single', 'powder', 'pdf' 
integer         , intent(in) :: ilist          ! SHELX LIST number
!
character(len=PREC_STRING) :: line     ! Dummy line
character(len=PREC_STRING) :: outfile     ! output file name(s)
!
! Create required folders
!
line = 'mkdir -p CELL'
call execute_command_line(line, wait=.false.)
line = 'mkdir -p DATA'
call execute_command_line(line, wait=.false.)
line = 'mkdir -p DIFFEV'
call execute_command_line(line, wait=.false.)
line = 'mkdir -p FINAL'
call execute_command_line(line, wait=.false.)
line = 'mkdir -p INDI'
call execute_command_line(line, wait=.false.)
line = 'mkdir -p LOGFILES'
call execute_command_line(line, wait=.false.)
line = 'mkdir -p TEMP'
call execute_command_line(line, wait=.false.)
!
if(ilist == 0) then
   hkl_file = substance(1:len_trim(substance))//'.hkl'
   fcf_file = ' '
   inquire(file=hkl_file, exist=lexist_hkl_file)
else
   hkl_file = substance(1:len_trim(substance))//'.hkl'
   fcf_file = substance(1:len_trim(substance))//'.fcf'
   inquire(file=hkl_file, exist=lexist_hkl_file)
   inquire(file=fcf_file, exist=lexist_fcf_file)
   if(lexist_fcf_file) then
      hkl_file = substance(1:len_trim(substance))//'_merged.hkl'    ! hkl file will be merged
      lexist_hkl_file = .false.          ! FCF exists, do not use original HKL file
   endif
endif
!
outfile = 'diffev_main_' // substance(1:len_trim(substance)) // '.mac'
open(IDF_SI, file=outfile, status='unknown')
open(IDF_NP, file='newpara.mac', status='unknown')
write(IDF_SI, '(a)') 'variable integer, diff_counter'
write(IDF_SI, '(a)') 'diffev'
write(IDF_SI, '(a)') '  reset'
if(c_style=='single') then
   if(lexist_hkl_file) then       ! An hkl file exists
      write(IDF_SI, '(4a)') '  data hklf4, DATA/', substance(1:len_trim(substance)), '.hkl', &
                            ' ! Load data file, populates F_XMIN etc.'
   elseif(lexist_fcf_file) then       ! An fcf file exists
      write(IDF_SI, '(4a)') '  data hklf4, DATA/', substance(1:len_trim(substance)), '_merged.hkl', &
                            ' ! Load data file, populates F_XMIN etc.'
   else
      write(IDF_SI, '(4a)') '  data h5, DATA/', substance(1:len_trim(substance)), '.h5', &
                            ' ! Load data file, populates F_XMIN etc.'
   endif
elseif(c_style=='powder') then
   write(IDF_SI, '(4a)') '  data xy, DATA/', substance(1:len_trim(substance)), '.tth', &
                         ' ! Load data file, populates F_XMIN etc.'
elseif(c_style=='pdf') then
   write(IDF_SI, '(4a)') '  data xy, DATA/', substance(1:len_trim(substance)), '.grobs', &
                         ' ! Load data file, populates F_XMIN etc.'
endif
write(IDF_SI, '(a)') '  @cleanup.mac       ! Remove old files'
write(IDF_SI, '(a)') '  @setup_rep.mac     ! Define number of individual repetitions'
write(IDF_SI, '(a)') '  @global.mac        ! Define global file and directory names'
write(IDF_SI, '(a)') '  @diffev_setup.mac  ! Define refinement details'
write(IDF_SI, '(a)') '  init               ! Initialize all parameters'
write(IDF_SI, '(a)') '  do diff_counter = 1, 10  ! Perform 10 loops!'
write(IDF_SI, '(a)') '    echo "In loop %4d ", diff_counter     ! Keep user posted'
write(IDF_SI, '(3a)') '    run_mpi discus, dis.diffev.mac  , repeat:REF_NINDIV, compute:', &
                     compute(1:len_trim(compute)),', logfile:LOGFILES/d'
if(compute=='parallel') then
  write(IDF_SI, '( a)') '    run_mpi kuplot, kup.diffev.mac  , repeat:1, compute:serial, logfile:LOGFILES/k'
endif
write(IDF_SI, '(a)') '    compare'
write(IDF_SI, '(a)') '  enddo'
write(IDF_SI, '(a)') 'exit   ! Back to SUITE'
write(IDF_SI, '(a)') 'exit   ! Terminate SUITE'
!
close(IDF_SI)
!
outfile = 'diffev_continue_' // substance(1:len_trim(substance)) // '.mac'
open(IDF_SI, file=outfile, status='unknown')
write(IDF_SI, '(a)') 'variable integer, diff_counter'
write(IDF_SI, '(a)') 'diffev'
write(IDF_SI, '(a)') '  reset'
if(c_style=='single') then
   if(lexist_hkl_file) then       ! An hkl file exists
      write(IDF_SI, '(4a)') '  data hklf4, DATA/', substance(1:len_trim(substance)), '.hkl', &
                            ' !  Load data file, populates F_XMIN etc.'
   elseif(lexist_fcf_file) then       ! An fcf file exists
      write(IDF_SI, '(4a)') '  data hklf4, DATA/', substance(1:len_trim(substance)), '_merged.hkl', &
                            ' !  Load data file, populates F_XMIN etc.'
   else
      write(IDF_SI, '(4a)') '  data h5, DATA/', substance(1:len_trim(substance)), '.h5', &
                            ' !  Load data file, populates F_XMIN etc.'
   endif
elseif(c_style=='powder') then
   write(IDF_SI, '(4a)') '  data xy, DATA/', substance(1:len_trim(substance)), '.tth', &
                         ' !  Load data file, populates F_XMIN etc.'
elseif(c_style=='pdf') then
   write(IDF_SI, '(4a)') '  data xy, DATA/', substance(1:len_trim(substance)), '.grobs', &
                         ' ! Load data file, populates F_XMIN etc.'
endif
write(IDF_SI, '(a)') '  @setup_rep.mac     ! Define number of individual repetitions'
write(IDF_SI, '(a)') '  @diffev_setup.mac  ! Define refinement details'
write(IDF_SI, '(a)') '  do diff_counter = 1, 10  ! Perform 10 loops!'
write(IDF_SI, '(a)') '    echo "In loop %4d ", diff_counter     ! Keep user posted'
write(IDF_SI, '(3a)') '    run_mpi discus, dis.diffev.mac  , repeat:REF_NINDIV, compute:', &
                     compute(1:len_trim(compute)),', logfile:LOGFILES/d'
if(compute=='parallel') then
  write(IDF_SI, '( a)') '    run_mpi kuplot, kup.diffev.mac  , repeat:1, compute:serial, logfile:LOGFILES/k'
endif
write(IDF_SI, '(a)') '    compare'
write(IDF_SI, '(a)') '  enddo'
write(IDF_SI, '(a)') 'exit   ! Back to SUITE'
write(IDF_SI, '(a)') 'exit   ! Terminate SUITE'
!
close(IDF_SI)
!
! Create a generic cleanup macro
!
outfile = 'cleanup.mac'
open(IDF_SI, file=outfile, status='unknown')
!
write(IDF_SI, '(a)') 'system rm -f DIFFEV/*'
write(IDF_SI, '(a)') 'system rm -f FINAL/*'
write(IDF_SI, '(a)') 'system rm -f INDI/*'
write(IDF_SI, '(a)') 'system rm -f LOGFILES/*'
write(IDF_SI, '(a)') 'system rm -f PLOT/*'
write(IDF_SI, '(a)') 'system rm -f STRU/*'
write(IDF_SI, '(a)') 'system rm -f TEMP/*'
!
close(IDF_SI)
!
!
! Create a generic setup_rep macro
!
outfile = 'setup_rep.mac'
open(IDF_SI, file=outfile, status='unknown')
!
if(compute=='serial') then
   write(IDF_SI, '(a)') 'REF_NINDIV = 1      ! No individual repetitions'
elseif(compute=='parallel') then
   write(IDF_SI, '(a)') 'REF_NINDIV = 2      ! Some individual repetitions'
endif
!
close(IDF_SI)
!
!
! Create the generic DIFFEV/DISCUS interface
!
outfile = 'dis.diffev.mac'
open(IDF_SI, file=outfile, status='unknown')
!
write(IDF_SI, '(a)') 'set error,exit'
write(IDF_SI, '(a)') '#'
write(IDF_SI, '(a)') '#@ HEADER'
write(IDF_SI, '(a)') '#@ NAME         dis.diffev.mac'
write(IDF_SI, '(a)') '#@'
write(IDF_SI, '(a)') '#@ KEYWORD      refinement, generic, powder'
write(IDF_SI, '(a)') '#@'
write(IDF_SI, '(a)') '#@ DESCRIPTION  Generic interface between DIFFEV and DISCUS'
write(IDF_SI, '(a)') '#@ DESCRIPTION  global.mac sets general file names'
write(IDF_SI, '(a)') '#@ DESCRIPTION  The loop over all individual repetitions is performed'
write(IDF_SI, '(a)') '#@ DESCRIPTION  serially within this calculation'
write(IDF_SI, '(a)') '#@ DESCRIPTION'
write(IDF_SI, '(a)') '#@ DESCRIPTION  The main calculation is performed by the macro'
if(c_style=='single') then
   if(lexist_hkl_file) then 
      write(IDF_SI, '(a)') '#@ DESCRIPTION  *_main_hkl.mac'
   elseif(lexist_fcf_file) then 
      write(IDF_SI, '(a)') '#@ DESCRIPTION  *_main_merged.mac'
   else
      write(IDF_SI, '(a)') '#@ DESCRIPTION  *_main_single.mac'
   endif
elseif(c_style=='powder') then
   write(IDF_SI, '(a)') '#@ DESCRIPTION  *_main_powder.mac'
elseif(c_style=='pdf') then
   write(IDF_SI, '(a)') '#@ DESCRIPTION  *_main_pdf.mac'
endif
write(IDF_SI, '(a)') '#@ DESCRIPTION  Averaging and possible scaling/background is done in '
write(IDF_SI, '(a)') '#@ DESCRIPTION  kup.diffev.mac'
write(IDF_SI, '(a)') ''
write(IDF_SI, '(a)') '#@'
write(IDF_SI, '(a)') '#@ PARAMETER    $0, 1'
write(IDF_SI, '(a)') '#@'
write(IDF_SI, '(a)') '#@ USAGE        @dis.diffev.mac'
write(IDF_SI, '(a)') '#@'
write(IDF_SI, '(a)') '#@ END'
write(IDF_SI, '(a)') '#'
write(IDF_SI, '(a)') '############################################################'
write(IDF_SI, '(a)') '#'
write(IDF_SI, '(a)') 'variable integer, indiv  ! Dummy variable for the local repetitions'
write(IDF_SI, '(a)') '!'
write(IDF_SI, '(a)') 'branch kuplot  ! Step sideways into KUPLOT'
write(IDF_SI, '(a)') '   reset       ! Ensure there are no data sets'
write(IDF_SI, '(a)') 'exit           ! Back to DISCUS section'
write(IDF_SI, '(a)') '#'
write(IDF_SI, '(a)') '@global.mac    ! Definitions of directories etc'
if(compute=='serial') write(IDF_SI, '(a)') 'do indiv = 1, REF_NINDIV'
if(c_style=='single') then
   if(lexist_hkl_file) then 
      write(IDF_SI, '(3a)') '  @',substance(1:len_trim(substance)), '_main_hkl.mac'
   elseif(lexist_fcf_file) then 
      write(IDF_SI, '(3a)') '  @',substance(1:len_trim(substance)), '_main_merged.mac'
   else
      write(IDF_SI, '(3a)') '  @',substance(1:len_trim(substance)), '_main_single.mac'
   endif
elseif(c_style=='powder') then
   write(IDF_SI, '(3a)') '  @',substance(1:len_trim(substance)), '_main_powder.mac'
elseif(c_style=='pdf') then
   write(IDF_SI, '(3a)') '  @',substance(1:len_trim(substance)), '_main_pdf.mac'
endif
if(compute=='serial') then
   write(IDF_SI, '(a)') 'enddo'
   write(IDF_SI, '(a)') 'branch   kuplot    !Switch to KUPLOT'
   write(IDF_SI, '(a)') '   @kup.diffev.mac ., REF_KID  ! contains a final exit'
   write(IDF_SI, '(a)') '#'
endif
write(IDF_SI, '(a)') 'exit    ! Finish DISCUS, back to DIFFEV'
!
close(IDF_SI)
!
! start diffev_setup.mac
!
open(unit=IDF_MA, file='diffev_setup.mac', status='unknown')
!
write(IDF_MA, '(a)') '#'
write(IDF_MA, '(a)') '#  diffev_setup.mac'
write(IDF_MA, '(a)') '#'
write(IDF_MA, '(a)') '################################################################################'
write(IDF_MA, '(a)') '#'
write(IDF_MA, '(a)') '#  Define population size, parameters and log files'
write(IDF_MA, '(a)') '#'
write(IDF_MA, '(a)') '#'
write(IDF_MA, '(a)') 'pop_gen[1]  =   0     ! We want to start the refinement in generation zero'
write(IDF_MA, '(a)') '#'
write(IDF_MA, '(a)') 'pop_n[1]    =  12    ! Define the population size should be 10 times the number of parameters'
write(IDF_MA, '(a)') 'pop_c[1]    =  12    ! Number of children in each generation usually identical to parent'
write(IDF_MA, '(a)') '#'
write(IDF_MA, '(a)') '#  For each parameter to be refined we define:'
write(IDF_MA, '(a)') '#     name, minimum/maximum allowed values , minimum/maximum for start range'
write(IDF_MA, '(a)') '#'
write(IDF_MA, '(a)') '#'
write(IDF_MA, '(a)') '@newpara.mac'
write(IDF_MA, '(a)') '#'
write(IDF_NP, '(a)') 'newpara  P_scale  ,  range:[0.001,  20.000], value:[ 0.9000, 1.100], status:free'
!
if(c_style/='single') then
   write(IDF_NP, '(a)') 'newpara  P_zero   ,  range:[0.000,   0.000], value:[ 0.0000, 0.000], status:fixed'
endif
!
!
! Create a 'diffev_fix_free macro
!
!outfile = 'diffev_fix_free.mac'
!open(IDF_FX, file=outfile, status='unknown')
!write(IDF_FX, '(a)') '#  '
!write(IDF_FX, '(a)') 'refine none  '
!write(IDF_FX, '(a)') 'refine P_scale  '
!
! Start the 'diffev_best macro
!
outfile = 'diffev_best.mac'
open(IDF_BS, file=outfile, status='unknown')
outfile = 'diffev_best_set_values.mac'
open(IDF_BS2, file=outfile, status='unknown')
!
write(IDF_BS, '(a)') '#  '
write(IDF_BS, '(a)') '#@ HEADER'
write(IDF_BS, '(a)') '#@ NAME         diffev_best.mac'
write(IDF_BS, '(a)') '#@ '
write(IDF_BS, '(a)') '#@ KEYWORD      diffev, best member'
write(IDF_BS, '(a)') '#@ '
write(IDF_BS, '(a)') '#@ DESCRIPTION  This macro contains the parameters for the current best'
write(IDF_BS, '(a)') '#@ DESCRIPTION  member. If run, the best member will be recreated.'
write(IDF_BS, '(a)') '#@ DESCRIPTION  As the random state is explicitely contained as well, the'
write(IDF_BS, '(a)') '#@ DESCRIPTION  best member will be recreated exactly.'
write(IDF_BS, '(a)') '#@'
write(IDF_BS, '(a)') '#@ PARAMETER    $0, 0'
write(IDF_BS, '(a)') '#@'
write(IDF_BS, '(a)') '#@ USAGE        @diffev_best.mac'
write(IDF_BS, '(a)') '#@'
write(IDF_BS, '(a)') '#@ END'
write(IDF_BS, '(a)') '#'
write(IDF_BS, '(a)') 'diffev'
if(c_style=='single') then
   if(lexist_hkl_file) then       ! An hkl file exists
      write(IDF_BS, '(4a)') '  data hklf4, DATA/', substance(1:len_trim(substance)), '.hkl', &
                            ' ! Load data file, populates F_XMIN etc.'
   elseif(lexist_fcf_file) then       ! An fcf file exists
      write(IDF_BS, '(4a)') '  data hklf4, DATA/', substance(1:len_trim(substance)), '_merged.hkl', &
                            ' ! Load data file, populates F_XMIN etc.'
   else
      write(IDF_BS, '(4a)') '  data h5, DATA/', substance(1:len_trim(substance)), '.h5', &
                            ' ! Load data file, populates F_XMIN etc.'
   endif
elseif(c_style=='powder') then
   write(IDF_BS, '(4a)') '  data xy, DATA/', substance(1:len_trim(substance)), '.tth', &
                         ' ! Load data file, populates F_XMIN etc.'
elseif(c_style=='pdf') then
   write(IDF_BS, '(4a)') '  data xy, DATA/', substance(1:len_trim(substance)), '.grobs', &
                         ' ! Load data file, populates F_XMIN etc.'
endif
!
write(IDF_BS, '(a)') 'exit'
write(IDF_BS, '(a)') 'discus'
write(IDF_BS, '(a)') '  reset'
write(IDF_BS, '(a)') 'exit'
write(IDF_BS, '(a)') 'kuplot'
write(IDF_BS, '(a)') '  reset'
write(IDF_BS, '(a)') 'exit'
write(IDF_BS, '(a)') '#'
write(IDF_BS, '(a)') 'discus'
write(IDF_BS, '(a)') 'REF_GENERATION =            0'
write(IDF_BS, '(a)') 'REF_MEMBER     =           12'
write(IDF_BS, '(a)') 'REF_CHILDREN   =           12'
write(IDF_BS, '(a)') 'REF_NINDIV     =            1'
write(IDF_BS, '(a)') 'REF_KID        =         9999'
write(IDF_BS, '(a)') 'REF_INDIV      =            1'
write(IDF_BS, '(a)') 'variable real, Rvalue'
diffev_dimension = 0
!
write(IDF_BS2, '(a)') 'Rvalue           =   1.0000'
!
end subroutine write_diffev_generic
!
!*****7**************************************************************** 
!
subroutine write_diffev_global(substance, datafile, user_name, compute, c_style)
!-
!  Writes "global.mac"  
!+
!
use precision_mod
!
implicit none
!
character(len=*), intent(in) :: substance
character(len=*), intent(in) :: datafile
character(len=*), intent(in) :: user_name
character(len=*), intent(in) :: compute       ! Individuals 'serial' or 'parallel'
character(len=*), intent(in) :: c_style       ! 'powder' or 'pdf'
!
character(len=PREC_STRING) :: outfile     ! output file name(s)
!
! Create a generic global macro
!
outfile = 'global.mac'
open(IDF_PO, file=outfile, status='unknown')
!
write(IDF_PO, '(a)') '#'
write(IDF_PO, '(a)') '# global.mac'
write(IDF_PO, '(a)') '#'
!
write(IDF_PO, '(a)') '#  define file names'
write(IDF_PO, '(a)') ''
write(IDF_PO, '(a)') 'variable character, TMPDIR    ! temporary DISCUS files   ''internal''     or ''.'''
write(IDF_PO, '(a)') 'variable character, INDIDIR   ! temporary indi directory  ''kuplot''      or ''.'''
write(IDF_PO, '(a)') 'variable character, DATADIR   ! input data directory     ''/tmp/user_name/'' or. ''.'''
write(IDF_PO, '(a)') 'variable character, DATAFILE  ! Input Data File  '
write(IDF_PO, '(a)') 'variable character, SUBSTANCE ! Substance name  '
write(IDF_PO, '(a)') '#'
write(IDF_PO, '(a)') '# LOCAL Version, serial'
write(IDF_PO, '(a)') '#'
if(compute=='serial') then
   write(IDF_PO, '(a )') 'DATADIR   = "%c",''.'''
   write(IDF_PO, '(3a)') 'DATAFILE  = "%c",''',datafile(1:len_trim(datafile)),''''
!  if(lexist_hkl_file) then
      write(IDF_PO, '(a )') 'INDIDIR   = "%c",''.'''
!  else
!     write(IDF_PO, '(a )') 'INDIDIR   = "%c",''kuplot'''
!  endif
   write(IDF_PO, '(a )') 'TMPDIR    = "%c",''internal'''
   write(IDF_PO, '(3a)') 'SUBSTANCE = "%c",''', substance(1:len_trim(substance)),''''
   write(IDF_PO, '(a )') '#'
   write(IDF_PO, '(a )') '# LOCAL Version, parallel'
   write(IDF_PO, '(a )') '#'
   write(IDF_PO, '(a )') '#DATADIR   = "%c",''.'''
   write(IDF_PO, '(3a)') '#DATAFILE  = "%c",''',datafile(1:len_trim(datafile)),''''
   write(IDF_PO, '(a )') '#INDIDIR   = "%c",''.'''
   write(IDF_PO, '(a )') '#TMPDIR    = "%c",''internal'''
   write(IDF_PO, '(3a)') '#SUBSTANCE = "%c",''', substance(1:len_trim(substance)),''''
elseif(compute=='parallel') then
   write(IDF_PO, '(a )') '#DATADIR   = "%c",''.'''
   write(IDF_PO, '(3a)') '#DATAFILE  = "%c",''',datafile(1:len_trim(datafile)),''''
!   if(lexist_hkl_file) then
      write(IDF_PO, '(a )') 'INDIDIR   = "%c",''#.'''
!   else
!      write(IDF_PO, '(a )') 'INDIDIR   = "%c",''#kuplot'''
!   endif
   write(IDF_PO, '(a )') '#TMPDIR    = "%c",''internal'''
   write(IDF_PO, '(3a)') '#SUBSTANCE = "%c",''', substance(1:len_trim(substance)),''''
   write(IDF_PO, '(a )') '#'
   write(IDF_PO, '(a )') '# LOCAL Version, parallel'
   write(IDF_PO, '(a )') '#'
   write(IDF_PO, '(a )') 'DATADIR   = "%c",''.'''
   write(IDF_PO, '(3a)') 'DATAFILE  = "%c",''',datafile(1:len_trim(datafile)),''''
   write(IDF_PO, '(a )') 'INDIDIR   = "%c",''.'''
   write(IDF_PO, '(a )') 'TMPDIR    = "%c",''internal'''
   write(IDF_PO, '(3a)') 'SUBSTANCE = "%c",''', substance(1:len_trim(substance)),''''
endif
write(IDF_PO, '(a)') '#'
write(IDF_PO, '(a)') '# High performace compute center individuals in parallel'
write(IDF_PO, '(a)') '#'
write(IDF_PO, '(3a)') '#DATADIR   = "%c",''/tmp/', user_name(1:len_trim(user_name)),'''      ! Make sure these are unique'
write(IDF_PO, '(3a)') '#DATAFILE  = "%c",''',datafile(1:len_trim(datafile)),''''
write(IDF_PO, '(3a)') '#INDIDIR   = "%c",''/tmp/', user_name(1:len_trim(user_name)),'''      ! Make sure these are unique'
write(IDF_PO, '(a )') '#TMPDIR    = "%c",''internal'''
write(IDF_PO, '(3a)') '#SUBSTANCE = "%c",''', substance(1:len_trim(substance)),''''
write(IDF_PO, '(a )') '#'
write(IDF_PO, '(a )') '# High performace compute center individuals serially'
write(IDF_PO, '(a )') '#'
write(IDF_PO, '(3a)') '#DATADIR   = "%c",''/tmp/', user_name(1:len_trim(user_name)),'''      ! Make sure these are unique'
write(IDF_PO, '(3a)') '#DATAFILE  = "%c",''',datafile(1:len_trim(datafile)),''''
write(IDF_PO, '(a )') '#INDIDIR   = "%c",''kuplot'''
write(IDF_PO, '(a )') '#TMPDIR    = "%c",''internal'''
write(IDF_PO, '(3a)') '#SUBSTANCE = "%c",''', substance(1:len_trim(substance)),''''
write(IDF_PO, '(a)') '#'
write(IDF_PO, '(a)') '#'
write(IDF_PO, '(a)') '# Create necessary directories'
write(IDF_PO, '(a)') '#'
write(IDF_PO, '(a)') 'if(''kuplot''/="%c",INDIDIR) then         ! Test for INDIDIR'
write(IDF_PO, '(a)') '  fexist "%c/INDI", INDIDIR'
write(IDF_PO, '(a)') '  if(res[1]==0) then'
write(IDF_PO, '(a)') '    set error, live'
write(IDF_PO, '(a)') '    system "mkdir -p %c/INDI", INDIDIR'
write(IDF_PO, '(a)') '    set error,exit'
write(IDF_PO, '(a)') '  endif'
write(IDF_PO, '(a)') 'endif'
write(IDF_PO, '(a)') '#'
write(IDF_PO, '(a)') 'if(''.''/="%c",DATADIR) then              ! Test for DATADIR'
write(IDF_PO, '(a)') '  fexist "%c/DATA", DATADIR'
write(IDF_PO, '(a)') '  if(res[1]==0) then'
write(IDF_PO, '(a)') '    set error, live'
write(IDF_PO, '(a)') '    system "mkdir -p %c/DATA", DATADIR'
write(IDF_PO, '(a)') '    set error, exit'
write(IDF_PO, '(a)') '  endif'
write(IDF_PO, '(a)') 'endif'
!
close(IDF_PO)
!
end subroutine write_diffev_global
!
!*****7**************************************************************** 
!
subroutine write_diffev_single_part1(substance, user_name,      &
           spcgr_syst, lattice_para, compute, c_style)
!-
!  Writes the first part for a DIFFEV powder refinement
!+
!
use precision_mod
!
implicit none
!
character(len=*), intent(in) :: substance
character(len=*), intent(in) :: user_name
integer, intent(in) :: spcgr_syst
real(kind=PREC_DP) , dimension(6), intent(in) :: lattice_para
character(len=*), intent(in) :: compute                         ! Individuals 'serial' or 'parallel'
character(len=*), intent(in) :: c_style                         ! 'powder' or 'pdf'
!
character(len=PREC_STRING) :: outfile     ! output file name(s)
!
if(lexist_hkl_file) then
   outfile = substance(1:len_trim(substance)) // '.hkl'
elseif(lexist_fcf_file) then
   outfile = substance(1:len_trim(substance)) // '_merged.hkl'
else
   outfile = substance(1:len_trim(substance)) // '.h5'
endif
!
call write_diffev_global(substance, outfile, user_name, compute, c_style)
!
! Start main discus macro
!
if(lexist_hkl_file) then 
   outfile = substance(1:len_trim(substance)) // '_main_hkl.mac'
elseif(lexist_fcf_file) then 
   outfile = substance(1:len_trim(substance)) // '_main_merged.mac'
else
   outfile = substance(1:len_trim(substance)) // '_main_single.mac'
endif
open(IDF_PO, file=outfile, status='unknown')
!
write(IDF_PO, '(a)') '#'
write(IDF_PO, '(a)') 'variable real, rlambda'
write(IDF_PO, '(a)') 'read'
write(IDF_PO, '(3a)') '  stru CELL/', substance(1:len_trim(substance)), '.cell'
write(IDF_PO, '(a)') '#'
!
!call write_powder_lattice(IDF_PO, spcgr_syst)
!call write_diffev_lattice_top(spcgr_syst, lattice_para)
!
end subroutine write_diffev_single_part1
!
!*****7**************************************************************** 
!
subroutine write_diffev_powder_part1(substance, user_name,      &
           spcgr_syst, lattice_para, compute, c_style)
!-
!  Writes the first part for a DIFFEV powder refinement
!+
!
use precision_mod
!
implicit none
!
character(len=*), intent(in) :: substance
character(len=*), intent(in) :: user_name
integer, intent(in) :: spcgr_syst
real(kind=PREC_DP) , dimension(6), intent(in) :: lattice_para
character(len=*), intent(in) :: compute                         ! Individuals 'serial' or 'parallel'
character(len=*), intent(in) :: c_style                         ! 'powder' or 'pdf'
!
character(len=PREC_STRING) :: outfile     ! output file name(s)
!
outfile = substance(1:len_trim(substance)) // '.tth'
call write_diffev_global(substance, outfile, user_name, compute, c_style)
!
! Start main discus macro
!
outfile = substance(1:len_trim(substance)) // '_main_powder.mac'
open(IDF_PO, file=outfile, status='unknown')
!
write(IDF_PO, '(a)') '#'
write(IDF_PO, '(a)') 'variable real, rlambda'
write(IDF_PO, '(a)') 'read'
write(IDF_PO, '(3a)') '  stru CELL/', substance(1:len_trim(substance)), '.cell'
write(IDF_PO, '(a)') '#'
!
call write_powder_lattice(IDF_PO, spcgr_syst)
call write_diffev_lattice_top(spcgr_syst, lattice_para)
!
end subroutine write_diffev_powder_part1
!
!*****7**************************************************************** 
!
subroutine write_diffev_pdf_part1(substance, user_name,      &
           spcgr_syst, lattice_para, compute, c_style)
!-
!  Writes the first part for a DIFFEV powder refinement
!+
!
use precision_mod
!
implicit none
!
character(len=*), intent(in) :: substance
character(len=*), intent(in) :: user_name
integer, intent(in) :: spcgr_syst
real(kind=PREC_DP) , dimension(6), intent(in) :: lattice_para
character(len=*), intent(in) :: compute                         ! Individuals 'serial' or 'parallel'
character(len=*), intent(in) :: c_style                         ! 'powder' or 'pdf'
!
character(len=PREC_STRING) :: outfile     ! output file name(s)
!
outfile = substance(1:len_trim(substance)) // '.grobs'
call write_diffev_global(substance, outfile, user_name, compute, c_style)
!
! Start main discus macro
!
outfile = substance(1:len_trim(substance)) // '_main_pdf.mac'
open(IDF_PO, file=outfile, status='unknown')
!
write(IDF_PO, '(a)') '#'
write(IDF_PO, '(a)') 'variable real, rlambda'
write(IDF_PO, '(a)') 'read'
write(IDF_PO, '(3a)') '  stru CELL/', substance(1:len_trim(substance)), '.cell'
write(IDF_PO, '(a)') '#'
!
call write_powder_lattice(IDF_PO, spcgr_syst)
call write_diffev_lattice_top(spcgr_syst, lattice_para)
!
end subroutine write_diffev_pdf_part1
!
!*****7**************************************************************** 
!
subroutine write_diffev_lattice_top(spcgr_syst, lattice_para)
!-
!  Top level for diffev_lattice parameters, used by powder and pdf
!+
!
use precision_mod
!
implicit none
!
integer, intent(in) :: spcgr_syst
real(kind=PREC_DP) , dimension(6), intent(in) :: lattice_para
!
if(spcgr_syst==1) then         ! Triclinic
   call write_diffev_lattice('P_alat' , lattice_para(1))
   call write_diffev_lattice('P_blat' , lattice_para(2))
   call write_diffev_lattice('P_clat' , lattice_para(3))
   call write_diffev_lattice('P_alpha', lattice_para(4))
   call write_diffev_lattice('P_beta' , lattice_para(5))
   call write_diffev_lattice('P_gamma', lattice_para(6))
elseif(spcgr_syst==2) then         ! Monoclinic b unique
   call write_diffev_lattice('P_alat' , lattice_para(1))
   call write_diffev_lattice('P_blat' , lattice_para(2))
   call write_diffev_lattice('P_clat' , lattice_para(3))
   call write_diffev_lattice('P_beta' , lattice_para(5))
elseif(spcgr_syst==3) then         ! Monoclinic c unique
   call write_diffev_lattice('P_alat' , lattice_para(1))
   call write_diffev_lattice('P_blat' , lattice_para(2))
   call write_diffev_lattice('P_clat' , lattice_para(3))
   call write_diffev_lattice('P_gamma', lattice_para(6))
elseif(spcgr_syst==4) then         ! Orthorhombic
   call write_diffev_lattice('P_alat' , lattice_para(1))
   call write_diffev_lattice('P_blat' , lattice_para(2))
   call write_diffev_lattice('P_clat' , lattice_para(3))
elseif(spcgr_syst==5) then         ! Tetragonal
   call write_diffev_lattice('P_alat' , lattice_para(1))
   call write_diffev_lattice('P_clat' , lattice_para(3))
elseif(spcgr_syst==6 .or. spcgr_syst==8) then         ! Trigonal, Hexagonal
   call write_diffev_lattice('P_alat' , lattice_para(1))
   call write_diffev_lattice('P_clat' , lattice_para(3))
elseif(spcgr_syst==7) then         ! Rhombohedral
   call write_diffev_lattice('P_alat' , lattice_para(1))
   call write_diffev_lattice('P_alpha', lattice_para(4))
elseif(spcgr_syst==9) then         ! Cubic
   call write_diffev_lattice('P_alat' , lattice_para(1))
endif
!
end subroutine write_diffev_lattice_top
!
!*****7**************************************************************** 
!
subroutine write_diffev_lattice(p_name, lp) 
!-
!  Write the diffev newpara line for lattice parameters
!  Values are adjusted to a small window of +-1 percent
!+
!
use precision_mod
!
implicit none
!
character(len=*)  , intent(in) :: p_name
real(kind=PREC_DP), intent(in) :: lp   ! lattice parameter
!
character(len=PREC_STRING) :: line 
real(kind=PREC_DP) :: l_low
real(kind=PREC_DP) :: l_hig
!
l_low = real((int(100*lp*0.99)   )*0.010D0, kind=PREC_DP)
l_hig = real((int(100*lp*0.99)+10)*0.010D0, kind=PREC_DP)
write(IDF_NP,'(a9, a7,4(a,f9.4),a)') 'newpara  ', p_name,      &
   ', range:[', l_low,    ', ', l_hig, '], value:[', lp*0.995, ', ', lp*1.00,'], status:free'
!
!write(IDF_FX, '(2a)') 'refine ', p_name
!
write(IDF_BS, '(2a)') 'variable real, ', p_name
diffev_dimension = diffev_dimension + 1
!
line = '                 = '
write(line(1:16), '(a)') p_name
write(IDF_BS2, '(a19,G18.10E3)') line(1:19), lp
!
end subroutine write_diffev_lattice
!
!*****7**************************************************************** 
!
subroutine write_diffev_single_part2(j, l, jj, iscat, lcontent, natoms, c_atom, uij_l, l_flag)
!
!  Write main discus macro; single crystal PART 2
!  Writes Uij or Ueqv
!  can be used for SINGLE and POWDER
!
!
use precision_mod
!
implicit none
!
integer            , intent(in) :: j
integer            , intent(in) :: l
integer            , intent(in) :: jj
integer            , intent(in) :: iscat
integer            , intent(in) :: lcontent
integer            , intent(in) :: natoms
character(len=4), dimension(lcontent), intent(in)  :: c_atom    ! Atom types 
real(kind=PREC_DP), dimension(6,natoms, natoms), intent(in) :: uij_l
logical            , intent(in) :: l_flag
!character(len=*)   , intent(in) :: c_flag
!
character(len=PREC_STRING) :: line 
character(len=5          ) :: c_flag 
integer :: m, k
integer :: ref_n
real(kind=PREC_DP), dimension(4,2) :: ref_window
real(kind=PREC_DP), dimension(4,2) :: ref_min
!
ref_window(:,1) = 0.00_PREC_DP
ref_window(1,2) = 0.10_PREC_DP
ref_window(2,2) = 0.10_PREC_DP
ref_window(3,2) = 0.05_PREC_DP
ref_window(4,2) = 0.05_PREC_DP
!
ref_min   (:,1) = 0.00000_PREC_DP
ref_min   (1,2) = 0.00010_PREC_DP
ref_min   (2,2) = 0.00010_PREC_DP
ref_min   (3,2) = 0.00005_PREC_DP
ref_min   (4,2) = 0.00005_PREC_DP

!
!if(c_flag=='fixed') then
if(l_flag         ) then
  ref_n = 2
  c_flag = 'free'
else
  ref_n = 1
  c_flag = 'fixed '
endif
if(l==1) then
   write(IDF_PO,'(a,i3,4a)') 'anis type:', jj,', values:[', 'U_',c_atom(iscat)(1:len_trim(c_atom(iscat))),'_1]'
   write(IDF_NP,'(3a,4(a,g15.8e3), 2a)') 'newpara U_',c_atom(iscat)(1:len_trim(c_atom(iscat))),'_1' &
        , ', range:[', uij_l(1,iscat, j)-max(ref_window(1,ref_n)*uij_l(1,iscat, j), ref_min(1,ref_n)) &
        , ', ', uij_l(1,iscat, j)+max(ref_window(2,ref_n)*uij_l(1,iscat, j), ref_min(2,ref_n)) &
        , '], value:[', uij_l(1,iscat, j)-max(ref_window(3,ref_n)*uij_l(1,iscat, j), ref_min(3,ref_n)) &
        , ', ', uij_l(1,iscat, j)+max(ref_window(4,ref_n)*uij_l(1,iscat, j), ref_min(4,ref_n))   &
        , '], status:',c_flag
   write(IDF_BS, '(4a)') 'variable real, ', 'U_',c_atom(iscat)(1:len_trim(c_atom(iscat))),'_1'
   line = '                 = '
   write(line(1:16), '(3a)') 'U_',c_atom(iscat)(1:len_trim(c_atom(iscat))),'_1'
   write(IDF_BS2,'(a19, G18.10E3)') line(1:19), uij_l(1,iscat, j)
   diffev_dimension = diffev_dimension + 1
!  if(l_flag) write(IDF_FX , '(3a,i1.1)') 'refine U',c_atom(iscat)(1:len_trim(c_atom(iscat))),'_',1
elseif(l==6) then
   write(IDF_PO,'(a,i3,a,5(3a,i1.1,a2),3a)') 'anis type:', jj,', values:[', &
        ('U_',c_atom(iscat)(1:len_trim(c_atom(iscat))),'_',k,', ',k=1,5), &
         'U_',c_atom(iscat)(1:len_trim(c_atom(iscat))),'_6]'
   do m=1,6
      write(IDF_NP,'(3a,i1.1,4(a,g15.8e3),2a)') 'newpara U_',c_atom(iscat)(1:len_trim(c_atom(iscat))),'_',m &
        , ', range:[', uij_l(m,iscat, j)-max(ref_window(1,ref_n)*uij_l(1,iscat, j), ref_min(1,ref_n)) &
        , ', ', uij_l(m,iscat, j)+max(ref_window(2,ref_n)*uij_l(1,iscat, j), ref_min(2,ref_n)) &
        , '], value:[', uij_l(m,iscat, j)-max(ref_window(3,ref_n)*uij_l(1,iscat, j), ref_min(3,ref_n)) &
        , ', ', uij_l(m,iscat, j)+max(ref_window(4,ref_n)*uij_l(1,iscat, j), ref_min(4,ref_n))         &
        , '], status:', c_flag
      write(IDF_BS, '(4a, i1.1)') 'variable real, ', 'U_',c_atom(iscat)(1:len_trim(c_atom(iscat))),'_',m
      line = '                 = '
      write(line(1:16), '(3a,i1.1)') 'U_',c_atom(iscat)(1:len_trim(c_atom(iscat))),'_',m
      write(IDF_BS2,'(a19, G18.10E3)') line(1:19), uij_l(m,iscat, j)
      diffev_dimension = diffev_dimension + 1
!     if(l_flag) write(IDF_FX , '(3a,i1.1)') 'refine U',c_atom(iscat)(1:len_trim(c_atom(iscat))),'_',m
   enddo
endif
!
end subroutine write_diffev_single_part2
!
!*****7**************************************************************** 
!
subroutine write_diffev_single_part3(&
   inline       , iatom, iscat, lcontent, natoms, c_atom, posit, &
   l_flag)
!
!  Write main discus macro; single crystal PART 3
!  can be used for SINGLE and POWDER
!
use spcgr_apply, only: get_wyckoff
use wyckoff_mod
use precision_mod
!
implicit none
!
character(len=*)   , intent(in) :: inline
integer            , intent(in) :: iatom
integer            , intent(in) :: iscat
integer            , intent(in) :: lcontent
integer            , intent(in) :: natoms
character(len=4), dimension(lcontent), intent(in)  :: c_atom    ! Atom types 
real(kind=PREC_DP), dimension(3,natoms), intent(in) :: posit
logical            , intent(in) :: l_flag
!
!character(len=5          ) :: c_flag 
character(len=1), dimension(3), parameter :: xyz_name = (/ 'x', 'y', 'z'/)
character(len=2), dimension(3), parameter :: xyz_appe = (/ '_x', '_y', '_z'/)
integer, parameter :: mode = 0
logical, parameter :: loutput = .false.
character(len=PREC_STRING) :: line 
character(len=PREC_STRING) :: string
integer :: i, j, k, l
real(kind=PREC_DP), dimension(3) :: coordinates
!
!if(l_flag         ) then
!  c_flag = 'fixed'
!else
!  c_flag = 'free '
!endif
!
coordinates = posit(:, iatom)
call get_wyckoff(coordinates, loutput, mode)
!
i = len_trim(inline)
loop_xyz: do j=1, 3         ! Loop over x, y, z
   if(wyc_fix(j)) then               ! Fixed or related
      if(all(wyc_fix_mat(j,1:3)==0)) cycle loop_xyz    ! Numerically fixed coordinate
!                                                      ! related coordinate
      write(string,'(2a1,i3,a4)') xyz_name(j),'[',iatom, '] = ' 
      l = 10   ! next free field in string
      do k =1, j-1
         if(wyc_fix_mat(j,k)/=0) then
            write(string(l:l+9), '(i2,3a1,i3,a1)') wyc_fix_mat(j,k), '*',xyz_name(k),'[',iatom, ']'
            l = l + 10
         endif
      enddo
      if(wyc_fix_mat(j,4)>0) then
         write(string(l:l+14), '(a2,f12.8)') ' +',    wyc_fix_mat(j,4)  / 24.0_prec_dp
      elseif(wyc_fix_mat(j,4)<0) then
         write(string(l:l+14), '(a2,f12.8)') ' -',abs(wyc_fix_mat(j,4)) / 24.0_prec_dp
      endif
      write(IDF_PO, '(a)') string(1:len_trim(string))
   else
!
!     Free variable
!
      write(IDF_PO,'(2a1,i3,3a)') xyz_name(j),'[',iatom, '] = P_',inline(1:i), xyz_appe(j)
!
      if(.not.l_flag) then                        ! Powder diffraction, initially fix the positions
         write(IDF_NP,'(3a,4(a,f12.8),a)') 'newpara P_', &
         c_atom(iscat)(1:len_trim(c_atom(iscat))),xyz_appe(j), &
          ', range:[', posit(j, iatom)     , ', ', posit(j, iatom)     , &
         '], value:[', posit(j, iatom)     , ', ', posit(j, iatom)     , &
         '], status:fixed'   
      else                                            ! Single crystal, free atom coordinates
         write(IDF_NP,'(3a,4(a,f12.8),a)') 'newpara P_', &
         c_atom(iscat)(1:len_trim(c_atom(iscat))),xyz_appe(j), &
          ', range:[', posit(j, iatom)-0.02, ', ', posit(j, iatom)+0.02, &
         '], value:[', posit(j, iatom)-0.01, ', ', posit(j, iatom)+0.01, &
         '], status:free'   
!        write(IDF_FX, '(3a)') 'refine P_', c_atom(iscat)(1:len_trim(c_atom(iscat))),xyz_appe(j)
      endif
!
!     Positions and variables in diffev_best*.mac
!
      write(IDF_BS,'(3a)') 'variable real, P_',c_atom(iscat)(1:len_trim(c_atom(iscat))),xyz_appe(j)
      line = '                 = '
      write(line(1:16), '(3a)') 'P_',c_atom(iscat)(1:len_trim(c_atom(iscat))),xyz_appe(j)
      write(IDF_BS2,'(a19, G18.10E3)') line(1:19), posit(j, iatom)
!
      diffev_dimension = diffev_dimension + 1       ! Adjust the number of members/children in DIFFEV       ! Adjust the number of members/children in DIFFEV
   endif
enddo loop_xyz
!
end subroutine write_diffev_single_part3
!
!*****7**************************************************************** 
!
subroutine write_diffev_single_part4(&
        substance,  ilist, c_form, lambda, slambda, compute, c_style,           &
        P_exti, hkl_max, P_ncycle  &
        )
!
!  Write main discus macro; single crystal PART 4
!  REFINE version
!
use blanks_mod
use precision_mod
!
implicit none
!
character(len=*)                , intent(in) :: substance
integer                         , intent(in) :: ilist
character(len=*)                , intent(in) :: c_form
character(len=*)                , intent(in) ::  lambda
character(len=*)                , intent(in) :: slambda
character(len=*)                , intent(in) :: compute
character(len=*)                , intent(in) :: c_style
real(kind=PREC_DP)              , intent(in) :: P_exti
integer           , dimension(3), intent(inout) :: hkl_max
integer                         , intent(in) :: P_ncycle
!
character(len=PREC_STRING) :: aspher_file
character(len=PREC_STRING) :: line
character(len=PREC_STRING) :: line1
character(len=PREC_STRING) :: line2
integer :: i,j, k
integer :: length
integer :: ios
logical :: lexist
real(kind=PREC_DP) :: rlambda
!
read(lambda, *) rlambda
!
j = 0
i = len_trim(substance)
call single_main_part4(IDF_PO, .false., substance, c_form, rlambda, slambda, P_exti)
close(IDF_PO)
!
if(P_exti> 0.00001) then
   write(IDF_NP,'(a,3(f8.5,a), f8.5, a)') 'newpara P_exti, range:[', &
                 P_exti*0.90, ', ', P_exti*1.10, '], value:[',      &
                 P_exti*0.95, ', ', P_exti*1.05, '], status:free'
else
   write(IDF_NP,'(a,3(f8.5,a), f8.5,a)') 'newpara P_exti, range:[', &
                 P_exti     , ', ', P_exti     ,'], value:[',       &
                 P_exti     , ', ', P_exti     ,'], status:fixed'   
endif
write(IDF_MA,'(a )') '#'
!
call diffev_setup_last
!write(IDF_MA,'(a )') 'exit  ! Back to SUITE'
close(IDF_MA)
!
open(IDF_MA, file='final_cell_cif.mac', status='unknown')
write(IDF_MA,'( a)') 'branch discus'
write(IDF_MA,'( a)') 'read'
write(IDF_MA,'(3a)') '  stru internal.' , substance(1:i), '.cell'
write(IDF_MA,'( a)') 'save'
write(IDF_MA,'(3a)') '  outf ', substance(1:i), '.cell'
write(IDF_MA,'( a)') '  write all'
write(IDF_MA,'( a)') '  run'
write(IDF_MA,'( a)') 'exit'
write(IDF_MA,'(3a)') 'export cif, ', substance(1:i), '.cif, spcgr:original'
write(IDF_MA,'( a)') 'exit  ! back to REFINE'
!
close(IDF_MA)
!
call write_prepare_hkl_full(substance, ilist, hkl_max, rlambda, slambda)
!
!  Interpret "aspher.json"
!
if(c_form=='discamb') then
   aspher_file = '../BUILD_TSC/aspher.json'  ! Use the special fixed name for aspher.json
   inquire(file=aspher_file, exist = lexist)
   cond_exist: if(lexist) then
      open(unit=IDF_MA, file=aspher_file, status='old')
      read(IDF_MA, '(a)', iostat=ios) line
      if(is_iostat_end(ios)) exit cond_exist
      if(ios/= 0           ) exit cond_exist
      if(line=='}'         ) exit cond_exist
      length=len_trim(line)
      call rem_leading_bl(line, length)
      loop_aspher: do
         i =index(line, 'structure')
         if(i > 0) then
            i = index(line, ':')
            j = index(line(i+1:len_trim(line)), '"')
            k = index(line, '"', .true.)
            line1 = line(i+j+1:k-1)
            line2 = 'cp ../ESSENTIAL_INPUT/' // line1(1:len_trim(line1)) // ' ../BUILD_TSC/'
            call execute_command_line(line2, wait=.false.)
            exit cond_exist
         endif
         read(IDF_MA, '(a)', iostat=ios) line
         if(is_iostat_end(ios)) exit cond_exist
         if(ios/= 0           ) exit cond_exist
         if(line=='}'         ) exit cond_exist
         length=len_trim(line)
         call rem_leading_bl(line, length)
      enddo loop_aspher
   endif cond_exist
   close(IDF_MA)
endif
!
call diffev_best_single_mac(.true., compute)
write(IDF_BS2,'(a)') 'P_scale          =  1.000'
!
close(IDF_BS)
close(IDF_BS2)
close(IDF_PO)
close(IDF_MA)
close(IDF_NP)
!
call write_kup_diffev_mac_hkl(substance, compute, c_style)
call single_k_inter_macro(IRE_MA, .false., substance)          ! Write "ksingle.mac"
!
write(line,'(3a)') 'cp ', substance(1:len_trim(substance)),'.cell CELL'
call execute_command_line(line, wait=.false.)
!
end subroutine write_diffev_single_part4
!
!*****7**************************************************************** 
!
subroutine write_diffev_powder_part4(&
           substance, c_form, lambda, slambda, compute, c_style )
!-
! Write main discus macro DIFFEV_VERSION, last part
! Write "diffev_setup.mac"              , last part
!+
!
use precision_mod
!
implicit none
!
character(len=*), intent(in) :: substance    ! current material
character(len=*), intent(in) :: c_form       ! Atom form factor
character(len=*), intent(in) ::  lambda      ! Wave length numerical value
character(len=*), intent(in) :: slambda      ! Wave length symbol
character(len=*), intent(in) :: compute      ! 'serial' or 'parallel' 
character(len=*), intent(in) :: c_style      ! 'powder' or 'pdf' 
!
character(len=PREC_STRING) :: line
real(kind=PREC_DP) :: rlambda
!
write(IDF_PO, '(a)') '#'
write(IDF_PO, '(a)') 'save'
write(IDF_PO, '(3a)') '  outfile internal.', substance(1:len_trim(substance)), '.cell'
write(IDF_PO, '(a)') '  write all'
write(IDF_PO, '(a)') '  run'
write(IDF_PO, '(a)') 'exit'
write(IDF_PO, '(a)') '#'
!
write(IDF_PO, '(a)') 'read'
write(IDF_PO, '(3a)') '  cell internal.', substance(1:len_trim(substance)), '.cell'
write(IDF_PO, '(a)') '#'
write(IDF_PO, '(a)') '@powder.mac'
!
read(lambda, *) rlambda
call write_macro_powder(c_form, hkl_file, rlambda, slambda)
!
write(IDF_PO, '(a)') 'output'
write(IDF_PO, '(a)') '  outf "%c/INDI/indi.%4D.%4D", INDIDIR, REF_KID, REF_INDIV'
write(IDF_PO, '(a)') '  value inte'
write(IDF_PO, '(a)') '  form powder, tth, F_XMIN, F_XMAX, F_XSTP'
write(IDF_PO, '(a)') '  run'
write(IDF_PO, '(a)') 'exit'
!
call diffev_main_powder_mac(.true.)
!
call diffev_best_powder_mac(.true., compute)
call diffev_setup_last
!
call diffev_best2_powder_mac
!
close(IDF_PO)
close(IDF_MA)
close(IDF_NP)
close(IDF_BS)
close(IDF_BS2)
!
call write_kup_diffev_mac(substance, compute, c_style)
!
write(line,'(3a)') 'cp ', substance(1:len_trim(substance)),'.cell CELL'
call execute_command_line(line, wait=.false.)
!
end subroutine write_diffev_powder_part4
!
!*****7**************************************************************** 
!
subroutine write_diffev_pdf_part4(&
           substance, c_form, lambda, slambda, compute, c_style )
!-
! Write main discus macro DIFFEV_VERSION, last part
! PDF version
! Write "diffev_setup.mac"              , last part
!+
!
use precision_mod
!
implicit none
!
character(len=*), intent(in) :: substance    ! current material
character(len=*), intent(in) :: c_form       ! Atom form factor
character(len=*), intent(in) ::  lambda      ! Wave length numerical value
character(len=*), intent(in) :: slambda      ! Wave length symbol
character(len=*), intent(in) :: compute      ! 'serial' or 'parallel' 
character(len=*), intent(in) :: c_style      ! 'powder' or 'pdf' 
!
character(len=PREC_STRING) :: line
real(kind=PREC_DP) :: rlambda
!
write(IDF_PO, '(a)') '#'
write(IDF_PO, '(a)') 'save'
write(IDF_PO, '(3a)') '  outfile internal.', substance(1:len_trim(substance)), '.cell'
write(IDF_PO, '(a)') '  write all'
write(IDF_PO, '(a)') '  run'
write(IDF_PO, '(a)') 'exit'
write(IDF_PO, '(a)') '#'
!
write(IDF_PO, '(a)') '#'
write(IDF_PO, '(a)') 'variable integer, ncellx'
write(IDF_PO, '(a)') 'variable integer, ncelly'
write(IDF_PO, '(a)') 'variable integer, ncellz'
write(IDF_PO, '(a)') 'ncellx = int(P_diam_a/lat[1]) + 4'
write(IDF_PO, '(a)') 'ncelly = int(P_diam_b/lat[2]) + 4'
write(IDF_PO, '(a)') 'ncellz = int(P_diam_c/lat[2]) + 4'
write(IDF_PO, '(a)') 'read'
write(IDF_PO, '(3a)') '  cell internal.', substance(1:len_trim(substance)), '.cell, ncellx, ncelly, ncellz' 
write(IDF_PO, '(a)') '#'
write(IDF_PO, '(a)') 'surface'
write(IDF_PO, '(a)') '  boundary ellipsoid, P_diam_a, P_diam_b, P_diam_c'
write(IDF_PO, '(a)') 'exit'
write(IDF_PO, '(a)') '#'
write(IDF_PO, '(a)') 'purge type:yes'
write(IDF_PO, '(a)') '#'
write(IDF_PO, '(a)') '#'
write(IDF_PO, '(a)') '@debye.mac'
!
read(lambda, *) rlambda
call write_macro_debye(c_form, hkl_file, rlambda, slambda)
!
write(IDF_PO, '(a)') 'output'
write(IDF_PO, '(a)') '  outf "%c/INDI/indi.%4D.%4D", INDIDIR, REF_KID, REF_INDIV'
write(IDF_PO, '(a)') '  value PDF'
write(IDF_PO, '(a)') '  form pdf, r, F_XMIN, F_XMAX, F_XSTP'
write(IDF_PO, '(a)') '  run'
write(IDF_PO, '(a)') 'exit'
!
!   PDF Additional parameters P_diam_a, P_diam_b, P_diam_c
!
write(IDF_NP, '(a)') 'newpara  P_diam_a , range:[20.0000, 100.0000], value:[ 20.0000,   30.0000], status:free     ! Diameter  a'
write(IDF_NP, '(a)') 'newpara  P_diam_b , range:[20.0000, 100.0000], value:[ 30.0000,   40.0000], status:free     ! Diameter  b'
write(IDF_NP, '(a)') 'newpara  P_diam_c , range:[20.0000, 100.0000], value:[ 50.0000,   60.0000], status:free     ! Diameter  c'
!
write(IDF_BS, '(a)') 'variable real, P_diam_a'
write(IDF_BS, '(a)') 'variable real, P_diam_b'
write(IDF_BS, '(a)') 'variable real, P_diam_c'
!
write(IDF_BS2, '(a)') 'P_diam_a         =   25.0000'
write(IDF_BS2, '(a)') 'P_diam_b         =   35.0000'
write(IDF_BS2, '(a)') 'P_diam_c         =   55.0000'
!
!write(IDF_FX , '(a)') 'refine P_diam_a'
!write(IDF_FX , '(a)') 'refine P_diam_b'
!write(IDF_FX , '(a)') 'refine P_diam_c'
!
diffev_dimension = diffev_dimension + 3
!
call diffev_main_powder_mac(.false.)
!
call diffev_best_powder_mac(.false., compute)
call diffev_setup_last
!
call diffev_best2_powder_mac
!
close(IDF_PO)
close(IDF_MA)
close(IDF_NP)
close(IDF_BS)
close(IDF_BS2)
!
call write_kup_diffev_mac(substance, compute, c_style)         ! Write "kup.diffev.mac"
!
close(IDF_PO)
!
write(line,'(3a)') 'cp ', substance(1:len_trim(substance)),'.cell CELL'
call execute_command_line(line, wait=.false.)
!
end subroutine write_diffev_pdf_part4
!
!*****7**************************************************************** 
!
subroutine write_kup_diffev_mac(substance, compute, c_style)
!-
!  Write the kup.diffev.mac
!+
!
implicit none
!
character(len=*), intent(in) :: substance    ! A the material
character(len=*), intent(in) :: compute      ! Atom form factor
character(len=*), intent(in) :: c_style      ! 'powder' or 'pdf' 
!
open(unit=IDF_PO, file='kup.diffev.mac', status='unknown')
!
write(IDF_PO, '(a)') '# kup.diffev.mac '
write(IDF_PO, '(a)') '# Merges data sets does background calculation and spline '
if(compute=='parallel') then
   write(IDF_PO, '(a)') '# Version for parallel computing of individual repetitions'
elseif(compute=='serial') then
   write(IDF_PO, '(a)') '# Version for serial computing of individual repetitions'
endif
if(c_style=='powder') then
   write(IDF_PO, '(a)') '# Version for fit to an experimental powder pattern'
elseif(c_style=='pdf') then
   write(IDF_PO, '(a)') '# Version for fit to an experimental PDF'
endif
!
write(IDF_PO, '(a)') '#'
write(IDF_PO, '(a)') '#   To avoid race conditions individual copies of the data dset are created'
write(IDF_PO, '(a)') 'fexist "%c/DATA/%c.%4D", DATADIR, DATAFILE, REF_KID'
write(IDF_PO, '(a)') 'if(res[1].eq.0) then  ! File does not exist load up to DATADIR'
write(IDF_PO, '(a)') '   system "cp DATA/%c %c/DATA/%c.%4D", DATAFILE, DATADIR, DATAFILE, REF_KID'
write(IDF_PO, '(a)') 'endif'
write(IDF_PO, '(a)') '#'
write(IDF_PO, '(a)') 'variable integer, idata   ! Data set number with experimental data'
write(IDF_PO, '(a)') 'variable integer, imerge  ! Data set number with merged calc data'
write(IDF_PO, '(a)') 'variable integer, ispline ! Data set number after spline'
write(IDF_PO, '(a)') 'variable integer, icalc   ! final calc data set number after background etc '
write(IDF_PO, '(a)') '#'
if(compute=='parallel') then
   write(IDF_PO, '(a)') 'reset'
   write(IDF_PO, '(a)') '@global.mac'
   write(IDF_PO, '(a)') '#'
   write(IDF_PO, '(a)') 'load xy, "%c/INDI/%c.%4D*", INDIDIR, SUBSTANCE, REF_KID'
elseif(compute=='serial') then
   write(IDF_PO, '(a)') 'if("%c", INDIDIR /= ''kuplot'') then ! Data are on disk'
   write(IDF_PO, '(a)') '  load xy, "%c/INDI/indi.%4D*", INDIDIR, REF_KID'
   write(IDF_PO, '(a)') 'endif'
endif
write(IDF_PO, '(a)') 'merge all'
write(IDF_PO, '(a)') 'imerge = n[1]             ! Data set number with merged calc data'
!
if(c_style=='powder') then
   write(IDF_PO, '(a)') 'load xy, "%c/DATA/%c.%4D", DATADIR, DATAFILE, REF_KID'
   write(IDF_PO, '(a)') 'idata = n[1]              ! Data set number with experimental data'
   write(IDF_PO, '(a)') 'spline imerge, idata      ! For safety spline ==> isplined'
   write(IDF_PO, '(a)') 'ispline = n[1]            ! Data set number after spline'
   write(IDF_PO, '(a)') 'fit idata                 ! Fit a background to experimental data'
   write(IDF_PO, '(a)') '  func back, ispline, 5   ! Use merged data and cubed polynomial'
   write(IDF_PO, '(a)') '  para 1, 0, 1.0          ! Scale factor'
   write(IDF_PO, '(a)') '  para 2, 1, ymin[2]      ! contant background'
   write(IDF_PO, '(a)') '  para 3, 1, 0.0          ! linear parameter'
   write(IDF_PO, '(a)') '  para 4, 1, 0.0          ! square parameter'
   write(IDF_PO, '(a)') '  para 5, 1, 0.0          ! cubed parameter'
   write(IDF_PO, '(a)') '  cycle 10'
   write(IDF_PO, '(a)') '  run                     ! create two data sets scaled+back and difference'      
   write(IDF_PO, '(a)') 'exit'
   write(IDF_PO, '(a)') 'icalc = n[1] - 1          ! Final calculated data set is scaled + back'
elseif(c_style=='pdf') then
   write(IDF_PO, '(a)') 'load xy, "%c/DATA/%c", DATADIR, DATAFILE'
   write(IDF_PO, '(a)') 'idata = n[1]              ! Data set number with experimental data'
   write(IDF_PO, '(a)') 'spline imerge, idata      ! For safety spline ==> isplined'
   write(IDF_PO, '(a)') 'ispline = n[1]            ! Data set number after spline'
   write(IDF_PO, '(a)') 'icalc = n[1]              ! Final calculated data set is merged'
endif
write(IDF_PO, '(a)') 'scale'
write(IDF_PO, '(a)') 'mark'
write(IDF_PO, '(a)') 'ksav icalc'
write(IDF_PO, '(a)') '  outfile "TEMP/calc.%4D", REF_KID'
write(IDF_PO, '(a)') '  form xy'
write(IDF_PO, '(a)') '  run'
write(IDF_PO, '(a)') '#'
write(IDF_PO, '(a)') 'rval idata, icalc, dat'
write(IDF_PO, '(a)') 'exit'
write(IDF_PO, '(a)') '#'
write(IDF_PO, '(a)') 'exit   ! Back to DISCUS'
!
close(IDF_PO)
!
!  Write ksingle.mac
open(unit=IDF_PO, file='ksingle.mac', status='unknown')
!
   write(IDF_PO, '(a)') 'reset'
   write(IDF_PO, '(a)') '@global.mac'
   write(IDF_PO, '(a)') 'load xy, "%c/DATA/%c", DATADIR, DATAFILE'
   write(IDF_PO, '(a)') 'load xy, "FINAL/final.%4D", $1'
   write(IDF_PO, '(a)') 'rval 1, 2, dat'
   write(IDF_PO, '(a)') 'r[101] = res[2]'
   write(IDF_PO, '(a)') 'kcal sub, 1, 2'
   write(IDF_PO, '(a)') 'lcol 1, blue'
   write(IDF_PO, '(a)') 'lcol 2, red'
   write(IDF_PO, '(a)') 'lcol 3, black'
   write(IDF_PO, '(a)') 'mtyp 1, 0'
   write(IDF_PO, '(a)') 'mtyp 2, 0'
   write(IDF_PO, '(a)') 'mtyp 3, 0'
   write(IDF_PO, '(a)') 'etyp 1, 0'
   write(IDF_PO, '(a)') 'etyp 2, 0'
   write(IDF_PO, '(a)') 'etyp 3, 0'
if(c_style=='powder') then
   write(IDF_PO, '(a)') 'r[1] = ymax[3]'
   write(IDF_PO, '(a)') 'ccal add, wy, 3, -r[1]'
   write(IDF_PO, '(a)') 'achx 2Theta'
   write(IDF_PO, '(a)') 'achy intensity'
   write(IDF_PO, '(a)') 'tit1 "Fit to powder: %c", SUBSTANCE'
elseif(c_style=='pdf') then
   write(IDF_PO, '(a)') 'r[1] = ymax[3] - min(ymin[1],ymin[2])'
   write(IDF_PO, '(a)') 'ccal add, wy, 3, -r[1]'
   write(IDF_PO, '(a)') 'achx Distance [\A]'
   write(IDF_PO, '(a)') 'achy G(r)'
   write(IDF_PO, '(a)') 'tit1 "Fit to PDF: %c", SUBSTANCE'
endif
   write(IDF_PO, '(a)') 'tit2 "R-value %8.4f", r[101]'
   write(IDF_PO, '(a)') 'scale xmin[1], xmax[1], ymin[3]*1.05, max(ymax[1], ymax[2])*1.05'
   write(IDF_PO, '(a)') 'mark'
   write(IDF_PO, '(a)') 'fnam off'
   write(IDF_PO, '(a)') 'plot'
   write(IDF_PO, '(a)') 'rval 1, 2, dat'
!
!
close(IDF_PO)
!
end subroutine write_kup_diffev_mac
!
!*****7**************************************************************** 
!
subroutine write_kup_diffev_mac_hkl(substance, compute, c_style)
!-
!  Write the kup.diffev.mac
!  hklf4 version
!+
!
implicit none
!
character(len=*), intent(in) :: substance    ! A the material
character(len=*), intent(in) :: compute      ! Atom form factor
character(len=*), intent(in) :: c_style      ! 'powder' or 'pdf' 
!
open(unit=IDF_PO, file='kup.diffev.mac', status='unknown')
!
write(IDF_PO, '(a)') '# kup.diffev.mac '
write(IDF_PO, '(a)') '# '
if(compute=='parallel') then
   write(IDF_PO, '(a)') '# Version for parallel computing of individual repetitions'
elseif(compute=='serial') then
   write(IDF_PO, '(a)') '# Version for serial computing of individual repetitions'
endif
write(IDF_PO, '(a)') '# Version for fit to an experimental single crystal pattern'
write(IDF_PO, '(a)') '#'
write(IDF_PO, '(a)') '#   To avoid race conditions individual copies of the data dset are created'
write(IDF_PO, '(a)') 'fexist "%c/DATA/%c.%4D", DATADIR, DATAFILE, REF_KID'
write(IDF_PO, '(a)') 'if(res[1].eq.0) then  ! File does not exist load up to DATADIR'
write(IDF_PO, '(a)') '   system "cp DATA/%c %c/DATA/%c.%4D", DATAFILE, DATADIR, DATAFILE, REF_KID'
write(IDF_PO, '(a)') 'endif'
write(IDF_PO, '(a)') '#'
!
!if(compute=='parallel') then
write(IDF_PO, '(a)') '@global.mac'
write(IDF_PO, '(a)') '#'
write(IDF_PO, '(a)') 'variable integer, idata   ! Data set number with experimental data'
write(IDF_PO, '(a)') 'variable integer, icalc   ! final calc data set number after background etc '
write(IDF_PO, '(a)') '#'
write(IDF_PO, '(a)') '#'
write(IDF_PO, '(a)') 'if("%c", INDIDIR /=''kuplot'') then'
write(IDF_PO, '(a)') '  reset'
if(lexist_hkl_file .or. lexist_fcf_file)  then
   write(IDF_PO, '(a)') '  load hklf4, "%c/INDI/indi.%4D", INDIDIR, REF_KID'
   write(IDF_PO, '(a)') 'endif'
   write(IDF_PO, '(a)') 'icalc = n[1]              ! Data set number with experimental data'
   write(IDF_PO, '(a)') 'load hklf4, "%c/DATA/%c.%4D", DATADIR, DATAFILE, REF_KID'
else
   write(IDF_PO, '(a)') '  load h5, "%c/INDI/indi.%4D", INDIDIR, REF_KID'
   write(IDF_PO, '(a)') 'endif'
   write(IDF_PO, '(a)') 'icalc = n[1]              ! Data set number with experimental data'
   write(IDF_PO, '(a)') 'load h5, "%c/DATA/%c.%4D", DATADIR, DATAFILE, REF_KID'
endif
write(IDF_PO, '(a)') 'idata = n[1]              ! Data set number with experimental data'
write(IDF_PO, '(a)') 'scale'
write(IDF_PO, '(a)') 'if("%c", INDIDIR /=''kuplot'') then'
write(IDF_PO, '(a)') '   system "cp %c/INDI/indi.%4D TEMP/calc.%4D", INDIDIR, REF_KID, REF_KID'
write(IDF_PO, '(a)') 'endif'
write(IDF_PO, '(a)') '#'
write(IDF_PO, '(a)') 'rval idata, icalc, dat'
write(IDF_PO, '(a)') 'exit'
write(IDF_PO, '(a)') '#'
write(IDF_PO, '(a)') 'exit   ! Back to DISCUS'
!
close(IDF_PO)
!
!
end subroutine write_kup_diffev_mac_hkl
!
!*****7**************************************************************** 
!
subroutine diffev_main_powder_mac(is_powder)
!-
!   Writes a generic diffev_main powder part, valid for POWDER and PDF
!+
!
implicit none
!
logical, intent(in) :: is_powder
!
write(IDF_NP, '(a)') 'newpara  P_eta    ,  range:[0.5000,   0.5000], value:[  0.5000,    0.5000], status:fixed     ! Profile eta'
write(IDF_NP, '(a)') 'newpara  P_eta_l  ,  range:[0.0000,   0.0000], value:[  0.0000,    0.0000], status:fixed     ! Profile eta_linear'
write(IDF_NP, '(a)') 'newpara  P_eta_q  ,  range:[0.0000,   0.0000], value:[  0.0000,    0.0000], status:fixed     ! Profile eta_quad'
write(IDF_NP, '(a)') 'newpara  P_u      ,  range:[0.0000,   0.0000], value:[  0.0000,    0.0000], status:fixed     ! Profile u'
write(IDF_NP, '(a)') 'newpara  P_v      ,  range:[0.0000,   0.0000], value:[  0.0000,    0.0000], status:fixed     ! Profile v'
write(IDF_NP, '(a)') 'newpara  P_w      ,  range:[0.0005,   0.0005], value:[  0.0005,    0.0005], status:fixed     ! Profile w'
write(IDF_NP, '(a)') 'newpara  P_as1_c  ,  range:[0.0000,   0.0000], value:[  0.0000,    0.0000], status:fixed     ! Asymmetry const 1'
write(IDF_NP, '(a)') 'newpara  P_as2_c  ,  range:[0.0000,   0.0000], value:[  0.0000,    0.0000], status:fixed     ! Asymmetry const 2'
write(IDF_NP, '(a)') 'newpara  P_as1_i  ,  range:[0.0000,   0.0000], value:[  0.0000,    0.0000], status:fixed     ! Asymmetry inv   1'
write(IDF_NP, '(a)') 'newpara  P_as2_i  ,  range:[0.0000,   0.0000], value:[  0.0000,    0.0000], status:fixed     ! Asymmetry inv   2'
write(IDF_NP, '(a)') 'newpara  P_as1_l  ,  range:[0.0000,   0.0000], value:[  0.0000,    0.0000], status:fixed     ! Asymmetry lin   1'
write(IDF_NP, '(a)') 'newpara  P_as2_l  ,  range:[0.0000,   0.0000], value:[  0.0000,    0.0000], status:fixed     ! Asymmetry lin   2'
write(IDF_NP, '(a)') 'newpara  P_as1_q  ,  range:[0.0000,   0.0000], value:[  0.0000,    0.0000], status:fixed     ! Asymmetry quad  1'
write(IDF_NP, '(a)') 'newpara  P_as2_q  ,  range:[0.0000,   0.0000], value:[  0.0000,    0.0000], status:fixed     ! Asymmetry quad  2'
if(is_powder) then
   write(IDF_NP, '(a)') '#'
   write(IDF_NP, '(a)') '#newpara  P_portion,  range:[1.000,   1.0000], value:[  1.000,    1.000], status:fixed     ! fixed '
   write(IDF_NP, '(a)') '#newpara  P_damping,  range:[0.000,   0.0000], value:[  0.000,    0.000], status:fixed     ! fixed '
endif
!
end subroutine diffev_main_powder_mac
!
!*****7**************************************************************** 
!
subroutine diffev_setup_last
!-
!  Write the last part of the diffev_setup.mac macro
!
implicit none
!
write(IDF_MA, '(a)') '#'
!write(IDF_MA, '(a)') '@diffev_fix_free.mac         ! ''refine none'' etc '
write(IDF_MA, '(a)') '#'
write(IDF_MA, '(a)') 'diff_cr[1]  = 0.9            ! Cross over probability !'
write(IDF_MA, '(a)') 'diff_f[1]   = 0.81           ! Scale vfactor for the difference vector '
write(IDF_MA, '(a)') 'diff_lo[1]  = 0.0            ! Probability for local search keep at 0.0'
write(IDF_MA, '(a)') 'diff_k[1]   = 1.0            ! Aligmen point of difference vector keep at 1.0'
write(IDF_MA, '(a)') '#'
write(IDF_MA, '(a)') 'donor      random            ! Choose donor randomly among parents'
write(IDF_MA, '(a)') '#selection  compare          ! Choose donor as the best parent'
write(IDF_MA, '(a)') 'selection  best,all          ! New parents are the better of all parents and all children'
write(IDF_MA, '(a)') '#'
write(IDF_MA, '(a)') 'lastfile   DIFFEV/Current    ! Parameter log fo the last generation'
write(IDF_MA, '(a)') 'logfile    DIFFEV/Parameter  ! Parameter log for all generations'
write(IDF_MA, '(a)') 'summary    DIFFEV/Summary    ! Short summary of the fir at each generation'
write(IDF_MA, '(a)') '#'
write(IDF_MA, '(a)') 'backup TEMP/calc, FINAL/final ! Back up files from TEMP/calc to FINAL/final'
!
end subroutine diffev_setup_last
!
!*****7**************************************************************** 
!
subroutine diffev_best_single_mac(is_powder, compute)
!
!  Writes a generic diffev_best part poder/ PDF
!+
use precision_mod
use random_state_mod
!
implicit none
!
logical , intent(in) :: is_powder
character(len=PREC_STRING), intent(in) :: compute
!
character(len=PREC_STRING) :: line
integer :: i, i1
integer :: nseeds
integer, dimension(:), allocatable :: seed_val
!
nseeds = random_nseeds()
allocate(seed_val(1:nseeds))
call random_current(nseeds, seed_val)
!
!
write(IDF_BS, '(a)') 'variable real, P_scale'
write(IDF_BS, '(a)') '#'
write(IDF_BS, '(a,i11)') 'REF_DIMENSION  = ', diffev_dimension
write(IDF_BS, '(a)') '#'
write(IDF_BS, '(a)') '@diffev_best_set_values.mac'
write(IDF_BS, '(a)') '#'
line = 'seed '
i1 = 6
do i=1, nseeds
   i1 = 6 + (i-1)*15
   write(line(i1:i1+16),'(I12,A1)') seed_val(i), ','
enddo
i= len_trim(LINE)
if(line(i:i)==',') line(i:i) = ' '
write(IDF_BS, '(a)') line(1:len_trim(line))
write(IDF_BS, '(a)') '#'
write(IDF_BS, '(a)') 'exit   ! Return to SUITE'
write(IDF_BS, '(a)') '#      ! Each macro on run_mpi command must have an exit to suite as well'
write(IDF_BS, '(a)') '#      ! Each macro is preceeded with a command that steps into the section'
write(IDF_BS, '(a)') '#'
if(compute=='parallel') then
   write(IDF_BS, '(a)') 'variable integer, indiv'
   write(IDF_BS, '(a)') 'do SLOW = 1, REF_NINDIV'
   write(IDF_BS, '(a)') '  discus'
   write(IDF_BS, '(a)') '  indiv     = SLOW'
   write(IDF_BS, '(a)') '  REF_INDIV = SLOW'
   write(IDF_BS, '(a)') '  @dis.diffev.mac'
   write(IDF_BS, '(a)') 'enddo '
   write(IDF_BS, '(a)') '#'
   write(IDF_BS, '(a)') 'kuplot'
   write(IDF_BS, '(a)') '@kup.diffev.mac'
   write(IDF_BS, '(a)') '#'
else
   write(IDF_BS, '(a)') 'discus'
   write(IDF_BS, '(a)') '@dis.diffev.mac'
endif
write(IDF_BS, '(a)') '#'
write(IDF_BS, '(a)') '#'
write(IDF_BS, '(a)') 'set error, continue'
!
deallocate(seed_val)
!
end subroutine diffev_best_single_mac
!
!*****7**************************************************************** 
!
subroutine diffev_best_powder_mac(is_powder, compute)
!
!  Writes a generic diffev_best part poder/ PDF
!+
use precision_mod
use random_state_mod
!
implicit none
!
logical , intent(in) :: is_powder
character(len=PREC_STRING), intent(in) :: compute
!
character(len=PREC_STRING) :: line
integer :: i, i1
integer :: nseeds
integer, dimension(:), allocatable :: seed_val
!
nseeds = random_nseeds()
allocate(seed_val(1:nseeds))
call random_current(nseeds, seed_val)
!
write(IDF_BS, '(a)') 'variable real, P_eta'
write(IDF_BS, '(a)') 'variable real, P_eta_l'
write(IDF_BS, '(a)') 'variable real, P_eta_q'
write(IDF_BS, '(a)') 'variable real, P_u'
write(IDF_BS, '(a)') 'variable real, P_v'
write(IDF_BS, '(a)') 'variable real, P_w'
write(IDF_BS, '(a)') 'variable real, P_as1_c'
write(IDF_BS, '(a)') 'variable real, P_as2_c'
write(IDF_BS, '(a)') 'variable real, P_as1_i'
write(IDF_BS, '(a)') 'variable real, P_as2_i'
write(IDF_BS, '(a)') 'variable real, P_as1_l'
write(IDF_BS, '(a)') 'variable real, P_as2_l'
write(IDF_BS, '(a)') 'variable real, P_as1_q'
write(IDF_BS, '(a)') 'variable real, P_as2_q'
!
if(is_powder) then
   write(IDF_BS, '(a)') '#variable real, P_portion'
   write(IDF_BS, '(a)') '#variable real, P_damping'
endif
write(IDF_BS, '(a)') 'variable real, P_scale'
write(IDF_BS, '(a)') 'variable real, P_zero'
write(IDF_BS, '(a)') '#'
write(IDF_BS, '(a,i11)') 'REF_DIMENSION  = ', diffev_dimension+ 16
write(IDF_BS, '(a)') '#'
write(IDF_BS, '(a)') '@diffev_best_set_values.mac'
write(IDF_BS, '(a)') '#'
line = 'seed '
i1 = 6
do i=1, nseeds
   i1 = 6 + (i-1)*15
   write(line(i1:i1+16),'(I12,A1)') seed_val(i), ','
enddo
i= len_trim(LINE)
if(line(i:i)==',') line(i:i) = ' '
write(IDF_BS, '(a)') line(1:len_trim(line))
write(IDF_BS, '(a)') '#'
write(IDF_BS, '(a)') 'exit   ! Return to SUITE'
write(IDF_BS, '(a)') '#      ! Each macro on run_mpi command must have an exit to suite as well'
write(IDF_BS, '(a)') '#      ! Each macro is preceeded with a command that steps into the section'
write(IDF_BS, '(a)') '#'
if(compute=='parallel') then
   write(IDF_BS, '(a)') 'variable integer, indiv'
   write(IDF_BS, '(a)') 'do SLOW = 1, REF_NINDIV'
   write(IDF_BS, '(a)') 'discus'
   write(IDF_BS, '(a)') '   indiv     = SLOW'
   write(IDF_BS, '(a)') '   REF_INDIV = SLOW'
   write(IDF_BS, '(a)') '   @dis.diffev.mac'
   write(IDF_BS, '(a)') 'enddo '
   write(IDF_BS, '(a)') '#'
   write(IDF_BS, '(a)') 'kuplot'
   write(IDF_BS, '(a)') '@kup.diffev.mac'
   write(IDF_BS, '(a)') '#'
else
   write(IDF_BS, '(a)') 'discus'
   write(IDF_BS, '(a)') '@dis.diffev.mac'
endif
write(IDF_BS, '(a)') '#'
write(IDF_BS, '(a)') '#'
write(IDF_BS, '(a)') 'set error, continue'
!
deallocate(seed_val)
!
end subroutine diffev_best_powder_mac
!
!*****7**************************************************************** 
!
subroutine diffev_best2_powder_mac
!-
!  Write a generic diffev_best_set_values part powder, PDF
!
write(IDF_BS2, '(a)') 'P_eta            =    0.5000'
write(IDF_BS2, '(a)') 'P_eta_l          =    0.0000'
write(IDF_BS2, '(a)') 'P_eta_q          =    0.0000'
write(IDF_BS2, '(a)') 'P_u              =    0.0000'
write(IDF_BS2, '(a)') 'P_v              =    0.0000'
write(IDF_BS2, '(a)') 'P_w              =    0.0005'
write(IDF_BS2, '(a)') 'P_as1_c          =    0.0000'
write(IDF_BS2, '(a)') 'P_as2_c          =    0.0000'
write(IDF_BS2, '(a)') 'P_as1_i          =    0.0000'
write(IDF_BS2, '(a)') 'P_as2_i          =    0.0000'
write(IDF_BS2, '(a)') 'P_as1_l          =    0.0000'
write(IDF_BS2, '(a)') 'P_as2_l          =    0.0000'
write(IDF_BS2, '(a)') 'P_as1_q          =    0.0000'
write(IDF_BS2, '(a)') 'P_as2_q          =    0.0000'
write(IDF_BS2, '(a)') '#P_portion       =    0.0000'
write(IDF_BS2, '(a)') '#P_damping       =    0.0000'
write(IDF_BS2, '(a)') 'P_scale          =    1.0000'
write(IDF_BS2, '(a)') 'P_zero           =    0.0000'
!
end subroutine diffev_best2_powder_mac
!
!*****7**************************************************************** 
!
subroutine write_macro_powder(c_form, hkl_file, rlambda, slambda)
!-
!  Write the powder diffraction macro
!+
use precision_mod
!
implicit none
!
character(len=*)  , intent(in) :: c_form       ! Atom form factor
character(len=*)  , intent(in) :: hkl_file     ! TSC file in case of DISCAMB
real(kind=PREC_DP), intent(in) :: rlambda
character(len=*)  , intent(in) :: slambda      ! Wave length symbol
!
integer, parameter :: IDI = 49
!
open(unit=IDI, file='powder.mac', status='unknown') 
!
write(IDI,'(a )') 'variable real, rlambda'
if(rlambda/=0.0D0) then
   write(IDI,'(a, f8.6)') 'rlambda = ', rlambda
else
   write(IDI,'(a )') 'rlambda = 1.54061'
endif
write(IDI,'(a )') 'powder'
if(c_form == 'discamb') then
   write(IDI,'(3a)') '  xray table:discamb, file:', hkl_file(1:len_trim(hkl_file)), '.tsc'
   write(IDI,'(a )') '  set axis, q'
   write(IDI,'(a )') '  set calc, grid'
else
   write(IDI,'(2a)') '  xray table:', c_form(1:len_trim(c_form))
   write(IDI,'(a )') '  set axis, q'
   write(IDI,'(a )') '  set calc, complete'
endif
write(IDI,'(a,f7.5,a )') '  set qmin, 4.*PI*sind(F_XMIN*0.5)/ rlambda -0.2'
write(IDI,'(a,f7.5,a )') '  set qmax, 4.*PI*sind(F_XMAX*0.5)/ rlambda +0.05'
write(IDI, '(a)') '  set dq  , 0.00025'
write(IDI,'(a )') '  set dh  , 1.0'
write(IDI,'(a )') '  set dk  , 1.0'
write(IDI,'(a )') '  set dl  , 1.0'
write(IDI,'(a )') '  set profile, pseudo'
write(IDI,'(a )') '  set profile, eta, P_eta, P_eta_l, P_eta_q'
write(IDI,'(a )') '  set profile, uvw, P_u, P_v, P_w'
write(IDI,'(a )') '  set profile, asym   , P_as1_c, P_as2_c'
write(IDI,'(a )') '  set profile, asym_i1, P_as1_i, P_as2_i'
write(IDI, '(a)') '# set profile, asym_l , P_as1_l, P_as2_l'
write(IDI, '(a)') '# set profile, asym_q , P_as1_q, P_as2_q'
write(IDI, '(a)') '  set profile, tthmax, 27.0'
write(IDI, '(a)') '# set preferred, march'
write(IDI, '(a)') '# set preferred, portion, P_portion'
write(IDI, '(a)') '# set preferred, damping, P_damping'
write(IDI, '(a)') '# set preferred, hkl, 0, 0, 1'
write(IDI,'(a )') '  set temp, use'
write(IDI,'(a )') '  set disp, off'
if(slambda/=' ') then
   write(IDI,'(2a)') '  set wvle, ', slambda(1:len_trim(slambda))
else
   write(IDI,'(a,f8.6)') '  set wvle, ', rlambda
endif
write(IDI,'(a )') '  set four, four'
write(IDI,'(a )') '  set lpcor, bragg, 0.00'
write(IDI,'(a )') '  set scale, P_scale'
write(IDI,'(a )') '  set tthzero, P_zero'
write(IDI,'(3a)') '  run'
write(IDI,'(a )') 'exit'
!
close(IDI)
!
end subroutine write_macro_powder
!
!*****7**************************************************************** 
!
subroutine write_macro_debye(c_form, hkl_file, rlambda, slambda)
!-
!  Write the powder diffraction macro   debye scattering version
!+
use precision_mod
!
implicit none
!
character(len=*)  , intent(in) :: c_form       ! Atom form factor
character(len=*)  , intent(in) :: hkl_file     ! TSC file in case of DISCAMB
real(kind=PREC_DP), intent(in) :: rlambda
character(len=*)  , intent(in) :: slambda      ! Wave length symbol
!
integer, parameter :: IDI = 49
!
open(unit=IDI, file='debye.mac', status='unknown') 
!
write(IDI,'(a )') 'variable real, rlambda'
if(rlambda/=0.0D0) then
   write(IDI,'(a, f8.6)') 'rlambda = ', rlambda
else
   write(IDI,'(a )') 'rlambda = 0.20000'
endif
write(IDI,'(a )') 'powder'
write(IDI,'(2a)') '  xray table:', c_form(1:len_trim(c_form))
write(IDI,'(a )') '  set axis, q'
write(IDI,'(a )') '  set calc, debye'
write(IDI,'(a,f7.5,a )') '  set qmin, 0.5000'
write(IDI,'(a,f7.5,a )') '  set qmax, 26.0000'
write(IDI, '(a)') '  set dq  , 0.00025'
write(IDI,'(a )') '  set profile, pseudo'
write(IDI,'(a )') '  set profile, eta, P_eta, P_eta_l, P_eta_q'
write(IDI,'(a )') '  set profile, uvw, P_u, P_v, P_w'
write(IDI,'(a )') '  set profile, asym   , P_as1_c, P_as2_c'
write(IDI,'(a )') '  set profile, asym_i1, P_as1_i, P_as2_i'
write(IDI, '(a)') '# set profile, asym_l , P_as1_l, P_as2_l'
write(IDI, '(a)') '# set profile, asym_q , P_as1_q, P_as2_q'
write(IDI, '(a)') '  set profile, tthmax, 27.0'
write(IDI,'(a )') '  set temp, use'
write(IDI,'(a )') '  set disp, off'
if(slambda/=' ') then
   write(IDI,'(2a)') '  set wvle, ', slambda(1:len_trim(slambda))
else
   write(IDI,'(a,f8.6)') '  set wvle, ', rlambda
endif
write(IDI,'(a )') '  set four, four'
write(IDI,'(a )') '  set lpcor, bragg, 0.00'
write(IDI,'(a )') '  set scale, P_scale'
write(IDI,'(a )') '  set tthzero, P_zero'
write(IDI,'(3a)') '  run'
write(IDI,'(a )') 'exit'
!
close(IDI)
!
end subroutine write_macro_debye
!
!*****7**************************************************************** 
!
subroutine write_prepare_hkl_full(substance, ilist, hkl_max, rlambda, slambda)
!-
!  Write the file 'prepare_hkl_full'
!+
!
use lib_conv_shelx_mod
!
implicit none
!
character(len=*)                , intent(in)    :: substance
integer                         , intent(in)    :: ilist
integer           , dimension(3), intent(inout) :: hkl_max
real(kind=PREC_DP)              , intent(in)    :: rlambda
character                       , intent(in)    :: slambda
!
integer :: i
!
if(ilist >  0) then
   i = len_trim(substance)! - 5
   if(lexist_fcf_file) then          ! !An fcf file exists
      call lib_convert_shelx(fcf_file, hkl_file, hkl_max)
   endif
   open(IRE_MA, file='prepare_hkl_full.mac', status='unknown')
   write(IRE_MA,'(a )') 'discus'
   write(IRE_MA,'(a )') '#'
   write(IRE_MA,'(a )') 'variable real, hhmax'
   write(IRE_MA,'(a )') 'variable real, kkmax'
   write(IRE_MA,'(a )') 'variable real, llmax'
   write(IRE_MA,'(a )') 'variable character, substance'
   write(IRE_MA,'(a )') '#'
   write(IRE_MA,'(a,i5)') 'hhmax = ', hkl_max(1)
   write(IRE_MA,'(a,i5)') 'kkmax = ', hkl_max(2)
   write(IRE_MA,'(a,i5)') 'llmax = ', hkl_max(3)
   write(IRE_MA,'(a )') '#'
   write(IRE_MA,'(3a)') 'substance = "%c", ''', substance(1:i), ''''
   write(IRE_MA,'(a )') '#'
   write(IRE_MA,'(a )') 'read'
   write(IRE_MA,'(a )') '  cell "%c.cell", substance'
   write(IRE_MA,'(a )') '#'
   write(IRE_MA,'(a )') 'fourier'
   write(IRE_MA,'(a )') '  xray table:waas'
   write(IRE_MA,'(a )') '  temp use'
   write(IRE_MA,'(a )') '  disp off'
   if(slambda==' ') then
      write(IRE_MA,'(a, f7.5 )') '  wvle ', rlambda
   else
      write(IRE_MA,'(2a      )') '  wvle ', slambda(1:len_trim(slambda))
   endif
   write(IRE_MA,'(a )') '  ll  -hhmax, -kkmax, -llmax'
   write(IRE_MA,'(a )') '  lr   hhmax, -kkmax, -llmax'
   write(IRE_MA,'(a )') '  ul  -hhmax,  kkmax, -llmax'
   write(IRE_MA,'(a )') '  tl  -hhmax, -kkmax,  llmax'
   write(IRE_MA,'(a )') '  na  2*hhmax + 1'
   write(IRE_MA,'(a )') '  no  2*kkmax + 1'
   write(IRE_MA,'(a )') '  nt  2*llmax + 1'
   write(IRE_MA,'(a )') '  set aver, 0.000'
   write(IRE_MA,'(a )') '  run'
   write(IRE_MA,'(a )') '  exit'
   write(IRE_MA,'(a )') '#'
   write(IRE_MA,'(a )') 'output'
   write(IRE_MA,'(3a)') '  outf  "%c_full.hkl", substance'
   write(IRE_MA,'(a )') '  form  hklf4'
   write(IRE_MA,'(a )') '  value inte'
   write(IRE_MA,'(a )') 'run'
   write(IRE_MA,'(a )') '#'
   write(IRE_MA,'(a )') 'exit  ! Back to main DISCUS menu'
   write(IRE_MA,'(a )') 'exit  ! Back to SUITE'
   close(IRE_MA)
endif
!
end subroutine write_prepare_hkl_full
!
!*****7**************************************************************** 
!
end module prep_refine_mod
