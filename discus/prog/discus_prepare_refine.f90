module prep_refine_mod
!
private
public write_refine_single_part1
public write_refine_single_part2
public write_refine_single_part3
public write_refine_single_part4
public write_refine_powder_part1
public write_refine_powder_part4
public write_diffev_generic
public write_diffev_powder_part1
public write_diffev_single_part2
public write_diffev_single_part3
public write_diffev_powder_part4
!
integer, parameter :: IDF_MA  = 38    ! DIFFEV_MAIN PART
integer, parameter :: IDF_FX  = 40    ! DIFFEV_FIX_FREE
integer, parameter :: IDF_SI  = 42    ! DIFFEV_SINGLE
integer, parameter :: IDF_PO  = 43    ! DIFFEV_POWDER
integer, parameter :: IDF_BS  = 44    ! diffev_best main
integer, parameter :: IDF_BS2 = 45    ! diffev_best part 2
integer :: diffev_dimension
!
contains
!
!*****7**************************************************************** 
!
subroutine write_refine_single_part1(IDI, IRE, ofile, infile, ilist, fv_1, &
           discus_file, hkl_file, fcf_file)
!
!  Write main discus macro; single crystal PART 1
!
!
use precision_mod
!
implicit none
!
integer            , intent(in) :: IDI
integer            , intent(in) :: IRE
character(len=*)   , intent(in) :: ofile
character(len=*)   , intent(in) :: infile
integer            , intent(in) :: ilist
real(kind=PREC_DP) , intent(in) :: fv_1
character(len=PREC_STRING), intent(out) :: discus_file
character(len=PREC_STRING), intent(out) ::    hkl_file
character(len=PREC_STRING), intent(out) ::    fcf_file
!
character(len=PREC_STRING) :: refine_file
integer :: i
!
i= len_trim(ofile)-5
discus_file = ofile(1:i)//'_main.mac'
refine_file = 'refine_'//ofile(1:i)//'_merged.mac'
if(ilist == 0) then
   hkl_file = infile(1:len_trim(infile)-4)//'.hkl'
   fcf_file = ' '
else
   hkl_file = infile(1:len_trim(infile)-4)//'_merged.hkl'
   fcf_file = infile(1:len_trim(infile)-4)//'.fcf'
endif
open(unit=IDI, file=discus_file, status='unknown')
open(unit=IRE, file=refine_file, status='unknown')
!
i= len_trim(ofile)
write(IDI,'(a )') 'branch discus'
write(IDI,'(a )') '#'
write(IDI,'(a )') 'read'
write(IDI,'(2a)') '  stru ', ofile(1:i)
!
write(IRE,'(a )') 'refine'
write(IRE,'(a )') 'rese'
write(IRE,'(2a)') 'data hklf4, ',hkl_file(1:len_trim(hkl_file))
write(IRE,'(a,f9.4,a )') 'newpara P_scale, value:', fv_1 , ', points:3, shift:0.001, status:free'
!
end subroutine write_refine_single_part1
!
!*****7**************************************************************** 
!
subroutine write_refine_single_part2(IDI, IRE, j, l, jj, iscat, lcontent, natoms, c_atom, uij_l, c_flag)
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
integer            , intent(in) :: IDI
integer            , intent(in) :: IRE
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
   write(IDI,'(a,i3,4a)') 'anis type:', jj,', values:[', 'U_',c_atom(iscat)(1:len_trim(c_atom(iscat))),'_1]'
   write(IRE,'(4a,g15.8e3,2a)') 'newpara U_',c_atom(iscat)(1:len_trim(c_atom(iscat))),'_1', &
                 ', value:', uij_l(1,iscat, j), ', points:3, shift:0.001, status:', c_flag 
elseif(l==6) then
   write(IDI,'(a,i3,a,5(3a,i1.1,a2),3a)') 'anis type:', jj,', values:[', &
        ('U_',c_atom(iscat)(1:len_trim(c_atom(iscat))),'_',k,', ',k=1,5), &
         'U_',c_atom(iscat)(1:len_trim(c_atom(iscat))),'_6]'
   do m=1,6
      write(IRE,'(3a,i1.1,a,g15.8e3,2a)') 'newpara U_',c_atom(iscat)(1:len_trim(c_atom(iscat))),'_',m, &
                ', value:', uij_l(m,iscat, j), ', points:3, shift:0.001, status:', c_flag 
   enddo
endif
!
end subroutine write_refine_single_part2
!
!*****7**************************************************************** 
!
subroutine write_refine_single_part3(IDI, IRE, line       , iatom, iscat, lcontent, natoms, c_atom, posit, c_flag)
!
!  Write main discus macro; single crystal PART 3
!  can be used for SINGLE and POWDER
!
use precision_mod
!
implicit none
!
integer            , intent(in) :: IDI
integer            , intent(in) :: IRE
character(len=*)   , intent(in) :: line
integer            , intent(in) :: iatom
integer            , intent(in) :: iscat
integer            , intent(in) :: lcontent
integer            , intent(in) :: natoms
character(len=4), dimension(lcontent), intent(in)  :: c_atom    ! Atom types 
real(kind=PREC_DP), dimension(3,natoms), intent(in) :: posit
character(len=*)   , intent(in) :: c_flag
!
integer :: i
!     line = content(jc)(1:4)
i = len_trim(line)
write(IDI,'(a,i3,3a)') 'x[',iatom, '] = P_',line(1:i), '_x'
write(IDI,'(a,i3,3a)') 'y[',iatom, '] = P_',line(1:i), '_y'
write(IDI,'(a,i3,3a)') 'z[',iatom, '] = P_',line(1:i), '_z'
!
write(IRE,'(3a,f12.8,2a)') 'newpara P_',c_atom(iscat)(1:len_trim(c_atom(iscat))),'_x, value:', &
   posit(1, iatom), ', points:3, shift:0.001, status:', c_flag
write(IRE,'(3a,f12.8,2a)') 'newpara P_',c_atom(iscat)(1:len_trim(c_atom(iscat))),'_y, value:', &
   posit(2, iatom), ', points:3, shift:0.001, status:', c_flag
write(IRE,'(3a,f12.8,2a)') 'newpara P_',c_atom(iscat)(1:len_trim(c_atom(iscat))),'_z, value:', &
   posit(3, iatom), ', points:3, shift:0.001, status:', c_flag
!
end subroutine write_refine_single_part3
!
!*****7**************************************************************** 
!
subroutine write_refine_single_part4(IDI, IRE, ofile, discus_file, hkl_file, fcf_file, &
        ilist, c_form, rlambda, P_exti, hkl_max)
!
!  Write main discus macro; single crystal PART 3
!
use blanks_mod
use lib_conv_shelx_mod
use precision_mod
!
implicit none
!
integer                         , intent(in) :: IDI
integer                         , intent(in) :: IRE
character(len=*)                , intent(in) :: ofile
character(len=*)                , intent(in) :: discus_file
character(len=*)                , intent(in) :: hkl_file
character(len=*)                , intent(in) :: fcf_file
integer                         , intent(in) :: ilist
character(len=*)                , intent(in) :: c_form
real(kind=PREC_DP)              , intent(in) :: rlambda
real(kind=PREC_DP)              , intent(in) :: P_exti
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
!
i= len_trim(ofile)
j= len_trim(hkl_file) - 4
write(IDI,'(a )') 'save'
write(IDI,'(2a)') '  outfile internal.',ofile(1:i)
write(IDI,'(a )') '  write all'
write(IDI,'(a )') '  run'
write(IDI,'(a )') 'exit'
write(IDI,'(a )') 'read'
write(IDI,'(2a)') '  cell internal.',ofile(1:i)
write(IDI,'(a )') 'fourier'
if(c_form == 'discamb') then
   write(IDI,'(3a)') '  xray table:discamb, file:', hkl_file(1:j), '.tsc'
else
   write(IDI,'(2a)') '  xray table:', c_form(1:len_trim(c_form))
endif
write(IDI,'(a )') '  temp use'
write(IDI,'(a )') '  disp off'
write(IDI,'(a, f7.5 )') '  wvle ', rlambda
write(IDI,'(a )') '  set aver, 0'
if(P_exti> 0.00001) then
   write(IDI,'(a )') '  set exti:P_exti'
endif
write(IDI,'(a )') '  set technique:turbo'
write(IDI,'(3a)') '  hkl in:',hkl_file(1:len_trim(hkl_file)), ', out:calc.hkl, scale:P_scale, style:hklf4'
write(IDI,'(a )') 'exit'
write(IDI,'(a )') 'branch kuplot'
write(IDI,'(a )') 'reset'
write(IDI,'(a )') 'load hklf4, calc.hkl'
write(IDI,'(a )') 'exit  ! Back to DISCUS'
write(IDI,'(a )') 'exit  ! Back to REFINE'
write(IDI,'(a )') 'finished'
close(IDI)
!
i= len_trim(ofile) - 5
open(IDI, file='k_fobs_fcalc.mac', status='unknown')
write(IDI,'(a )') 'reset'
write(IDI,'(3a)') 'load csv, ',hkl_file(1:len_trim(hkl_file)), ', colx:4, coly:5, separator:[4,4,4,8,8], skip:0'
write(IDI,'(a )') 'load csv, calc.hkl, colx:5, coly:4, separator:[4,4,4,8,8], skip:0'
write(IDI,'(a )') 'ccal mul, wy, 1, 0.0'
write(IDI,'(a )') 'kcal add, 1, 2'
write(IDI,'(a )') 'kfra 1, 3'
write(IDI,'(a )') 'scale 0, max(xmax[3], ymax[3]), 0.0, max(xmax[3], ymax[3])'
write(IDI,'(a )') 'aver 1'
write(IDI,'(a )') 'mark'
write(IDI,'(a )') 'ltyp 3, 0'
write(IDI,'(a )') 'mtyp 3, 3'
write(IDI,'(a )') 'fnam off'
write(IDI,'(a )') 'fset 2'
write(IDI,'(a )') 'grid on'
write(IDI,'(a )') 'achx Iobs'
write(IDI,'(a )') 'achy Icalc'
write(IDI,'(a )') 'mcol 3, black'
write(IDI,'(3a)') 'load csv, ',hkl_file(1:len_trim(hkl_file)), ', colx:5, coly:4, separator:[4,4,4,8,8], skip:0'
write(IDI,'(a )') 'load csv, calc.hkl, colx:5, coly:4, separator:[4,4,4,8,8], skip:0'
write(IDI,'(a )') 'do LOOP = 1, np[4]'
write(IDI,'(a )') '  dy[4,LOOP] = x[4,LOOP]'
write(IDI,'(a )') 'enddo'
write(IDI,'(a )') 'rval 4,5, dat'
write(IDI,'(2a)') 'tit1 Substance: ', ofile(1:i)
write(IDI,'(a )') 'tit2 "wR-value %8.5f",res[2]'
write(IDI,'(a )') 'plot'
write(IDI,'(a )') 'exit '
close(IDI)
!
if(P_exti> 0.00001) then
   write(IRE,'(a,f8.5,a)') 'newpara P_exti, value:', P_exti, ', points:3, shift:0.001, status:free'
else
   write(IRE,'(a,f8.5,a)') 'newpara P_exti, value:', P_exti, ', points:3, shift:0.001, status:fixed'
endif
write(IRE,'(a )') 'set cycle,   5'
write(IRE,'(a )') 'set conver, status:on, dchi:0.050, chisq:1.10, pshift:2.0, conf:1.0, lambda:65000.'
write(IRE,'(2a)') '@', discus_file(1:len_trim(discus_file))
write(IRE,'( a)') 'branch kuplot'
write(IRE,'( a)') '  @k_fobs_fcalc.mac'
write(IRE,'(3a)') 'run ', discus_file(1:len_trim(discus_file)), ', plot:k_fobs_fcalc.mac'
write(IRE,'(a )') '#'
write(IRE,'(a )') '@final_cell_cif.mac'
write(IRE,'(a )') '#'
write(IRE,'(a )') 'exit  ! Back to SUITE'
close(IRE)
!
open(IRE, file='final_cell_cif.mac', status='unknown')
write(IRE,'( a)') 'branch discus'
write(IRE,'( a)') 'read'
write(IRE,'(3a)') '  stru internal.' , ofile(1:i), '.cell'
write(IRE,'( a)') 'save'
write(IRE,'(3a)') '  outf ', ofile(1:i), '.cell'
write(IRE,'( a)') '  write all'
write(IRE,'( a)') '  run'
write(IRE,'( a)') 'exit'
write(IRE,'(3a)') 'export cif, ', ofile(1:i), '.cif, spcgr:original'
write(IRE,'( a)') 'exit  ! back to REFINE'
close(IRE)
!
!
if(ilist >  0) then
   i = len_trim(ofile) - 5
   call lib_convert_shelx(fcf_file, hkl_file, hkl_max)
   open(IRE, file='prepare_hkl_full.mac', status='unknown')
   write(IRE,'(a )') 'discus'
   write(IRE,'(a )') '#'
   write(IRE,'(a )') 'variable real, hhmax'
   write(IRE,'(a )') 'variable real, kkmax'
   write(IRE,'(a )') 'variable real, llmax'
   write(IRE,'(a )') 'variable character, substance'
   write(IRE,'(a )') '#'
   write(IRE,'(a,i5)') 'hhmax = ', hkl_max(1)
   write(IRE,'(a,i5)') 'kkmax = ', hkl_max(2)
   write(IRE,'(a,i5)') 'llmax = ', hkl_max(3)
   write(IRE,'(a )') '#'
   write(IRE,'(3a)') 'substance = "%c", ''', ofile(1:i), ''''
   write(IRE,'(a )') '#'
   write(IRE,'(a )') 'read'
   write(IRE,'(a )') '  cell "%c.cell", substance'
   write(IRE,'(a )') '#'
   write(IRE,'(a )') 'fourier'
   write(IRE,'(a )') '  xray table:waas'
   write(IRE,'(a )') '  temp use'
   write(IRE,'(a )') '  disp off'
   write(IRE,'(a, f7.5 )') '  wvle ', rlambda
   write(IRE,'(a )') '  ll  -hhmax, -kkmax, -llmax'
   write(IRE,'(a )') '  lr   hhmax, -kkmax, -llmax'
   write(IRE,'(a )') '  ul  -hhmax,  kkmax, -llmax'
   write(IRE,'(a )') '  tl  -hhmax, -kkmax,  llmax'
   write(IRE,'(a )') '  na  2*hhmax + 1'
   write(IRE,'(a )') '  no  2*kkmax + 1'
   write(IRE,'(a )') '  nt  2*llmax + 1'
   write(IRE,'(a )') '  set aver, 0.000'
   write(IRE,'(a )') '  run'
   write(IRE,'(a )') '  exit'
   write(IRE,'(a )') '#'
   write(IRE,'(a )') 'output'
   write(IRE,'(3a)') '  outf  "%c_full.hkl", substance'
   write(IRE,'(a )') '  form  hklf4'
   write(IRE,'(a )') '  value inte'
   write(IRE,'(a )') 'run'
   write(IRE,'(a )') '#'
   write(IRE,'(a )') 'exit  ! Back to main DISCUS menu'
   write(IRE,'(a )') 'exit  ! Back to SUITE'
   close(IRE)
endif
!
!  Interpret "aspher.json"
!
if(c_form=='discamb') then
   aspher_file = '../BUILD_TSC/aspher.json'  ! Use the special fixed name for aspher.json
   inquire(file=aspher_file, exist = lexist)
   cond_exist: if(lexist) then
      open(unit=IRE, file=aspher_file, status='old')
      read(IRE, '(a)', iostat=ios) line
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
         read(IRE, '(a)', iostat=ios) line
         if(is_iostat_end(ios)) exit cond_exist
         if(ios/= 0           ) exit cond_exist
         if(line=='}'         ) exit cond_exist
         length=len_trim(line)
         call rem_leading_bl(line, length)
      enddo loop_aspher
   endif cond_exist
   close(IRE)
endif
!
end subroutine write_refine_single_part4
!
!*****7**************************************************************** 
!
subroutine write_refine_powder_part1(IDI, IRE, substance, infile, ilist, P_scale, &
           spcgr_syst, lattice_para, discus_file, hkl_file)
!
!  Write main discus macro; powder diffraction PART 1
!
!
use precision_mod
!
implicit none
!
integer            , intent(in) :: IDI
integer            , intent(in) :: IRE
character(len=*)   , intent(in) :: substance
character(len=*)   , intent(in) :: infile
integer            , intent(in) :: ilist
real(kind=PREC_DP) , intent(in) :: P_scale
integer            , intent(in) :: spcgr_syst
real(kind=PREC_DP) , dimension(6), intent(in) :: lattice_para
character(len=PREC_STRING), intent(out) :: discus_file
character(len=PREC_STRING), intent(out) ::    hkl_file
!character(len=PREC_STRING), intent(out) ::    fcf_file
!
character(len=PREC_STRING) :: refine_file
integer :: i
!
i= len_trim(substance)
discus_file = substance(1:i)//'_main_powder.mac'
refine_file = 'refine_'//substance(1:i)//'_powder.mac'
if(ilist == 0) then
   hkl_file = infile(1:len_trim(infile)-4)//'.tth'
else
   hkl_file = infile(1:len_trim(infile)-4)//'.Q'
endif
open(unit=IDI, file=discus_file, status='unknown')
open(unit=IRE, file=refine_file, status='unknown')
!
i= len_trim(substance)
write(IDI,'(a )') 'branch discus'
write(IDI,'(a )') '#'
write(IDI,'(a )') 'read'
write(IDI,'(3a)') '  stru CELL/', substance(1:i), '.cell'
!
write(IRE,'(a )') 'refine'
write(IRE,'(a )') 'rese'
write(IRE,'(2a)') 'data xy, DATA/',hkl_file(1:len_trim(hkl_file))
write(IRE,'(a,f9.4,a )') 'newpara P_scale, value:', P_scale , ', points:3, shift:0.001, status:free'
write(IRE,'(a,f9.4,a )') 'newpara P_zero , value:', 0.00000 , ', points:3, shift:0.001, status:free'
write(IRE,'(a )') '#'
write(IRE,'(a,f9.4,a )') 'newpara P_eta  , value:', 0.50000 , ', points:3, shift:0.001, status:fixed'
write(IRE,'(a,f9.4,a )') 'newpara P_eta_l, value:', 0.00000 , ', points:3, shift:0.001, status:fixed'
write(IRE,'(a,f9.4,a )') 'newpara P_eta_q, value:', 0.00000 , ', points:3, shift:0.001, status:fixed'
write(IRE,'(a,f9.4,a )') 'newpara P_u    , value:', 0.00000 , ', points:3, shift:0.001, status:fixed'
write(IRE,'(a,f9.4,a )') 'newpara P_v    , value:', 0.00000 , ', points:3, shift:0.001, status:fixed'
write(IRE,'(a,f9.4,a )') 'newpara P_w    , value:', 0.00700 , ', points:3, shift:0.001, status:fixed'
write(IRE,'(a )') '#'
write(IRE,'(a,f9.4,a )') 'newpara P_as1_c, value:', 0.00000 , ', points:3, shift:0.001, status:fixed'
write(IRE,'(a,f9.4,a )') 'newpara P_as2_c, value:', 0.00000 , ', points:3, shift:0.001, status:fixed'
write(IRE,'(a,f9.4,a )') 'newpara P_as1_i, value:', 0.00000 , ', points:3, shift:0.001, status:fixed'
write(IRE,'(a,f9.4,a )') 'newpara P_as2_i, value:', 0.00000 , ', points:3, shift:0.001, status:fixed'
write(IRE,'(a,f9.4,a )') 'newpara P_as1_l, value:', 0.00000 , ', points:3, shift:0.001, status:fixed'
write(IRE,'(a,f9.4,a )') 'newpara P_as2_l, value:', 0.00000 , ', points:3, shift:0.001, status:fixed'
write(IRE,'(a,f9.4,a )') 'newpara P_as1_q, value:', 0.00000 , ', points:3, shift:0.001, status:fixed'
write(IRE,'(a,f9.4,a )') 'newpara P_as2_q, value:', 0.00000 , ', points:3, shift:0.001, status:fixed'
write(IRE,'(a )') '#'
!
if(spcgr_syst==1) then         ! Triclinic
   write(IRE,'(a,f9.4,a )') 'newpara P_alat , value:', lattice_para(1) , ', points:3, shift:0.001, status:free'
   write(IRE,'(a,f9.4,a )') 'newpara P_blat , value:', lattice_para(2) , ', points:3, shift:0.001, status:free'
   write(IRE,'(a,f9.4,a )') 'newpara P_clat , value:', lattice_para(3) , ', points:3, shift:0.001, status:free'
   write(IRE,'(a,f9.4,a )') 'newpara P_alpha, value:', lattice_para(4) , ', points:3, shift:0.001, status:free'
   write(IRE,'(a,f9.4,a )') 'newpara P_beta , value:', lattice_para(5) , ', points:3, shift:0.001, status:free'
   write(IRE,'(a,f9.4,a )') 'newpara P_gamma, value:', lattice_para(6) , ', points:3, shift:0.001, status:free'
elseif(spcgr_syst==2) then         ! Monoclinic-B
   write(IRE,'(a,f9.4,a )') 'newpara P_alat , value:', lattice_para(1) , ', points:3, shift:0.001, status:free'
   write(IRE,'(a,f9.4,a )') 'newpara P_blat , value:', lattice_para(2) , ', points:3, shift:0.001, status:free'
   write(IRE,'(a,f9.4,a )') 'newpara P_clat , value:', lattice_para(3) , ', points:3, shift:0.001, status:free'
   write(IRE,'(a,f9.4,a )') 'newpara P_beta , value:', lattice_para(5) , ', points:3, shift:0.001, status:free'
elseif(spcgr_syst==3) then         ! Monoclinic-C
   write(IRE,'(a,f9.4,a )') 'newpara P_alat , value:', lattice_para(1) , ', points:3, shift:0.001, status:free'
   write(IRE,'(a,f9.4,a )') 'newpara P_blat , value:', lattice_para(2) , ', points:3, shift:0.001, status:free'
   write(IRE,'(a,f9.4,a )') 'newpara P_clat , value:', lattice_para(3) , ', points:3, shift:0.001, status:free'
   write(IRE,'(a,f9.4,a )') 'newpara P_gamma, value:', lattice_para(5) , ', points:3, shift:0.001, status:free'
elseif(spcgr_syst==4) then         ! Orthorhombic
   write(IRE,'(a,f9.4,a )') 'newpara P_alat , value:', lattice_para(1) , ', points:3, shift:0.001, status:free'
   write(IRE,'(a,f9.4,a )') 'newpara P_blat , value:', lattice_para(2) , ', points:3, shift:0.001, status:free'
   write(IRE,'(a,f9.4,a )') 'newpara P_clat , value:', lattice_para(3) , ', points:3, shift:0.001, status:free'
elseif(spcgr_syst==5) then         ! tetragonal  
   write(IRE,'(a,f9.4,a )') 'newpara P_alat , value:', lattice_para(1) , ', points:3, shift:0.001, status:free'
   write(IRE,'(a,f9.4,a )') 'newpara P_clat , value:', lattice_para(3) , ', points:3, shift:0.001, status:free'
elseif(spcgr_syst==6) then         ! trigonal  
   write(IRE,'(a,f9.4,a )') 'newpara P_alat , value:', lattice_para(1) , ', points:3, shift:0.001, status:free'
   write(IRE,'(a,f9.4,a )') 'newpara P_clat , value:', lattice_para(3) , ', points:3, shift:0.001, status:free'
elseif(spcgr_syst==7) then         ! trigonal  
   write(IRE,'(a,f9.4,a )') 'newpara P_alat , value:', lattice_para(1) , ', points:3, shift:0.001, status:free'
   write(IRE,'(a,f9.4,a )') 'newpara P_alpha, value:', lattice_para(4) , ', points:3, shift:0.001, status:free'
elseif(spcgr_syst==8) then         ! hexagonal  
   write(IRE,'(a,f9.4,a )') 'newpara P_alat , value:', lattice_para(1) , ', points:3, shift:0.001, status:free'
   write(IRE,'(a,f9.4,a )') 'newpara P_clat , value:', lattice_para(3) , ', points:3, shift:0.001, status:free'
elseif(spcgr_syst==9) then         ! cubic      
   write(IRE,'(a,f9.4,a )') 'newpara P_alat , value:', lattice_para(1) , ', points:3, shift:0.001, status:free'
endif
!
end subroutine write_refine_powder_part1
!
!*****7**************************************************************** 
!
subroutine write_refine_powder_part4(IDI, IRE, substance, discus_file, hkl_file, fcf_file, &
        ilist, c_form, rlambda, P_exti, spcgr_syst, lattice_para, hkl_max)
!
!  Write main discus macro; single crystal PART 3
!
use blanks_mod
use lib_conv_shelx_mod
use precision_mod
!
implicit none
!
integer                         , intent(in) :: IDI
integer                         , intent(in) :: IRE
character(len=*)                , intent(in) :: substance
character(len=*)                , intent(in) :: discus_file
character(len=*)                , intent(in) :: hkl_file
character(len=*)                , intent(in) :: fcf_file
integer                         , intent(in) :: ilist
character(len=*)                , intent(in) :: c_form
real(kind=PREC_DP)              , intent(in) :: rlambda
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
!
i= len_trim(substance)
j= len_trim(hkl_file) - 4
!
call write_powder_lattice(IDI, spcgr_syst)
!
write(IDI,'(a )') 'save'
write(IDI,'(3a)') '  outfile internal.',substance(1:i), '.cell'
write(IDI,'(a )') '  write all'
write(IDI,'(a )') '  run'
write(IDI,'(a )') 'exit'
write(IDI,'(a )') 'read'
write(IDI,'(3a)') '  cell internal.',substance(1:i), '.cell'
write(IDI,'(a )') 'variable real, rlambda'
write(IDI,'(a )') 'rlambda = 1.5409'
write(IDI,'(a )') '@powder.mac'
write(IDI,'(a )') 'powder'
if(c_form == 'discamb') then
   write(IDI,'(3a)') '  xray table:discamb, file:', hkl_file(1:j), '.tsc'
   write(IDI,'(a )') '  set axis, q'
   write(IDI,'(a )') '  set calc, grid'
else
   write(IDI,'(2a)') '  xray table:', c_form(1:len_trim(c_form))
   write(IDI,'(a )') '  set axis, q'
   write(IDI,'(a )') '  set calc, complete'
endif
write(IDI,'(a,f7.5,a )') '  set qmin, 4.*PI*sind(F_XMIN*0.5)/ rlambda -0.2'
write(IDI,'(a,f7.5,a )') '  set qmax, 4.*PI*sind(F_XMAX*0.5)/ rlambda +0.05'
write(IDI,'(a )') '  set dq  , 0.005'
write(IDI,'(a )') '  set dh  , 1.0'
write(IDI,'(a )') '  set dk  , 1.0'
write(IDI,'(a )') '  set dl  , 1.0'
write(IDI,'(a )') '  set profile, pseudo'
write(IDI,'(a )') '  set profile, eta, P_eta, P_eta_l, P_eta_q'
write(IDI,'(a )') '  set profile, uvw, P_u, P_v, P_w'
write(IDI,'(a )') '  set profile, asym   , P_as1_c, P_as2_c'
write(IDI,'(a )') '  set profile, asym_i1, P_as1_i, P_as2_i'
write(IDI,'(a )') '  set profile, asym_l , P_as1_l, P_as2_l'
write(IDI,'(a )') '  set profile, asym_q , P_as1_q, P_as2_q'
write(IDI,'(a )') '  set temp, use'
write(IDI,'(a )') '  set disp, off'
write(IDI,'(a )') '  set wvle, CU12'
write(IDI,'(a )') '  set four, four'
write(IDI,'(a )') '  set lpcor, bragg, 0.00'
write(IDI,'(a )') '  set scale, P_scale'
write(IDI,'(a )') '  set tthzero, P_zero'
write(IDI,'(3a)') '  run'
write(IDI,'(a )') 'exit'
write(IDI,'(a )') 'output'
write(IDI,'(3a)') '  outf POWDER/', substance(1:i), '.raw'
write(IDI,'(a )') '  value inte'
write(IDI,'(a )') '  form powder, tth, F_XMIN, F_XMAX, F_XSTP'
write(IDI,'(a )') '  run'
write(IDI,'(a )') 'exit'
write(IDI,'(a )') 'branch kuplot'
write(IDI,'(a )') 'reset'
write(IDI,'(3a)') 'load xy, POWDER/', substance(1:i), '.raw'
write(IDI,'(3a)') 'load xy, DATA/', substance(1:i), '.tth'
write(IDI,'(a )') 'fit 2'
write(IDI,'(a )') '  func back, 1, 5'
write(IDI,'(a )') '  para 1, 0, 1.0'
write(IDI,'(a )') '  para 2, 1, ymin[2]'
write(IDI,'(a )') '  para 3, 1, 0.0'
write(IDI,'(a )') '  para 4, 1, 0.0'
write(IDI,'(a )') '  para 5, 1, 0.0'
write(IDI,'(a )') '  cycle 10'
write(IDI,'(a )') '  run'
write(IDI,'(a )') 'exit'
write(IDI,'(a )') 'ksav n[1]-1'
write(IDI,'(3a)') '  outf POWDER/', substance(1:i),'.tth'
write(IDI,'(a )') '  form xy'
write(IDI,'(a )') '  run'
write(IDI,'(a )') 'n[1] = n[1]-1'
write(IDI,'(a )') 'exit  ! Back to DISCUS'
write(IDI,'(a )') 'exit  ! Back to REFINE'
write(IDI,'(a )') 'finished'
close(IDI)
!
i= len_trim(substance)
open(IDI, file='k_inter.mac', status='unknown')
write(IDI,'(a )') 'reset'
write(IDI,'(3a)') 'load xy, DATA/', substance(1:i), '.tth'
write(IDI,'(3a)') 'load xy, POWDER/', substance(1:i), '.tth'
write(IDI,'(a )') 'kcal sub, 1, 2'
write(IDI,'(a )') 'r[0] = min(ymin[1], ymin[2]) - ymax[3]'
write(IDI,'(a )') 'ccal add, wy, 3, r[0]'
write(IDI,'(a )') 'aver'
write(IDI,'(a )') 'scale'
write(IDI,'(a )') 'mark'
write(IDI,'(a )') 'ltyp 1, 1'
write(IDI,'(a )') 'ltyp 2, 1'
write(IDI,'(a )') 'ltyp 3, 1'
write(IDI,'(a )') 'lcol 1, blue'
write(IDI,'(a )') 'lcol 2, red'
write(IDI,'(a )') 'lcol 3, black'
write(IDI,'(a )') 'mtyp 1, 0'
write(IDI,'(a )') 'mtyp 2, 0'
write(IDI,'(a )') 'mtyp 3, 0'
write(IDI,'(a )') 'fnam off'
write(IDI,'(a )') 'fset 2'
write(IDI,'(a )') 'grid on'
write(IDI,'(a,f9.7 )') 'achx 2Theta @ ', rlambda
write(IDI,'(a )') 'achy Intensity'
write(IDI,'(a )') 'rval 1,2, dat'
write(IDI,'(2a)') 'tit1 Substance: ', substance(1:i)
write(IDI,'(a )') 'tit2 "wR-value %8.5f",res[2]'
write(IDI,'(a )') 'plot'
write(IDI,'(a )') 'exit '
close(IDI)
!
write(IRE,'(a )') 'set cycle,   5'
write(IRE,'(a )') 'set conver, status:on, dchi:0.050, chisq:1.10, pshift:2.0, conf:1.0, lambda:65000.'
write(IRE,'(a )') 'set relax, start:0.25'
write(IRE,'(2a)') '@', discus_file(1:len_trim(discus_file))
write(IRE,'( a)') 'branch kuplot'
write(IRE,'( a)') '  @k_inter.mac'
write(IRE,'(3a)') 'run ', discus_file(1:len_trim(discus_file)), ', plot:k_inter.mac'
write(IRE,'(a )') '#'
write(IRE,'(a )') '@final_cell_cif.mac'
write(IRE,'(a )') '#'
write(IRE,'(a )') 'exit  ! Back to SUITE'
close(IRE)
!
open(IRE, file='final_cell_cif.mac', status='unknown')
write(IRE,'( a)') 'branch discus'
write(IRE,'( a)') 'read'
write(IRE,'(3a)') '  stru internal.' , substance(1:i), '.cell'
write(IRE,'( a)') 'save'
write(IRE,'(3a)') '  outf ', substance(1:i), '.cell'
write(IRE,'( a)') '  write all'
write(IRE,'( a)') '  run'
write(IRE,'( a)') 'exit'
write(IRE,'(3a)') 'export cif, ', substance(1:i), '.cif, spcgr:original'
write(IRE,'( a)') 'exit  ! back to REFINE'
close(IRE)
!
!  Interpret "aspher.json"
!
if(c_form=='discamb') then
   aspher_file = '../BUILD_TSC/aspher.json'  ! Use the special fixed name for aspher.json
   inquire(file=aspher_file, exist = lexist)
   cond_exist: if(lexist) then
      open(unit=IRE, file=aspher_file, status='old')
      read(IRE, '(a)', iostat=ios) line
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
         read(IRE, '(a)', iostat=ios) line
         if(is_iostat_end(ios)) exit cond_exist
         if(ios/= 0           ) exit cond_exist
         if(line=='}'         ) exit cond_exist
         length=len_trim(line)
         call rem_leading_bl(line, length)
      enddo loop_aspher
   endif cond_exist
   close(IRE)
endif
!
end subroutine write_refine_powder_part4
!
!*****7**************************************************************** 
!
subroutine write_powder_lattice(IDI, spcgr_syst)
!-
! Write the attice parameters into the discus_main file 
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
subroutine write_diffev_generic(substance )
!-
!  Writes the general part for a DIFFEV refinement
!+
!
use precision_mod
!
implicit none
!
character(len=*), intent(in) :: substance
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
outfile = 'diffev_main_' // substance(1:len_trim(substance)) // '.mac'
open(IDF_SI, file=outfile, status='unknown')
write(IDF_SI, '(a)') 'variable integer, diff_counter'
write(IDF_SI, '(a)') 'diffev'
write(IDF_SI, '(a)') '  reset'
write(IDF_SI, '(a)') '  @cleanup.mac       ! Remove old files'
write(IDF_SI, '(a)') '  @setup_rep.mac     ! Define number of individual repetitions'
write(IDF_SI, '(a)') '  @global.mac        ! Define global file and directory names'
write(IDF_SI, '(a)') '  @diffev_setup.mac  ! Define refinement details'
write(IDF_SI, '(a)') '  init               ! Initialize all parameters'
write(IDF_SI, '(a)') '  do diff_counter = 1, 10  ! Perform 10 loops!'
write(IDF_SI, '(a)') '    echo "In loop %4d ", diff_counter     ! Keep user posted'
write(IDF_SI, '(a)') '    run_mpi discus, dis.diffev.mac  , repeat:REF_NINDIV, compute:serial, logfile:LOGFILES/d'
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
write(IDF_SI, '(a)') '  @setup_rep.mac     ! Define number of individual repetitions'
write(IDF_SI, '(a)') '  @diffev_setup.mac  ! Define refinement details'
write(IDF_SI, '(a)') '  do diff_counter = 1, 10  ! Perform 10 loops!'
write(IDF_SI, '(a)') '    echo "In loop %4d ", diff_counter     ! Keep user posted'
write(IDF_SI, '(a)') '    run_mpi discus, dis.diffev.mac  , repeat:REF_NINDIV, compute:serial, logfile:LOGFILES/d'
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
write(IDF_SI, '(a)') 'system rm -f GRCALC/*'
write(IDF_SI, '(a)') 'system rm -f INDI/*'
write(IDF_SI, '(a)') 'system rm -f LOGFILES/*'
write(IDF_SI, '(a)') 'system rm -f PLOT/*'
write(IDF_SI, '(a)') 'system rm -f POWDER/*'
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
write(IDF_SI, '(a)') 'REF_NINDIV = 1      ! No individual repetitions'
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
write(IDF_SI, '(a)') '#@ KEYWORD      refinemetn, generic, nanoparticle'
write(IDF_SI, '(a)') '#@'
write(IDF_SI, '(a)') '#@ DESCRIPTION  Main macro for the refinement of ellipsoidal'
write(IDF_SI, '(a)') '#@ DESCRIPTION  ZnSe nanoparticles with growth type stacking faults.'
write(IDF_SI, '(a)') '#@ DESCRIPTION  In this set of macros the structure is refined versus'
write(IDF_SI, '(a)') '#@ DESCRIPTION  the PDF'
write(IDF_SI, '(a)') '#@ DESCRIPTION'
write(IDF_SI, '(a)') '#@ DESCRIPTION  As this macro is run as redirected input file on the'
write(IDF_SI, '(a)') '#@ DESCRIPTION  call to DISCUS, the first line MUST be "set prompt,redirect"'
write(IDF_SI, '(a)') ''
write(IDF_SI, '(a)') '#@'
write(IDF_SI, '(a)') '#@ PARAMETER    $0, 1'
write(IDF_SI, '(a)') '#@'
write(IDF_SI, '(a)') '#@ USAGE        @dis.diffev.mac <cwd>, <kid>'
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
write(IDF_SI, '(a)') 'do indiv = 1, REF_NINDIV'
write(IDF_SI, '(3a)') '  @',substance(1:len_trim(substance)), '_main_powder.mac REF_KID'
write(IDF_SI, '(a)') 'enddo'
write(IDF_SI, '(a)') 'branch   kuplot    !Switch to KUPLOT'
write(IDF_SI, '(a)') '   @kup.diffev.mac ., REF_KID  ! contains a final exit'
write(IDF_SI, '(a)') '#'
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
write(IDF_MA, '(a)') 'pop_n[1]    = 102    ! Define the population size should be 10 times the number of parameters'
write(IDF_MA, '(a)') 'pop_c[1]    = 102    ! Number of children in each generation usually identical to parent'
write(IDF_MA, '(a)') '#'
write(IDF_MA, '(a)') '#  For each parameter to be refined we define:'
write(IDF_MA, '(a)') '#     name, minimum/maximum allowed values , minimum/maximum for start range'
write(IDF_MA, '(a)') '#'
write(IDF_MA, '(a)') '#  '
write(IDF_MA, '(a)') 'newpara  P_scale  ,  0.001,  20.000,  0.9000, 1.100'
write(IDF_MA, '(a)') 'newpara  P_zero   ,  0.000,   0.000,  0.0000, 0.000   ! fixed'
!
!
! Create a 'diffev_fix_free macro
!
outfile = 'diffev_fix_free.mac'
open(IDF_FX, file=outfile, status='unknown')
write(IDF_FX, '(a)') '#  '
write(IDF_FX, '(a)') 'refine none  '
write(IDF_FX, '(a)') 'refine P_scale  '
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
write(IDF_BS, '(a)') 'discus'
write(IDF_BS, '(a)') 'reset'
write(IDF_BS, '(a)') 'exit'
write(IDF_BS, '(a)') 'kuplot'
write(IDF_BS, '(a)') 'reset'
write(IDF_BS, '(a)') 'exit'
write(IDF_BS, '(a)') '#'
write(IDF_BS, '(a)') 'discus'
write(IDF_BS, '(a)') 'REF_GENERATION =            0'
write(IDF_BS, '(a)') 'REF_MEMBER     =          102'
write(IDF_BS, '(a)') 'REF_CHILDREN   =          102'
write(IDF_BS, '(a)') 'REF_DIMENSION  =           33'
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
subroutine write_diffev_global(substance, datafile, user_name)
!-
!  Writes the first part for a DIFFEV powder refinement
!+
!
use precision_mod
!
implicit none
!
character(len=*), intent(in) :: substance
character(len=*), intent(in) :: datafile
character(len=*), intent(in) :: user_name
!
character(len=PREC_STRING) :: line     ! Dummy line
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
write(IDF_PO, '(a)') 'variable character, INDIDIR   ! temporary PDF directory  ''kuplot''       or ''TEMP'''
write(IDF_PO, '(a)') 'variable character, DATADIR   ! input data directory     ''/tmp/user_name/'' or. ''.'''
write(IDF_PO, '(a)') 'variable character, DATAFILE  ! Input Data File  '
write(IDF_PO, '(a)') 'variable character, SUBSTANCE ! Substance name  '
write(IDF_PO, '(a)') '#'
write(IDF_PO, '(a)') '# LOCAL Version, serial'
write(IDF_PO, '(a)') '#'
write(IDF_PO, '(a)') '#DATADIR   = "%c",''.'''
write(IDF_PO, '(3a)') '#DATAFILE  = "%c",''',datafile(1:len_trim(datafile)),''''
write(IDF_PO, '(a)') '#INDIDIR   = "%c",''kuplot'''
write(IDF_PO, '(a)') '#TMPDIR    = "%c",''internal'''
write(IDF_PO, '(3a)') '#SUBSTANCE = "%c",''', substance(1:len_trim(substance)),''''
write(IDF_PO, '(a)') '#'
write(IDF_PO, '(a)') '# LOCAL Version, parallel'
write(IDF_PO, '(a)') '#'
write(IDF_PO, '(a)') 'DATADIR   = "%c",''.'''
write(IDF_PO, '(3a)') 'DATAFILE  = "%c",''',datafile(1:len_trim(datafile)),''''
write(IDF_PO, '(a)') 'INDIDIR   = "%c",''.'''
write(IDF_PO, '(a)') 'TMPDIR    = "%c",''.'''
write(IDF_PO, '(3a)') 'SUBSTANCE = "%c",''', substance(1:len_trim(substance)),''''
write(IDF_PO, '(a)') '#'
write(IDF_PO, '(a)') '# High performace compute center individuals in parallel'
write(IDF_PO, '(a)') '#'
write(IDF_PO, '(3a)') '#DATADIR   = "%c",''/tmp/', user_name(1:len_trim(user_name)),'''      ! Make sure these are unique'
write(IDF_PO, '(3a)') '#DATAFILE  = "%c",''',datafile(1:len_trim(datafile)),''''
write(IDF_PO, '(3a)') '#INDIDIR   = "%c",''/tmp/', user_name(1:len_trim(user_name)),'''      ! Make sure these are unique'
write(IDF_PO, '(a)') '#TMPDIR    = "%c",''internal'''
write(IDF_PO, '(3a)') '#SUBSTANCE = "%c",''', substance(1:len_trim(substance)),''''
write(IDF_PO, '(a)') '#'
write(IDF_PO, '(a)') '# High performace compute center individuals serially'
write(IDF_PO, '(a)') '#'
write(IDF_PO, '(3a)') '#DATADIR   = "%c",''/tmp/', user_name(1:len_trim(user_name)),'''      ! Make sure these are unique'
write(IDF_PO, '(3a)') '#DATAFILE  = "%c",''',datafile(1:len_trim(datafile)),''''
write(IDF_PO, '(a)') '#INDIDIR   = "%c",''kuplot'''
write(IDF_PO, '(a)') '#TMPDIR    = "%c",''internal'''
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
subroutine write_diffev_powder_part1(substance, user_name,      &
           spcgr_syst, lattice_para)
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
!
character(len=PREC_STRING) :: line     ! Dummy line
character(len=PREC_STRING) :: outfile     ! output file name(s)
!
line = 'mkdir -p POWDER'
call execute_command_line(line, wait=.false.)
!
outfile = substance(1:len_trim(substance)) // '.tth'
call write_diffev_global(substance, outfile, user_name)
!
! Start main discus macro
!
outfile = substance(1:len_trim(substance)) // '_main_powder.mac'
open(IDF_PO, file=outfile, status='unknown')
!
write(IDF_PO, '(a)') '#'
write(IDF_PO, '(a)') 'variable real, rlambda'
write(IDF_PO, '(a)') 'variable real, F_XMIN'
write(IDF_PO, '(a)') 'variable real, F_XMAX'
write(IDF_PO, '(a)') 'variable real, F_XSTP'
write(IDF_PO, '(a)') 'read'
write(IDF_PO, '(3a)') '  stru CELL/', substance(1:len_trim(substance)), '.cell'
write(IDF_PO, '(a)') '#'
!
call write_powder_lattice(IDF_PO, spcgr_syst)
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
end subroutine write_diffev_powder_part1
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
write(IDF_MA,'(a9, a7,4(a,f9.4)  )') 'newpara  ', p_name,      &
   ', ', l_low,    ', ', l_hig, ', ', lp*0.995, ', ', lp*1.00
!
write(IDF_FX, '(2a)') 'refine ', p_name
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
subroutine write_diffev_single_part2(j, l, jj, iscat, lcontent, natoms, c_atom, uij_l, c_flag)
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
character(len=PREC_STRING) :: line 
integer :: m, k
!
if(c_flag=='fixed') then
if(l==1) then
   write(IDF_PO,'(a,i3,4a)') 'anis type:', jj,', values:[', 'U_',c_atom(iscat)(1:len_trim(c_atom(iscat))),'_1]'
   write(IDF_MA,'(3a,4(a2,g15.8e3))') 'newpara U_',c_atom(iscat)(1:len_trim(c_atom(iscat))),'_1' &
        , ', ', uij_l(1,iscat, j)*1.000 &
        , ', ', uij_l(1,iscat, j)*1.000 &
        , ', ', uij_l(1,iscat, j)*1.000 &
        , ', ', uij_l(1,iscat, j)*1.000  
   write(IDF_BS, '(4a)') 'variable real, ', 'U_',c_atom(iscat)(1:len_trim(c_atom(iscat))),'_1'
   line = '                 = '
   write(line(1:16), '(3a)') 'U_',c_atom(iscat)(1:len_trim(c_atom(iscat))),'_1'
   write(IDF_BS2,'(a19, G18.10E3)') line(1:19), uij_l(1,iscat, j)
   diffev_dimension = diffev_dimension + 1
elseif(l==6) then
   write(IDF_PO,'(a,i3,a,5(3a,i1.1,a2),3a)') 'anis type:', jj,', values:[', &
        ('U_',c_atom(iscat)(1:len_trim(c_atom(iscat))),'_',k,', ',k=1,5), &
         'U_',c_atom(iscat)(1:len_trim(c_atom(iscat))),'_6]'
   do m=1,6
      write(IDF_MA,'(3a,i1.1,4(a2,g15.8e3))') 'newpara U_',c_atom(iscat)(1:len_trim(c_atom(iscat))),'_',m &
        , ', ', uij_l(m,iscat, j)*1.000 &
        , ', ', uij_l(m,iscat, j)*1.000 &
        , ', ', uij_l(m,iscat, j)*1.000 &
        , ', ', uij_l(m,iscat, j)*1.000  
      write(IDF_BS, '(4a, i1.1)') 'variable real, ', 'U_',c_atom(iscat)(1:len_trim(c_atom(iscat))),'_',m
      line = '                 = '
      write(line(1:16), '(3a,i1.1)') 'U_',c_atom(iscat)(1:len_trim(c_atom(iscat))),'_',m
      write(IDF_BS2,'(a19, G18.10E3)') line(1:19), uij_l(m,iscat, j)
      diffev_dimension = diffev_dimension + 1
   enddo
endif
endif
!
end subroutine write_diffev_single_part2
!
!*****7**************************************************************** 
!
subroutine write_diffev_single_part3(&
   inline       , iatom, iscat, lcontent, natoms, c_atom, posit, c_flag)
!
!  Write main discus macro; single crystal PART 3
!  can be used for SINGLE and POWDER
!
use spcgr_apply, only: get_wyckoff
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
character(len=*)   , intent(in) :: c_flag
!
integer, parameter :: mode = 0
logical, parameter :: loutput = .true.
character(len=PREC_STRING) :: line 
integer :: i
real(kind=PREC_DP), dimension(3) :: coordinates
!
!     line = content(jc)(1:4)
coordinates = posit(:, iatom)
call get_wyckoff(coordinates, loutput, mode)
i = len_trim(inline)
write(IDF_PO,'(a,i3,3a)') 'x[',iatom, '] = P_',inline(1:i), '_x'
write(IDF_PO,'(a,i3,3a)') 'y[',iatom, '] = P_',inline(1:i), '_y'
write(IDF_PO,'(a,i3,3a)') 'z[',iatom, '] = P_',inline(1:i), '_z'
!
if(c_flag=='fixed') then
write(IDF_MA,'(3a,4(a2,f12.8)  )') 'newpara P_',c_atom(iscat)(1:len_trim(c_atom(iscat))),'_x', &
   ', ', posit(1, iatom)-0.02, ', ', posit(1, iatom)+0.02, &
   ', ', posit(1, iatom)-0.01, ', ', posit(1, iatom)+0.01
write(IDF_MA,'(3a,4(a2,f12.8)  )') 'newpara P_',c_atom(iscat)(1:len_trim(c_atom(iscat))),'_y', &
   ', ', posit(2, iatom)-0.02, ', ', posit(2, iatom)+0.02, &
   ', ', posit(2, iatom)-0.01, ', ', posit(2, iatom)+0.01
write(IDF_MA,'(3a,4(a2,f12.8)  )') 'newpara P_',c_atom(iscat)(1:len_trim(c_atom(iscat))),'_z', &
   ', ', posit(3, iatom)-0.02, ', ', posit(3, iatom)+0.02, &
   ', ', posit(3, iatom)-0.01, ', ', posit(3, iatom)+0.01
endif
!
write(IDF_FX, '(3a)') 'refine P_', c_atom(iscat)(1:len_trim(c_atom(iscat))),'_x'
write(IDF_FX, '(3a)') 'refine P_', c_atom(iscat)(1:len_trim(c_atom(iscat))),'_y'
write(IDF_FX, '(3a)') 'refine P_', c_atom(iscat)(1:len_trim(c_atom(iscat))),'_z'
!
write(IDF_BS,'(3a)') 'variable real, P_',c_atom(iscat)(1:len_trim(c_atom(iscat))),'_x'
write(IDF_BS,'(3a)') 'variable real, P_',c_atom(iscat)(1:len_trim(c_atom(iscat))),'_y'
write(IDF_BS,'(3a)') 'variable real, P_',c_atom(iscat)(1:len_trim(c_atom(iscat))),'_z'
line = '                 = '
write(line(1:16), '(3a)') 'P_',c_atom(iscat)(1:len_trim(c_atom(iscat))),'_x'
write(IDF_BS2,'(a19, G18.10E3)') line(1:19), posit(1, iatom)
line = '                 = '
write(line(1:16), '(3a)') 'P_',c_atom(iscat)(1:len_trim(c_atom(iscat))),'_y'
write(IDF_BS2,'(a19, G18.10E3)') line(1:19), posit(2, iatom)
line = '                 = '
write(line(1:16), '(3a)') 'P_',c_atom(iscat)(1:len_trim(c_atom(iscat))),'_z'
write(IDF_BS2,'(a19, G18.10E3)') line(1:19), posit(3, iatom)
diffev_dimension = diffev_dimension + 3
!
end subroutine write_diffev_single_part3
!
!*****7**************************************************************** 
!
subroutine write_diffev_powder_part4(&
           substance, c_form, lambda, slambda, hkl_file)
!-
! Write main discus macro DIFFEV_VERSION, last part
! Write "diffev_setup.mac"              , last part
!+
!
use precision_mod
use random_state_mod
!
implicit none
!
character(len=*), intent(in) :: substance    ! current material
character(len=*), intent(in) :: c_form       ! Atom form factor
character(len=*), intent(in) :: hkl_file     ! TSC file in case of DISCAMB
character(len=*), intent(in) ::  lambda      ! Wave length numerical value
character(len=*), intent(in) :: slambda      ! Wave length symbol
!integer, intent(in) :: spcgr_syst
!
character(len=PREC_STRING) :: line
integer :: i, i1
integer :: nseeds
integer, dimension(:), allocatable :: seed_val
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
write(*,*) ' RLAMBDA ', lambda(1:10), slambda
read(lambda, *) rlambda
call write_macro_powder(c_form, hkl_file, rlambda, slambda)
!
write(IDF_PO, '(a)') 'output'
write(IDF_PO, '(a)') '  outf "POWDER/%c.%4D.%4D", SUBSTANCE, REF_KID, REF_INDIV'
write(IDF_PO, '(a)') '  value inte'
write(IDF_PO, '(a)') '  form powder, tth, F_XMIN, F_XMAX, F_XSTP'
write(IDF_PO, '(a)') '  run'
write(IDF_PO, '(a)') 'exit'
write(IDF_PO, '(a)') 'branch kuplot'
write(IDF_PO, '(a)') 'reset'
write(IDF_PO, '(a)') 'load xy, "POWDER/%c.%4D.%4D", SUBSTANCE, REF_KID, REF_INDIV'
write(IDF_PO, '(a)') 'load xy, "%c/DATA/%c.tth", DATADIR, SUBSTANCE'
write(IDF_PO, '(a)') 'fit 2'
write(IDF_PO, '(a)') '  func back, 1, 5'
write(IDF_PO, '(a)') '  para 1, 0, 1.0'
write(IDF_PO, '(a)') '  para 2, 1, ymin[2]'
write(IDF_PO, '(a)') '  para 3, 1, 0.0'
write(IDF_PO, '(a)') '  para 4, 1, 0.0'
write(IDF_PO, '(a)') '  para 5, 1, 0.0'
write(IDF_PO, '(a)') '  cycle 10'
write(IDF_PO, '(a)') '  run'
write(IDF_PO, '(a)') 'exit'
write(IDF_PO, '(a)') 'ksav n[1]-1'
write(IDF_PO, '(a)') '  outf "INDI/indi.%4D.%4D.calc", REF_KID, REF_INDIV'
write(IDF_PO, '(a)') '  form xy'
write(IDF_PO, '(a)') '  run'
write(IDF_PO, '(a)') 'n[1] = n[1]-1'
write(IDF_PO, '(a)') 'exit  ! Back to DISCUS'
!
write(IDF_MA, '(a)') 'newpara  P_eta    ,  0.5000,   0.5000,   0.5000,    0.5000     ! Profile eta'
write(IDF_MA, '(a)') 'newpara  P_eta_l  ,  0.0000,   0.0000,   0.0000,    0.0000     ! Profile eta_linear'
write(IDF_MA, '(a)') 'newpara  P_eta_q  ,  0.0000,   0.0000,   0.0000,    0.0000     ! Profile eta_quad'
write(IDF_MA, '(a)') 'newpara  P_u      ,  0.0000,   0.0000,   0.0000,    0.0000     ! Profile u'
write(IDF_MA, '(a)') 'newpara  P_v      ,  0.0000,   0.0000,   0.0000,    0.0000     ! Profile v'
write(IDF_MA, '(a)') 'newpara  P_w      ,  0.0005,   0.0005,   0.0005,    0.0005     ! Profile w'
write(IDF_MA, '(a)') 'newpara  P_as1_c  ,  0.0000,   0.0000,   0.0000,    0.0000     ! Asymmetry const 1'
write(IDF_MA, '(a)') 'newpara  P_as2_c  ,  0.0000,   0.0000,   0.0000,    0.0000     ! Asymmetry const 2'
write(IDF_MA, '(a)') 'newpara  P_as1_i  ,  0.0000,   0.0000,   0.0000,    0.0000     ! Asymmetry inv   1'
write(IDF_MA, '(a)') 'newpara  P_as2_i  ,  0.0000,   0.0000,   0.0000,    0.0000     ! Asymmetry inv   2'
write(IDF_MA, '(a)') 'newpara  P_as1_l  ,  0.0000,   0.0000,   0.0000,    0.0000     ! Asymmetry lin   1'
write(IDF_MA, '(a)') 'newpara  P_as2_l  ,  0.0000,   0.0000,   0.0000,    0.0000     ! Asymmetry lin   2'
write(IDF_MA, '(a)') 'newpara  P_as1_q  ,  0.0000,   0.0000,   0.0000,    0.0000     ! Asymmetry quad  1'
write(IDF_MA, '(a)') 'newpara  P_as2_q  ,  0.0000,   0.0000,   0.0000,    0.0000     ! Asymmetry quad  2'
write(IDF_MA, '(a)') '#'
write(IDF_MA, '(a)') '#newpara  P_portion,  1.000,   1.0000,   1.000,    1.000     ! fixed '
write(IDF_MA, '(a)') '#newpara  P_damping,  0.000,   0.0000,   0.000,    0.000     ! fixed '
write(IDF_MA, '(a)') '#'
write(IDF_MA, '(a)') '@diffev_fix_free.mac         ! ''refine none'' etc '
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
write(IDF_BS, '(a)') '#variable real, P_portion'
write(IDF_BS, '(a)') '#variable real, P_damping'
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
write(IDF_BS, '(a)') 'discus'
write(IDF_BS, '(a)') '@dis.diffev.mac  ., REF_KID, REF_INDIV'
write(IDF_BS, '(a)') '#'
write(IDF_BS, '(a)') '#'
write(IDF_BS, '(a)') 'set error, continue'
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
close(IDF_PO)
close(IDF_MA)
close(IDF_FX)
close(IDF_BS)
close(IDF_BS2)
!
open(unit=IDF_PO, file='kup.diffev.mac', status='unknown')
write(IDF_PO, '(a)') 'reset'
write(IDF_PO, '(a)') 'load xy, "INDI/indi.%4D.*.calc", REF_KID'
write(IDF_PO, '(a)') 'merge all'
write(IDF_PO, '(a)') 'scale'
write(IDF_PO, '(a)') 'mark'
write(IDF_PO, '(a)') 'ksav n[1]'
write(IDF_PO, '(a)') '  outfile "TEMP/calc.%4D", REF_KID'
write(IDF_PO, '(a)') '  form xy'
write(IDF_PO, '(a)') '  run'
write(IDF_PO, '(a)') '#'
write(IDF_PO, '(3a)') 'load xy, DATA/',substance(1:len_trim(substance)),'.tth'
write(IDF_PO, '(a)') 'rval n[1], n[1]-1, dat'
write(IDF_PO, '(a)') 'exit'
write(IDF_PO, '(a)') '#'
write(IDF_PO, '(a)') 'exit   ! Back to DISCUS'
!
close(IDF_PO)

!
deallocate(seed_val)
!
end subroutine write_diffev_powder_part4
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
write(IDI,'(a )') 'variable real, F_XMIN'
write(IDI,'(a )') 'variable real, F_XMAX'
write(IDI,'(a )') 'variable real, F_XSTP'
if(rlambda/=0.0D0) then
   write(IDI,'(a, f8.6)') 'rlambda = ', rlambda
else
   write(IDI,'(a )') 'rlambda = 1.54061'
endif
write(IDI,'(a )') 'F_XMIN  =  10.0000'
write(IDI,'(a )') 'F_XMAX  = 100.0000'
write(IDI,'(a )') 'F_XSTP  =   0.0050'
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
end module prep_refine_mod
