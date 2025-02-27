module prep_refine_mod
!
private
public write_refine_single_part1
public write_refine_single_part2
public write_refine_single_part3
public write_refine_single_part4
public write_refine_powder_part1
public write_refine_powder_part4
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
elseif(spcgr_syst==8) then         ! cubic      
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
elseif(spcgr_syst==7) then         ! Triclinic
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
write(IDI,'(a )') 'powder'
if(c_form == 'discamb') then
   write(IDI,'(3a)') '  xray table:discamb, file:', hkl_file(1:j), '.tsc'
else
   write(IDI,'(2a)') '  xray table:', c_form(1:len_trim(c_form))
endif
write(IDI,'(a )') '  set axis, q'
write(IDI,'(a )') '  set calc, complete'
write(IDI,'(a,f7.5,a )') '  set qmin, 4.*PI*sind(F_XMIN*0.5)/ rlambda -0.2'
write(IDI,'(a,f7.5,a )') '  set qmax, 4.*PI*sind(F_XMAX*0.5)/ rlambda +0.2'
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
end module prep_refine_mod
