module exp2pdf_write_mod
!-
!  Writing routines for exp2pdf
!+
!
contains
!
!*******************************************************************************
!
subroutine exp2pdf_write_unified(  datafile, ncols, cvalue, exp_data, exp_back, exp_sigma, &
                exp_bscale, &
                exp_radiation, exp_rlambda, exp_renergy, &
                exp_natom, exp_atname, exp_atocc, exp_comp_current, &
                exp_npoly, exp_qmin_i, exp_qmax_u, &
                exp_qmin_f, exp_qmax_f, qmin, qmax, &
                exp_npoints ,   exp_x, exp_values, &
                unit_cell_lengths, unit_cell_angles,               &
                symmetry_H_M, symmetry_origin, symmetry_abc, symmetry_n_mat, &
                symmetry_mat)
!-
!  Write the file in the Diffuse developers common file format
!+
!
use blanks_mod
use envir_mod
use errlist_mod
use unified_write_mod
use precision_mod
use prompt_mod
!
implicit none
!
character(len=*)                          , intent(in) ::   datafile
integer                                   , intent(in) :: ncols
character(len=*)  , dimension(2)          , intent(in) :: cvalue
character(len=*)                          , intent(in) :: exp_data
character(len=*)                          , intent(in) :: exp_back
character(len=*)                          , intent(in) :: exp_sigma
real(kind=PREC_DP)                        , intent(in) :: exp_bscale
character(len=*)                          , intent(in) :: exp_radiation
real(kind=PREC_DP)                        , intent(in) :: exp_rlambda
real(kind=PREC_DP)                        , intent(in) :: exp_renergy
integer                                   , intent(in) :: exp_natom
character(len=4)  , dimension(exp_natom)  , intent(in) :: exp_atname
real(kind=PREC_DP), dimension(exp_natom)  , intent(in) :: exp_atocc
logical                                   , intent(in) :: exp_comp_current
integer                                   , intent(in) :: exp_npoly
real(kind=PREC_DP)                        , intent(in) :: exp_qmin_i
real(kind=PREC_DP)                        , intent(in) :: exp_qmax_u
real(kind=PREC_DP)                        , intent(in) :: exp_qmin_f
real(kind=PREC_DP)                        , intent(in) :: exp_qmax_f
real(kind=PREC_DP)                        , intent(in) :: qmin
real(kind=PREC_DP)                        , intent(in) :: qmax
integer                                   , intent(in) :: exp_npoints
real(kind=PREC_DP), dimension(0:exp_npoints), intent(in) :: exp_x
real(kind=PREC_DP), dimension(exp_npoints), intent(in) :: exp_values
!
real(kind=PREC_DP), dimension(3)            , intent(in)  :: unit_cell_lengths
real(kind=PREC_DP), dimension(3)            , intent(in)  :: unit_cell_angles
character(len=32)                           , intent(in)  :: symmetry_H_M       ! Hermann-Mauguin Symbol
integer                                     , intent(in)  :: symmetry_origin    ! Space group origin 1 / 2
character(len=3)                            , intent(in)  :: symmetry_abc       ! Permutation for orthorhombic space groups
integer                                     , intent(in)  :: symmetry_n_mat     ! Number of symmetry elements
real(kind=PREC_DP)        , dimension(:,:,:), allocatable, intent(in) :: symmetry_mat       ! Actual Symmetry matrices
!
character(len=32)                                         :: data_type_experiment ! 'experimental', 'calculated'
character(len=32)                                         :: data_type_style      ! 'powder', 'powder_pdf', 'single_diffraction', 'single_pdf' ...
character(len=32)                                         :: data_type_axes       ! 'hkl', 'Q', 'theta', 'theta', 'uvw', 'dstar', 'sinR(theta)/lambda'
character(len=32)                                         :: data_type_content    ! 'intensity', '3d-delta-pdf', 'amplitide', 'density' ...
character(len=32)                                         :: data_type_reciprocal ! 'reciprocal', 'patterson', 'direct'
character(len=32)                                         :: data_type_with_bragg ! With Bragg   = 0, No Bragg   = 1
character(len=32)                                         :: data_type_symmetrized! No = 0; Point=1; Laue=2; Space=3;SymElem=4
character(len=32)                                         :: data_type_number     ! 'real; 'complex'
character(len=32)                                         :: data_rad_radiation   ! Xray=0; Neutron=1; Electron=2
character(len=16)                                         :: data_rad_symbol      ! CU, CUA1, CU12 ...
real(kind=PREC_DP)        , dimension(3)                  :: data_rad_length      ! CU, CUA1, CU12 ...
integer                                                   :: data_abs_is_hkl      ! Abscissa is 1=h 2=k 3=l
integer                                                   :: data_ord_is_hkl      ! Ordinate is 1=h 2=k 3=l
integer                                                   :: data_top_is_hkl      ! top-axis is 1=h 2=k 3=l
integer                   , dimension(3   )               :: data_dimension       ! Data have dimensions along (HKL) / (UVW)
real(kind=PREC_DP)        , dimension(3   )               :: data_corner          ! Lower left bottom corner in fractional coordinates
character(len=32)                                         :: coordinate_unit      ! Axis units 'basecell_fractional' ...
real(kind=PREC_DP)        , dimension(3, 3)               :: data_vector          ! Increment vectors abs: (:,1); ord: (:,2); top: (:,3)
! PDF info
character(len=PREC_STRING)           :: sample_composition     ! string like "Zr1.0, O2.0"
character(len=PREC_STRING)           :: source_conversion      ! "DISCUS_SINGLE 6.22.3"
character(len=PREC_STRING), dimension(3) :: source_file        ! "DISCUS_CALCULATION"
!
real(kind=PREC_DP), dimension(1)     :: source_back_scale      ! 0.0
integer, dimension(3)                :: source_dimension       ! Dimensions in reciprocal space
!integer                              :: source_abs_is_hkl      ! Source Abscissa is 1=h 2=k 3 = k
!integer                              :: source_ord_is_hkl      ! Source Ordinate is 1=h 2=k 3=l
!integer                              :: source_top_is_hkl      ! Source top-axis is 1=h 2=k 3=l
!character(len=32)                    :: source_coordinate_unit ! 'basecell_fractional'
!real(kind=PREC_DP), dimension(3)     :: source_corner          ! Reciprocal lower left corner
!real(kind=PREC_DP), dimension(3, 3)  :: source_vector          ! Reciprocal increment vecors abs, ord, top
real(kind=PREC_DP), dimension(2)     :: source_inst_qmin_qmax  ! Powder qmin, qmax Fit range
real(kind=PREC_DP), dimension(2)     :: source_fourier_qmin_qmax  ! Powder qmin, qmax Fourier range
real(kind=PREC_DP), dimension(1)     :: source_correction         ! Powder poly correction parameters
!
character(len=PREC_STRING), dimension(  5) :: crystal_meta       ! Metadata, see "c_meta"
!
real(kind=PREC_DP), dimension(:,:,:), allocatable :: data_values
!
character(len=8)   :: date    ! Date string
integer :: NMSG, i, i1, i2
!
NMSG = ubound(ier_msg,1)
do i=1, exp_natom
   i1 = (i-1)*13 +  1
   i2 = (i-1)*13 + 13
   write(sample_composition(i1:i2),'(a4,a1,f7.4,a1)') exp_atname(i), ':',exp_atocc(i), ','
enddo
sample_composition(i2:i2) = ' '
call rem_bl(sample_composition, i2)
!
if(cvalue(2)(1:4)=='G(r)') then    ! Powder PDF  
   data_type_experiment  = 'experimental'
   data_type_style       = 'powder_pdf'
   data_type_axes        = 'r'
   data_type_content     = 'powder_pdf'
   data_type_reciprocal  = 'direct'
   data_type_with_bragg  = 'bragg_present'
   data_type_symmetrized = 'none'
   data_type_number      = 'real'
   data_rad_radiation    = exp_radiation
   data_rad_symbol       = ' '
   data_rad_length       = exp_rlambda
   data_dimension        = (/exp_npoints, 1, 1/)
   data_abs_is_hkl       = 1
   data_ord_is_hkl       = 2
   data_top_is_hkl       = 3
   coordinate_unit       = 'angstroem'
   data_corner           = 0.0_PREC_DP
   data_corner(1)        = exp_x(0)
   data_vector           = 0.0_PREC_DP
   data_vector(1, 1)     = (exp_x(exp_npoints-1)-exp_x(0))/(real(exp_npoints-1))
   allocate(data_values(exp_npoints, 1, 1))
   data_values(1:exp_npoints, 1, 1) = exp_values
!
!  crystal_meta
   call date_and_time(date=date)
   crystal_meta(1) = 'Disorder unified data'
   crystal_meta(2) = '0.0.0'
   crystal_meta(3) = date(1:4) // '-' // date(5:6) // '-' // date(7:8)
   crystal_meta(4) = version_discus
   crystal_meta(5) = author
   source_conversion     = 'discus_exp2pdf ' // version
!
   source_file(1)        = exp_data
   source_file(2)        = exp_back
   source_file(3)        = exp_sigma
!
   source_back_scale(1)     = exp_bscale
   source_inst_qmin_qmax(1) = exp_qmin_i
   source_inst_qmin_qmax(2) = exp_qmax_u
   source_fourier_qmin_qmax(1) = exp_qmin_f
   source_fourier_qmin_qmax(2) = exp_qmax_f
   source_correction(1)  = exp_npoly
!
   call unified_write_data(datafile, unit_cell_lengths, unit_cell_angles, &
                           symmetry_H_M, symmetry_origin, symmetry_abc, symmetry_n_mat, &
                           symmetry_mat,                                                &
                           data_type_experiment, &
                           data_type_style     , &
                           data_type_axes      , &
                           data_type_content   , &
                           data_type_reciprocal, &
                           data_type_with_bragg, &
                           data_type_symmetrized, &
                           data_type_number     , &
                           data_rad_radiation , data_rad_symbol, data_rad_length, &
                           data_dimension  , &
                           data_abs_is_hkl     , &
                           data_ord_is_hkl     , &
                           data_top_is_hkl     , &
                           coordinate_unit , &
                           data_corner     , &
                           data_vector     , &
                           data_values     , &
                           NMSG, ier_num, ier_msg,                                      &
                           crystal_meta,                  &
                           sample_composition = sample_composition, &
                           source_conversion  = source_conversion , &
                           source_file        = source_file       , &
                           source_back_scale  = source_back_scale , &
                           source_inst_qmin_qmax   = source_inst_qmin_qmax, &
                           source_fourier_qmin_qmax   = source_fourier_qmin_qmax, &
                           source_correction   = source_correction &
                           )
endif
!
deallocate(data_values)
!
end subroutine exp2pdf_write_unified
!
!*******************************************************************************
!
end module exp2pdf_write_mod
