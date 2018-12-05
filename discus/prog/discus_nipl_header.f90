MODULE discus_nipl_header
!
CONTAINS
!
SUBROUTINE write_discus_nipl_header(header_lines, nheader, layer)
!
!USE crystal_mod
USE diffuse_mod
IMPLICIT NONE
!
CHARACTER (LEN=160), DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: header_lines
INTEGER,                                        INTENT(OUT) :: nheader! number of lines in header
INTEGER,                                        INTENT(IN ) :: layer  ! Current layer along the top axis
!
CHARACTER (LEN=160), DIMENSION(:), ALLOCATABLE :: lines_fourier
CHARACTER (LEN=160), DIMENSION(:), ALLOCATABLE :: lines_crystal
!NTEGER                         :: nheader
INTEGER                         :: nsection
INTEGER                         :: ncrystal
INTEGER                         :: nfourier
INTEGER                         :: i, j
!
nheader  = 1
nsection = 0
ncrystal = 0
nfourier = 0
!
!  Always write a basic crystal header, this is the first section
!
CALL write_discus_nipl_crystal(ncrystal, lines_crystal)
nheader  = nheader + 1
nsection = nsection + 1
!
!  Test the various Fourier types (FOURIER, INVERSE, PATTERSON)
!
IF(IABS(four_last)< 3) THEN               ! A Fourier was calculated
   CALL write_discus_nipl_fourier(nfourier, lines_fourier, layer)
   nheader  = nheader + 1
   nsection = nsection + 1
ENDIF
!
nheader = nheader + ncrystal + nfourier
ALLOCATE(header_lines(nheader))
header_lines(:) = ' '
!
j = 1
WRITE(header_lines(j), '(''# HEADER  : lines, sec  ;2I::'' , i8,'', '',i8)') nheader, nsection
j = 2
WRITE(header_lines(j), '(''# HEADER  : section     ;1C::'' , a )') 'CRYSTAL'
!
IF(IABS(four_last)< 3) THEN               ! A Fourier was calculated
   j = j + 1
   WRITE(header_lines(j), '(''# HEADER  : section     ;1C::'' , a )') 'FOURIER'
ENDIF
!
DO i=1, ncrystal
!  WRITE(iff,'(a)') lines_crystal(i)(1:LEN_TRIM(lines_crystal(i)))
   j = j + 1
   header_lines(j) = lines_crystal(i)(1:LEN_TRIM(lines_crystal(i)))
ENDDO
!
!  There is always a "FOURIER" / "INVERSE" / "PATTERSON" section 
!
DO i=1, nfourier
!  WRITE(iff,'(a)') lines_fourier(i)(1:LEN_TRIM(lines_fourier(i)))
   j = j + 1
   header_lines(j) = lines_fourier(i)(1:LEN_TRIM(lines_fourier(i)))
ENDDO
!
IF(ALLOCATED(lines_crystal)) DEALLOCATE(lines_crystal)
IF(ALLOCATED(lines_fourier)) DEALLOCATE(lines_fourier)
!
END SUBROUTINE write_discus_nipl_header
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
SUBROUTINE write_discus_nipl_crystal(ncrystal, lines_crystal)
!
USE crystal_mod
USE chem_aver_mod
!
IMPLICIT NONE
!
INTEGER,                                        INTENT(OUT) :: ncrystal
CHARACTER (LEN=160), DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: lines_crystal
!
INTEGER  :: i
!
CALL chem_elem(.FALSE.)
ncrystal = 2 + cr_nscat + 1
ALLOCATE(lines_crystal(1:ncrystal))
WRITE(lines_crystal(1), '(''# CRYSTAL : lines       ;1I:: '', i8)') ncrystal
WRITE(lines_crystal(2), '(''# CRYSTAL : no. atoms   ;3I:: '', 2(i12,'', ''),I12)') cr_natoms, cr_nscat, cr_ncatoms
DO i=0, cr_nscat
   WRITE(lines_crystal(3+i), '(''# CRYSTAL : type atom   ;1C 2I 1R:: '',a5,'', '',2(i12,'', ''),f8.5)') &
         cr_at_lis(i), i, cr_amount(i), cr_dw(i)
ENDDO
!
END SUBROUTINE write_discus_nipl_crystal
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
SUBROUTINE write_discus_nipl_fourier(nfourier, lines_fourier, layer)
!
USE crystal_mod
USE diffuse_mod
USE four_angles_mod
USE output_mod
!
IMPLICIT NONE
!
INTEGER,                                        INTENT(OUT) :: nfourier
CHARACTER (LEN=160), DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: lines_fourier
INTEGER,                                        INTENT(IN ) :: layer  ! Current layer along the top axis
!
CHARACTER (LEN=8), DIMENSION(3) :: cradiation
CHARACTER (LEN=1), DIMENSION(3) :: caxis
INTEGER                         :: i
LOGICAL                         :: ltop
REAL                            ::  angle_vh
REAL                            ::  ratio_vh
REAL                            ::   aver_vh
REAL                            ::  angle_ht
REAL                            ::  ratio_ht
REAL                            ::   aver_ht
REAL                            ::  angle_tv
REAL                            ::  ratio_tv
REAL                            ::   aver_tv
REAL            , DIMENSION(3)  ::  length = (/0.0, 0.0, 0.0/)
!
DATA cradiation /'xray    ','neutron ', 'electron'/
DATA caxis      /'h','k','l'/
!
nfourier = 18
ALLOCATE(lines_fourier(1:nfourier))
!
ltop = .NOT. (eck(1,1) == eck(1,4) .AND. eck(2,1) == eck(2,4) .AND. &
              eck(3,1) == eck(3,4) )                                &
       .OR. inc(3) > 1
CALL four_angles(ltop, length, angle_vh, ratio_vh, aver_vh, &
                               angle_ht, ratio_ht, aver_ht, &
                               angle_tv, ratio_tv, aver_tv)
!
!     Calculate lengths in Ang-1
!
WRITE(lines_fourier( 1), '(''# FOURIER : lines       ;1I:: '', i8)') nfourier
if(four_last<0) THEN
   WRITE(lines_fourier( 2), '(''# FOURIER : type        ;1C:: '', i1,a)') IABS(four_last),'D Fourier with lots'
ELSE
   WRITE(lines_fourier( 2), '(''# FOURIER : type        ;1C:: '', i1,a)') IABS(four_last),'D Fourier complete crystal'
ENDIF
WRITE(lines_fourier( 3), '(''# FOURIER : radiation   ;1C:: '', a )') cradiation(diff_radiation)
WRITE(lines_fourier( 4), '(''# FOURIER : wave length ;1R::'' ,   f12.5               )') rlambda
i=1
IF(lambda==' ') i=0
WRITE(lines_fourier( 5), '(''# FOURIER : wave name   ;'',i1,''C::'' ,   a           )') i,lambda
IF(ldbw) THEN
   WRITE(lines_fourier( 6), '(''# FOURIER : adp         ;1C:: '', a )') 'used'
ELSE
   WRITE(lines_fourier( 6), '(''# FOURIER : adp         ;1C:: '', a )') 'ignored'
ENDIF
IF(ano) THEN
   WRITE(lines_fourier( 7), '(''# FOURIER : dispersion  ;1C:: '', a )') 'anomalous'
ELSE
   WRITE(lines_fourier( 7), '(''# FOURIER : dispersion  ;1C:: '', a )') 'ignored'
ENDIF
WRITE(lines_fourier( 8), '(''# FOURIER : lower left  ;3R::'' , 2(f12.5,'', ''),f12.5 )') eck(:,1)
WRITE(lines_fourier( 9), '(''# FOURIER : lower right ;3R::'' , 2(f12.5,'', ''),f12.5 )') eck(:,2)
WRITE(lines_fourier(10), '(''# FOURIER : upper left  ;3R::'' , 2(f12.5,'', ''),f12.5 )') eck(:,3)
WRITE(lines_fourier(11), '(''# FOURIER : top   left  ;3R::'' , 2(f12.5,'', ''),f12.5 )') eck(:,4)
WRITE(lines_fourier(12), '(''# FOURIER : no. points  ;3I::'' , 2(i12  ,'', ''),i12   )') inc(:)
WRITE(lines_fourier(13), '(''# FOURIER : abscissa    ;1C:: '', a )') caxis(extr_abs)
WRITE(lines_fourier(14), '(''# FOURIER : ordinate    ;1C:: '', a )') caxis(extr_ord)
WRITE(lines_fourier(15), '(''# FOURIER : top         ;1C:: '', a )') caxis(extr_top)
WRITE(lines_fourier(16), '(''# FOURIER : current top ;1I:: '', i12 )') layer
WRITE(lines_fourier(17), '(''# FOURIER : angles      ;3R::'' , 2(f12.5,'', ''),f12.5 )') angle_vh, angle_ht, angle_tv
WRITE(lines_fourier(18), '(''# FOURIER : aver        ;3R::'' , 2(f12.5,'', ''),f12.5 )')  aver_vh,  aver_ht,  aver_tv
!
END SUBROUTINE write_discus_nipl_fourier
!
END MODULE discus_nipl_header
