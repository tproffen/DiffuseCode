!*****7****************************************************************
!
      subroutine errlist_appl
!-
!     Displays error Messages for the error type APPLication
!+
      implicit      none
!
      include      'errlist.inc'
!
      integer       iu,io
      PARAMETER    (IU=-112,IO=4)
!
      CHARACTER*41  ERROR(IU:IO)
!
      DATA ERROR (-112: -101) /                     &
     &  'Error reading ADP  instruction',           & !-112 ! discus
     &  'Error reading SCAT instruction',           & !-111 ! discus
     &  'Connectivity definition does not exist',   & !-110 ! discus
     &  'Connectivity definition does not exist',   & !-109 ! discus
     &  'Q limits or step width are illegal',       & !-108 ! discus
     &  '2Theta limits or step width are illegal',  & !-107 ! discus
     &  'HKL steps for complete powder must be >0', & !-106 ! discus
     &  'Atom number outside limits',               & !-105 ! discus
     &  'Powder output not defined as TTH or Q',    & !-104 ! discus
     &  'Too many atoms per unit cell',             & !-103 ! discus
     &  'Property value outside defined range',     & !-102 ! discus
     &  'Dimension of lots < 0 or > than crystal'   & !-101 ! discus
     &  /
      DATA ERROR (-100:  -81) /                     &
     &  'Space group symbol missing in cell file',  & !-100 ! discus
     &  'No wavelength has been set ',              & !-99  ! discus
     &  'No atom types exist at present ',          & !-98  ! discus
     &  'Atom type outside proper limits',          & !-97  ! discus
     &  'Illegal keyword in domain input file',     & !-96  ! discus
     &  'First domain keyword has parameters',      & !-95  ! discus
     &  'Invalid domain descriptor in input file',  & !-94  ! discus
     &  'Unit cell constants <= zero',              & !-93  ! discus
     &  'Error reading generators from structure',  & !-92  ! discus
     &  'Unexpected pseudoatom name read',          & !-91  ! discus
     &  'Unknown diffractometer geometry',          & !-90  ! discus
     &  'Unknown keyword in unit cell file',        & !-89  ! discus
     &  'Bravais types differ',                     & !-88  ! discus
     &  'Different lattice constants',              & !-87  ! discus
     &  'Unknown import format',                    & !-86  ! discus
     &  'Sharpened patternson requires HKLF4 file', & !-85  ! discus
     &  'Error reading molecule parameters',        & !-84  ! discus
     &  'Molecule buildup failed',                  & !-83  ! discus
     &  'Invalid molecule character',               & !-82  ! discus
     &  'Invalid flag for space type'               & !-81  ! discus
     &  /
      DATA ERROR ( -80:  -61) /                     &
     &  'Angle between Normal and Abszissa is 0',   & !-80  ! discus
     &  'Too many atoms in result array',           & !-79  ! discus
     &  'Wrong optional parameter for HKLF4 format',& !-78  ! discus
     &  'Unknown wave length symbol used          ',& !-77  ! discus
     &  'Function not available for molecules     ',& !-76  ! discus
     &  'No bond valence parameters for atom pair ',& !-75  ! discus
     &  'Molecule atom-number outside limits      ',& !-74  ! discus
     &  'Transform. requires too many generators  ',& !-73  ! discus
     &  'Too many different atom types in file    ',& !-72  ! discus
     &  'Too many points in direct space layer    ',& !-71  ! discus
     &  'Delta value must be in interval 0 -> 1   ',& !-70  ! discus
     &  'Microdomains overlap                     ',& !-69  ! discus
     &  'Mode only available for molecules        ',& !-68  ! discus
     &  'Molecules have different number of atoms ',& !-67  ! discus
     &  'Too many molecule types created          ',& !-66  ! discus
     &  'Too many molecules created               ',& !-65  ! discus
     &  'Molecule type outside limits             ',& !-64  ! discus
     &  'Molecule number outside limits           ',& !-63  ! discus
     &  'Too many additional symmetry operators   ',& !-62  ! discus
     &  'Too many additional generators           ' & !-61  ! discus
     &  /
      DATA ERROR ( -60:  -41) /                     &
     &  'Output value NOT allowed using lots      ',& !-60  ! discus
     &  'Invalid color or typ selected for atom   ',& !-59  ! discus
     &  'No atoms written to file                 ',& !-58  ! discus
     &  'Av. Translation is zero                  ',& !-57  ! discus
     &  'Av. Transl. in plane of moduli vectors   ',& !-56  ! discus
     &  'No layers created at all                 ',& !-55  ! discus
     &  'Index outside limits                     ',& !-54  ! discus
     &  'Too many different layer types           ',& !-53  ! discus
     &  'Unsuitable input value for SHELXL format ',& !-52  ! discus
     &  'Unsuitable file type for SHELXL format   ',& !-51  ! discus
     &  'Wrong format for 1-dimensional data      ',& !-50  ! discus
     &  'Error reading atom coordinates           ',& !-49  ! discus
     &  'Error reading lattice constants          ',& !-48  ! discus
     &  'Error reading space group symbol         ',& !-47  ! discus
     &  'Error reading title of structure         ',& !-46  ! discus
     &  'Too many atoms in environment',            & !-45  ! discus
     &  'Right quotation mark missing in format',   & !-44  ! discus
     &  'Not enough parameter for filename format ',& !-43  ! discus
     &  'Type must be: inten,ampl,phase,real,imag ',& !-42  ! discus
     &  'File specifier must be "a" or "b"'         & !-41  ! discus
     &  /
      DATA ERROR ( -40:  -21) /                     &
     &  'Microdomain type cannot be removed',       & !-40  ! discus
     &  'All elements of correlation matrix zero',  & !-39  ! discus
     &  'Unsuitable file types for Patterson',      & !-38  ! discus
     &  'No filename defined yet            ',      & !-37  ! discus
     &  'Unsuitable file types for inverse Fourier',& !-36  ! discus
     &  'Volume of unit cell <= zero',              & !-35  ! discus
     &  'Form does not appear to be closed',        & !-34  ! discus
     &  'No microdomain types defined yet',         & !-33  ! discus
     &  'Length of vector is zero',                 & !-32  ! discus
     &  'Unknown distribution mode',                & !-31  ! discus
     &  'Unknown boundary type',                    & !-30  ! discus
     &  'Too many different microdomain types',     & !-29  ! discus
     &  'Input parameters must be > zero',          & !-28  ! discus
     &  'No atom of this type present in crystal',  & !-27  ! discus
     &  'Too many different atoms in crystal',      & !-26  ! discus
     &  'Different coordinates in 1-dim files',     & !-25  ! discus
     &  'Different coordinates in GNUPLOT files',   & !-24  ! discus
     &  'Incompatible GNUPLOT file sizes',          & !-23  ! discus
     &  'Incompatible standard file sizes',         & !-22  ! discus
     &  'No element present, no Fourier calculated' & !-21  ! discus
     &  /
      DATA ERROR ( -20:   -1) /                     &
     &  'Unknown element, no Fourier calculated',   & !-20  ! discus
     &  'Atom number outside limits',               & !-19  ! discus
     &  'No orientation with this number exists',   & !-18  ! discus
     &  'Index of matrix outside limits',           & !-17  ! discus
     &  'Status for log must be : "on" or "off"',   & !-16  ! discus
     &  'No microdomain input file name defined',   & !-15  ! discus
     &  'Invalid space group & lattice constants',  & !-14  ! discus
     &  'Correlation matrix index outside limits',  & !-13  ! discus
     &  'Number of points must be > zero ',         & !-12  ! discus
     &  'Unknown threshold type',                   & !-11  ! discus
     &  'Too many Atoms in crystal',                & !-10  ! discus
     &  'Unknown Output Format   ',                 & ! -9  ! discus
     &  'Too many points in reciprocal layer',      & ! -8  ! discus
     &  'Unknown space group symbol         ',      & ! -7  ! discus
     &  'Unknown microdomain type',                 & ! -6  ! discus
     &  'Too many stacking layers within crystal',  & ! -5  ! discus
     &  'Extend of plotspace is zero',              & ! -4  ! discus
     &  'No atoms selected yet',                    & ! -3  ! discus
     &  'Improper limits for atom number',          & ! -2  ! discus
     &  'Maximum number of Orient. matrizes read'   & ! -1  ! discus
     &  /
      DATA ERROR (   0:   io) /                     &
     &  ' ',                                        & !  0  ! discus
     &  'Large distance between micro and host',    & ! +1  ! discus
     &  'Molecule generator is obsolete >help data',& ! +2  ! discus
     &  'Molecule symmetry  is obsolete >help data',& ! +3  ! discus
     &  'Monte Carlo level is obsolete >help mmc'   & ! +4  ! discus
     &           /
!
      if (ier_typ.eq.ER_RMC) then
        call errlist_rmc
      ELSEIF (ier_typ.eq.ER_MMC) then
        call errlist_mmc
      ELSEIF (ier_typ.eq.ER_CHEM) then
        call errlist_chem
      ELSEIF (ier_typ.eq.ER_FOUR) then
        call errlist_four
      ELSEIF (ier_typ.eq.ER_PDF) then
        call errlist_pdf
      ELSE
        call disp_error ('APPL',error,iu,io)
      endif
      end
!*****7****************************************************************
