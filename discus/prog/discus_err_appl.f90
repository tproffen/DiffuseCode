!*****7****************************************************************
!
      SUBROUTINE discus_errlist_appl
!-
!     Displays error Messages for the error type APPLication
!+
      USE errlist_mod 
      implicit      none
!
!
      integer       iu,io
      PARAMETER    (IU=-154,IO=6)
!
      CHARACTER(LEN=45) ::  ERROR(IU:IO)
!
      DATA ERROR ( IU : -141) /                          &
     &  'Angle Ordinate to (Normal/Abscissa) is 0',      & !-154 ! discus
     &  'Angle between Normal and Ordinate is 0',        & !-153 ! discus
     &  'Interactive plot for JMOL only',                & !-152 ! discus
     &  'AND and DEFAULT simultaneously in isprop',      & !-151 ! discus
     &  'Occupancy outside [0:1]',                       & !-150 ! discus
     &  'Error reading SCAT instruction',                & !-149 ! discus
     &  '1bar not at origin',                            & !-148 ! discus
     &  'Unknown export format',                         & !-147 ! discus
     &  'No. of atoms not an integer multiple of sites', & !-146 ! discus
     &  'Decoration name not recognized      ',          & !-145 ! discus
     &  'Too few atoms in ligand molecule    ',          & !-144 ! discus
     &  'Wrong parameters for this bond type ',          & !-143 ! discus
     &  'Cubeoctahedron only allowed in cubic systems.', & !-142 ! discus
     &  'Atoms are too close to each other  '            & !-141 ! discus
     &  /
      DATA ERROR (-140: -121) /                          &
     &  'CSSR file not allowed for read cell',           & !-140 ! discus
     &  'Layer type outside limits',                     & !-139 ! discus
     &  'Atom is already inside a molecule',             & !-138 ! discus
     &  'Atom number is outside crystal',                & !-137 ! discus
     &  'Connectivity name does not match',              & !-136 ! discus
     &  'Did not find a connectivity for this atom',     & !-135 ! discus
     &  'Connectivity list has not been created ',       & !-134 ! discus
     &  'Refinement param index outside limits  ',       & !-133 ! discus
     &  'Mismatch between corners and increments'  ,     & !-132 ! discus
     &  'No surface sites for ligands found'       ,     & !-131 ! discus
     &  'No surface atoms found                   ',     & !-130 ! discus
     &  'No decoration definition exists yet      ',     & !-129 ! discus
     &  'Could not add the new decoration         ',     & !-128 ! discus
     &  'Empty content file                       ',     & !-127 ! discus
     &  'H-M symbol in CIF file is a question mark',     & !-126 ! discus
     &  'S(Q), F(Q) require Q-axis                ',     & !-125 ! discus
     &  'Powder output type wrong /= I, S(Q), F(Q)',     & !-124 ! discus
     &  'Atoms are at identical positions',              & !-123 ! discus
     &  'Atom type number outside limits',               & !-122 ! discus
     &  'Error calculating x-position for powder'        & !-121 ! discus
     &  /
      DATA ERROR (-120: -101) /                          &
     &  'Conn. Name is equal to variable name',          & !-120 ! discus
     &  'Error reading atom number from RMCPROFILE',     & !-119 ! discus
     &  'No Fourier calculated yet, no output',          & !-118 ! discus
     &  'This DISCUS has NeXus suppport disabled  ',     & !-117 ! discus
     &  'Could not find definition',                     & !-116 ! discus
     &  'Different atom no on SCAT and ADP',             & !-115 ! discus
     &  'Error allocating              ',                & !-114 ! discus
     &  'Could not find internal storage',               & !-113 ! discus
     &  'Error reading ADP  instruction',                & !-112 ! discus
     &  'Error reading SCAT instruction',                & !-111 ! discus
     &  'No connectivity definitions exist at all',      & !-110 ! discus
     &  'Connectivity definition does not exist',        & !-109 ! discus
     &  'Q limits or step width are illegal',            & !-108 ! discus
     &  '2Theta limits or step width are illegal',       & !-107 ! discus
     &  'HKL steps for complete powder must be >0',      & !-106 ! discus
     &  'Atom number outside limits',                    & !-105 ! discus
     &  'Powder output not defined as TTH or Q',         & !-104 ! discus
     &  'Too many atoms per unit cell',                  & !-103 ! discus
     &  'Property value outside defined range',          & !-102 ! discus
     &  'Dimension of lots < 0 or > than crystal'        & !-101 ! discus
     &  /
      DATA ERROR (-100:  -81) /                          &
     &  'Space group symbol missing in cell file',       & !-100 ! discus
     &  'No wavelength has been set ',                   & !-99  ! discus
     &  'No atom types exist at present ',               & !-98  ! discus
     &  'Atom type outside proper limits',               & !-97  ! discus
     &  'Illegal keyword in domain input file',          & !-96  ! discus
     &  'First domain keyword has parameters',           & !-95  ! discus
     &  'Invalid domain descriptor in input file',       & !-94  ! discus
     &  'Unit cell constants <= zero',                   & !-93  ! discus
     &  'Error reading generators from structure',       & !-92  ! discus
     &  'Unexpected pseudoatom name read',               & !-91  ! discus
     &  'Unknown diffractometer geometry',               & !-90  ! discus
     &  'Unknown keyword in unit cell file',             & !-89  ! discus
     &  'Bravais types differ',                          & !-88  ! discus
     &  'Different lattice constants',                   & !-87  ! discus
     &  'Unknown import format',                         & !-86  ! discus
     &  'Sharpened patternson requires HKLF4 file',      & !-85  ! discus
     &  'Error reading molecule parameters',             & !-84  ! discus
     &  'Molecule buildup failed',                       & !-83  ! discus
     &  'Invalid molecule character',                    & !-82  ! discus
     &  'Invalid flag for space type'                    & !-81  ! discus
     &  /
      DATA ERROR ( -80:  -61) /                          &
     &  'Angle between Normal and Abszissa is 0',        & !-80  ! discus
     &  'Too many atoms in result array',                & !-79  ! discus
     &  'Wrong optional parameter for HKLF4 format',     & !-78  ! discus
     &  'Unknown wave length symbol used          ',     & !-77  ! discus
     &  'Function not available for molecules     ',     & ! !-76  ! discus
     &  'No bond valence parameters for atom pair ',     & ! !-75  ! discus
     &  'Molecule atom-number outside limits      ',     & ! !-74  ! discus
     &  'Transform. requires too many generators  ',     & ! !-73  ! discus
     &  'Too many different atom types in file    ',     & ! !-72  ! discus
     &  'Too many points in direct space layer    ',     & ! !-71  ! discus
     &  'Delta value must be in interval 0 -> 1   ',     & ! !-70  ! discus
     &  'Microdomains overlap                     ',     & ! !-69  ! discus
     &  'Mode only available for molecules        ',     & ! !-68  ! discus
     &  'Molecules have different number of atoms ',     & ! !-67  ! discus
     &  'Too many molecule types created          ',     & ! !-66  ! discus
     &  'Too many molecules created               ',     & ! !-65  ! discus
     &  'Molecule type outside limits             ',     & ! !-64  ! discus
     &  'Molecule number outside limits           ',     & ! !-63  ! discus
     &  'Too many additional symmetry operators   ',     & ! !-62  ! discus
     &  'Too many additional generators           '      & ! !-61  ! discus
     &  /
      DATA ERROR ( -60:  -41) /                     &
     &  'Output value NOT allowed using lots      ',     & ! !-60  ! discus
     &  'Invalid color or typ selected for atom   ',     & ! !-59  ! discus
     &  'No atoms written to file                 ',     & ! !-58  ! discus
     &  'Av. Translation is zero                  ',     & ! !-57  ! discus
     &  'Av. Transl. in plane of moduli vectors   ',     & ! !-56  ! discus
     &  'No layers created at all                 ',     & ! !-55  ! discus
     &  'Index outside limits                     ',     & ! !-54  ! discus
     &  'Too many different layer types           ',     & ! !-53  ! discus
     &  'Unsuitable input value for SHELXL format ',     & ! !-52  ! discus
     &  'Unsuitable file type for SHELXL format   ',     & ! !-51  ! discus
     &  'Wrong format for 1-dimensional data      ',     & ! !-50  ! discus
     &  'Error reading atom coordinates           ',     & ! !-49  ! discus
     &  'Error reading lattice constants          ',     & ! !-48  ! discus
     &  'Error reading space group symbol         ',     & ! !-47  ! discus
     &  'Error reading title of structure         ',     & ! !-46  ! discus
     &  'Too many atoms in environment',                 & ! !-45  ! discus
     &  'Right quotation mark missing in format',        & ! !-44  ! discus
     &  'Not enough parameter for filename format ',     & ! !-43  ! discus
     &  'Type must be: inten,ampl,phase,real,imag ',     & ! !-42  ! discus
     &  'File specifier must be "a" or "b"'              & ! !-41  ! discus
     &  /
      DATA ERROR ( -40:  -21) /                     &
     &  'Microdomain type cannot be removed',            & ! !-40  ! discus
     &  'All elements of correlation matrix zero',       & ! !-39  ! discus
     &  'Unsuitable file types for Patterson',           & ! !-38  ! discus
     &  'No filename defined yet            ',           & ! !-37  ! discus
     &  'Unsuitable file types for inverse Fourier',     & ! !-36  ! discus
     &  'Volume of unit cell <= zero',                   & ! !-35  ! discus
     &  'Form does not appear to be closed',             & ! !-34  ! discus
     &  'No microdomain types defined yet',              & ! !-33  ! discus
     &  'Length of vector is zero',                      & ! !-32  ! discus
     &  'Unknown distribution mode',                     & ! !-31  ! discus
     &  'Unknown boundary type',                         & ! !-30  ! discus
     &  'Too many different microdomain types',          & ! !-29  ! discus
     &  'Input parameters must be > zero',               & ! !-28  ! discus
     &  'No atom of this type present in crystal',       & ! !-27  ! discus
     &  'Too many different atoms in crystal',           & ! !-26  ! discus
     &  'Different coordinates in 1-dim files',          & ! !-25  ! discus
     &  'Different coordinates in GNUPLOT files',        & ! !-24  ! discus
     &  'Incompatible GNUPLOT file sizes',               & ! !-23  ! discus
     &  'Incompatible standard file sizes',              & ! !-22  ! discus
     &  'No element present, no Fourier calculated'      & ! !-21  ! discus
     &  /
      DATA ERROR ( -20:   -1) /                     &
     &  'Unknown element, no Fourier calculated',        & ! !-20  ! discus
     &  'Atom number outside limits',                    & ! !-19  ! discus
     &  'No orientation with this number exists',        & ! !-18  ! discus
     &  'Index of matrix outside limits',                & ! !-17  ! discus
     &  'Status for log must be : "on" or "off"',        & ! !-16  ! discus
     &  'No microdomain input file name defined',        & ! !-15  ! discus
     &  'Invalid space group & lattice constants',       & ! !-14  ! discus
     &  'Correlation matrix index outside limits',       & ! !-13  ! discus
     &  'Number of points must be > zero ',              & ! !-12  ! discus
     &  'Unknown threshold type',                        & ! !-11  ! discus
     &  'Too many Atoms in crystal',                     & ! !-10  ! discus
     &  'Unknown Output Format   ',                      & ! ! -9  ! discus
     &  'Too many points in reciprocal layer',           & ! ! -8  ! discus
     &  'Unknown space group symbol         ',           & ! ! -7  ! discus
     &  'Unknown microdomain type',                      & ! ! -6  ! discus
     &  'Too many stacking layers within crystal',       & ! ! -5  ! discus
     &  'Extend of plotspace is zero',                   & ! ! -4  ! discus
     &  'No atoms selected yet',                         & ! ! -3  ! discus
     &  'Improper limits for atom number',               & ! ! -2  ! discus
     &  'Maximum number of Orient. matrizes read'        & ! ! -1  ! discus
     &  /
      DATA ERROR (   0:   io) /                     &
     &  ' ',                                             & ! !  0  ! discus
     &  'Large distance between micro and host',         & ! ! +1  ! discus
     &  'Molecule generator is obsolete >help data',     & ! ! +2  ! discus
     &  'Molecule symmetry  is obsolete >help data',     & ! ! +3  ! discus
     &  'Monte Carlo level is obsolete >help mmc',       & ! ! +4  ! discus
     &  'Element charge was dropped',                    & ! ! +5  ! discus
     &  'Element name is unknown'                        & ! ! +6  ! discus
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
      END SUBROUTINE discus_errlist_appl
!*****7*****************************************************************
      SUBROUTINE errlist_four 
!-                                                                      
!     Displays error Messages for the error type FOUR                    
!+                                                                      
      USE errlist_mod 
      IMPLICIT none 
!                                                                       
!                                                                       
      INTEGER iu, io 
      PARAMETER (IU = -15, IO = 1) 
!                                                                       
      CHARACTER(LEN=45) :: ERROR (IU:IO) 
!                                                                       
      DATA ERROR( IU:-1) /                               &
        ' ',                                             &   ! -15  ! FOUR
        ' ',                                             &   ! -14  ! FOUR
        ' ',                                             &   ! -13  ! FOUR
        ' ',                                             &   ! -12  ! FOUR
        ' ',                                             &   ! -11  ! FOUR
        ' ',                                             &   ! -10  ! FOUR
        ' ',                                             &   ! -9   ! FOUR
        ' ',                                             &   ! -8   ! FOUR
        ' ',                                             &   ! -7   ! FOUR
        ' ',                                             &   ! -6   ! FOUR
        ' ',                                             &   ! -5   ! FOUR
        'Component of increment vector is zero    ',     &   ! -4   ! FOUR
        'SIN(THETA)/LAMBDA > lookup table limits  ',     &   ! -3   ! FOUR
        'Invalid lot shape selected               ',     &   ! -2   ! FOUR
        'Invalid Fourier mode selected            '      &   ! -1   ! FOUR
          /
      DATA ERROR (   0:   io) /                          &
     &  ' ',                                             &   !  0  ! FOUR
     &  'The average intensity sampled zero cells!'      &   ! +1  ! FOUR
          /
!                                                                       
      CALL disp_error ('FOUR', error, iu, io) 
      END SUBROUTINE errlist_four                   
!*****7*****************************************************************
      SUBROUTINE errlist_pdf 
!-                                                                      
!     Displays error Messages for the error type PDF                    
!+                                                                      
      USE errlist_mod 
      IMPLICIT none 
!                                                                       
!                                                                       
      INTEGER iu, io 
      PARAMETER (IU = -11, IO = 0) 
!                                                                       
      CHARACTER(LEN=45) ERROR (IU:IO) 
!                                                                       
      DATA ERROR ( IU: IO) /                             &
          'User Fit maximum outside data range ',        & ! !  -11
          'User Fit minimum outside data range ',        & ! !  -10
          'No atoms in asymmetric unit',                 & ! !  -9
          'Atom type ALL not allowed              ',     & ! !  -8
          'Disable Gaussian mode and recalculate  ',     & ! !  -7
          'PDF range fixed with data loaded ',           & ! !  -6
          'PDF data must start with r=Dr          ',     & ! !  -5
          'No structure defined yet (>= 1 atoms)  ',     & ! !  -4
          'Crystal too large for peridic bound .   ',    & ! !  -3
          'Cannot extend r-range for corr. convol.',     & ! !  -2
          'Too many points in PDF                 ',     & ! !  -1
          ' '                                            & ! !   0
          /                                  
!                                                                       
      CALL disp_error ('PDF ', error, iu, io) 
      END SUBROUTINE errlist_pdf                    
!********************************************************************** 
      SUBROUTINE errlist_rmc 
!-                                                                      
!     Displays error Messages for the error type RMC                    
!+                                                                      
      USE errlist_mod 
      IMPLICIT none 
!                                                                       
!                                                                       
      INTEGER iu, io 
      PARAMETER (iu = -21, io = 0) 
!                                                                       
      CHARACTER(LEN=45) :: ERROR (IU:IO) 
!                                                                       
      DATA ERROR / &
      'Number of LOTS exceeds maximum           ',       &  ! -21
      'No valid move after 1000 disp. intervalls',       &  ! -20
      'Invalid constrain entered                ',       &  ! -19
      'Data file is not an ASCII PGM file       ',       &  ! -18
      'Data and weight file have different sizes',       &  ! -17
      'Invalid weighting scheme / weighting file',       &  ! -16
      'Invalid data type selected               ',       &  ! -15
      'No experimental data within given q limit',       &  ! -14
      'Too many symmetrically equivalent planes ',       &  ! -13
      'Displacements too small for SWDISP mode  ',       &  ! -12
      'Only ONE atom type present in SWCHEM mode',       &  ! -11
      'No atom types selected for RMC run       ',       &  ! -10
      'Invalid RMC/MC mode selected             ',       &  !  -9
      'Invalid symmetry number selected         ',       &  !  -8
      'Invalid plane selected                   ',       &  !  -7
      'Too many atoms per molecule for RMC      ',       &  !  -6
      'Invalid method (x,n) selected            ',       &  !  -5
      'No experimental data present             ',       &  !  -4
      'No atoms present in model crystal        ',       &  !  -3
      'Too many experimental data points        ',       &  !  -2
      'Too many experimental data planes        ',       &  !   1
      ' ' /                                                 !   0
!                                                                       
      CALL disp_error ('RMC ', error, iu, io) 
      END SUBROUTINE errlist_rmc                    
!*****7*****************************************************************
      SUBROUTINE errlist_mmc 
!-                                                                      
!     Displays error Messages for the error type MC                     
!+                                                                      
      USE errlist_mod 
      IMPLICIT none 
!                                                                       
!                                                                       
      INTEGER iu, io 
      PARAMETER (IU = -5, IO = 0) 
!                                                                       
      CHARACTER(LEN=45) ::ERROR (IU:IO) 
!                                                                       
      DATA ERROR /                                       &
      'Number of feedback intervalls is zero    ',       & ! !  -5
      'Number of MC cycles is zero              ',       & ! !  -4
      'Invalid mode selected for COCC MC run    ',       & ! !  -3
      'No valid move after 1000 cycles          ',       & ! !  -2
      'Invalid or no energy type selected       ',       & ! !  -1
      ' ' /                                                !   0
!                                                                       
      CALL disp_error ('MMC ', error, iu, io) 
      END SUBROUTINE errlist_mmc                    
!*****7**************************************************************** 
      SUBROUTINE errlist_chem 
!-                                                                      
!     Displays error messages for the error type CHEM                   
!+                                                                      
      USE errlist_mod 
      IMPLICIT none 
!                                                                       
!                                                                       
      INTEGER iu, io 
      PARAMETER (IU = - 31, IO = 0) 
!                                                                       
      CHARACTER(LEN=45) :: ERROR (IU:IO) 
!                                                                       
      DATA ERROR ( IU:-21) /                             &
      'Periodic boundries disabled, PURGE was done ',    & ! !-31
      'Multiple identical sites in unit cell   ',        & ! !-30
      'Atom type outside valid range           ',        & ! !-29
      'Invalid correlation conn   index given  ',        & ! !-28
      'No atoms present in crystal             ',        & ! !-27
      'Invalid correlation environment index   ',        & ! !-26
      'Invalid range for bond-angle histogramm ',        & ! !-25
      'Invalid correlation angle index given   ',        & ! !-24
      'Too many neighbouring atoms/molecules   ',        & ! !-23
      'Command not available in molecule mode  ',        & ! !-22
      'Molecule types need to be different     '         & ! !-21
      /                                 
      DATA ERROR (-20: -1) /                             &
      'No molecules present in crystal         ',        & ! !-20
      'No neighbouring molecules found         ',        & ! !-19
      'Correlation fields require same # vector',        & ! !-18
      'Correlation fields require same mode    ',        & ! !-17
      'Failed to apply periodic boundaries     ',        & ! !-16
      'No displacement directions selected     ',        & ! !-15
      'Invalid neighbour definition selected   ',        & ! !-14
      'Correlation direction invalid           ',        & ! !-13
      'Too many neighbour definitions          ',        & ! !-12
      'No neighbours defined                   ',        & ! !-11
      'Invalid crystal site or atom index given',        & ! !-10
      'Invalid correlation vector index given  ',        & ! !- 9
      'No neighbouring atoms found             ',        & ! !- 8
      'Atoms need to be different              ',        & ! !- 7
      'Atom name ALL not allowed for command   ',        & ! !- 6
      'Too many different atoms found          ',        & ! !- 5
      'Invalid SIGMA entered                   ',        & ! !- 4
      'Invalid range for bond-length histogramm',        & ! !- 3
      'Not enough space for all result in res[]',        & ! !- 2
      'Too many points for histogramm          '         & ! !- 1
      /                                 
      DATA ERROR (  0:  0) /                             &
      ' '                                                & ! !  0
      /                                 
!                                                                       
      CALL disp_error ('CHEM', error, iu, io) 
      END SUBROUTINE errlist_chem                   
