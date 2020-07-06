!*****7****************************************************************
      SUBROUTINE kuplot_errlist_appl
!-
!     Displays error Messages for the error type APPLication
!+
      USE errlist_mod 
USE lib_errlist_func
      implicit      none
!
!
      integer       iu,io
      parameter    (iu=-72,io=0)
!
      character(LEN=45) ::  error(iu:io)
!
      data error ( iu: -61) /                      &
     &  'HDF5 is not supported in this version',   &! -72  ! kupl
     &  'Optional layer number invalid',           &! -71  ! kupl
     &  'Trouble reading dataset list   ',         &! -70  ! kupl
     &  'Not a Yell/DISCUS h5 file      ',         &! -69  ! kupl
     &  'Sigma is zero with weight = DAT',         &! -68  ! kupl
     &  'Error reading GENERATION file',           &! -67  ! kupl
     &  'DIFFEV Rvalue file not found',            &! -66  ! kupl
     &  'DIFFEV LOG  or Current file not found',   &! -65  ! kupl
     &  'DIFFEV Summary file not found',           &! -64  ! kupl
     &  'Parameter no outside range',              &! -63  ! kupl
     &  'Data set is not 2D',                      &! -62  ! kupl
     &  'Singular matrix in Savitzky calculation ' &! -61  ! kupl
     &  /
      data error (-60: -41) /                      &
     &  'Multiple point with same x, no spline   ',&! -60  ! kupl
     &  'Invalid column number specified         ',&! -59  ! kupl
     &  'SDS section value out of range          ',&! -58  ! kupl
     &  'Invalid SDS name specified (try nxdir)  ',&! -57  ! kupl
     &  'NeXus SDS has too many dimensions       ',&! -56  ! kupl
     &  'No NeXus file currently open            ',&! -55  ! kupl
     &  'Close current NeXus file first          ',&! -54  ! kupl
     &  'This KUPLOT has NeXus support disabled  ',&! -53  ! kupl
     &  'Unsupported incident spectr. function   ',&! -52  ! kupl
     &  'Unit for GSAS TOF conversion invalid    ',&! -51  ! kupl
     &  'Requested bank not found in GSAS file   ',&! -50  ! kupl
     &  'No TIME_MAP entry found in GSAS file    ',&! -49  ! kupl
     &  'Unsupported GSAS binning type found     ',&! -48  ! kupl
     &  'Error reading instrument parameter file ',&! -47  ! kupl
     &  'Invalid bank number found in iparm. file',&! -46  ! kupl
     &  'Step size of data sets is different     ',&! -45  ! kupl
     &  'Data sets must have same x-range        ',&! -44  ! kupl
     &  'Invalid value/range for log. axis       ',&! -43  ! kupl
     &  'Invalid window ID selected              ',&! -42  ! kupl
     &  'MCA scan not found                      ' &! -41  ! kupl
     &  /
      data error (-40: -21) /                      &
     &  'Maximum derivative exceeded             ',&! -40  ! kupl
     &  'Invalid smoothing size specified        ',&! -39  ! kupl
     &  'Invalid column value specified          ',&! -38  ! kupl
     &  'Scan not found in SPEC file             ',&! -37  ! kupl
     &  'Configuration mismatch file & KUPLOT    ',&! -36  ! kupl
     &  'No valid frame selected                 ',&! -35  ! kupl
     &  'Data set to large for bitmap drawing    ',&! -34  ! kupl
     &  'Invalid bond definition selected        ',&! -33  ! kupl
     &  'Error in KUPLOT contour routines        ',&! -32  ! kupl
     &  'To many fit parameters                  ',&! -31  ! kupl
     &  'Not enough maxima found for start values',&! -30  ! kupl
     &  'x or y-range of data set is zero        ',&! -29  ! kupl
     &  'Invalid annotation number selected      ',&! -28  ! kupl
     &  'Invalid weighting scheme selected       ',&! -27  ! kupl
     &  'Invalid fit parameter selected          ',&! -26  ! kupl
     &  'Invalid or no fit function selected     ',&! -25  ! kupl
     &  'Size of res[] array exceeded            ',&! -24  ! kupl
     &  'Incompatible data sets for KCAL         ',&! -23  ! kupl
     &  'No maxima found                         ',&! -22  ! kupl
     &  'Too many maxima found at search         ' &! -21  ! kupl
     &  /
      data error (-20:  -1) /                      &
     &  'No data points in selected area         ',&! -20  ! kupl
     &  'Invalid peak number found               ',&! -19  ! kupl
     &  'Invalid RGB color found                 ',&! -18  ! kupl
     &  'Maximum number of frames exceeded       ',&! -17  ! kupl
     &  'Invalid X11 screen width entered        ',&! -16  ! kupl
     &  'Invalid frame selected                  ',&! -15  ! kupl
     &  'Invalid contour line set selected       ',&! -14  ! kupl
     &  'Too many points for spline / steps      ',&! -13  ! kupl
     &  'No data present to plot                 ',&! -12  ! kupl
     &  'Invalid plot device selected            ',&! -11  ! kupl
     &  'Too many major tick marks               ',&! -10  ! kupl
     &  'Too many excluded regions               ',&!  -9  ! kupl
     &  'Use only if NO data set is loaded       ',&!  -8  ! kupl
     &  'Invalid parameter value entered         ',&!  -7  ! kupl
     &  'Number of data points exceeds limit     ',&!  -6  ! kupl
     &  'Not an ASCII PGM file                   ',&!  -5  ! kupl
     &  'Invalid data set selected               ',&!  -4  ! kupl
     &  'Maximum number of data sets exceeded    ',&!  -3  ! kupl
     &  'Unkown file format                      ',&!  -2  ! kupl
     &  'Maximum number of data sets exceeded    ' &!  -1  ! kupl
     &  /
      data error (  0:   0) /                      &
     &  ' '                                        &!   0  ! kupl
     &           /
!
      call disp_error ('KUPL',error,iu,io)
      END SUBROUTINE kuplot_errlist_appl
!*****7****************************************************************
