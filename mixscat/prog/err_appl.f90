!*****7****************************************************************
!
      subroutine mixscat_errlist_appl
!-
!     Displays error Messages for the error type APPLication
!+
      USE errlist_mod 
      implicit      none
!
!
      integer       iud,iod,iup,iop
      parameter    (iud=-54,iod=0)
      parameter    (iup=-18,iop=0)
!
      character(LEN=45) ::  error_discus(iud:iod)
      character(LEN=45) ::  error_pdffit(iup:iop)
!
      data error_pdffit (-18: -1) /                 &
     &  'Gaussian too broad for computation       ',&!-18  ! pdffit
     &  'Can not appy periodic boundaries         ',&!-17  ! pdffit
     &  'Invalid r-range specified                ',&!-16  ! pdffit
     &  'Observed and reference PDF do not match  ',&!-15  ! pdffit
     &  'Parameter index in p[i] out of range     ',&!-14  ! pdffit
     &  'R-value out of range                     ',&!-13  ! pdffit
     &  'Cannot extend r-range for convolution    ',&!-12  ! pdffit
     &  'Too many data points in experimental PDF ',&!-11  ! pdffit
     &  'Invalid parameter selected               ',&!-10  ! pdffit
     &  'Inconsistency between NCELL and # atoms  ',&! -9  ! pdffit
     &  'Too many different parameters in constr. ',&! -8  ! pdffit
     &  'Invalid parameter constraint specified   ',&! -7  ! pdffit
     &  'Invalid data set specified               ',&! -6  ! pdffit
     &  'Invalid radiation type selected          ',&! -5  ! pdffit
     &  'Maximum number of data sets exceeded     ',&! -4  ! pdffit
     &  'Invalid occupancy specified              ',&! -3  ! pdffit
     &  'Number of phases outside limits          ',&! -2  ! pdffit
     &  'Invalid structure phase selected         ' &! -1  ! pdffit
        /
      data error_pdffit (  0:  0) /                 &
     &  ' '                                         &!  0  ! pdffit
     &       /
!
      data error_discus (-54:-41)/                  &
     &  'Index outside limits                     ',&!-54  ! discus
     &  ' ',                                        &!-53  ! unused
     &  ' ',                                        &!-52  ! unused
     &  ' ',                                        &!-51  ! unused
     &  ' ',                                        &!-50  ! unused
     &  'Error reading atom coordinates           ',&!-49  ! discus
     &  'Error reading lattice constants          ',&!-48  ! discus
     &  'Error reading space group symbol         ',&!-47  ! discus
     &  'Error reading title of structure         ',&!-46  ! discus
     &  'Too many atoms in environment',            &!-45  ! discus
     &  'Right quotation mark missing in format',   &!-44  ! discus
     &  'Not enough parameter for filename format ',&!-43  ! discus
     &  ' ',                                        &!-42  ! unused
     &  ' '                                         &!-41  ! unused
        /
      data error_discus (-40:-21)/                  &
     &  ' ',                                        &!-40  ! unused
     &  ' ',                                        &!-39  ! unused
     &  ' ',                                        &!-38  ! unused
     &  'No filename defined yet            ',      &!-37  ! discus
     &  ' ',                                        &!-36  ! unused
     &  'Volume of unit cell <= zero',              &!-35  ! discus
     &  ' ',                                        &!-34  ! unused
     &  ' ',                                        &!-33  ! unused
     &  'Length of vector is zero',                 &!-32  ! discus
     &  ' ',                                        &!-31  ! unused
     &  ' ',                                        &!-31  ! unused
     &  ' ',                                        &!-29  ! unused
     &  'Input parameters must be > zero',          &!-28  ! discus
     &  'No atom of this type present in crystal',  &!-27  ! discus
     &  'Too many different atoms in crystal',      &!-26  ! discus
     &  ' ',                                        &!-25  ! unused
     &  ' ',                                        &!-24  ! unused
     &  ' ',                                        &!-23  ! unused
     &  ' ',                                        &!-22  ! unused
     &  'No element present, no weights calculated' &!-21  ! discus
        /
      data error_discus (-20: -1)/                  &
     &  'Unknown element, no weights calculated',   &!-20  ! discus
     &  'Atom number outside limits',               &!-19  ! discus
     &  ' ',                                        &!-18  ! unused
     &  ' ',                                        &!-17  ! unused
     &  ' ',                                        &!-16  ! unused
     &  ' ',                                        &!-15  ! unused
     &  'Invalid space group & lattice constants',  &!-14  ! discus
     &  ' ',                                        &!-13  ! unused
     &  'Number of points must be > zero ',         &!-12  ! discus
     &  'Run command calc first',                   &!-11  ! mixsca
     &  'Too many Atoms in crystal',                &!-10  ! discus
     &  'Rmin or Dr of datasets do not match',      &! -9  ! mixsca
     &  'At least 2 data sets needed',              &! -8  ! mixsca
     &  'Unknown space group symbol         ',      &! -7  ! discus
     &  'No atom pair to remove defined',           &! -6  ! mixsca
     &  'No sample composition defined',            &! -5  ! mixsca
     &  'Element not present in sample',            &! -4  ! mixsca
     &  'No atoms selected yet',                    &! -3  ! discus
     &  'Improper limits for atom number',          &! -2  ! discus
     &  'Too many elements in composition'          &! -1  ! mixsca  
        /
      data error_discus (  0:  0)/                  &
     &  ' '                                         &!  0  ! discus
     &       /
!
      if (ier_num.le.-100) then
        ier_num = ier_num + 100
        call disp_error ('APPL',error_pdffit,iup,iop)
      else
        call disp_error ('APPL',error_discus,iud,iod)
      endif
!
      end subroutine mixscat_errlist_appl
