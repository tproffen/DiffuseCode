MODULE kuplot_words_mod
!
IMPLICIT NONE
!
!INTEGER, PARAMETER :: user         = -6
!INTEGER, PARAMETER :: kupl         = -5
!INTEGER, PARAMETER :: thermal      = -4
!INTEGER, PARAMETER :: pdf          = -3
!INTEGER, PARAMETER :: ice          = -2
!INTEGER, PARAMETER :: fire         = -1
INTEGER, PARAMETER :: map          = -1
INTEGER, PARAMETER :: none         =  0
INTEGER, PARAMETER :: red          =  1
INTEGER, PARAMETER :: green        =  2
INTEGER, PARAMETER :: blue         =  3
INTEGER, PARAMETER :: magenta      =  4
INTEGER, PARAMETER :: yellow       =  5
INTEGER, PARAMETER :: black        =  6
INTEGER, PARAMETER :: dark_red     =  7
INTEGER, PARAMETER :: dark_green   =  8
INTEGER, PARAMETER :: dark_blue    =  9
INTEGER, PARAMETER :: dark_magenta = 10
INTEGER, PARAMETER :: dark_yellow  = 11
INTEGER, PARAMETER :: gray         = 12
INTEGER, PARAMETER :: cyan         = 13
INTEGER, PARAMETER :: dark_cyan    = 14
INTEGER, PARAMETER :: white        = 15
!
INTEGER, PARAMETER :: nomarker     =  0
INTEGER, PARAMETER :: dot          =  1
INTEGER, PARAMETER :: circle       =  2
INTEGER, PARAMETER :: filled_circle=  3
INTEGER, PARAMETER :: square       =  4
INTEGER, PARAMETER :: filled_square=  5
INTEGER, PARAMETER :: triangle     =  6
INTEGER, PARAMETER :: times        =  7
INTEGER, PARAMETER :: plus         =  8
INTEGER, PARAMETER :: bar          =  9
INTEGER, PARAMETER :: slash        = 10
INTEGER, PARAMETER :: backslash    = 11
INTEGER, PARAMETER :: minus        = 12
INTEGER, PARAMETER :: vertical     = 13
INTEGER, PARAMETER :: x_coord      = -1
INTEGER, PARAMETER :: y_coord      = -2
INTEGER, PARAMETER :: xy_coord     = -3
!
INTEGER, PARAMETER :: noline       =  0
INTEGER, PARAMETER :: solid        =  1
INTEGER, PARAMETER :: broken       =  2
INTEGER, PARAMETER :: dashed       =  3
INTEGER, PARAMETER :: short_broken =  4
INTEGER, PARAMETER :: short_dashed =  5
!
INTEGER, PARAMETER :: noerror      =  0
INTEGER, PARAMETER :: xbar         =  1
INTEGER, PARAMETER :: ybar         =  2
INTEGER, PARAMETER :: cross        =  3
!
INTEGER, PARAMETER :: fill_solid   =  1
INTEGER, PARAMETER :: fill_slash   =  2
INTEGER, PARAMETER :: fill_back    =  3
INTEGER, PARAMETER :: fill_cross   =  4
INTEGER, PARAMETER :: b_fill_solid =  5
INTEGER, PARAMETER :: b_fill_slash =  6
INTEGER, PARAMETER :: b_fill_back  =  7
INTEGER, PARAMETER :: b_fill_cross =  8
!
CONTAINS
   SUBROUTINE get_words (ianz, cpara, lpara, maxw, item, bef)
!
!  REPLACE the string in cpara by the appropriate number 
!  If a word is not recognized no error is specified as its likely 
!  a numerical value
USE str_comp_mod
   IMPLICIT NONE
!
   INTEGER                           , INTENT(IN)    :: maxw
   INTEGER                           , INTENT(IN)    :: ianz
   CHARACTER (LEN=*), DIMENSION(maxw), INTENT(INOUT) :: cpara
   INTEGER          , DIMENSION(maxw), INTENT(INOUT) :: lpara
   INTEGER                           , INTENT(IN)    :: item
   CHARACTER (LEN=*)                 , INTENT(IN)    :: bef
!
   INTEGER :: lbef
!
   lbef = LEN_TRIM(bef)
   IF(str_comp(bef, 'lcolor', 3, lbef, 6 ) .OR.  &
      str_comp(bef, 'hcolor', 3, lbef, 6 ) .OR.  &
      str_comp(bef, 'mcolor', 3, lbef, 6 ) .OR.  &
      str_comp(bef, 'fcolor', 3, lbef, 6 ) .OR.  &
      str_comp(bef, 'ecolor', 3, lbef, 6 )     ) THEN
      CALL get_colors (ianz, cpara, lpara, maxw, item)
   ELSEIF(str_comp(bef, 'mtype', 3, lbef, 5 ) .OR. &
          str_comp(bef, 'ptype', 3, lbef, 5 ))  THEN
      CALL get_marker (ianz, cpara, lpara, maxw, item)
   ELSEIF(str_comp(bef, 'ltype', 3, lbef, 5 ) .OR. &
          str_comp(bef, 'htype', 3, lbef, 5 ))  THEN
      CALL get_linetype(ianz, cpara, lpara, maxw, item)
   ELSEIF(str_comp(bef, 'etype', 3, lbef, 5 ))  THEN 
      CALL get_err_type(ianz, cpara, lpara, maxw, item)
   ELSEIF(str_comp(bef, 'fill', 3, lbef, 4 ))  THEN 
      CALL get_filltype(ianz, cpara, lpara, maxw, item)
   ENDIF
!
   END SUBROUTINE get_words
!
   SUBROUTINE get_colors (ianz, cpara, lpara, maxw, item)
!
!  REPLACE the string in cpara by the appropriate number 
!  If a word is not recognized no error is specified as its likely 
!  a numerical value
USE str_comp_mod
      USE string_convert_mod
   IMPLICIT NONE
!
   INTEGER                           , INTENT(IN)    :: maxw
   INTEGER                           , INTENT(IN)    :: ianz
   CHARACTER (LEN=*), DIMENSION(maxw), INTENT(INOUT) :: cpara
   INTEGER          , DIMENSION(maxw), INTENT(INOUT) :: lpara
   INTEGER                           , INTENT(IN)    :: item
!
!
   CALL do_cap(cpara(item))   ! Lets not worry about capitalization
!  IF(str_comp (cpara(item),     'USER', 4, lpara(item) , 4)) THEN
!      WRITE(cpara(item),'(i3)') user
!      lpara(item) = 3
!  ELSEIF(str_comp (cpara(item),     'KUPL', 4, lpara(item) , 4)) THEN
!      WRITE(cpara(item),'(i3)') kupl
!      lpara(item) = 3
!  ELSEIF(str_comp (cpara(item),     'THERMAL', 4, lpara(item) , 7)) THEN
!      WRITE(cpara(item),'(i3)') thermal
!      lpara(item) = 3
!  ELSEIF(str_comp (cpara(item),     'PDF' , 2, lpara(item) , 4)) THEN
!      WRITE(cpara(item),'(i3)') pdf 
!      lpara(item) = 3
!  ELSEIF(str_comp (cpara(item),     'ICE' , 2, lpara(item) , 4)) THEN
!      WRITE(cpara(item),'(i3)') ice 
!      lpara(item) = 3
!  ELSEIF(str_comp (cpara(item),     'FIRE', 2, lpara(item) , 4)) THEN
!      WRITE(cpara(item),'(i3)') fire
!      lpara(item) = 3
       IF(str_comp (cpara(item),     'MAP' , 3, lpara(item) , 3)) THEN
       WRITE(cpara(item),'(i3)') map
       lpara(item) = 3
   ELSEIF(str_comp (cpara(item),     'NONE', 2, lpara(item) , 4)) THEN
       WRITE(cpara(item),'(i3)') none
       lpara(item) = 3
   ELSEIF(str_comp (cpara(item), 'RED', 2, lpara(item) , 3)) THEN
       WRITE(cpara(item),'(i3)') red
       lpara(item) = 3
   ELSEIF(str_comp (cpara(item), 'GREEN', 3, lpara(item) , 5)) THEN
       WRITE(cpara(item),'(i3)') green
       lpara(item) = 3
   ELSEIF(str_comp (cpara(item), 'BLUE' , 3, lpara(item) , 4)) THEN
       WRITE(cpara(item),'(i3)') blue
       lpara(item) = 3
   ELSEIF(str_comp (cpara(item), 'MAGENTA', 2, lpara(item) , 7)) THEN
       WRITE(cpara(item),'(i3)') magenta
       lpara(item) = 3
   ELSEIF(str_comp (cpara(item), 'YELLOW', 2, lpara(item) , 6)) THEN
       WRITE(cpara(item),'(i3)') yellow
       lpara(item) = 3
   ELSEIF(str_comp (cpara(item), 'BLACK', 3, lpara(item) , 5)) THEN
       WRITE(cpara(item),'(i3)') black
       lpara(item) = 3
   ELSEIF(str_comp (cpara(item), 'DARKRED', 5, lpara(item) , 7)) THEN
       WRITE(cpara(item),'(i3)') dark_red
       lpara(item) = 3
   ELSEIF(str_comp (cpara(item), 'DARKGREEN', 5, lpara(item) , 9)) THEN
       WRITE(cpara(item),'(i3)') dark_green
       lpara(item) = 3
   ELSEIF(str_comp (cpara(item), 'DARKBLUE', 5, lpara(item) , 8)) THEN
       WRITE(cpara(item),'(i3)') dark_blue
       lpara(item) = 3
   ELSEIF(str_comp (cpara(item), 'DARKMAGENTA', 5, lpara(item) ,11)) THEN
       WRITE(cpara(item),'(i3)') dark_magenta
       lpara(item) = 3
   ELSEIF(str_comp (cpara(item), 'DARKYELLOW', 5, lpara(item) ,10)) THEN
       WRITE(cpara(item),'(i3)') dark_yellow
       lpara(item) = 3
   ELSEIF(str_comp (cpara(item), 'GRAY', 3, lpara(item) , 4)) THEN
       WRITE(cpara(item),'(i3)') gray
       lpara(item) = 3
   ELSEIF(str_comp (cpara(item), 'CYAN', 2, lpara(item) , 4)) THEN
       WRITE(cpara(item),'(i3)') cyan
       lpara(item) = 3
   ELSEIF(str_comp (cpara(item), 'DARKCYAN', 5, lpara(item) , 8)) THEN
       WRITE(cpara(item),'(i3)') dark_cyan
       lpara(item) = 3
   ELSEIF(str_comp (cpara(item), 'WHITE', 2, lpara(item) , 5)) THEN
       WRITE(cpara(item),'(i3)') white
       lpara(item) = 3
   ENDIF
   END SUBROUTINE get_colors
!
   SUBROUTINE get_marker (ianz, cpara, lpara, maxw, item)
!
!  REPLACE the string in cpara by the appropriate number 
!  If a word is not recognized no error is specified as its likely 
!  a numerical value
USE str_comp_mod
      USE string_convert_mod
   IMPLICIT NONE
!
   INTEGER                           , INTENT(IN)    :: maxw
   INTEGER                           , INTENT(IN)    :: ianz
   CHARACTER (LEN=*), DIMENSION(maxw), INTENT(INOUT) :: cpara
   INTEGER          , DIMENSION(maxw), INTENT(INOUT) :: lpara
   INTEGER                           , INTENT(IN)    :: item
!
!
   CALL do_cap(cpara(item))   ! Lets not worry about capitalization
   IF(str_comp (cpara(item),     'NONE', 2, lpara(item) , 4) .OR.  &
      str_comp (cpara(item), 'NOMARKER', 2, lpara(item) , 8)) THEN
       WRITE(cpara(item),'(i3)') nomarker
       lpara(item) = 3
   ELSEIF(str_comp (cpara(item), 'DOT'   , 2, lpara(item) , 3)) THEN
       WRITE(cpara(item),'(i3)') dot
       lpara(item) = 3
   ELSEIF(str_comp (cpara(item), 'CIRCLE', 3, lpara(item) , 6)) THEN
       WRITE(cpara(item),'(i3)') circle
       lpara(item) = 3
   ELSEIF(str_comp (cpara(item), 'FILLEDCIRCLE' , 7, lpara(item) ,12)) THEN
       WRITE(cpara(item),'(i3)') filled_circle
       lpara(item) = 3
   ELSEIF(str_comp (cpara(item), 'SQUARE', 2, lpara(item) , 6)) THEN
       WRITE(cpara(item),'(i3)') square 
       lpara(item) = 3
   ELSEIF(str_comp (cpara(item), 'FILLEDSQUARE', 7, lpara(item) ,12)) THEN
       WRITE(cpara(item),'(i3)') filled_square
       lpara(item) = 3
   ELSEIF(str_comp (cpara(item), 'TRIANGLE', 2, lpara(item) , 8)) THEN
       WRITE(cpara(item),'(i3)') triangle
       lpara(item) = 3
   ELSEIF(str_comp (cpara(item), 'TIMES', 2, lpara(item) , 5)) THEN
       WRITE(cpara(item),'(i3)') times
       lpara(item) = 3
   ELSEIF(str_comp (cpara(item), 'PLUS', 2, lpara(item) , 4)) THEN
       WRITE(cpara(item),'(i3)') plus
       lpara(item) = 3
   ELSEIF(str_comp (cpara(item), 'BAR', 3, lpara(item) , 3)) THEN
       WRITE(cpara(item),'(i3)') bar
       lpara(item) = 3
   ELSEIF(str_comp (cpara(item), 'SLASH', 2, lpara(item) , 5)) THEN
       WRITE(cpara(item),'(i3)') slash
       lpara(item) = 3
   ELSEIF(str_comp (cpara(item), 'BACKSLASH', 3, lpara(item) , 9)) THEN
       WRITE(cpara(item),'(i3)') backslash
       lpara(item) = 3
   ELSEIF(str_comp (cpara(item), 'MINUS', 2, lpara(item) , 5)) THEN
       WRITE(cpara(item),'(i3)') minus
       lpara(item) = 3
   ELSEIF(str_comp (cpara(item), 'VERTICAL', 2, lpara(item) , 8)) THEN
       WRITE(cpara(item),'(i3)') vertical
       lpara(item) = 3
   ELSEIF(str_comp (cpara(item), 'X_COORD' , 2, lpara(item) , 7)) THEN
       WRITE(cpara(item),'(i2)') x_coord 
       lpara(item) = 2
   ELSEIF(str_comp (cpara(item), 'Y_COORD' , 2, lpara(item) , 7)) THEN
       WRITE(cpara(item),'(i2)') y_coord 
       lpara(item) = 2
   ELSEIF(str_comp (cpara(item), 'XY'      , 2, lpara(item) , 2)) THEN
       WRITE(cpara(item),'(i2)') xy_coord 
       lpara(item) = 2
   ENDIF
   END SUBROUTINE get_marker
!
   SUBROUTINE get_linetype (ianz, cpara, lpara, maxw, item)
!
!  REPLACE the string in cpara by the appropriate number 
!  If a word is not recognized no error is specified as its likely 
!  a numerical value
USE str_comp_mod
      USE string_convert_mod
   IMPLICIT NONE
!
   INTEGER                           , INTENT(IN)    :: maxw
   INTEGER                           , INTENT(IN)    :: ianz
   CHARACTER (LEN=*), DIMENSION(maxw), INTENT(INOUT) :: cpara
   INTEGER          , DIMENSION(maxw), INTENT(INOUT) :: lpara
   INTEGER                           , INTENT(IN)    :: item
!
!
   CALL do_cap(cpara(item))   ! Lets not worry about capitalization
   IF(str_comp (cpara(item),     'NONE', 2, lpara(item) , 4) .OR.  &
      str_comp (cpara(item),   'NOLINE', 2, lpara(item) , 6)) THEN
       WRITE(cpara(item),'(i3)') noline
       lpara(item) = 3
   ELSEIF(str_comp (cpara(item), 'SOLID'   , 2, lpara(item) , 5)) THEN
       WRITE(cpara(item),'(i3)') solid
       lpara(item) = 3
   ELSEIF(str_comp (cpara(item), 'BROKEN', 2, lpara(item) , 6)) THEN
       WRITE(cpara(item),'(i3)') broken
       lpara(item) = 3
   ELSEIF(str_comp (cpara(item), 'DASHED', 2, lpara(item) , 6)) THEN
       WRITE(cpara(item),'(i3)') dashed
       lpara(item) = 3
   ELSEIF(str_comp (cpara(item), 'SHORTBROKEN' , 7, lpara(item) ,11)) THEN
       WRITE(cpara(item),'(i3)') short_broken
       lpara(item) = 3
   ELSEIF(str_comp (cpara(item), 'SHORTDASHED', 7, lpara(item) ,11)) THEN
       WRITE(cpara(item),'(i3)') short_dashed
       lpara(item) = 3
   ENDIF
   END SUBROUTINE get_linetype
!
   SUBROUTINE get_err_type (ianz, cpara, lpara, maxw, item)
!
!  REPLACE the string in cpara by the appropriate number 
!  If a word is not recognized no error is specified as its likely 
!  a numerical value
USE str_comp_mod
      USE string_convert_mod
   IMPLICIT NONE
!
   INTEGER                           , INTENT(IN)    :: maxw
   INTEGER                           , INTENT(IN)    :: ianz
   CHARACTER (LEN=*), DIMENSION(maxw), INTENT(INOUT) :: cpara
   INTEGER          , DIMENSION(maxw), INTENT(INOUT) :: lpara
   INTEGER                           , INTENT(IN)    :: item
!
!
   CALL do_cap(cpara(item))   ! Lets not worry about capitalization
   IF(str_comp (cpara(item),     'NONE', 2, lpara(item) , 4) .OR.  &
      str_comp (cpara(item),  'NOERROR', 2, lpara(item) , 7)) THEN
       WRITE(cpara(item),'(i3)') noline
       lpara(item) = 3
   ELSEIF(str_comp (cpara(item), 'XBAR', 2, lpara(item) , 4)) THEN
       WRITE(cpara(item),'(i3)') xbar
       lpara(item) = 3
   ELSEIF(str_comp (cpara(item), 'YBAR', 2, lpara(item) , 4)) THEN
       WRITE(cpara(item),'(i3)') ybar
       lpara(item) = 3
   ELSEIF(str_comp (cpara(item), 'CROSS', 2, lpara(item) , 5)) THEN
       WRITE(cpara(item),'(i3)') cross
       lpara(item) = 3
   ENDIF
   END SUBROUTINE get_err_type
!
   SUBROUTINE get_filltype (ianz, cpara, lpara, maxw, item)
!
!  REPLACE the string in cpara by the appropriate number 
!  If a word is not recognized no error is specified as its likely 
!  a numerical value
USE str_comp_mod
      USE string_convert_mod
   IMPLICIT NONE
!
   INTEGER                           , INTENT(IN)    :: maxw
   INTEGER                           , INTENT(IN)    :: ianz
   CHARACTER (LEN=*), DIMENSION(maxw), INTENT(INOUT) :: cpara
   INTEGER          , DIMENSION(maxw), INTENT(INOUT) :: lpara
   INTEGER                           , INTENT(IN)    :: item
!
!
   CALL do_cap(cpara(item))   ! Lets not worry about capitalization
       IF(str_comp (cpara(item), 'SOLID', 2, lpara(item) , 5)) THEN
       WRITE(cpara(item),'(i3)') fill_solid
       lpara(item) = 3
   ELSEIF(str_comp (cpara(item), 'SLASH', 2, lpara(item) , 5)) THEN
       WRITE(cpara(item),'(i3)') fill_slash
       lpara(item) = 3
   ELSEIF(str_comp (cpara(item), 'BACKSLASH', 2, lpara(item) , 9)) THEN
       WRITE(cpara(item),'(i3)') fill_back
       lpara(item) = 3
   ELSEIF(str_comp (cpara(item), 'CROSS', 2, lpara(item) , 5)) THEN
       WRITE(cpara(item),'(i3)') fill_cross
       lpara(item) = 3
   ELSEIF(str_comp (cpara(item), 'BSOLID', 3, lpara(item) , 6)) THEN
       WRITE(cpara(item),'(i3)') b_fill_solid
       lpara(item) = 3
   ELSEIF(str_comp (cpara(item), 'BSLASH', 3, lpara(item) , 6)) THEN
       WRITE(cpara(item),'(i3)') b_fill_slash
       lpara(item) = 3
   ELSEIF(str_comp (cpara(item), 'BBACKSLASH', 2, lpara(item) ,10)) THEN
       WRITE(cpara(item),'(i3)') b_fill_back
       lpara(item) = 3
   ELSEIF(str_comp (cpara(item), 'BCROSS', 2, lpara(item) , 6)) THEN
       WRITE(cpara(item),'(i3)') b_fill_cross
       lpara(item) = 3
   ENDIF
   END SUBROUTINE get_filltype
!
END MODULE kuplot_words_mod
