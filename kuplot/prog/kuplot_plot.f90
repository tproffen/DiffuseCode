module kuplot_plot_mod
!
!********************************************************************   
!     These are the main plot routines for 2D and 3D data. 3D data      
!     can be plotted as contour plots or bitmaps or both. All           
!     settings are defined in 'kuplot.inc'.                             
!
!     Modified to just do_hardcopy
!*********************************************************************  
!
contains
!
!*******************************************************************************
!
SUBROUTINE do_plot (lmenu) 
!                                                                       
!     Main plotting routine                                             
!                                                                       
      USE errlist_mod 
      USE kuplot_config 
      USE kuplot_mod 
use kuplot_draw_low_mod
use kuplot_plot_low_mod
!                                                                       
      IMPLICIT none 
!                                                                       
      REAL x1, x2, y1, y2, width, ratio, dev_sf_old 
      INTEGER ii
      LOGICAL tfr, lmenu 
!                                                                       
      INTEGER PGOPEN 
!                                                                       
      DATA dev_sf_old / 0.7 / 
!                                                                       
!------ see if there are data at all                                    
!                                                                       
      tfr = .false. 
      DO ii = 1, iaf (iwin) 
      tfr = tfr.or. (infra (iwin, ii, 1) .eq. - 1) 
      ENDDO 
!                                                                       
      IF (iz.le.1.and..not.tfr) then 
         ier_num = - 12 
         ier_typ = ER_APPL 
         RETURN 
      ENDIF 
!                                                                       
!------ Open screen graphics device                                     
!                                                                       
      IF (lmenu) then 
         dev_draw (iwin, 1) = 0.20 
         dev_draw (iwin, 2) = 0.15 
      ELSE 
         dev_draw (iwin, 1) = 0.00 
         dev_draw (iwin, 2) = 0.00 
      ENDIF 
!                                                                       
      IF (dev_id (iwin, x11) .eq. - 1) then 
         dev_id (iwin, x11) = PGOPEN (dev_name (x11) ) 
         IF (dev_id (iwin, x11) .le.0) stop 
         width = dev_sf (iwin, x11) * dev_width (iwin) 
      ELSE 
         CALL PGSLCT (dev_id (iwin, x11) ) 
         IF (dev_name (x11) (2:3) .ne.'GW') then 
            CALL PGCLOS 
            dev_id (iwin, x11) = PGOPEN (dev_name (x11) ) 
            CALL PGQVSZ (1, x1, x2, y1, y2) 
            width = x2 - x1 
            IF (dev_sf (iwin, x11) .ne.dev_sf_old) then 
               width = (x2 - x1) * dev_sf (iwin, x11) / dev_sf_old 
               dev_sf_old = dev_sf (iwin, x11) 
            ENDIF 
         ELSE 
            width = dev_sf (iwin, x11) * dev_width (iwin) 
         ENDIF 
      ENDIF 
!                                                                       
!------ Getting ready to plot                                           
!------ Colour_setup is done twice due to some strange                  
!------ difference between the GW and X11 drivers ????                  
!                                                                       
      CALL colour_setup 
      ratio = dev_height (iwin) / dev_width (iwin) 
      CALL PGPAP (width, ratio) 
      CALL PGASK (.false.) 
      CALL PGPAGE 
      CALL colour_setup 
!                                                                       
!------ here starts the plotting                                        
!                                                                       
      DO ii = 1, iaf (iwin) 
      CALL draw_frame (ii, lmenu) 
      ENDDO 
!                                                                       
!------ draw menu parts                                                 
!                                                                       
      IF (lmenu) call draw_menu 
!                                                                       
!------ draw identification                                             
!                                                                       
      IF (iden (iwin) ) then 
         CALL PGSCI (foncol (iwin, iframe, 3) ) 
         CALL PGIDEN 
      ENDIF 
!                                                                       
      END SUBROUTINE do_plot                        
!
!*********************************************************************  
!
SUBROUTINE do_hardcopy (befehl, zeile, lbef, lp) 
!+                                                                      
!     This is the plotting routine for hardcopies                       
!-                                                                      
USE build_name_mod
USE errlist_mod 
USE get_params_mod
USE prompt_mod 
USE kuplot_config 
USE kuplot_mod 
use kuplot_plot_low_mod
USE lib_do_operating_mod
USE lib_length
USE precision_mod
USE string_convert_mod
USE support_mod
use take_param_mod
!                                                                       
IMPLICIT none 
!                                                                       
INTEGER, PARAMETER :: maxw = 20
!                                                                       
CHARACTER(len=*), intent(inout) :: befehl
CHARACTER(len=*), intent(inout) :: zeile 
integer         , intent(inout) :: lbef
integer         , intent(inout) :: lp
!
CHARACTER(LEN=PREC_STRING) :: cpara (maxw), prnbef , line
CHARACTER(len=256) :: filname, uname , outname
CHARACTER(len=256), dimension(2) :: conv_opt
REAL(KIND=PREC_DP) :: werte (maxw) 
REAL(kind=PREC_SP) :: width, ratio 
INTEGER :: lpara (maxw) 
INTEGER :: ianz, ii, idev, i
integer :: udev        ! User device definition
LOGICAL :: tfr, lrena, lmenu 
!                                                                       
INTEGER :: PGOPEN 
!
!
integer, parameter :: NOPTIONAL = 2
integer, parameter :: O_RES     = 1
integer, parameter :: O_TRANS   = 2
character(LEN=   5), dimension(NOPTIONAL) :: oname   !Optional parameter names
character(LEN=PREC_STRING), dimension(NOPTIONAL) :: opara   !Optional parameter strings returned
integer            , dimension(NOPTIONAL) :: loname  !Lenght opt. para name
integer            , dimension(NOPTIONAL) :: lopara  !Lenght opt. para name returned
logical            , dimension(NOPTIONAL) :: lpresent!opt. para is present
real(kind=PREC_DP) , dimension(NOPTIONAL) :: owerte   ! Calculated values
integer, parameter                        :: ncalc = 1 ! Number of values to calculate 
!
data oname  / 'res', 'trans'  /
data loname /  3,     5       /
opara  =  (/ '300 ', 'none' /)   ! Always provide fresh default values
lopara =  (/  3,    4       /)
owerte =  (/  0.0,  0.0     /)
!
!
!                                                                       
lmenu = .false. 
lrena = .false. 
idev  = png
udev  = png
filname = ' '
uname = ' '
outname = ' '
!                                                                       
!------ see if there are data at all                                    
!                                                                       
tfr = .false. 
DO ii = 1, iaf (iwin) 
   tfr = tfr.or. (infra (iwin, ii, 1) .eq. - 1) 
ENDDO 
!                                                                       
IF (iz.le.1.and..not.tfr) then 
      ier_num = - 12 
      ier_typ = ER_APPL 
      RETURN 
ENDIF 
!                                                                       
CALL do_cap (befehl) 
CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
IF (ier_num.ne.0) return 
call get_optional(ianz, MAXW, cpara, lpara, NOPTIONAL,  ncalc, &
                  oname, loname, opara, lopara, lpresent, owerte)
!                                                                       
!------ Command save or print                                           
!                                                                       
IF (befehl (1:2) .eq.'SA') then 
   IF (ianz.ge.1) then 
      lrena = .false. 
      CALL do_cap (cpara (1) ) 
      IF (cpara (1) (1:2) .eq.'PS') then      ! Postscript
         filname = 'kuplot_temp.ps' 
           uname = 'kuplot.ps' 
         lrena   = .FALSE.
         IF (orient (iwin) ) then 
            idev = ps 
         ELSE 
            idev = vps 
         ENDIF 
         udev = idev 
      elseIF (cpara(1)(1:3) .eq.'EPS') then      ! Postscript
         filname = 'kuplot_temp.ps' 
           uname = 'kuplot.eps' 
         lrena   = .FALSE.
         IF (orient (iwin) ) then 
            idev = ps 
         ELSE 
            idev = vps 
         ENDIF 
         udev = eps 
      ELSEIF (cpara (1) (1:2) .eq.'PI') then  ! GIF, 
         filname = 'kuplot_temp.ps' 
           uname = 'kuplot.gif' 
         udev = pic 
         IF (orient (iwin) ) then 
            idev = ps 
         ELSE 
            idev = vps 
         ENDIF 
      ELSEIF (cpara (1) (1:2) .eq.'PN') then  ! PNG, make postscript and convert
         filname = 'kuplot_temp.ps' 
           uname = 'kuplot.png' 
         udev = png 
         if (orient (iwin) ) then
            idev = ps
         ELSE
            idev = vps
         ENDIF
      ELSEIF (cpara (1) (1:2) .eq.'LA') then  ! Latex
         filname = 'kuplot.tex' 
           uname = 'kuplot.tex' 
         idev = lat 
         udev = idev
      ELSEIF (cpara (1) (1:2) .eq.'PD') then  ! PDF, make postscript and convert
         filname = 'kuplot_temp.ps' 
           uname = 'kuplot.pdf' 
         udev = pdf 
         IF (orient (iwin) ) then 
            idev = ps 
         ELSE 
            idev = vps 
         ENDIF 
         udev = pdf
      ELSE                                   ! Unknown output format
         ier_num = - 11 
         ier_typ = ER_APPL 
         RETURN 
      ENDIF 
      IF (ianz.ge.2) then                    ! User specified file name
         CALL del_params (1, ianz, cpara, lpara, maxw) 
         CALL do_build_name (ianz, cpara, lpara, werte, maxw, 1) 
         uname = cpara (1) (1:MIN(256,LEN_TRIM(cpara(1))))
         lrena = .true. 
         if(udev==ps) then                ! User specified Postscript without filename
            filname = uname               ! No need for different file names
            lrena = .FALSE. 
         endif
      ELSE
         if(udev==ps) then                ! User specified Postscript without filename
            filname = 'kuplot.ps'         ! Default file name
            uname = filname 
            lrena = .FALSE. 
         endif
      ENDIF 
   ELSE 
      ier_num = - 6 
      ier_typ = ER_COMM 
      RETURN 
   ENDIF 
!                                                                       
   ii = len_str (uname) 
   WRITE (output_io, 1000) iwin, uname (1:ii) 
!                                                                       
!------ Command prin                                                    
!                                                                       
ELSEIF (befehl (1:2) .eq.'PR') then 
   filname = 'kuplot.plt' 
   idev = ps 
   IF (orient (iwin) ) then 
      idev = ps 
   ELSE 
      idev = vps 
   ENDIF 
   udev = idev
   lrena = .false. 
!                                                                       
   IF (ianz.ge.1) then 
      CALL do_cap (cpara (1) ) 
      IF (cpara (1) (1:2) .ne.'PS') then 
         ier_num = - 11 
         ier_typ = ER_APPL 
         RETURN 
      ENDIF 
   ENDIF 
!                                                                       
   IF (ianz.gt.1) then 
      lp = -lp
      CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
      CALL del_params (1, ianz, cpara, lpara, maxw) 
      CALL do_build_name (ianz, cpara, lpara, werte, maxw, 1) 
      prnbef = cpara (1) 
   ELSE 
      prnbef = dev_prn (idev) 
   ENDIF 
ENDIF 
!
IF (ier_num.ne.0) return 
!                                                                       
!------ Opening device                                                  
!                                                                       
dev_draw (iwin, 1) = 0.00 
dev_draw (iwin, 2) = 0.00 
!                                                                       
dev_id (iwin, idev) = PGOPEN (filname(1:len_str(filname))//dev_name (idev) ) 
IF (dev_id (iwin, idev) .le.0) then
   ier_num = -75
   ier_typ =ER_APPL
   return
endif
!                                                                       
width = dev_sf (iwin, idev) * dev_width (iwin) 
ratio = dev_height (iwin) / dev_width (iwin) 
!                                                                       
CALL PGPAP (width, ratio) 
CALL colour_setup 
!                                                                       
!------ here starts the plotting                                        
!                                                                       
DO ii = 1, iaf (iwin) 
   CALL draw_frame (ii, lmenu) 
ENDDO 
!                                                                       
!------ draw identification                                             
!                                                                       
IF (iden (iwin) ) then 
   CALL PGSCI (foncol (iwin, iframe, 3) ) 
   CALL PGIDEN 
ENDIF 
!                                                                       
CALL PGCLOS 
!                                                                       
!------ check if we need to rename ouputfile                            
if(udev/=lat)  then           ! No need if latex
   if(udev==ps) then
      if (lrena) call do_rename_file (filname, uname) 
   elseif(udev==eps) then
      line = 'ps2eps '//filname(1:len(filname))// ' -f ' //uname(1:len(uname))
      call system(line, ier_num)
   else                       ! All other picture types (pdf, png, gif)
!     ! Write the options for convert: density=Pixel per inch and 
!     !                                transparency
      write(conv_opt(1),'(a,i6,a)') ' -density ', nint(owerte(O_RES)),          &
           ' -units PixelsPerInch '
      conv_opt(2) = ' -background white -alpha off  -alpha remove '
      if(opara(O_TRANS)=='none'.or. opara(O_TRANS)=='solid') then
         conv_opt(2) = ' -background white -alpha off  -alpha remove '
      elseif(opara(O_TRANS)=='trans') then
         conv_opt(2) = ' '
      endif
      i = LEN_TRIM(uname)
      outname = uname
      if(udev==png) then            ! Create png file out of Postscript
         IF(uname(i-3:i)=='.png' .or. uname(i-3:i)=='.PNG') THEN
            uname(i-3:i) = '.ps '
            outname(i-3:i) = '.png'   ! 
         ENDIF
      elseif(udev==pic) then            ! Create GIF file out of Postscript
         IF(uname(i-3:i)=='.gif' .or. uname(i-3:i)=='.GIF') THEN
            uname(i-3:i) = '.ps '
            outname(i-3:i) = '.gif'   ! 
         ENDIF
      elseif(udev==pdf) then
         IF(uname(i-3:i)=='.pdf' .or. uname(i-3:i)=='.PDF') THEN
            uname(i-3:i) = '.ps '
            outname(i-3:i+1) = '.pdf'
         endif
      endif
!
      write(line,'(8a)') 'convert ', &
          conv_opt(1)(1:len_trim(conv_opt(1))), ' ',&
          filname(1:len_trim(filname)),       &
          ' -rotate 90 ',     &
          conv_opt(2)(1:len_trim(conv_opt(2))), ' ',&
          outname(1:len_trim(outname))
!         ' -density 300 -units PixelsPerInch ',                                 &
!         ' -background white -alpha off  -alpha remove ',     &
      CALL system(line, ier_num)
      line = 'rm -r kuplot_temp.ps'
      CALL system(line, ier_num)
   endif
endif
!                                                                       
!------ if the command was print we print                               
!                                                                       
IF (befehl (1:2) .eq.'PR') then 
   ii = len_str (prnbef) 
   prnbef = prnbef (1:ii) //' '//filname 
   ii = len_str (prnbef) 
   WRITE (output_io, 2000) iwin, prnbef (1:ii) 
   CALL do_operating (prnbef (1:ii), ii) 
   CALL do_del_file (filname) 
ENDIF 
!                                                                       
 1000 FORMAT     (' ------ > Saving window ',i3,' to file : ',a,' ..') 
 2000 FORMAT     (' ------ > Printing window ',i3,' using : ',a,' ..') 
!
END SUBROUTINE do_hardcopy                    
!
!*****7*****************************************************************
!
!*****7*****************************************************************
!
end module kuplot_plot_mod
