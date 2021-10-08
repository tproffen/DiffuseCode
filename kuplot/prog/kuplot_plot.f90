!********************************************************************   
!     These are the main plot routines for 2D and 3D data. 3D data      
!     can be plotted as contour plots or bitmaps or both. All           
!     settings are defined in 'kuplot.inc'.                             
!*********************************************************************  
SUBROUTINE do_plot (lmenu) 
!                                                                       
!     Main plotting routine                                             
!                                                                       
      USE errlist_mod 
      USE kuplot_config 
      USE kuplot_mod 
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
      SUBROUTINE colour_setup 
!+                                                                      
!     This sets the default KUPLOT colours ...                          
!-                                                                      
      USE kuplot_config 
      USE kuplot_mod 
!                                                                       
      IMPLICIT none 
!                                                                       
      INTEGER ic 
!                                                                       
      DO ic = 0, 15 
      CALL PGSCR (ic, colour (iwin, ic, 1), colour (iwin, ic, 2),       &
      colour (iwin, ic, 3) )                                            
      ENDDO 
!                                                                       
      END SUBROUTINE colour_setup                   
!*****7*****************************************************************
      SUBROUTINE open_frame (ii, lwhite) 
!+                                                                      
!     This routine sets background and border for frame ii              
!-                                                                      
      USE kuplot_config 
      USE kuplot_mod 
!                                                                       
      IMPLICIT none 
!                                                                       
      REAL x1, x2, y1, y2 
      INTEGER ic, ii 
      LOGICAL lwhite 
!                                                                       
!------ Plot border and background colour if required                   
!                                                                       
      x1 = frame (iwin, ii, 1) * (1.0 - dev_draw (iwin, 1) ) 
      x2 = frame (iwin, ii, 3) * (1.0 - dev_draw (iwin, 1) ) 
      y1 = frame (iwin, ii, 2) * (1.0 - dev_draw (iwin, 2) ) + dev_draw &
      (iwin, 2)                                                         
      y2 = frame (iwin, ii, 4) * (1.0 - dev_draw (iwin, 2) ) + dev_draw &
      (iwin, 2)                                                         
!                                                                       
      IF (y2.ge.1.0) y2 = y2 - 0.001 
      IF (x2.ge.1.0) x2 = x2 - 0.001 
!                                                                       
      CALL PGSVP (x1, x2, y1, y2) 
      CALL PGSWIN (0.0, 1.0, 0.0, 1.0) 
!                                                                       
      IF (frback (iwin, ii, 1) .ne.1.0.or.frback (iwin, ii, 2)          &
      .ne.1.0.or.frback (iwin, ii, 3) .ne.1.0.or.lwhite) then           
         ic = ii + 15 
         CALL PGSCR (ic, frback (iwin, ii, 1), frback (iwin, ii, 2),    &
         frback (iwin, ii, 3) )                                         
         IF (lwhite) then 
            CALL PGSCI (0) 
         ELSE 
            CALL PGSCI (ic) 
         ENDIF 
         CALL PGSFS (1) 
         CALL PGRECT (0.0, 1.0, 0.0, 1.0) 
      ENDIF 
!                                                                       
      IF (tot_frame (iwin) ) call border_frame (ii, 1, 0.0) 
!                                                                       
      END SUBROUTINE open_frame                     
!*****7*****************************************************************
      SUBROUTINE border_frame (ii, icol, delta) 
!+                                                                      
!     This routine draws border around frame ii                         
!-                                                                      
      USE kuplot_config 
      USE kuplot_mod 
!                                                                       
      IMPLICIT none 
!                                                                       
      REAL x1, x2, y1, y2, delta, deltay 
      INTEGER ii, icol 
!                                                                       
!------ Plot border and background colour if required                   
!                                                                       
      x1 = frame (iwin, ii, 1) * (1.0 - dev_draw (iwin, 1) ) 
      x2 = frame (iwin, ii, 3) * (1.0 - dev_draw (iwin, 1) ) 
      y1 = frame (iwin, ii, 2) * (1.0 - dev_draw (iwin, 2) ) + dev_draw &
      (iwin, 2)                                                         
      y2 = frame (iwin, ii, 4) * (1.0 - dev_draw (iwin, 2) ) + dev_draw &
      (iwin, 2)                                                         
!                                                                       
      IF (y2.ge.1.0) y2 = y2 - 0.001 
      IF (x2.ge.1.0) x2 = x2 - 0.001 
!                                                                       
      deltay = (y2 - y1) / (x2 - x1) 
      deltay = delta / deltay 
!                                                                       
      CALL PGSVP (x1, x2, y1, y2) 
      CALL PGSWIN (0.0, 1.0, 0.0, 1.0) 
!                                                                       
      CALL PGSCI (icol) 
      CALL PGSFS (2) 
      CALL PGRECT (0.0 + delta, 1.0 - delta, 0.0 + deltay, 1.0 - deltay) 
!                                                                       
      END SUBROUTINE border_frame                   
!*****7*****************************************************************
      SUBROUTINE open_viewport 
!+                                                                      
!     This routine sets the viewport and window for the current         
!     frame (-> iframe)                                                 
!-                                                                      
      USE koordinate_mod
      USE wink_mod
      USE kuplot_config 
      USE kuplot_mod 
!                                                                       
      IMPLICIT none 
!                                                                       
!     REAL wxmi, wxma, wymi, wyma 
      REAL xx, yy, x1, x2, y1, y2, vl, vr, vb, vt, yf, off 
!                                                                       
!------ Set viewport for current frame                                  
!                                                                       
      vl = min (1.0, (frame (iwin, iframe, 1) + ibuf (iwin, iframe, 1) )&
      )                                                                 
      vr = max (0.0, (frame (iwin, iframe, 3) - ibuf (iwin, iframe, 2) )&
      )                                                                 
      vb = min (1.0, (frame (iwin, iframe, 2) + ibuf (iwin, iframe, 3) )&
      )                                                                 
      vt = max (0.0, (frame (iwin, iframe, 4) - ibuf (iwin, iframe, 4) )&
      )                                                                 
!                                                                       
      vl = vl * (1.0 - dev_draw (iwin, 1) ) 
      vr = vr * (1.0 - dev_draw (iwin, 1) ) 
      vb = vb * (1.0 - dev_draw (iwin, 2) ) + dev_draw (iwin, 2) 
      vt = vt * (1.0 - dev_draw (iwin, 2) ) + dev_draw (iwin, 2) 
!                                                                       
      CALL PGSVP (vl, vr, vb, vt) 
!                                                                       
!------ Force user defined aspect ratio if required                     
!                                                                       
      IF (.not.lyskal (iwin, iframe) .and.shear (iwin, iframe) .ne.90.0)&
      then                                                              
         lyskal (iwin, iframe) = .true. 
         yskal_u (iwin, iframe) = 1.0 
      ENDIF 
!                                                                       
      yf = sin (REAL(rad) * shear (iwin, iframe) ) 
      IF (lyskal (iwin, iframe) ) then 
         yskal (iwin, iframe) = yskal_u (iwin, iframe) * yf 
         xx = pex (iwin, iframe, 2) - pex (iwin, iframe, 1) 
         yy = pey (iwin, iframe, 2) - pey (iwin, iframe, 1) 
         x1 = pex (iwin, iframe, 2) 
         y1 = pey (iwin, iframe, 2) 
         CALL koor_shear (1, x1, y1) 
         off = x1 - pex (iwin, iframe, 2) 
         CALL PGWNAD (0.0, 1.0, 0.0, (yy / (xx + abs (off) ) ) * yskal (&
         iwin, iframe) )                                                
      ELSE 
         CALL PGQVP (2, x1, x2, y1, y2) 
         yskal (iwin, iframe) = (y2 - y1) / (x2 - x1) 
         yskal_u (iwin, iframe) = yskal (iwin, iframe) / yf 
      ENDIF 
!                                                                       
!------ Set coordinates for current frame                               
!                                                                       
      x1 = pex (iwin, iframe, 2) 
      y1 = pey (iwin, iframe, 2) 
      CALL koor_shear (1, x1, y1) 
      off = x1 - pex (iwin, iframe, 2) 
!                                                                       
      IF (shear (iwin, iframe) .le.90.0) then 
         CALL PGSWIN (pex (iwin, iframe, 1), pex (iwin, iframe, 2)      &
         + off, pey (iwin, iframe, 1), pey (iwin, iframe, 2) )          
      ELSE 
         CALL PGSWIN (pex (iwin, iframe, 1) + off, pex (iwin, iframe, 2)&
         , pey (iwin, iframe, 1), pey (iwin, iframe, 2) )               
      ENDIF 
!                                                                       
      END SUBROUTINE open_viewport                  
!*********************************************************************  
      SUBROUTINE draw_frame (ii, lmenu) 
!+                                                                      
!     Here we actually do the plotting ...                              
!-                                                                      
      USE errlist_mod 
      USE kuplot_config 
      USE kuplot_mod 
USE kuplot_3dm_mod
USE kuplot_3dm_draw
!                                                                       
      IMPLICIT none 
!                                                                       
      INTEGER ik, ii, iframe_old 
      LOGICAL k_in_f, lmenu 
!                                                                       
      CALL PGBBUF 
!                                                                       
      iframe_old = iframe 
      iframe = ii 
!                                                                       
      CALL open_frame (iframe, lmenu) 
main: IF (infra (iwin, iframe, 1) .eq. - 1) then 
   CALL draw_text_frame 
ELSE main
         CALL skalieren 
         CALL open_viewport 
   DO ik = 1, iz - 1 
      IF(k3dm_ik==ik) THEN
         CALL kuplot_draw_3d_static(ik)
         EXIT main
      ELSE
      IF ( (hlineart (iwin, iframe, ik) .eq.2.or.hlineart (iwin,     &
         iframe, ik) .eq.3) .and.k_in_f (ik) .and.lni (ik) ) call       &
         draw_bitmap (ik)                                               
      ENDIF
         ENDDO 
         CALL draw_werte 
         CALL draw_rahmen 
      ENDIF main 
!                                                                       
      iframe = iframe_old 
!                                                                       
      CALL PGEBUF 
!                                                                       
      END SUBROUTINE draw_frame                     
!****7******************************************************************
SUBROUTINE draw_text_frame 
!+                                                                      
!       plot text frame                                                 
!-                                                                      
USE errlist_mod 
USE wink_mod
USE kuplot_config 
USE kuplot_mod 
USE kuplot_fit6
USE lib_length
use precision_mod
USE string_convert_mod
USE support_mod
!                                                                       
IMPLICIT none 
!                                                                       
CHARACTER(len=PREC_STRING) :: zeile, filname 
REAL :: xh, yh, xt, yt 
INTEGER :: i 
!                                                                       
!                                                                       
      filname = ftext (iwin, iframe) 
      CALL do_cap (filname) 
!                                                                       
      IF (filname.eq.'PARA'.or.filname.eq.'KUPL.PAR') then 
         CALL write_para 
         filname = 'kupl.par' 
      ELSEIF (filname.eq.'FIT'.or.filname.eq.'KUPL.FIT') then 
         CALL write_fit 
         filname = 'kupl.fit' 
      ELSE 
         filname = ftext (iwin, iframe) 
      ENDIF 
!                                                                       
      CALL oeffne (21, filname, 'unknown') 
      IF (ier_num.ne.0) return 
!                                                                       
      CALL PGSCI (foncol (iwin, iframe, 5) ) 
      CALL PGSLW (nint (linewid (iwin, iframe, 0) / 0.13) ) 
      CALL PGSCF (fon_id (iwin, iframe, 5) ) 
      CALL PGSCH (fonscal (iwin, iframe) * fonsize (iwin, iframe, 5)    &
      / 120.0)                                                          
      CALL PGQCS (4, xh, yh) 
!                                                                       
      i = 0 
   10 CONTINUE 
      READ (21, '(a)', end = 20) zeile 
      IF (frjust (iwin, iframe) .eq.if_left) then 
         xt = 0.0 
      ELSE 
         xt = 0.5 
      ENDIF 
      yt = 1.0 - (i + 1) * 1.05 * yh 
      IF (len_str (zeile) .gt.0) then 
         CALL PGPTXT (xt, yt, 0.0, xt, zeile (1:len_str (zeile) ) ) 
      ENDIF 
      i = i + 1 
      GOTO 10 
   20 CONTINUE 
      CLOSE (21) 
!                                                                       
      END SUBROUTINE draw_text_frame                
!*****7*****************************************************************
      SUBROUTINE draw_rahmen 
!+                                                                      
!     This routine draws the labeled frame around the plot.             
!-                                                                      
      USE koordinate_mod
      USE kuplot_config 
      USE kuplot_mod 
!                                                                       
      IMPLICIT none 
!                                                                       
      CHARACTER(2) cal, car, cax 
      CHARACTER(1) clx, cly 
      REAL ax1, ay1, off, off0, adis 
      REAL rj, xt, yt, ltx, lty 
      REAL xpl (2), ypl (2) 
      REAL xh (5), yh (5) 
      INTEGER  :: ninterv
      INTEGER  :: i
!                                                                       
      ax1 = pex (iwin, iframe, 2) 
      ay1 = pey (iwin, iframe, 2) 
      CALL koor_shear (1, ax1, ay1) 
      off = ax1 - pex (iwin, iframe, 2) 
!                                                                       
      CALL PGSLS (1) 
      CALL PGSLW (nint (linewid (iwin, iframe, 0) / 0.13) ) 
      CALL PGSCI (6) 
!                                                                       
!------ IBOX = 0 : nothing                                              
!                                                                       
      IF (ibox (iwin, iframe) .eq.0) return 
!                                                                       
!------ IBOX = 1 : just a box                                           
!                                                                       
      IF (abs (ibox (iwin, iframe) ) .eq.1) then 
         xh (1) = pex (iwin, iframe, 1) 
         yh (1) = pey (iwin, iframe, 1) 
         xh (2) = pex (iwin, iframe, 2) 
         yh (2) = pey (iwin, iframe, 1) 
         xh (3) = pex (iwin, iframe, 2) 
         yh (3) = pey (iwin, iframe, 2) 
         xh (4) = pex (iwin, iframe, 1) 
         yh (4) = pey (iwin, iframe, 2) 
         xh (5) = pex (iwin, iframe, 1) 
         yh (5) = pey (iwin, iframe, 1) 
         CALL koor_shear (5, xh, yh) 
         CALL PGSCLP (0) 
         CALL PGLINE (5, xh, yh) 
         CALL PGSCLP (1) 
!                                                                       
!------ IBOX >= 2 : tickmarks, numbers, ..                              
!                                                                       
      ELSE 
!                                                                       
!------ - Grid lines                                                    
!                                                                       
         IF (igrid (iwin, iframe) ) then 
            CALL PGSCI (ilinecol (iwin, iframe, 0) ) 
            CALL PGSLS (4) 
!                                                                       
            IF (pex (iwin, iframe, 1) .lt.0.0) then 
               ninterv = IABS(NINT((pex (iwin, iframe, 1) - 0.0 )/pt (iwin, iframe, 1)))
               DO i = 0, ninterv !xt = 0.0, pex (iwin, iframe, 1), - pt (iwin, iframe, 1)
                  xt = 0.0 - i*pt (iwin, iframe, 1)
               IF (xt.lt.pex (iwin, iframe, 2) ) then 
                  xt = 0.0 - i* pt (iwin, iframe, 1)
                  xpl (1) = xt 
                  ypl (1) = pey (iwin, iframe, 1) 
                  xpl (2) = xt 
                  ypl (2) = pey (iwin, iframe, 2) 
                  CALL koor_shear (2, xpl, ypl) 
                  CALL PGLINE (2, xpl, ypl) 
               ENDIF 
               ENDDO 
            ENDIF 
            IF (pex (iwin, iframe, 2) .gt.0.0) then 
               ninterv = IABS(NINT((pex (iwin, iframe, 2) - 0.0 )/pt (iwin, iframe, 1)))
               DO i = 0, ninterv !xt = 0.0, pex (iwin, iframe, 2), pt (iwin, iframe, 1) 
                  xt = 0.0 + i*pt (iwin, iframe, 1)
               IF (xt.gt.pex (iwin, iframe, 1) ) then 
                  xt = 0.0 + i* pt (iwin, iframe, 1)
                  xpl (1) = xt 
                  ypl (1) = pey (iwin, iframe, 1) 
                  xpl (2) = xt 
                  ypl (2) = pey (iwin, iframe, 2) 
                  CALL koor_shear (2, xpl, ypl) 
                  CALL PGLINE (2, xpl, ypl) 
               ENDIF 
               ENDDO 
            ENDIF 
!                                                                       
            IF (pey (iwin, iframe, 1) .lt.0.0) then 
               ninterv = IABS(NINT((pey (iwin, iframe, 1) - 0.0 )/pt (iwin, iframe, 2)))
               DO i = 0, ninterv !yt = 0.0, pey (iwin, iframe, 1), - pt (iwin, iframe, 2)
                  yt = 0.0 - i* pt (iwin, iframe, 2)
               xpl (1) = pex (iwin, iframe, 1) 
               ypl (1) = yt 
               xpl (2) = pex (iwin, iframe, 2) 
               ypl (2) = yt 
               CALL koor_shear (2, xpl, ypl) 
               CALL PGLINE (2, xpl, ypl) 
               ENDDO 
            ENDIF 
            IF (pey (iwin, iframe, 2) .gt.0.0) then 
               ninterv = IABS(NINT((pey (iwin, iframe, 2) - 0.0 )/pt (iwin, iframe, 2)))
               DO i  = 0,  ninterv !yt = 0.0, pey (iwin, iframe, 2), pt (iwin, iframe, 2) 
                  yt = 0.0 + i* pt (iwin, iframe, 2)
               xpl (1) = pex (iwin, iframe, 1) 
               ypl (1) = yt 
               xpl (2) = pex (iwin, iframe, 2) 
               ypl (2) = yt 
               CALL koor_shear (2, xpl, ypl) 
               CALL PGLINE (2, xpl, ypl) 
               ENDDO 
            ENDIF 
            CALL PGSCI (6) 
         ENDIF 
!                                                                       
!------ - Axes                                                          
!                                                                       
         clx = ' ' 
         cly = ' ' 
         IF (lachse (iwin, iframe, 1) ) clx = 'L' 
         IF (lachse (iwin, iframe, 2) ) cly = 'L' 
!                                                                       
         xpl (1) = pex (iwin, iframe, 1) 
         xpl (2) = pex (iwin, iframe, 2) 
         ypl (1) = pey (iwin, iframe, 1) 
         ypl (2) = pey (iwin, iframe, 2) 
!                                                                       
         ltx = pt (iwin, iframe, 1) 
         lty = pt (iwin, iframe, 2) 
!                                                                       
         IF (ibox (iwin, iframe) .gt.1) then 
            cal = cly//'N' 
            car = cly//' ' 
            adis = - lab_d (iwin, iframe, 2) 
         ELSE 
            cal = cly//' ' 
            car = cly//'N' 
            adis = lab_d (iwin, iframe, 2) 
         ENDIF 
!                                                                       
         cax = clx//'N' 
      IF (achse (iwin, iframe, 1)  (1:3) .eq.'OFF') cax = '  ' 
      IF (achse (iwin, iframe, 2)  (1:3) .eq.'OFF') cal = '  ' 
      IF (achse (iwin, iframe, 2)  (1:3) .eq.'OFF') car = '  ' 
!                                                                       
         CALL PGSLS (1) 
         CALL PGSLW (nint (linewid (iwin, iframe, 0) / 0.13) ) 
         CALL PGSCI (foncol (iwin, iframe, 4) ) 
         CALL PGSCF (fon_id (iwin, iframe, 4) ) 
         CALL PGSCH (fonscal (iwin, iframe) * fonsize (iwin, iframe, 4) &
         / 120.0)                                                       
         CALL PGQCS (4, xh (4), yh (4) ) 
         CALL PGAXIS (cax, xpl (1), ypl (1), xpl (2), ypl (1), xpl (1), &
         xpl (2), ltx, tick_nsub (iwin, iframe, 1), tick_ma_h (iwin,    &
         iframe, 1), 0.0, tick_mi_h (iwin, iframe, 1), lab_d (iwin,     &
         iframe, 1), lab_angle (iwin, iframe, 1) )                      
         CALL PGAXIS (clx, xpl (1) + off, ypl (2), xpl (2) + off, ypl ( &
         2), xpl (1), xpl (2), ltx, tick_nsub (iwin, iframe, 1),        &
         0.0, tick_ma_h (iwin, iframe, 1), tick_mi_h (iwin, iframe, 1), &
         lab_d (iwin, iframe, 1), lab_angle (iwin, iframe, 1) )         
!                                                                       
         CALL PGAXIS (cal, xpl (1), ypl (1), xpl (1) + off, ypl (2),    &
         ypl (1), ypl (2), lty, tick_nsub (iwin, iframe, 2), 0.0,       &
         tick_ma_h (iwin, iframe, 2), tick_mi_h (iwin, iframe, 2),      &
         adis, lab_angle (iwin, iframe, 2) )                            
         CALL PGAXIS (car, xpl (2), ypl (1), xpl (2) + off, ypl (2),    &
         ypl (1), ypl (2), lty, tick_nsub (iwin, iframe, 2), tick_ma_h (&
         iwin, iframe, 2), 0.0, tick_mi_h (iwin, iframe, 2), adis,      &
         lab_angle (iwin, iframe, 2) )                                  
!                                                                       
!------ - Axis at x=0.0                                                 
!                                                                       
         IF (abs (ibox (iwin, iframe) ) .ge.3.and.pex (iwin, iframe, 1) &
         .lt.0.0.and.pex (iwin, iframe, 2) .gt.0.0) then                
            ax1 = 0.0 
            ay1 = pey (iwin, iframe, 2) 
            CALL koor_shear (1, ax1, ay1) 
            off0 = ax1 
            CALL PGAXIS (cly, 0.0, ypl (1), off0, ypl (2), ypl (1),     &
            ypl (2), lty, tick_nsub (iwin, iframe, 1), tick_ma_h (iwin, &
            iframe, 2), tick_ma_h (iwin, iframe, 2), tick_mi_h (iwin,   &
            iframe, 2), 0.0, 0.0)                                       
         ENDIF 
!                                                                       
!------ - Axis at y=0.0                                                 
!                                                                       
         IF (abs (ibox (iwin, iframe) ) .ge.3.and.ey (iwin, iframe, 1)  &
         .lt.0.0.and.ey (iwin, iframe, 2) .gt.0.0) then                 
            ax1 = pex (iwin, iframe, 2) 
            ay1 = 0.0 
            CALL koor_shear (1, ax1, ay1) 
            off0 = ax1 - pex (iwin, iframe, 2) 
      CALL PGAXIS (clx, xpl (1)  + off0, 0.0, xpl (2)  + off0, 0.0, xpl &
     &(1),  xpl (2),  ltx, tick_nsub (iwin, iframe, 2),  tick_ma_h (iwin&
     &, iframe, 1),  tick_ma_h (iwin, iframe, 1),  tick_mi_h (iwin, ifra&
     &me, 1),  0.0, 0.0)                                                
         ENDIF 
!                                                                       
!------ - Axes labels                                                   
!                                                                       
         CALL PGSCI (foncol (iwin, iframe, 3) ) 
         CALL PGSCF (fon_id (iwin, iframe, 3) ) 
         CALL PGSCH (fonscal (iwin, iframe) * fonsize (iwin, iframe, 3) &
         / 120.0)                                                       
         CALL PGQCS (4, xh (3), yh (3) ) 
!                                                                       
         IF (achse (iwin, iframe, 1) (1:3) .ne.'OFF') then 
            xt = pex (iwin, iframe, 1) + 0.5 * (pex (iwin, iframe, 2)   &
            - pex (iwin, iframe, 1) )                                   
            yt = pey (iwin, iframe, 1) - ax_d (iwin, iframe, 1) * yh (4)&
            - 1.2 * yh (3)                                              
            CALL PGPTXT (xt, yt, 0.0, 0.5, achse (iwin, iframe, 1) ) 
         ENDIF 
!                                                                       
         IF (achse (iwin, iframe, 2) (1:3) .ne.'OFF') then 
            IF (ibox (iwin, iframe) .gt.1) then 
               yt = pey (iwin, iframe, 1) + 0.5 * (pey (iwin, iframe, 2)&
               - pey (iwin, iframe, 1) )                                
               xt = pex (iwin, iframe, 1) - ax_d (iwin, iframe, 2)      &
               * xh (4) - 0.2 * xh (3)                                  
               CALL koor_shear (1, xt, yt) 
               CALL PGPTXT (xt, yt, shear (iwin, iframe), 0.5, achse (  &
               iwin, iframe, 2) )                                       
            ELSE 
               yt = pey (iwin, iframe, 1) + 0.5 * (pey (iwin, iframe, 2)&
               - pey (iwin, iframe, 1) )                                
               xt = pex (iwin, iframe, 2) + (ax_d (iwin, iframe, 2)     &
               + 1) * xh (4) + 0.2 * xh (3)                             
               CALL koor_shear (1, xt, yt) 
               CALL PGPTXT (xt, yt, shear (iwin, iframe), 0.5, achse (  &
               iwin, iframe, 2) )                                       
            ENDIF 
         ENDIF 
                                                                        
      ENDIF 
!                                                                       
!------ Title                                                           
!                                                                       
      ax1 = pex (iwin, iframe, 1) 
      ay1 = pey (iwin, iframe, 2) 
      CALL koor_shear (1, ax1, ay1) 
      off = ax1 - pex (iwin, iframe, 1) 
!                                                                       
      IF (frjust (iwin, iframe) .eq.if_left) then 
         xt = pex (iwin, iframe, 1) + off 
         yt = pey (iwin, iframe, 2) 
         rj = 0.0 
      ELSE 
         xt = pex (iwin, iframe, 1) + 0.5 * (pex (iwin, iframe, 2)      &
         - pex (iwin, iframe, 1) ) + off                                
         yt = pey (iwin, iframe, 2) 
         rj = 0.5 
      ENDIF 
!                                                                       
      CALL PGSCI (foncol (iwin, iframe, 2) ) 
      CALL PGSCF (fon_id (iwin, iframe, 2) ) 
      CALL PGSCH (fonscal (iwin, iframe) * fonsize (iwin, iframe, 2)    &
      / 120.0)                                                          
      CALL PGQCS (4, xh (2), yh (2) ) 
      yt = yt + 0.5 * yh (2) 
      CALL PGPTXT (xt, yt, 0.0, rj, titel (iwin, iframe, 2) ) 
!                                                                       
      CALL PGSCI (foncol (iwin, iframe, 1) ) 
      CALL PGSCF (fon_id (iwin, iframe, 1) ) 
      CALL PGSCH (fonscal (iwin, iframe) * fonsize (iwin, iframe, 1)    &
      / 120.0)                                                          
      CALL PGQCS (4, xh (1), yh (1) ) 
      yt = yt + 1.1 * yh (1) 
      CALL PGPTXT (xt, yt, 0.0, rj, titel (iwin, iframe, 1) ) 
!                                                                       
      END SUBROUTINE draw_rahmen                    
!*****7*****************************************************************
      SUBROUTINE draw_werte 
!+                                                                      
!     This routine plots the datasets ...                               
!-                                                                      
      USE koordinate_mod
      USE kuplot_config 
      USE kuplot_mod 
USE lib_length
!                                                                       
      IMPLICIT none 
!                                                                       
      INTEGER maxmax 
      PARAMETER (maxmax = 50) 
!                                                                       
      CHARACTER(80) text 
      REAL xpl (maxarray), ypl (maxarray) 
      REAL wmax (maxmax) 
      REAL eex (2), eey (2) 
      REAL xleg, yleg, xfehl (2), yfehl (2) 
      REAL rj, xh, yh, xp, yp, xt, yt, sx, sy, xma, yma, delx 
      INTEGER ixm (maxmax), iym (maxmax) 
      INTEGER ikurv, ipkt, npkt 
      INTEGER il, ima, ikr, iityp, iicol 
      LOGICAL k_in_f 
      INTEGER :: jkurv
      REAL    :: max_value, min_value, range, factor, sma
!                                                                       
      LOGICAL inrect 
!                                                                       
!------ Loop over all data sets                                         
!                                                                       
      ikr = 0 
      DO ikurv = 1, iz - 1 
      IF (.not.k_in_f (ikurv) ) goto 999 
!                                                                       
!-------- draw lines                                                    
!                                                                       
      IF (lni (ikurv) ) then 
         IF (hlineart (iwin, iframe, ikurv) .eq.1.or.hlineart (iwin,    &
         iframe, ikurv) .eq.3) then                                     
            CALL draw_contour (ikurv) 
         ENDIF 
      ELSE 
         IF (fform (ikurv) .ne.'CR') then 
            DO ipkt = 2, lenc(ikurv) 
            xpl (ipkt - 1) = x (offxy (ikurv - 1) + (ipkt - 1) ) 
            ypl (ipkt - 1) = y (offxy (ikurv - 1) + (ipkt - 1) ) 
            xpl (ipkt) = x (offxy (ikurv - 1) + ipkt) 
            ypl (ipkt) = y (offxy (ikurv - 1) + ipkt) 
            ENDDO 
            npkt = lenc(ikurv) 
            IF (npkt.gt.1) then 
               CALL draw_points (xpl, ypl, npkt, ikurv) 
            ENDIF 
         ENDIF 
      ENDIF 
!                                                                       
!-------- Markers                                                       
!                                                                       
      IF (fform (ikurv) .eq.'CR') then 
         CALL draw_bonds (ikurv) 
         DO ipkt = 1, lenc(ikurv) 
         delx = dy (offxy (ikurv - 1) + ipkt) * sizemark (iwin, iframe, &
         ikurv)                                                         
         xma = x (offxy (ikurv - 1) + ipkt) 
         yma = y (offxy (ikurv - 1) + ipkt) 
         iicol = int (dx (offxy (ikurv - 1) + ipkt) / 1000) 
         iityp = int (dx (offxy (ikurv - 1) + ipkt) - 1000.0 * iicol) 
         CALL draw_marker (xma, yma, iityp, iicol, delx) 
         ENDDO 
      ELSE 
         IF (imarktyp (iwin, iframe, ikurv) .ne.0.and..not.lni (ikurv) )&
         then                                                           
            IF(rel_mark(iwin, iframe, ikurv) == 0 ) THEN
            CALL draw_bonds (ikurv) 
            DO ipkt = 1, lenc(ikurv) 
            xma = x (offxy (ikurv - 1) + ipkt) 
            yma = y (offxy (ikurv - 1) + ipkt) 
            CALL draw_marker (xma, yma, imarktyp (iwin, iframe, ikurv), &
            imarkcol (iwin, iframe, ikurv), sizemark (iwin, iframe,     &
            ikurv) )                                                    
            ENDDO 
            ELSE    ! proportional marker sizes
               IF(rel_mark(iwin, iframe, ikurv) == -1 ) THEN  ! x-scale
                  max_value = ABS(x (offxy (ikurv - 1) + 1))
                  min_value = ABS(x (offxy (ikurv - 1) + 1))
                  range = 0.0
                  DO ipkt = 1, lenc(ikurv)
                     max_value = MAX(max_value, ABS(x (offxy (ikurv - 1) + ipkt)))
                     min_value = MIN(min_value, ABS(x (offxy (ikurv - 1) + ipkt)))
                  ENDDO
                  IF(max_value==min_value) THEN
                      range = 1.0
                  ELSE
                      range = (max_value-min_value)
                  ENDIF
                  DO ipkt = 1, lenc(ikurv)
                     xma = x (offxy (ikurv - 1) + ipkt) 
                     yma = y (offxy (ikurv - 1) + ipkt) 
                     factor = 0.1+0.9*(ABS(x (offxy (ikurv - 1) + ipkt))-min_value)/range
                     IF(factor==0.0) factor = 1.0
                     sma = MAX(0.01,sizemark (iwin, iframe,ikurv)*factor)
                     CALL draw_marker (xma, yma, imarktyp (iwin, iframe, ikurv), &
                     imarkcol (iwin, iframe, ikurv), sma )
                  ENDDO
               ELSEIF(rel_mark(iwin, iframe, ikurv) == -2 ) THEN  ! y-scale
                  max_value = ABS(y (offxy (ikurv - 1) + 1))
                  min_value = ABS(y (offxy (ikurv - 1) + 1))
                  range = 0.0
                  DO ipkt = 1, lenc(ikurv)
                     max_value = MAX(max_value, ABS(y (offxy (ikurv - 1) + ipkt)))
                     min_value = MIN(min_value, ABS(y (offxy (ikurv - 1) + ipkt)))
                  ENDDO
                  IF(max_value==min_value) THEN
                      range = 1.0
                  ELSE
                      range = (max_value-min_value)
                  ENDIF
                  DO ipkt = 1, lenc(ikurv)
                     xma = x (offxy (ikurv - 1) + ipkt) 
                     yma = y (offxy (ikurv - 1) + ipkt) 
                     factor = 0.1+0.9*(ABS(y (offxy (ikurv - 1) + ipkt))-min_value)/range
                     IF(factor==0.0) factor = 1.0
                     sma = MAX(0.01,sizemark (iwin, iframe,ikurv)*factor)
                     CALL draw_marker (xma, yma, imarktyp (iwin, iframe, ikurv), &
                     imarkcol (iwin, iframe, ikurv), sma )
                  ENDDO
               ELSEIF(rel_mark(iwin, iframe, ikurv) >   0 ) THEN  ! scale by data set
                  jkurv = rel_mark(iwin, iframe, ikurv)
                  max_value = ABS(y (offxy (jkurv - 1) + 1))
                  min_value = ABS(y (offxy (jkurv - 1) + 1))
                  range = 0.0
                  DO ipkt = 1, lenc(jkurv)
                     max_value = MAX(max_value, ABS(y (offxy (jkurv - 1) + ipkt)))
                     min_value = MIN(min_value, ABS(y (offxy (jkurv - 1) + ipkt)))
                  ENDDO
                  IF(max_value==min_value) THEN
                      range = 1.0
                  ELSE
                      range = (max_value-min_value)
                  ENDIF
                  DO ipkt = 1, lenc(ikurv)
                     xma = x (offxy (ikurv - 1) + ipkt) 
                     yma = y (offxy (ikurv - 1) + ipkt) 
                     factor = 0.1+0.9*(ABS(y (offxy (jkurv - 1) + ipkt))-min_value)/range
                     IF(factor==0.0) factor = 1.0
                     sma = MAX(0.01,sizemark (iwin, iframe,ikurv)*factor)
                     CALL draw_marker (xma, yma, imarktyp (iwin, iframe, ikurv), &
                     imarkcol (iwin, iframe, ikurv), sma )
                  ENDDO
               ENDIF 
            ENDIF 
         ENDIF 
      ENDIF 
!                                                                       
!-------- Marker for maxima                                             
!                                                                       
      IF (imarkmax (iwin, iframe, ikurv) .ne.0) then 
         IF (lni (ikurv) ) then 
            CALL do_fmax_z (ikurv, wmax, ixm, iym, maxmax, ima) 
         ELSE 
            CALL do_fmax_xy (ikurv, wmax, ixm, maxmax, ima) 
         ENDIF 
         DO ipkt = 1, ima 
         IF (.not.lni (ikurv) ) iym (ipkt) = ixm (ipkt) 
         xma = x (ixm (ipkt) ) 
         yma = y (iym (ipkt) ) 
         CALL draw_marker (xma, yma, imarkmax (iwin, iframe, ikurv),    &
         imarkcol (iwin, iframe, ikurv), sizemark (iwin, iframe, ikurv) &
         )                                                              
         ENDDO 
      ENDIF 
!                                                                       
!-------- Error bars                                                    
!                                                                       
      IF (ierr (iwin, iframe, ikurv) .ne.0.and..not.lni (ikurv) ) then 
         CALL PGSLS (1) 
         CALL PGSLW (nint (linewid (iwin, iframe, ikurv) / 0.13) ) 
         CALL PGSCI (ierrcol (iwin, iframe, ikurv) ) 
         DO ipkt = 1, lenc(ikurv) 
!                                                                       
         IF (ierr (iwin, iframe, ikurv) .eq.1.or.ierr (iwin, iframe,    &
         ikurv) .eq.3) then                                             
            xfehl (1) = x (offxy (ikurv - 1) + ipkt) - dx (offxy (ikurv &
            - 1) + ipkt)                                                
            yfehl (1) = y (offxy (ikurv - 1) + ipkt) 
            xfehl (2) = x (offxy (ikurv - 1) + ipkt) + dx (offxy (ikurv &
            - 1) + ipkt)                                                
            yfehl (2) = y (offxy (ikurv - 1) + ipkt) 
            CALL koor_shear (2, xfehl, yfehl) 
            CALL koor_log (2, xfehl, yfehl) 
            CALL PGERRX (1, xfehl (1), xfehl (2), yfehl (1), 1.0) 
         ENDIF 
!                                                                       
         IF (ierr (iwin, iframe, ikurv) .eq.2.or.ierr (iwin, iframe,    &
         ikurv) .eq.3) then                                             
            xfehl (1) = x (offxy (ikurv - 1) + ipkt) 
            yfehl (1) = y (offxy (ikurv - 1) + ipkt) - dy (offxy (ikurv &
            - 1) + ipkt)                                                
            xfehl (2) = x (offxy (ikurv - 1) + ipkt) 
            yfehl (2) = y (offxy (ikurv - 1) + ipkt) + dy (offxy (ikurv &
            - 1) + ipkt)                                                
            CALL koor_shear (2, xfehl, yfehl) 
            CALL koor_log (2, xfehl, yfehl) 
            CALL PGERRY (1, xfehl (1), yfehl (1), yfehl (2), 1.0) 
         ENDIF 
         ENDDO 
      ENDIF 
!                                                                       
!-------- Filenames                                                     
!                                                                       
      IF (ifname (iwin, iframe) ) then 
         text = fname (ikurv) (1:MIN(80,LEN_TRIM(fname(ikurv))))
         il = len_str (text) 
         IF (il.gt.0) then 
            ikr = ikr + 1 
            CALL PGSCI (ilinecol (iwin, iframe, ikurv) ) 
            CALL PGSCF (fon_id (iwin, iframe, 6) ) 
            CALL PGSLW (nint (linewid (iwin, iframe, 0) / 0.13) ) 
            CALL PGSCH (fonscal (iwin, iframe) * fonsize (iwin, iframe, &
            6) / 120.0)                                                 
            CALL PGQCS (4, xh, yh) 
            IF (info_orig (iwin, iframe, 1) .ne. - 9999.) then 
               xt = info_orig (iwin, iframe, 1) 
               yt = info_orig (iwin, iframe, 2) 
            ELSE 
               xt = pex (iwin, iframe, 2) 
               yt = pey (iwin, iframe, 2) 
            ENDIF 
            xt = xt - xh 
            yt = yt - ikr * yh * 1.2 - 0.3 * yh 
            CALL koor_shear (1, xt, yt) 
            CALL PGPTXT (xt, yt, 0.0, 1.0, text (1:il) ) 
         ENDIF 
      ENDIF 
!                                                                       
!-------- Possible legend                                               
!                                                                       
      IF (ilegend (iwin, iframe, ikurv) .eq.1) then 
         ikr = ikr + 1 
         CALL PGSCF (fon_id (iwin, iframe, 6) ) 
         CALL PGSCH (fonscal (iwin, iframe) * fonsize (iwin, iframe, 6) &
         / 120.0)                                                       
         CALL PGQCS (4, xh, yh) 
         text = clegend (iwin, iframe, ikurv) 
         IF (info_orig (iwin, iframe, 1) .ne. - 9999.) then 
            xleg = info_orig (iwin, iframe, 1) 
            yleg = info_orig (iwin, iframe, 2) 
         ELSE 
            xleg = pex (iwin, iframe, 2) 
            yleg = ey (iwin, iframe, 2) 
         ENDIF 
         xleg = xleg - xh 
         yleg = yleg - ikr * yh * 1.2 
         CALL koor_shear (1, xleg, yleg) 
!                                                                       
         IF (ilinetyp (iwin, iframe, ikurv) .ne.0) then 
            CALL PGSCI (ilinecol (iwin, iframe, ikurv) ) 
            CALL PGSLS (ilinetyp (iwin, iframe, ikurv) ) 
            CALL PGSLW (nint (linewid (iwin, iframe, ikurv) / 0.13) ) 
            CALL PGERRX (1, xleg, xleg - 4.0 * xh, yleg, 1.0) 
         ENDIF 
!                                                                       
         CALL draw_marker (xleg - 2.0 * xh, yleg, imarktyp (iwin,       &
         iframe, ikurv), imarkcol (iwin, iframe, ikurv), sizemark (iwin,&
         iframe, ikurv) )                                               
!                                                                       
         CALL PGSCI (foncol (iwin, iframe, 6) ) 
         CALL PGSCF (fon_id (iwin, iframe, 6) ) 
         CALL PGSLW (nint (linewid (iwin, iframe, 0) / 0.13) ) 
         CALL PGSCH (fonscal (iwin, iframe) * fonsize (iwin, iframe, 6) &
         / 120.0)                                                       
         xt = xleg - 5.0 * xh 
         yt = yleg - 0.3 * yh 
         il = len_str (text) 
         IF (il.gt.0) then 
            CALL PGPTXT (xt, yt, 0.0, 1.0, text (1:il) ) 
         ENDIF 
      ENDIF 
!                                                                       
  999 CONTINUE 
      ENDDO 
!                                                                       
!------ plot possible annotations                                       
!                                                                       
      eex (1) = pex (iwin, iframe, 1) 
      eex (2) = pex (iwin, iframe, 2) 
      eey (1) = pey (iwin, iframe, 1) 
      eey (2) = pey (iwin, iframe, 2) 
!                                                                       
      CALL PGSCI (foncol (iwin, iframe, 6) ) 
      CALL PGSLW (nint (linewid (iwin, iframe, 0) / 0.13) ) 
      CALL PGSCF (fon_id (iwin, iframe, 6) ) 
      CALL PGSCH (fonscal (iwin, iframe) * fonsize (iwin, iframe, 6)    &
      / 120.0)                                                          
      CALL PGQCS (4, xh, yh) 
      DO ipkt = 1, maxan 
      IF (antext (iwin, iframe, ipkt) (1:3) .ne.'OFF') then 
         il = len_str (antext (iwin, iframe, ipkt) ) 
         xt = antx (iwin, iframe, ipkt) 
         yt = anty (iwin, iframe, ipkt) 
         xp = anx (iwin, iframe, ipkt) 
         yp = any (iwin, iframe, ipkt) 
         CALL koor_shear (1, xp, yp) 
         CALL koor_shear (1, xt, yt) 
         CALL koor_log (1, xp, yp) 
         CALL koor_log (1, xt, yt) 
         IF (il.gt.0.and.inrect (eex, eey, xt, yt) ) then 
            sx = 0.05 * (xt - xp) 
            sy = 0.05 * (yt - yp) 
            IF (xt.ne.xp.or.yt.ne.yp) call PGARRO (xt - sx, yt - sy, xp,&
            yp)                                                         
            IF (anjust (iwin, iframe, ipkt) .eq.IF_LEFT) rj = 0.0 
            IF (anjust (iwin, iframe, ipkt) .eq.IF_RIGHT) rj = 1.0 
            IF (anjust (iwin, iframe, ipkt) .eq.IF_CENTRE) rj = 0.5 
            CALL PGPTXT (xt, yt, anangle (iwin, iframe, ipkt), rj,      &
            antext (iwin, iframe, ipkt) (1:il) )                        
         ENDIF 
      ENDIF 
      ENDDO 
!                                                                       
      END SUBROUTINE draw_werte                     
!****7******************************************************************
      SUBROUTINE draw_line (xpl, ypl, npkt, ikurv) 
!+                                                                      
!     The drawing and filling itself                                    
!-                                                                      
      USE koordinate_mod
!
      USE debug_mod 
      USE errlist_mod 
      USE kuplot_config 
      USE kuplot_mod 
!                                                                       
      IMPLICIT none 
!                                                                       
      INTEGER npkt, ikurv 
      REAL xpl (npkt), ypl (npkt) 
!                                                                       
      REAL xfill (maxarray), yfill (maxarray) 
      REAL x1, x2, y1, y2, yb 
      INTEGER ip, np 
!                                                                       
      CALL koor_shear (npkt, xpl, ypl) 
      CALL koor_log (npkt, xpl, ypl) 
!                                                                       
!------ Fill area if required                                           
!                                                                       
      IF (ifilltyp (iwin, iframe, ikurv) .gt.0) then 
         IF (fillrange (iwin, iframe, ikurv, 1) .eq. - 9999.) then 
            x1 = pex (iwin, iframe, 1) 
            x2 = pex (iwin, iframe, 2) 
         ELSE 
            x1 = fillrange (iwin, iframe, ikurv, 1) 
            x2 = fillrange (iwin, iframe, ikurv, 2) 
         ENDIF 
         IF (fillrange (iwin, iframe, ikurv, 3) .eq. - 9999.) then 
            y1 = pey (iwin, iframe, 1) 
            y2 = pey (iwin, iframe, 2) 
            yb = max (0.0, ey (iwin, iframe, 1) ) 
         ELSE 
            y1 = fillrange (iwin, iframe, ikurv, 3) 
            y2 = fillrange (iwin, iframe, ikurv, 4) 
            yb = fillrange (iwin, iframe, ikurv, 3) 
         ENDIF 
!                                                                       
         np = 1 
         xfill (np) = x1 
         yfill (np) = yb 
!                                                                       
         DO ip = 1, npkt 
         IF (xpl (ip) .ge.x1.and.xpl (ip) .le.x2) then 
            np = np + 1 
            xfill (np) = xpl (ip) 
            IF (ypl (ip) .ge.y1.and.ypl (ip) .le.y2) then 
               yfill (np) = ypl (ip) 
            ELSE 
               IF (ypl (ip) .lt.y1) yfill (np) = y1 
               IF (ypl (ip) .gt.y2) yfill (np) = y2 
            ENDIF 
         ENDIF 
         ENDDO 
!                                                                       
         np = np + 1 
         xfill (np) = x2 
         yfill (np) = yb 
!                                                                       
         IF (dbg) then 
            WRITE ( *, 1000) 1, xfill (1), yfill (1) 
            WRITE ( *, 1000) 2, xfill (2), yfill (2) 
            WRITE ( *, 1000) np - 1, xfill (np - 1), yfill (np - 1) 
            WRITE ( *, 1000) np, xfill (np), yfill (np) 
 1000 FORMAT         (' debug  > Fill point ',i5,' at ',2(g12.4,1x)) 
         ENDIF 
!                                                                       
         IF (ifilltyp (iwin, iframe, ikurv) .eq.1.or.ifilltyp (iwin,    &
         iframe, ikurv) .eq.5) then                                     
            CALL PGSFS (1) 
            CALL PGSHS (45.0, 1.0, 0.0) 
         ELSEIF (ifilltyp (iwin, iframe, ikurv) .eq.2.or.ifilltyp (iwin,&
         iframe, ikurv) .eq.6) then                                     
            CALL PGSFS (3) 
            CALL PGSHS (45.0, 1.0, 0.0) 
         ELSEIF (ifilltyp (iwin, iframe, ikurv) .eq.3.or.ifilltyp (iwin,&
         iframe, ikurv) .eq.7) then                                     
            CALL PGSFS (3) 
            CALL PGSHS ( - 45.0, 1.0, 0.0) 
         ELSEIF (ifilltyp (iwin, iframe, ikurv) .eq.4.or.ifilltyp (iwin,&
         iframe, ikurv) .eq.8) then                                     
            CALL PGSFS (4) 
            CALL PGSHS (45.0, 1.0, 0.0) 
         ENDIF 
!                                                                       
         CALL PGSCI (ifillcol (iwin, iframe, ikurv) ) 
         CALL PGPOLY (np, xfill, yfill) 
!                                                                       
         IF (ifilltyp (iwin, iframe, ikurv) .gt.4) then 
            CALL PGSCI (6) 
            CALL PGSFS (2) 
            CALL PGPOLY (np, xfill, yfill) 
         ENDIF 
      ENDIF 
!                                                                       
!------ Draw the line                                                   
!                                                                       
      IF (abs(ilinetyp(iwin, iframe, ikurv)) .gt.0) then 
         CALL PGSLS (abs(ilinetyp (iwin, iframe, ikurv)) ) 
         CALL PGSLW (nint (linewid (iwin, iframe, ikurv) / 0.13) ) 
         CALL PGSCI (ilinecol (iwin, iframe, ikurv) ) 
         CALL PGLINE (npkt, xpl, ypl) 
      ENDIF 
!                                                                       
      END SUBROUTINE draw_line                      
!****7******************************************************************
      SUBROUTINE draw_poly (xpl, ypl, npkt, ikurv) 
!+                                                                      
!     The drawing and filling itself                                    
!-                                                                      
      USE koordinate_mod
!
      USE debug_mod 
      USE errlist_mod 
      USE kuplot_config 
      USE kuplot_mod 
!                                                                       
      IMPLICIT none 
!                                                                       
      INTEGER npkt, ikurv 
      REAL xpl (npkt), ypl (npkt) 
!                                                                       
      REAL xfill (maxarray), yfill (maxarray) 
      REAL x1, x2, y1, y2, yb 
      INTEGER ip, np 
!                                                                       
      CALL koor_shear (npkt, xpl, ypl) 
      CALL koor_log (npkt, xpl, ypl) 
!                                                                       
!------ Fill area if required                                           
!                                                                       
      IF (ifilltyp (iwin, iframe, ikurv) .gt.0) then 
         IF (fillrange (iwin, iframe, ikurv, 1) .eq. - 9999.) then 
            x1 = pex (iwin, iframe, 1) 
            x2 = pex (iwin, iframe, 2) 
         ELSE 
            x1 = fillrange (iwin, iframe, ikurv, 1) 
            x2 = fillrange (iwin, iframe, ikurv, 2) 
         ENDIF 
         IF (fillrange (iwin, iframe, ikurv, 3) .eq. - 9999.) then 
            y1 = pey (iwin, iframe, 1) 
            y2 = pey (iwin, iframe, 2) 
            yb = max (0.0, ey (iwin, iframe, 1) ) 
         ELSE 
            y1 = fillrange (iwin, iframe, ikurv, 3) 
            y2 = fillrange (iwin, iframe, ikurv, 4) 
            yb = fillrange (iwin, iframe, ikurv, 3) 
         ENDIF 
!                                                                       
         np = 0 
!        xfill (np) = x1 
!        yfill (np) = yb 
!                                                                       
         DO ip = 1, npkt 
         IF (xpl (ip) .ge.x1.and.xpl (ip) .le.x2) then 
            np = np + 1 
            xfill (np) = xpl (ip) 
            IF (ypl (ip) .ge.y1.and.ypl (ip) .le.y2) then 
               yfill (np) = ypl (ip) 
            ELSE 
               IF (ypl (ip) .lt.y1) yfill (np) = y1 
               IF (ypl (ip) .gt.y2) yfill (np) = y2 
            ENDIF 
         ENDIF 
         ENDDO 
!                                                                       
!        np = np + 1 
!        xfill (np) = x2 
!        yfill (np) = yb 
!                                                                       
         IF (dbg) then 
            WRITE ( *, 1000) 1, xfill (1), yfill (1) 
            WRITE ( *, 1000) 2, xfill (2), yfill (2) 
            WRITE ( *, 1000) np - 1, xfill (np - 1), yfill (np - 1) 
            WRITE ( *, 1000) np, xfill (np), yfill (np) 
 1000 FORMAT         (' debug  > Fill point ',i5,' at ',2(g12.4,1x)) 
         ENDIF 
!                                                                       
         IF (ifilltyp (iwin, iframe, ikurv) .eq.1.or.ifilltyp (iwin,    &
         iframe, ikurv) .eq.5) then                                     
            CALL PGSFS (1) 
            CALL PGSHS (45.0, 1.0, 0.0) 
         ELSEIF (ifilltyp (iwin, iframe, ikurv) .eq.2.or.ifilltyp (iwin,&
         iframe, ikurv) .eq.6) then                                     
            CALL PGSFS (3) 
            CALL PGSHS (45.0, 1.0, 0.0) 
         ELSEIF (ifilltyp (iwin, iframe, ikurv) .eq.3.or.ifilltyp (iwin,&
         iframe, ikurv) .eq.7) then                                     
            CALL PGSFS (3) 
            CALL PGSHS ( - 45.0, 1.0, 0.0) 
         ELSEIF (ifilltyp (iwin, iframe, ikurv) .eq.4.or.ifilltyp (iwin,&
         iframe, ikurv) .eq.8) then                                     
            CALL PGSFS (4) 
            CALL PGSHS (45.0, 1.0, 0.0) 
         ENDIF 
!                                                                       
         CALL PGSCI (ifillcol (iwin, iframe, ikurv) ) 
         CALL PGPOLY (np, xfill, yfill) 
!                                                                       
         IF (ifilltyp (iwin, iframe, ikurv) .gt.4) then 
            CALL PGSCI (6) 
            CALL PGSFS (2) 
            CALL PGPOLY (np, xfill, yfill) 
         ENDIF 
      ENDIF 
!                                                                       
!------ Draw the line                                                   
!                                                                       
      IF (abs(ilinetyp (iwin, iframe, ikurv)) .gt.0) then 
         CALL PGSLS (abs(ilinetyp (iwin, iframe, ikurv))) 
         CALL PGSLW (nint (linewid (iwin, iframe, ikurv) / 0.13) ) 
         CALL PGSCI (ilinecol (iwin, iframe, ikurv) ) 
         CALL PGLINE (npkt, xpl, ypl) 
      ENDIF 
!                                                                       
      END SUBROUTINE draw_poly                      
!****7******************************************************************
      SUBROUTINE draw_points (xpl, ypl, nnpkt, ikurv) 
!+                                                                      
!     Draws line as simple line, step or spline                         
!-                                                                      
      USE errlist_mod 
      USE kuplot_config 
      USE kuplot_mod 
!                                                                       
      IMPLICIT none 
!                                                                       
      INTEGER nnpkt, npkt, ikurv, is, ipkt 
      REAL xpl (nnpkt), ypl (nnpkt) 
      REAL eex (2), eey (2) 
      REAL y2a (maxarray), xhe (maxsp), yhe (maxsp) 
      REAL xxx, yyy, xst, xen, dxx 
      INTEGER :: i
      INTEGER :: ninterv
!                                                                       
      eex (1) = pex (iwin, iframe, 1) 
      eex (2) = pex (iwin, iframe, 2) 
      eey (1) = pey (iwin, iframe, 1) 
      eey (2) = pey (iwin, iframe, 2) 
!                                                                       
      npkt = nnpkt 
!                                                                       
!------ Normal line between points                                      
!                                                                       
      IF (abs(ilineart (iwin, iframe, ikurv)) .eq.1) then 
         IF (shear (iwin, iframe) .ne.90.0) then 
            CALL clip_array (eex, eey, xpl, ypl, npkt) 
         ENDIF 
         if(ilinetyp(iwin, iframe, ikurv)>=0) then
         CALL draw_line (xpl, ypl, npkt, ikurv) 
         else
         CALL draw_poly (xpl, ypl, npkt, ikurv) 
         endif
!                                                                       
!------ Histogram style steps                                           
!                                                                       
      ELSEIF (abs(ilineart (iwin, iframe, ikurv)) .eq.2) then 
         IF (4 * npkt.gt.maxsp) then 
            ier_num = - 13 
            ier_typ = ER_APPL 
         ELSE 
            is = 0 
            DO ipkt = 1, npkt - 1 
            xhe (is + 1) = xpl (ipkt) 
            yhe (is + 1) = ypl (ipkt) 
            xhe (is + 2) = xpl (ipkt) + (xpl (ipkt + 1) - xpl (ipkt) )  &
            / 2.0                                                       
            yhe (is + 2) = ypl (ipkt) 
            xhe (is + 3) = xpl (ipkt) + (xpl (ipkt + 1) - xpl (ipkt) )  &
            / 2.0                                                       
            yhe (is + 3) = ypl (ipkt + 1) 
            xhe (is + 4) = xpl (ipkt + 1) 
            yhe (is + 4) = ypl (ipkt + 1) 
            is = is + 4 
            ENDDO 
            IF (shear (iwin, iframe) .ne.90.0) then 
               CALL clip_array (eex, eey, xhe, yhe, npkt) 
            ENDIF 
            CALL draw_line (xhe, yhe, is, ikurv) 
         ENDIF 
!                                                                       
!------ Spline interpolation between points                             
!                                                                       
      ELSEIF (abs(ilineart (iwin, iframe, ikurv)) .eq.3) then 
         IF (npkt.gt.maxsp) then 
            ier_num = - 13 
            ier_typ = ER_APPL 
         ELSE 
            CALL spline (xpl, ypl, npkt, 1e30, 1e30, y2a) 
            xst = max (xmin (ikurv), pex (iwin, iframe, 1) ) 
            xen = min (xmax (ikurv), pex (iwin, iframe, 2) ) 
            dxx = (xen - xst) / maxsp 
            is = 1
            ninterv = IABS( NINT( ((xen-dxx)-(xst+dxx))/dxx )) 
!           DO xxx = xst + dxx, xen - dxx, dxx 
            DO i = 0, ninterv
               xxx = (xst+dxx) + i * dxx
            CALL splint (xpl, ypl, y2a, npkt, xxx, yyy, ier_num)
            IF(ier_num /= 0) THEN
               ier_typ =ER_APPL
               RETURN 
            ENDIF
            xhe (is) = xxx 
            yhe (is) = yyy 
            is = is + 1 
            ENDDO 
            IF (shear (iwin, iframe) .ne.90.0) then 
               CALL clip_array (eex, eey, xhe, yhe, is - 1) 
            ENDIF 
            CALL draw_line (xhe, yhe, is - 1, ikurv) 
         ENDIF 
      ENDIF 
!                                                                       
      END SUBROUTINE draw_points                    
!****7******************************************************************
      SUBROUTINE draw_marker (px, py, ityp, icol, size) 
!+                                                                      
!     Draws marker of typ ityp and color icol at px,py                  
!-                                                                      
      USE koordinate_mod
      USE kuplot_config 
      USE kuplot_mod 
!                                                                       
      IMPLICIT none 
!                                                                       
      CHARACTER(40) xstr, ystr 
      REAL px, py, size 
      REAL eex (2), eey (2) 
      REAL xh, yh, tx, ty 
      INTEGER ix, iy, ityp, icol, iptyp 
      LOGICAL inrect 
!                                                                       
      eex (1) = pex (iwin, iframe, 1) 
      eex (2) = pex (iwin, iframe, 2) 
      eey (1) = pey (iwin, iframe, 1) 
      eey (2) = pey (iwin, iframe, 2) 
      IF (ityp.eq.0.or..not.inrect (eex, eey, px, py) ) return 
!                                                                       
      CALL koor_shear (1, px, py) 
      CALL koor_log (1, px, py) 
!                                                                       
      CALL PGSLW (nint (linewid (iwin, iframe, 0) / 0.13) ) 
      CALL PGSCI (icol) 
!                                                                       
!------ just a dot                                                      
!                                                                       
      IF (ityp.eq.1) then 
         CALL PGPT (1, px, py, - 1) 
!                                                                       
!------ empty/filled circle                                             
!                                                                       
      ELSEIF (ityp.eq.2) then 
         CALL PGSCH (7.143 * size) 
         CALL PGPT (1, px, py, 21) 
      ELSEIF (ityp.eq.3) then 
         CALL PGSCH (8.333 * size) 
         CALL PGPT (1, px, py, 17) 
!                                                                       
!------ empty/filled square                                             
!                                                                       
      ELSEIF (ityp.eq.4) then 
         CALL PGSCH (5.000 * size) 
         CALL PGPT (1, px, py, 6) 
      ELSEIF (ityp.eq.5) then 
         CALL PGSCH (8.333 * size) 
         CALL PGPT (1, px, py, 16) 
!                                                                       
!------ emtyp triangle                                                  
!                                                                       
      ELSEIF (ityp.eq.6) then 
         CALL PGSCH (5.000 * size) 
         CALL PGPT (1, px, py, 7) 
!                                                                       
!------ cross 'x'                                                       
!                                                                       
      ELSEIF (ityp.eq.7) then 
         CALL PGSCH (5.556 * size) 
         CALL PGPT (1, px, py, 5) 
!                                                                       
!------ cross '+'                                                       
!                                                                       
      ELSEIF (ityp.eq.8) then 
         CALL PGSCH (4.545 * size) 
         CALL PGPT (1, px, py, 2) 
!                                                                       
!------ line '|'                                                        
!                                                                       
      ELSEIF (ityp.eq.9) then 
         CALL PGSCH (2.083 * size) 
         CALL PGPT (1, px, py, ichar ('|') ) 
!                                                                       
!------ line '/'                                                        
!                                                                       
      ELSEIF (ityp.eq.10) then 
         CALL PGSCH (1.852 * size) 
         CALL PGPT (1, px, py, ichar ('/') ) 
!                                                                       
!------ line '\'                                                        
!                                                                       
      ELSEIF (ityp.eq.11) then 
         CALL PGSCH (2.500 * size) 
         CALL PGPT (1, px, py, ichar ('\') ) 
!                                                                       
!------ line '-'                                                        
!                                                                       
      ELSEIF (ityp.eq.12) then 
         CALL PGSCH (3.571 * size) 
         CALL PGPT (1, px, py, ichar ('-') ) 
!                                                                       
!------ line from y-axis                                                
!                                                                       
      ELSEIF (ityp.eq.13) then 
         CALL PGMOVE (px, ey (iwin, iframe, 1) ) 
         CALL PGDRAW (px, py) 
!                                                                       
!------ marker types 100+PGPLOT marker number                           
!                                                                       
      ELSEIF (ityp.ge.92) then 
         iptyp = ityp - 100 
         CALL PGSCH (5.000 * size) 
         CALL PGPT (1, px, py, iptyp) 
!                                                                       
!------ cross '+' with coordinates                                      
!                                                                       
      ELSEIF (ityp.lt.0) then 
         CALL PGSCI (foncol (iwin, iframe, 6) ) 
         CALL PGSCH (4.545 * size) 
         CALL PGPT (1, px, py, 2) 
!                                                                       
         CALL PGSLW (nint (linewid (iwin, iframe, 0) / 0.13) ) 
         CALL PGSCF (fon_id (iwin, iframe, 6) ) 
         CALL PGSCH (fonscal (iwin, iframe) * fonsize (iwin, iframe, 6) &
         / 120.0)                                                       
         CALL PGQCS (4, xh, yh) 
         CALL realtostr (px, xstr, ix) 
         CALL realtostr (py, ystr, iy) 
         IF (ityp.eq. - 1.or.ityp.eq. - 3) then 
            tx = px + 0.3 * xh 
            ty = py + 0.6 * yh 
            CALL PGPTXT (tx, ty, 90.0, 0.0, xstr (1:ix) ) 
         ENDIF 
         IF (ityp.eq. - 2.or.ityp.eq. - 3) then 
            tx = px + 0.4 * xh 
            ty = py - 0.3 * yh 
            CALL PGPTXT (tx, ty, 0.0, 0.0, ystr (1:iy) ) 
         ENDIF 
      ENDIF 
!                                                                       
      END SUBROUTINE draw_marker                    
!****7******************************************************************
      SUBROUTINE draw_bonds (ikurv) 
!+                                                                      
!     Draws bonds                                                       
!-                                                                      
      USE errlist_mod 
      USE koordinate_mod
      USE kuplot_config 
      USE kuplot_mod 
!                                                                       
      IMPLICIT none 
!                                                                       
      INTEGER ikurv 
      INTEGER ib, iii, jjj, ipkt, jpkt 
      REAL xpl (2), ypl (2), dist 
      REAL eex (2), eey (2) 
      REAL blow, bhigh 
!                                                                       
      LOGICAL inrect 
!                                                                       
      eex (1) = pex (iwin, iframe, 1) 
      eex (2) = pex (iwin, iframe, 2) 
      eey (1) = pey (iwin, iframe, 1) 
      eey (2) = pey (iwin, iframe, 2) 
!                                                                       
      DO ib = 1, maxbond 
      IF (bond_rad (iwin, iframe, ib) .gt.0.0) then 
         blow = bond_rad (iwin, iframe, ib) * (1 - bond_sig (iwin,      &
         iframe, ib) )                                                  
         bhigh = bond_rad (iwin, iframe, ib) * (1 + bond_sig (iwin,     &
         iframe, ib) )                                                  
         CALL PGSLS (bond_ltyp (iwin, iframe, ib) ) 
         CALL PGSLW (nint (bond_lwid (iwin, iframe, ib) / 0.13) ) 
         CALL PGSCI (bond_lcol (iwin, iframe, ib) ) 
         DO ipkt = 1, lenc(ikurv) 
         iii = offxy (ikurv - 1) + ipkt 
         IF (inrect (eex, eey, x (iii), y (iii) ) ) then 
            DO jpkt = 1, lenc(ikurv) 
            jjj = offxy (ikurv - 1) + jpkt 
            xpl (1) = x (iii) 
            ypl (1) = y (iii) 
            xpl (2) = x (jjj) 
            ypl (2) = y (jjj) 
            CALL koor_shear (2, xpl, ypl) 
            CALL koor_log (2, xpl, ypl) 
            ypl (1) = yskal (iwin, iframe) * ypl (1) 
            ypl (2) = yskal (iwin, iframe) * ypl (2) 
            dist = sqrt ( (xpl (2) - xpl (1) ) **2 + (ypl (2) - ypl (1) &
            ) **2)                                                      
            IF (dist.ge.blow.and.dist.le.bhigh) then 
               xpl (1) = x (iii) 
               ypl (1) = y (iii) 
               xpl (2) = x (jjj) 
               ypl (2) = y (jjj) 
               CALL koor_shear (2, xpl, ypl) 
               CALL koor_log (2, xpl, ypl) 
               CALL PGLINE (2, xpl, ypl) 
            ENDIF 
            ENDDO 
         ENDIF 
         ENDDO 
      ENDIF 
      ENDDO 
!                                                                       
      END SUBROUTINE draw_bonds                     
!*****7*****************************************************************
      SUBROUTINE draw_bitmap (ik) 
!+                                                                      
!     This routine draws the bitmaps                                    
!-                                                                      
      USE debug_mod 
      USE errlist_mod 
      USE prompt_mod 
      USE wink_mod
      USE kuplot_config 
      USE kuplot_mod 
!                                                                       
      IMPLICIT none 
!                                                                       
      REAL zpl (maxz, maxz), tr (6) 
      REAL rc, rdx, rdy, yf 
      REAL zzmin, zzmax, x1, x2 
      INTEGER i, ic, ik, ix, iy, ikk 
      INTEGER nx_min, nx_max, ny_min, ny_max 
!                                                                       
      rdx = (xmax (ik) - xmin (ik) ) / REAL(nx (ik) - 1) 
      rdy = (ymax (ik) - ymin (ik) ) / REAL(ny (ik) - 1) 
      nx_min = max (1, nint ( (pex (iwin, iframe, 1) - xmin (ik) )      &
      / rdx) + 1)                                                       
      nx_max = min (nx (ik), nint ( (pex (iwin, iframe, 2) - xmin (ik) )&
      / rdx) + 1)                                                       
      ny_min = max (1, nint ( (pey (iwin, iframe, 1) - ymin (ik) )      &
      / rdy) + 1)                                                       
      ny_max = min (ny (ik), nint ( (pey (iwin, iframe, 2) - ymin (ik) )&
      / rdy) + 1)                                                       
!                                                                       
      IF (nx_min.eq.nx_max.or.ny_min.eq.ny_max) return 
!                                                                       
      IF (nx_max.gt.maxz.or.ny_max.gt.maxz) then 
         ier_num = - 34 
         ier_typ = ER_APPL 
         RETURN 
      ENDIF 
!                                                                       
!------ Copy array and call drawing routine                             
!                                                                       
      zzmin = z_min (iwin, iframe, 1) 
      zzmax = nz (iwin, iframe, 1) * z_inc (iwin, iframe, 1) + z_min (  &
      iwin, iframe, 1)                                                  
      DO iy = ny_min, ny_max 
      DO ix = nx_min, nx_max 
      ikk = offz (ik - 1) + (ix - 1) * ny (ik) + iy 
      zpl (ix, iy) = z (ikk) 
      ENDDO 
      ENDDO 
!                                                                       
!------ Transformation matrix                                           
!                                                                       
      tr (1) = xmin (ik) - rdx 
      tr (2) = rdx 
      tr (3) = 0.0 
      tr (4) = ymin (ik) - rdy 
      tr (5) = 0.0 
      tr (6) = rdy 
!                                                                       
      IF (shear (iwin, iframe) .ne.90.0) then 
         yf = yskal (iwin, iframe) / tan (REAL(rad) * shear (iwin, iframe) ) 
         tr (1) = tr (1) + yf * (tr (4) - ey (iwin, iframe, 1) ) 
         tr (3) = rdy * yf 
!                                                                       
!------ - Account for pixel overlap for ANGL not 90.                    
!                                                                       
         x1 = tr (1) + tr (2) * (REAL(nx_min) - 0.5) + tr (3) *       &
         (REAL(ny_min) - 0.5)                                         
         x2 = tr (1) + tr (2) * (REAL(nx_max) + 0.5) + tr (3) *       &
         (REAL(ny_min) + 0.5)                                         
         pex (iwin, iframe, 1) = x1 
         pex (iwin, iframe, 2) = x2 
      ENDIF 
!                                                                       
      IF (dbg) then 
         WRITE (output_io, 8000) rdx, rdy 
         WRITE (output_io, 8010) nx_min, nx_max, ny_min, ny_max 
         WRITE (output_io, 8020) tr 
 8000 FORMAT      (' debug > dx,dy               = ',2(g13.6,2x)) 
 8010 FORMAT      (' debug > nx(mi,ma),ny(mi,ma) = ',2i5,2x,2i5) 
 8020 FORMAT      (' debug > TR matrix           = ',3(g13.6,2x),/,     &
     &                    '                               ',3(g13.6,2x))
      ENDIF 
!                                                                       
!------ Setting colours                                                 
!                                                                       
      IF(ABS(col_map_type(iwin,iframe))==COL_MAP_GRAY) THEN
         CALL cmap_gray(.FALSE.)
      ELSEIF(ABS(col_map_type(iwin,iframe))==COL_MAP_FIRE) THEN
         CALL cmap_fire(.FALSE.)
      ELSEIF(ABS(col_map_type(iwin,iframe))==COL_MAP_ICE ) THEN
         CALL cmap_ice(.FALSE.)
      ELSEIF(ABS(col_map_type(iwin,iframe))==COL_MAP_KUPL) THEN
         CALL cmap_kupl(.FALSE.)
      ELSEIF(ABS(col_map_type(iwin,iframe))==COL_MAP_THER) THEN
         CALL cmap_thermal(zzmin, zzmax, .FALSE.)
      ELSEIF(ABS(col_map_type(iwin,iframe))==COL_MAP_PDF) THEN
         CALL cmap_pdf(zzmin, zzmax, .FALSE.)
      ENDIF
      IF(col_map_type(iwin,iframe)<0) THEN
         CALL cmap_invert(.FALSE.)
      ENDIF
      CALL PGQCIR (ix, iy) 
      IF (ix.lt. (iaf (iwin) + 19) ) ix = iaf (iwin) + 19 
      CALL PGSCIR (ix, iy) 
      rc = 254.0 / REAL(iy - ix) 
      DO i = ix, iy 
      ic = int ( (i - ix) * rc) + 1 
      CALL PGSCR (i, col_map (iwin, ic, 1), col_map (iwin, ic, 2),      &
      col_map (iwin, ic, 3) )                                           
      ENDDO 
!                                                                       
      CALL PGSITF (0) 
      CALL PGIMAG (zpl, maxz, maxz, nx_min, nx_max, ny_min, ny_max,     &
      zzmin, zzmax, tr)                                                 
!                                                                       
!------ Plot wedge next to plot if required                             
!                                                                       
      IF (achse (iwin, iframe, 3) (1:3) .ne.'OFF') then 
         CALL PGSCI (foncol (iwin, iframe, 4) ) 
         CALL PGSCF (fon_id (iwin, iframe, 4) ) 
         CALL PGSCH (fonscal (iwin, iframe) * fonsize (iwin, iframe, 4) &
         / 120.0)                                                       
         CALL PGWEDG ('RI', 0.5, 2.5, zzmin, zzmax, achse (iwin, iframe,&
         3) )                                                           
      ENDIF 
!                                                                       
      END SUBROUTINE draw_bitmap                    
!*****7*****************************************************************
      SUBROUTINE draw_contour (ik) 
!+                                                                      
!     This routine draws the contour lines                              
!-                                                                      
      USE debug_mod 
      USE errlist_mod 
      USE prompt_mod 
      USE wink_mod
      USE kuplot_config 
      USE kuplot_mod 
!                                                                       
      IMPLICIT none 
!                                                                       
      CHARACTER(25) label 
      REAL zpl (maxz, maxz), tr (6) 
      REAL rdx, rdy, zm, zi, h, log10, yf 
      INTEGER il, ic, ik, ix, iy, ikk, ihp, ihl, lmi, lin 
      INTEGER nx_min, nx_max, ny_min, ny_max 
!                                                                       
      rdx = (xmax (ik) - xmin (ik) ) / REAL(nx (ik) - 1) 
      rdy = (ymax (ik) - ymin (ik) ) / REAL(ny (ik) - 1) 
      nx_min = max (1, nint ( (pex (iwin, iframe, 1) - xmin (ik) )      &
      / rdx) + 1)                                                       
      nx_max = min (nx (ik), nint ( (pex (iwin, iframe, 2) - xmin (ik) )&
      / rdx) + 1)                                                       
      ny_min = max (1, nint ( (ey (iwin, iframe, 1) - ymin (ik) )       &
      / rdy) + 1)                                                       
      ny_max = min (ny (ik), nint ( (ey (iwin, iframe, 2) - ymin (ik) ) &
      / rdy) + 1)                                                       
!                                                                       
      IF (nx_min.eq.nx_max.or.ny_min.eq.ny_max) return 
!                                                                       
      lmi = nint (0.2 * min ( (nx_max - nx_min), (ny_max - ny_min) ) ) 
      lin = 3 * lmi 
      ic = iframe+20 
!                                                                       
      log10 = log (10.0) 
!                                                                       
      IF (nx_max.gt.maxz.or.ny_max.gt.maxz) then 
         ier_num = - 34 
         ier_typ = ER_APPL 
         RETURN 
      ENDIF 
!                                                                       
!------ Copy array and call drawing routine                             
!                                                                       
      DO iy = ny_min, ny_max 
      DO ix = nx_min, nx_max 
      ikk = offz (ik - 1) + (ix - 1) * ny (ik) + iy 
      zpl (ix, iy) = z (ikk) 
      ENDDO 
      ENDDO 
!                                                                       
!------ Transformation matrix                                           
!                                                                       
      tr (1) = xmin (ik) - rdx 
      tr (2) = rdx 
      tr (3) = 0.0 
      tr (4) = ymin (ik) - rdy 
      tr (5) = 0.0 
      tr (6) = rdy 
!                                                                       
      IF (shear (iwin, iframe) .ne.90.0) then 
         yf = yskal (iwin, iframe) / tan (REAL(rad) * shear (iwin, iframe) ) 
         tr (1) = tr (1) + yf * (tr (4) - ey (iwin, iframe, 1) ) 
         tr (3) = rdy * yf 
      ENDIF 
!                                                                       
      IF (dbg) then 
         WRITE (output_io, 8000) rdx, rdy 
         WRITE (output_io, 8010) nx_min, nx_max, ny_min, ny_max 
         WRITE (output_io, 8020) tr 
 8000 FORMAT      (' debug > dx,dy               = ',2(g13.6,2x)) 
 8010 FORMAT      (' debug > nx(mi,ma),ny(mi,ma) = ',2i5,2x,2i5) 
 8020 FORMAT      (' debug > TR matrix           = ',3(g13.6,2x),/,     &
     &                    '                               ',3(g13.6,2x))
      ENDIF 
!                                                                       
!------ Drawing the contours                                            
!                                                                       
      DO ihp = 1, iho (iwin, iframe) 
      CALL PGSLS (hlinetyp (iwin, iframe, ik, ihp) ) 
      CALL PGSLW (nint (linewid (iwin, iframe, ik) / 0.13) ) 
      zm = z_min (iwin, iframe, ihp) 
      zi = z_inc (iwin, iframe, ihp) 
      IF (hlinetyp (iwin, iframe, ik, ihp) .ne.0) then 
         DO ihl = 1, nz (iwin, iframe, ihp) 
         h = zm + (ihl - 1) * zi 
         CALL PGSCI (hlinecol (iwin, iframe, ik, ihp) ) 
         CALL PGCONT (zpl, maxz, maxz, nx_min, nx_max, ny_min, ny_max,  &
         h, - 1, tr)                                                    
!                                                                       
         IF (hlabel (iwin, iframe, ik) .ne.0) then 
            IF (mod (ihl, hlabel (iwin, iframe, ik) ) .eq.0) then 
               CALL realtostr (h, label, il) 
               CALL PGSLS (1) 
               CALL PGSLW (nint (linewid (iwin, iframe, 0) / 0.13) ) 
               CALL PGSCI (foncol (iwin, iframe, 4) ) 
               CALL PGSCF (fon_id (iwin, iframe, 4) ) 
               CALL PGSCH (fonscal (iwin, iframe) * fonsize (iwin,      &
               iframe, 4) / 160.0)                                      
               CALL PGSTBG (ic) 
               CALL PGCONL (zpl, maxz, maxz, nx_min, nx_max, ny_min,    &
               ny_max, h, tr, label, lin, lmi)                          
               CALL PGSTBG ( - 1) 
            ENDIF 
         ENDIF 
         ENDDO 
      ENDIF 
      ENDDO 
!                                                                       
      END SUBROUTINE draw_contour                   
!*****7*****************************************************************
      SUBROUTINE realtostr (r, str, ll) 
!+                                                                      
!     Converts real number to string                                    
!-                                                                      
      IMPLICIT none 
!                                                                       
      CHARACTER ( * ) str 
      REAL r 
      INTEGER ll, pp, mm 
!                                                                       
      IF (r.eq.0) then 
         pp = 1 
         mm = 0 
      ELSE 
         pp = int (log (abs (r) ) / log (10.0) ) - 4 
         IF (pp.lt.0) then 
            mm = nint (r * 10** (abs (pp) ) ) 
         ELSE 
            mm = nint (r / 10** (abs (pp) ) ) 
         ENDIF 
      ENDIF 
!                                                                       
      CALL PGNUMB (mm, pp, 0, str, ll) 
!                                                                       
      END SUBROUTINE realtostr                      
!*****7*****************************************************************
      LOGICAL function k_in_f (ik) 
!+                                                                      
!     Check if data set ik should be plotted in actual frame            
!-                                                                      
      USE kuplot_config 
      USE kuplot_mod 
!                                                                       
      IMPLICIT none 
!                                                                       
      INTEGER ik, i 
!                                                                       
      k_in_f = .false. 
      DO i = 1, maxkurvtot 
      k_in_f = k_in_f.or. (infra (iwin, iframe, i) .eq.ik) 
      ENDDO 
!                                                                       
      END FUNCTION k_in_f                           
!*****7*****************************************************************
      LOGICAL function inrect (ex, ey, x, y) 
!-                                                                      
!       Checks if point x,y is within rectangle ex(2),ey(2).            
!+                                                                      
      REAL x, y, ex (2), ey (2) 
!                                                                       
      inrect = (x.ge.ex (1) .and.x.le.ex (2) ) .and. (y.ge.ey (1)       &
      .and.y.le.ey (2) )                                                
      END FUNCTION inrect                           
!***********************************************************************
      SUBROUTINE clip_array (ex, ey, ax, ay, npkt) 
!                                                                       
      IMPLICIT none 
!                                                                       
      INTEGER i, j, npkt 
      LOGICAL inrect, in1, in2, in3, sf 
      REAL ex (2), ey (2), ax (npkt), ay (npkt) 
      REAL xs, ys, axp, ayp 
!                                                                       
      IF (npkt.lt.2) return 
!                                                                       
      j = 1 
!                                                                       
!------ First point                                                     
!                                                                       
      in2 = inrect (ex, ey, ax (1), ay (1) ) 
      in3 = inrect (ex, ey, ax (2), ay (2) ) 
      IF (in2) then 
         ax (j) = ax (1) 
         ay (j) = ay (1) 
         j = j + 1 
      ELSEIF (.not.in2.and.in3) then 
         CALL schnitt (ex, ey, ax (1), ay (1), ax (2), ay (2), xs, ys,  &
         sf)                                                            
         ax (j) = xs 
         ay (j) = ys 
         j = j + 1 
      ENDIF 
!                                                                       
      axp = ax (1) 
      ayp = ay (1) 
!                                                                       
!------ Check supsequent points                                         
!                                                                       
      DO i = 2, npkt - 1 
      in1 = inrect (ex, ey, axp, ayp) 
      in2 = inrect (ex, ey, ax (i), ay (i) ) 
      in3 = inrect (ex, ey, ax (i + 1), ay (i + 1) ) 
      IF (in2) then 
         axp = ax (i) 
         ayp = ay (i) 
         ax (j) = ax (i) 
         ay (j) = ay (i) 
         j = j + 1 
      ELSEIF (.not.in2.and.in1) then 
         CALL schnitt (ex, ey, axp, ayp, ax (i), ay (i), xs, ys, sf) 
         axp = ax (i) 
         ayp = ay (i) 
         ax (j) = xs 
         ay (j) = ys 
         j = j + 1 
      ELSEIF (.not.in2.and.in3) then 
         CALL schnitt (ex, ey, ax (i), ay (i), ax (i + 1), ay (i + 1),  &
         xs, ys, sf)                                                    
         axp = ax (i) 
         ayp = ay (i) 
         ax (j) = xs 
         ay (j) = ys 
         j = j + 1 
      ELSE 
         axp = ax (i) 
         ayp = ay (i) 
         CONTINUE 
      ENDIF 
      ENDDO 
!                                                                       
!------ Last point                                                      
!                                                                       
      in1 = inrect (ex, ey, axp, ayp) 
      in2 = inrect (ex, ey, ax (npkt), ay (npkt) ) 
      IF (in2) then 
         ax (j) = ax (npkt) 
         ay (j) = ay (npkt) 
         j = j + 1 
      ELSEIF (.not.in2.and.in1) then 
         CALL schnitt (ex, ey, axp, ayp, ax (npkt), ay (npkt), xs, ys,  &
         sf)                                                            
         ax (j) = xs 
         ay (j) = ys 
         j = j + 1 
      ENDIF 
!                                                                       
      npkt = j - 1 
!                                                                       
      END SUBROUTINE clip_array                     
!***********************************************************************
      SUBROUTINE schnitt (ex, ey, x1, y1, x2, y2, xs, ys, sflag) 
!-                                                                      
!           Calculates the intersection of the line from x1,y1 to x2,y2 
!           with rectangle given by ex(2),ey(2)                         
!+                                                                      
      IMPLICIT none 
!                                                                       
      REAL ex (2), ey (2) 
      REAL xs, ys, x1, y1, x2, y2 
      REAL lam, mue 
      LOGICAL sflag 
!                                                                       
      sflag = .false. 
      CALL cross (lam, mue, ex (1), ex (2), x1, x2, ey (1), ey (1),     &
      y1, y2, xs, ys)                                                   
      IF (0..le.lam.and.lam.le.1..and.0..le.mue.and.mue.le.1.) then 
         sflag = .true. 
         RETURN 
      ENDIF 
      CALL cross (lam, mue, ex (2), ex (2), x1, x2, ey (1), ey (2),     &
      y1, y2, xs, ys)                                                   
      IF (0..le.lam.and.lam.le.1..and.0..le.mue.and.mue.le.1.) then 
         sflag = .true. 
         RETURN 
      ENDIF 
      CALL cross (lam, mue, ex (2), ex (1), x1, x2, ey (2), ey (2),     &
      y1, y2, xs, ys)                                                   
      IF (0..le.lam.and.lam.le.1..and.0..le.mue.and.mue.le.1.) then 
         sflag = .true. 
         RETURN 
      ENDIF 
      CALL cross (lam, mue, ex (1), ex (1), x1, x2, ey (2), ey (1),     &
      y1, y2, xs, ys)                                                   
      IF (0..le.lam.and.lam.le.1..and.0..le.mue.and.mue.le.1.) then 
         sflag = .true. 
         RETURN 
      ENDIF 
      END SUBROUTINE schnitt                        
!*********************************************************************  
      SUBROUTINE cross (lam, mue, x1, x2, x3, x4, y1, y2, y3, y4, xs,   &
      ys)                                                               
!-                                                                      
!           Calculates intersection of two vectors : 12 and 34          
!           if singulaer lam and mue are set to  -20.                   
!+                                                                      
      IMPLICIT none 
!                                                                       
      REAL x1, x2, x3, x4, x21, x43, x13 
      REAL y1, y2, y3, y4, y21, y43, y13 
      REAL xs, ys, det 
      REAL lam, mue 
!                                                                       
      x21 = x2 - x1 
      x43 = x4 - x3 
      x13 = x1 - x3 
      y21 = y2 - y1 
      y43 = y4 - y3 
      y13 = y1 - y3 
      det = y43 * x21 - x43 * y21 
      IF (abs (det) .eq.0.0) then 
         lam = - 20. 
         mue = - 20. 
      ELSE 
         lam = (x21 * y13 - x13 * y21) / det 
         mue = (x43 * y13 - x13 * y43) / det 
         xs = x3 + lam * x43 
         ys = y3 + lam * y43 
      ENDIF 
      END SUBROUTINE cross                          
