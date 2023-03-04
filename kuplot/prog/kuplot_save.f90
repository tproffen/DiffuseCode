module kuplot_save_mod
!
!*****7**************************************************************** 
!     This file contains subroutines saving files                       
!*****7**************************************************************** 
!
contains
!
!*******************************************************************************
!
SUBROUTINE do_save_data (zeile, lp) 
!-                                                                      
!       Main menu for save command.                                     
!+                                                                      
USE build_name_mod
USE errlist_mod 
USE get_params_mod
USE kuplot_config 
use kuplot_math_mod
USE lib_help
USE precision_mod
USE str_comp_mod
!                                                                       
IMPLICIT none 
!                                                                       
INTEGER maxw 
PARAMETER (maxw = 10) 
!                                                                       
CHARACTER(len=*), intent(inout) :: zeile 
INTEGER         , intent(inout) :: lp 
!                                                                       
CHARACTER(LEN=PREC_STRING) :: cpara (maxw) 
INTEGER :: lpara (maxw) 
integer :: ianz 
real(KIND=PREC_DP) :: werte (maxw) 
!                                                                       
!                                                                       
CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
!                                                                       
IF (ianz.lt.2) then 
   ier_num = - 6 
   ier_typ = ER_COMM 
ELSE 
!                                                                       
!------ - 'gsas': Save data as GSAS file                                
!                                                                       
   IF (str_comp (cpara (1) , 'gsas', 2, lpara (1) , 4) ) then 
      CALL del_params (1, ianz, cpara, lpara, maxw) 
      CALL do_build_name (ianz, cpara, lpara, werte, maxw, 1) 
      IF (ier_num.ne.0) return 
      IF (ianz.lt.2) cpara (2) = 'none' 
      CALL do_save_gsas (cpara (1), cpara (2) ) 
!                                                                       
!------ - 'merge': Saves all current data as zz file                    
!                                                                       
   ELSEIF (str_comp (cpara (1) , 'merge', 2, lpara (1) , 5) ) then
      CALL del_params (1, ianz, cpara, lpara, maxw) 
      CALL do_build_name (ianz, cpara, lpara, werte, maxw, 1) 
      IF (ier_num.ne.0) return 
      CALL do_save_merge (cpara (1) ) 
!                                                                       
!------ - 'xml': Save current datasets as PlotML (XML) file             
!                                                                       
   ELSEIF (str_comp (cpara (1) , 'xml', 2, lpara (1) , 3) ) then 
      CALL del_params (1, ianz, cpara, lpara, maxw) 
      CALL do_build_name (ianz, cpara, lpara, werte, maxw, 1) 
      IF (ier_num.ne.0) return 
      IF (ianz.lt.2) cpara (2) = '/ptplot' 
      CALL do_save_xml (cpara (1), cpara (2) ) 
!                                                                       
!------ - Unknown command                                               
!                                                                       
   ELSE 
      ier_num = - 6 
      ier_typ = ER_COMM 
   ENDIF 
ENDIF 
!                                                                       
END SUBROUTINE do_save_data                   
!
!*****7**************************************************************** 
!
SUBROUTINE do_save_merge (filname) 
!-                                                                      
!       Saves all current data sets as 'x dataset y' file               
!+                                                                      
USE debug_mod 
USE errlist_mod 
USE prompt_mod 
USE kuplot_config 
USE kuplot_mod 
use kuplot_low_mod
USE support_mod
!                                                                       
IMPLICIT none 
!                                                                       
INTEGER ifil 
PARAMETER (ifil = 33) 
!                                                                       
CHARACTER(len= *), intent(in) :: filname
! 
INTEGER :: ikurv, ip 
!                                                                       
      !LOGICAL k_in_f 
!                                                                       
CALL oeffne (ifil, filname, 'unknown') 
IF (ier_num.ne.0) return 
!                                                                       
DO ikurv = 1, iz - 1 
   IF (k_in_f (ikurv) .and..not.lni (ikurv) ) then 
      DO ip = 1, lenc(ikurv) 
         WRITE (ifil, 1000) x (offxy (ikurv - 1) + ip), REAL(ikurv),  &
         y (offxy (ikurv - 1) + ip)                                     
      ENDDO 
   ENDIF 
ENDDO 
!                                                                       
 1000 FORMAT    (3(g14.6,1x)) 
!
END SUBROUTINE do_save_merge                  
!
!*****7**************************************************************** 
!
SUBROUTINE do_save_xml (filname, path) 
!-                                                                      
!       Saves all current data sets as PlotML file                      
!+                                                                      
USE debug_mod 
USE errlist_mod 
USE prompt_mod 
USE kuplot_config 
USE kuplot_mod 
use kuplot_low_mod
USE lib_length
USE support_mod
!                                                                       
IMPLICIT none 
!                                                                       
INTEGER ifil 
PARAMETER (ifil = 33) 
!                                                                       
CHARACTER(len= *), intent(in) :: filname
CHARACTER(len= *), intent(in) :: path
!
CHARACTER(len=20) :: connect, marks 
INTEGER :: ikurv, ip 
!                                                                       
!      LOGICAL k_in_f 
!                                                                       
CALL oeffne (ifil, filname, 'unknown') 
IF (ier_num.ne.0) return 
!                                                                       
WRITE (ifil, 1000) path (1:len_str (path) ) 
WRITE (ifil, 1100) 
WRITE (ifil, 1200) titel (iwin, iframe, 1) (1:len_str (titel (    &
iwin, iframe, 1) ) )                                              
WRITE (ifil, 1210) achse (iwin, iframe, 1) (1:len_str (achse (    &
iwin, iframe, 1) ) )                                              
WRITE (ifil, 1220) achse (iwin, iframe, 2) (1:len_str (achse (    &
iwin, iframe, 2) ) )                                              
WRITE (ifil, 1300) 
!                                                                       
DO ikurv = 1, iz - 1 
   IF (k_in_f (ikurv) .and..not.lni (ikurv) ) then 
      IF (ilinetyp (iwin, iframe, ikurv) .eq.0) then 
         connect = 'no' 
      ELSE 
         connect = 'yes' 
      ENDIF 
!                                                                       
      IF (imarktyp (iwin, iframe, ikurv) .ne.0) then 
            marks = 'marks="dots"' 
      ENDIF 
!                                                                       
      IF (ilegend (iwin, iframe, ikurv) .eq.1) then 
         WRITE (ifil, 1400) clegend (iwin, iframe, ikurv) (1:len_str &
         (clegend (iwin, iframe, ikurv) ) ), connect (1:len_str (    &
         connect) ), marks (1:len_str (marks) )                      
      ELSE 
         WRITE (ifil, 1400) fname (ikurv) (1:len_str (fname (ikurv))), &
            connect (1:len_str (connect) ), marks (1:len_str (marks))
      ENDIF 
      DO ip = 1, lenc(ikurv) 
         WRITE(ifil, 1500) x(offxy(ikurv - 1) + ip), y(offxy(ikurv - 1) + ip)
      ENDDO 
      WRITE (ifil, 2400) 
   ENDIF 
ENDDO 
!                                                                       
WRITE (ifil, 2100) 
CLOSE (ifil) 
!                                                                       
 1000 FORMAT('<?xml version="1.0" standalone="no"?>',/,                 &
     &       '<!DOCTYPE plot SYSTEM "',a,'/PlotML_1.dtd">',/,           &
     &       '<!-- Ptolemy plot, version 3.1, PlotML format. -->')      
 1100 FORMAT('<plot>') 
 1200 FORMAT('<title>',a,'</title>') 
 1210 FORMAT('<xLabel>',a,'</xLabel>') 
 1220 FORMAT('<yLabel>',a,'</yLabel>') 
 1300 FORMAT('<default connected="yes" marks="none"/>') 
 1400 FORMAT('<dataset name="',a,'" connected="',a,                     &
     &       '" ',a,'>')                                                
 1500 FORMAT('<p x="',e13.6,'" y="',e13.6,'"/>') 
 2100 FORMAT('</plot>') 
 2400 FORMAT('</dataset>') 
!                                                                       
END SUBROUTINE do_save_xml                    
!
!*****7**************************************************************** 
!
SUBROUTINE do_save_gsas (filname, iname) 
!-                                                                      
!       Saves all current data sets as GSAS file                        
!+                                                                      
USE debug_mod 
USE errlist_mod 
USE prompt_mod 
USE kuplot_config 
USE kuplot_mod 
USE lib_length
USE support_mod
!                                                                       
IMPLICIT none 
!                                                                       
INTEGER ifil 
PARAMETER (ifil = 33) 
!                                                                       
CHARACTER(len= *), intent(inout) :: filname
CHARACTER(len= *), intent(in) :: iname 
!                                                                       
CHARACTER(len=200) :: outstr 
CHARACTER(len=2) :: eol 
REAL   , dimension(:), allocatable :: gdat  !(MAXGSAS) 
REAL   , dimension(:), allocatable :: gsig  !(MAXGSAS) 
REAL :: delt 
INTEGER, dimension(:), allocatable :: tmap  !(MAXGSAS) 
INTEGER, dimension(:), allocatable :: tof  !(MAXGSAS) 
INTEGER :: ig, il, ixx 
INTEGER :: j, k, ibeg, ifin, imap, ioff, itof 
INTEGER :: ntmap, nrec, irec
LOGICAL::  lfile, lnew 
!                                                                       
allocate(gdat(MAXGSAS))
allocate(gsig(MAXGSAS))
allocate(tmap(MAXGSAS))
allocate( tof(MAXGSAS))
!                                                                       
eol = char (13) //char (10) 
irec = 1 
!                                                                       
INQUIRE (file = filname, exist = lfile) 
IF (lfile) call do_del_file (filname) 
OPEN (ifil, file = filname, status = 'new', access = 'direct',    &
      form = 'formatted', recl = 82, err = 999)                         
!                                                                       
WRITE (output_io, 5000) filname (1:len_str (filname) ) 
!                                                                       
!------ Write GSAS file header                                          
!                                                                       
il = len_str (titel (iwin, iframe, 1) ) 
WRITE (outstr, 1000) titel (iwin, iframe, 1) (1:il) 
WRITE (ifil, 1000, rec = irec) outstr (1:80) //eol 
irec = irec + 1 
IF (iname.ne.'none') then 
   WRITE (outstr, 1040) iname (1:len_str (iname) ) 
   WRITE (ifil, 1000, rec = irec) outstr (1:80) //eol 
   irec = irec + 1 
ENDIF 
WRITE (outstr, 1050) nint (1.0) 
WRITE (ifil, 1000, rec = irec) outstr (1:80) //eol 
irec = irec + 1 
!                                                                       
!------ Comment section                                                 
!                                                                       
WRITE (outstr, 1100) 
WRITE (ifil, 1000, rec = irec) outstr (1:80) //eol 
irec = irec + 1 
!                                                                       
WRITE (outstr, 1110) 'Created by KUPLOT ..' 
WRITE (ifil, 1000, rec = irec) outstr (1:80) //eol 
irec = irec + 1 
!                                                                       
WRITE (outstr, 1100) 
WRITE (ifil, 1000, rec = irec) outstr (1:80) //eol 
irec = irec + 1 
!                                                                       
!------ Loop over all data sets                                         
!                                                                       
itof = 0 
DO ig = 1, iz - 1 
!                                                                       
!------ - Make TOF array - convert back to bin boundaries               
!                                                                       
   ioff = offxy (ig - 1) 
   tof (1) = nint (10. * (x (ioff + 1) - tof_offset) ) 
   tof (2) = nint (10. * (x (ioff + 1) + tof_offset) ) 
   DO k = 3, lenc(ig) + 1 
      ixx = nint (20. * x (ioff + k - 1) ) 
      tof (k) = ixx - tof (k - 1) 
   ENDDO 
!                                                                       
!------ - Check if we need to write new time map entry                  
!                                                                       
   CALL make_timemap (ig, tof, tmap, ntmap, lnew) 
   IF (lnew) then 
      itof = itof + 1 
      nrec = (ntmap - 1) / 10 + 1 
      WRITE (outstr, 2000) 'TIME_MAP', itof, ntmap, nrec, ' TIME_MAP', 100
      WRITE (ifil, 1000, rec = irec) outstr (1:80) //eol 
      irec = irec + 1 
      ibeg = 1 
      DO j = 1, nrec 
         ifin = min (ibeg + 9, ntmap) 
         WRITE (outstr, 2100) (tmap (k), k = ibeg, ifin) 
         WRITE (ifil, 1000, rec = irec) outstr (1:80) //eol 
         irec = irec + 1 
         ibeg = ibeg + 10 
      ENDDO 
      imap = itof 
   ENDIF 
!                                                                       
!     - Convert to Counts                                               
!                                                                       
   ioff = offxy (ig - 1) 
   DO k = 1, lenc(ig) 
      delt = REAL(tof (k + 1) - tof (k) ) 
      gdat (k) = y (ioff + k) * delt 
      gsig (k) = dy (ioff + k) * delt 
   ENDDO 
!                                                                       
!     - Write data here                                                 
!                                                                       
   nrec = (lenc(ig) - 1) / 5 + 1 
   WRITE (outstr, 2200) 'BANK ', ig, lenc(ig) , nrec, ' TIME_MAP', imap, ' ESD'
   WRITE (ifil, 1000, rec = irec) outstr (1:80) //eol 
   irec = irec + 1 
!                                                                       
   ibeg = 1 
   DO j = 1, nrec 
      ifin = min (lenc(ig), ibeg + 4) 
      WRITE (outstr, 2300) (gdat (k), gsig (k), k = ibeg, ifin) 
      WRITE (ifil, 1000, rec = irec) outstr (1:80) //eol 
      irec = irec + 1 
      ibeg = ibeg + 5 
   ENDDO 
ENDDO 
!                                                                       
CLOSE (ifil) 
deallocate(gdat)
deallocate(gsig)
deallocate(tmap)
deallocate( tof)
RETURN 
!                                                                       
  999 CONTINUE 
      ier_num = - 2 
      ier_typ = ER_IO 
deallocate(gdat)
deallocate(gsig)
deallocate(tmap)
deallocate( tof)
!                                                                       
 1000 FORMAT    (a) 
 1040 FORMAT    ('Instrument parameter file: ',a) 
 1050 FORMAT    ('Monitor: ',i15) 
 1100 FORMAT    ('#',79('-')) 
 1110 FORMAT    ('# ',a) 
 2000 FORMAT    (a,3i5,a,i5) 
 2100 FORMAT    (10i8) 
 2200 FORMAT    (a,i3,i6,i5,a,i3,a) 
 2300 FORMAT    (10f8.1) 
 5000 FORMAT    (1x,'Saving GSAS file ',a,' ..') 
!                                                                       
END SUBROUTINE do_save_gsas                   
!
!*****7**************************************************************** 
!
SUBROUTINE make_timemap (idat, tof, tmap, nt, lnew) 
!-                                                                      
!     Creates timemap for GSAS file                                     
!     NOTE: We need to go back from musec to 100ns for this !!!         
!+                                                                      
USE errlist_mod 
USE prompt_mod 
USE kuplot_config 
USE kuplot_mod 
!                                                                       
IMPLICIT none 
!                                                                       
INTEGER, intent(in)  :: idat 
INTEGER, intent(in)  :: tof (MAXGSAS) 
INTEGER, intent(out)  :: tmap (MAXGSAS) 
integer, intent(out)  :: nt
LOGICAL, intent(out)  :: lnew 
!
INTEGER ::  dt, j 
!                                                                       
IF (idat.eq.1) then 
   lnew = .true. 
ELSE 
   IF (lenc(idat) .ne.lenc(idat - 1) ) then 
      lnew = .true. 
   ELSE 
      lnew = .false. 
      DO j = 1, lenc(idat) 
         lnew = lnew.or. (x (offxy (idat - 1) + j) .ne.x (offxy (    &
         idat - 2) + j) )                                            
      ENDDO 
   ENDIF 
ENDIF 
!                                                                       
IF (.not.lnew) return 
!                                                                       
tmap (1) = 1 
tmap (2) = tof (1) 
tmap (3) = tof (2) - tof (1) 
!                                                                       
nt = 1 
dt = tmap (3) 
!                                                                       
DO j = 2, lenc(idat) 
   IF ( (tof (j + 1) - tof (j) ) .ne.dt) then 
      nt = nt + 3 
      tmap (nt) = j 
      tmap (nt + 1) = tof (j) 
      dt = tof (j + 1) - tof (j) 
      tmap (nt + 2) = dt 
   ENDIF 
ENDDO 
!                                                                       
nt = nt + 3 
tmap (nt) = tof (lenc(idat) + 1) 
!                                                                       
END SUBROUTINE make_timemap                   
!
!*****7**************************************************************** 
!                                                                       
!     Here are the routines necessary to save data sets in              
!     various formats .. old style                                      
!                                                                       
!*****7**************************************************************** 
!
SUBROUTINE do_ksav (zei, lp) 
!                                                                       
!     Main save menu                                                    
!                                                                       
USE ber_params_mod
USE build_name_mod
USE doact_mod 
USE do_wait_mod
USE errlist_mod 
USE get_params_mod
USE learn_mod 
USE class_macro_internal
USE param_mod 
USE prompt_mod 
USE kuplot_config 
USE kuplot_mod 
use kuplot_math_mod
USE sup_mod
USE lib_errlist_func
USE lib_help
USE lib_macro_func
USE precision_mod
USE str_comp_mod
!                                                                       
IMPLICIT none 
!                                                                       
INTEGER maxw 
PARAMETER (maxw = 10) 
!                                                                       
CHARACTER(len=*), intent(inout) :: zei 
integer         , intent(inout) :: lp
!
CHARACTER(LEN=PREC_STRING) :: cpara (maxw) 
CHARACTER(LEN=PREC_STRING) :: zeile, line 
CHARACTER(LEN=PREC_STRING) :: filname 
CHARACTER(LEN=40         ) :: orig_prompt
CHARACTER(len=11) :: befehl, czeile 
CHARACTER(len=2) :: form, cdummy 
INTEGER :: lpara (maxw), lbef 
INTEGER :: i, ianz, iianz, ik, length 
LOGICAL :: l_m999 
REAL(KIND=PREC_DP) :: werte (maxw), wwerte (maxw) 
!                                                                       
!                                                                       
CALL no_error 
l_m999 = .false. 
!                                                                       
CALL get_params (zei, ianz, cpara, lpara, maxw, lp) 
IF (ier_num.ne.0) return 
!                                                                       
!------ get data set number                                             
!                                                                       
cond_ianz: IF(ianz.eq.1) then 
   CALL ber_params (1, cpara, lpara, werte, maxw) 
   IF (ier_num.ne.0) return 
   ik = nint (werte (1) ) 
   cond_ik: IF (ik.ge.1.and.ik.lt.iz) then 
      filname = 'kuplot.sav' 
      form = fform (ik) 
      IF(form.eq.'Z1'.or.form.eq.'Z2'.or.form.eq.'T1'.or.form.eq.'T2'.or.       &
         form.eq.'C1'.or.form.eq.'C2'.or.form.eq.'P1'.or.form.eq.'P2'.or.       &
         form.eq.'O1'.or.form.eq.'O2'.or.form.eq.'H1'.or.form.eq.'H2'.or.       &
         form.eq.'K1'.or.form.eq.'K2'.or.form.eq.'L1'.or.form.eq.'L2'.or.       &
         form.eq.'DM'.or.form.eq.'TE'.or.form.eq.'ST'.or.form.eq.'SC'.or.       &
         form.eq.'XX'.or.form.eq.'MC'.or.form.eq.'ZE'.or.form.eq.'YX'.or.       &
         form.eq.'GS'.or.form.eq.'SM') form = 'XY'                                 
      IF(form.eq.'ZZ'.or.form.eq.'DA'.or.form.eq.'AS'.or.form.eq.'MP') form = 'PG'
      iianz = 0 
      CALL check_form (ik, form, iianz, wwerte, maxw) 
      pgmlow = zmin (ik) 
      pgmhigh = zmax (ik) 
!
      orig_prompt = prompt
      prompt = pname//'/ksav' 
!                                                                       
!------ --- here starts the sublevel loop                               
!                                                                       
!  10       CONTINUE 
      loop_main: do
         CALL get_cmd (line, length, befehl, lbef, zeile, lp, prompt) 
         cond_ier: IF (ier_num.eq.0) then 
            IF (line.eq.' '.or.line (1:1) .eq.'#') cycle loop_main
!           IF (line.eq.' '.or.line (1:1) .eq.'#') goto 10 
!                                                                       
!------ --- subcommand 'run'                                            
!                                                                       
           IF (str_comp (befehl, 'run', 2, lbef, 3) ) then 
               CALL do_save(ik, filname, form, iianz, wwerte, maxw, l_m999)
!           GOTO 20 
               exit loop_main
!                                                                       
!     --- help 'help','?'                                               
!                                                                       
           ELSEIF(str_comp (befehl, 'help', 2, lbef, 4) .or.                    &
                  str_comp (befehl, '?   ', 1, lbef, 4) ) then
               IF (str_comp (zeile, 'errors', 2, lp, 6) ) then 
                  lp = lp + 7 
                  CALL do_hel ('kuplot '//zeile, lp) 
               ELSE 
                  lp = lp + 12 
                  CALL do_hel ('kuplot ksav '//zeile, lp) 
               ENDIF 
!                                                                       
!     --- continues a macro 'continue'                                  
!                                                                       
            ELSEIF (str_comp (befehl, 'continue', 3, lbef, 8) ) then 
               CALL macro_continue (zeile, lp) 
!                                                                       
!------ --- subcommand 'outfile'                                        
!                                                                       
            ELSEIF (str_comp (befehl, 'outfile', 2, lbef, 7) ) then 
               CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
               IF (ianz.ge.1) then 
                  CALL do_build_name (ianz, cpara, lpara, werte, maxw, 1)
                  IF (ier_num.eq.0) then 
                     filname = cpara (1) 
                  ENDIF 
               ELSE 
                  ier_num = - 6 
                  ier_typ = ER_COMM 
               ENDIF 
!                                                                       
!------ --- subcommand 'form'                                           
!                                                                       
            ELSEIF (str_comp (befehl, 'format', 2, lbef, 6) ) then 
               CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
               cdummy = cpara (1) (1:2)
               IF(cdummy=='xy') THEN
                  l_two_col = .false.
                  IF(str_comp(cpara(ianz), 'two',3, lpara(ianz), 3)) THEN
                     l_two_col = .true.
                     ianz = ianz - 1
                  ELSEIF(str_comp(cpara(ianz), 'four',3, lpara(ianz), 4)) THEN
                     l_two_col = .false.
                     ianz = ianz - 1
                  ENDIF
               ENDIF
               IF (ianz.gt.1) then 
                  CALL del_params (1, ianz, cpara, lpara, maxw) 
                  CALL ber_params (ianz, cpara, lpara, werte, maxw) 
               ELSE 
                  ianz = 0 
               ENDIF 
               CALL check_form (ik, cdummy, ianz, werte, maxw) 
               IF (ier_num.eq.0) then 
                  form = cdummy 
                  iianz = ianz 
                  DO i = 1, iianz 
                     wwerte (i) = werte (i) 
                  ENDDO 
               ENDIF 
!                                                                       
!------ --- subcommand 'm999'                                           
!                                                                       
            ELSEIF (str_comp (befehl, 'm999', 2, lbef, 4) ) then 
               CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
               IF (ianz.eq.1) then 
                  IF (str_comp(cpara(1), 'data', 2, lpara(1) , 4)) then
                     l_m999 = .true. 
                  ELSEIF(str_comp(cpara(1) , 'excl', 2, lpara (1), 4) ) then
                     l_m999 = .false. 
                  ENDIF 
               ELSE 
                  ier_num = - 6 
               ENDIF 
!                                                                       
!------ --- subcommand 'show'                                           
!                                                                       
            ELSEIF (str_comp (befehl, 'show', 2, lbef, 4) ) then 
                  WRITE (output_io, 1000) ik, filname 
                  WRITE (output_io, 1100) form 
                  IF (iianz.ne.0) write (output_io, 1200) (wwerte (i),  &
                  i = 1, iianz)                                         
                  IF (form.eq.'PG') write (output_io, 1300) pgmlow,     &
                  pgmhigh                                               
!                                                                       
!     --- Set threshold for PGM files (0..255)                          
!                                                                       
            ELSEIF (str_comp (befehl, 'threshold', 2, lbef, 9) ) then 
               CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
               IF (ier_num.eq.0) then 
                  IF (ianz.eq.2) then 
                     IF (str_comp (cpara (1) , 'high', 1, lpara (1) , 4) ) then                                       
                        CALL del_params (1, ianz, cpara, lpara, maxw) 
                        CALL ber_params (ianz, cpara, lpara, werte, maxw)                                        
                        IF (ier_num.eq.0) then 
                           pgmhigh = werte (1) * 0.01 * zmax (ik) 
                        ENDIF 
                     ELSEIF (str_comp (cpara (1) , 'low', 1, lpara ( &
                        1) , 3) ) then                                  
                        CALL del_params (1, ianz, cpara, lpara, maxw) 
                        CALL ber_params (ianz, cpara, lpara, werte, maxw)                                        
                        IF (ier_num.eq.0) then 
                           pgmlow = werte (1) * 0.01 * zmin (ik) 
                        ENDIF 
                     ELSEIF (str_comp (cpara (1) , 'sigma', 1, lpara &
                        (1) , 5) ) then                                 
                        CALL del_params (1, ianz, cpara, lpara, maxw) 
                        CALL ber_params (ianz, cpara, lpara, werte, maxw)                                        
                        IF (ier_num.eq.0) then 
                           WRITE (czeile, '(i4)') ik 
                           length = 4
                           CALL do_mean (czeile, length, .false.) 
                           pgmhigh = MIN(REAL(zmax(ik),KIND(1.D0)), res_para (3)    &
                                     + werte (1) * res_para (6) )              
                           pgmlow = MAX(REAL(zmin(ik),KIND(1.D0)), res_para (3)     &
                                     - werte (1) * res_para (6) )              
                        ENDIF 
                     ELSEIF (str_comp (cpara (1) , 'zmax', 3, lpara (&
                        1) , 4) ) then                                  
                        CALL del_params (1, ianz, cpara, lpara, maxw) 
                        CALL ber_params (ianz, cpara, lpara, werte, maxw)                                        
                           IF (ier_num.eq.0) then 
                              pgmhigh = werte (1) 
                           ENDIF 
                     ELSEIF (str_comp (cpara (1) , 'zmin', 3, lpara( 1) , 4) ) then                                  
                        CALL del_params (1, ianz, cpara, lpara, maxw) 
                        CALL ber_params (ianz, cpara, lpara, werte, maxw)                                        
                        IF (ier_num.eq.0) then 
                           pgmlow = werte (1) 
                        ENDIF 
                     ELSE 
                        ier_num = - 6 
                        ier_typ = ER_COMM 
                     ENDIF 
                  ELSE 
                     ier_num = - 6 
                     ier_typ = ER_COMM 
                  ENDIF 
               ENDIF 
!                                                                       
!------ --- subcommand 'exit'                                           
!                                                                    
            ELSEIF (str_comp (befehl, 'exit', 2, lbef, 4) ) then 
!                 GOTO 20 
                  exit loop_main
!                                                                       
!------ --- subcommand 'wait'                                           
!                                                                       
            ELSEIF (str_comp (befehl, 'wait', 2, lbef, 4) ) then 
               CALL do_input (zeile, lp) 
!                                                                       
!------ --- no command found                                            
!                                                                       
            ELSE 
               ier_num = - 8 
               ier_typ = ER_COMM 
            ENDIF 
         ENDIF cond_ier
!                                                                       
         IF (ier_num.ne.0) THEN 
            CALL errlist 
            IF (ier_sta.ne.ER_S_LIVE) THEN 
               IF (lmakro .OR. lmakro_error) THEN 
                  IF(sprompt /= prompt .OR. lmakro_error) THEN
                     ier_num = -10
                     ier_typ = ER_COMM
                     ier_msg(1) = ' Error occured in save menu'
                     prompt_status = PROMPT_ON 
                     prompt = orig_prompt
                     RETURN 
                  ELSE
                     IF(lmacro_close) THEN
                        CALL macro_close 
                        prompt_status = PROMPT_ON 
                     ENDIF 
                  ENDIF 
               ENDIF 
               IF (lblock) THEN 
                  ier_num = - 11 
                  ier_typ = ER_COMM 
                  prompt_status = PROMPT_ON 
                  prompt = orig_prompt
                  RETURN 
               ENDIF 
               CALL no_error 
               lmakro_error = .FALSE.
               sprompt = ' '
            ENDIF 
         ENDIF
!           GOTO 10 
      enddo loop_main
!                                                                       
!------ --- here ends the sublevel loop                                 
!                                                                       
!  20 CONTINUE 
      prompt = orig_prompt
   ELSE  cond_ik
      ier_num = - 4 
      ier_typ = ER_APPL 
   ENDIF  cond_ik
ELSE  cond_ianz
   ier_num = - 6 
   ier_typ = ER_COMM 
ENDIF  cond_ianz
!                                                                       
 1000 FORMAT    (' Current KSAV settings (data set ',i3,') :',//,       &
     &                  3x,'Output filename   : ',a40)                  
 1100 FORMAT    (3x,'Format            : ',a2) 
 1200 FORMAT    (3x,'Parameters        : ',5(g12.5,2x)) 
 1300 FORMAT    (3x,'Z range -> 0..255 : ',2(g12.5,2x)) 
!
END SUBROUTINE do_ksav                        
!
!*****7************************************************************     
!
SUBROUTINE check_form (ikk, form, ianz, werte, maxw) 
!+                                                                      
!     This routine checks if file format and parameters                 
!     are valid.                                                        
!-                                                                      
USE errlist_mod 
USE kuplot_config 
USE kuplot_mod 
use kuplot_para_mod
!
USE precision_mod
USE string_convert_mod
!                                                                       
IMPLICIT none 
!                                                                       
integer         , intent(in)    :: ikk          ! ksav ikk
CHARACTER(len=*), intent(inout) :: form 
integer         , intent(inout) :: ianz
integer         , intent(in) :: MAXW
REAL(KIND=PREC_DP), dimension(MAXW), intent(out) :: werte !(MAXW) 
!
INTEGER :: ik
!                                                                       
CALL do_cap (form) 
CALL skalieren 
!                                                                       
!------ xy files : parameters: xmin, xmax                               
!                                                                       
IF (form.eq.'XY'.or.form.eq.'DY') then 
   IF (ianz.eq.0) then 
      werte (1) = ex (iwin, iframe, 1) 
      werte (2) = ex (iwin, iframe, 2) 
      ianz = 2 
   ELSEIF (ianz.eq.2) then 
      CONTINUE 
   ELSE 
      ier_num = - 6 
      ier_typ = ER_COMM 
   ENDIF 
elseif (form=='H5') then
   if (ianz.eq.0) then 
      werte (1) = ex (iwin, iframe, 1) 
      werte (2) = ex (iwin, iframe, 2) 
      werte (3) = ey (iwin, iframe, 1) 
      werte (4) = ey (iwin, iframe, 2) 
      ianz = 4 
   elseif(ianz==2 .or. ianz==4) then
      continue
   endif
!                                                                       
!------ gn,ni,pg files : parameters: xmin, xmax, ymin, ymax             
!                                                                       
ELSEIF (form.eq.'GN'.or.form.eq.'NI'.or.form.eq.'PG') then 
   IF (ianz.eq.0) then 
      werte (1) = ex (iwin, iframe, 1) 
      werte (2) = ex (iwin, iframe, 2) 
      werte (3) = ey (iwin, iframe, 1) 
      werte (4) = ey (iwin, iframe, 2) 
      ianz = 4 
   ELSEIF (ianz.eq.4) then 
      CONTINUE 
   ELSE 
      ier_num = - 6 
      ier_typ = ER_COMM 
   ENDIF 
!                                                                       
!------ sx, sy files : parameter: cutting at x/y                        
!                                                                       
ELSEIF (form.eq.'SX'.or.form.eq.'SY' .or. form=='SZ') then 
   if(lh5(ikk)) then
      if (ianz /= 2) then 
         ier_num = -6 
         ier_typ = ER_COMM 
      endif 
   else
      if (ianz.ne.1) then 
         ier_num = -6 
         ier_typ = ER_COMM 
      endif 
   endif
!                                                                       
!------ mx, my files : parameter: cutting through max no.               
!                                                                       
ELSEIF (form.eq.'MX'.or.form.eq.'MY') then 
   IF (ianz.ne.1) then 
      ier_num = - 6 
      ier_typ = ER_COMM 
   ENDIF 
!                                                                       
!------ sl files : parameter: cutting from x1,y1 to x2,y2, n points     
!                                                                       
ELSEIF (form.eq.'SL') then 
   IF (ianz.ne.5) then 
      ier_num = - 6 
      ier_typ = ER_COMM 
   ENDIF 
!                                                                       
!------ sk files : parameter: cut along xy of data set ik               
!                                                                       
ELSEIF (form.eq.'SK') then 
   IF (ianz.eq.1) then 
      ik = nint (werte (1) ) 
      IF (ik.lt.1.or.ik.gt. (iz - 1) ) then 
         ier_num = - 4 
         ier_typ = ER_APPL 
      ENDIF 
   ELSE 
      ier_num = - 6 
      ier_typ = ER_COMM 
   ENDIF 
!                                                                       
ELSE 
   ier_num = - 7 
   ier_typ = ER_APPL 
ENDIF 
!                                                                       
END SUBROUTINE check_form                     
!
!*****7************************************************************     
!
SUBROUTINE do_save (ik, filname, form, ianz, werte, maxw, l_m999) 
!+                                                                      
!     Here is the main save routine                                     
!-                                                                      
USE errlist_mod 
USE kuplot_config 
USE kuplot_mod 
use kuplot_extrema_mod
use kuplot_math_mod
!
use gen_hdf_write_mod
USE lib_length
USE support_mod
use precision_mod
!                                                                       
IMPLICIT none 
!                                                                       
integer                            , intent(in) :: ik
CHARACTER(len=*)                   , intent(in) :: filname
CHARACTER(len=*)                   , intent(in) :: form 
integer                            , intent(in) :: ianz
integer                            , intent(in) :: MAXW
REAL(KIND=PREC_DP), dimension(MAXW), intent(in) :: werte (maxw)
LOGICAL                            , intent(in) :: l_m999 
!
INTEGER isa, maxmax, maxf 
PARAMETER (isa = 66, maxmax = 50, maxf = 3) 
!                                                                       
REAL :: xf (maxf), yf (maxf), zf (maxf, maxf) 
REAL(KIND=PREC_DP) :: wmax (maxmax) 
REAL :: rdx, rdy 
REAL :: xsteig, xabsch, xanf, xend, xdel 
REAL :: xxx, yyy, zzz, dzzz 
REAL :: wx_min, wx_max, wy_min, wy_max 
real(kind=PREC_DP) :: TOL    ! Tolerance at xmin,xmax... set to 0.1*x(2)-x(1)
logical                              :: hh5_laver
real(kind=PREC_DP)                   :: hh5_valmax
integer                              :: hh5_value
integer                              :: HH5_VAL_PDF
integer                              :: HH5_VAL_3DPDF
INTEGER              , DIMENSION(3)  :: hh5_out_inc
REAL(kind=PREC_DP)   , DIMENSION(3,4):: hh5_out_eck ! (3,4)
REAL(kind=PREC_DP)   , DIMENSION(3,3):: hh5_out_vi
REAL(kind=PREC_DP)   , DIMENSION(3)  :: hh5_cr_a0
REAL(kind=PREC_DP)   , DIMENSION(3)  :: hh5_cr_win
INTEGER :: ixm (maxmax), iym (maxmax) 
INTEGER, dimension(:), allocatable :: ipg ! (maxarray) 
INTEGER :: i, ix, iy, ie, k 
integer :: iix, iiy
INTEGER :: ninterv  ! number of (points-1)==no of intervals writen to file
INTEGER :: ispk, ixxx, ikk, ima, nma 
INTEGER :: nx_min, nx_max, ny_min, ny_max, nx_s, ny_s 
real(kind=PREC_DP), dimension(:,:,:), allocatable :: qvalues
!                                                                       
if(lh5(ik)) then
   call do_save_global(ik, filname, form, MAXW, werte)
   return
endif
!                                                                       
if(form(1:2) /= 'H5') CALL oeffne (isa, filname, 'unknown') 
IF (ier_num.ne.0) return 
!
! Tolerance as 1/10 of step width
TOL = abs(x(offxy(ik - 1) + 2) - x(offxy(ik - 1) + 1))*0.1D0
!                                                                       
!-------xy-kurve abspeichern                                            
!                                                                       
!write(*,*) ' KSAV ', form(1:2), lni(ik), lh5(ik), epsilon(x), epsilon(werte)
!write(*,*) ' LENC ', lenc(ik), werte(1), werte(2)
!i = 1
!write(*,*) ' x(1) ', x(offxy(ik - 1) + i), x(offxy(ik - 1) + i)-werte(1), x(offxy(ik - 1) + i) .ge. werte(1), &
!x(offxy(ik - 1) + i)-werte(1)>-TOL
!i = lenc(ik)
!write(*,*) ' x(1) ', x(offxy(ik - 1) + i), werte(2)-x(offxy(ik - 1) + i), x(offxy(ik - 1) + i) .le. werte(2), &
!x(offxy(ik - 1) + i)-werte(2)<TOL
IF (form (1:2) .eq.'XY'.and..not.lni (ik) ) then 
   WRITE (isa, 3000) titel (iwin, iframe, 1) (1:len_str (titel (  &
         iwin, iframe, 1) ) )                                           
   WRITE (isa, 3000) titel (iwin, iframe, 2) (1:len_str (titel (  &
         iwin, iframe, 2) ) )                                           
   IF(l_two_col) THEN
      DO i = 1, lenc(ik) 
!        IF(x(offxy(ik - 1) + i) .ge. werte(1) .and. x(offxy(ik - 1) + i) .le.werte(2) ) then
         if(x(offxy(ik - 1) + i)-werte(1)>-TOL .and. x(offxy(ik - 1) + i)-werte(2)<TOL) then
            WRITE(isa, 4000) x(offxy(ik - 1) + i), y(offxy(ik - 1) + i)
         ENDIF 
      ENDDO 
   ELSE
      DO i = 1, lenc(ik) 
!        IF(x(offxy(ik - 1) + i) .ge. werte(1) .and. x(offxy(ik - 1) + i) .le. werte(2)) then
         if(x(offxy(ik - 1) + i)-werte(1)>-TOL .and. x(offxy(ik - 1) + i)-werte(2)<TOL) then
            WRITE(isa, 4000) x (offxy(ik - 1) + i), y (offxy(ik - 1) + i), & 
                             dx(offxy(ik - 1) + i), dy(offxy(ik - 1) + i)      
         ENDIF 
      ENDDO 
   ENDIF 
!                                                                       
ELSEIF (form (1:2) .eq.'DY'.and..not.lni (ik) ) then 
   WRITE(isa, 3000) titel(iwin, iframe, 1)(1:len_str(titel(iwin, iframe, 1)))
   WRITE(isa, 3000) titel(iwin, iframe, 2)(1:len_str(titel(iwin, iframe, 2)))
   DO i = 1, lenc(ik) 
!     IF(x(offxy(ik - 1) + i) .ge. werte(1) .and. x(offxy(ik - 1) + i) .le. werte(2)) then
      if(x(offxy(ik - 1) + i)-werte(1)>-TOL .and. x(offxy(ik - 1) + i)-werte(2)<TOL) then
         WRITE(isa, 4000) x(offxy (ik - 1) + i), y(offxy(ik - 1) + i),          &
                          dy(offxy (ik - 1) + i)
      ENDIF 
   ENDDO 
!                                                                       
!-------nipl-file fuer gnuplot abspeichern                              
!                                                                       
ELSEIF (form (1:2) .eq.'GN'.and.lni (ik) ) then 
   rdx = (xmax (ik) - xmin (ik) ) / REAL(nx (ik) - 1) 
   rdy = (ymax (ik) - ymin (ik) ) / REAL(ny (ik) - 1) 
   nx_min = max (1, nint ( ( (werte (1) - xmin (ik) ) / rdx) ) + 1)
   nx_max = min (nx (ik), nint ( ( (werte (2) - xmin (ik) ) / rdx) ) + 1) 
   ny_min = max (1, nint ( ( (werte (3) - ymin (ik) ) / rdy) ) + 1)
   ny_max = min (ny (ik), nint ( ( (werte (4) - ymin (ik) )  / rdy) ) + 1)
   DO ix = nx_min, nx_max 
      DO iy = ny_min, ny_max 
         WRITE (isa, 4000) x (offxy (ik - 1) + ix), y (offxy (ik - 1)   &
         + iy), z (offz (ik - 1) + (ix - 1) * ny (ik) + iy)             
      ENDDO 
      WRITE (isa, 2000) ' ' 
   ENDDO 
!                                                                       
!-------nipl-file-schnitt abspeichern (beliebige richtung)              
!                                                                       
ELSEIF (form (1:2) .eq.'SL'.and.lni (ik) ) then 
   xsteig = (werte (4) - werte (2) ) / (werte (3) - werte (1) ) 
   xabsch = werte (2) - xsteig * werte (1) 
   xanf = ex (iwin, iframe, 1) 
   xend = ex (iwin, iframe, 2) 
   xdel = (werte (3) - werte (1) ) / werte (5) 
   ninterv = NINT((xend-xanf)/xdel)
!         DO xxx = xanf, xend, xdel
   DO i = 0, ninterv
      xxx = xanf + i* xdel
      yyy = xsteig * xxx + xabsch 
      CALL extract_subarray (xf, yf, zf, xxx, yyy, maxf, ik, ie) 
      IF (ie.eq.0) then 
         CALL polin2 (xf, yf, zf, maxf, maxf, xxx, yyy, zzz, dzzz, ier_num)
         IF(ier_num /= 0) THEN
            ier_typ = ER_APPL
            RETURN
         ENDIF 
         WRITE (isa, 4000) xxx, zzz, 0.0, abs (dzzz), yyy 
      ENDIF 
   ENDDO 
!                                                                       
!-------nipl-file-schnitt abspeichern (spur aus kurve ik)               
!                                                                       
ELSEIF (form (1:2) .eq.'SK'.and.lni (ik) ) then 
   ispk = nint (werte (1) ) 
   DO ixxx = 1, lenc(ispk) 
      xxx = x (offxy (ispk - 1) + ixxx) 
      yyy = y (offxy (ispk - 1) + ixxx) 
      CALL extract_subarray (xf, yf, zf, xxx, yyy, maxf, ik, ie) 
      IF (ie.eq.0) then 
         CALL polin2 (xf, yf, zf, maxf, maxf, xxx, yyy, zzz, dzzz,ier_num) 
         IF(ier_num /= 0) THEN
            ier_typ = ER_APPL
            RETURN
         ENDIF 
         WRITE(isa, 4100) REAL(ixxx), zzz, 0.0, abs (dzzz), xxx, yyy
      ENDIF 
   ENDDO 
!                                                                       
!-------nipl-file-schnitt abspeichern                                   
!                                                                       
ELSEIF (form (1:2) .eq.'SX'.and.lni (ik) ) then 
   rdy = (ymax (ik) - ymin (ik) ) / REAL(ny (ik) - 1) 
   ny_s = nint ( ( (werte (1) - ymin (ik) ) / rdy) ) + 1 
   DO i = 1, nx (ik) 
      ikk = offz (ik - 1) + (i - 1) * ny (ik) + ny_s 
      IF(z(ikk) .ne. - 9999.0 .or. z(ikk) .eq. -9999.0 .and. l_m999) then
         WRITE (isa, 4000) x (offxy (ik - 1) + i), z (ikk) 
      ENDIF 
   ENDDO 
!                                                                       
!-------nipl-file-schnitt durch i-tes maximum abspeichern               
!                                                                       
ELSEIF (form (1:2) .eq.'MX'.and.lni (ik) ) then 
   CALL do_fmax_z (ik, wmax, ixm, iym, maxmax, ima) 
   nma = nint (werte (1) ) 
   IF (ianz.eq.0) nma = 1 
   IF (nma.gt.ima) then 
      ier_num = - 19 
      ier_typ = ER_APPL 
      RETURN 
   ENDIF 
   rdy = (ymax (ik) - ymin (ik) ) / REAL(ny (ik) - 1) 
   ny_s = nint ( ( (y (iym (nma) ) - ymin (ik) ) / rdy) ) + 1 
   DO i = 1, nx (ik) 
      ikk = offz (ik - 1) + (i - 1) * ny (ik) + ny_s 
      IF (z (ikk) .ne. - 9999.0.or.z (ikk) .eq. - 9999.0.and.l_m999) then                                                           
         WRITE (isa, 4000) x (offxy (ik - 1) + i), z (ikk) 
      ENDIF 
   ENDDO 
!                                                                       
!-------nipl-file-schnitt abspeichern                                   
!                                                                       
ELSEIF (form (1:2) .eq.'SY'.and.lni (ik) ) then 
   rdx = (xmax (ik) - xmin (ik) ) / REAL(nx (ik) - 1) 
   nx_s = nint ( ( (werte (1) - xmin (ik) ) / rdx) ) + 1 
   DO i = 1, ny (ik) 
      ikk = offz (ik - 1) + (nx_s - 1) * ny (ik) + i 
      IF (z (ikk) .ne. - 9999.0.or.z (ikk) .eq. - 9999.0.and.l_m999) then
         WRITE (isa, 4000) y (offxy (ik - 1) + i), z (ikk) 
      ENDIF 
   ENDDO 
!                                                                       
!-------nipl-file-schnitt abspeichern durch maximum                     
!                                                                       
ELSEIF (form (1:2) .eq.'MY'.and.lni (ik) ) then 
   CALL do_fmax_z (ik, wmax, ixm, iym, maxmax, ima) 
   nma = nint (werte (1) ) 
   IF (ianz.eq.0) nma = 1 
   IF (nma.gt.ima) then 
      ier_num = - 19 
      ier_typ = ER_APPL 
      RETURN 
   ENDIF 
   rdx = (xmax (ik) - xmin (ik) ) / REAL(nx (ik) - 1) 
   nx_s = nint ( ( (x (ixm (nma) ) - xmin (ik) ) / rdx) ) + 1 
   DO i = 1, ny (ik) 
      ikk = offz (ik - 1) + (nx_s - 1) * ny (ik) + i 
      IF (z (ikk) .ne. - 9999.0.or.z (ikk) .eq. - 9999.0.and.l_m999) then
         WRITE (isa, 4000) y (offxy (ik - 1) + i), z (ikk) 
      ENDIF 
   ENDDO 
!                                                                       
!-------pgm-file abspeichern (aktueller ausschnitt)                     
!                                                                       
ELSEIF (form (1:2) .eq.'PG'.and.lni (ik) ) then 
   rdx = (xmax (ik) - xmin (ik) ) / REAL(nx (ik) - 1) 
   rdy = (ymax (ik) - ymin (ik) ) / REAL(ny (ik) - 1) 
   nx_min = max (1, nint ( ( (werte (1) - xmin (ik) ) / rdx) ) + 1)
   nx_max = min (nx (ik), nint ( ( (werte (2) - xmin (ik) )  / rdx) ) + 1)
   ny_min = max (1, nint ( ( (werte (3) - ymin (ik) ) / rdy) ) + 1)
   ny_max = min (ny (ik), nint ( ( (werte (4) - ymin (ik) )  / rdy) ) + 1)
   WRITE (isa, 2000) 'P2' 
   WRITE (isa, 2000) '# PGM-file created by KUPLOT' 
   WRITE (isa, * ) nx_max - nx_min + 1, ny_max - ny_min + 1 
   WRITE (isa, * ) 255 
   allocate(ipg(maxarray))
   DO iy = ny_max, ny_min, - 1 
      DO ix = nx_min, nx_max 
         k = offz (ik - 1) + (ix - 1) * ny (ik) + iy 
         IF (z (k) .lt.pgmlow) then 
            ipg (ix - nx_min + 1) = 0 
         ELSEIF (z (k) .gt.pgmhigh) then 
            ipg (ix - nx_min + 1) = 255 
         ELSEIF (z (k) .eq. - 9999.0) then 
            ipg (ix - nx_min + 1) = 0 
         ELSE 
            ipg (ix - nx_min + 1) = nint (255. * (z (k) - pgmlow)       &
            / (pgmhigh - pgmlow) )                                      
         ENDIF 
      ENDDO 
      WRITE (isa, 4500) (ipg (ix), ix = 1, nx_max - nx_min + 1) 
   ENDDO 
   deallocate(ipg)
!                                                                       
!-------nipl-file abspeichern (aktueller ausschnitt)                    
!                                                                       
ELSEIF (form (1:2) .eq.'NI'.and.lni (ik) ) then 
   rdx = (xmax (ik) - xmin (ik) ) / REAL(nx (ik) - 1) 
   rdy = (ymax (ik) - ymin (ik) ) / REAL(ny (ik) - 1) 
   nx_min = max (1, nint ( ( (werte (1) - xmin (ik) ) / rdx) ) + 1)
   nx_max = min (nx (ik), nint ( ( (werte (2) - xmin (ik) ) / rdx) ) + 1)
   ny_min = max (1, nint ( ( (werte (3) - ymin (ik) ) / rdy) ) + 1)
   ny_max = min (ny (ik), nint ( ( (werte (4) - ymin (ik) ) / rdy) ) + 1)
   wx_min = xmin (ik) + (nx_min - 1) * rdx 
   wx_max = xmin (ik) + (nx_max - 1) * rdx 
   wy_min = ymin (ik) + (ny_min - 1) * rdy 
   wy_max = ymin (ik) + (ny_max - 1) * rdy 
   WRITE (isa, * ) nx_max - nx_min + 1, ny_max - ny_min + 1 
   WRITE (isa, 4000) wx_min, wx_max, wy_min, wy_max 
   DO iy = ny_min, ny_max 
      WRITE (isa, 4000) (z (offz (ik - 1) + (ix - 1) * ny (ik) + iy), ix = nx_min, nx_max)
      ENDDO 
elseif(form(1:2)=='H5') then
   if(lh5(ik)) then
      call do_save_global(ik, filname, form, MAXW, werte)
   else
   hh5_value     = 1              ! regular data
   HH5_VAL_PDF   = 1              ! regular data
   HH5_VAL_3DPDF = 1              ! regular data
   hh5_laver     = .false.        ! No averaged intensities
  
   if(allocated(qvalues)) deallocate(qvalues)
   if(lni(ik)) then               ! 2D- Nipl file
!
      rdx = (xmax (ik) - xmin (ik) ) / REAL(nx (ik) - 1) 
      rdy = (ymax (ik) - ymin (ik) ) / REAL(ny (ik) - 1) 
      nx_min = max (1, nint ( ( (werte (1) - xmin (ik) ) / rdx) ) + 1)
      nx_max = min (nx (ik), nint ( ( (werte (2) - xmin (ik) ) / rdx) ) + 1)
      ny_min = max (1, nint ( ( (werte (3) - ymin (ik) ) / rdy) ) + 1)
      ny_max = min (ny (ik), nint ( ( (werte (4) - ymin (ik) ) / rdy) ) + 1)
      wx_min = xmin (ik) + (nx_min - 1) * rdx 
      wx_max = xmin (ik) + (nx_max - 1) * rdx 
      wy_min = ymin (ik) + (ny_min - 1) * rdy 
      wy_max = ymin (ik) + (ny_max - 1) * rdy 
!
      hh5_out_inc(1) = nx_max - nx_min + 1
      hh5_out_inc(2) = ny_max - ny_min + 1
      hh5_out_inc(3) = 1
      hh5_out_eck      = 0.0D0
      hh5_out_eck(1,1) = wx_min
      hh5_out_eck(2,1) = wy_min
!
      hh5_out_vi       = 0.0D0
      hh5_out_vi (1,1) = rdx
      hh5_out_vi (2,2) = rdy
!
      hh5_cr_a0        =  1.0D0    ! Make cartesian space
      hh5_cr_win       = 90.0D0
      hh5_valmax       =  0.0D0
!
      allocate(qvalues(hh5_out_inc(1), hh5_out_inc(2), hh5_out_inc(3)))
      iiy = 0
      do iy = ny_min, ny_max
        iiy = iiy + 1
        iix = 0
        do ix = nx_min, nx_max
           iix = iix + 1
          qvalues(iix, iiy,1) = (z (offz (ik - 1) + (ix - 1) * ny (ik) + iy))
         enddo
      enddo
   else                               ! 1D-      file
!
      rdx = (xmax (ik) - xmin (ik) ) / REAL(lenc (ik) - 1) 
      nx_min = max (1, nint ( ( (werte (1) - xmin (ik) ) / rdx) ) + 1)
      nx_max = min (lenc(ik), nint ( ( (werte (2) - xmin (ik) ) / rdx) ) + 1)
      wx_min = xmin (ik) + (nx_min - 1) * rdx 
      wx_max = xmin (ik) + (nx_max - 1) * rdx 
!
      hh5_out_inc(1) = nx_max - nx_min + 1
      hh5_out_inc(2) = 1
      hh5_out_inc(3) = 1
      hh5_out_eck      = 0.0D0
      hh5_out_eck(1,1) = wx_min
!
      hh5_out_vi       = 0.0D0
      hh5_out_vi (1,1) = rdx
!
      hh5_cr_a0        =  1.0D0    ! Make cartesian space
      hh5_cr_win       = 90.0D0
      hh5_valmax       =  0.0D0
!
      allocate(qvalues(hh5_out_inc(1), hh5_out_inc(2), hh5_out_inc(3)))
      iix = 0
      do ix = nx_min, nx_max
         iix = iix + 1
         qvalues(iix,  1, 1) = y(offxy(ik - 1) + ix)
      enddo
   endif
   call gen_hdf5_write (hh5_value, hh5_laver, filname, hh5_out_inc, hh5_out_eck, hh5_out_vi, &
                        hh5_cr_a0, hh5_cr_win, qvalues, HH5_VAL_PDF, HH5_VAL_3DPDF, hh5_valmax, &
                        ier_num, ier_typ, ER_IO, ER_APPL)
   deallocate(qvalues)
   endif
!
ENDIF 
if(form(1:2) /= 'H5') CLOSE (isa) 
!                                                                       
 2000 FORMAT     (a) 
 3000 FORMAT     ('#',a) 
 4000 FORMAT     (5(g15.8,2x)) 
 4100 FORMAT     (6(g13.6,1x)) 
 4500 FORMAT     (20(i3,1x)) 
!
END SUBROUTINE do_save                        
!
!*******************************************************************************
!
subroutine do_save_global(ik, filname, form, MAXW, werte)
!-
!  Save a data set from the global storage into an HDF5 File
!+
!
!use kuplot_global
!
use lib_data_struc_h5
use gen_hdf_write_mod
use errlist_mod
USE support_mod       , only:oeffne
use precision_mod
!
implicit none
!
integer         , intent(in) :: ik        ! Data set number
character(len=*), intent(in) :: filname   ! Output file name
character(len=*), intent(in) :: form      ! Output file format
integer         , intent(in) :: MAXW      ! Werte dimension
real(kind=PREC_DP), dimension(MAXW), intent(in) :: werte
!
integer, PARAMETER :: IWR = 43
logical                              :: hh5_laver
real(kind=PREC_DP)                   :: hh5_valmax
integer                              :: hh5_value
integer                              :: HH5_VAL_PDF
integer                              :: HH5_VAL_3DPDF
real(kind=PREC_DP), dimension(:)    , allocatable :: c_x
real(kind=PREC_DP), dimension(:)    , allocatable :: c_y
real(kind=PREC_DP), dimension(:)    , allocatable :: c_z
real(kind=PREC_DP), dimension(:)    , allocatable :: c_dx
real(kind=PREC_DP), dimension(:)    , allocatable :: c_dy
real(kind=PREC_DP), dimension(:)    , allocatable :: c_dz
real(kind=PREC_DP), dimension(:,:,:), allocatable :: qvalues
real(kind=PREC_DP), dimension(:,:,:), allocatable :: sigma
!
integer :: i,j,k ! Dummy loop indices
!
character(len=PREC_STRING)                         :: infile
integer                                            :: node_number  ! Node in global data
integer                                            :: nlayer       ! Current layer (3-D only
logical                                            :: is_direct    ! date in direct / reciprocal space
integer                                            :: ndims        ! Number of dimensions
integer, dimension(3)                              :: dims         ! Dimensions global array
logical                                            :: is_grid      ! date in direct / reciprocal space
logical                                            :: has_dxyz     ! date in direct / reciprocal space
logical                                            :: has_dval     ! date in direct / reciprocal space
REAL(kind=PREC_DP)   , DIMENSION(3,4)              :: corners      ! steps along each axis
REAL(kind=PREC_DP)   , DIMENSION(3,3)              :: vectors      ! steps along each axis
REAL(kind=PREC_DP)   , DIMENSION(3)                :: a0           ! Lower limits
REAL(kind=PREC_DP)   , DIMENSION(3)                :: win          ! Lower limits
REAL(kind=PREC_DP)   , DIMENSION(3)                :: llims        ! Lower limits
REAL(kind=PREC_DP)   , DIMENSION(3)                :: steps        ! steps along each axis
REAL(kind=PREC_DP)   , DIMENSION(3,3)              :: steps_full   ! steps along each axis
REAL(kind=PREC_DP)   , DIMENSION(2)                :: minmaxval    ! steps along each axis
REAL(kind=PREC_DP)   , DIMENSION(3,2)              :: minmaxcoor   ! steps along each axis
!
call data2local(ik, ier_num, ier_typ, node_number, infile, nlayer, is_direct,   &
                ndims, dims, is_grid, has_dxyz, has_dval, corners, vectors,     &
                a0, win, c_x, c_y, c_z, c_dx, c_dy,c_dz, qvalues, sigma, llims, &
                steps, steps_full, minmaxval, minmaxcoor)
!
if(form/='H5') then
   call oeffne(IWR, filname, 'unknown')
   if(ier_num/=0) then
      return
   endif
endif
!
if(form=='H5')then
!
   if(is_direct) then
      hh5_value     = 1              ! Direct space data
      HH5_VAL_PDF   = 1              ! Direct space data
      HH5_VAL_3DPDF = 1              ! Direct space data
      hh5_laver     = .false.        ! No averaged intensities
   else
      hh5_value     = 0              ! reciprocal data
      HH5_VAL_PDF   = 0              ! reciprocal data
      HH5_VAL_3DPDF = 0              ! reciprocal data
      hh5_laver     = .false.        ! No averaged intensities
   endif
!
   hh5_valmax      = 0.0D0
!
!  
   call gen_hdf5_write(hh5_value, hh5_laver, filname, dims, corners, vectors,      &
                       a0, win, qvalues, HH5_VAL_PDF, HH5_VAL_3DPDF, hh5_valmax,   &
                       ier_num, ier_typ, ER_IO, ER_APPL)
elseif(form=='SX') then
   j = 1
   k = 1
   if(dims(2)>1) then
      j = nint( (werte(1)-c_y(1))/steps(2)) + 1
   endif
   if(dims(3)>1) then
      k = nint( (werte(2)-c_z(1))/steps(3)) + 1
   endif
   write(IWR, '(a)') '#'
   write(IWR, '(a)') '#'
   do i=1, dims(1)
      write(IWR, 4000) c_x(i), qvalues(i,j,k), 0.0, 0.0
   enddo
elseif(form=='SY') then
   i = 1
   k = 1
   if(dims(1)>1) then
      i = nint( (werte(1)-c_x(1))/steps(1)) + 1
   endif
   if(dims(3)>1) then
      k = nint( (werte(2)-c_z(1))/steps(3)) + 1
   endif
   write(IWR, '(a)') '#'
   write(IWR, '(a)') '#'
   do j=1, dims(2)
      write(IWR, 4000) c_y(j), qvalues(i,j,k), 0.0, 0.0
   enddo
elseif(form=='SZ') then
   i = 1
   j = 1
   if(dims(1)>1) then
      i = nint( (werte(1)-c_x(1))/steps(1)) + 1
   endif
   if(dims(2)>1) then
      j = nint( (werte(2)-c_y(1))/steps(2)) + 1
   endif
write(*,*) ' SCHNITT AT ', i,j, werte(1:2)
   write(IWR, '(a)') '#'
   write(IWR, '(a)') '#'
   do k=1, dims(3)
      write(IWR, 4000) c_z(k), qvalues(i,j,k), 0.0, 0.0
   enddo
endif
!
deallocate(qvalues)
if(allocated(c_x)) deallocate( c_x)
if(allocated(c_y)) deallocate( c_y)
if(allocated(c_z)) deallocate( c_z)
if(allocated(c_dx)) deallocate( c_dx)
if(allocated(c_dy)) deallocate( c_dy)
if(allocated(c_dz)) deallocate( c_dz)
if(allocated(sigma)) deallocate( sigma)
!
if(form/='H5') then
  close(IWR)
endif
!
4000 format(4(g15.8e3,2x)) 
!
end subroutine do_save_global
!
!*******************************************************************************
!
end module kuplot_save_mod
