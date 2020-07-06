!*****7**************************************************************** 
!     This file contains subroutines saving files                       
!*****7**************************************************************** 
      SUBROUTINE do_save_data (zeile, lp) 
!-                                                                      
!       Main menu for save command.                                     
!+                                                                      
      USE build_name_mod
      USE errlist_mod 
      USE get_params_mod
      USE kuplot_config 
USE lib_help
USE precision_mod
USE str_comp_mod
!                                                                       
      IMPLICIT none 
!                                                                       
      INTEGER maxw 
      PARAMETER (maxw = 10) 
!                                                                       
      CHARACTER ( * ) zeile 
      INTEGER lp 
!                                                                       
      CHARACTER(LEN=PREC_STRING) :: cpara (maxw) 
      INTEGER lpara (maxw) 
      integer ianz 
      rEAL(KIND=PREC_DP) :: werte (maxw) 
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
         ELSEIF (str_comp (cpara (1) , 'merge', 2, lpara (1) , 5) )     &
         then                                                           
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
!*****7**************************************************************** 
      SUBROUTINE do_save_merge (filname) 
!-                                                                      
!       Saves all current data sets as 'x dataset y' file               
!+                                                                      
      USE debug_mod 
      USE errlist_mod 
      USE prompt_mod 
      USE kuplot_config 
      USE kuplot_mod 
USE support_mod
!                                                                       
      IMPLICIT none 
!                                                                       
      INTEGER ifil 
      PARAMETER (ifil = 33) 
!                                                                       
      CHARACTER ( * ) filname 
      INTEGER ikurv, ip 
!                                                                       
      LOGICAL k_in_f 
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
      END SUBROUTINE do_save_merge                  
!*****7**************************************************************** 
      SUBROUTINE do_save_xml (filname, path) 
!-                                                                      
!       Saves all current data sets as PlotML file                      
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
      CHARACTER ( * ) filname, path 
      CHARACTER(20) connect, marks 
      INTEGER ikurv, ip 
!                                                                       
      LOGICAL k_in_f 
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
            WRITE (ifil, 1400) fname (ikurv) (1:len_str (fname (ikurv) )&
            ), connect (1:len_str (connect) ), marks (1:len_str (marks) &
            )                                                           
         ENDIF 
         DO ip = 1, lenc(ikurv) 
         WRITE (ifil, 1500) x (offxy (ikurv - 1) + ip), y (offxy (ikurv &
         - 1) + ip)                                                     
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
!*****7**************************************************************** 
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
      CHARACTER ( * ) filname 
      CHARACTER ( * ) iname 
!                                                                       
      CHARACTER(200) outstr 
      CHARACTER(2) eol 
      REAL gdat (MAXGSAS) 
      REAL gsig (MAXGSAS) 
      REAL delt 
      INTEGER tmap (MAXGSAS) 
      INTEGER tof (MAXGSAS) 
      INTEGER ig, il, ixx 
      INTEGER j, k, ibeg, ifin, imap, ioff, itof 
      INTEGER ntmap, nrec, irec
      LOGICAL lfile, lnew 
!                                                                       
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
      WRITE (outstr, 2000) 'TIME_MAP', itof, ntmap, nrec, ' TIME_MAP', 1&
     &00                                                                
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
      WRITE (outstr, 2200) 'BANK ', ig, lenc(ig) , nrec, ' TIME_MAP',   &
      imap, ' ESD'                                                      
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
      RETURN 
!                                                                       
  999 CONTINUE 
      ier_num = - 2 
      ier_typ = ER_IO 
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
!*****7**************************************************************** 
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
      INTEGER tof (MAXGSAS) 
      INTEGER tmap (MAXGSAS) 
      INTEGER idat 
      INTEGER nt, dt, j 
      LOGICAL lnew 
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
!*****7**************************************************************** 
!                                                                       
!     Here are the routines necessary to save data sets in              
!     various formats .. old style                                      
!                                                                       
!*****7**************************************************************** 
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
      CHARACTER ( * ) zei 
      CHARACTER(LEN=PREC_STRING) :: cpara (maxw) 
      CHARACTER(LEN=PREC_STRING) :: zeile, line 
      CHARACTER(LEN=PREC_STRING) :: filname 
      CHARACTER(LEN=40         ) :: orig_prompt
      CHARACTER(4) befehl, czeile 
      CHARACTER(2) form, cdummy 
      INTEGER lpara (maxw), lp, lbef 
      INTEGER i, ianz, iianz, ik, length 
      LOGICAL l_m999 
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
      IF (ianz.eq.1) then 
         CALL ber_params (1, cpara, lpara, werte, maxw) 
         IF (ier_num.ne.0) return 
         ik = nint (werte (1) ) 
         IF (ik.ge.1.and.ik.lt.iz) then 
            filname = 'kuplot.sav' 
            form = fform (ik) 
      IF (form.eq.'Z1'.or.form.eq.'Z2'.or.form.eq.'T1'.or.form.eq.'T2'.o&
     &r.form.eq.'C1'.or.form.eq.'C2'.or.form.eq.'P1'.or.form.eq.'P2'.or.&
     &form.eq.'O1'.or.form.eq.'O2'.or.form.eq.'H1'.or.form.eq.'H2'.or.fo&
     &rm.eq.'K1'.or.form.eq.'K2'.or.form.eq.'L1'.or.form.eq.'L2'.or.form&
     &.eq.'DM'.or.form.eq.'TE'.or.form.eq.'ST'.or.form.eq.'SC'.or.form.e&
     &q.'XX'.or.form.eq.'MC'.or.form.eq.'ZE'.or.form.eq.'YX'.or.form.eq.&
     &'GS'.or.form.eq.'SM') form = 'XY'                                 
      IF (form.eq.'ZZ'.or.form.eq.'DA'.or.form.eq.'AS'.or.form.eq.'MP') &
     &form = 'PG'                                                       
            iianz = 0 
            CALL check_form (form, iianz, wwerte, maxw) 
            pgmlow = zmin (ik) 
            pgmhigh = zmax (ik) 
!
            orig_prompt = prompt
            prompt = pname//'/ksav' 
!                                                                       
!------ --- here starts the sublevel loop                               
!                                                                       
   10       CONTINUE 
            CALL get_cmd (line, length, befehl, lbef, zeile, lp, prompt) 
            IF (ier_num.eq.0) then 
               IF (line.eq.' '.or.line (1:1) .eq.'#') goto 10 
!                                                                       
!------ --- subcommand 'run'                                            
!                                                                       
               IF (str_comp (befehl, 'run', 2, lbef, 3) ) then 
                  CALL do_save (ik, filname, form, iianz, wwerte, maxw, &
                  l_m999)                                               
                  GOTO 20 
!                                                                       
!     --- help 'help','?'                                               
!                                                                       
      ELSEIF (str_comp (befehl, 'help', 2, lbef, 4) .or.str_comp (befehl&
     &, '?   ', 1, lbef, 4) ) then                                      
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
               ELSEIF (str_comp (befehl, 'outf', 2, lbef, 4) ) then 
                  CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
                  IF (ianz.ge.1) then 
                     CALL do_build_name (ianz, cpara, lpara, werte,     &
                     maxw, 1)                                           
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
               ELSEIF (str_comp (befehl, 'form', 2, lbef, 4) ) then 
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
                  CALL check_form (cdummy, ianz, werte, maxw) 
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
                     IF (str_comp (cpara (1) , 'data', 2, lpara (1) , 4)&
                     ) then                                             
                        l_m999 = .true. 
                     ELSEIF (str_comp (cpara (1) , 'excl', 2, lpara (1) &
                     , 4) ) then                                        
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
               ELSEIF (str_comp (befehl, 'thre', 2, lbef, 4) ) then 
                  CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
                  IF (ier_num.eq.0) then 
                     IF (ianz.eq.2) then 
                        IF (str_comp (cpara (1) , 'high', 1, lpara (1) ,&
                        4) ) then                                       
                           CALL del_params (1, ianz, cpara, lpara, maxw) 
                           CALL ber_params (ianz, cpara, lpara, werte,  &
                           maxw)                                        
                           IF (ier_num.eq.0) then 
                              pgmhigh = werte (1) * 0.01 * zmax (ik) 
                           ENDIF 
                        ELSEIF (str_comp (cpara (1) , 'low', 1, lpara ( &
                        1) , 3) ) then                                  
                           CALL del_params (1, ianz, cpara, lpara, maxw) 
                           CALL ber_params (ianz, cpara, lpara, werte,  &
                           maxw)                                        
                           IF (ier_num.eq.0) then 
                              pgmlow = werte (1) * 0.01 * zmin (ik) 
                           ENDIF 
                        ELSEIF (str_comp (cpara (1) , 'sigma', 1, lpara &
                        (1) , 5) ) then                                 
                           CALL del_params (1, ianz, cpara, lpara, maxw) 
                           CALL ber_params (ianz, cpara, lpara, werte,  &
                           maxw)                                        
                           IF (ier_num.eq.0) then 
                              WRITE (czeile, '(i4)') ik 
                              CALL do_mean (czeile, 4, .false.) 
                              pgmhigh = MIN(REAL(zmax(ik),KIND(1.D0)), res_para (3)    &
                              + werte (1) * res_para (6) )              
                              pgmlow = MAX(REAL(zmin(ik),KIND(1.D0)), res_para (3)     &
                              - werte (1) * res_para (6) )              
                           ENDIF 
                        ELSEIF (str_comp (cpara (1) , 'zmax', 3, lpara (&
                        1) , 4) ) then                                  
                           CALL del_params (1, ianz, cpara, lpara, maxw) 
                           CALL ber_params (ianz, cpara, lpara, werte,  &
                           maxw)                                        
                           IF (ier_num.eq.0) then 
                              pgmhigh = werte (1) 
                           ENDIF 
                        ELSEIF (str_comp (cpara (1) , 'zmin', 3, lpara (&
                        1) , 4) ) then                                  
                           CALL del_params (1, ianz, cpara, lpara, maxw) 
                           CALL ber_params (ianz, cpara, lpara, werte,  &
                           maxw)                                        
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
                  GOTO 20 
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
            ENDIF 
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
            GOTO 10 
!                                                                       
!------ --- here ends the sublevel loop                                 
!                                                                       
   20       CONTINUE 
            prompt = orig_prompt
         ELSE 
            ier_num = - 4 
            ier_typ = ER_APPL 
         ENDIF 
      ELSE 
         ier_num = - 6 
         ier_typ = ER_COMM 
      ENDIF 
!                                                                       
 1000 FORMAT    (' Current KSAV settings (data set ',i3,') :',//,       &
     &                  3x,'Output filename   : ',a40)                  
 1100 FORMAT    (3x,'Format            : ',a2) 
 1200 FORMAT    (3x,'Parameters        : ',5(g12.5,2x)) 
 1300 FORMAT    (3x,'Z range -> 0..255 : ',2(g12.5,2x)) 
      END SUBROUTINE do_ksav                        
!*****7************************************************************     
      SUBROUTINE check_form (form, ianz, werte, maxw) 
!+                                                                      
!     This routine checks if file format and parameters                 
!     are valid.                                                        
!-                                                                      
      USE errlist_mod 
      USE kuplot_config 
      USE kuplot_mod 
USE precision_mod
      USE string_convert_mod
!                                                                       
      IMPLICIT none 
!                                                                       
      CHARACTER ( * ) form 
      INTEGER ik, ianz, maxw 
      REAL(KIND=PREC_DP) :: werte (maxw) 
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
      ELSEIF (form.eq.'SX'.or.form.eq.'SY') then 
         IF (ianz.ne.1) then 
            ier_num = - 6 
            ier_typ = ER_COMM 
         ENDIF 
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
!*****7************************************************************     
      SUBROUTINE do_save (ik, filname, form, ianz, werte, maxw, l_m999) 
!+                                                                      
!     Here is the main save routine                                     
!-                                                                      
      USE errlist_mod 
      USE kuplot_config 
      USE kuplot_mod 
USE lib_length
USE precision_mod
USE support_mod
!                                                                       
      IMPLICIT none 
!                                                                       
      INTEGER isa, maxw, maxmax, maxf 
      PARAMETER (isa = 66, maxmax = 50, maxf = 3) 
!                                                                       
      CHARACTER ( * ) filname, form 
      REAL xf (maxf), yf (maxf), zf (maxf, maxf) 
      REAL(KIND=PREC_DP) :: werte (maxw), wmax (maxmax) 
      REAL rdx, rdy 
      REAL xsteig, xabsch, xanf, xend, xdel 
      REAL xxx, yyy, zzz, dzzz 
      REAL wx_min, wx_max, wy_min, wy_max 
      INTEGER ixm (maxmax), iym (maxmax) 
      INTEGER ipg (maxarray) 
      INTEGER ianz, ik, i, ix, iy, ie, k 
      INTEGER :: ninterv  ! number of (points-1)==no of intervals writen to file
      INTEGER ispk, ixxx, ikk, ima, nma 
      INTEGER nx_min, nx_max, ny_min, ny_max, nx_s, ny_s 
      LOGICAL l_m999 
!                                                                       
!                                                                       
      CALL oeffne (isa, filname, 'unknown') 
      IF (ier_num.ne.0) return 
!                                                                       
!-------xy-kurve abspeichern                                            
!                                                                       
      IF (form (1:2) .eq.'XY'.and..not.lni (ik) ) then 
         WRITE (isa, 3000) titel (iwin, iframe, 1) (1:len_str (titel (  &
         iwin, iframe, 1) ) )                                           
         WRITE (isa, 3000) titel (iwin, iframe, 2) (1:len_str (titel (  &
         iwin, iframe, 2) ) )                                           
         IF(l_two_col) THEN
         DO i = 1, lenc(ik) 
         IF (x (offxy (ik - 1) + i) .ge.werte (1) .and.x (offxy (ik - 1)&
         + i) .le.werte (2) ) then                                      
            WRITE (isa, 4000) x (offxy (ik - 1) + i), y (offxy (ik - 1) &
            + i)
         ENDIF 
         ENDDO 
         ELSE
         DO i = 1, lenc(ik) 
         IF (x (offxy (ik - 1) + i) .ge.werte (1) .and.x (offxy (ik - 1)&
         + i) .le.werte (2) ) then                                      
            WRITE (isa, 4000) x (offxy (ik - 1) + i), y (offxy (ik - 1) &
            + i), dx (offxy (ik - 1) + i), dy (offxy (ik - 1) + i)      
         ENDIF 
         ENDDO 
         ENDIF 
!                                                                       
      ELSEIF (form (1:2) .eq.'DY'.and..not.lni (ik) ) then 
         WRITE (isa, 3000) titel (iwin, iframe, 1) (1:len_str (titel (  &
         iwin, iframe, 1) ) )                                           
         WRITE (isa, 3000) titel (iwin, iframe, 2) (1:len_str (titel (  &
         iwin, iframe, 2) ) )                                           
         DO i = 1, lenc(ik) 
         IF (x (offxy (ik - 1) + i) .ge.werte (1) .and.x (offxy (ik - 1)&
         + i) .le.werte (2) ) then                                      
            WRITE (isa, 4000) x (offxy (ik - 1) + i), y (offxy (ik - 1) &
            + i), dy (offxy (ik - 1) + i)                               
         ENDIF 
         ENDDO 
!                                                                       
!-------nipl-file fuer gnuplot abspeichern                              
!                                                                       
      ELSEIF (form (1:2) .eq.'GN'.and.lni (ik) ) then 
         rdx = (xmax (ik) - xmin (ik) ) / REAL(nx (ik) - 1) 
         rdy = (ymax (ik) - ymin (ik) ) / REAL(ny (ik) - 1) 
         nx_min = max (1, nint ( ( (werte (1) - xmin (ik) ) / rdx) )    &
         + 1)                                                           
         nx_max = min (nx (ik), nint ( ( (werte (2) - xmin (ik) )       &
         / rdx) ) + 1)                                                  
         ny_min = max (1, nint ( ( (werte (3) - ymin (ik) ) / rdy) )    &
         + 1)                                                           
         ny_max = min (ny (ik), nint ( ( (werte (4) - ymin (ik) )       &
         / rdy) ) + 1)                                                  
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
            WRITE (isa, 4100) REAL(ixxx), zzz, 0.0, abs (dzzz),       &
            xxx, yyy                                                    
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
         IF (z (ikk) .ne. - 9999.0.or.z (ikk) .eq. - 9999.0.and.l_m999) &
         then                                                           
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
         IF (z (ikk) .ne. - 9999.0.or.z (ikk) .eq. - 9999.0.and.l_m999) &
         then                                                           
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
         IF (z (ikk) .ne. - 9999.0.or.z (ikk) .eq. - 9999.0.and.l_m999) &
         then                                                           
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
         IF (z (ikk) .ne. - 9999.0.or.z (ikk) .eq. - 9999.0.and.l_m999) &
         then                                                           
            WRITE (isa, 4000) y (offxy (ik - 1) + i), z (ikk) 
         ENDIF 
         ENDDO 
!                                                                       
!-------pgm-file abspeichern (aktueller ausschnitt)                     
!                                                                       
      ELSEIF (form (1:2) .eq.'PG'.and.lni (ik) ) then 
         rdx = (xmax (ik) - xmin (ik) ) / REAL(nx (ik) - 1) 
         rdy = (ymax (ik) - ymin (ik) ) / REAL(ny (ik) - 1) 
         nx_min = max (1, nint ( ( (werte (1) - xmin (ik) ) / rdx) )    &
         + 1)                                                           
         nx_max = min (nx (ik), nint ( ( (werte (2) - xmin (ik) )       &
         / rdx) ) + 1)                                                  
         ny_min = max (1, nint ( ( (werte (3) - ymin (ik) ) / rdy) )    &
         + 1)                                                           
         ny_max = min (ny (ik), nint ( ( (werte (4) - ymin (ik) )       &
         / rdy) ) + 1)                                                  
         WRITE (isa, 2000) 'P2' 
         WRITE (isa, 2000) '# PGM-file created by KUPLOT' 
         WRITE (isa, * ) nx_max - nx_min + 1, ny_max - ny_min + 1 
         WRITE (isa, * ) 255 
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
!                                                                       
!-------nipl-file abspeichern (aktueller ausschnitt)                    
!                                                                       
      ELSEIF (form (1:2) .eq.'NI'.and.lni (ik) ) then 
         rdx = (xmax (ik) - xmin (ik) ) / REAL(nx (ik) - 1) 
         rdy = (ymax (ik) - ymin (ik) ) / REAL(ny (ik) - 1) 
         nx_min = max (1, nint ( ( (werte (1) - xmin (ik) ) / rdx) )    &
         + 1)                                                           
         nx_max = min (nx (ik), nint ( ( (werte (2) - xmin (ik) )       &
         / rdx) ) + 1)                                                  
         ny_min = max (1, nint ( ( (werte (3) - ymin (ik) ) / rdy) )    &
         + 1)                                                           
         ny_max = min (ny (ik), nint ( ( (werte (4) - ymin (ik) )       &
         / rdy) ) + 1)                                                  
         wx_min = xmin (ik) + (nx_min - 1) * rdx 
         wx_max = xmin (ik) + (nx_max - 1) * rdx 
         wy_min = ymin (ik) + (ny_min - 1) * rdy 
         wy_max = ymin (ik) + (ny_max - 1) * rdy 
         WRITE (isa, * ) nx_max - nx_min + 1, ny_max - ny_min + 1 
         WRITE (isa, 4000) wx_min, wx_max, wy_min, wy_max 
         DO iy = ny_min, ny_max 
         WRITE (isa, 4000) (z (offz (ik - 1) + (ix - 1) * ny (ik)       &
         + iy), ix = nx_min, nx_max)                                    
         ENDDO 
      ENDIF 
      CLOSE (isa) 
!                                                                       
 2000 FORMAT     (a) 
 3000 FORMAT     ('#',a) 
 4000 FORMAT     (5(g15.8,2x)) 
 4100 FORMAT     (6(g13.6,1x)) 
 4500 FORMAT     (20(i3,1x)) 
      END SUBROUTINE do_save                        
