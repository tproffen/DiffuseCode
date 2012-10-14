!*****7**************************************************************** 
!     Routines for NeXus support ..                                     
!*****7**************************************************************** 
      SUBROUTINE do_nxinit 
!                                                                       
      include'napif.inc' 
      include'nexus.inc' 
      include'prompt.inc' 
!                                                                       
      nxs_open = .false. 
      WRITE ( * , * ) 'NeXus support enabled ..' 
      END SUBROUTINE do_nxinit                      
!*****7**************************************************************** 
      SUBROUTINE do_nxopen (line, ll) 
!                                                                       
      include'napif.inc' 
      include'nexus.inc' 
      include'prompt.inc' 
      include'debug.inc' 
      include'errlist.inc' 
!                                                                       
      INTEGER maxw 
      PARAMETER (maxw = 10) 
!                                                                       
      CHARACTER ( * ) line 
      INTEGER ll 
!                                                                       
      CHARACTER(1024) cpara (maxw) 
      CHARACTER(1024) ent_name, ent_class 
      INTEGER lpara (maxw), ianz, itype, stat 
      REAL werte (maxw) 
      LOGICAL lexist 
!                                                                       
      INTEGER len_str 
!                                                                       
      IF (nxs_open) then 
         ier_num = - 54 
         ier_typ = ER_APPL 
         RETURN 
      ENDIF 
!                                                                       
      CALL get_params (line, ianz, cpara, lpara, maxw, ll) 
      IF (ier_num.ne.0) return 
      CALL do_build_name (ianz, cpara, lpara, werte, maxw, 1) 
      IF (ier_num.ne.0) return 
!                                                                       
      INQUIRE (file = cpara (1), exist = lexist) 
      IF (lexist) stat = NXopen (cpara (1), NXACC_READ, nxs_id) 
      IF (stat.eq.NX_ERROR.or..not.lexist) then 
         ier_num = - 2 
         ier_typ = ER_IO 
         RETURN 
      ENDIF 
!                                                                       
      nxs_open = .true. 
      nxs_fname = cpara (1) (1:lpara (1) ) 
!                                                                       
      stat = NXgetnextentry (nxs_id, ent_name, ent_class, itype) 
      stat = NXopengroup (nxs_id, ent_name, ent_class) 
      IF (dbg) write ( *, 1000) ent_name (1:len_str (ent_name) ),       &
      ent_class (1:len_str (ent_class) )                                
!                                                                       
      stat = NXopendata (nxs_id, 'title') 
      stat = NXgetchardata (nxs_id, run_title) 
      stat = NXclosedata (nxs_id) 
      stat = NXopendata (nxs_id, 'start_time') 
      stat = NXgetchardata (nxs_id, run_stime) 
      stat = NXclosedata (nxs_id) 
      stat = NXopendata (nxs_id, 'end_time') 
      stat = NXgetchardata (nxs_id, run_etime) 
      stat = NXclosedata (nxs_id) 
      stat = NXopendata (nxs_id, 'experiment_number') 
      stat = NXgetdata (nxs_id, run_iexp, NX_UINT32) 
      stat = NXclosedata (nxs_id) 
      stat = NXopendata (nxs_id, 'run_number') 
      stat = NXgetdata (nxs_id, run_irun, NX_UINT32) 
      stat = NXclosedata (nxs_id) 
!                                                                       
      WRITE (output_io, 2000) cpara (1) (1:lpara (1) ), run_iexp,       &
      run_irun, run_title (1:len_str (run_title) ), run_stime (1:       &
      len_str (run_stime) ), run_etime (1:len_str (run_etime) )         
!                                                                       
 1000 FORMAT    (' debug > NXentry ',a,'(',a,')') 
 2000 FORMAT    (' Reading file ',a,/                                   &
     &                  '   Proposal  : ',i8,/                          &
     &             '   Run       : ',i8,/                               &
     &         '   Title     : ',a,/                                    &
     &         '   Run start : ',a,/                                    &
     &         '   Run end   : ',a)                                     
      END SUBROUTINE do_nxopen                      
!*****7**************************************************************** 
      SUBROUTINE do_nxdir 
!                                                                       
      include'napif.inc' 
      include'nexus.inc' 
      include'prompt.inc' 
      include'errlist.inc' 
!                                                                       
      INTEGER nmax 
      PARAMETER (nmax = 3) 
!                                                                       
      CHARACTER(40) nxs_dname, nxs_aname (nmax) 
      CHARACTER(40) ent_name, ent_class 
      INTEGER nxs_dim (nmax) 
      INTEGER stat, itype, irank, i 
!                                                                       
      INTEGER len_str 
!                                                                       
      IF (.not.nxs_open) then 
         ier_num = - 55 
         ier_typ = ER_APPL 
         RETURN 
      ENDIF 
!                                                                       
      WRITE (output_io, 1000) 
      stat = NX_OK 
      DO while (stat.eq.NX_OK) 
      ent_name = '                                    ' 
      ent_class = 'None' 
      stat = NXgetnextentry (nxs_id, ent_name, ent_class, itype) 
      IF (ent_class.eq.'NXdata'.and.stat.eq.NX_OK) then 
         stat = NXopengroup (nxs_id, ent_name, ent_class) 
         CALL nxs_datainfo (nxs_dname, nxs_aname, nxs_dim, irank, nmax) 
         stat = NXclosegroup (nxs_id) 
         IF (irank.ge.1.and.irank.le.3) then 
            WRITE (output_io, 1500) ent_name 
            DO i = irank, 1, - 1 
      WRITE (output_io, 1510) irank - i + 1, nxs_dim (i),  nxs_aname (i)&
     &  (1:len_str (nxs_aname (i) ) )                                   
            ENDDO 
         ELSE 
            WRITE (output_io, 1530) ent_name 
         ENDIF 
      ENDIF 
      ENDDO 
!                                                                       
      stat = NXinitgroupdir (nxs_id) 
!                                                                       
 1000 FORMAT    (' NXdata blocks in current file :') 
 1500 FORMAT    (1x,79('-'),/,' Name     : ',a) 
 1510 FORMAT    ('   Axis ',i1,' : Dimension : ',i6,'  - Name : ',a) 
 1530 FORMAT    (1x,79('-'),/,' Name     : ',a,' (unsupported rank)') 
      END SUBROUTINE do_nxdir                       
!*****7**************************************************************** 
      SUBROUTINE do_nxclose 
!                                                                       
      include'napif.inc' 
      include'nexus.inc' 
      include'errlist.inc' 
!                                                                       
      INTEGER stat 
!                                                                       
      IF (.not.nxs_open) then 
         ier_num = - 55 
         ier_typ = ER_APPL 
      ELSE 
         stat = NXclosegroup (nxs_id) 
         stat = NXclose (nxs_id) 
         nxs_open = .false. 
      ENDIF 
!                                                                       
      END SUBROUTINE do_nxclose                     
!*****7**************************************************************** 
      SUBROUTINE do_nxload (line, ll) 
!                                                                       
      include'napif.inc' 
      include'nexus.inc' 
      include'errlist.inc' 
!                                                                       
      INTEGER nmax, maxw 
      PARAMETER (nmax = 3) 
      PARAMETER (maxw = 6) 
!                                                                       
      CHARACTER ( * ) line 
      INTEGER ll 
!                                                                       
      CHARACTER(1024) cpara (maxw), clabel, csds 
      INTEGER lpara (maxw), ianz 
      REAL werte (maxw) 
!                                                                       
      CHARACTER(40) nxs_dname, nxs_aname (nmax) 
      INTEGER nxs_dim (nmax) 
      INTEGER istart (3), isize (3) 
      INTEGER inx, iny 
      INTEGER llabel, irank, stat, icol, irow, ioth, i 
!                                                                       
      IF (.not.nxs_open) then 
         ier_num = - 55 
         ier_typ = ER_APPL 
         RETURN 
      ENDIF 
!                                                                       
      CALL get_params (line, ianz, cpara, lpara, maxw, ll) 
      IF (ier_num.ne.0) return 
!                                                                       
      csds = cpara (1) (1:lpara (1) ) 
      stat = NXopengroup (nxs_id, cpara (1) (1:lpara (1) ) , 'NXdata') 
      IF (stat.ne.NX_OK) then 
         ier_num = - 57 
         ier_typ = ER_APPL 
         RETURN 
      ENDIF 
!                                                                       
      CALL nxs_datainfo (nxs_dname, nxs_aname, nxs_dim, irank, nmax) 
!                                                                       
      DO i = 1, 3 
      istart (i) = 0 
      isize (i) = 0 
      ENDDO 
!                                                                       
!------ Replace * with 0.0                                              
!                                                                       
      IF (ianz.gt.1) then 
         cpara (1) = '0.0' 
         lpara (1) = 3 
         DO i = 2, ianz 
         IF (cpara (i) .eq." * ") then 
            cpara (i) = "0.0" 
            lpara (i) = 3 
         ENDIF 
         ENDDO 
         CALL ber_params (ianz, cpara, lpara, werte, maxw) 
         IF (ier_num.ne.0) return 
      ENDIF 
!                                                                       
!------ SDS is 1D                                                       
!                                                                       
      IF (irank.eq.1) then 
         IF (ianz.eq.1) then 
            icol = 1 
            istart (icol) = 1 
            isize (icol) = nxs_dim (icol) 
            clabel = nxs_fname (1:len_str (nxs_fname) ) //':'//csds (1: &
            len_str (csds) )                                            
            CALL nxs_read1d (nxs_dname, nxs_aname (icol), clabel,       &
            istart, isize, icol, irank)                                 
         ELSE 
            ier_num = - 6 
            ier_typ = ER_COMM 
         ENDIF 
!                                                                       
!------ SDS is 2D                                                       
!                                                                       
      ELSEIF (irank.eq.2) then 
!                                                                       
!------ - 2D complete                                                   
!                                                                       
         IF (ianz.eq.1.or.ianz.eq.5.or. (ianz.eq.3.and.werte (2)        &
         .gt.0.0.and.werte (3) .gt.0.0) ) then                          
!                                                                       
            IF (ianz.eq.1) then 
               inx = 1 
               iny = 1 
               istart (1) = 1 
               istart (2) = 1 
               isize (1) = nxs_dim (1) 
               isize (2) = nxs_dim (2) 
            ELSEIF (ianz.eq.3) then 
               inx = nint (werte (2) ) 
               iny = nint (werte (3) ) 
               istart (1) = 1 
               istart (2) = 1 
               isize (1) = nxs_dim (1) 
               isize (2) = nxs_dim (2) 
            ELSEIF (ianz.eq.5) then 
               inx = 1 
               iny = 1 
               istart (1) = nint (werte (2) ) 
               istart (2) = nint (werte (4) ) 
               isize (1) = nint (werte (3) ) - istart (1) + 1 
               isize (2) = nint (werte (5) ) - istart (2) + 1 
            ENDIF 
!                                                                       
            icol = 1 
            irow = 2 
!                                                                       
            clabel = nxs_fname (1:len_str (nxs_fname) ) //':'//csds (1: &
            len_str (csds) )                                            
            CALL nxs_read2d (nxs_dname, nxs_aname (1), nxs_aname (2),   &
            clabel, istart, isize, nxs_dim (1), nxs_dim (2), inx, iny,  &
            icol, irow, irank)                                          
!                                                                       
!------ - 2D sections (x or y)                                          
!                                                                       
         ELSEIF (ianz.eq.3.and. (werte (2) .eq.0.0.xor.werte (3)        &
         .eq.0.0) ) then                                                
            IF (werte (2) .eq.0.0) then 
               icol = 2 
               irow = 1 
               istart (irow) = nint (werte (3) ) 
            ELSE 
               icol = 1 
               irow = 2 
               istart (irow) = nint (werte (2) ) 
            ENDIF 
            istart (icol) = 1 
            isize (irow) = 1 
            isize (icol) = nxs_dim (icol) 
!                                                                       
            IF (istart (irow) .le.0.or.istart (irow) .gt.nxs_dim (irow) &
            ) then                                                      
               ier_num = - 58 
               ier_typ = ER_APPL 
            ELSE 
               WRITE (clabel, 1000) nxs_fname (1:len_str (nxs_fname) ), &
               csds (1:len_str (csds) ), nxs_aname (irow), istart (irow)
               llabel = len_str (clabel) 
               CALL rem_bl (clabel, llabel) 
               CALL nxs_read1d (nxs_dname, nxs_aname (icol), clabel,    &
               istart, isize, icol, irank)                              
            ENDIF 
         ENDIF 
!                                                                       
!------ SDS is 3D                                                       
!                                                                       
      ELSEIF (irank.eq.3) then 
!                                                                       
!------ - 3D sections (x, y or z)                                       
!                                                                       
         IF (ianz.eq.4.and. (werte (2) .eq.0.0.xor.werte (3)            &
         .eq.0.0.xor.werte (4) .eq.0.0) ) then                          
            IF (werte (2) .eq.0.0) then 
               icol = 3 
               irow = 2 
               ioth = 1 
               istart (irow) = nint (werte (3) ) 
               istart (ioth) = nint (werte (4) ) 
            ELSEIF (werte (3) .eq.0.0) then 
               icol = 2 
               irow = 1 
               ioth = 3 
               istart (irow) = nint (werte (4) ) 
               istart (ioth) = nint (werte (2) ) 
            ELSEIF (werte (4) .eq.0.0) then 
               icol = 1 
               irow = 2 
               ioth = 3 
               istart (irow) = nint (werte (3) ) 
               istart (ioth) = nint (werte (2) ) 
            ENDIF 
!                                                                       
            istart (icol) = 1 
!                                                                       
            isize (icol) = nxs_dim (icol) 
            isize (irow) = 1 
            isize (ioth) = 1 
!                                                                       
            IF (istart (irow) .le.0.or.istart (irow) .gt.nxs_dim (irow) &
            .or.istart (ioth) .le.0.or.istart (ioth) .gt.nxs_dim (ioth) &
            ) then                                                      
               ier_num = - 58 
               ier_typ = ER_APPL 
            ELSE 
               WRITE (clabel, 1100) nxs_fname (1:len_str (nxs_fname) ), &
               csds (1:len_str (csds) ), nxs_aname (irow), istart (irow)&
               , nxs_aname (ioth), istart (ioth)                        
               llabel = len_str (clabel) 
               CALL rem_bl (clabel, llabel) 
               CALL nxs_read1d (nxs_dname, nxs_aname (icol), clabel,    &
               istart, isize, icol, irank)                              
            ENDIF 
         ENDIF 
!                                                                       
      ELSE 
         ier_num = - 56 
         ier_typ = ER_APPL 
      ENDIF 
!                                                                       
      stat = NXclosegroup (nxs_id) 
      stat = NXinitgroupdir (nxs_id) 
!                                                                       
 1000 FORMAT    (a,':',a,':',a,'=',i8) 
 1100 FORMAT    (a,':',a,':',a,'=',i8,':',a,'=',i8) 
      END SUBROUTINE do_nxload                      
!*****7**************************************************************** 
!     Helper routines                                                   
!*****7**************************************************************** 
      SUBROUTINE nxs_read2d (dname, xname, yname, label, istart, isize, &
      dimx, dimy, inx, iny, icol, irow, irank)                          
!                                                                       
      include'config.inc' 
      include'kuplot.inc' 
      include'napif.inc' 
      include'nexus.inc' 
      include'debug.inc' 
      include'errlist.inc' 
!                                                                       
      CHARACTER ( * ) dname, xname, yname, label 
      INTEGER istart (3), isize (3) 
      INTEGER iistart (3), iisize (3) 
      INTEGER irank, irow, inx, iny, icol, dimx, dimy 
!                                                                       
      REAL rval (maxarray) 
      INTEGER ival (maxarray) 
      INTEGER stat, i, j, k, l, m, ii, jj, kk, itype, nr 
      INTEGER maxpkt, maxzz 
!                                                                       
      INTEGER len_str 
!                                                                       
      IF (dbg) then 
      WRITE ( * , 8000) 'Data   : ', dname (1:len_str (dname) ) 
         WRITE ( * , 8000) 'Axis X : ', xname (1:len_str (xname) ) 
         WRITE ( * , 8000) 'Axis Y : ', yname (1:len_str (yname) ) 
         WRITE ( * , 8100) 'Istart : ', istart 
      WRITE ( * , 8100) 'Isize  : ', isize 
      WRITE ( * , 8200) 'In X   : ', inx 
      WRITE ( * , 8200) 'In Y   : ', iny 
      WRITE ( * , 8200) 'Dim X  : ', dimx 
      WRITE ( * , 8200) 'Dim Y  : ', dimy 
      WRITE ( * , 8200) 'Irow   : ', irow 
      WRITE ( * , 8200) 'Icol   : ', icol 
      WRITE ( * , 8200) 'Irank  : ', irank 
!                                                                       
 8000 FORMAT      (' debug > ',a,a) 
 8100 FORMAT      (' debug > ',a,3(2x,i5)) 
 8200 FORMAT      (' debug > ',a,2x,i5) 
      ENDIF 
!                                                                       
!------ Space for desired data set ?                                    
!                                                                       
      maxpkt = maxarray - offxy (iz - 1) 
      maxzz = maxarray - offz (iz - 1) 
!                                                                       
      nx (iz) = isize (icol) / inx 
      ny (iz) = isize (irow) / iny 
!                                                                       
      IF (max (nx (iz), ny (iz) ) .gt.maxzz.or.max (nx (iz), ny (iz) )  &
      .gt.maxpkt) then                                                  
         ier_num = - 6 
         ier_typ = ER_APPL 
         GOTO 999 
      ENDIF 
!                                                                       
!------ Read data                                                       
!                                                                       
      stat = NXopendata (nxs_id, dname) 
      IF (stat.ne.NX_OK) goto 55 
!                                                                       
      jj = 1 
      iistart (icol) = istart (icol) 
      iistart (irow) = istart (irow) 
      iisize (icol) = isize (icol) 
      iisize (irow) = 1 
!                                                                       
      DO j = 1, isize (irow), iny 
      kk = 1 
      DO k = 1, iny 
      stat = NXgetslab (nxs_id, ival, iistart, iisize) 
      IF (stat.ne.NX_OK) goto 55 
      DO l = 1, isize (icol), inx 
      DO m = 2, inx 
      ival (kk) = ival (l) + ival (l + m - 1) 
      ENDDO 
      kk = kk + 1 
      ENDDO 
      DO i = 1, nx (iz) 
      z (offz (iz - 1) + (i - 1) * ny (iz) + jj) = float (ival (i) )    &
      / inx / iny                                                       
      ENDDO 
      jj = jj + 1 
      iistart (irow) = iistart (irow) + 1 
      ENDDO 
      ENDDO 
!                                                                       
      stat = NXclosedata (nxs_id) 
!                                                                       
!------ Read values for X and Y                                         
!                                                                       
      kk = 1 
!                                                                       
      stat = NXopendata (nxs_id, xname) 
      IF (stat.ne.NX_OK) goto 55 
      stat = NXgetslab (nxs_id, rval, istart (icol), isize (icol) ) 
      IF (stat.ne.NX_OK) goto 55 
      DO l = 1, isize (icol), inx 
      DO m = 2, inx 
      rval (kk) = rval (l) + rval (l + m - 1) 
      ENDDO 
      kk = kk + 1 
      ENDDO 
      DO i = 1, nx (iz) 
      x (offxy (iz - 1) + i) = rval (i) / float (inx) 
      ENDDO 
      stat = NXclosedata (nxs_id) 
!                                                                       
      kk = 1 
!                                                                       
      stat = NXopendata (nxs_id, yname) 
      IF (stat.ne.NX_OK) goto 55 
      stat = NXgetslab (nxs_id, rval, istart (irow), isize (irow) ) 
      IF (stat.ne.NX_OK) goto 55 
      DO l = 1, isize (irow), iny 
      DO m = 2, iny 
      rval (kk) = rval (l) + rval (l + m - 1) 
      ENDDO 
      kk = kk + 1 
      ENDDO 
      DO i = 1, ny (iz) 
      y (offxy (iz - 1) + i) = rval (i) / float (iny) 
      ENDDO 
      stat = NXclosedata (nxs_id) 
!                                                                       
!------ set remaining parameters                                        
!                                                                       
      lni (iz) = .true. 
      len (iz) = max (nx (iz), ny (iz) ) 
      offxy (iz) = offxy (iz - 1) + len (iz) 
      offz (iz) = offz (iz - 1) + nx (iz) * ny (iz) 
      fname (iz) = label (1:len_str (label) ) 
      iz = iz + 1 
!                                                                       
      titel (iwin, iframe, 1) = run_title (1:len_str (run_title) ) 
      achse (iwin, iframe, 1) = xname (1:len_str (xname) ) 
      achse (iwin, iframe, 2) = yname (1:len_str (yname) ) 
      achse (iwin, iframe, 3) = dname (1:len_str (dname) ) 
!                                                                       
      CALL show_data (iz - 1) 
      GOTO 999 
!                                                                       
   50 CONTINUE 
      ier_num = - 6 
      ier_typ = ER_IO 
      GOTO 999 
!                                                                       
   55 CONTINUE 
      ier_num = - 3 
      ier_typ = ER_IO 
      GOTO 999 
!                                                                       
  999 CONTINUE 
      END SUBROUTINE nxs_read2d                     
!*****7**************************************************************** 
      SUBROUTINE nxs_read1d (dname, aname, label, istart, isize, irow,  &
      irank)                                                            
!                                                                       
      include'config.inc' 
      include'kuplot.inc' 
      include'napif.inc' 
      include'nexus.inc' 
      include'debug.inc' 
      include'errlist.inc' 
!                                                                       
      CHARACTER ( * ) dname, aname, label 
      INTEGER istart (3), isize (3) 
      INTEGER irank, irow 
!                                                                       
      REAL rval (maxarray) 
      INTEGER ival (maxarray) 
      INTEGER stat, i, itype, nr, maxpp, ntot 
!                                                                       
      INTEGER len_str 
!                                                                       
      maxpp = maxarray - offxy (iz - 1) 
      IF (isize (irow) .gt.maxpp) then 
         ier_num = - 6 
         ier_typ = ER_APPL 
         RETURN 
      ENDIF 
!                                                                       
      IF (dbg) then 
      WRITE ( * , 8000) 'Data   : ', dname (1:len_str (dname) ) 
      WRITE ( * , 8000) 'Axis   : ', aname (1:len_str (aname) ) 
         WRITE ( * , 8100) 'Istart : ', istart 
      WRITE ( * , 8100) 'Isize  : ', isize 
      WRITE ( * , 8200) 'Irow   : ', irow 
      WRITE ( * , 8200) 'Irank  : ', irank 
!                                                                       
 8000 FORMAT      (' debug > ',a,a) 
 8100 FORMAT      (' debug > ',a,3(2x,i5)) 
 8200 FORMAT      (' debug > ',a,2x,i5) 
      ENDIF 
!                                                                       
!------ read data SDS                                                   
!                                                                       
      stat = NXopendata (nxs_id, dname) 
      IF (stat.ne.NX_OK) goto 90 
      stat = NXgetslab (nxs_id, ival, istart, isize) 
      IF (stat.ne.NX_OK) goto 90 
      stat = NXclosedata (nxs_id) 
!                                                                       
!------ read x-axis data SDS                                            
!                                                                       
      stat = NXopendata (nxs_id, aname) 
      IF (stat.ne.NX_OK) goto 90 
      stat = NXgetdata (nxs_id, rval) 
      IF (stat.ne.NX_OK) goto 90 
      stat = NXclosedata (nxs_id) 
!                                                                       
      DO nr = 1, isize (irow) 
      x (offxy (iz - 1) + nr) = rval (nr) 
      y (offxy (iz - 1) + nr) = float (ival (nr) ) 
      dx (offxy (iz - 1) + nr) = 0.0 
      dy (offxy (iz - 1) + nr) = 0.0 
      ENDDO 
!                                                                       
      len (iz) = isize (irow) 
      offxy (iz) = offxy (iz - 1) + len (iz) 
      offz (iz) = offz (iz - 1) 
      fname (iz) = label (1:len_str (label) ) 
      fform (iz) = 'XY' 
      lni (iz) = .false. 
      iz = iz + 1 
!                                                                       
      titel (iwin, iframe, 1) = run_title (1:len_str (run_title) ) 
      achse (iwin, iframe, 1) = aname (1:len_str (aname) ) 
      achse (iwin, iframe, 2) = dname (1:len_str (dname) ) 
!                                                                       
      CALL show_data (iz - 1) 
      RETURN 
!                                                                       
   90 CONTINUE 
      ier_num = - 3 
      ier_typ = ER_IO 
!                                                                       
      END SUBROUTINE nxs_read1d                     
!*****7**************************************************************** 
      SUBROUTINE nxs_datainfo (nxs_dname, nxs_aname, nxs_dim, irank,    &
      nmax)                                                             
!                                                                       
      include'napif.inc' 
      include'nexus.inc' 
      include'errlist.inc' 
!                                                                       
      INTEGER nmax 
!                                                                       
      CHARACTER(40) nxs_dname, nxs_aname (nmax) 
      CHARACTER(40) ent_name, ent_class, attr_name 
      CHARACTER(40) cdummy 
      INTEGER idimension (NX_MAXRANK) 
      INTEGER nxs_dim (nmax) 
      INTEGER estat, astat, stat 
      INTEGER itype, ilen, irank, ival, i, idummy 
!                                                                       
      irank = 0 
      estat = NX_OK 
!                                                                       
      DO while (estat.eq.NX_OK) 
      ent_name = '                                    ' 
      ent_class = 'None' 
      estat = NXgetnextentry (nxs_id, ent_name, ent_class, itype) 
      IF (ent_class.eq.'SDS'.and.estat.eq.NX_OK) then 
         astat = NXopendata (nxs_id, ent_name) 
         DO while (astat.eq.NX_OK) 
         attr_name = 'None' 
         astat = NXgetnextattr (nxs_id, attr_name, ilen, itype) 
!                                                                       
         IF (astat.eq.NX_OK) then 
            IF (attr_name.eq.'signal') then 
               ilen = 1 
               stat = NXgetinfo (nxs_id, irank, idimension, itype) 
               stat = NXgetattr (nxs_id, attr_name, ival, ilen, itype) 
               IF (irank.gt.nmax) then 
                  ier_num = - 56 
                  ier_typ = ER_APPL 
               ELSE 
                  nxs_dname = ent_name 
                  DO i = 1, irank 
                  nxs_dim (i) = idimension (i) 
                  ENDDO 
               ENDIF 
!                                                                       
            ELSEIF (attr_name.eq.'axis') then 
               ilen = 1 
               stat = NXgetattr (nxs_id, attr_name, ival, ilen, itype) 
               IF (ival.le.nmax) then 
                  nxs_aname (ival) = ent_name 
               ELSE 
                  ier_num = - 56 
                  ier_typ = ER_APPL 
               ENDIF 
            ENDIF 
         ENDIF 
         ENDDO 
         stat = NXclosedata (nxs_id) 
      ENDIF 
      ENDDO 
!                                                                       
!------ Sort the names backwards                                        
!                                                                       
      IF (irank.eq.2) then 
         cdummy = nxs_aname (1) 
         nxs_aname (1) = nxs_aname (2) 
         nxs_aname (2) = cdummy 
      ELSEIF (irank.eq.3) then 
         cdummy = nxs_aname (1) 
         nxs_aname (1) = nxs_aname (3) 
         nxs_aname (3) = cdummy 
      ENDIF 
!                                                                       
      END SUBROUTINE nxs_datainfo                   
