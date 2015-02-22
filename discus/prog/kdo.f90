!*****7*****************************************************************
      SUBROUTINE discus_mache_kdo (line, lend, length) 
!+                                                                      
!     This is the main routine for command interpretation, each         
!     command is identified here and the corresponding subroutine       
!     executed. A leading @ indicates a macro.                          
!-                                                                      
      USE addfile_mod
      USE discus_allocate_appl_mod
      USE charact_mod 
      USE chem_menu
      USE conn_mod
      USE domain_menu
      USE fourier_menu
      USE insert_menu
      USE metric_mod
      USE mmc_menu
      USE interpret_menu
      USE inverse_mod 
      USE modify_mod
      USE mole_surf_mod
      USE output_menu
      USE patters_menu
      USE pdf_menu
      USE plot_menu
      USE powder
      USE rmc_menu
      USE save_menu
      USE show_menu
      USE stack_menu
      USE symm_menu
      USE shear
      USE structur
      USE spcgr_apply, ONLY: wyckoff_main
      USE thermal_mod
      USE transform_menu
      USE waves_do_menu
!
      USE errlist_mod 
      USE learn_mod 
      IMPLICIT none 
!                                                                       
      CHARACTER (LEN= * ), INTENT(INOUT) :: line 
      LOGICAL            , INTENT(OUT)   :: lend
      INTEGER            , INTENT(INOUT) :: length 
!                                                                       
      CHARACTER(1024) zeile 
      CHARACTER(5) befehl 
      INTEGER indxb, indxg, lcomm, lbef 
      INTEGER                  :: indxt ! position of a TAB
      INTEGER                  ::  inverse_type
      LOGICAL lout_rho, lkick 
      LOGICAL str_comp 
!                                                                       
      DATA lout_rho / .false. / 
!                                                                       
      CALL no_error 
!                                                                       
!-------If a commentary return immediately                              
!                                                                       
      IF (line (1:1)  == ' '.or.line (1:1)  == '#' .or.   & 
          line == char(13) .or. line(1:1) == '!'  ) RETURN
!                                                                       
!     Only the first 5 characters are significant. The command consists 
!     of the four nonblank characters                                   
!                                                                       
      befehl = '    ' 
      indxt = INDEX (line, tab)       ! find a tabulator
      IF(indxt==0) indxt = length + 1
      indxb = index (line, ' ')       ! find a blank
      IF(indxb==0) indxb = length + 1
      indxb = MIN(indxb,indxt)
      lbef = min (indxb - 1, 5) 
      befehl = line (1:lbef) 
!                                                                       
!------ command parameters start at the first character following       
!       the blank                                                       
!                                                                       
      zeile = ' ' 
      lcomm = 0 
      IF (indxb + 1.le.length) then 
         zeile = line (indxb + 1:length) 
         lcomm = length - indxb 
         CALL rem_leading_bl ( zeile, lcomm )
      ENDIF 
!                                                                       
!-------Suche nach einem "="                                            
!                                                                       
      indxg = index (line, '=') 
      IF (indxg.ne.0.and..not. (str_comp (befehl, 'echo', 2, lbef, 4) ) &
     &.and..not. (str_comp (befehl, 'syst', 2, lbef, 4) ) .and..not. (st&
     &r_comp (befehl, 'fput', 2, lbef, 4) ) .and..not. (str_comp (befehl&
     &, 'socket', 2, lbef, 5) ) .and..not. (str_comp (befehl, 'help', 2,&
     & lbef, 4) .or.str_comp (befehl, '?   ', 2, lbef, 4) ) ) then      
!                                                                       
!-------Zuweisung eines Funktionswertes                                 
!                                                                       
         CALL do_math (line, indxg, length) 
      ELSE 
!                                                                       
!     --execute a macro file                                            
!                                                                       
         IF (befehl (1:1) .eq.'@') then 
            IF (length.ge.2) then 
               CALL file_kdo (line (2:length), length - 1) 
            ELSE 
               ier_num = - 13 
               ier_typ = ER_MAC 
            ENDIF 
!
!     -- Allocate array sizes
!
         ELSEIF (str_comp (befehl, 'allocate', 3, lbef,  8) ) then
            CALL discus_do_allocate_appl (zeile, lcomm)
!
!     -- Deallocate array sizes
!
         ELSEIF (str_comp (befehl, 'deallocate', 3, lbef, 10) ) then
            CALL discus_do_deallocate_appl (zeile, lcomm)
!                                                                       
!-------add two files 'addf'                                            
!                                                                       
         ELSEIF (str_comp (befehl, 'addf', 2, lbef, 4) ) then 
            CALL do_addfile (zeile, lcomm) 
!                                                                       
!     append a new atom 'appe'                                          
!                                                                       
         ELSEIF (str_comp (befehl, 'appe', 2, lbef, 4) ) then 
            lkick = .false. 
            CALL do_app (zeile, lcomm, lkick) 
!                                                                       
!-------show the asymmetric unit 'asym'                                 
!                                                                       
         ELSEIF (str_comp (befehl, 'asym', 2, lbef, 4) ) then 
            CALL show_asym 
!                                                                       
!-------remove atomes outside a boundary 'boundary'                     
!                                                                       
         ELSEIF (str_comp (befehl, 'boundary', 2, lbef, 8) ) then 
            CALL boundary (zeile, lcomm) 
!                                                                       
!-------change some properties in the crystal 'change'                  
!                                                                       
         ELSEIF (str_comp (befehl, 'change', 3, lbef, 6) ) then 
            CALL do_change (zeile, lcomm) 
!                                                                       
!-------show the atoms present in the crystal 'chem'                    
!                                                                       
         ELSEIF (str_comp (befehl, 'chem', 3, lbef, 4) ) then 
            CALL chem 
!                                                                       
!-------go to the connectivity menue  'connectivity'                    
!                                                                       
         ELSEIF (str_comp (befehl, 'connectivity', 3, lbef, 12) ) then 
            CALL conn_menu 
!                                                                       
!-------Define some DISCUS parameters 'define'                          
!                                                                       
         ELSEIF (str_comp (befehl, 'define', 3, lbef, 6) ) then 
            CALL do_define (zeile, lcomm) 
!                                                                       
!-------Decorate a surface by molecules 'decorate'                      
!                                                                       
         ELSEIF (str_comp (befehl, 'decorate', 3, lbef, 8) ) then 
            CALL do_place_molecule
!                                                                       
!-------Handling of domains within the host structure 'domain'          
!                                                                       
         ELSEIF (str_comp (befehl, 'domain', 3, lbef, 6) ) then 
            CALL do_domain (zeile, lcomm) 
!                                                                       
!-------copy an atom 'copy'                                             
!                                                                       
         ELSEIF (str_comp (befehl, 'copy', 3, lbef, 4) ) then 
            CALL do_copy (zeile, lcomm) 
!                                                                       
!     Transform vector from direct to reciprocal space                  
!                                                                       
         ELSEIF (str_comp (befehl, 'd2r', 2, lbef, 3) ) then 
            CALL d2r (zeile, lcomm, .true.) 
!                                                                       
!     Difference  Fourier 'diff'                                        
!                                                                       
         ELSEIF (str_comp (befehl, 'diff', 2, lbef, 4) ) then 
            inverse_type = INV_DIFF 
            lout_rho = .true. 
            CALL patterson (inverse_type) 
!                                                                       
!-------Terminate DISCUS 'exit'                                         
!                                                                       
         ELSEIF (str_comp (befehl, 'exit', 2, lbef, 4) ) then 
            LEND = .TRUE. 
!                                                                       
!-------Find properties   'find'                                        
!                                                                       
         ELSEIF (str_comp (befehl, 'find', 2, lbef, 4) ) then 
            CALL do_find (zeile, lcomm) 
!                                                                       
!-------Fourier transform 'four'                                        
!                                                                       
         ELSEIF (str_comp (befehl, 'four', 2, lbef, 4) ) then 
            lout_rho = .false. 
            CALL fourier 
!                                                                       
!-------import a file into discus.cell format 'import'                  
!                                                                       
         ELSEIF (str_comp (befehl, 'import', 3, lbef, 6) ) then 
            CALL do_import (zeile, lcomm) 
!                                                                       
!-------calculate cell and site <==> atomindex
!                                                                       
         ELSEIF (str_comp (befehl, 'index2cell', 3, lbef, 10) ) then 
            CALL chem_trans (zeile, lcomm) 
!                                                                       
!     inserte a new atom 'inse'                                         
!                                                                       
         ELSEIF (str_comp (befehl, 'inse', 3, lbef, 4) ) then 
            IF (str_comp (zeile, 'domain', 3, lcomm, 6) ) then 
               CALL insert ( - 1) 
            ELSEIF (str_comp (zeile, 'molecule', 3, lcomm, 8) ) then 
               CALL insert (0) 
            ELSEIF (str_comp (zeile, 'object', 3, lcomm, 6) ) then 
               CALL insert (1) 
            ELSE 
               CALL do_ins (zeile, lcomm) 
            ENDIF 
!                                                                       
!     interpret an electron density 'interpret'                         
!                                                                       
         ELSEIF (str_comp (befehl, 'inte', 3, lbef, 4) ) then 
            CALL interpret 
!                                                                       
!     inverse Fourier 'inverse'                                         
!                                                                       
         ELSEIF (str_comp (befehl, 'inve', 3, lbef, 4) ) then 
            inverse_type = INV_INV 
            lout_rho = .true. 
            CALL patterson (inverse_type) 
!                                                                       
!     append a new atom while kicking out old atoms 'kick'              
!                                                                       
         ELSEIF (str_comp (befehl, 'kick', 2, lbef, 4) ) then 
            lkick = .true. 
            CALL do_app (zeile, lcomm, lkick) 
!                                                                       
!     Output routines 'outp'                                            
!                                                                       
         ELSEIF (str_comp (befehl, 'outp', 1, lbef, 4) ) then 
            CALL do_niplps (lout_rho) 
!                                                                       
!     Patterson 'patterson'                                             
!                                                                       
         ELSEIF (str_comp (befehl, 'patt', 2, lbef, 4) ) then 
            inverse_type = INV_PATT 
            lout_rho = .true. 
            CALL patterson (inverse_type) 
!                                                                       
!     PDF calculation and analysis                                      
!                                                                       
         ELSEIF (str_comp (befehl, 'pdf', 2, lbef, 3) ) then 
            CALL pdf 
!                                                                       
!     Plot the crystal 'plot'                                           
!                                                                       
         ELSEIF (str_comp (befehl, 'plot', 2, lbef, 4) ) then 
            CALL plot 
!                                                                       
!     Calculate Powder pattern                                          
!                                                                       
         ELSEIF (str_comp (befehl, 'powder', 2, lbef, 6) ) then 
            CALL do_powder 
!                                                                       
!     Project a vector onto another and onto a plane 'proj'             
!                                                                       
         ELSEIF (str_comp (befehl, 'proj', 4, lbef, 4) ) then 
            CALL do_proj (zeile, lcomm) 
!                                                                       
!     Go to property menu 'property'                                    
!                                                                       
         ELSEIF (str_comp (befehl, 'property', 4, lbef, 8) ) then 
            CALL property_menu 
!                                                                       
!     Purge the list of atoms 'purg'                                    
!                                                                       
         ELSEIF (str_comp (befehl, 'purg', 2, lbef, 4) ) then 
            CALL do_purge 
!                                                                       
!     Transform vector from reciprocal to direct space                  
!                                                                       
         ELSEIF (str_comp (befehl, 'r2d', 2, lbef, 3) ) then 
            CALL d2r (zeile, lcomm, .false.) 
!                                                                       
!-------Einlesen einer Struktur/Zelle 'read'                            
!                                                                       
         ELSEIF (str_comp (befehl, 'read', 3, lbef, 4) ) then 
            CALL read_struc 
!                                                                       
!     Remove a single atom 'remo'                                       
!                                                                       
         ELSEIF (str_comp (befehl, 'remo', 3, lbef, 4) ) then 
            CALL do_remove (zeile, lcomm) 
!                                                                       
!     Replace atom(s) or molecule(s)                                    
!                                                                       
         ELSEIF (str_comp (befehl, 'replace', 3, lbef, 7) ) then 
            CALL do_replace (zeile, lcomm) 
!                                                                       
!     Do Monte Carlo simulations 'rmc', 'mc', 'mmc' and 'amc'           
!                                                                       
         ELSEIF (str_comp (befehl, 'rmc', 2, lbef, 3) ) then 
            CALL rmc 
!        ELSEIF (str_comp (befehl, 'mc', 2, lbef, 2) ) then 
!           CALL mc 
         ELSEIF (str_comp (befehl, 'mmc', 2, lbef, 3) ) then 
            CALL mmc 
!                                                                       
!     save structure to file 'save'                                     
!                                                                       
         ELSEIF (str_comp (befehl, 'save', 2, lbef, 4) ) then 
            CALL save_struc (zeile, lcomm) 
!                                                                       
!     generalized shear operation 'shear'                               
!                                                                       
         ELSEIF (str_comp (befehl, 'shear', 3, lbef, 4) ) then 
            CALL shear_menue
!                                                                       
!     Show something                                 'show'             
!                                                                       
         ELSEIF (str_comp (befehl, 'show', 2, lbef, 4) ) then 
            CALL discus_do_show (zeile, lcomm) 
!                                                                       
!     Go to stacking fault menu 'stack'                                 
!                                                                       
         ELSEIF (str_comp (befehl, 'stac', 2, lbef, 4) ) then 
            CALL stack 
!                                                                       
!     Go to surface menu 'surface'                                      
!                                                                       
         ELSEIF (str_comp (befehl, 'surface', 2, lbef, 8) ) then 
            CALL surface_menu 
!                                                                       
!     switch two atoms 'swit'                                           
!                                                                       
         ELSEIF (str_comp (befehl, 'swit', 2, lbef, 4) ) then 
            CALL do_switch (zeile, lcomm) 
!                                                                       
!     generalized symmetry operation 'symm'                             
!                                                                       
         ELSEIF (str_comp (befehl, 'symm', 3, lbef, 4) ) then 
            CALL symm 
!                                                                       
!     Thermal displacement of all atoms 'ther'                          
!                                                                       
         ELSEIF (str_comp (befehl, 'ther', 2, lbef, 4) ) then 
            CALL ther_displ (zeile, lcomm) 
!                                                                       
!     unit cell transformations      'tran'                             
!                                                                       
         ELSEIF (str_comp (befehl, 'tran', 2, lbef, 4) ) then 
            CALL transform 
!                                                                       
!       Vector Product 'vprod'                                          
!                                                                       
         ELSEIF (str_comp (befehl, 'vprod', 1, lbef, 5) ) then 
            CALL vprod (zeile, lcomm) 
!                                                                       
!------   Waves traveling through the crystal 'wave'                    
!                                                                       
         ELSEIF (str_comp (befehl, 'wave', 3, lbef, 4) ) then 
            CALL waves_menu
!                                                                       
!       Determine Wyckoff symmetry                                      
!                                                                       
         ELSEIF (str_comp (befehl, 'wyckoff', 1, lbef, 7) ) then 
            CALL wyckoff_main (zeile, lcomm) 
!                                                                       
!------   Try general commands                                          
!                                                                       
         ELSE 
            CALL kdo_all (befehl, lbef, zeile, lcomm) 
         ENDIF 
      ENDIF 
!                                                                       
      END SUBROUTINE discus_mache_kdo                      
!*****7**************************************************************** 
      SUBROUTINE do_define (zeile, lp) 
!-                                                                      
!     Sets the value of status variables, DISCUS specific routine       
!+                                                                      
      USE unitcell_mod 
      USE errlist_mod 
      IMPLICIT none 
!                                                                       
!                                                                       
      INTEGER maxw 
      PARAMETER (maxw = 5) 
!                                                                       
      CHARACTER ( * ) zeile 
      CHARACTER(1024) cpara (maxw) 
      INTEGER lpara (maxw), lp 
      INTEGER ianz 
!                                                                       
      LOGICAL str_comp 
!                                                                       
      IF (zeile.ne.' ') then 
         CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
         IF (ier_num.eq.0) then 
!                                                                       
!----- ---- define error                                                
!                                                                       
            IF (str_comp (cpara (1) , 'gener', 1, lpara (1) , 5) ) then 
               IF (ianz.eq.2) then 
                  IF (str_comp (cpara (2) , 'center', 2, lpara (2) , 5) &
                  ) then                                                
                     gen_sta = GEN_CENTER 
                  ELSEIF (str_comp (cpara (2) , 'symmetry', 2, lpara (2)&
                  , 8) ) then                                           
                     gen_sta = GEN_SYMM 
                  ELSE 
                     ier_num = - 6 
                     ier_typ = ER_COMM 
                  ENDIF 
               ELSE 
                  ier_num = - 6 
                  ier_typ = ER_COMM 
               ENDIF 
!                                                                       
            ELSE 
               ier_num = - 6 
               ier_typ = ER_COMM 
            ENDIF 
         ENDIF 
      ELSE 
         ier_num = - 6 
         ier_typ = ER_COMM 
      ENDIF 
!                                                                       
      END SUBROUTINE do_define                      
