!module discus_mache_kdo_mod
!
!contains
!
!*****7*****************************************************************
SUBROUTINE discus_mache_kdo (line, lend, length) 
!+                                                                      
!     This is the main routine for command interpretation, each         
!     command is identified here and the corresponding subroutine       
!     executed. A leading @ indicates a macro.                          
!-                                                                      
USE addfile_mod
USE discus_allocate_appl_mod
USE discus_reset_all_mod
use bragg_mod
USE charact_mod 
USE chem_menu
USE chem_aver_mod, ONLY: get_displacement
USE conn_mod
USE demolec
use discus_3dpdf_mod
USE do_find_top
USE domain_menu
use exp2pdf_menu
USE fourier_menu
USE insert_menu
USE metric_mod
USE mmc_menu
USE interpret_menu
USE inverse_mod 
USE modify_mod
USE molecule_func_mod
USE mole_surf_mod
USE output_menu
USE patters_menu
use perioditize_mod
USE pdf_menu
use prep_anis_mod
USE discus_plot_menu
USE powder_top_mod
USE prop_para_func
USE rmc_menu
USE save_menu
USE discus_show_menu
USE stack_menu
USE symm_menu
USE shear
USE structur
USE spcgr_apply, ONLY: wyckoff_main
USE surface_func_mod
USE thermal_mod
USE transform_menu
USE waves_do_menu
USE discus_init_mod
USE discus_export
USE storage_menu_mod
!
use private_mod
!
USE blanks_mod
USE calc_expr_mod
USE doact_mod
USE errlist_mod 
use exit_para_mod
USE get_params_mod
USE class_macro_internal
USE kdo_all_mod
USE learn_mod 
USE lib_errlist_func
USE lib_do_operating_mod
USE lib_macro_func
USE precision_mod
USE prompt_mod
USE variable_mod
USE set_sub_generic_mod
USE str_comp_mod
!
IMPLICIT none 
!                                                                       
CHARACTER (LEN= * ), INTENT(INOUT) :: line 
LOGICAL            , INTENT(OUT)   :: lend
INTEGER            , INTENT(INOUT) :: length 
!                                                                       
CHARACTER(LEN=MAX(PREC_STRING, LEN(line))) :: zeile 
CHARACTER(len=14) :: befehl 
INTEGER :: indxb, indxg, lcomm, lbef 
INTEGER                  :: indxt ! position of a TAB
INTEGER                  ::  inverse_type
LOGICAL :: lout_rho, lkick 
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
      lbef = min (indxb - 1, len(befehl)) 
      befehl = line (1:lbef) 
!                                                                       
!------ command parameters start at the first character following       
!       the blank                                                       
!                                                                       
      zeile = ' ' 
      lcomm = 0 
      IF (indxb + 1.le.length) THEN 
         zeile = line (indxb + 1:length) 
         lcomm = length - indxb 
         CALL rem_leading_bl ( zeile, lcomm )
      ENDIF 
!                                                                       
!-------Suche nach einem "="                                            
!                                                                       
indxg = index (line, '=') 
IF(indxg /= 0.AND. .NOT. (str_comp (befehl, 'echo',   2, lbef, 4) )       &
             .AND. .NOT. (str_comp (befehl, 'system', 2, lbef,6) )       &
             .AND. .NOT. (str_comp (befehl, 'fput',   2, lbef, 4) )       &
             .AND. .NOT. (str_comp (befehl, 'help',   2, lbef, 4) .OR.    &
                          str_comp (befehl, '?   ',   2, lbef, 4) )       &
             .AND. INDEX(line,'==') == 0                            ) THEN      
!                                                                       
!             .AND. .NOT. (str_comp (befehl, 'socket', 2, lbef, 5) )     &
!-------Zuweisung eines Funktionswertes                                 
!                                                                       
         CALL do_math (line, indxg, length) 
      ELSE 
!                                                                       
!     --execute a macro file                                            
!                                                                       
         IF (befehl (1:1) .eq.'@') THEN 
            IF (length.ge.2) THEN 
               line(1:length-1) = line(2:length)
               length = length - 1
               CALL file_kdo(line, length)
            ELSE 
               ier_num = - 13 
               ier_typ = ER_MAC 
            ENDIF 
!
!     -- Allocate array sizes
!
         ELSEIF (str_comp (befehl, 'allocate', 3, lbef,  8) ) THEN
            CALL discus_do_allocate_appl (zeile, lcomm)
!
!     -- Deallocate array sizes
!
         ELSEIF (str_comp (befehl, 'deallocate', 3, lbef, 10) ) THEN
            CALL discus_do_deallocate_appl (zeile, lcomm)
!                                                                       
!-------add two files 'addf'                                            
!                                                                       
         ELSEIF (str_comp (befehl, 'addfile', 2, lbef, 7) ) THEN 
            CALL do_addfile (zeile, lcomm) 
!                                                                       
!     Set anisotropic APD 'anis'                                        
!                                                                       
         ELSEIF(str_comp(befehl, 'anis', 3, lbef, 4) ) THEN 
            CALL do_anis(zeile, lcomm) 
!                                                                       
!     append a new atom 'appe'                                          
!                                                                       
         ELSEIF (str_comp (befehl, 'append', 2, lbef, 6) ) THEN 
            lkick = .false. 
            CALL do_app (zeile, lcomm, lkick) 
!                                                                       
!-------show the asymmetric unit 'asym'                                 
!                                                                       
         ELSEIF (str_comp (befehl, 'asym', 2, lbef, 4) ) THEN 
            CALL show_asym 
!                                                                       
!-------do Bragg's equation 'bragg'
!                                                                       
         ELSEIF (str_comp (befehl, 'bragg', 2, lbef, 5) ) THEN 
            CALL bragg_eq(zeile, lcomm) 
!                                                                       
!-------remove atomes outside a boundary 'boundary'                     
!                                                                       
         ELSEIF (str_comp (befehl, 'boundary', 2, lbef, 8) ) THEN 
            CALL boundary (zeile, lcomm) 
!                                                                       
!-------change some properties in the crystal 'change'                  
!                                                                       
         ELSEIF (str_comp (befehl, 'change', 3, lbef, 6) ) THEN 
            CALL do_change (zeile, lcomm) 
!                                                                       
!-------show the atoms present in the crystal 'chem'                    
!                                                                       
         ELSEIF((linteractive.OR.lblock.OR.lmakro) .AND.  &
                str_comp (befehl, 'chemistry', 3, lbef, 8) ) THEN 
            CALL chem 
!                                                                       
!-------go to the connectivity menue  'connectivity'                    
!                                                                       
         ELSEIF ((linteractive.OR.lblock.OR.lmakro) .AND.  &
                str_comp (befehl, 'connectivity', 3, lbef, 12) ) THEN 
            CALL conn_menu 
!                                                                       
!-------Define some DISCUS parameters 'define'                          
!                                                                       
         ELSEIF (str_comp (befehl, 'define', 3, lbef, 6) ) THEN 
            CALL do_define (zeile, lcomm) 
!                                                                       
!-------Decorate a surface by molecules 'decorate'                      
!                                                                       
         ELSEIF ((linteractive.OR.lblock.OR.lmakro) .AND.  &
                str_comp (befehl, 'decorate', 3, lbef, 8) ) THEN 
            CALL do_place_molecule
!
!-------Demolecularize
!
        ELSEIF ((linteractive.OR.lblock.OR.lmakro) .AND.  &
                str_comp (befehl, 'demolecularize', 3, lbef, 14) ) THEN
           CALL demolecularize
!                                                                       
!-------Handling of domains within the host structure 'domain'          
!                                                                       
         ELSEIF ((linteractive.OR.lblock.OR.lmakro) .AND. str_comp (befehl, 'domain', 3, lbef, 6) ) THEN 
            CALL do_domain
!                                                                       
!-------copy an atom 'copy'                                             
!                                                                       
         ELSEIF (str_comp (befehl, 'copy', 3, lbef, 4) ) THEN 
            CALL do_copy (zeile, lcomm) 
!                                                                       
!-------Determine displacement of an atom 'displacement'                
!                                                                       
         ELSEIF (str_comp (befehl, 'displacent', 3, lbef, 10) ) THEN 
            CALL get_displacement (zeile, lcomm) 
!                                                                       
!     Transform vector from direct to reciprocal space                  
!                                                                       
         ELSEIF (str_comp (befehl, 'd2r', 2, lbef, 3) ) THEN 
            CALL d2r (zeile, lcomm, .true.) 
!                                                                       
!     Difference  Fourier 'diff'                                        
!                                                                       
         ELSEIF ((linteractive.OR.lblock.OR.lmakro) .AND. str_comp (befehl, 'diff', 2, lbef, 4) ) THEN 
            inverse_type = INV_DIFF 
            lout_rho = .true. 
            CALL patterson (inverse_type) 
!                                                                       
!-------Transform epsrimental powder pattern to PDF 'exp2pdf'
!                                                                       
         ELSEIF (str_comp (befehl, 'exp2pdf', 5, lbef, 7) ) THEN 
            CALL exp2pdf
!                                                                       
!-------export a file from discus format 'export'                  
!                                                                       
         ELSEIF (str_comp (befehl, 'export', 6, lbef, 6) ) THEN 
            CALL do_export (zeile, lcomm) 
!                                                                       
!-------Terminate DISCUS 'exit'                                         
!                                                                       
         ELSEIF (str_comp (befehl, 'exit', 2, lbef, 4) ) THEN 
            LEND = .TRUE. 
!                                                                       
!-------Find properties   'find'                                        
!                                                                       
         ELSEIF (str_comp (befehl, 'find', 3, lbef, 4) ) THEN 
            CALL do_find (zeile, lcomm) 
!                                                                       
!-------Fit something     'fit'                                        
!                                                                       
         ELSEIF (str_comp (befehl, 'fit', 3, lbef, 3) ) THEN 
            CALL discus_do_fit (zeile, lcomm) 
!                                                                       
!-------Fourier transform 'four'                                        
!                                                                       
         ELSEIF ((linteractive.OR.lblock.OR.lmakro) .AND. str_comp (befehl, 'fourier', 2, lbef, 7) ) THEN 
            lout_rho = .false. 
            CALL fourier 
!                                                                       
!-------import a file into discus.cell format 'import'                  
!                                                                       
         ELSEIF (str_comp (befehl, 'import', 3, lbef, 6) ) THEN 
            CALL do_import (zeile, lcomm) 
!                                                                       
!-------calculate cell and site <==> atomindex
!                                                                       
         ELSEIF (str_comp (befehl, 'index2cell', 3, lbef, 10) ) THEN 
            CALL chem_trans (zeile, lcomm) 
!                                                                       
!     inserte a new atom 'inse'                                         
!                                                                       
         ELSEIF (str_comp (befehl, 'insert', 3, lbef, 6) ) THEN 
            IF (str_comp (zeile, 'domain', 3, lcomm, 6) ) THEN 
               CALL insert ( - 1) 
            ELSEIF (str_comp (zeile, 'molecule', 3, lcomm, 8) ) THEN 
               CALL insert (0) 
            ELSEIF (str_comp (zeile, 'object', 3, lcomm, 6) ) THEN 
               CALL insert (1) 
            ELSE 
               CALL do_ins (zeile, lcomm) 
            ENDIF 
!                                                                       
!     interpret an electron density 'interpret'                         
!                                                                       
         ELSEIF (str_comp (befehl, 'interpret', 3, lbef, 9) ) THEN 
            CALL interpret 
!                                                                       
!     inverse Fourier 'inverse'                                         
!                                                                       
         ELSEIF ((linteractive.OR.lblock.OR.lmakro) .AND. str_comp (befehl, 'inverse', 3, lbef, 7) ) THEN 
            inverse_type = INV_INV 
            lout_rho = .true. 
            CALL patterson (inverse_type) 
!                                                                       
!     append a new atom while kicking out old atoms 'kick'              
!                                                                       
         ELSEIF (str_comp (befehl, 'kick', 2, lbef, 4) ) THEN 
            lkick = .true. 
            CALL do_app (zeile, lcomm, lkick) 
!                                                                       
!     group atoms into a molecule 'molecularize'
!                                                                       
         ELSEIF (str_comp (befehl, 'molecularize', 2, lbef, 12) ) THEN 
            CALL do_molecularize (zeile, lcomm)
!                                                                       
!     Output routines 'outp'                                            
!                                                                       
         ELSEIF ((linteractive.OR.lblock.OR.lmakro) .AND. str_comp (befehl, 'output', 1, lbef, 6) ) THEN 
            CALL do_niplps (lout_rho) 
!                                                                       
!     Patterson 'patterson'                                             
!                                                                       
         ELSEIF ((linteractive.OR.lblock.OR.lmakro) .AND. str_comp (befehl, 'patterson', 2, lbef, 9) ) THEN 
            inverse_type = INV_PATT 
            lout_rho = .true. 
            CALL patterson (inverse_type) 
!                                                                       
!     PDF calculation and analysis                                      
!                                                                       
         ELSEIF ((linteractive.OR.lblock.OR.lmakro) .AND. str_comp (befehl, 'pdf', 2, lbef, 3) ) THEN 
            CALL pdf 
!                                                                       
!     Create peridic crystal
!                                                                       
         ELSEIF ((linteractive.OR.lblock.OR.lmakro) .AND. str_comp (befehl, 'perioditize', 3, lbef, 11) ) THEN 
            CALL perioditize_menu 
!                                                                       
!     Plot the crystal 'plot'                                           
!                                                                       
         ELSEIF ((linteractive.OR.lblock.OR.lmakro) .AND. str_comp (befehl, 'plot', 2, lbef, 4) ) THEN 
            CALL plot 
!                                                                       
!     Calculate Powder pattern                                          
!                                                                       
         ELSEIF ((linteractive.OR.lblock.OR.lmakro) .AND. str_comp (befehl, 'powder', 2, lbef, 6) ) THEN 
            CALL do_powder 
!                                                                       
!     Project a vector onto another and onto a plane 'proj'             
!                                                                       
         ELSEIF (str_comp (befehl, 'project', 4, lbef, 7) ) THEN 
            CALL do_proj (zeile, lcomm) 
!                                                                       
!     Pivat user subroutine 'privat'
!                                                                       
         ELSEIF(str_comp(befehl, 'private', 6, lbef, 7) ) THEN 
            CALL do_private(zeile)
!                                                                       
!     Go to property menu 'property'                                    
!                                                                       
         ELSEIF ((linteractive.OR.lblock.OR.lmakro) .AND. str_comp (befehl, 'property', 4, lbef, 8) ) THEN 
            CALL property_menu 
!                                                                       
!     Purge the list of atoms 'purg'                                    
!                                                                       
         ELSEIF (str_comp (befehl, 'purge', 2, lbef, 5) ) THEN 
            CALL do_purge (zeile, lcomm)
!                                                                       
!     Transform vector from reciprocal to direct space                  
!                                                                       
         ELSEIF (str_comp (befehl, 'r2d', 2, lbef, 3) ) THEN 
            CALL d2r (zeile, lcomm, .false.) 
!                                                                       
!-------Einlesen einer Struktur/Zelle 'read'                            
!                                                                       
         ELSEIF ((linteractive.OR.lblock.OR.lmakro) .AND. str_comp (befehl, 'read', 3, lbef, 4) ) THEN 
            CALL read_struc 
!                                                                       
!     Remove a single atom 'remo'                                       
!                                                                       
         ELSEIF (str_comp (befehl, 'remove', 3, lbef, 5) ) THEN 
            CALL do_remove (zeile, lcomm) 
!                                                                       
!     Replace atom(s) or molecule(s)                                    
!                                                                       
         ELSEIF (str_comp (befehl, 'replace', 3, lbef, 7) ) THEN 
            CALL do_replace (zeile, lcomm) 
!                                                                       
!     Do Monte Carlo simulations 'rmc', 'mmc'
!                                                                       
         ELSEIF ((linteractive.OR.lblock.OR.lmakro) .AND. str_comp (befehl, 'rmc', 2, lbef, 3) ) THEN 
            CALL rmc 
         ELSEIF ((linteractive.OR.lblock.OR.lmakro) .AND. str_comp (befehl, 'mmc', 2, lbef, 3) ) THEN 
            CALL mmc 
!
!     reset discus to system start
!
         ELSEIF (str_comp(befehl, 'reset', 3, lbef, 5)) THEN
            CALL discus_initarrays
            CALL discus_reset_all
!                                                                       
!     save structure to file 'save'                                     
!                                                                       
         ELSEIF ((linteractive.OR.lblock.OR.lmakro) .AND. str_comp (befehl, 'save', 2, lbef, 4) ) THEN 
            CALL save_struc (zeile, lcomm) 
!                                                                       
!     generalized shear operation 'shear'                               
!                                                                       
         ELSEIF ((linteractive.OR.lblock.OR.lmakro) .AND. str_comp (befehl, 'shear', 3, lbef, 5) ) THEN 
            CALL shear_menue
!                                                                       
!     Show something                                 'show'             
!                                                                       
         ELSEIF (str_comp (befehl, 'show', 2, lbef, 4) ) THEN 
            CALL discus_do_show (zeile, lcomm) 
!                                                                       
!     set or show space group                                           
!                                                                       
         ELSEIF (str_comp (befehl, 'spacegroup', 2, lbef, 10) ) THEN 
            CALL set_spcgr(zeile, lcomm) 
!                                                                       
!     Go to stacking fault menu 'stack'                                 
!                                                                       
         ELSEIF ((linteractive.OR.lblock.OR.lmakro) .AND. str_comp (befehl, 'stack', 2, lbef, 5) ) THEN 
            CALL stack 
!
!     Go into storage menu 'storage'
!
         ELSEIF ((linteractive.OR.lblock.OR.lmakro) .AND. str_comp (befehl, 'storage', 3, lbef, 7) ) THEN
            CALL storage
!                                                                       
!     Go to surface menu 'surface'                                      
!                                                                       
         ELSEIF ((linteractive.OR.lblock.OR.lmakro) .AND. str_comp (befehl, 'surface', 2, lbef, 8) ) THEN 
            CALL surface_menu 
!                                                                       
!     switch two atoms 'swit'                                           
!                                                                       
         ELSEIF (str_comp (befehl, 'switch', 2, lbef, 6) ) THEN 
            CALL do_switch (zeile, lcomm) 
!                                                                       
!     generalized symmetry operation 'symm'                             
!                                                                       
         ELSEIF ((linteractive.OR.lblock.OR.lmakro) .AND. str_comp (befehl, 'symmetry', 3, lbef, 8) ) THEN 
            CALL symm 
!                                                                       
!     Thermal displacement of all atoms 'ther'                          
!                                                                       
         ELSEIF (str_comp (befehl, 'thermal', 2, lbef, 7) ) THEN 
            CALL ther_displ (zeile, lcomm) 
!                                                                       
!     unit cell transformations      'tran'                             
!                                                                       
         ELSEIF ((linteractive.OR.lblock.OR.lmakro) .AND. str_comp (befehl, 'transform', 2, lbef, 9) ) THEN 
            CALL transform 
!                                                                       
!     Three-D-PDF interpretation 'three'
!                                                                       
         ELSEIF (str_comp (befehl, 'three',2, lbef, 5) ) THEN 
            CALL three_main
!                                                                       
!       Vector Product 'vprod'                                          
!                                                                       
         ELSEIF (str_comp (befehl, 'vprod', 1, lbef, 5) ) THEN 
            CALL vprod (zeile, lcomm) 
!                                                                       
!------   Waves traveling through the crystal 'wave'                    
!                                                                       
         ELSEIF ((linteractive.OR.lblock.OR.lmakro) .AND. str_comp (befehl, 'waves', 3, lbef, 5) ) THEN 
            CALL waves_menu
!                                                                       
!       Determine Wyckoff symmetry                                      
!                                                                       
         ELSEIF (str_comp (befehl, 'wyckoff', 1, lbef, 7) ) THEN 
            CALL wyckoff_main (zeile, lcomm) 
!                                                                       
!       Branch to KUPLOT (standalone call system, suite do branch)
!                                                                       
         ELSEIF ((linteractive.OR.lblock.OR.lmakro) .AND. str_comp (befehl, 'branch', 2, lbef, 6) ) THEN 
            CALL p_branch (zeile, lcomm, .FALSE., 0     ) 
!                                                                       
!------   Try general commands                                          
!                                                                       
         ELSE 
            CALL kdo_all (befehl, lbef, zeile, lcomm) 
            IF(zeile == 'EXIT') THEN ! kdo_all detected "continue suite"
               lend = .TRUE. 
            ENDIF
         ENDIF 
      ENDIF 
!
if(ex_do_exit) lend = .true.   ! A global exit was flagged
!                                                                       
END SUBROUTINE discus_mache_kdo                      
!
!*****7**************************************************************** 
!
      SUBROUTINE do_define (zeile, lp) 
!-                                                                      
!     Sets the value of status variables, DISCUS specific routine       
!+                                                                      
      USE unitcell_mod 
      USE errlist_mod 
      USE get_params_mod
USE precision_mod
USE str_comp_mod
      IMPLICIT none 
!                                                                       
!                                                                       
      INTEGER maxw 
      PARAMETER (maxw = 5) 
!                                                                       
      CHARACTER(LEN=*), INTENT(INOUT) ::  zeile 
      CHARACTER(LEN=MAX(PREC_STRING, LEN(zeile))) :: cpara (maxw) 
      INTEGER lpara (maxw), lp 
      INTEGER ianz 
!                                                                       
!                                                                       
      IF (zeile.ne.' ') THEN 
         CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
         IF (ier_num.eq.0) THEN 
!                                                                       
!----- ---- define error                                                
!                                                                       
            IF (str_comp (cpara (1) , 'generator', 1, lpara (1) , 9) ) THEN 
               IF (ianz.eq.2) THEN 
                  IF (str_comp (cpara (2) , 'center', 2, lpara (2) , 5) &
                  ) THEN                                                
                     gen_sta = GEN_CENTER 
                  ELSEIF (str_comp (cpara (2) , 'symmetry', 2, lpara (2)&
                  , 8) ) THEN                                           
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
subroutine discus_do_fit(line, length)
use crystal_mod
use fit_mod
use precision_mod
!
character(len=*), INTENT(INOUT) :: Line
integer         , INTENT(INOUT) :: length
!
INTEGER, DIMENSION(:), ALLOCATABLE :: list
INTEGER               :: i
REAL(kind=PREC_DP)   , DIMENSION(3) :: hkl
REAL(kind=PREC_DP)                  :: dist
!
ALLOCATE(list(1:cr_natoms))
         DO i=1,cr_natoms
            list(i) = i
         ENDDO
         CALL dis_fit_plane(cr_natoms, list, hkl, dist)
         DEALLOCATE(list)
!
end subroutine discus_do_fit
!
!end module discus_mache_kdo_mod
