!*****7*****************************************************************
!     This routine interprets the commands and executes the             
!     corresponding function.                                           
!*****7*****************************************************************
recursive SUBROUTINE kuplot_mache_kdo (line, lend, length) !, previous) 
!                                                                       
!     Main menu for KUPLOT                                              
!                                                                       
USE nexus_kuplot
USE kuplot_diffev_mod
USE kuplot_2dm
!
!USE kuplot_fit_old_mod
use kuplot_adt_mod
use kuplot_blk_mod
use kuplot_calc_mod
use kuplot_color_mod
use kuplot_exit_mod
USE kuplot_fit6
use kuplot_frame_mod
use kuplot_load_mod
use kuplot_math_mod
use kuplot_para_mod
use kuplot_plot_mod
use kuplot_reset_mod
use kuplot_save_mod
use kuplot_show_mod
USE kuplot_toglobal
!
USE ber_params_mod
USE blanks_mod
USE calc_expr_mod
USE doact_mod
USE errlist_mod 
use exit_para_mod
USE get_params_mod
USE class_macro_internal
USE kdo_all_mod
USE learn_mod 
USE prompt_mod 
USE param_mod 
!                                                                       
USE kuplot_config 
USE kuplot_mod 
!     USE kuplot_fit_mod
USE lib_errlist_func
USE lib_macro_func
USE precision_mod
USE set_sub_generic_mod
USE str_comp_mod
use take_param_mod
USE variable_mod
!
      IMPLICIT none 
!                                                                       
      CHARACTER (LEN= * ), INTENT(INOUT) :: line 
      LOGICAL            , INTENT(OUT)   :: lend
      INTEGER            , INTENT(INOUT) :: length 
!     CHARACTER (LEN= * ), DIMENSION(3), INTENT(INOUT) :: previous
!                                                                       
integer, parameter :: maxw  = 2
!      PARAMETER (maxw = 2) 
!                                                                       
CHARACTER(LEN=PREC_STRING) :: zei
CHARACTER(LEN=PREC_STRING) :: cpara (maxw) 
CHARACTER(LEN=11)          :: bef 
      REAL(KIND=PREC_DP), DIMENSION(MAXW) :: werte
      REAL dummy 
      INTEGER lpara (maxw) 
      INTEGER lc, lbef, i, j
      INTEGER ianz, indxg, indxb , indxt
      LOGICAL ldummy 
!                                                                       
integer, parameter :: NOPTIONAL = 1
integer, parameter :: O_DATA    = 1
character(len=4), dimension(NOPTIONAL) :: oname
character(len=MAX(PREC_STRING,len(line))), dimension(NOPTIONAL) :: opara   !Optional parameter strings returned
integer            , dimension(NOPTIONAL) :: loname  !Lenght opt. para name
integer            , dimension(NOPTIONAL) :: lopara  !Lenght opt. para name returned
logical            , dimension(NOPTIONAL) :: lpresent!opt. para is present
real(kind=PREC_DP) , dimension(NOPTIONAL) :: owerte   ! Calculated values
integer, parameter                        :: ncalc = 1 ! Number of values to calculate
!
data  oname / 'data' /
data loname / 1      /
!
!REFDIF
opara  = (/'1'    /)
lopara = (/ 1     /)
owerte = (/ 1.0D0 /)

!                                                                       
      CALL no_error 
!                                                                       
!------ If a commentary return immediately                              
!                                                                       
      IF (line (1:1) .eq.' '.or.line (1:1) .eq.'#'.or.line (1:1) .eq.'!'&
     &.or.length.eq.0) return                                           
!                                                                       
!     Only the first 4 characters are significant. The command consists 
!     of the four nonblank characters                                   
!                                                                       
      bef = '    ' 
      indxt = INDEX (line, achar(9))       ! find a tabulator
      IF(indxt==0) indxt = length + 1
      indxb = index (line, ' ')       ! find a blank
      IF(indxb==0) indxb = length + 1
      indxb = MIN(indxb,indxt)
      lbef = min (indxb - 1, 4)
      bef = line (1:lbef) 
!                                                                       
!     command parameters start at the first character                   
!     following the blank                                               
!                                                                       
      zei = ' ' 
      lc = 0 
      IF (indxb + 1.le.length) then 
         zei = line (indxb + 1:length) 
         lc = length - indxb 
         call rem_leading_bl(zei, lc)
      ENDIF 
!                                                                       
!-------Suche nach einem "="                                            
!                                                                       
      indxg = index (line, '=') 
      IF (indxg.ne.0  &
         .and..not. (str_comp (bef, 'echo', 2, lbef, 4) )   &
         .and..not. (str_comp (bef, 'system', 2, lbef, 6) )   &
         .and..not. (str_comp (bef, 'achx', 2, lbef, 4) )   &
         .and..not. (str_comp (bef, 'achy', 2, lbef, 4) )   &
         .and..not. (str_comp (bef, 'tit1', 2, lbef, 4) )   &
         .and..not. (str_comp (bef, 'tit2', 2, lbef, 4) )   &
         .and..not. (str_comp (bef, 'sleg', 2, lbef, 4) )   &
         .and..not. (str_comp (bef, 'sann', 2, lbef, 4) )   &
         .and..not. (str_comp (bef, 'fput', 2, lbef, 4) )   &
         .and..not. (str_comp (bef, 'help', 2, lbef, 4) .or.&
                     str_comp (bef, '?   ', 2, lbef, 4) )   &
         .AND. INDEX(line,'==') == 0 ) then
!                                                                       
!-------Zuweisung eines Funktionswertes                                 
!                                                                       
         CALL do_math (line, indxg, length) 
      ELSE 
!                                                                       
!     --execute a macro file                                            
!                                                                       
         IF (bef (1:1) .eq.'@') then 
            IF (length.ge.2) then 
               line =  line(2:length)
               length = length - 1 
               CALL file_kdo (line, length)
            ELSE 
               ier_num = - 13 
               ier_typ = ER_MAC 
            ENDIF 
!
!------ FINISHED AN LOW level Macro return to specified place
!
!        ELSEIF(str_comp(bef, 'finished', 8, lbef, 8) ) THEN
!           IF(str_comp (previous(1), 'fit', 3, LEN_TRIM(previous), 3) ) THEN 
!              zei = bef            ! flag down the 'finished' command
!              lc  = LEN_TRIM(bef)
!           write(*,*) ' RETURNED TO FIT  with ', zei(1:lc)
!              CALL do_fit (zei, lc, previous) 
!           ENDIF
!           previous(:) = ' '
!                                                                       
!------ Terminate KUPLOT 'exit'                                         
!                                                                       
         ELSEIF (str_comp (bef, 'exit', 2, lbef, 4) ) then 
            call kuplot_exit
            lend = .true. 
         ELSEIF (str_comp (bef, 'finished', 2, lbef, 8) ) then 
            call kuplot_exit
            lend = .true. 
!                                                                       
!-------  Start of KUPLOT specific commands                             
!-------  Allocate memory for new dataset                               
!                                                                       
         elseif(str_comp(bef, 'adt', 3, lbef, 3)) then
            call adt_menu()
         ELSEIF (str_comp (bef, 'allocate', 2, lbef, 8) ) then 
            CALL do_allocate (zei, lc, .TRUE.) 
!                                                                       
!-------  Set bond drawing options                                      
!                                                                       
         ELSEIF (str_comp (bef, 'bond', 2, lbef, 4) ) then 
            CALL set_bond (zei, lc) 
!                                                                       
!-------  Perform calculations, merge and rebin                         
!                                                                       
         elseif(str_comp (bef, 'mergeg',6, lbef, 6)) then
            call do_merge (zei, lc, .true.)
         ELSEIF (str_comp (bef, 'ccalc', 3, lbef, 5) ) then 
            CALL do_calc (zei, lc) 
         ELSEIF (str_comp (bef, 'convolute', 3, lbef, 9) ) then 
            CALL do_convolute (zei, lc) 
         ELSEIF (str_comp (bef, 'exclude', 3, lbef, 7) ) then 
            CALL do_exclude (zei, lc) 
         ELSEIF (str_comp (bef, 'kcalc', 3, lbef, 5) ) then 
            CALL do_kmath (zei, lc) 
         ELSEIF (str_comp (bef, 'merge', 3, lbef, 5) ) then 
            CALL do_merge (zei, lc, .false.) 
         ELSEIF (str_comp (bef, 'rebin', 3, lbef, 5) ) then 
            CALL do_rebin (zei, lc) 
         ELSEIF(str_comp (bef, 'close', 3, lbef, 5)) THEN
            CALL get_params (zei, ianz, cpara, lpara, maxw, lc) 
            IF(ianz>= 1) THEN
               CALL ber_params(ianz, cpara, lpara, werte, maxw)
               IF(ier_num == 0 ) THEN
                  DO j=1, ianz
                     i = NINT(werte(j))
                     IF(j<=MAXWIN) THEN
                        IF(dev_id(i,x11)>0) THEN
                           CALL PGSLCT (dev_id (i, x11) )
                           CALL PGCLOS
                           dev_id(i,x11) = -1
                        ENDIF
                     ELSE
                        ier_num = -9999
                        ier_typ = ER_APPL
                     ENDIF
                  ENDDO
               ENDIF
            ELSE
               DO i=1, MAXWIN
                  IF(dev_id(i,x11)>0) THEN
                     CALL PGSLCT (dev_id (i, x11) )
                     CALL PGCLOS
                     dev_id(i,x11) = -1
                  ENDIF
               ENDDO
            ENDIF
!                                                                       
!-------  Set colours                                                   
!                                                                       
         ELSEIF (str_comp (bef, 'color', 3, lbef, 5) ) then 
            CALL set_color (zei, lc) 
!                                                                       
!------- define a cost function value to be returned to diffev
!                                                                       
         ELSEIF (str_comp (bef, 'costvalue', 3, lbef, 9) ) then 
            CALL get_params (zei, ianz, cpara, lpara, maxw, lc) 
            IF(ianz>= 1) then
               CALL ber_params(ianz, cpara, lpara, werte, maxw)
               IF(ier_num == 0 ) THEN
                  DO i=0, ianz-1
                     rvalues(2, i) = werte(1+i)
                  ENDDO
                  IF(ianz>1) THEN
                     nrvalues = MAX(nrvalues, ianz-1)
                  ELSEIF(ianz==1) THEN
                     nrvalues = MAX(nrvalues, 1)
                  ENDIF
                  rvalue_yes = .true.
               ENDIF
            ENDIF
!                                                                       
!-------  Set bitmap colormap                                           
!                                                                       
         ELSEIF (str_comp (bef, 'cmap', 3, lbef, 4) ) then 
            CALL set_cmap (zei, lc) 
!                                                                       
!-------  Calculate derivative                                          
!                                                                       
         ELSEIF (str_comp (bef, 'derivative', 3, lbef, 10) ) then 
            CALL do_derivative (zei, lc) 
!
!-------  Make kpara or kpar_par plot
!
         ELSEIF (str_comp (bef, 'kpara', 3, lbef, 6) ) then 
            CALL do_diffev_plot (zei, lc) 
!                                                                       
!-------  Set fill parameters                                           
!                                                                       
         ELSEIF (str_comp (bef, 'fill', 3, lbef, 4) ) then 
            CALL set_fill (zei, lc) 
!                                                                       
!-------  Enter fit sublevel                                            
!                                                                       
!        ELSEIF ((linteractive.OR.lblock.OR.lmakro) .AND. &
!                str_comp (bef, 'f55', 3, lbef, 3)       ) THEN
!           CALL do_fit (zei, lc) !, previous) 
         ELSEIF ((linteractive.OR.lblock.OR.lmakro) .AND. &
                 str_comp (bef, 'fit', 3, lbef, 3)       ) THEN
            CALL do_f66 (zei, lc) !, previous) 
!           CALL do_fit (zei, lc) !, previous) 
!           previous(1) = 'fit 1'
!           previous(2) = 'run'
!        ELSEIF ((linteractive.OR.lblock.OR.lmakro) .AND. &
!                str_comp (bef, 'f66', 3, lbef, 3)       ) THEN
!           CALL do_f66 (zei, lc) !, previous) 
!                                                                       
!-------  Plot filenames on plot                                        
!                                                                       
         ELSEIF (str_comp (bef, 'fname', 3, lbef, 5) ) then 
            CALL get_params (zei, ianz, cpara, lpara, maxw, lc) 
            IF (ier_num.ne.0) return 
            ifname (iwin, iframe) = str_comp (cpara (1) , 'on', 2,      &
            lpara (1) , 2)                                              
!                                                                       
!-------  Set typ of box & axis                                         
!                                                                       
         ELSEIF (str_comp (bef, 'fset', 3, lbef, 4) ) then 
            CALL set_ibox (zei, lc) 
!                                                                       
!-------  Set font size, type and color                                 
!                                                                       
         ELSEIF (str_comp (bef, 'font', 3, lbef, 4) ) then 
            CALL set_font (zei, lc) 
!                                                                       
!-------  Plot frame around plot/frames                                 
!                                                                       
         ELSEIF (str_comp (bef, 'frame', 3, lbef, 5) ) then 
            CALL get_params (zei, ianz, cpara, lpara, maxw, lc) 
            IF (ier_num.ne.0) return 
            tot_frame (iwin) = str_comp (cpara (1) , 'on', 2, lpara (1) &
            , 2)                                                        
!                                                                       
!-------  Calculate Fourier Transform from data set                     
!                                                                       
         ELSEIF (str_comp (bef, 'fourier', 2, lbef, 7) ) then 
            CALL do_four (zei, lc) 
         ELSEIF (str_comp (bef, 'fft' , 2, lbef, 3) ) then 
            CALL kuplot_do_fft (zei, lc) 
!                                                                       
!-------  Create data set from function                                 
!                                                                       
         ELSEIF (str_comp (bef, 'function', 2, lbef, 8) ) then 
            CALL do_func (zei, lc) 
!                                                                       
!-------  smooth data 'glat' and 'smooth'                               
!                                                                       
         ELSEIF (str_comp (bef, 'glat', 2, lbef, 4) ) then 
            CALL do_glat (zei, lc, .false.) 
         ELSEIF (str_comp (bef, 'smooth', 3, lbef, 6) ) then 
            CALL do_glat (zei, lc, .true.) 
!                                                                       
!-------  grid on/off                                                   
!                                                                       
         ELSEIF (str_comp (bef, 'grid', 3, lbef, 4) ) then 
            CALL get_params (zei, ianz, cpara, lpara, maxw, lc) 
            IF (ier_num.ne.0) return 
            igrid (iwin, iframe) = str_comp (cpara (1) , 'on', 2, lpara &
            (1) , 2)                                                    
!                                                                       
!-------  ident on/off                                                  
!                                                                       
         ELSEIF (str_comp (bef, 'ident', 3, lbef, 5) ) then 
            CALL get_params (zei, ianz, cpara, lpara, maxw, lc) 
            IF (ier_num.ne.0) return 
            iden (iwin) = str_comp (cpara (1) , 'on', 2, lpara (1) , 2) 
!                                                                       
!-------  set contour line parameters                                   
!                                                                       
         ELSEIF (str_comp (bef, 'hline', 3, lbef, 5) ) then 
            CALL set_hlin (zei, lc) 
         ELSEIF (str_comp (bef, 'hcolor', 3, lbef, 6) ) then 
            CALL para_setii (zei, lc, hlinecol, maxkurvtot, maxhl, bef, &
            1, 15)                                                      
         ELSEIF (str_comp (bef, 'htype', 3, lbef, 5) ) then 
            CALL para_setii (zei, lc, hlinetyp, maxkurvtot, maxhl, bef, &
            0, 5)                                                       
         ELSEIF (str_comp (bef, 'hpak', 3, lbef, 4)    .or. &
                 str_comp (bef, 'hpackage', 3, lbef, 8) ) then 
            CALL set_hpak (zei, lc) 
!                                                                       
!-------  calculate integral                                            
!                                                                       
         ELSEIF (str_comp (bef, 'integral', 2, lbef, 8) ) then 
            CALL do_inte (zei, lc) 
!                                                                       
!-------  Interpolate function on grid of other function                
!                                                                       
         ELSEIF (str_comp (bef, 'spline', 3, lbef, 6) ) then 
            CALL do_interpolate (zei, lc) 
!                                                                       
!-------  save data set                                                 
!                                                                       
         ELSEIF (str_comp (bef, 'ksave', 3, lbef, 5) ) then 
            CALL do_ksav (zei, lc) 
         ELSEIF (str_comp (bef, 'dsave', 3, lbef, 5) ) then 
            CALL do_save_data (zei, lc) 
!                                                                       
!-------  Load files 'load'                                             
!                                                                       
         ELSEIF (str_comp (bef, 'load', 2, lbef, 4) ) then 
            CALL do_load (zei, lc, .TRUE.) 
!                                                                       
!-------  Load 2d files into a map '2dm'                                             
!                                                                       
         ELSEIF (str_comp (bef, '2dm', 2, lbef, 3) ) then 
            CALL kuplot_2dm_menu
!                                                                       
!-------  Set line width                                                
!                                                                       
         ELSEIF (str_comp (bef, 'lwidth', 2, lbef, 6) ) then 
            CALL set_linewidth (zei, lc) 
!                                                                       
!-------  Set tick marks                                                
!                                                                       
         ELSEIF (str_comp (bef, 'mark', 3, lbef, 4) ) then 
            CALL set_mark (zei, lc) 
!                                                                       
!-------  Match data sets                                               
!                                                                       
         ELSEIF (str_comp (bef, 'match', 3, lbef, 5) ) then 
            CALL do_match (zei, lc) 
!                                                                       
!-------  Activate mouse menu                                           
!                                                                       
!         ELSEIF ((linteractive.OR.lblock.OR.lmakro) .AND. str_comp (bef, 'mouse', 3, lbef, 5) ) then 
!            CALL do_mouse (zei, lc) 
!                                                                       
!-------  Calculate averages, ..                                        
!                                                                       
         ELSEIF (str_comp (bef, 'mean', 3, lbef, 4) ) then 
            CALL do_mean (zei, lc, .true.) 
!                                                                       
!-------  NeXus commands                                                
!                                                                       
         ELSEIF (str_comp (bef, 'nxopen', 3, lbef, 6) ) then 
            CALL do_nxopen (zei, lc) 
         ELSEIF (str_comp (bef, 'nxclose', 3, lbef, 7) ) then 
            CALL do_nxclose 
         ELSEIF (str_comp (bef, 'nxdir', 3, lbef, 5) ) then 
            CALL do_nxdir 
         ELSEIF (str_comp (bef, 'nxload', 3, lbef, 6) ) then 
            CALL do_nxload (zei, lc) 
!                                                                       
!-------  Set marker size                                               
!                                                                       
         ELSEIF (str_comp (bef, 'msize', 2, lbef, 5) ) then 
            CALL set_sizemark (zei, lc) 
!                                                                       
!-------  set landscape / portrait orientation                          
!                                                                       
         ELSEIF (str_comp (bef, 'orientation', 3, lbef, 11) ) then 
            CALL get_params (zei, ianz, cpara, lpara, maxw, lc) 
            IF (ier_num.ne.0) return 
            ldummy = str_comp(cpara(1), 'landscape', 1, lpara(1), 9)
            IF (ldummy.neqv.orient (iwin) ) then 
               orient (iwin) = ldummy 
               dummy = dev_width (iwin) 
               dev_width (iwin) = dev_height (iwin) 
               dev_height (iwin) = dummy 
            ENDIF 
!                                                                       
!-------  Calculate the residual of two data sets                       
!                                                                       
         ELSEIF (str_comp (bef, 'rvalue', 3, lbef, 6) ) then 
            CALL do_rvalue (zei, lc, .TRUE.) 
!                                                                       
!-------  Setting plotting window                                       
!                                                                       
         ELSEIF (str_comp (bef, 'window', 2, lbef, 6) ) then 
            CALL set_window (zei, lc) 
!                                                                       
!-------  Here are all frame related commands: nfra,afra,kfra,...       
!                                                                       
         ELSEIF (str_comp (bef, 'nframe', 2, lbef, 6) ) then 
            CALL set_nfra (zei, lc) 
         ELSEIF (str_comp (bef, 'kframe', 2, lbef, 6) ) then 
            CALL set_kfra (zei, lc) 
         ELSEIF (str_comp (bef, 'sframe', 2, lbef, 6) ) then 
            CALL set_sfra (zei, lc) 
         ELSEIF (str_comp (bef, 'aframe', 2, lbef, 6) ) then 
            CALL set_afra (zei, lc) 
         ELSEIF (str_comp (bef, 'cframe', 2, lbef, 6) ) then 
            CALL set_cfra (zei, lc) 
         ELSEIF (str_comp (bef, 'bframe', 2, lbef, 6) ) then 
            CALL set_bfra (zei, lc) 
         ELSEIF (str_comp (bef, 'mass', 3, lbef, 4) ) then 
            WRITE (output_io, 2000) 
!                                                                       
!-------  Plot graphic                                                  
!                                                                       
         ELSEIF (str_comp (bef, 'plot', 2, lbef, 4) ) then 
            if(l_plot_status) CALL do_plot (.false.) 
!                                                                       
!-------  Save or print graphic                                         
!                                                                       
         ELSEIF (str_comp (bef, 'save', 2, lbef, 4) .or. &
                 str_comp (bef, 'print', 2, lbef, 5) ) then                                     
            CALL do_hardcopy (bef, zei, lbef, lc) 
!                                                                       
!-------  Commands 'achx','achy','tit1' and 'tit2'                      
!                                                                       
         ELSEIF (str_comp (bef, 'achx', 4, lbef, 4) ) then 
            CALL para_set_achse (zei, lc, achse (iwin, iframe, 1),      &
            lachse (iwin, iframe, 1) )                                  
         ELSEIF (str_comp (bef, 'achy', 4, lbef, 4) ) then 
            CALL para_set_achse (zei, lc, achse (iwin, iframe, 2),      &
            lachse (iwin, iframe, 2) )                                  
         ELSEIF (str_comp (bef, 'achz', 4, lbef, 4) ) then 
            CALL para_set_achse (zei, lc, achse (iwin, iframe, 3),      &
            lachse (iwin, iframe, 3) )                                  
!                                                                       
         ELSEIF (str_comp (bef, 'tit1', 4, lbef, 4) ) then 
            CALL para_set_title (zei, lc, titel (iwin, iframe, 1) ) 
         ELSEIF (str_comp (bef, 'tit2', 4, lbef, 4) ) then 
            CALL para_set_title (zei, lc, titel (iwin, iframe, 2) ) 
!
!
         ELSEIF (str_comp (bef, 'torefine', 4, lbef, 8)) THEN
            CALL kuplot_to_global(zei)
!                                                                       
!-------  Commands 'mass','aver',...                                    
!                                                                       
         ELSEIF (str_comp (bef, 'aver', 3, lbef, 4) ) then 
            call get_params (zei, ianz, cpara, lpara, maxw, lc) 
            call get_optional(ianz, MAXW, cpara, lpara, NOPTIONAL,  ncalc, &
                              oname, loname, opara, lopara, lpresent, owerte)
            if(lpresent(O_DATA)) then       ! Optional parameter data: 
               j = nint(owerte(O_DATA))
               call get_aver_angle(MAXW, j, werte)
               yskal_u(iwin, iframe) = real(werte(1), kind=PREC_SP)
            else
               CALL para_setr(zei, lc, yskal_u(iwin, iframe), bef, -9999.)
            endif
            lyskal(iwin, iframe) = (yskal_u(iwin, iframe) /=  -9999.)
         ELSEIF (str_comp (bef, 'angle', 3, lbef, 5) ) then 
            call get_params (zei, ianz, cpara, lpara, maxw, lc) 
            call get_optional(ianz, MAXW, cpara, lpara, NOPTIONAL,  ncalc, &
                              oname, loname, opara, lopara, lpresent, owerte)
            if(lpresent(O_DATA)) then       ! Optional parameter data: 
               j = nint(owerte(O_DATA))
               call get_aver_angle(MAXW, j, werte)
               shear(iwin, iframe) = real(werte(2), kind=PREC_SP)
            else
               CALL para_setr (zei, lc, shear (iwin, iframe), bef, 90.0) 
            endif
!                                                                       
!-------  Commands 'ltyp','lcol','mtyp','mcol',...                      
!                                                                       
         ELSEIF (str_comp (bef, 'ltype', 3, lbef, 5) ) then 
            CALL para_seti (zei, lc, ilinetyp, 1, maxkurvtot, bef,-5, 5,&
            .false.)                                                    
         ELSEIF (str_comp (bef, 'mtype', 3, lbef, 5) ) then 
            CALL para_seti (zei, lc, imarktyp, 1, maxkurvtot, bef,      &
            - 3, 5000, .false.)                                         
         ELSEIF (str_comp (bef, 'lcolor', 3, lbef, 6) ) then 
            CALL para_seti (zei, lc, ilinecol, 0, maxkurvtot, bef,-1,   &
            15, .true.)                                                 
         ELSEIF (str_comp (bef, 'mcolor', 3, lbef, 6) ) then 
            CALL para_seti (zei, lc, imarkcol, 1, maxkurvtot, bef, 1,   &
            15, .false.)                                                
         ELSEIF (str_comp (bef, 'ecolor', 3, lbef, 6) ) then 
            CALL para_seti (zei, lc, ierrcol, 1, maxkurvtot, bef, 1, 15,&
            .false.)                                                    
         ELSEIF (str_comp (bef, 'lart', 3, lbef, 4) ) then 
            CALL para_seti (zei, lc, ilineart, 1, maxkurvtot, bef, 1, 3,&
            .false.)                                                    
         ELSEIF (str_comp (bef, 'hart', 3, lbef, 4) ) then 
            CALL para_seti (zei, lc, hlineart, 1, maxkurvtot, bef, 1, 4,&
            .false.)                                                    
         ELSEIF (str_comp (bef, 'hlab', 3, lbef, 4) ) then 
            CALL para_seti (zei, lc, hlabel, 1, maxkurvtot, bef, 0, 99, &
            .false.)                                                    
         ELSEIF (str_comp (bef, 'ptype', 3, lbef, 5) ) then 
            CALL para_seti (zei, lc, imarkmax, 1, maxkurvtot, bef,      &
            - 3, 5000, .false.)                                         
         ELSEIF (str_comp (bef, 'etype', 3, lbef, 5) ) then 
            CALL para_seti (zei, lc, ierr, 1, maxkurvtot, bef, 0, 3,    &
            .false.)                                                    
!                                                                       
!-------  Reset                                                         
!                                                                       
         ELSEIF (str_comp (bef, 'reset', 3, lbef, 5) ) then 
            CALL kuplot_do_reset (zei, lc) 
!                                                                       
!-------  Set annotations                                               
!                                                                       
         ELSEIF (str_comp (bef, 'sann', 3, lbef, 4) ) then 
            CALL set_annotation (zei, lc) 
!                                                                       
!-------  Read/Write defaults                                           
!                                                                       
         ELSEIF (str_comp (bef, 'sdef', 3, lbef, 4) ) then 
            CALL do_write_def (zei, lc) 
         ELSEIF (str_comp (bef, 'rdef', 3, lbef, 4) ) then 
            CALL do_read_def (zei, lc) 
!                                                                       
!-------  Show parameters 'show'                                        
!                                                                       
         ELSEIF (str_comp (bef, 'show', 3, lbef, 4) ) then 
            CALL kuplot_do_show (zei, lc) 
!                                                                       
!-------  Set data set caption                                          
!                                                                       
         ELSEIF (str_comp (bef, 'sleg', 3, lbef, 4) ) then 
            CALL set_legend (zei, lc) 
!                                                                       
!-------  Set buffer around plotting window                             
!                                                                       
         ELSEIF (str_comp (bef, 'buffer', 3, lbef, 6) ) then 
            CALL set_buff (zei, lc) 
!                                                                       
!-------  Set plotting window                                           
!                                                                       
         ELSEIF (str_comp (bef, 'skal', 3, lbef, 4)   .or. &
                 str_comp (bef, 'scale', 3, lbef, 5) ) then 
            CALL set_skal (zei, lc) 
!                                                                       
!-------  Solve zero points in a polynomial                             
!                                                                       
         ELSEIF (str_comp (bef, 'solve', 3, lbef, 5) ) then 
            CALL do_solve(zei, lc) 
!                                                                       
!-------  Sort data set (up in x)                                       
!                                                                       
         ELSEIF (str_comp (bef, 'sort', 3, lbef, 4) ) then 
            CALL do_sort (zei, lc) 
!                                                                       
!-------  Search maxima/minima                                          
!                                                                       
         ELSEIF (str_comp (bef, 'smax', 3, lbef, 4) ) then 
            CALL do_smax (zei, lc) 
!                                                                       
!       Branch to DISCUS (standalone call system, suite do branch)
!                                                                       
         ELSEIF ((linteractive.OR.lblock.OR.lmakro) .AND. &
            str_comp (bef, 'branch', 2, lbef, 6) ) then 
            CALL p_branch (zei, lc, .FALSE., 0     ) 
!                                                                       
!-------  Check for generic command                                     
!                                                                       
         ELSE 
            CALL kdo_all (bef, lbef, zei, lc) 
            IF(zei == 'EXIT') THEN ! kdo_all detected "continue suite"
               lend = .TRUE.
            ENDIF 
         ENDIF 
      ENDIF 
!
if(ex_do_exit) lend = .true.   ! A global exit was flagged
!                                                                       
 2000 FORMAT     (' ****WARN**** Redundant command, do not use',        &
     &          13x,'**** -99 ****')                                    
      END SUBROUTINE kuplot_mache_kdo                      
