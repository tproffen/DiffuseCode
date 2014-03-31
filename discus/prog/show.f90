MODULE show_menu
!
CONTAINS
!*****7*****************************************************************
!                                                                       
      SUBROUTINE do_show (line, laenge) 
!-                                                                      
!     These subroutine is the main routine for showing various          
!     parameters.                                                       
!+                                                                      
      USE config_mod 
      USE allocate_appl_mod
      USE crystal_mod 
      USE conn_mod
!                                                                       
      USE errlist_mod 
      USE param_mod 
      USE prompt_mod 
      IMPLICIT none 
       
!                                                                       
      INTEGER maxw 
      PARAMETER (maxw = 7) 
!                                                                       
      INTEGER mode 
!                                                                       
      INTEGER FULL, SYMBOL, XYZ, MATRIX 
      PARAMETER (FULL = 0) 
      PARAMETER (SYMBOL = 1) 
      PARAMETER (XYZ = 2) 
      PARAMETER (MATRIX = 3) 
!                                                                       
      CHARACTER ( * ) line 
      CHARACTER(1024) cpara (maxw) 
      INTEGER lpara (maxw) 
      INTEGER ianz , iianz
      INTEGER i, j 
      INTEGER laenge 
      REAL werte (maxw) 
      CHARACTER (LEN=256)  :: c_name   ! Connectivity name
      INTEGER              :: c_name_l ! connectivity name length
      INTEGER              :: is1, ino ! Connectivity id, connectivity no
      INTEGER              :: iatom    ! atoms no for show
      LOGICAL              :: long
!                                                                       
      LOGICAL str_comp 
!                                                                       
      CALL get_params (line, ianz, cpara, lpara, maxw, laenge) 
      IF (ier_num.eq.0) then 
!                                                                       
!     --Interprete first parameter as command                           
!                                                                       
!                                                                       
!     ----Show composition of asymmetric unit 'asym'                    
!                                                                       
         IF (str_comp (cpara (1) , 'asym', 2, lpara (1) , 4) ) then 
            CALL show_asym 
!                                                                       
!     ----Show an atom                     'atom'                       
!                                                                       
         ELSEIF (str_comp (cpara (1) , 'atom', 2, lpara (1) , 4) ) then 
            CALL do_show_atom (ianz, cpara, lpara, werte, maxw) 
!                                                                       
!     ----Show bond valence parameters     'bval'                       
!                                                                       
         ELSEIF (str_comp (cpara (1) , 'bval', 2, lpara (1) , 4) ) then 
            CALL do_show_bval (ianz, cpara, maxw) 
!                                                                       
!     ----Show the chemistry               'chem'                       
!                                                                       
         ELSEIF (str_comp (cpara (1) , 'chem', 2, lpara (1) , 4) ) then 
            CALL show_chem 
!                                                                       
!     ----Show current configuration       'config'                     
!                                                                       
         ELSEIF (str_comp (cpara(1), 'config', 4, lpara(1), 6) )then
            CALL show_config 
!                                                                       
!     ----Show connectivity around an atom 'connect'                     
!                                                                       
         ELSEIF (str_comp (cpara(1), 'connect', 4, lpara(1), 7) ) then
            CALL del_params (1, ianz, cpara, lpara, maxw)
            IF(str_comp (cpara(ianz), 'long',3, lpara(ianz), 4)) THEN
               long = .true.
               ianz = ianz - 1
            ELSE
               long = .false.
            ENDIF
            iianz = 1
            CALL ber_params (iianz, cpara, lpara, werte, maxw) 
            iatom = NINT(werte(1))
            CALL del_params (1, ianz, cpara, lpara, maxw)
            CALL ber_params (ianz, cpara, lpara, werte, maxw) 
            IF(ier_num/=0) THEN
               c_name   = cpara(1)
               c_name_l = lpara(1)
               ino      = 0
               CALL no_error
            ELSE                                               ! Success set to value
               ino = nint (werte (1) ) 
               c_name   = ' '
               c_name_l = 1
            ENDIF
            CALL get_connectivity_identity( cr_iscat(iatom), ino, c_name, c_name_l)
            CALL do_show_connectivity ( iatom, ino, long)
!                                                                       
!     ----Show the dimensions              'cdim'                       
!                                                                       
         ELSEIF (str_comp (cpara (1) , 'cdim', 2, lpara (1) , 4) ) then 
            WRITE (output_io, 2000) ( (cr_dim (j, i), i = 1, 2),        &
            j = 1, 3)                                                   
!                                                                       
!     ----Show a domain                    'domain'                     
!                                                                       
         ELSEIF (str_comp (cpara (1) , 'domain', 2, lpara (1) , 6) )    &
         then                                                           
            CALL do_show_molecule (ianz, cpara, lpara, werte, maxw) 
!                                                                       
!     ----Show the atom environment        'envi'                       
!                                                                       
         ELSEIF (str_comp (cpara (1) , 'envi', 2, lpara (1) , 4) ) then 
            CALL do_show_env 
            cpara (2) = 'envi' 
            lpara (2) = 4 
            CALL do_show_atom (ianz, cpara, lpara, werte, maxw) 
!                                                                       
!     ----Show the molecular environment        'menvi'                 
!                                                                       
         ELSEIF (str_comp (cpara (1) , 'menvi', 3, lpara (1) , 5) )     &
         then                                                           
            CALL do_show_menv 
            IF (ianz.eq.2) then 
               IF (str_comp (cpara (2) , 'full', 1, lpara (2) , 4) )    &
               then                                                     
                  cpara (2) = 'envi' 
                  lpara (2) = 4 
                  CALL do_show_molecule (ianz, cpara, lpara, werte,     &
                  maxw)                                                 
               ENDIF 
            ENDIF 
!                                                                       
!     ----Show current unit cell metrics   'metric'                     
!                                                                       
         ELSEIF (str_comp (cpara (1) , 'metric', 3, lpara (1) , 6) )    &
         then                                                           
            CALL do_show_metric 
!                                                                       
!     ----Show a molecule                  'molecule'                   
!                                                                       
         ELSEIF (str_comp (cpara (1) , 'molecule', 2, lpara (1) , 8) )  &
         then                                                           
            CALL do_show_molecule (ianz, cpara, lpara, werte, maxw) 
!                                                                       
!     ----Show an object                   'object'                     
!                                                                       
         ELSEIF (str_comp (cpara (1) , 'object', 2, lpara (1) , 6) )    &
         then                                                           
            CALL do_show_molecule (ianz, cpara, lpara, werte, maxw) 
!                                                                       
!     ----Show scattering curve            'scat'                       
!                                                                       
         ELSEIF (str_comp (cpara (1) , 'scat', 1, lpara (1) , 4) ) then 
            CALL do_show_scat (ianz, cpara, lpara, werte, maxw) 
!                                                                       
!     ----Show Symmetry matrices           'symmetry'                   
!                                                                       
         ELSEIF (str_comp (cpara (1) , 'symmetry', 1, lpara (1) , 8) )  &
         then                                                           
            IF (ianz.eq.1.or.str_comp (cpara (2) , 'full', 2, lpara (2) &
            , 4) ) then                                                 
               mode = FULL 
            ELSEIF (str_comp (cpara (2) , 'symbol', 2, lpara (2) , 6) ) &
            then                                                        
               mode = SYMBOL 
            ELSEIF (str_comp (cpara (2) , 'xyz', 2, lpara (2) , 3) )    &
            then                                                        
               mode = XYZ 
            ELSEIF (str_comp (cpara (2) , 'matrix', 2, lpara (2) , 6) ) &
            then                                                        
               mode = MATRIX 
            ELSE 
               ier_num = - 6 
               ier_typ = ER_COMM 
               RETURN 
            ENDIF 
            CALL do_show_symmetry (mode) 
         ELSE 
            CALL do_show_generic (cpara, lpara, maxw) 
         ENDIF 
      ENDIF 
!                                                                       
 2000 FORMAT    (' Current crystal dimensions:'/                        &
     &                  '          minimum       maximum'/              &
     &                  '  X ',2(2x,f12.4)/                             &
     &                  '  Y ',2(2x,f12.4)/                             &
     &                  '  Z ',2(2x,f12.4))                             
!                                                                       
      END SUBROUTINE do_show                        
!*****7*****************************************************************
      SUBROUTINE show_asym 
!-                                                                      
!     This subroutine shows the content of the asymmetric unit          
!+                                                                      
      USE config_mod 
      USE crystal_mod 
!
      USE envir_mod 
      USE macro_mod 
      USE prompt_mod 
      IMPLICIT none 
!                                                                       
      CHARACTER(1) cdummy 
      INTEGER i, j, k, l 
!                                                                       
      DO j = 1, as_natoms / (LINES - 2) 
      WRITE (output_io, 2000) 
      DO i = 1, LINES - 2 
      k = (j - 1) * (LINES - 2) + i 
      WRITE (output_io, 2010) as_at_lis (k), as_iscat (k), (as_pos (l,  &
      k), l = 1, 3), as_dw (as_iscat (k) )                              
      ENDDO 
      IF (.not.lmakro) then 
         WRITE (output_io, 2020) 
         READ ( *, 1000) cdummy 
      ENDIF 
      ENDDO 
      WRITE (output_io, 2000) 
      j = as_natoms / (LINES - 2) + 1 
      DO i = 1, as_natoms - (as_natoms / (LINES - 2) ) * (lines - 2) 
      k = (j - 1) * (LINES - 2) + i 
      WRITE (output_io, 2010) as_at_lis (k), as_iscat (k), (as_pos (l,  &
      k), l = 1, 3), as_dw (as_iscat (k) )                              
      ENDDO 
!                                                                       
 1000 FORMAT    (a) 
 2000 FORMAT    (' Asymmetric unit'/ ' Name    Type       x',           &
     &           '             y             z               B')        
 2010 FORMAT    (1x,a4,2x,i4,3(2x,f12.5),4x,f12.5) 
 2020 FORMAT    (' Press RETURN for more ...') 
      END SUBROUTINE show_asym                      
!*****7*****************************************************************
      SUBROUTINE show_chem 
!-                                                                      
!     This subroutine shows the different atom types present            
!     in the crystal                                                    
!+                                                                      
      USE config_mod 
      USE crystal_mod 
      USE atom_name
!
      USE envir_mod 
      USE macro_mod 
      USE prompt_mod 
      IMPLICIT none 
       
!                                                                       
      CHARACTER(1) cdummy 
      INTEGER i, j, k 
      INTEGER kvoid 
!                                                                       
      CHARACTER(9) at_name_d 
!     CHARACTER(9) at_name 
!                                                                       
      kvoid = 0 
!                                                                       
      DO j = 1, (cr_nscat - kvoid) / (LINES - 2) 
      WRITE (output_io, 2000) 
      DO i = 1, LINES - 2 
      k = (j - 1) * (LINES - 2) + i + kvoid 
      at_name_d = at_name (k) 
      WRITE (output_io, 2010) at_name_d, k, cr_dw (k) 
      ENDDO 
      IF (.not.lmakro) then 
         WRITE (output_io, 2020) 
         READ ( *, 1000) cdummy 
      ENDIF 
      ENDDO 
      WRITE (output_io, 2000) 
      j = (cr_nscat - kvoid) / (LINES - 2) + 1 
      DO i = 1, cr_nscat - ( (cr_nscat - kvoid) / (LINES - 2) ) *       &
      (LINES - 2) - kvoid                                               
      k = (j - 1) * (LINES - 2) + i + kvoid 
      WRITE (output_io, 2010) cr_at_lis (k), k, cr_dw (k) 
      ENDDO 
 1000 FORMAT    (a)
 2000 FORMAT    (' Atoms present  '/ ' Name         Type       B')
 2010 FORMAT    (1x,a9,2x,i4,2x,f12.5)
 2020 FORMAT    (' Press RETURN for more ...')
      END SUBROUTINE show_chem
!*****7*****************************************************************
      SUBROUTINE do_show_env 
!-                                                                      
!     Shows the atom_environment array                                  
!+                                                                      
      USE config_mod 
      USE atom_env_mod 
!
      USE errlist_mod 
      USE prompt_mod 
      IMPLICIT none 
!                                                                       
!                                                                       
      INTEGER i, j, k, k1, k2 
!                                                                       
      IF (atom_env (0) .eq.0) then 
         WRITE (output_io, * ) 'Atom environment array is empty' 
      ELSE 
         WRITE (output_io, 2000) atom_env (0) 
         j = atom_env (0) / 5 
         IF (mod (atom_env (0), 5) .eq.0) then 
            j = j - 1 
         ENDIF 
         DO i = 0, j 
         k1 = 5 * i + 1 
         k2 = min (5 * i + 5, atom_env (0) ) 
         WRITE (output_io, 2001) (atom_env (k), k = k1, k2) 
         ENDDO 
      ENDIF 
!                                                                       
 2000 FORMAT    ( i8,' Atom numbers in the environment array') 
 2001 FORMAT    (5(2x,i13)) 
      END SUBROUTINE do_show_env                    
!*****7*****************************************************************
      SUBROUTINE do_show_menv 
!-                                                                      
!     Shows the molecular _environment array                            
!+                                                                      
      USE config_mod 
      USE mole_env_mod 
!                                                                       
      USE errlist_mod 
      USE prompt_mod 
      IMPLICIT none 
!                                                                       
      INTEGER i, j, k, k1, k2 
!                                                                       
      IF (mole_env (0) .eq.0) then 
         WRITE (output_io, * ) 'Molecular environment array is empty' 
      ELSE 
         WRITE (output_io, 2000) mole_env (0) 
         j = mole_env (0) / 5 
         IF (mod (mole_env (0), 5) .eq.0) then 
            j = j - 1 
         ENDIF 
         DO i = 0, j 
         k1 = 5 * i + 1 
         k2 = min (5 * i + 5, mole_env (0) ) 
         WRITE (output_io, 2001) (mole_env (k), k = k1, k2) 
         ENDDO 
      ENDIF 
!                                                                       
 2000 FORMAT    ( i8,' Molecule numbers in the environment array') 
 2001 FORMAT    (5(2x,i13)) 
      END SUBROUTINE do_show_menv                   
!*****7*****************************************************************
      SUBROUTINE do_show_atom (ianz, cpara, lpara, werte, maxw) 
!-                                                                      
!     Shows information about an atom.                                  
!+                                                                      
      USE config_mod 
      USE crystal_mod 
      USE atom_env_mod 
      USE atom_name
      USE modify_mod
!                                                                       
      USE errlist_mod 
      USE prompt_mod 
      IMPLICIT none 
       
!                                                                       
      INTEGER maxw 
!                                                                       
      CHARACTER ( * ) cpara (maxw) 
      INTEGER ianz 
      INTEGER lpara (maxw) 
      REAL werte (maxw) 
!                                                                       
      CHARACTER(32) c_property 
      INTEGER i, istart, iend, l ,k
      INTEGER, DIMENSION(1:3) :: ioffset
      INTEGER length 
!                                                                       
      CHARACTER(9) at_name_d 
!     CHARACTER(9) at_name 
!                                                                       
      LOGICAL str_comp 
!                                                                       
      cpara (1) = '0.0' 
      lpara (1) = 3 
!                                                                       
      IF (str_comp (cpara (2) , 'all', 1, lpara (2) , 3) ) then 
         werte (2) = 1 
         werte (3) = cr_natoms 
         ianz = 3 
      ELSEIF (str_comp (cpara (2) , 'last', 1, lpara (2) , 4) ) then 
         werte (2) = cr_natoms 
         werte (3) = cr_natoms 
         ianz = 3 
      ELSE 
         CALL ber_params (ianz, cpara, lpara, werte, maxw) 
      ENDIF 
      IF (ier_num.eq.0) then 
!                                                                       
!     --List atom environment                                           
!                                                                       
         IF (str_comp (cpara (2) , 'envi', 1, lpara (2) , 4) ) then 
            WRITE (output_io, 3000) 
            DO l = 1, atom_env (0) 
            i = atom_env (l) 
            at_name_d = at_name (cr_iscat (i) ) 
            CALL char_prop_1 (c_property, cr_prop (i), length) 
            WRITE (output_io, 3010) at_name_d, cr_pos (1, i), cr_pos (2,&
            i), cr_pos (3, i), cr_dw (cr_iscat (i) ), c_property (1:    &
            length)
            DO k=1,3
               ioffset(k) = NINT(atom_pos(k,l)-cr_pos(k,i))
            ENDDO
            WRITE (output_io, 3020) ioffset
            ENDDO 
!                                                                       
!     --List sequence of atoms                                          
!                                                                       
         ELSE 
            IF (ianz.eq.2.or.ianz.eq.3) then 
               istart = nint (werte (2) ) 
               IF (ianz.eq.2) then 
                  iend = istart 
               ELSE 
                  iend = nint (werte (3) ) 
               ENDIF 
               iend = min (iend, cr_natoms) 
               WRITE (output_io, 3000) 
               DO i = istart, iend 
               at_name_d = at_name (cr_iscat (i) ) 
               CALL char_prop_1 (c_property, cr_prop (i), length) 
               WRITE (output_io, 3010) at_name_d, cr_pos (1, i),        &
               cr_pos (2, i), cr_pos (3, i), cr_dw (cr_iscat (i) ),     &
               c_property (1:length)                                    
               ENDDO 
            ELSE 
               ier_num = - 6 
               ier_typ = ER_COMM 
            ENDIF 
         ENDIF 
      ENDIF 
!                                                                       
 3000 FORMAT    (' Name',11x,'x',13x,'y',13x,'z',13x,'B',9x,'Property') 
 3010 FORMAT    (1x,a9,3(2x,f12.6),4x,f10.6,2x,a) 
 3020 FORMAT    ( 3x  ,3(8x,i6   )              ) 
      END SUBROUTINE do_show_atom                   
!*****7*****************************************************************
      SUBROUTINE do_show_bval (ianz, cpara, maxw) 
!-                                                                      
!     Shows bond valence parameters for an atom.                        
!+                                                                      
      USE config_mod 
      USE crystal_mod 
      USE chem_mod 
      USE bv_data_mod
      USE element_data_mod
!                                                                       
      USE errlist_mod 
      USE prompt_mod 
      IMPLICIT none 
       
!                                                                       
      INTEGER maxw 
!                                                                       
      CHARACTER ( * ) cpara (maxw) 
      CHARACTER (LEN=4) el_name
      INTEGER ianz 
!                                                                       
      INTEGER ie_1, ie_2, i 
      INTEGER  :: ilook  ! lookup dummy for bv_index
!                                                                       
      IF (ianz.eq.3) then 
         CALL do_cap (cpara (2) ) 
         CALL do_cap (cpara (3) ) 
!         READ (cpara (2) (1:4), 1000) elem_1 
!         READ (cpara (3) (1:4), 1000) elem_2 
!         CALL symbf (elem_1, ie_1) 
!         CALL symbf (elem_2, ie_2) 
         el_name = cpara(2)
         CALL symbf ( el_name, ie_1 )
         el_name = cpara(3)
         CALL symbf ( el_name, ie_2 )
         CALL bv_lookup ( ie_1, ie_2, ilook )
         IF (ilook /=  0) then 
            i = bv_target(ilook)
            WRITE (output_io, 2000) cpara (2) (1:4), cpara (3) (1:4),   &
            bv_r0 (i), bv_b (i)                                         
         ELSE 
            ier_num = - 75 
            ier_typ = ER_APPL 
         ENDIF 
      ELSE 
         ier_num = - 6 
         ier_typ = ER_COMM 
      ENDIF 
!                                                                       
 2000 FORMAT    (' Bond valence parameters ',a4,' - ',a4,' : r0 = ',    &
     &                   f7.4,' b = ',f7.4)                             
      END SUBROUTINE do_show_bval                   
!*****7*****************************************************************
      SUBROUTINE do_show_molecule (ianz, cpara, lpara, werte, maxw) 
!-                                                                      
!     Shows information about a molecule.                               
!+                                                                      
      USE config_mod 
      USE crystal_mod 
      USE atom_name
      USE molecule_mod 
      USE mole_env_mod 
!                                                                       
      USE errlist_mod 
      USE prompt_mod 
      IMPLICIT none 
       
!                                                                       
      INTEGER maxw 
!                                                                       
      CHARACTER ( * ) cpara (maxw) 
      CHARACTER(15) C_MOLE ( - 4:4) 
      INTEGER ianz 
      INTEGER lpara (maxw) 
      REAL werte (maxw) 
!                                                                       
      INTEGER i, istart, iend, j, k 
!                                                                       
      CHARACTER(9) at_name_d 
!     CHARACTER(9) at_name 
!                                                                       
      INTEGER len_str 
      LOGICAL str_comp 
!                                                                       
      DATA C_MOLE / 'Domain    Fuzzy', 'Domain   Sphere', 'Domain Cylind&
     &er', 'Domain     Cube', '   Atoms       ', '    Cube       ', 'Cyl&
     &inder       ', '  Sphere       ', '    Edge       ' /             
!                                                                       
      cpara (1) = '0.0' 
      lpara (1) = 3 
!                                                                       
      IF (ianz.eq.2.or.ianz.eq.3) then 
         IF (str_comp (cpara (2) , 'all', 1, lpara (2) , 3) ) then 
            istart = 1 
            iend = mole_num_mole 
         ELSEIF (str_comp (cpara (2) , 'last', 1, lpara (2) , 5) ) then 
            istart = mole_num_mole 
            iend = mole_num_mole 
!                                                                       
!     --List molecular environment                                      
!                                                                       
         ELSEIF (str_comp (cpara (2) , 'envi', 1, lpara (2) , 4) ) then 
            WRITE (output_io, 3000) 
            DO i = 1, mole_env (0) 
            WRITE (output_io, 3000) mole_env (i), mole_type (mole_env ( &
            i) )                                                        
            DO j = 1, mole_len (i) 
            k = mole_cont (mole_off (i) + j) 
            at_name_d = at_name (cr_iscat (k) ) 
            WRITE (output_io, 3010) at_name_d, k, cr_pos (1, k),        &
            cr_pos (2, k), cr_pos (3, k), cr_dw (cr_iscat (k) )         
            ENDDO 
            ENDDO 
            RETURN 
         ELSE 
            CALL ber_params (ianz, cpara, lpara, werte, maxw) 
            IF (ier_num.ne.0) return 
            IF (ianz.eq.2) then 
               istart = nint (werte (2) ) 
               iend = istart 
            ELSE 
               istart = nint (werte (2) ) 
               iend = nint (werte (3) ) 
            ENDIF 
            istart = max (istart, 1) 
            iend = min (iend, mole_num_mole) 
         ENDIF 
      ELSE 
         ier_num = - 6 
         ier_typ = ER_COMM 
         RETURN 
      ENDIF 
!                                                                       
      DO i = istart, iend 
      IF (mole_char (i) .eq.MOLE_ATOM) then 
         WRITE (output_io, 3000) i, mole_type (i), C_MOLE (mole_char (i)&
         )                                                              
      ELSEIF (mole_char (i) .gt.MOLE_ATOM) then 
         WRITE (output_io, 4000) i, mole_type (i), C_MOLE (mole_char (i)&
         ), mole_dens (i)                                               
      ELSEIF (mole_char (i) .lt.MOLE_ATOM) then 
         j = max (1, len_str (mole_file (i) ) ) 
         WRITE (output_io, 5000) i, mole_type (i), C_MOLE (mole_char (i)&
         ), mole_file (i) (1:j), mole_fuzzy (i)                         
      ENDIF 
      DO j = 1, mole_len (i) 
      k = mole_cont (mole_off (i) + j) 
      at_name_d = at_name (cr_iscat (k) ) 
      WRITE (output_io, 3010) at_name_d, k, cr_pos (1, k), cr_pos (2, k)&
      , cr_pos (3, k), cr_dw (cr_iscat (k) )                            
      ENDDO 
      ENDDO 
!                                                                       
 3000 FORMAT(/' Molecule Number    : ',i11/                             &
     &       10x,'Type      : ',i11/                                    &
     &       10x,'Character : ',a8/                                     &
     &       ' Name',11x,'Number',6x,'x',13x,'y',13x,'z',15x,'B')       
 4000 FORMAT(/' Object   Number    : ',i11,/                            &
     &       10x,'Type      : ',i11/                                    &
     &       10x,'Character : ',a15/                                    &
     &       10x,'Density   : ',f13.4/                                  &
     &       ' Name',11x,'Number',6x,'x',13x,'y',13x,'z',15x,'B')       
 5000 FORMAT(/' Domain   Number    : ',i11/                             &
     &       10x,'Type      : ',i11/                                    &
     &       10x,'Character : ',a15/                                    &
     &       10x,'File      : ',a/                                      &
     &       10x,'Fuzzy     : ',f12.4/                                  &
     &       ' Name',11x,'Number',6x,'x',13x,'y',13x,'z',15x,'B')       
 3010 FORMAT(1x,a9,i11,1x,3(2x,f12.6),4x,f10.6) 
      END SUBROUTINE do_show_molecule               
!*****7*****************************************************************
      SUBROUTINE do_show_scat (ianz, cpara, lpara, werte, maxw) 
!-                                                                      
!     Shows the scattering curve of an atom type.                       
!+                                                                      
      USE config_mod 
      USE allocate_generic
      USE atom_name
      USE crystal_mod 
      USE diffuse_mod 
      USE fourier_sup
      USE modify_mod
      USE molecule_mod 
!                                                                       
      USE errlist_mod 
      USE prompt_mod 
      IMPLICIT none 
       
!                                                                       
      INTEGER maxw 
!                                                                       
      CHARACTER ( * ) cpara (maxw) 
      INTEGER ianz 
      INTEGER lpara (maxw) 
      REAL werte (maxw) 
!                                                                       
      LOGICAL lold 
      PARAMETER (lold = .false.) 
!                                                                       
      INTEGER i, j 
!                                                                       
      CHARACTER(9) at_name_d 
!     CHARACTER(9) at_name 
!                                                                       
      LOGICAL, DIMENSION(:),ALLOCATABLE :: latom ! (0:MAXSCAT) 
      INTEGER                           :: all_status
      INTEGER                           :: size_of
!
!-----allocate local logical array
!
      CALL alloc_arr ( latom  ,0,MAXSCAT,  all_status, .false.  , size_of )
!                                                                       
      DO i = 2, ianz 
      cpara (i - 1) = cpara (i) 
      lpara (i - 1) = lpara (i) 
      ENDDO 
      ianz = ianz - 1 
!                                                                       
      CALL dlink (ano, lambda, rlambda, diff_radiation, &
                  diff_power) 
!                                                                       
      DO i = 0, cr_nscat 
      latom (i) = .false. 
      ENDDO 
      CALL get_iscat (ianz, cpara, lpara, werte, maxw, lold) 
      IF (ier_num.eq.0) then 
         IF (werte (1) .eq. - 1) then 
            DO i = 0, cr_nscat 
            latom (i) = .true. 
            ENDDO 
         ELSE 
            DO i = 1, ianz 
            latom (nint (werte (i) ) ) = .true. 
            ENDDO 
         ENDIF 
      ENDIF 
!     IF (lxray) then 
!        WRITE (output_io, 3000) 'Xray', lambda, rlambda 
!     ELSE 
!        WRITE (output_io, 3000) 'Neutronen', lambda, rlambda 
!     ENDIF 
      SELECTCASE(diff_radiation)
         CASE(RAD_XRAY)
            WRITE (output_io, 3000) 'Xray', lambda, rlambda 
         CASE(RAD_NEUT)
            WRITE (output_io, 3000) 'Neutron', lambda, rlambda 
         CASE(RAD_ELEC)
            WRITE (output_io, 3000) 'Electron', lambda, rlambda 
      END SELECT
      IF (ano) then 
         WRITE (output_io, 3005) 'calculated' 
      ELSE 
         WRITE (output_io, 3005) 'ignored' 
      ENDIF 
      DO i = 0, cr_nscat 
      IF (latom (i) ) then 
         at_name_d = at_name ( (i) ) 
         SELECTCASE(diff_power)
           CASE(0)
             WRITE (output_io, 3030) at_name_d,  &
             cr_scat (1, i)
           CASE(4)
             WRITE (output_io, 3010) at_name_d,  &
             (cr_scat (2 * j    , i),j = 1, diff_power),      &
             (cr_scat (2 * j + 1, i),j = 1, diff_power),      &
             cr_scat (1, i), cr_delfr (i), cr_delfi (i)
           CASE(5)
             WRITE (output_io, 3020) at_name_d,  &
             (cr_scat (2 * j    , i),j = 1, diff_power),      &
             (cr_scat (2 * j + 1, i),j = 1, diff_power),      &
             cr_scat (1, i), cr_delfr (i), cr_delfi (i)
         END SELECT
      ENDIF 
      ENDDO 
!
      DEALLOCATE (latom)
!                                                                       
 3000 FORMAT    (' Scattering technique : ',a/                          &
     &                  ' Wave length symbol   : ',a/                   &
     &                  ' Wave length          : ',f7.5)                
 3005 FORMAT    (' Anomalous scat.      : ',a/) 
 3030 FORMAT    (1x,a9,    'bcoh        : ', f12.7)
 3010 FORMAT    (1x,a9,    'a(i)        : ',4(f12.7,2x)/                &
     &                  '          b(i)        : ',4(f12.7,2x)/         &
     &                  '          c           : ',  f12.7/             &
     &               '          f'',f''''      : ',2(f12.7,2x)/)        
 3020 FORMAT    (1x,a9,    'a(i)        : ',5(f12.7,2x)/                &
     &                  '          b(i)        : ',5(f12.7,2x)/         &
     &                  '          c           : ',  f12.7/             &
     &               '          f'',f''''      : ',2(f12.7,2x)/)        
!                                                                       
      END SUBROUTINE do_show_scat                   
!*****7*****************************************************************
      SUBROUTINE do_show_metric 
!-                                                                      
!     Shows the current unit cell dimensions                            
!+                                                                      
      USE config_mod 
      USE crystal_mod 
!                                                                       
      USE prompt_mod 
      IMPLICIT none 
!                                                                       
      INTEGER i, j 
!                                                                       
      WRITE (output_io, 2001) (cr_a0 (i), i = 1, 3), (cr_win (i),       &
      i = 1, 3), cr_v                                                   
      WRITE (output_io, 2002) ( (cr_gten (i, j), j = 1, 3), i = 1, 3) 
      WRITE (output_io, 2003) (cr_ar (i), i = 1, 3), (cr_wrez (i),      &
      i = 1, 3), cr_vr                                                  
      WRITE (output_io, 2004) ( (cr_rten (i, j), j = 1, 3), i = 1, 3) 
!                                                                       
 2001 FORMAT     ( ' Lattice constants :'                               &
     &           ,/,4x,'a',10x,'b',10x,'c', 9x,                         &
     &           'alpha',6x,'beta',7x,'gamma', 6x,'volume',             &
     &           /,6(2X,F9.5),2X,G12.6)                                 
 2002 FORMAT     (/' Metric Tensor     :'/(3(' ',3(2X,F11.5)/))) 
 2003 FORMAT     ( ' Reciprocal Lattice constants :'                    &
     &           ,/,4x,'a*', 9x,'b*', 9x,'c*', 8x,                      &
     &           'alpha*',5x,'beta*',6x,'gamma*', 5x,'volume',          &
     &           /,6(2X,F9.5),2X,G12.6)                                 
 2004 FORMAT     (/' Reciprocal metric tensor     : '/                  &
     &            (3(' ',3(2X,F11.5)/)))                                
!                                                                       
      END SUBROUTINE do_show_metric                 
!*****7*****************************************************************
      SUBROUTINE do_show_symmetry (mode) 
!-                                                                      
!     Shows all symmetry matrices for current space group               
!+                                                                      
      USE config_mod 
      USE crystal_mod 
      USE wyckoff_mod 
      USE unitcell_mod 
!                                                                       
      USE prompt_mod 
      IMPLICIT none 
!                                                                       
      INTEGER mode 
!                                                                       
      INTEGER FULL, SYMBOL, XYZ, MATRIX 
      PARAMETER (FULL = 0) 
      PARAMETER (SYMBOL = 1) 
      PARAMETER (XYZ = 2) 
      PARAMETER (MATRIX = 3) 
!                                                                       
      CHARACTER(26) cr_system (1:9) 
      INTEGER is, j 
      INTEGER igroup 
      INTEGER :: block = 1
      INTEGER n_center 
!                                                                       
      DATA cr_system / 'triclinic', 'monoclinic b-unique', 'monoclinic c&
     &-unique', 'orthorhombic', 'tetragonal', 'trigonal hexagonal axes',&
     & 'trigonal rhombohedral axes', 'hexagonal', 'cubic' /             
!                                                                       
      n_center = 1 
      IF (cr_spcgr (1:1) .eq.'P') then 
         n_center = 1 
      ELSEIF (cr_spcgr (1:1) .eq.'A') then 
         n_center = 2 
      ELSEIF (cr_spcgr (1:1) .eq.'B') then 
         n_center = 2 
      ELSEIF (cr_spcgr (1:1) .eq.'C') then 
         n_center = 2 
      ELSEIF (cr_spcgr (1:1) .eq.'I') then 
         n_center = 2 
      ELSEIF (cr_spcgr (1:1) .eq.'F') then 
         n_center = 4 
      ELSEIF (cr_spcgr (1:1) .eq.'R'.and.cr_syst.eq.6) then 
         n_center = 3 
      ENDIF 
      IF (gen_sta.eq.GEN_SYMM) then 
         block = spc_n / n_center 
      ELSEIF (gen_sta.eq.GEN_CENTER) then 
         block = n_center 
      ENDIF 
!                                                                       
      WRITE (output_io, 2100) cr_spcgr, cr_spcgrno, cr_system (cr_syst) 
!                                                                       
      DO is = 1, spc_n 
      IF (gen_sta.eq.GEN_SYMM) then 
         igroup = mod (is - 1, block) + 1 
      ELSEIF (gen_sta.eq.GEN_CENTER) then 
         igroup = (is - 1) / block + 1 
      ENDIF 
      IF (mode.eq.FULL) then 
         WRITE (output_io, 2200) is, igroup 
         WRITE (output_io, 2300) (spc_mat (1, j, is), j = 1, 4), spc_char (is),&
         (spc_mat (2, j, is), j = 1, 4), (spc_mat (3, j, is), j = 1, 4),&
         spc_xyz (is)                                                   
      ELSEIF (mode.eq.SYMBOL) then 
         WRITE (output_io, 3200) is, igroup, spc_char (is) 
      ELSEIF (mode.eq.XYZ) then 
         WRITE (output_io, 4200) is, igroup, spc_xyz (is) 
      ELSEIF (mode.eq.MATRIX) then 
         WRITE (output_io, 5200) is, igroup, (spc_mat (1, j, is),       &
         j = 1, 4), (spc_mat (2, j, is), j = 1, 4), (spc_mat (3, j, is),&
         j = 1, 4)                                                      
      ENDIF 
      IF (gen_sta.eq.GEN_SYMM.and.n_center.gt.1.and.mod (is - 1, spc_n /&
      n_center) + 1.eq.block) then                                      
         WRITE (output_io, * ) 
      ENDIF 
      ENDDO 
!                                                                       
 2100 FORMAT    (/,' Space group ',a16,' No.: ',i3,2x,a26) 
 2200 FORMAT    ('Symmetry No.      [',i3,']  (',i3,')') 
 2300 FORMAT    (  ' ( ',3(f4.1,', '),f8.5,' )','  ',a65,/,             &
     &                    ' ( ',3(f4.1,', '),f8.5,' )',/,               &
     &                    ' ( ',3(f4.1,', '),f8.5,' )','  ',a87,/)      
 3200 FORMAT    ('Symmetry No.      [',i3,']  (',i3,')  ',a65) 
 4200 FORMAT    ('Symmetry No.      [',i3,']  (',i3,')  ',a87) 
 5200 FORMAT    ('Symmetry No.      [',i3,']  (',i3,')  ',              &
     &                    ' ( ',3(f4.1,', '),f8.5,' )',/,               &
     &                32x,' ( ',3(f4.1,', '),f8.5,' )',/,               &
     &                32x,' ( ',3(f4.1,', '),f8.5,' )'   )              
!                                                                       
!                                                                       
      END SUBROUTINE do_show_symmetry               
END MODULE show_menu
