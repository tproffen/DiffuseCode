MODULE reserved_mod
!
IMPLICIT NONE
SAVE
!
INTEGER, PARAMETER :: lib_reserved_n    = 92
INTEGER, PARAMETER :: discus_reserved_n = 42 
INTEGER, PARAMETER :: diffev_reserved_n = 24
INTEGER, PARAMETER :: refine_reserved_n =  1
INTEGER, PARAMETER :: kuplot_reserved_n = 26
INTEGER, PARAMETER ::  suite_reserved_n =  1
                                                                        
CHARACTER(LEN=16), DIMENSION(1: suite_reserved_n) ::  suite_reserved 
CHARACTER(LEN=16), DIMENSION(1:kuplot_reserved_n) :: kuplot_reserved 
CHARACTER(LEN=16), DIMENSION(1:diffev_reserved_n) :: diffev_reserved 
CHARACTER(LEN=16), DIMENSION(1:refine_reserved_n) :: refine_reserved
CHARACTER(LEN=16), DIMENSION(1:discus_reserved_n) :: discus_reserved
CHARACTER(LEN=16), DIMENSION(1:lib_reserved_n)    :: lib_reserved
!                                                                       
DATA lib_reserved /                                                         &
'REF_GENERATION', 'REF_DIMENSION ', 'REF_CHILDREN  ', 'REF_MEMBER    ', &
'REF_COMPUTE   ', 'REF_NINDIV    ', 'REF_INDIV     ', 'variable   '   , &
'continue      ', 'REF_KID       ', 'fformat      ' , 'elseif   '     , &
'fclose'        , 'fexist'        , 'socket'        , 'system'        , &
'getcwd'        , 'getenv'        , 'endif'         , 'asind'         , &
'acosd'         , 'atand'         , 'enddo'         , 'while'         , &
'until'         , 'break'         , 'fopen'         , 'learn'         , &
'sleep'         , 'gskew'         , 'fdate'         , 'fmodt'         , &
'asin'          , 'acos'          , 'atan'          , 'sind'          , &
'cosd'          , 'tand'          , 'sinh'          , 'cosh'          , &
'tanh'          , 'sqrt'          , 'nint'          , 'frac'          , &
'gran'          , 'logn'          , 'else'          , 'then'          , &
'stop'          , 'exit'          , 'echo'          , 'eval'          , &
'fend'          , 'fput'          , 'fget'          , 'fsub'          , &
'help'          , 'lend'          , 'seed'          , 'show'          , &
'wait'          , 'pois'          , 'gbox'          , 'date'          , &
'exp'           , 'sin'           , 'cos'           , 'tan'           , &
'abs'           , 'mod'           , 'max'           , 'min'           , &
'int'           , 'and'           , 'xor'           , 'ran'           , &
'set'           , 'res'           , 'cd'            , 'PI'            , &
'ln'            , 'do'            , 'if'            , 'lt'            , &
'le'            , 'gt'            , 'ge'            , 'eq'            , &
'or'            , 'by'            , 'i'             , 'r'               &
/
!                                                                       
DATA discus_reserved /                                                  &
'mol_biso'      , 'mol_cont'      , 'mol_dens'      , 'mol_test'      , &
'mol_type'      , 'pdf_dens'      , 'pdf_scal'      , 'at_name'       , &
'at_type'       , 'mc_type'       , 'in_mole'       , 'mc_orig'       , &
'md_dist'       , 'md_next'       , 'md_test'       , 'mol_len'       , &
'scalpro'       , 'mc_rad'        , 'md_num'        , 'md_cre'        , &
'mc_num'        , 'md_rad'        , 'mr_run'        , 'dstar'         , &
'sym_n'         , 'bang'          , 'blen'          , 'cdim'          , &
'menv'          , 'rang'          , 'rlat'          , 'rvol'          , &
'occ'           , 'env'           , 'lat'           , 'vol'           , &
'b'             , 'm'             , 'n'             , 'x'             , &
'y'             , 'z'                                                   &
/         
!
DATA diffev_reserved /                                                  &
'child_val     ', 'pop_xmin      ', 'pop_xmax      ', 'pop_smin      ', &
'pop_smax '     , 'pop_lsig '     , 'pop_dimx '     , 'diff_sel '     , &
'pop_sig  '     , 'pop_gen  '     , 'diff_lo  '     , 'diff_cr  '     , &
'worstr   '     , 'worstm   '     , 'rvalue   '     , 'diff_k   '     , &
'diff_f   '     , 'pop_v    '     , 'pop_t    '     , 'pop_n    '     , &
'pop_c    '     , 'bestr    '     , 'bestm    '     , 'p        '       &
/
!
DATA refine_reserved /                                                  &
'ref_nonsense  '                                                        &
/
!
DATA kuplot_reserved /                                                  &
'F_DERIV       ', 'F_VALUE       ', 'F_PARA        ',                                     &
'zmax          ', 'zmin          ', 'ymax          ', 'ymin          ', 'xmax          ', &
'xmin          ', 'size          ', 'pwin          ', 'cmax          ', 'cmap          ', &
'axis          ', 'ny            ', 'nx            ', 'np            ', 'ni            ', &
'dy            ', 'dx            ', 'z             ', 'y             ', 'x             ', &
's             ', 'p             ', 'n             '                          &
/
DATA  suite_reserved /                                                  &
'              '                                                        &
/
!
END MODULE reserved_mod
