MODULE suite
!
!  Module to interface python with the discus_suite
!
!  reinhard.neder@fau.de
!
!  Version initial test not ready for release
!
USE suite_python_support
!
PRIVATE
PUBLIC initialize_suite    ! Initialize the discus_suite as if started directly
PUBLIC interactive         ! start an interactive suite session
PUBLIC execute_macro       ! Execute a macro at suite, discus, diffev, kuplot
PUBLIC execute_help        ! Execute the help
PUBLIC execute_command     ! Execute a single command
PUBLIC test_macro_param    ! Test a macro for the number of parameters required by the macro
PUBLIC gui_do_init         ! Initialize a Do/IF block
PUBLIC gui_do_insert       ! Read Commands into a Do/If block and execute once finished
PUBLIC send_i              ! Send an integer array to the suite
PUBLIC send_r              ! Send a real valued array to the suite
PUBLIC get_i               ! Get an integer valued array from the suite
PUBLIC get_r               ! Get a real valued array from the suite
PUBLIC suite_learn         ! Transfer a line for learning
PUBLIC discus_get_wave_number  ! Interface to get number of wave length symbols
PUBLIC discus_get_wave_symbol  ! Interface to get wave length symbols and values
PUBLIC discus_get_spcgr_number ! Interface to get number of wave length symbols
PUBLIC discus_get_spcgr_symbol ! Interface to get wave length symbols and values
PUBLIC discus_get_nscat_number  ! Interface to get number of wave length symbols
PUBLIC discus_get_scat_symbol  ! Interface to get wave length symbols and values
PUBLIC discus_get_natoms       ! Interface to get number of atoms in the structure
PUBLIC discus_read_structure   ! Use discus/read to read a structure or unit cell
PUBLIC discus_calc_fourier     ! Use discus/fourier to calculate a Fourier
PUBLIC discus_get_fourier      ! Interface to get Fourier menu items from DISCUS
PUBLIC discus_get_four_last    ! Interface to get last Fourier type
PUBLIC discus_set_fourier      ! Interface to set Fourier menu items to   DISCUS
PUBLIC discus_get_powder       ! Interface to get Powder  menu items from DISCUS
PUBLIC discus_calc_powder      ! Use discus/do_powder to calculate a Powder
PUBLIC discus_get_pdf          ! Interface to get PDF     menu items from DISCUS
PUBLIC discus_calc_pdf         ! Use discus/pdf to calculate a PDF
PUBLIC discus_get_save         ! Use discus/pdf to save a crystal structure
PUBLIC discus_run_save         ! Use discus/pdf to save a crystal structure
PUBLIC discus_output           ! Interface to run OUTPUT
PUBLIC discus_run_sro          ! Interface to run SRO
PUBLIC kuplot_load             ! Use kuplot/load to load a data set
!
CONTAINS
!
SUBROUTINE initialize_suite()
!
!   Initialization of the discus_suite, to be run only once
!
USE suite_setup_mod
USE suite_loop_mod
USE suite_init_mod
USE discus_setup_mod
USE kuplot_setup_mod
USE diffev_setup_mod
USE diffev_mpi_mod
USE run_mpi_mod
!
USE charact_mod
USE prompt_mod
USE envir_mod
USE terminal_mod
!
IMPLICIT NONE
!
INTEGER, PARAMETER :: master = 0 ! Master ID for MPI
EXTERNAL :: suite_sigint
!
run_mpi_myid      = 0
lstandalone       = .false.      ! No standalone for DIFFEV, DISCUS, KUPLOT
!
IF( .NOT. lsetup_done ) THEN    ! ! If necessary do initial setup
   CALL run_mpi_init    ! Do the initial MPI configuration
   CALL setup_suite     ! Define initial parameter, array values
   CALL suite_set_sub
   CALL discus_setup(lstandalone)
   CALL kuplot_setup(lstandalone)
   CALL diffev_setup(lstandalone)
   suite_discus_init = .TRUE.
   suite_kuplot_init = .TRUE.
   suite_diffev_init = .TRUE.
   pname     = 'suite'
   pname_cap = 'SUITE'
   prompt    = pname
   hlpfile   = hlpdir(1:hlp_dir_l)//pname(1:LEN(TRIM(pname)))//'.hlp'
   hlpfile_l = LEN(TRIM(hlpfile))
   IF(.NOT.run_mpi_active) THEN
      CALL suite_set_sub_cost ()
   ENDIF
   lsetup_done = .TRUE.
ELSE
   CALL suite_set_sub
ENDIF
lstandalone = .false.
WRITE(output_io,'(a5,a,a5)') COLOR_HIGH,'Control turned to GUI ...',COLOR_FG_DEFAULT
!
END SUBROUTINE initialize_suite
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
SUBROUTINE interactive(prog)
!
!  Generic interface routine to start an interactive discus_suite session
!  from the python host
!
USE suite_loop_mod
USE discus_loop_mod
USE diffev_loop_mod
USE kuplot_loop_mod
USE prompt_mod
!
IMPLICIT NONE
!
CHARACTER(LEN=*), INTENT(IN)    :: prog
!
IF( .NOT. lsetup_done) CALL initialize_suite
linteractive = .TRUE.
section: SELECT CASE (prog)
   CASE ('suite')
      CALL suite_prae
      CALL suite_loop      ! Perform the normal main loop
   CASE ('discus')
      CALL discus_prae
      CALL discus_loop     ! Perform the normal discus loop
   CASE ('diffev')
      CALL diffev_prae
      CALL diffev_loop     ! Perform the normal discus loop
   CASE ('kuplot')
      CALL kuplot_prae
      CALL kuplot_loop     ! Perform the normal discus loop
END SELECT section
lsetup_done = .TRUE.
WRITE(output_io,'(a)') 'Control returned to GUI ...'
CALL back_to_suite      ! Go back to the suite
!
END SUBROUTINE interactive
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
SUBROUTINE execute_macro(prog, line)
!
!  Execute the macro given on line for the section defined by prog
!  The macro name and all parameters must be specified, parameters
!  must be separated from each other by a comma.
!
USE suite_loop_mod
USE discus_loop_mod
USE diffev_loop_mod
USE kuplot_loop_mod
!
USE charact_mod
USE prompt_mod
USE terminal_mod
!
IMPLICIT NONE
!
CHARACTER(LEN=*), INTENT(IN)    :: prog
CHARACTER(LEN=*), INTENT(INOUT) :: line
!
INTEGER :: length
!
linteractive = .FALSE.
section: SELECT CASE (prog)
   CASE ('suite')
      CALL suite_prae
   CASE ('discus')
      CALL discus_prae
   CASE ('diffev')
      CALL diffev_prae
   CASE ('kuplot')
      CALL kuplot_prae
END SELECT section
!
length = LEN_TRIM(line)
IF(line(1:1) == '@' ) THEN
   line = line(2:length)
   length = length - 1
ENDIF
!
WRITE(output_io,'(a5,''@''a,a5)') COLOR_INFO,line(1:length),COLOR_FG_DEFAULT
CALL file_kdo(line,length)
exec: SELECT CASE (prog)
   CASE ('suite')
      CALL suite_loop
   CASE ('discus')
      CALL discus_loop
   CASE ('diffev')
      CALL diffev_loop
   CASE ('kuplot')
      CALL kuplot_loop
END SELECT exec
!
linteractive = .TRUE.
WRITE(output_io,'(a)') 'Control returned to GUI ...'
CALL back_to_suite      ! Go back to the suite
!
END SUBROUTINE execute_macro
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
SUBROUTINE execute_help(prog)
!
!  Execute the help for the section defined by prog
!
!
USE charact_mod
USE prompt_mod
USE terminal_mod
!
IMPLICIT NONE
!
CHARACTER(LEN=*), INTENT(IN)    :: prog
!
INTEGER :: length
!
linteractive = .TRUE.
section: SELECT CASE (prog)
   CASE ('suite')
      CALL suite_prae
   CASE ('discus')
      CALL discus_prae
   CASE ('diffev')
      CALL diffev_prae
   CASE ('kuplot')
      CALL kuplot_prae
END SELECT section
!
length = LEN_TRIM(prog)
CALL do_hel(prog, length)
!
linteractive = .TRUE.
WRITE(output_io,'(a5,a,a5)') COLOR_HIGH,'Control returned to GUI ...',COLOR_FG_DEFAULT
CALL back_to_suite      ! Go back to the suite
!
END SUBROUTINE execute_help
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
SUBROUTINE execute_command(prog,line)
!
!   Executes a single command as given on line
!   Returns to the GUI immediately
!
USE discus_loop_mod
USE diffev_loop_mod
USE kuplot_loop_mod
USE suite_loop_mod
USE charact_mod
USE prompt_mod
USE terminal_mod
IMPLICIT NONE
!
CHARACTER(LEN=*), INTENT(IN) :: prog
CHARACTER(LEN=*), INTENT(IN) :: line
INTEGER                      :: length
LOGICAL                      :: lend
!
IF( .NOT. lsetup_done) CALL initialize_suite   ! Do we need to initialize?
!
length = LEN_TRIM(line)
lend   = .FALSE.
linteractive = .FALSE.
input_gui = line
WRITE(output_io,'(a5,a,a5)') COLOR_INFO,line(1:length),COLOR_FG_DEFAULT
section: SELECT CASE (prog)
   CASE ('suite')
      CALL suite_prae
      CALL suite_loop
   CASE ('discus')
      CALL discus_prae
      CALL discus_loop
   CASE ('diffev')
      CALL diffev_prae
      CALL diffev_loop
   CASE ('kuplot')
      CALL kuplot_prae
      CALL kuplot_loop
END SELECT section
CALL back_to_suite      ! Go back to the suite
linteractive=.TRUE.     ! Tell get_cmd to read input from standard I/O
CALL back_to_suite      ! Go back to the suite
!
END SUBROUTINE execute_command
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
INTEGER FUNCTION test_macro_param(line)
!
CHARACTER(LEN=*), INTENT(IN) :: line
!
INTEGER :: length
INTEGER :: numpar
!
length = LEN_TRIM(line)
CALL test_macro(line,length, numpar)
test_macro_param = numpar
!
END FUNCTION test_macro_param
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
SUBROUTINE gui_do_init(prog,line)
!
USE charact_mod
USE doact_mod
USE prompt_mod
USE terminal_mod
!
IMPLICIT NONE
!
CHARACTER(LEN=*), INTENT(IN) :: prog
CHARACTER(LEN=*), INTENT(IN) :: line
INTEGER                      :: length
LOGICAL                      :: lend
!
linteractive = .FALSE.
section: SELECT CASE (prog)
   CASE ('suite')
      CALL suite_prae
   CASE ('discus')
      CALL discus_prae
   CASE ('diffev')
      CALL diffev_prae
   CASE ('kuplot')
      CALL kuplot_prae
END SELECT section
!
length = LEN_TRIM(line)
WRITE(output_io,'(a5,a,a5)') COLOR_INFO,line(1:length), COLOR_FG_DEFAULT
IF(.NOT.lblock_read) THEN   ! This is the first DO/IF statement
   CALL do_do_init (line, lend, length)
ELSE
   input_gui = line
   CALL do_insert_line
ENDIF
CALL back_to_suite      ! Go back to the suite
!
length = LEN_TRIM(prog)
END SUBROUTINE gui_do_init
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
SUBROUTINE gui_do_insert(prog,line)
!
USE charact_mod
USE doact_mod
USE doexec_mod
USE errlist_mod
USE class_macro_internal
USE prompt_mod
USE terminal_mod
!
IMPLICIT NONE
!
CHARACTER(LEN=*), INTENT(IN) :: prog
CHARACTER(LEN=*), INTENT(IN) :: line
LOGICAL                      :: lend
!
IF(lblock_read) THEN    ! Only if we are reading into a Do/If block
   linteractive = .FALSE.
   section: SELECT CASE (prog)
   CASE ('suite')
         CALL suite_prae
      CASE ('discus')
         CALL discus_prae
      CASE ('diffev')
         CALL diffev_prae
      CASE ('kuplot')
         CALL kuplot_prae
   END SELECT section
!
   insert: DO
      IF(level > -1) THEN
         WRITE(output_io,'(a5,a,a5)') COLOR_INFO,line(1:LEN_TRIM(line)), COLOR_FG_DEFAULT
         input_gui = line
         CALL do_insert_line
      ELSE
         EXIT insert
      ENDIF
      IF(.NOT.lmakro) EXIT insert
   ENDDO insert
   IF(level < 0) THEN   ! Reached last enddo/endif, execute block
      lblock_read = .FALSE.
      CALL do_execute_block(lend)
      lblock = .false.
   ENDIF
!
   IF(ier_num /= 0) THEN
      CALL errlist
   ENDIF
!
   CALL back_to_suite      ! Go back to the suite
!
ENDIF
!
!
END SUBROUTINE gui_do_insert
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!
! These are generic send and get routines that are used by 
! the programs in the interface to python. They allow
! python to send / get parts of the arrays i[] and r[]
! to / from (discus, diffev, kuplot)
!
! This section used to be an independent file in lib_f90.
! As the pythion interface has moved to an explcit 
! directorty it is no longer needed in lib_f90
!
SUBROUTINE send_i (iin, lower, upper )
!
! The outer routine sends integer valued numbers for i[:]
!
USE param_mod
USE prompt_mod
IMPLICIT NONE
!
INTEGER,                         INTENT(IN) :: lower
INTEGER,                         INTENT(IN) :: upper
INTEGER, DIMENSION(lower:upper), INTENT(IN) :: iin
!
IF( .not. lsetup_done ) THEN    ! If necessary do initial setup
   CALL initialize_suite
ENDIF
!
IF(lower>0 .and. upper<UBOUND(inpara,1) .and. lower <= upper) THEN
   inpara(lower:upper) = iin(lower:upper)
ENDIF
!
END SUBROUTINE send_i
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
SUBROUTINE send_r (rin, lower, upper )
!
! The outer routine sends real valued numbers for r[:]
!
USE param_mod
USE prompt_mod
IMPLICIT NONE
!
INTEGER,                      INTENT(IN) :: lower
INTEGER,                      INTENT(IN) :: upper
REAL, DIMENSION(lower:upper), INTENT(IN) :: rin
!
IF( .not. lsetup_done ) THEN    ! If necessary do initial setup
   CALL initialize_suite
ENDIF
!
IF(lower>0 .and. upper<UBOUND(rpara,1) .and. lower <= upper) THEN
   rpara(lower:upper) = rin(lower:upper)
ENDIF
!
END SUBROUTINE send_r
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
SUBROUTINE get_i (iout, lower, upper )
!
! The outer routine gets integer valued numbers from i[:]
!
USE param_mod
USE prompt_mod
IMPLICIT NONE
!
INTEGER,                         INTENT(IN ) :: lower
INTEGER,                         INTENT(IN ) :: upper
INTEGER, DIMENSION(lower:upper), INTENT(OUT) :: iout
!
IF( .not. lsetup_done ) THEN    ! If necessary do initial setup
   CALL initialize_suite
ENDIF
!
IF(lower>0 .and. upper<UBOUND(inpara,1) .and. lower <= upper) THEN
   iout(lower:upper) = inpara(lower:upper) 
ENDIF
!
END SUBROUTINE get_i
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
SUBROUTINE get_r (rout, lower, upper )
!
! The outer routine gets real valued numbers from r[:]
!
USE param_mod
USE prompt_mod
IMPLICIT NONE
!
INTEGER,                      INTENT(IN ) :: lower
INTEGER,                      INTENT(IN ) :: upper
REAL, DIMENSION(lower:upper), INTENT(OUT) :: rout
!
IF( .not. lsetup_done ) THEN    ! If necessary do initial setup
   CALL initialize_suite
ENDIF
!
IF(lower>0 .and. upper<UBOUND(rpara,1) .and. lower <= upper) THEN
   rout(lower:upper) = rpara(lower:upper) 
ENDIF
!
END SUBROUTINE get_r
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
SUBROUTINE suite_learn(line)
!
!  This subroutine is called every time the GUI starts a menu in order 
!  to save the relevant menu changing command to a learn sequence. 
!
USE learn_mod
!
CHARACTER(LEN=*), INTENT(IN) :: line
!
IF(llearn) THEN
   IF(LEN_TRIM(line)>0) THEN
      WRITE(33, '(a)') line(1:LEN_TRIM(LINE))
   ENDIF
ENDIF
END SUBROUTINE suite_learn
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
SUBROUTINE discus_get_wave_number(nwave)
!
USE element_data_mod
IMPLICIT NONE
!
INTEGER, INTENT(OUT) :: nwave
!
nwave = PER_MAX_WAVE
!
END SUBROUTINE discus_get_wave_number
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
SUBROUTINE discus_get_wave_symbol(i,symbols, wavelengths)
!
USE element_data_mod
IMPLICIT NONE
!
INTEGER, PARAMETER :: MAX_WAVE = 40
INTEGER           , INTENT(IN ) :: i
CHARACTER  (LEN=4), INTENT(OUT) :: symbols
REAL              , INTENT(OUT) :: wavelengths
!
CALL get_sym_length(i,symbols, wavelengths)  ! Get names and values from element module
!
END SUBROUTINE discus_get_wave_symbol
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
SUBROUTINE discus_get_spcgr_number(nspcgr)
!
USE spcgr_mod
IMPLICIT NONE
!
INTEGER, INTENT(OUT) :: nspcgr
!
nspcgr = SPCGR_MAX
!
END SUBROUTINE discus_get_spcgr_number
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
SUBROUTINE discus_get_spcgr_symbol(i,get_spcgr_name, get_spcgr_syst)
!
USE spcgr_mod
IMPLICIT NONE
INTEGER          , INTENT(IN)  :: i
CHARACTER(LEN=16), INTENT(OUT) :: get_spcgr_name
INTEGER          , INTENT(OUT) :: get_spcgr_syst
!
get_spcgr_name = spcgr_name(i)
get_spcgr_syst = spcgr_syst(i)
!
END SUBROUTINE discus_get_spcgr_symbol
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
SUBROUTINE discus_get_nscat_number(nscat)
!
USE crystal_mod
IMPLICIT NONE
!
INTEGER, INTENT(OUT) :: nscat
!
nscat = cr_nscat
!
END SUBROUTINE discus_get_nscat_number
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
SUBROUTINE discus_get_scat_symbol(i,at_list)
!
USE crystal_mod
IMPLICIT NONE
INTEGER          , INTENT(IN)  :: i
CHARACTER(LEN= 4), INTENT(OUT) :: at_list
!
at_list = cr_at_lis(i)
!
END SUBROUTINE discus_get_scat_symbol
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
SUBROUTINE discus_get_natoms(natoms)
!
USE crystal_mod
IMPLICIT NONE
!
INTEGER, INTENT(OUT) :: natoms
!
natoms = cr_natoms
!
END SUBROUTINE discus_get_natoms
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
SUBROUTINE discus_read_structure(line)
!
!  A first interface that allows to read a structre from python via
!  suite.read_cell( python_string )
!  where python string is any of: 
!       cell      crystal_structure.cell, nx, ny, nz
!       lcell     crystal_structure.cell, nx, ny, nz
!       structure crystal_structure.cell
!       free      [optional parameters]
!
USE structur
USE prompt_mod
IMPLICIT NONE
!
CHARACTER(LEN=*), INTENT(IN) :: line
!
IF( .NOT. lsetup_done) CALL initialize_suite   ! Do we need to initialize?
!
linteractive=.FALSE.    ! Tell get_cmd to get input from input_gui
CALL discus_prae        ! Switch to discus section
input_gui = line        ! copy the input line to the automatic command line
CALL read_struc         ! Call the actual task at hand
CALL back_to_suite      ! Go back to the suite
linteractive=.TRUE.     ! Tell get_cmd to read input from standard I/O
!
END SUBROUTINE discus_read_structure
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
SUBROUTINE discus_calc_fourier(line)
!
!  A first interface that allows to calculate a Fourier from python via
!  suite.fourier( python_string )
!  where python string is any Fourier command.
!
USE fourier_menu
USE prompt_mod
IMPLICIT NONE
!
CHARACTER(LEN=*), INTENT(IN) :: line
!
IF( .NOT. lsetup_done) CALL initialize_suite   ! Do we need to initialize?
!
linteractive=.FALSE.    ! Tell get_cmd to get input from input_gui
CALL discus_prae        ! Switch to discus section
input_gui = line        ! copy the input line to the automatic command line
CALL fourier            ! Call the actual task at hand
CALL back_to_suite      ! Go back to the suite
linteractive=.TRUE.     ! Tell get_cmd to read input from standard I/O
!
END SUBROUTINE discus_calc_fourier
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
SUBROUTINE discus_get_fourier(corners, increment, radiation, element, wavelength, &
           energy, use_adp, use_ano, percent, lot_type, lot_numb, lot_dim, lot_per)
!
!  Interface to get all Fourier Menu items from DISCUS to Python
!
USE ISO_C_BINDING
USE element_data_mod
USE diffuse_mod
USE fourier_sup
IMPLICIT NONE
!
REAL (8), DIMENSION(1:4, 1:3), INTENT(OUT) :: corners
INTEGER,  DIMENSION(1:3)     , INTENT(OUT) :: increment
INTEGER                      , INTENT(OUT) :: radiation
CHARACTER (LEN=4)            , INTENT(OUT) :: element
REAL (8)                     , INTENT(OUT) :: wavelength
REAL (8)                     , INTENT(OUT) :: energy
INTEGER                      , INTENT(OUT) :: use_adp
INTEGER                      , INTENT(OUT) :: use_ano
REAL (8)                     , INTENT(OUT) :: percent
INTEGER                      , INTENT(OUT) :: lot_type
INTEGER                      , INTENT(OUT) :: lot_numb
INTEGER, DIMENSION(1:3)      , INTENT(OUT) :: lot_dim
INTEGER                      , INTENT(OUT) :: lot_per
!
CALL dlink (ano, lambda, rlambda, renergy, l_energy, &
            diff_radiation, diff_power) 
corners    = TRANSPOSE(eck)
increment  = inc
radiation  = diff_radiation - 1   ! Python counts 0 to 2
element    = lambda
wavelength = rlambda
energy     = renergy
IF(ldbw) THEN
   use_adp = 1
ELSE
   use_adp = 0
ENDIF
IF(ano) THEN
   use_ano = 1
ELSE
   use_ano = 0
ENDIF
percent  = fave
lot_type = ilots - 1
lot_numb = nlots
lot_dim(:)  = ls_xyz(:)
IF(lperiod) THEN
   lot_per = 1
ELSE
   lot_per = 0
ENDIF

!
END SUBROUTINE discus_get_fourier
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
SUBROUTINE discus_get_four_last(last_type)
!
!  Interface to get the last Fourier type calculated
!
USE ISO_C_BINDING
USE diffuse_mod
IMPLICIT NONE
!
INTEGER, INTENT(OUT) :: last_type
!
last_type = four_last
!
END SUBROUTINE discus_get_four_last
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
SUBROUTINE discus_set_fourier(corners, increment, radiation, wavelength, &
           energy, use_adp, use_ano, lot_type, lot_numb, lot_dim, lot_per)
!
!  Interface to set all Fourier Menu items in DISCUS from PYTHON
!
USE ISO_C_BINDING
USE diffuse_mod
USE fourier_sup
IMPLICIT NONE
!
REAL (8), DIMENSION(1:4, 1:3), INTENT(IN) :: corners
INTEGER,  DIMENSION(1:3)     , INTENT(IN) :: increment
INTEGER                      , INTENT(IN) :: radiation
REAL (8)                     , INTENT(IN) :: wavelength
REAL (8)                     , INTENT(IN) :: energy
INTEGER                      , INTENT(IN) :: use_adp
INTEGER                      , INTENT(IN) :: use_ano
INTEGER                      , INTENT(IN) :: lot_type
INTEGER                      , INTENT(IN) :: lot_numb
INTEGER, DIMENSION(1:3)      , INTENT(IN) :: lot_dim
INTEGER                      , INTENT(IN) :: lot_per
!
!CALL dlink (ano, lambda, rlambda, renergy, l_energy, &
!            diff_radiation, diff_power) 
eck        = REAL(TRANSPOSE(corners))
inc        = increment(:)
diff_radiation  = radiation + 1   ! Python counts 0 to 2
rlambda    = REAL(wavelength)
renergy    = REAL(energy)
ldbw       = use_adp==1
ano        = use_ano==1
ilots      = lot_type +1          ! Python counts 0 to 2
nlots      = lot_numb
ls_xyz     = lot_dim(:)
lperiod    = lot_per==1

!
END SUBROUTINE discus_set_fourier
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
SUBROUTINE discus_get_powder(calc_mode, axis, radiation, element, wavelength, &
           energy, use_adp, use_ano, theta, qvalues, dhkl, shkl,  &
           profile, profile_eta, profile_uvw, profile_asy, profile_width, &
           profile_delta, preferred, pref_dp, pref_pt, pref_hkl, &
           lp, lp_ang, lp_fac )
!
!  Interface to get all Powder Menu items from DISCUS to Python
!
USE ISO_C_BINDING
USE element_data_mod
USE diffuse_mod
USE powder_mod
USE fourier_sup
IMPLICIT NONE
!
INTEGER                      , INTENT(OUT) :: calc_mode
INTEGER                      , INTENT(OUT) :: axis
INTEGER                      , INTENT(OUT) :: radiation
CHARACTER (LEN=4)            , INTENT(OUT) :: element
REAL (8)                     , INTENT(OUT) :: wavelength
REAL (8)                     , INTENT(OUT) :: energy
INTEGER                      , INTENT(OUT) :: use_adp
INTEGER                      , INTENT(OUT) :: use_ano
REAL (8), DIMENSION(1:3)     , INTENT(OUT) :: theta
REAL (8), DIMENSION(1:3)     , INTENT(OUT) :: qvalues
REAL (8), DIMENSION(1:3)     , INTENT(OUT) :: dhkl
REAL (8), DIMENSION(1:3)     , INTENT(OUT) :: shkl
INTEGER                      , INTENT(OUT) :: profile
REAL (8)                     , INTENT(OUT) :: profile_eta
REAL (8), DIMENSION(1:3)     , INTENT(OUT) :: profile_uvw
REAL (8), DIMENSION(1:4)     , INTENT(OUT) :: profile_asy
REAL (8)                     , INTENT(OUT) :: profile_width
REAL (8)                     , INTENT(OUT) :: profile_delta
INTEGER                      , INTENT(OUT) :: preferred
REAL (8)                     , INTENT(OUT) :: pref_dp
REAL (8)                     , INTENT(OUT) :: pref_pt
REAL (8), DIMENSION(1:3)     , INTENT(OUT) :: pref_hkl
INTEGER                      , INTENT(OUT) :: lp
REAL (8)                     , INTENT(OUT) :: lp_fac
REAL (8)                     , INTENT(OUT) :: lp_ang
!
IF(pow_four_type==POW_HIST) THEN
   calc_mode  = 0
ELSE
   calc_mode = 1
ENDIF
IF(pow_axis==POW_AXIS_Q) THEN
   axis       = 0
ELSEIF(pow_axis==POW_AXIS_TTH) THEN
   axis       = 1
ENDIF
CALL dlink (ano, lambda, rlambda, renergy, l_energy, &
            diff_radiation, diff_power) 
radiation  = diff_radiation - 1   ! Python counts 0 to 2
element    = lambda
wavelength = rlambda
energy     = renergy
IF(ldbw) THEN
   use_adp = 1
ELSE
   use_adp = 0
ENDIF
IF(ano) THEN
   use_ano = 1
ELSE
   use_ano = 0
ENDIF
theta(1)   = pow_tthmin
theta(2)   = pow_tthmax
theta(3)   = pow_deltatth
qvalues(1) = pow_qmin
qvalues(2) = pow_qmax
qvalues(3) = pow_deltaq
!
dhkl(:)        = pow_hkl_del(:)
shkl(:)        = pow_hkl_shift(:)
!
profile        = pow_profile
profile_eta    = pow_eta
profile_uvw(1) = pow_u
profile_uvw(2) = pow_v
profile_uvw(3) = pow_w
profile_asy(1) = pow_p1
profile_asy(2) = pow_p2
profile_asy(3) = pow_p3
profile_asy(4) = pow_p4
profile_width  = pow_width
profile_delta  = pow_delta
IF(pow_pref) THEN
   preferred   = pow_pref_type
ELSE
   preferred   = 0
ENDIF
pref_dp        = pow_pref_g1
pref_pt        = pow_pref_g2
pref_hkl(:)    = pow_pref_hkl(:)
!
lp             = pow_lp
lp_ang         = pow_lp_ang
lp_fac         = pow_lp_fac
!
END SUBROUTINE discus_get_powder
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
SUBROUTINE discus_get_pdf(radiation, use_adp, r_max, r_step, lrho0, rho0, &
           corr_lin, corr_quad, period, exact, weight, finite, diameter,  &
           qmax, qbroad, &
           qdamp )
!
USE ISO_C_BINDING
USE diffuse_mod
USE pdf_mod
USE fourier_sup
USE chem_mod
!
INTEGER                      , INTENT(OUT) :: radiation
INTEGER                      , INTENT(OUT) :: use_adp
REAL (8)                     , INTENT(OUT) :: r_max
REAL (8)                     , INTENT(OUT) :: r_step
INTEGER                      , INTENT(OUT) :: lrho0
REAL (8)                     , INTENT(OUT) :: rho0
REAL (8)                     , INTENT(OUT) :: corr_lin
REAL (8)                     , INTENT(OUT) :: corr_quad
INTEGER                      , INTENT(OUT) :: period
INTEGER                      , INTENT(OUT) :: exact
REAL (8)                     , INTENT(OUT) :: weight
INTEGER                      , INTENT(OUT) :: finite
REAL (8)                     , INTENT(OUT) :: diameter
REAL (8)                     , INTENT(OUT) :: qmax
REAL (8)                     , INTENT(OUT) :: qbroad
REAL (8)                     , INTENT(OUT) :: qdamp
!
CALL dlink (ano, lambda, rlambda, renergy, l_energy, &
            diff_radiation, diff_power) 
radiation  = pdf_radiation - 1   ! Python counts 0 to 2
IF(pdf_gauss) THEN
   use_adp = 1
ELSE
   use_adp = 0
ENDIF
r_max  = pdf_rmax
r_step = pdf_deltaru
lrho0   = 0
IF(pdf_lrho0) lrho0 = 1
rho0   = pdf_rho0
corr_lin  = pdf_gamma
corr_quad = pdf_delta
IF(chem_period(1)) THEN
   period = 1
ELSE
   period = 0
ENDIF
IF(pdf_lexact) THEN
   exact = 1
ELSE
   exact = 0
ENDIF
weight   = pdf_scale
finite   = pdf_finite
diameter = pdf_sphere
!
qmax = pdf_qmax
qbroad = pdf_qalp
qdamp  = pdf_sigmaq
!
END SUBROUTINE discus_get_pdf
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
SUBROUTINE discus_calc_powder(line)
!
!  A first interface that allows to calculate a Powder from python via
!  suite.fourier( python_string )
!  where python string is any Powder command.
!
USE powder
USE prompt_mod
IMPLICIT NONE
!
CHARACTER(LEN=*), INTENT(IN) :: line
!
IF( .NOT. lsetup_done) CALL initialize_suite   ! Do we need to initialize?
!
linteractive=.FALSE.    ! Tell get_cmd to get input from input_gui
CALL discus_prae        ! Switch to discus section
input_gui = line        ! copy the input line to the automatic command line
CALL do_powder          ! Call the actual task at hand
CALL back_to_suite      ! Go back to the suite
linteractive=.TRUE.     ! Tell get_cmd to read input from standard I/O
!
END SUBROUTINE discus_calc_powder
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
SUBROUTINE discus_calc_pdf(line)
!
!  A first interface that allows to calculate a Powder from python via
!  suite.fourier( python_string )
!  where python string is any Powder command.
!
USE pdf_menu
USE prompt_mod
IMPLICIT NONE
!
CHARACTER(LEN=*), INTENT(IN) :: line
!
IF( .NOT. lsetup_done) CALL initialize_suite   ! Do we need to initialize?
!
linteractive=.FALSE.    ! Tell get_cmd to get input from input_gui
CALL discus_prae        ! Switch to discus section
input_gui = line        ! copy the input line to the automatic command line
CALL pdf                ! Call the actual task at hand
CALL back_to_suite      ! Go back to the suite
linteractive=.TRUE.     ! Tell get_cmd to read input from standard I/O
!
END SUBROUTINE discus_calc_pdf
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
SUBROUTINE discus_run_save(line)
!
!  A first interface that allows to save a structure   from python via
!  suite.save( python_string )
!  where python string is any Save command.
!
USE save_menu
USE prompt_mod
IMPLICIT NONE
!
CHARACTER(LEN=*), INTENT(IN) :: line
CHARACTER(LEN=1) :: zeile
INTEGER          :: length 
!
IF( .NOT. lsetup_done) CALL initialize_suite   ! Do we need to initialize?
!
linteractive=.FALSE.    ! Tell get_cmd to get input from input_gui
CALL discus_prae        ! Switch to discus section
input_gui = line        ! copy the input line to the automatic command line
zeile = ' '
length = 0
CALL save_struc(zeile,length)! Call the actual task at hand
CALL back_to_suite      ! Go back to the suite
linteractive=.TRUE.     ! Tell get_cmd to read input from standard I/O
!
END SUBROUTINE discus_run_save
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
SUBROUTINE discus_get_save(savefile, file_length, sel_atom, w_scat, w_ncell, &
     w_mol, &
     w_dom, w_obj, w_gen, w_sym, w_start, w_end, p_all, p_nor, p_mol,  &
     p_dom, p_out, p_ext, p_int, p_lig)
!
USE ISO_C_BINDING
USE discus_save_mod
!
IMPLICIT NONE
!
CHARACTER (LEN=200), INTENT(OUT) :: savefile
INTEGER          , INTENT(OUT) :: file_length
INTEGER          , INTENT(OUT) :: sel_atom
INTEGER          , INTENT(OUT) :: w_scat
INTEGER          , INTENT(OUT) :: w_ncell
INTEGER          , INTENT(OUT) :: w_mol
INTEGER          , INTENT(OUT) :: w_dom
INTEGER          , INTENT(OUT) :: w_obj
INTEGER          , INTENT(OUT) :: w_gen
INTEGER          , INTENT(OUT) :: w_sym
INTEGER          , INTENT(OUT) :: w_start
INTEGER          , INTENT(OUT) :: w_end
INTEGER          , INTENT(OUT) :: p_all
INTEGER          , INTENT(OUT) :: p_nor
INTEGER          , INTENT(OUT) :: p_mol
INTEGER          , INTENT(OUT) :: p_dom
INTEGER          , INTENT(OUT) :: p_out
INTEGER          , INTENT(OUT) :: p_ext
INTEGER          , INTENT(OUT) :: p_int
INTEGER          , INTENT(OUT) :: p_lig
!
savefile = sav_file
file_length = LEN_TRIM(savefile)
sel_atom = log2int(sav_sel_atom)
w_scat   = log2int(sav_w_scat)
w_ncell  = log2int(sav_w_ncell)
w_mol    = log2int(sav_w_mole )
w_dom    = log2int(sav_w_doma )
w_obj    = log2int(sav_w_obje )
w_gen    = log2int(sav_w_gene )
w_sym    = log2int(sav_w_symm )
w_start  = sav_start
w_end    = sav_end
p_nor    = prop2int(sav_sel_prop,0)
p_mol    = prop2int(sav_sel_prop,0)
p_dom    = prop2int(sav_sel_prop,0)
p_out    = prop2int(sav_sel_prop,0)
p_ext    = prop2int(sav_sel_prop,0)
p_int    = prop2int(sav_sel_prop,0)
p_lig    = prop2int(sav_sel_prop,0)
!
END SUBROUTINE discus_get_save
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
SUBROUTINE discus_output(line, four_last_type)
!
!  A first interface that allows to write the Fourier output 
!  suite.discus_output( python_string )
!  where python string is any Output command.
!
USE diffuse_mod
USE output_menu
USE prompt_mod
IMPLICIT NONE
!
CHARACTER(LEN=*), INTENT(IN) :: line
INTEGER         , INTENT(IN) :: four_last_type
!
LOGICAL :: linverse
!
IF( .NOT. lsetup_done) CALL initialize_suite   ! Do we need to initialize?
!
linteractive=.FALSE.    ! Tell get_cmd to get input from input_gui
CALL discus_prae        ! Switch to discus section
input_gui = line        ! copy the input line to the automatic command line
linverse=.FALSE.
IF(four_last_type==PATT_1D .or. &
   four_last_type==PATT_2D .or. &
   four_last_type==PATT_3D     ) THEN
   linverse=.TRUE.
ENDIF
CALL do_niplps(linverse)! Call the actual task at hand
CALL back_to_suite      ! Go back to the suite
linteractive=.TRUE.     ! Tell get_cmd to read input from standard I/O
!
END SUBROUTINE discus_output
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
SUBROUTINE discus_run_sro(line)
!
!  A first interface that allows to save a structure   from python via
!  suite.save( python_string )
!  where python string is any Save command.
!
USE mmc_menu
USE prompt_mod
IMPLICIT NONE
!
CHARACTER(LEN=*), INTENT(IN) :: line
!
IF( .NOT. lsetup_done) CALL initialize_suite   ! Do we need to initialize?
!
linteractive=.FALSE.    ! Tell get_cmd to get input from input_gui
CALL discus_prae        ! Switch to discus section
input_gui = line        ! copy the input line to the automatic command line
CALL mmc()              ! Call the actual task at hand
CALL back_to_suite      ! Go back to the suite
linteractive=.TRUE.     ! Tell get_cmd to read input from standard I/O
!
END SUBROUTINE discus_run_sro
!
!________KUPLOT_________________________________________________________________
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
SUBROUTINE kuplot_load(line)
!
!  A first interface that allows to load a data file from python via
!  suite.kuplot_load( python_string )
!  where python string is any of: 
!       load xy, filename
!
USE prompt_mod
IMPLICIT NONE
!
CHARACTER(LEN=*), INTENT(IN) :: line
INTEGER                      :: length
!
IF( .NOT. lsetup_done) CALL initialize_suite   ! Do we need to initialize?
!
linteractive=.FALSE.    ! Tell get_cmd to get input from input_gui
CALL kuplot_prae        ! Switch to discus section
length = LEN_TRIM(line)
CALL do_load(line, length) ! Call the actual task at hand
CALL back_to_suite      ! Go back to the suite
linteractive=.TRUE.     ! Tell get_cmd to read input from standard I/O
!
END SUBROUTINE kuplot_load
END MODULE suite
