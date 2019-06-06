SUBROUTINE file_kdo(line, ilen)
!
USE blanks_mod
USE class_macro_internal
USE envir_mod
USE errlist_mod
USE doact_mod
USE macro_mod
USE precision_mod
USE prompt_mod
!
IMPLICIT NONE
!
INTEGER, PARAMETER ::  maxw = MAC_MAX_PARA + 1
!
CHARACTER (LEN=*),  INTENT(INOUT) :: line
INTEGER          ,  INTENT(INOUT) :: ilen
!
CHARACTER (LEN=1024), DIMENSION(1:MAXW) :: cpara
CHARACTER (LEN=1024)                    :: filename, string
!
INTEGER             , DIMENSION(1:MAXW) :: lpara
INTEGER                                 :: ianz, i
LOGICAL                                 :: fileda
REAL(KIND=PREC_DP)  , DIMENSION(1:MAXW) :: werte
!
!
INTEGER, PARAMETER :: imc   = 63
!
CHARACTER (LEN=1024), DIMENSION(:), ALLOCATABLE  :: content
CHARACTER (LEN=1024)  :: macrofile
INTEGER               :: length, file_length
INTEGER               :: iline
INTEGER               :: iseof
INTEGER               :: istatus
LOGICAL               :: is_stored
INTEGER               :: lslash          ! position os slash in filename
!
INTEGER   len_str
!
is_stored = .false.                     ! assume macro does not exist in storage
IF(macro_level==0 .AND. .NOT.lmakro) THEN
   sprompt = prompt                     ! Store prompt at top macro level start
ENDIF
macro_level = macro_level + 1
CALL build_macro_name(line, ilen, filename, MAXW, ianz, cpara, lpara, werte)
!
!  Copy filename, if no '/' within prepend with current directory
!
lslash = INDEX ( filename , '/' )
macrofile = ' '
IF ( lslash == 0 ) THEN       ! No slash in filename
   macrofile = current_dir(1:current_dir_l) // filename(1:LEN_TRIM(filename)) ! Copy filename
ELSE
   macrofile = filename       ! Copy filename
ENDIF
!
!  Let's first test if macro is stored internally
!
IF(ASSOCIATED(macro_root)) THEN       ! We do have macros in storage, test for existence
   CALL macro_find_node(macro_root, macrofile, macro_temp, ier_num)
   IF(ier_num == 0 ) THEN             ! Found the macro file in storage
      is_stored = .true.
   ELSE
      length = len_str(macrofile)
!
!              File not found, test with appended mac
!
      IF(macrofile(length-3:length) /= '.mac' ) THEN
!                 macrofile = macro_temp%macrofile(1:length) // '.mac'
         macrofile =            macrofile(1:length) // '.mac'
         CALL macro_find_node(macro_root, macrofile, macro_temp, ier_num)
         IF(ier_num == 0 ) THEN             ! Found the macro file in storage
            is_stored = .true.
         ELSE
            ALLOCATE(macro_temp, STAT=istatus)    ! Allocate a temporary macro storage
            macro_counter = macro_counter + 1
            macro_temp%number = macro_counter
            NULLIFY(macro_temp%before)            ! None before and after
            NULLIFY(macro_temp%after)
            macro_temp%macrofile = macrofile
            CALL macro_add_node(macro_root, macro_temp)  ! Add to storage
            ALLOCATE(macro_temp%macros,STAT=istatus)
            macro_temp%macros%macro_length = 0
            macro_temp%macros%lmacro       = .false.
         ENDIF
      ELSE    ! Macro ends on '.mac' but was not found
         ALLOCATE(macro_temp, STAT=istatus)    ! Allocate a temporary macro storage
            macro_counter = macro_counter + 1
            macro_temp%number = macro_counter
         NULLIFY(macro_temp%before)            ! None before and after
         NULLIFY(macro_temp%after)
         macro_temp%macrofile = macrofile
         CALL macro_add_node(macro_root, macro_temp)     ! Add to storage
         ALLOCATE(macro_temp%macros,STAT=istatus)
         macro_temp%macros%macro_length = 0
         macro_temp%macros%lmacro       = .false.
      ENDIF
   ENDIF
ELSE           ! No internal storage yet, make new storage, and add
   CALL inquire_macro_name(fileda, filename)  ! We need to locate the macro on the disk
   file_length = len_str(filename)
   IF(fileda) THEN   ! FILE EXISTS make storage
!
      ALLOCATE(macro_root, STAT=istatus)    ! Allocate a temporary macro storage
      macro_counter = macro_counter + 1
      macro_root%number = macro_counter
      NULLIFY(macro_root%before)            ! None before and after
      NULLIFY(macro_root%after)
      macro_root%macrofile = macrofile
      ALLOCATE(macro_root%macros,STAT=istatus)
      macro_root%macros%macro_length = 0
      macro_root%macros%lmacro       = .false.
   ELSE              ! File does not exist
!              MACRO not found
      ier_num = - 12
      ier_typ = ER_MAC
      oprompt = prompt
      CALL macro_close
      IF(lblock) THEN                ! If inside do/if terminate the block
         lblock_dbg = .false.
         lblock = .false.
      ENDIF
      RETURN
   ENDIF
ENDIF
CALL no_error
!
is_new: IF(.NOT. is_stored ) THEN             ! This is a new macro
   CALL inquire_macro_name(fileda, filename)  ! We need to locate the macro on the disk
   file_length = len_str(filename)
!
   file_exist: IF (fileda) THEN               ! File exist on disk
      CALL oeffne(imc, filename, 'old')
      ALLOCATE(content(1:5000) )
      content = ' '
      iline   = 0
      readcont: DO
         READ(IMC,'(a)',IOSTAT=iseof) string
         IF ( IS_IOSTAT_END(iseof) ) EXIT readcont
         length = len_str(string)
         CALL rem_leading_bl(string,length)
         iline          = iline + 1
         content(iline) = string
      ENDDO readcont
!
!     Store macro in the structure
!
      IF(ASSOCIATED(macro_temp)) THEN
         CALL macro_temp%macros%alloc_arrays(iline)
         CALL macro_temp%macros%set_macro   (iline,content)
      ELSE
         CALL macro_root%macros%alloc_arrays(iline)
         CALL macro_root%macros%set_macro   (iline,content)
      ENDIF
      CLOSE(imc)
!
      DEALLOCATE(content)
   ELSE file_exist
!       MACRO not found
      ier_num = - 12
      ier_typ = ER_MAC
      oprompt = prompt
      CALL macro_close
      IF(lblock) THEN                ! If inside do/if terminate the block
         lblock_dbg = .false.
         lblock = .false.
      ENDIF
      RETURN
   ENDIF file_exist
ELSE is_new
!
ENDIF is_new
!
IF(macro_level == 1 ) THEN                 ! Top level, start execution tree
!
   ALLOCATE(mac_tree_root, STAT=istatus)      ! Allocate next node
   macro_tree_co = macro_tree_co + 1
   mac_tree_root%number = macro_tree_co
   NULLIFY(mac_tree_root%kid)                 ! Currently no kid
   NULLIFY(mac_tree_root%sib)                 ! Root will never have siblings
   mac_tree_root%params   = ' '               ! Initialize parameters
   mac_tree_root%lparams  = 0  
   DO I=1,ianz-1
      mac_tree_root%params(i)  = cpara(i+1)(1:lpara(i+1))   ! store parameters
      mac_tree_root%lparams(i) = lpara(i+1)
   ENDDO
   WRITE(mac_tree_root%params(0), '(I5)') ianz-1 ! Store "$0"
   mac_tree_root%lparams(0) = 5               ! Length of "$0"
   mac_tree_root%nparams = ianz - 1           ! Store number of parameters
   mac_tree_root%current = 0                  ! Currently in line 0
   mac_tree_root%level   = macro_level        ! Currently at depth macro_level
   IF(ASSOCIATED(macro_temp)) THEN
      mac_tree_root%active => macro_temp         ! active macro is currently loaded macro
   ELSE
      mac_tree_root%active => macro_root         ! active macro is currently loaded macro
   ENDIF
!
   NULLIFY(mac_tree_root%parent)           ! This one has no parent, as top level
   mac_tree_active  => mac_tree_root       ! Point to currently active macro
   mac_tree_tail    => mac_tree_root       ! Point to last macro
   lmakro = .true.
   lmakro_error = .FALSE.                  ! Start with macro termination error off
ELSE
!
   ALLOCATE(mac_tree_temp, STAT=istatus)      ! Allocate next node
   macro_tree_co = macro_tree_co + 1
   mac_tree_temp%number = macro_tree_co
   NULLIFY(mac_tree_temp%kid)                 ! Currently no kid
   NULLIFY(mac_tree_temp%sib)                 ! Currently no siblings
   mac_tree_temp%params   = ' '               ! Initialize parameters
   mac_tree_temp%lparams  = 0  
   DO I=1,ianz-1
      mac_tree_temp%params(i)  = cpara(i+1)(1:lpara(i+1))   ! store parameters
      mac_tree_temp%lparams(i) = lpara(i+1)
   ENDDO
   WRITE(mac_tree_temp%params(0), '(I5)') ianz-1 ! Store "$0"
   mac_tree_temp%lparams(0) = 5               ! Length of "$0"
   mac_tree_temp%nparams = ianz - 1           ! Store number of parameters
   mac_tree_temp%current = 0                  ! Currently in line 0
   mac_tree_temp%level   = macro_level        ! Currently at depth macro_level
!
   IF(ASSOCIATED(macro_temp)) THEN
      mac_tree_temp%active => macro_temp         ! active macro is currently loaded macro
   ELSE
      mac_tree_temp%active => macro_root         ! active macro is currently loaded macro
   ENDIF
!
   IF(associated(mac_tree_active%kid)) THEN
      mac_tree_srch => mac_tree_active%kid
      DO WHILE(ASSOCIATED(mac_tree_srch%sib))
         mac_tree_srch => mac_tree_srch%sib
      ENDDO
      mac_tree_srch%sib => mac_tree_temp   ! Store new macro as sibling of active macro
   ELSE
      mac_tree_active%kid  => mac_tree_temp   ! Store new macro as kid of active macro
   ENDIF
   mac_tree_temp%parent => mac_tree_active ! Store parent of current macro
   mac_tree_active      => mac_tree_temp   ! Point to currently active macro
   mac_tree_tail        => mac_tree_temp   ! Point to last macro
ENDIF
!
IF (prompt.ne.'macro ') oprompt = prompt
!
END SUBROUTINE file_kdo
!*****7*****************************************************************
SUBROUTINE build_macro_name(line, ilen, filename, MAXW, &
           ianz, cpara, lpara, werte)
!
USE build_name_mod
USE envir_mod
USE errlist_mod
USE get_params_mod
USE precision_mod
USE prompt_mod
!
IMPLICIT NONE
!
CHARACTER (LEN=*), INTENT(IN   )  :: line
INTEGER          , INTENT(IN   )  :: ilen
CHARACTER (LEN=*), INTENT(OUT  )  :: filename
INTEGER          , INTENT(IN   )  :: MAXW
INTEGER          , INTENT(OUT  )  :: ianz
CHARACTER (LEN=*), DIMENSION(1:MAXW),INTENT(OUT  )  :: cpara
INTEGER          , DIMENSION(1:MAXW),INTENT(OUT  )  :: lpara
REAL(KIND=PREC_DP),DIMENSION(1:MAXW),INTENT(OUT  )  :: werte
!
CHARACTER(LEN=1024)      :: string
INTEGER                  :: i
INTEGER                  :: ip
INTEGER                  :: length
!
!     --Get filename from command line and string for parameters
!
string = line
ip = INDEX (string, ' ')
string (ip:ip) = ','
length = -IABS(ilen)
CALL get_params (string, ianz, cpara, lpara, maxw, length)
!
!     --Try to build filename
!
CALL do_build_name (ianz, cpara, lpara, werte, maxw, 1)
build: IF (ier_num == 0) THEN
   filename = cpara(1)(1:lpara(1))
   DO i=1, ianz
      lpara(i) = LEN_TRIM(cpara(i))   ! Remove trailing blanks from parameters
   ENDDO
ELSE build
!       ungueltiges MAKRO
   ier_num = - 13
   ier_typ = ER_MAC
ENDIF build
!
END SUBROUTINE build_macro_name
!
!*****7*****************************************************************
!
SUBROUTINE inquire_macro_name(fileda, infile)
!
USE envir_mod
USE errlist_mod
USE prompt_mod
!
IMPLICIT NONE
!
LOGICAL          , INTENT(OUT  )  :: fileda
CHARACTER (LEN=*), INTENT(INOUT)  :: infile
!
CHARACTER(LEN=1024)      :: ldir            ! local directory
CHARACTER(LEN=1024)      :: string
CHARACTER(LEN=1024)      :: filename
INTEGER                  :: ldir_length
INTEGER                  :: filename_length ! length of filename string
INTEGER                  :: infile_length   ! length of filename string
INTEGER                  :: lc
INTEGER                  :: lslash          ! position os slash in filename
!
INTEGER, EXTERNAL :: len_str
!
ldir        = current_dir       ! Build current file name
ldir_length = current_dir_l
!
build: IF (ier_num == 0) THEN
!
!       Try to open the macro file
!
!
!-----   Try filename as is
!
!           filename = cpara (1) (1:lpara (1) )
   infile_length = len_str(infile)
   filename      = infile(1:infile_length)
   INQUIRE (file = filename, exist = fileda)
!
!-----   Try filename as is with appended '.mac'
!
   IF (.NOT.fileda) THEN
!              filename = cpara (1) (1:lpara (1) ) //'.mac'
      filename = infile(1:infile_length) // '.mac'
      INQUIRE (file = filename, exist = fileda)
   ENDIF
   filename_length = len_str(filename)
!
!-----   If file is found, convert to absolute path
!
   IF ( fileda ) THEN
      lslash = INDEX ( filename , '/' )
      IF ( lslash == 0 ) THEN       ! No slash in filename
         IF(ldir(ldir_length:ldir_length) == '/') THEN
            filename = ldir(1:ldir_length) // filename(1:filename_length)
         ELSE
            filename = ldir(1:ldir_length) // '/' // filename(1:filename_length)
         ENDIF
      ENDIF
   ENDIF
!
!-----   Try local directory with appended '.mac'
!
   IF (.NOT.fileda) THEN
      string = ' '
      lc = 1
      CALL do_chdir (string, lc, .false.)
      IF (string (lc:lc) .ne.'/') THEN
         lc = lc + 1
         string (lc:lc) = '/'
      ENDIF
!           filename = string (1:lc) //cpara (1) (1:lpara (1) ) //&
!           '.mac'
         filename = string (1:lc) //infile(1:infile_length) // '.mac'
         INQUIRE (file = filename, exist = fileda)
!
!-----   Try local directory
!
         IF (.NOT.fileda) THEN
!              filename = string (1:lc) //cpara (1) (1:lpara (1) )
            filename = string (1:lc) //infile(1:infile_length)
            INQUIRE (file = filename, exist = fileda)
         ENDIF
      ENDIF
!
!-----   Try users system directory with appended '.mac'
!
      IF (.NOT.fileda) THEN
!           filename = umac_dir (1:umac_dir_l) //cpara (1)        &
!           (1:lpara (1) ) //'.mac'
         filename = umac_dir (1:umac_dir_l) //infile(1:infile_length) // '.mac'
         INQUIRE (file = filename, exist = fileda)
      ENDIF
!
!-----   Try users system directory
!
      IF (.NOT.fileda) THEN
!           filename = umac_dir (1:umac_dir_l) //cpara (1)        &
!           (1:lpara (1) )
         filename = umac_dir (1:umac_dir_l) //infile(1:infile_length)
         INQUIRE (file = filename, exist = fileda)
      ENDIF
!
!-----   Try system directory with appended '.mac'
!
      IF (.NOT.fileda) THEN
!           filename = mac_dir (1:mac_dir_l) //cpara (1) (1:lpara &
!           (1) ) //'.mac'
         filename = mac_dir (1:mac_dir_l) //infile(1:infile_length) // '.mac'
         INQUIRE (file = filename, exist = fileda)
      ENDIF
!
!-----   Try system directory
!
      IF (.NOT.fileda) THEN
!           filename = mac_dir (1:mac_dir_l) //cpara (1) (1:lpara &
!           (1) )
         filename = mac_dir (1:mac_dir_l) //infile(1:infile_length)
         INQUIRE (file = filename, exist = fileda)
      ENDIF
   infile = filename(1:len_str(filename))
ELSE build
!       ungueltiges MAKRO
   ier_num = - 13
   ier_typ = ER_MAC
ENDIF build
!
END SUBROUTINE inquire_macro_name
!
!*****7*****************************************************************
!
SUBROUTINE macro_read (line, laenge)
!-
!     Reads a single line from the current macro storage
!+
USE ber_params_mod
USE charact_mod
USE doact_mod
USE do_eval_mod
USE errlist_mod
USE get_params_mod
USE macro_mod
USE class_macro_internal
USE precision_mod
USE prompt_mod
IMPLICIT none
!
!
CHARACTER(LEN=*), INTENT(INOUT) :: line
INTEGER         , INTENT(INOUT) :: laenge
!
CHARACTER(LEN=1024) :: zeile
CHARACTER(LEN=1024), DIMENSION(1) :: string
INTEGER :: ndol, nexcl, nquote1, nquote2
INTEGER :: lpar
INTEGER :: n_par
INTEGER :: lx, nx, x, lll, sdol
INTEGER :: il
INTEGER, DIMENSION(1) :: lstring
LOGICAL :: lnum
REAL(KIND=PREC_DP)   , DIMENSION(1) :: r_par
!
INTEGER, EXTERNAL :: len_str
LOGICAL, EXTERNAL :: str_comp
!
ier_num = 0
ier_typ = ER_NONE
mac_tree_active%current = mac_tree_active%current + 1
IF(mac_tree_active%current > mac_tree_active%active%macros%macro_length) THEN
   IF(.NOT. ASSOCIATED(mac_tree_active%parent)) THEN  ! Got back to the top 

      lmakro = .false.
      macro_level = 0
      CALL macro_close
   ELSE
      mac_tree_active => mac_tree_active%parent
!            DEALLOCATE(mac_tree_active%kid)
      macro_level = macro_level - 1
   ENDIF
   RETURN
ENDIF
line                  = mac_tree_active%active%macros%macro_line(mac_tree_active%current)
laenge = len_str (line)
IF (laenge == 0) THEN
   line = ' '
   laenge = 0
   RETURN
ENDIF
!
nocomment: IF (line(1:1) /= '#' .AND. line (1:1) /=  '!' .AND. laenge /= 0) THEN
   ndol    = 0
   nexcl   = 0
   nquote1 = 0
   nquote2 = 0
   zeile = line
   nexcl   = INDEX(zeile (1:laenge) , '!', .TRUE.)   ! Find last '!'
   nquote1 = INDEX(zeile (1:laenge) , '"', .TRUE.)   ! Find last '"'
   IF(nexcl > 0 .AND. nexcl > nquote1) THEN
      laenge = nexcl -1
   ENDIF
!
!     If inside a macro, check for parameter substitution
!
   ndol = INDEX (zeile (1:laenge) , '$')
   IF (ndol>   laenge) THEN
      ndol = 0
   ENDIF
!
!     Replace all '$'-parameters
!
   sub_param: DO WHILE (ndol.ne.0)
!
!------ --Determine length of the numerical parameter string
!       i.e.: '1', '12'
!
      lx   = 0
      nx   = ndol + lx + 1
      x    = IACHAR (zeile (nx:nx) )
      lnum = zero <= x.AND.x <= nine
      DO WHILE (lnum.AND.nx <  laenge)
         lx   = lx + 1
         nx   = ndol + lx + 1
         x    = IACHAR (zeile (nx:nx) )
         lnum = zero <= x.AND.x <= nine
      ENDDO
      IF (nx <  laenge) THEN
         nx = nx - 1
      ELSEIF (nx == laenge.AND..NOT.lnum) THEN
         nx = nx - 1
      ENDIF
!
!     --Read parameter number and substitute rest of string
!
      IF (nx >= ndol + 1) THEN
         string = zeile (ndol + 1:nx)
         lstring = nx - (ndol + 1) + 1
         CALL ber_params (1, string, lstring, r_par, 1)
      ELSE
         ier_num = - 12
         ier_typ = ER_FORT
      ENDIF
!DBG        read(zeile(ndol+1:nx),*) n_par
      no_errorA: IF (ier_num == 0) THEN
         n_par = nint (r_par(1))
         IF(n_par <= mac_tree_active%nparams) THEN
            lpar = mac_tree_active%lparams(n_par)
            line = ' '
            IF (ndol>   1) THEN
               line (1:ndol - 1) = zeile (1:ndol - 1)
            ENDIF
!
            IF(.NOT.lblock_read .AND. mac_tree_active%params(n_par)(1:6)=='value(') THEN   ! Evaluate the parametere
               zeile = mac_tree_active%params(n_par)(7:lpar-1)
               lpar   = lpar-2
               CALL do_eval(zeile, lpar, .FALSE.)
               IF(ier_num/=0) RETURN
               line (ndol:ndol+lpar-1) = zeile(1:lpar)
            ELSE
               line (ndol:ndol+lpar-1) = mac_tree_active%params(n_par)(1:lpar)
            ENDIF
            lll = ndol + lpar - 1
            IF (nx <  laenge) THEN
               line (ndol + lpar:ndol + lpar + laenge-nx) = zeile (  &
               nx + 1:laenge)
               lll = lll + laenge-nx
            ENDIF
            zeile = ' '
            zeile = line
            laenge = lll
            sdol = ndol + 1
            IF (sdol>   laenge) THEN
               ndol = 0
            ELSE
               ndol = INDEX (zeile (sdol:laenge) , '$')
               IF (ndol>   laenge) THEN
                  ndol = 0
               ELSEIF (ndol>   0) THEN
                  ndol = ndol + sdol - 1
               ENDIF
            ENDIF
         ELSE
            lmakro = .false.
            il     = len_str (line)
            WRITE (output_io, 1000) line (1:il)
            ier_num = - 41
            ier_typ = ER_MAC
            CALL macro_close
            RETURN
         ENDIF
      ELSE no_errorA
         ier_num = 0
         ier_typ = ER_NONE
         sdol = ndol + 1
         IF (sdol>   laenge) THEN
            ndol = 0
         ELSE
            ndol = INDEX (zeile (sdol:laenge) , '$')
            IF (ndol>   laenge) THEN
               ndol = 0
            ELSEIF (ndol>   0) THEN
               ndol = ndol + sdol - 1
            ENDIF
         ENDIF
      ENDIF no_errorA
   ENDDO sub_param
ENDIF nocomment
!
!     line read, return to calling routine
!
il = len_str (line)
IF (prompt_status.ne.PROMPT_OFF) THEN
   IF (il>   0) THEN
      WRITE (output_io, 1000) line (1:il)
   ELSE
      WRITE (output_io, 1000) ' '
   ENDIF
ENDIF
!
!     Only if the string has length longer than zero
!
IF(il>   0) THEN
!
!     Check for 'stop' command, unless we are reading a block structure
!
   IF(str_comp(line, 'stop', 4, il, 4) .AND..NOT.lblock_read) THEN
      WRITE ( *, 2000) char (7)
      lmakro = .false.
      lmakro_error = .false.    ! Macro termination error off
      line = '#'
      il = 1
   ENDIF
ENDIF
!
!
 1000 FORMAT     (a)
 2000 FORMAT     (' ------ > Macro halted, continue with cont ...',a1)
END SUBROUTINE macro_read
!
!*****7****************************************************************
!
SUBROUTINE macro_terminate
!
USE class_macro_internal
USE macro_mod
!
IMPLICIT NONE
!
!IF(mac_tree_active%current > mac_tree_active%active%macros%macro_length) THEN
   IF(.NOT. ASSOCIATED(mac_tree_active%parent)) THEN  ! Got back to the top 
      lmakro = .false.
      macro_level = 0
      CALL macro_close
   ELSE
      mac_tree_active => mac_tree_active%parent
      macro_level = macro_level - 1
   ENDIF
!   RETURN
!ENDIF
!
END SUBROUTINE macro_terminate
!
!*****7****************************************************************
!
SUBROUTINE macro_close
!-
!     Closes the macro file, switches macro status off and sets the
!     macro level back to zero.
!     The macro tree is deallocated.
!     In an interactive session (promot /= redirect ) the stored 
!     macros are deallocated. This allows the user to modify a macro
!     and run the modified version
!+
USE class_macro_internal
USE macro_mod
USE mpi_slave_mod
USE errlist_mod
USE prompt_mod
IMPLICIT none
!
IF(mpi_is_slave) THEN
   RETURN
ENDIF
!
IF(prompt_status/=PROMPT_OFF) THEN
   WRITE(output_io,*) ' '
ENDIF
!     
lmakro = .false.
IF(ier_num/=0) lmakro_error = .TRUE.     ! Macro terminated with error
macro_level = 0
IF(ASSOCIATED(mac_tree_root)) THEN                ! We have stored macros
   CALL macro_rem_tree(mac_tree_root)
   NULLIFY(mac_tree_root)                      ! Clear pointer status
   NULLIFY(mac_tree_temp)
   NULLIFY(mac_tree_tail)
   NULLIFY(mac_tree_srch)
   macro_tree_co = 0
ENDIF
!
IF(ASSOCIATED(macro_root)) THEN                ! We have stored macros
   CALL macro_rem_all(macro_root)
   macro_counter = 0
ENDIF
!
NULLIFY(macro_root)                      ! Clear pointer status
NULLIFY(macro_temp)
NULLIFY(mac_tree_tail)
!
END SUBROUTINE macro_close
!
!*****7*****************************************************************
!
SUBROUTINE macro_close_mpi(first_mac, mac_l)
!
! close the tree associated with currently running mpi slave
!
USE class_macro_internal
USE macro_mod
USE errlist_mod
USE prompt_mod
!
IMPLICIT NONE
!
CHARACTER (LEN=*), INTENT(IN) :: first_mac
INTEGER          , INTENT(IN) :: mac_l
!
IF(ASSOCIATED(mac_tree_root)) THEN                ! We have stored macros
   CALL macro_rem_tree(mac_tree_root)
   NULLIFY(mac_tree_root)                      ! Clear pointer status
   NULLIFY(mac_tree_temp)
   NULLIFY(mac_tree_tail)
   NULLIFY(mac_tree_srch)
   lmakro = .false.
   macro_tree_co = 0
   macro_level = 0
   IF(ier_num/=0) lmakro_error = .TRUE.     ! Macro terminated with error
ENDIF
!
END SUBROUTINE macro_close_mpi
!
!*****7*****************************************************************
!
      SUBROUTINE macro_continue (zeile, lcomm)
!-
!     Continues the macro file, that had been interupted for
!     debugging purposes
!+
      USE doact_mod
      USE errlist_mod
      USE get_params_mod
      USE class_macro_internal
      USE macro_mod
      USE prompt_mod
      IMPLICIT none
!
!
      CHARACTER (LEN= * ), INTENT(INOUT) :: zeile
      INTEGER            , INTENT(INOUT) :: lcomm
      INTEGER maxw
      PARAMETER (maxw = 1)
      CHARACTER(1024) cpara (maxw)
      CHARACTER(LEN=40)  :: cprompt
      INTEGER lpara (maxw)
      INTEGER ianz
!
      INTEGER len_str
      LOGICAL str_comp
!
      CALL get_params (zeile, ianz, cpara, lpara, maxw, lcomm)
!
!     --No parameter for macro or block structures
!
      IF (ianz == 0) THEN
!
!     --Only active inside macros stopped by a previous 'stop' command
!
!        IF (mac_level>   0) THEN
         IF (macro_level>   0) THEN
            lmakro = .true.
            lmakro_error = .FALSE.
            IF (prompt.ne.'macro ') oprompt = prompt
            prompt = 'macro '
         ENDIF
!
!     Only active inside macros stopped by a previous 'stop' command
!
         IF (lblock_dbg) THEN
            lblock_dbg = .false.
            lblock = .true.
         ENDIF
!
!     --One parameter and '<pname>' command
!
      ELSEIF (ianz == 1) THEN
         IF (str_comp (cpara (1), pname, 3, lpara (1), len_str (pname) )&
         ) THEN
            cprompt = prompt   ! Remember current prompt, as macro close 
                               ! goes all the way back....
            CALL macro_close
            prompt = cprompt   ! Set current prompt as active prompt
            lblock_dbg = .false.
            lblock = .false.
         ELSEIF(str_comp(cpara(1), 'suite', 3, lpara (1), 5)) THEN
            CALL macro_close
            lblock_dbg = .false.
            lblock = .false.
            IF(pname /= 'suite') THEN
               zeile = 'EXIT'
               lcomm = 4
            ENDIF
         ELSE
            ier_num = - 6
            ier_typ = ER_COMM
            CALL macro_close
         ENDIF
      ELSE
         ier_num = - 6
         ier_typ = ER_COMM
         CALL macro_close
      ENDIF
!
      END SUBROUTINE macro_continue
!
!*****7*****************************************************************
!
SUBROUTINE test_macro(line,ilen, numpar)
!
!  Tests how many macro parameters a macro needs
!
USE ber_params_mod
USE charact_mod
USE errlist_mod
USE macro_mod
USE precision_mod
!
IMPLICIT NONE
!
CHARACTER(LEN=*), INTENT(IN)  :: line
INTEGER         , INTENT(IN)  :: ilen
INTEGER         , INTENT(OUT) :: numpar
!
INTEGER, PARAMETER ::  maxw = MAC_MAX_PARA + 1
INTEGER, PARAMETER ::  imc  = 63
!
!
CHARACTER (LEN=1024), DIMENSION(1:MAXW) :: cpara
CHARACTER (LEN=1024)                    :: filename, zeile
CHARACTER (LEN=1024), DIMENSION(1:MAXW) :: string
!
INTEGER             , DIMENSION(1:MAXW) :: lpara
INTEGER                                 :: ianz, length
INTEGER                                 :: ndol, nexcl, nquote1,nquote2, nx, lx, x
INTEGER             , DIMENSION(1:MAXW) :: lstring
LOGICAL                                 :: lnum
INTEGER                                 :: iseof
LOGICAL                                 :: fileda
REAL(KIND=PREC_DP)  , DIMENSION(1:MAXW) :: werte
REAL(KIND=PREC_DP)  , DIMENSION(1:MAXW) :: r_par
!
numpar = 0                                 ! Assume no parameters are required
CALL build_macro_name(line, ilen, filename, MAXW, ianz, cpara, lpara, werte)
CALL inquire_macro_name(fileda, filename)  ! We need to locate the macro on the disk
IF(fileda) THEN
   CALL oeffne(imc, filename, 'old')
   readcont: DO                            ! Read all lines from macro
      READ(IMC,'(a)',IOSTAT=iseof) string(1)
      IF ( IS_IOSTAT_END(iseof) ) EXIT readcont
      length = LEN_TRIM(string(1))
!
      nocomment: IF (string(1)(1:1) /= '#' .AND. string(1) (1:1) /=  '!' .AND. length /= 0) THEN
         ndol    = 0
         nexcl   = 0
         nquote1 = 0
         nquote2 = 0
         zeile   = string(1)
         nexcl   = INDEX(zeile (1:length) , '!', .TRUE.)   ! Find last '!'
         nquote1 = INDEX(zeile (1:length) , '"', .TRUE.)   ! Find last '"'
         IF(nexcl > 0 .AND. nexcl > nquote1) THEN
            length = nexcl -1
         ENDIF
!
!     If inside a macro, check for parameter substitution
!
         ndol = INDEX (zeile (1:length) , '$')
         IF (ndol>   length) THEN
            ndol = 0
         ENDIF
!
!     Replace all '$'-parameters
!
         sub_param: DO WHILE (ndol.ne.0)
!
!------ --Determine length of the numerical parameter string
!       i.e.: '1', '12'
!
            lx   = 0
            nx   = ndol + lx + 1
            x    = IACHAR (zeile (nx:nx) )
            lnum = zero <= x.AND.x <= nine
            DO WHILE (lnum.AND.nx <  length)
               lx   = lx + 1
               nx   = ndol + lx + 1
               x    = IACHAR (zeile (nx:nx) )
               lnum = zero <= x.AND.x <= nine
            ENDDO
            IF (nx <  length) THEN
               nx = nx - 1
            ELSEIF (nx == length.AND..NOT.lnum) THEN
               nx = nx - 1
            ENDIF
!
!     --Read parameter number and substitute rest of string
!
            IF (nx >= ndol + 1) THEN
               string(1) = zeile (ndol + 1:nx)
               lstring(1) = nx - (ndol + 1) + 1
               CALL ber_params (1, string, lstring, r_par, 1)
            ELSE
               ier_num = - 12
               ier_typ = ER_FORT
            ENDIF
!DBG        read(zeile(ndol+1:nx),*) n_par
            IF (ier_num == 0) THEN
               numpar = MAX(numpar, NINT(r_par(1)))
            ENDIF
            IF(nx+1 > length) EXIT sub_param
            zeile = zeile(nx+1:length)
            length=LEN_TRIM(zeile)
            ndol = INDEX (zeile (1:length) , '$')
         ENDDO sub_param
      ENDIF  nocomment
   ENDDO readcont                            ! Read all lines from macro
ELSE
   ier_num = - 12
   ier_typ = ER_MAC
ENDIF
CLOSE(IMC)
!
!*****7*****************************************************************
!
END SUBROUTINE test_macro
