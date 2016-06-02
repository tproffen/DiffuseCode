SUBROUTINE file_kdo(line, ilen)
!
      USE class_macro_internal
      USE envir_mod
      USE errlist_mod
      USE doact_mod
      USE macro_mod
      USE prompt_mod
!
      IMPLICIT NONE
!
      INTEGER, PARAMETER ::  maxw = MAC_MAX_PARA + 1
!
      CHARACTER (LEN=*),  INTENT(INOUT) ::  line
      INTEGER          ,  INTENT(INOUT) :: ilen
!
      CHARACTER (LEN=1024), DIMENSION(1:MAXW) :: cpara
      CHARACTER (LEN=1024)                    :: filename, string
!
      INTEGER             , DIMENSION(1:MAXW) :: lpara
      INTEGER                                 :: ianz, i
      LOGICAL                                 :: fileda
      REAL                , DIMENSION(1:MAXW) :: werte
!
!
      LOGICAL, PARAMETER :: lread = .true.
      INTEGER, PARAMETER :: imc   = 63
!
      CHARACTER (LEN=1024), DIMENSION(:), ALLOCATABLE  :: content
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
      macro_level = macro_level + 1
      CALL build_macro_name(line, ilen, fileda, filename, MAXW, ianz, cpara, lpara, werte)
!
!     Build a temporary storage
!
         ALLOCATE(macro_temp, STAT=istatus)    ! Allocate a temporary macro storage
         NULLIFY(macro_temp%before)            ! None before and after
         NULLIFY(macro_temp%after)
         IF ( istatus /= 0) THEN
            ier_num = -114
            ier_typ = ER_MAC
            macro_level = 0
!           THIS should be a routine "destroy tree"
!           IF(ASSOCIATED(mac_tree_root)) DEALLOCATE(mac_tree_root) !WARNING MEMORY LEAK
!           IF(ASSOCIATED(mac_tree_active)) DEALLOCATE(mac_tree_active) !WARNING MEMORY LEAK
            lmakro = .false.
            oprompt = prompt
            CALL macro_close
            RETURN
         ENDIF
!
!  Copy filename, if no '/' within prepend with current directory
!
         lslash = index ( filename , '/' )
         IF ( lslash == 0 ) THEN       ! No slash in filename
            macro_temp%macrofile = current_dir(1:current_dir_l) // filename       ! Copy filename
         ELSE
            macro_temp%macrofile = filename       ! Copy filename
         ENDIF
!
!  Let's first test if macro is stored internally
!
         IF(ASSOCIATED(macro_root)) THEN       ! We do have macros in storage, test for existence
            CALL macro_find_node(macro_root, macro_temp, ier_num)
            IF(ier_num == 0 ) THEN             ! Found the macro file in storage
               is_stored = .true.
            ELSE
               length = len_str(macro_temp%macrofile)
!
!              File not found, test with appended mac
!
               IF(macro_temp%macrofile(length-3:length) /= '.mac' ) THEN
                  macro_temp%macrofile = macro_temp%macrofile(1:length) // '.mac'
                  CALL macro_find_node(macro_root, macro_temp, ier_num)
                  IF(ier_num == 0 ) THEN             ! Found the macro file in storage
                     is_stored = .true.
                  ELSE
                     CALL macro_add_node(macro_root, macro_temp)  ! Add to storage
                     ALLOCATE(macro_temp%macros,STAT=istatus)
                     macro_temp%macros%macro_length = 0
                     macro_temp%macros%lmacro       = .false.
                  ENDIF
               ELSE    ! Macro ends on '.mac' but was not found
                  CALL macro_add_node(macro_root, macro_temp)     ! Add to storage
                  ALLOCATE(macro_temp%macros,STAT=istatus)
                  macro_temp%macros%macro_length = 0
                  macro_temp%macros%lmacro       = .false.
               ENDIF
            ENDIF
         ELSE           ! No internal storage yes, make new storage, and add
            CALL inquire_macro_name(fileda, filename)  ! We need to locate the macro on the disk
            file_length = len_str(filename)
            IF(fileda) THEN   ! FILE EXISTS make storage
               CALL macro_add_node(macro_root, macro_temp)
               ALLOCATE(macro_temp%macros,STAT=istatus)
               macro_temp%macros%macro_length = 0
               macro_temp%macros%lmacro       = .false.
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
         is_new: IF(.not. is_stored ) THEN             ! This is a new macro
            CALL inquire_macro_name(fileda, filename)  ! We need to locate the macro on the disk
            file_length = len_str(filename)
!
            file_exist: IF (fileda) then               ! File exist on disk
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
            CALL macro_temp%macros%alloc_arrays(iline)
            CALL macro_temp%macros%set_macro   (iline,content)
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
         ENDIF is_new
!
      ALLOCATE(mac_tree_temp, STAT=istatus)      ! Allocate next node
      NULLIFY(mac_tree_temp%kid)                 ! Currently no kid
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
      mac_tree_temp%active => macro_temp         ! active macro is currently loaded macro
!
      IF(macro_level == 1 ) THEN                 ! Top level, start execution tree
         NULLIFY(mac_tree_temp%parent)           ! This one has no parent, as top level
         mac_tree_root    => mac_tree_temp       ! root points to current
         mac_tree_active  => mac_tree_temp       ! Point to currently active macro
         mac_tree_tail    => mac_tree_temp       ! Point to last macro
         lmakro = .true.
      ELSE
         mac_tree_active%kid  => mac_tree_temp   ! Store new macro as kid of active macro
         mac_tree_temp%parent => mac_tree_active ! Store parent of current macro
         mac_tree_active      => mac_tree_temp   ! Point to currently active macro
         mac_tree_tail        => mac_tree_temp   ! Point to last macro
      ENDIF
!
!
!
      IF (prompt.ne.'macro ') oprompt = prompt
!
      END SUBROUTINE file_kdo
!*****7*****************************************************************
      SUBROUTINE build_macro_name(line, ilen, fileda, filename, MAXW, &
                 ianz, cpara, lpara, werte)
!
      USE envir_mod
      USE errlist_mod
      USE prompt_mod
!
      IMPLICIT NONE
!
      CHARACTER (LEN=*), INTENT(IN   )  :: line
      INTEGER          , INTENT(IN   )  :: ilen
      LOGICAL          , INTENT(OUT  )  :: fileda
      CHARACTER (LEN=*), INTENT(OUT  )  :: filename
      INTEGER          , INTENT(IN   )  :: MAXW
      INTEGER          , INTENT(OUT  )  :: ianz
      CHARACTER (LEN=*), DIMENSION(1:MAXW),INTENT(OUT  )  :: cpara
      INTEGER          , DIMENSION(1:MAXW),INTENT(OUT  )  :: lpara
      REAL             , DIMENSION(1:MAXW),INTENT(OUT  )  :: werte
!
      CHARACTER(LEN=1024)      :: string
      INTEGER                  :: ip
!
!     --Get filename from command line and string for parameters
!
      string = line
      ip = index (string, ' ')
      string (ip:ip) = ','
      CALL get_params (string, ianz, cpara, lpara, maxw, ilen)
!
!     --Try to build filename
!
      CALL do_build_name (ianz, cpara, lpara, werte, maxw, 1)
      build: IF (ier_num == 0) then
         filename = cpara(1)(1:lpara(1))
      ELSE build
!       ungueltiges MAKRO
         ier_num = - 13
         ier_typ = ER_MAC
      ENDIF build
!
      END SUBROUTINE build_macro_name
!*****7*****************************************************************
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
      INTEGER :: len_str
!
      ldir        = current_dir       ! Build current file name
      ldir_length = current_dir_l
!
      build: IF (ier_num.eq.0) then
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
            IF (.not.fileda) then
!              filename = cpara (1) (1:lpara (1) ) //'.mac'
               filename = infile(1:infile_length) // '.mac'
               INQUIRE (file = filename, exist = fileda)
            ENDIF
            filename_length = len_str(filename)
!
!-----   If file is found, convert to absolute path
!
            IF ( fileda ) THEN
              lslash = index ( filename , '/' )
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
            IF (.not.fileda) then
               string = ' '
               lc = 1
               CALL do_chdir (string, lc, .false.)
               IF (string (lc:lc) .ne.'/') then
                  lc = lc + 1
                  string (lc:lc) = '/'
               ENDIF
!              filename = string (1:lc) //cpara (1) (1:lpara (1) ) //&
!              '.mac'
               filename = string (1:lc) //infile(1:infile_length) // '.mac'
               INQUIRE (file = filename, exist = fileda)
!
!-----   Try local directory
!
               IF (.not.fileda) then
!                 filename = string (1:lc) //cpara (1) (1:lpara (1) )
                  filename = string (1:lc) //infile(1:infile_length)
                  INQUIRE (file = filename, exist = fileda)
               ENDIF
            ENDIF
!
!-----   Try users system directory with appended '.mac'
!
            IF (.not.fileda) then
!              filename = umac_dir (1:umac_dir_l) //cpara (1)        &
!              (1:lpara (1) ) //'.mac'
               filename = umac_dir (1:umac_dir_l) //infile(1:infile_length) // '.mac'
               INQUIRE (file = filename, exist = fileda)
            ENDIF
!
!-----   Try users system directory
!
            IF (.not.fileda) then
!              filename = umac_dir (1:umac_dir_l) //cpara (1)        &
!              (1:lpara (1) )
               filename = umac_dir (1:umac_dir_l) //infile(1:infile_length)
               INQUIRE (file = filename, exist = fileda)
            ENDIF
!
!-----   Try system directory with appended '.mac'
!
            IF (.not.fileda) then
!              filename = mac_dir (1:mac_dir_l) //cpara (1) (1:lpara &
!              (1) ) //'.mac'
               filename = mac_dir (1:mac_dir_l) //infile(1:infile_length) // '.mac'
               INQUIRE (file = filename, exist = fileda)
            ENDIF
!
!-----   Try system directory
!
            IF (.not.fileda) then
!              filename = mac_dir (1:mac_dir_l) //cpara (1) (1:lpara &
!              (1) )
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
!*****7*****************************************************************
      SUBROUTINE macro_read (line, laenge)
!-
!     Reads a single line from the current macro storage
!+
      USE charact_mod
      USE doact_mod
      USE errlist_mod
      USE macro_mod
      USE class_macro_internal
      USE prompt_mod
      IMPLICIT none
!
!
      CHARACTER ( * ) line
!
      CHARACTER(1024) zeile
      CHARACTER(1024) string
      INTEGER ndol
      INTEGER lpar
      INTEGER n_par
      INTEGER laenge
      INTEGER lx, nx, x, lll, sdol
      INTEGER il
      INTEGER lstring
      LOGICAL lnum
      REAL r_par
!
      INTEGER len_str
      LOGICAL str_comp
!
      ier_num = 0
      ier_typ = ER_NONE
      mac_tree_active%current = mac_tree_active%current + 1
      IF(mac_tree_active%current > mac_tree_active%active%macros%macro_length) THEN
         IF(.not. ASSOCIATED(mac_tree_active%parent)) THEN  ! Got back to the top 

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
      IF (laenge.eq.0) then
         line = ' '
         laenge = 0
         RETURN
      ENDIF
!
      nocomment: IF (line(1:1) /= '#' .and. line (1:1) /=  '!' .and. laenge /= 0) THEN
!
!     If inside a macro, check for parameter substitution
!
         zeile = line
         ndol = index (zeile (1:laenge) , '$')
         IF (ndol.gt.laenge) then
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
         x    = iachar (zeile (nx:nx) )
         lnum = zero.le.x.and.x.le.nine
         DO WHILE (lnum.and.nx.lt.laenge)
            lx   = lx + 1
            nx   = ndol + lx + 1
            x    = iachar (zeile (nx:nx) )
            lnum = zero.le.x.and.x.le.nine
         ENDDO
         IF (nx.lt.laenge) then
            nx = nx - 1
         ELSEIF (nx.eq.laenge.and..not.lnum) then
            nx = nx - 1
         ENDIF
!
!     --Read parameter number and substitute rest of string
!
         IF (nx.ge.ndol + 1) then
            string = zeile (ndol + 1:nx)
            lstring = nx - (ndol + 1) + 1
            CALL ber_params (1, string, lstring, r_par, 1)
         ELSE
            ier_num = - 12
            ier_typ = ER_FORT
         ENDIF
!DBG        read(zeile(ndol+1:nx),*) n_par
         IF (ier_num.eq.0) then
            n_par = nint (r_par)
               IF(n_par <= mac_tree_active%nparams) THEN
               lpar = mac_tree_active%lparams(n_par)
               line = ' '
               IF (ndol.gt.1) then
                  line (1:ndol - 1) = zeile (1:ndol - 1)
               ENDIF
!
               line (ndol:ndol+lpar-1) = mac_tree_active%params(n_par)(1:lpar)
               lll = ndol + lpar - 1
               IF (nx.lt.laenge) then
                  line (ndol + lpar:ndol + lpar + laenge-nx) = zeile (  &
                  nx + 1:laenge)
                  lll = lll + laenge-nx
               ENDIF
               zeile = ' '
               zeile = line
               laenge = lll
               sdol = ndol + 1
               IF (sdol.gt.laenge) then
                  ndol = 0
               ELSE
                  ndol = index (zeile (sdol:laenge) , '$')
                  IF (ndol.gt.laenge) then
                     ndol = 0
                  ELSEIF (ndol.gt.0) then
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
         ELSE
            ier_num = 0
            ier_typ = ER_NONE
            sdol = ndol + 1
            IF (sdol.gt.laenge) then
               ndol = 0
            ELSE
               ndol = index (zeile (sdol:laenge) , '$')
               IF (ndol.gt.laenge) then
                  ndol = 0
               ELSEIF (ndol.gt.0) then
                  ndol = ndol + sdol - 1
               ENDIF
            ENDIF
         ENDIF
         ENDDO sub_param
      ENDIF nocomment
!
!     line read, return to calling routine
!
      il = len_str (line)
      IF (prompt_status.ne.PROMPT_OFF) then
         IF (il.gt.0) then
            WRITE (output_io, 1000) line (1:il)
         ELSE
            WRITE (output_io, 1000) ' '
         ENDIF
      ENDIF
!
!     Only if the string has length longer than zero
!
      IF (il.gt.0) then
!
!     Check for 'stop' command, unless we are reading a block structure
!
         IF (str_comp (line, 'stop', 4, il, 4) .and..not.lblock_read)   &
         then
            WRITE ( *, 2000) char (7)
            lmakro = .false.
            line = '#'
            il = 1
         ENDIF
      ENDIF
      RETURN
!
!
 1000 FORMAT     (a)
 2000 FORMAT     (' ------ > Macro halted, continue with cont ...',a1)
      END SUBROUTINE macro_read
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
      USE prompt_mod
      IMPLICIT none
!
      INTEGER :: all_status
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
      prompt = oprompt
!     mac_level = 0
      macro_level = 0
IF(ASSOCIATED(mac_tree_root)) THEN
!
   mac_tree_active => mac_tree_root
!
   DO WHILE(ASSOCIATED(mac_tree_active%kid))   ! There are more macros in the tree
      mac_tree_active => mac_tree_active%kid   ! Point to kid
      CALL          mac_tree_active%parent%active%macros%finalize_macro ! Deallocate macro lines
      IF(ASSOCIATED(mac_tree_active%parent%active%macros)) &
      DEALLOCATE(mac_tree_active%parent%active%macros,STAT=all_status)         ! Remove macro storage
      DEALLOCATE(mac_tree_active%parent,stat=all_status)       ! Remove previous node
   ENDDO
   CALL mac_tree_active%active%macros%finalize_macro ! Deallocate macro lines
   DEALLOCATE(mac_tree_active%active%macros)         ! Remove macro storage
   DEALLOCATE(mac_tree_active)                 ! Finally remove current node
   NULLIFY(mac_tree_root)                      ! Clear pointer status
   NULLIFY(mac_tree_temp)
ENDIF
!
!!!IF(prompt_status/=PROMPT_REDIRECT) THEN        ! Assume an interactive session
   IF(ASSOCIATED(macro_root)) THEN             ! We have stored macros
      macro_temp => macro_root
   ENDIF
!
   IF(ASSOCIATED(macro_temp)) THEN             ! There are some macros
      DO WHILE(ASSOCIATED(macro_temp%after))   ! There are more macros in the list
      macro_temp => macro_temp%after           ! Point to next macro
         IF(ASSOCIATED(macro_temp%before)) THEN
            DEALLOCATE(macro_temp%before)      ! Remove previous macro
         ENDIF
      ENDDO
      DEALLOCATE(macro_temp)                   ! Finally remove current node
      NULLIFY(macro_root)                      ! Clear pointer status
      NULLIFY(macro_temp)
   ENDIF
!!!ENDIF
!
END SUBROUTINE macro_close
!
!*****7*****************************************************************
SUBROUTINE macro_close_mpi(first_mac, mac_l)
!
! close the tree associated with currently running mpi slave
!
USE class_macro_internal
USE macro_mod
USE prompt_mod
!
IMPLICIT NONE
!
CHARACTER (LEN=*), INTENT(IN) :: first_mac
INTEGER          , INTENT(IN) :: mac_l
INTEGER :: indx
!
mac_tree_active => mac_tree_tail   ! start at tail 
indx = INDEX(mac_tree_active%active%macrofile,first_mac(1:mac_l))

DO WHILE(indx == 0)
   mac_tree_active => mac_tree_active%parent
   DEALLOCATE(mac_tree_active%kid)
   indx = INDEX(mac_tree_active%active%macrofile,first_mac(1:mac_l))
ENDDO
IF(ASSOCIATED(mac_tree_active%parent)) THEN
   mac_tree_active => mac_tree_active%parent
   DEALLOCATE(mac_tree_active%kid)
   macro_level = mac_tree_active%level
ELSE
   DEALLOCATE(mac_tree_active)
   lmakro = .false.
   macro_level = 0
ENDIF
END SUBROUTINE macro_close_mpi
!*****7*****************************************************************
      SUBROUTINE macro_continue (zeile, lcomm)
!-
!     Continues the macro file, that had been interupted for
!     debugging purposes
!+
      USE doact_mod
      USE errlist_mod
      USE class_macro_internal
      USE macro_mod
      USE prompt_mod
      IMPLICIT none
!
!
      INTEGER maxw
      PARAMETER (maxw = 1)
      CHARACTER ( * ) zeile
      CHARACTER(1024) cpara (maxw)
      INTEGER lpara (maxw)
      INTEGER ianz, lcomm
!
      INTEGER len_str
      LOGICAL str_comp
!
      CALL get_params (zeile, ianz, cpara, lpara, maxw, lcomm)
!
!     --No parameter for macro or block structures
!
      IF (ianz.eq.0) then
!
!     --Only active inside macros stopped by a previous 'stop' command
!
!        IF (mac_level.gt.0) then
         IF (macro_level.gt.0) then
            lmakro = .true.
            IF (prompt.ne.'macro ') oprompt = prompt
            prompt = 'macro '
         ENDIF
!
!     Only active inside macros stopped by a previous 'stop' command
!
         IF (lblock_dbg) then
            lblock_dbg = .false.
            lblock = .true.
         ENDIF
!
!     --One parameter and '<pname>' command
!
      ELSEIF (ianz.eq.1) then
         IF (str_comp (cpara (1), pname, 3, lpara (1), len_str (pname) )&
         ) then
            CALL macro_close
            lblock_dbg = .false.
            lblock = .false.
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
