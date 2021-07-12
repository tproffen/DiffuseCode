MODULE do_variable_mod
!
CONTAINS
!*****7**************************************************************** 
SUBROUTINE ersetz_variable (line, length, lmask, amask) 
!-                                                                      
!       replaces a substring in an expression by the value of the       
!       appropriate user defined variable.                              
!     This is needed if the parameter is read.                          
!                                                                       
!       Version : 1.0                                                   
!       Date    : 22 Oct 03                                             
!                                                                       
!       Author  : R.B. Neder  (reinhard.neder@fau.de)    
!+                                                                      
USE blanks_mod
USE constants_mod
USE errlist_mod 
USE lib_length
USE precision_mod
USE string_extract_mod
USE variable_mod
!
IMPLICIT none 
!                                                                       
CHARACTER (LEN=*), INTENT(INOUT) :: line 
INTEGER          , INTENT(INOUT) :: length 
LOGICAL  , DIMENSION(   (PREC_STRING          ),0:1), INTENT(INOUT) :: lmask 
INTEGER          , INTENT(INOUT) :: amask     ! active mask line
CHARACTER(LEN=   (PREC_STRING          )) :: string    = ' '
CHARACTER(LEN=   (PREC_STRING          )) :: substring = ' '
CHARACTER(LEN=   (PREC_STRING          )) :: zeile     = ' '
CHARACTER(LEN=   (PREC_STRING          )) :: dummy     = ' '
!                                                                       
INTEGER :: i, ianf, iend, ll
INTEGER :: linsert 
INTEGER :: laenge
INTEGER :: istart, istop
INTEGER :: s1,s2,s3
INTEGER :: omask   = 0         ! old mask 
INTEGER :: nmask   = 1         ! new mask 
LOGICAL :: success = .FALSE.
logical :: lres    = .FALSE.   ! Found a reserved word => true
!                                                                       
!                                                                       
string = line
istart = 1
istop  = length
s1 = 0
s2 = 0
s3 = 0
line = ' '
lmask(:,:)= .TRUE.
omask = 0
nmask = 1
substring = ' '
!
main: DO WHILE(s2<istop)     ! Loop over all non-quoted section of string
   CALL string_extract(string,istart, istop, substring, s1,s2,s3)
!write(*,*) ' STRING >', string(1:LEN_TRIM(string)),'<'
!write(*,*) ' SUB    >', substring(1:len_trim(substring)),'<', len_trim(substring),s1, s2, s3
!write(*,'(1x,a,80L1)') ' MASK   >', lmask    (1:len_trim(substring),omask)
!write(*,*) '         123456789 123456789 1234567890'
!
   linsert = 0
   laenge  = s2 - s1 + 1     ! Length of current string
   ll = laenge 
!                                                                       
!     Loop over all defined variables and search for coincidences       
!     of the variable name. Works since names like "cos" etc are        
!     forbidden, and the names have been sorted in alphabetically        
!     descending order.                                                 
!
!
!     --If an apostrophe is found, ignore the string                    
!
!IF (max (INDEX (string, '''') , INDEX (string, '"') ) .gt.0) THEN 
!   RETURN 
!ENDIF 
   names:DO i = 1, var_num 
!lmask(:,:)= .TRUE.
   success = .FALSE.
   ianf = INDEX_MASK (substring, var_name (i) (1:var_l (i) ), lmask(1:LEN_TRIM(substring),omask)) !, .TRUE. ) 
   DO while (ianf.ne.0) 
!write(*,*) ' zeile Z>', substring(1:len_trim(substring)),'<',var_name (i) (1:var_l (i)), ianf
!write(*,'(1x,a,80L1)') ' MASK Z >', lmask    (1:len_trim(substring),omask)
!write(*,'(1x,a,80L1)') ' MASK Z >', lmask    (1:len_trim(substring),nmask)
!write(*,*) '         123456789 123456789 1234567890'
!write(*,'(1x,a,3i8)') '        ', ianf, var_l(i), var_entry(i)
      zeile = ' ' 
      iend = ianf + var_l (i) - 1 
      call test_reserved(substring, ianf, iend, var_val(VAR_PROGRAM), &    ! Test if inside a reserved string
                         lmask, omask, nmask, lres)
      if_lres:if(lres) then         ! Inside a reserved string, skip this occurence
!write(*,*) ' MASKED RESERVED ??? '
!write(*,*) ' zeile Z>', substring(1:len_trim(substring)),'<',var_name (i) (1:var_l (i)), ianf
!write(*,'(1x,a,80L1)') ' MASK Z >', lmask    (1:len_trim(substring),omask)
!write(*,'(1x,a,80L1)') ' MASK Z >', lmask    (1:len_trim(substring),nmask)
!write(*,*) '         123456789 123456789 123456789 123456789 123456789 '
!write(*,'(1x,a,3i8)') '        ', ianf, var_l(i), var_entry(i)
      else   if_lres
         IF(var_entry(i)>0) CYCLE names        ! This is a variable field
         IF (ianf.gt.1) THEN
            zeile (1:ianf - 1) = substring (1:ianf - 1) 
            lmask (1:ianf-1,nmask) = lmask (1:ianf-1,omask)
!write(*,*) ' zeile A>', zeile    (1:len_trim(zeile    )),'<', ianf
!write(*,'(1x,a,80L1)') ' MASK A >', lmask    (1:len_trim(zeile),omask)
!write(*,'(1x,a,80L1)') ' MASK A >', lmask    (1:len_trim(zeile),nmask)
!write(*,*) '         123456789 123456789 1234567890'
         ENDIF
         IF (var_type (i) .eq.      IS_REAL) THEN
            WRITE (dummy (1:PREC_WIDTH) , PREC_F_REAL) var_val (i) 
            dummy (PREC_MANTIS+1:PREC_MANTIS+1) = 'd' 
            ll = PREC_WIDTH
            CALL rem_bl (dummy, ll) 
            zeile (ianf:ianf + ll - 1) = dummy (1:ll) 
            lmask (ianf:ianf+ll-1,nmask) = .FALSE.
!write(*,*) ' zeile R>', zeile    (1:len_trim(zeile    )),'<'
!write(*,'(1x,a,80L1)') ' MASK R >', lmask    (1:len_trim(zeile),omask)
!write(*,'(1x,a,80L1)') ' MASK R >', lmask    (1:len_trim(zeile),nmask)
!write(*,*) '         123456789 123456789 1234567890'
!write(*,*) ' ianf, ll ', ianf, ll
            linsert = ll 
         ELSEIF (var_type (i) .eq.      IS_INTE) THEN 
            WRITE (dummy (1:PREC_WIDTH) , PREC_F_INTE) nint (var_val (i) ) 
            ll = PREC_WIDTH 
            CALL rem_bl (dummy, ll) 
            zeile (ianf:ianf + ll - 1) = dummy (1:ll) 
            lmask (ianf:ianf+ll-1,nmask) = .FALSE.
!write(*,*) ' zeile I>', zeile    (1:len_trim(zeile    )),'<'
!write(*,'(1x,a,80L1)') ' MASK I >', lmask    (1:len_trim(zeile),omask)
!write(*,'(1x,a,80L1)') ' MASK I >', lmask    (1:len_trim(zeile),nmask)
!write(*,*) '         123456789 123456789 1234567890'
!write(*,*) ' ianf, ll ', ianf, ll
            linsert = ll 
         ELSEIF (var_type (i) .eq.      IS_CHAR) THEN 
            ll = len_str (var_char (i) ) 
            zeile (ianf:ianf) = '''' 
            zeile (ianf + 1:ianf + ll) = var_char (i) (1:ll) 
            zeile (ianf + ll + 1:ianf + ll + 1) = '''' 
            lmask (ianf:ianf+ll+1,nmask) = .FALSE.
!write(*,*) ' zeile Q>', zeile    (1:len_trim(zeile    )),'<'
!write(*,'(1x,a,80L1)') ' MASK Q >', lmask    (1:len_trim(zeile),omask)
!write(*,'(1x,a,80L1)') ' MASK Q >', lmask    (1:len_trim(zeile),nmask)
!write(*,*) '         123456789 123456789 1234567890'
!write(*,*) ' ianf, ll ', ianf, ll
            linsert = ll + 2 
         ENDIF 
         ll = laenge+linsert - (iend-ianf + 1) 
!write(*,*) ' ll neu ', ll, laenge, linsert, ianf, iend, iend.lt.laenge
         IF(iend.lt.laenge) THEN
!write(*,*) 'ADD   ',ianf+linsert, ll, iend+1, laenge, '>',substring(iend+1:laenge),'<'
            zeile(ianf + linsert:ll)       =substring(iend+1:laenge)
            lmask(ianf + linsert:ll,nmask) =lmask(    iend+1:laenge    ,omask)
!        lmask(ianf + linsert:ll,nmask) =lmask(    ianf + linsert:ll,omask)
         ENDIF
         substring = zeile 
!write(*,*) ' zeile X>', zeile    (1:len_trim(zeile    )),'<'
!write(*,*) ' SUB   X>', substring(1:len_trim(substring)),'<'
!write(*,'(1x,a,90L1)') ' MASK  X>', lmask    (1:len_trim(zeile),omask)
!write(*,'(1x,a,900L1)') ' MASK  X>', lmask    (1:len_trim(zeile),nmask)
!write(*,*) '         123456789 123456789 1234567890'
         istart = ianf+linsert                  ! Start further search after insert
      endif if_lres                          ! If block reserved words
      laenge = ll 
      success = .TRUE.
      omask = MOD(omask+1,2)
      nmask = MOD(nmask+1,2)
!     ianf = INDEX      (substring(istart:len_trim(substring)), var_name (i) (1:var_l (i) ) ) 
      ianf = INDEX_MASK (substring, var_name (i) (1:var_l (i) ), lmask(1:LEN_TRIM(substring),omask)) !, .TRUE. ) 
!write(*,*) ' REPEATED ? ', ianf, var_name(i)(1:var_l(i)), len_trim(substring)
!     IF(ianf>0) ianf = ianf + istart - 1    ! if found, correct the offset
   ENDDO 
   ENDDO names
!
   line = line(1:LEN_TRIM(line))//substring(1:LEN_TRIM(substring))  &
          //string(s2+1:s3-1)
   istart = s3
!
ENDDO main
length = LEN_TRIM(line)
amask  = omask
!write(*,*) ' LINE  X>', line(1:len_trim(line)),'<'
!write(*,'(1x,a,90L1)') ' MASK  X>', lmask    (1:len_trim(zeile),omask)
!write(*,*) ' OMASK NMASK ', omask, nmask, amask
!rite(*,*) ' AFTER ERSETZ_VARIABLE >',line(1:len_trim(line)),'<' 
!
END SUBROUTINE ersetz_variable                
!
!*****7**************************************************************** 
!
INTEGER FUNCTION index_mask(string, find, lmask)
!
IMPLICIT NONE
!
CHARACTER(LEN=*), INTENT(IN) :: string   ! The string to be searched
CHARACTER(LEN=*), INTENT(IN) :: find     ! The string to be found
LOGICAL, DIMENSION(LEN_TRIM(string)), INTENT(IN) :: lmask
!
INTEGER :: ls   ! Length of string
INTEGER :: lf   ! Length of find
integer :: i
!
index_mask = 0
ls = LEN_TRIM(string)
lf = LEN_TRIM(find)
!
main: DO i=1,ls-lf + 1
   IF(string(i:i+lf-1)==find .AND. ALL(lmask(i:i+lf-1))) THEN
      index_mask = i
      EXIT main
   ENDIF
ENDDO main
!write(*,*) 'INDMSK ', index_mask, INDEX(string, find)
!write(*,'(a,a)'    ) ' STRINGS ', string(1:len_trim(string))
!write(*,'(a,200L1)') ' lmask   ', lmask (1:len_trim(string))
!write(*,'(a,a)'    ) ' FIND    ', find  (1:len_trim(find ))
!if(index_mask==0 .and. INDEX(string, find)>0) then
!write(*,*) ' LS lF   ', LEN_TRIM(string), len(string), len_trim(find),len(find)
!write(*,'(a,a)'    ) ' search  ', string(INDEX(string, find):INDEX(string, find)+LEN_TRIM(find)-1)
!write(*,'(a,200L1)') ' lmask   ', lmask (INDEX(string, find):INDEX(string, find)+LEN_TRIM(find)-1)
!endif
!
END FUNCTION index_mask
!
!*****7**************************************************************** 
!
subroutine test_reserved(substring, ianf, iend,           prog_name,            &
        lmask, omask, nmask, lres)
!-
!  Test if the user variable is inside a reserved string, if so, the mask is
!  adapted to ignore this section of the substring
!+
!
use reserved_mod
use precision_mod
!
implicit none
!
character(len=*)                                , intent(in)    :: substring   ! String to be tested
integer                                         , intent(in)    :: ianf        ! Start position of variable
integer                                         , intent(in)    :: iend        ! End  position of variable
real(kind=PREC_DP)                              , intent(in)    :: prog_name   ! Suite, discus, diffev etc
logical           , dimension((PREC_STRING),0:1), intent(inout) :: lmask       ! Mask to flag areas to ignore
integer                                         , intent(in)    :: omask       ! Old mask
integer                                         , intent(in)    :: nmask       ! new mask
logical                                         , intent(inout) :: lres        ! true if reserved word was found
!
lres = .FALSE.
call test_reserved_group(substring, ianf, iend,    lib_reserved_n,              &
                         lib_reserved, lmask, omask, nmask, lres)
!
if(nint(prog_name)==0) then         ! SUITE  section
   call test_reserved_group(substring, ianf, iend,  suite_reserved_n,           &
                             suite_reserved, lmask, omask, nmask, lres)
elseif(nint(prog_name)==1) then     ! DISCUS section
   call test_reserved_group(substring, ianf, iend, discus_reserved_n,           &
                            discus_reserved, lmask, omask, nmask, lres)
elseif(nint(prog_name)==2) then     ! DIFFEV section
   call test_reserved_group(substring, ianf, iend, diffev_reserved_n,           &
                            diffev_reserved, lmask, omask, nmask, lres)
elseif(nint(prog_name)==3) then     ! KUPLOT section
   call test_reserved_group(substring, ianf, iend, kuplot_reserved_n,           &
                            kuplot_reserved, lmask, omask, nmask, lres)
elseif(nint(prog_name)==4) then     ! REFINE section
   call test_reserved_group(substring, ianf, iend, refine_reserved_n,           &
                            refine_reserved, lmask, omask, nmask, lres)
endif
!
end subroutine test_reserved
!
!*****7**************************************************************** 
!
subroutine test_reserved_group(substring, ianf, iend,                           &
        reserved_n, reserved, lmask, omask, nmask, lres)
!-
!  Test if the user variable is inside a reserved string, ignore
!  Test the group lib, suite, discus, diffev, kuplot, refine
!+
!
use precision_mod
!
implicit none
!
character(len=*)                              , intent(in)    :: substring   ! String to be tested
integer                                       , intent(in)    :: ianf        ! Start position of variable 
integer                                       , intent(in)    :: iend        ! End position of variable
integer                                       , intent(in)    :: reserved_n  ! Number of reserved words
character(len=*), dimension(1:reserved_n)     , intent(in)    :: reserved    ! The reserved words
logical         , dimension((PREC_STRING),0:1), intent(inout) :: lmask       ! Mask to flag areas to ignore
integer                                       , intent(in)    :: omask       ! Old mask index
integer                                       , intent(in)    :: nmask       ! New mask index
logical                                       , intent(inout) :: lres        ! true if reserved word was found
!
integer :: i      ! Dummy loop index
integer :: ires   ! Location of reserved string
!
loop_res: do i=1,reserved_n
   ires = index(substring, reserved(i)(1:len_trim(reserved(i))) )
   if(ires>0) then
      if(ires<=ianf .and. ires+len_trim(reserved(i))-1>=iend) then
         lres = .TRUE.
         lmask(ires:ires+len_trim(reserved(i))-1,omask) = .FALSE.
         lmask(1:ires+len_trim(reserved(i))-1,nmask) = &
         lmask(1:ires+len_trim(reserved(i))-1,omask)
         exit loop_res
      endif
   endif
enddo loop_res
!
end subroutine test_reserved_group
!
!*****7**************************************************************** 
!
SUBROUTINE upd_variable (string, laenge, wert, dummy, length, substr) 
!-                                                                      
!       updates the user defined variable by the value wert             
!                                                                       
!       Version : 1.0                                                   
!       Date    : 22 Oct 03                                             
!                                                                       
!       Author  : R.B. Neder  (reinhard.neder@mail.uni-wuerzburg.de)    
!+                                                                      
USE constants_mod
USE lib_f90_allocate_mod
USE errlist_mod 
USE param_mod
USE variable_mod
USE precision_mod
!
IMPLICIT none 
!                                                                       
!                                                                       
CHARACTER(LEN=*)     , INTENT(IN) :: string 
INTEGER              , INTENT(IN) :: laenge 
REAL(KIND=PREC_DP)   , INTENT(IN) :: wert 
CHARACTER(LEN=*)     , INTENT(IN) :: dummy 
INTEGER              , INTENT(IN) :: length 
INTEGER, DIMENSION(2), INTENT(IN) :: substr 
!                                                                       
INTEGER :: i 
!                                                                       
      ier_num = - 24 
      ier_typ = ER_FORT 
!                                                                       
!     loop over all variable names. Strict equality is required         
!                                                                       
      DO i = 1, var_num 
      IF (string (1:laenge) .eq.var_name (i) (1:var_l (i) ) ) THEN 
         IF(var_entry(i)>0) RETURN     ! This is a variable array
         IF (var_type (i) .eq.      IS_REAL) THEN 
            var_val (i) = wert 
         ELSEIF (var_type (i) .eq.      IS_INTE) THEN 
            var_val (i) = nint (wert) 
         ELSEIF (var_type (i) .eq.      IS_CHAR) THEN 
            var_char (i)(substr(1):substr(2)) = dummy (1:length) 
         ENDIF 
         IF(string (1:laenge) == 'REF_DIMENSION') THEN
            IF(var_val (i)>UBOUND(ref_para,1)) THEN
               CALL alloc_ref_para(NINT(wert))
            ENDIF
         ENDIF
         ier_num = 0 
         ier_typ = ER_NONE 
      ENDIF 
      ENDDO 
END SUBROUTINE upd_variable                   
!
!*****7**************************************************************** 
!
SUBROUTINE rese_variables(is_diffev)
!-
!     Removes all user variables
!+
USE variable_mod
IMPLICIT none 
!
LOGICAL, INTENT(IN)    :: is_diffev
!
INTEGER :: i, j
INTEGER :: ndel
!
IF(var_num==var_sys) RETURN         ! No variables were defined
!
DO i=var_sys+1, var_num
   IF(var_diff(i) .EQV. is_diffev) THEN
      var_name( i) = ' '
      var_char( i) = ' '
      var_l   ( i) = 1
      var_type( i) = 0
      var_diff( i) = .FALSE.
      var_val ( i) = 0
      IF(var_entry(i)/=0) THEN      ! GOT An array
         IF(ALLOCATED(var_field(var_entry(i))%var_value)) &
          DEALLOCATE(var_field(var_entry(i))%var_value)
         IF(ALLOCATED(var_field(var_entry(i))%var_char)) &
          DEALLOCATE(var_field(var_entry(i))%var_char)
         var_field(var_entry(i))%var_shape(:) = 0
         var_entry(i) = 0
      ENDIF
   ENDIF
ENDDO
ndel = 0
i = var_sys+1
remove: DO 
   IF(var_name(i) == ' ') THEN
      ndel = ndel +1
      IF(i>var_num-ndel) EXIT remove
      DO j=i, var_num-ndel
         var_name (j) = var_name (j+1)
         var_char (j) = var_char (j+1)
         var_l    (j) = var_l    (j+1)
         var_type (j) = var_type (j+1)
         var_diff (j) = var_diff (j+1)
         var_val  (j) = var_val  (j+1)
         var_entry(j) = var_entry(j+1)
      ENDDO
   ELSE
      i=i+1
   ENDIF
ENDDO remove
var_num = var_num - ndel
!
END SUBROUTINE rese_variables 
!
!*****7**************************************************************** 
 SUBROUTINE del_variables (MAXW, ianz, cpara, lpara, is_diffev)
!-
!     Removes all or specified user variables
!+
USE errlist_mod
USE str_comp_mod
USE variable_mod
!
IMPLICIT none 
!
INTEGER                               , INTENT(IN)    :: MAXW
INTEGER                               , INTENT(IN)    :: ianz
CHARACTER(LEN=*   ), DIMENSION(1:MAXW), INTENT(INOUT) :: cpara
INTEGER            , DIMENSION(1:MAXW), INTENT(INOUT) :: lpara
LOGICAL                               , INTENT(IN)    :: is_diffev
!
INTEGER :: i, j, k, upper
!
!
IF(ianz==1) THEN
   CALL rese_variables (is_diffev)
ELSEIF(str_comp (cpara (2) , 'all', 3, lpara(1), 3) ) THEN 
   CALL rese_variables (is_diffev)
ELSE
   main: DO i=2, ianz
      upper = var_num
      loop: DO j=var_sys+1, upper
         IF(cpara(i) == var_name(j)) THEN 
            IF(var_diff(j) .EQV. is_diffev) THEN
               IF(var_entry(j)/=0) THEN      ! GOT An array
                  IF(ALLOCATED(var_field(var_entry(j))%var_value)) &
                    DEALLOCATE(var_field(var_entry(j))%var_value)
                  IF(ALLOCATED(var_field(var_entry(j))%var_char )) &
                    DEALLOCATE(var_field(var_entry(j))%var_char )
                  var_field(var_entry(j))%var_shape(:) = 0
                  var_entry(j) = 0
               ENDIF
               DO k=j+1, var_num
                  var_name (k-1) = var_name(k)
                  var_char (k-1) = var_char(k)
                  var_l    (k-1) = var_l(k)
                  var_type (k-1) = var_type(k)
                  var_diff (k-1) = var_diff(k)
                  var_val  (k-1) = var_val(k)
                  var_entry(k-1) = var_entry(k)
               ENDDO
               var_num = var_num - 1
               CYCLE main
            ELSE
               ier_num =  -39
               ier_typ = ER_FORT
               ier_msg(1) = 'Offending variable: '//cpara(i)(1:LEN_TRIM(cpara(i)))
               RETURN
            ENDIF
         ENDIF
      ENDDO loop
      ier_num = -24
      ier_typ = ER_FORT
      ier_msg(1) = 'Offending variable: '//cpara(i)(1:LEN_TRIM(cpara(i)))
      RETURN
   ENDDO main
ENDIF
!
END SUBROUTINE del_variables 
!
!*****7**************************************************************** 
!
SUBROUTINE show_variables 
!-                                                                      
!     shows the variable definitions                                    
!+                                                                      
USE constants_mod
USE blanks_mod
USE lib_length
USE precision_mod
USE prompt_mod 
USE variable_mod
IMPLICIT none 
!
!
CHARACTER(LEN=6), DIMENSION(0:1) :: cdiffev
CHARACTER(LEN=PREC_STRING) :: string
INTEGER :: i , j
!
INTEGER :: length
DATA cdiffev/'DIFFEV','      '/
!
IF (var_num.gt.0) THEN 
   WRITE (output_io, 2000) var_num, VAR_MAX 
   DO i = 1, var_num 
      j=1
      IF(var_diff(i)) j=0
      string = ' '
      length = 0
      IF(var_entry(i)>0) THEN
         IF(var_field(var_entry(i))%var_shape(2)>1) THEN
            WRITE(string,'(a1,i8,a1,i8,a1)') '[',var_field(var_entry(i))%var_shape(1), &
            ',',var_field(var_entry(i))%var_shape(2),']'
            length = 19
            CALL rem_bl(string, length)
         ELSE
            WRITE(string,'(a1,i8,a1)') '[',var_field(var_entry(i))%var_shape(1),']'
            length = 10
            CALL rem_bl(string, length)
         ENDIF
      ENDIF
      IF (var_type (i) .eq.      IS_REAL) THEN 
         IF (ABS(var_val(i)) .lt.1e-5.or. ABS(var_val(i)).ge.1e6) THEN
            WRITE(output_io, 2100) var_name (i), var_val (i) , cdiffev(j), string(1:length)
         ELSE 
            WRITE(output_io, 2150) var_name (i), var_val (i) , cdiffev(j), string(1:length)
         ENDIF 
      ELSEIF (var_type (i) .eq.      IS_INTE) THEN 
         WRITE(output_io, 2200) var_name (i), nint (var_val (i) ) , cdiffev(j), string(1:length)
      ELSEIF (var_type (i) .eq.      IS_CHAR) THEN 
         WRITE(output_io, 2300) var_name(i), var_char(i)(1:len_str(var_char(i))), cdiffev(j)
      ENDIF 
   ENDDO 
ELSE 
   WRITE (output_io, 2050) VAR_MAX 
ENDIF 
WRITE (output_io, * ) 
 2000 FORMAT    (' User defined variables',i4,                          &
     &                  ' Maximum number: ',i3/)                        
 2050 FORMAT    (' No user defined variables',                          &
     &                  ' Maximum number: ',i3/)                        
 2100 FORMAT    (' Name: ',a30,' = ',8x,e20.8e2,' Real      ',a,a) 
 2150 FORMAT    (' Name: ',a30,' = ',8x,f16.8  ,'     Real      ',a,a) 
 2200 FORMAT    (' Name: ',a30,' = ',i15,13x,   ' Integer   ',a,a) 
 2300 FORMAT    (' Name: ',a30,' = ''',a,     ''' Character ',a) 
END SUBROUTINE show_variables                 
!
!*****7**************************************************************** 
!
SUBROUTINE validate_variable (zeile, lp) 
!-                                                                      
!       Checks whether the variable name is legal                       
!     All intrinsic function names, variable names and parts of         
!     these names are illegal                                           
!                                                                       
!       Version : 1.0                                                   
!       Date    : 22 Oct 03                                             
!                                                                       
!       Author  : R.B. Neder  (reinhard.neder@mail.uni-wuerzburg.de)    
!+                                                                      
USE reserved_mod
USE charact_mod
USE errlist_mod 
USE set_sub_generic_mod
USE prompt_mod
!
IMPLICIT none 
!                                                                       
CHARACTER (LEN=*), INTENT(IN) :: zeile 
INTEGER,           INTENT(IN) :: lp 
!                                                                       
INTEGER :: i, ii 
INTEGER :: length
LOGICAL :: lok 
!                                                                       
!                                                                       
ier_num = 0 
ier_typ = ER_NONE 
!                                                                       
main:DO i = 1, lib_reserved_n 
!        IF (INDEX (lib_reserved (i), zeile (1:lp) ) .NE.0) THEN 
   length = MAX(LEN_TRIM(lib_reserved(i)), LEN_TRIM(zeile(1:lp)))
   length = MIN(length, LEN(lib_reserved), LEN(zeile))
   IF ((lib_reserved (i)(1:length)== zeile (1:length) ) ) THEN 
      ier_num = - 25 
      ier_typ = ER_FORT 
      RETURN
   ENDIF 
ENDDO main 
!                                                                       
!     Check against variables/functions of the main program             
!                                                                       
IF(lstandalone) THEN
   IF (ier_num.eq.0) THEN 
      CALL p_validate_var_spec (zeile, lp) 
      IF(ier_num /= 0) RETURN
   ENDIF 
ELSE        ! Part of the suite, check all parts 
   CALL general_validate_var_spec (zeile, lp, diffev_reserved_n, diffev_reserved)
   IF(ier_num /= 0) RETURN
   CALL general_validate_var_spec (zeile, lp, discus_reserved_n, discus_reserved)
   IF(ier_num /= 0) RETURN
   CALL general_validate_var_spec (zeile, lp, kuplot_reserved_n, kuplot_reserved)
   IF(ier_num /= 0) RETURN
   CALL general_validate_var_spec (zeile, lp, refine_reserved_n, refine_reserved)
   IF(ier_num /= 0) RETURN
   CALL general_validate_var_spec (zeile, lp,  suite_reserved_n,  suite_reserved)
   IF(ier_num /= 0) RETURN
ENDIF 
!                                                                       
!     Check that the name contains only letters, Numbers and the "_"    
!                                                                       
IF (ier_num.eq.0) THEN 
   lok = .true. 
   DO i = 1, lp 
      ii = iachar (zeile (i:i) ) 
      lok = lok.and. (zero.le.ii.and.ii.le.nine.or.aa.le.ii.and. & 
                      ii.le.zz .or.u.eq.ii.or.a.le.ii.and.ii.le.z)                               
   ENDDO 
   IF (.not.lok) THEN 
      ier_num = - 26 
      ier_typ = ER_FORT 
   ENDIF 
ENDIF 
!                                                                       
END SUBROUTINE validate_variable              
!
!*******************************************************************************
!
SUBROUTINE general_validate_var_spec (zeile, lp, general_reserved_n, general_reserved)
!-                                                                      
!  checks whether the variable name is legal, GENERAL specific part 
!                                                                       
!  Version : 1.0                                                   
!  Date    : 20 Jun 19                                             
!                                                                       
!  Author  : R.B. Neder  (reinhard.neder@fau.de)    
!+                                                                      
USE errlist_mod
!                                                                       
IMPLICIT none
!                                                                       
CHARACTER (LEN=*), INTENT(IN) :: zeile
INTEGER,           INTENT(IN) :: lp
INTEGER,           INTENT(IN) :: general_reserved_n
CHARACTER(LEN=16), DIMENSION(1:general_reserved_n), INTENT(IN) :: general_reserved
!
INTEGER  :: i , length
!                                                                       
ier_num = 0
ier_typ = ER_NONE
!                                                                       
main: DO i = 1, general_reserved_n
!
   length = MAX(LEN_TRIM(general_reserved(i)), LEN_TRIM(zeile(1:lp)))
   length = MIN(length, LEN(general_reserved), LEN(zeile))
   IF(general_reserved (i)(1:length)== zeile(1:length) ) THEN
      ier_num = -25
      ier_typ = ER_FORT
      EXIT main
   ENDIF
ENDDO main
!                                                                       
END SUBROUTINE general_validate_var_spec
!
!*******************************************************************************
!
END MODULE do_variable_mod
