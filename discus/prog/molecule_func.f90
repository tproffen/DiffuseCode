module molecule_func_mod
!
USE conn_mod
USE conn_sup_mod
USE molecule_mod
!
IMPLICIT NONE
PRIVATE
PUBLIC  :: do_molecularize
PUBLIC  :: molecularize_sub
PUBLIC  :: molecularize_numbers
!
!
LOGICAL, DIMENSION(:  ), ALLOCATABLE :: t_list    ! Temporary list of all atoms in the molecule
!
CONTAINS
!
!*******************************************************************************
!
SUBROUTINE do_molecularize (line, length)
!-
!  Groups atoms into a molecule
!  Two styles are available:
!  molecularize conn, <no>, <excl1> [, <excl2>...] [,type:<mol_type>],
!               [, start:<number>]
!               [, finish:<number>]
!               [, exclude:<number>]  [, exclude:[ex1, ex2,...]]
!               [, conn:{<name>| <number}]
!               [, type:<mol_type>]
!               [, biso:<biso>]
!               [, corrlin:<clin>]
!               [, corrquad:<cquad>]
!
!  molecularize range, <start> , <finish>, <mol_type>, <biso>
!
USE crystal_mod
use get_iscat_mod
!
USE ber_params_mod
USE errlist_mod
USE get_params_mod
USE precision_mod
USE str_comp_mod
use take_param_mod
!
IMPLICIT NONE
!
CHARACTER(LEN=*), INTENT(IN) :: line
INTEGER         , INTENT(INOUT) :: length
!
INTEGER            , PARAMETER       :: MAXW  = 20
INTEGER            , PARAMETER       :: MAXWW = 20
logical            , parameter       :: lnew = .FALSE.
!
CHARACTER(LEN=MAX(PREC_STRING,LEN(line))), DIMENSION(MAXW) :: cpara ! 
CHARACTER(LEN=MAX(PREC_STRING,LEN(line))), DIMENSION(MAXW) :: ccpara ! 
INTEGER            , DIMENSION(MAXW) :: lpara ! 
INTEGER            , DIMENSION(MAXW) :: llpara ! 
REAL(KIND=PREC_DP) , DIMENSION(MAXW) :: werte ! 
REAL(KIND=PREC_DP) , DIMENSION(MAXW) :: wwerte ! 
INTEGER                              :: ianz
INTEGER                              :: iianz
INTEGER                              :: istart, ifinish
INTEGER                              :: natom
INTEGER, DIMENSION(:), ALLOCATABLE   :: iatoms
INTEGER                              :: itype
INTEGER                              :: i
REAL(kind=PREC_DP)                   :: biso
REAL(kind=PREC_DP)                   :: clin
REAL(kind=PREC_DP)                   :: cqua
character(len=PREC_STRING), dimension(:), allocatable :: con_name
integer                   , dimension(:), allocatable :: con_ino
integer :: con_nnums    ! Number of connectivities stated with numbers
integer :: con_nnames   ! Number of connectivities stated with names
!
logical                   , dimension(:), allocatable :: cen_ino
integer :: cen_nnums    ! Number of central atom types stated with numbers
integer :: cen_nnames   ! Number of central atom types stated with names
!
integer, parameter :: NOPTIONAL = 9
integer, parameter :: O_CLIN    = 1
integer, parameter :: O_CQUAD   = 2
integer, parameter :: O_BISO    = 3
integer, parameter :: O_TYPE    = 4
integer, parameter :: O_FINISH  = 5
integer, parameter :: O_START   = 6
integer, parameter :: O_EXCLUDE = 7
integer, parameter :: O_CONN    = 8
integer, parameter :: O_CENTRAL = 9
character(LEN=   8), dimension(NOPTIONAL) :: oname   !Optional parameter names
character(LEN=PREC_STRING), dimension(NOPTIONAL) :: opara   !Optional parameter strings returned
integer            , dimension(NOPTIONAL) :: loname  !Lenght opt. para name
integer            , dimension(NOPTIONAL) :: lopara  !Lenght opt. para name returned
logical            , dimension(NOPTIONAL) :: lpresent!opt. para is present
real(kind=PREC_DP) , dimension(NOPTIONAL) :: owerte   ! Calculated values
integer, parameter                        :: ncalc = 4 ! Number of values to calculate 
!
data oname  / 'corrlin', 'corrquad', 'biso'   , 'moletype',                     &
              'finish' , 'start'   , 'exclude', 'conn'    , 'central'  /
data loname /  7,         8,          4       ,  8        ,  6       ,          &
               5       ,  7        ,  4       ,  7         /
opara  =  (/ '0.0000'  , '0.0000'  , '0.0000' , '1.0000'  , '1.0000' ,          &
             '-1.000'  , '0.0000'  , '1.0000' , 'all   '             /)   ! Always provide fresh default values
lopara =  (/  6,          6,          6       ,  6        ,  6       ,          &
              6,          6        ,  6       ,  6                   /)
owerte =  (/  0.0,        0.0,        0.0     ,  1.0      , 0.0,                &
              0.0,        0.0      ,  0.0     , -1.0                 /)
!
!
werte = 0.0D0
!
CALL get_params (line, ianz, cpara, lpara, MAXW, length)
!      CALL ber_params (ianz, cpara, lpara, werte, maxw)
if(ier_num /=0 ) return
!
call get_optional(ianz, MAXW, cpara, lpara, NOPTIONAL,  ncalc, &
                  oname, loname, opara, lopara, lpresent, owerte)
if(ier_num /=0 ) return
!
clin  = owerte(O_CLIN)
cqua  = owerte(O_CQUAD)
biso  = owerte(O_BISO)
itype = owerte(O_TYPE)
!
! Evaluate optional START
!
istart  = 1
ifinish = cr_natoms        ! Default is full range
!
if(lpresent(O_START)) then
   if(opara(O_START)=='all') then
      istart  = 1
      ifinish = cr_natoms
   else
      ccpara(1) = opara(O_START)
      llpara(1) = lopara(O_START)
      iianz = 1
      call ber_params(iianz, ccpara, llpara, wwerte, MAXWW)
      if(ier_num /=0 ) return
      istart = NINT(wwerte(1))
   endif
else
   istart = 1
endif
!
! Evaluate optional FINISH
!
if(lpresent(O_FINISH)) then
   if(opara(O_FINISH)=='all') then
      istart  = 1
      ifinish = cr_natoms
   elseif(opara(O_FINISH)=='last') then
      ifinish = cr_natoms
   else
      ccpara(1) = opara(O_FINISH)
      llpara(1) = lopara(O_FINISH)
      iianz = 1
      call ber_params(iianz, ccpara, llpara, wwerte, MAXWW)
      if(ier_num /=0 ) return
      ifinish = NINT(wwerte(1))
   endif
endif
!
! Evaluate optional EXCLUDE
!
if(lpresent(O_EXCLUDE)) then 
   if(istart==ifinish) then    ! A single atom is checked
      call get_optional_multi(MAXWW, opara(O_EXCLUDE), lopara(O_EXCLUDE), wwerte, iianz)
      natom = iianz + 1
      allocate(iatoms(1:natom))
      do i=1, iianz
         iatoms(1+i) = nint(wwerte(i))
      enddo
   else
      ier_num = -181
      ier_typ = ER_APPL
      return
   endif
endif
!
! Evaluate optional CONNectivity
!
con_nnums = 0
con_nnames = 0
if(lpresent(O_CONN)) then
   call eva_optional_multi(opara(O_CONN), lopara(O_CONN), MAXW, ccpara, wwerte, con_nnums, con_nnames)
   if(allocated(con_ino)) deallocate(con_ino)
   if(allocated(con_name)) deallocate(con_name)
   allocate(con_ino(1:con_nnums))
   allocate(con_name(1:con_nnames))
   do i=1, con_nnums
      con_ino(i) = nint(wwerte(i))
   enddo
   do i=1, con_nnames
      con_name(i) = ccpara(i)
   enddo
else
   con_nnums = 1
   con_nnames = 0
   if(allocated(con_ino)) deallocate(con_ino)
   if(allocated(con_name)) deallocate(con_name)
   allocate(con_ino(1:1))
   allocate(con_name(1:1))
   con_ino(1) = 1
   con_name(1) = ' '
endif
!
! Evaluate optional CENTRAL atom type
!
if(lpresent(O_CENTRAL)) then
   call eva_optional_multi(opara(O_CENTRAL), lopara(O_CENTRAL), MAXWW, ccpara, wwerte, cen_nnums, cen_nnames)
   if(allocated(cen_ino)) deallocate(cen_ino)
   allocate(cen_ino(0:MAXSCAT) )
   cen_ino = .FALSE.
   do i=1, cen_nnums
      cen_ino(nint(wwerte(i))) = .TRUE.
   enddo
   do i=1, cen_nnames
      llpara(i) = len_trim(ccpara(i))
   enddo
   if(cen_nnames>0) then
      call get_iscat(cen_nnames, ccpara, llpara, wwerte, MAXWW, lnew)
      if(ier_num/=0) then
         if(allocated(cen_ino)) deallocate(cen_ino)
         ier_msg(1) = 'Central atom types are wrong'
         ier_msg(2) = 'Wrong atoms: ' // opara(O_CENTRAL)(1:lopara(O_CENTRAL))
         return
      endif
      if(nint(wwerte(1))==-1) then    ! all atoms allowed
         cen_ino = .TRUE.
      else
         do i=1, cen_nnames
            cen_nnums = cen_nnums + 1
            cen_ino(nint(wwerte(i))) = .TRUE.
         enddo
      endif
   endif
else
   cen_nnums = 1
   cen_nnames = 0
   if(allocated(cen_ino)) deallocate(cen_ino)
   allocate(cen_ino(0:MAXSCAT))
   cen_ino = .TRUE.
endif
!
!     Start actual work
!
IF(IANZ == 1 ) THEN
   IF(str_comp(cpara(1), 'connectivity', 1, lpara(1), 12)) THEN
      if(allocated(iatoms)) then        ! singel atom with exclusion
         iatoms(1) = istart
         CALL molecularize(natom,iatoms, cen_ino, con_nnums, con_ino, &
                           con_nnames, con_name, itype, biso, clin, cqua)
      else
         natom = 1
         ALLOCATE(iatoms(1:natom))
         do i=istart,ifinish
            iatoms(1) = i
            CALL molecularize(natom,iatoms, cen_ino, con_nnums, con_ino, &
                              con_nnames, con_name, itype, biso, clin, cqua)
         ENDDO
      endif
      DEALLOCATE(iatoms)
   ELSEIF(str_comp(cpara (1) , 'range', 1, lpara (1) , 5)) THEN
      CALL molecularize_numbers(istart, ifinish, itype, biso, clin, cqua)
   endif
ELSEIF(IANZ > 0 ) THEN
   IF(str_comp(cpara (1) , 'connectivity', 1, lpara (1) , 12)) THEN
      CALL del_params (1, ianz, cpara, lpara, maxw)
      CALL ber_params (ianz, cpara, lpara, werte, maxw)
      IF(ier_num == 0) THEN
         IF(NINT(werte(1)) == -1 .AND. ianz==1) THEN
            istart  = 1
            ifinish = cr_natoms
         ELSEIF(NINT(werte(1)) > 0) THEN
            istart  = NINT(werte(1))
            ifinish = istart
         ELSE
            ier_num = -6
            ier_typ = ER_FORT
         ENDIF
         IF(ier_num == 0) THEN
            natom = ianz
            ALLOCATE(iatoms(1:natom))
            iatoms(1:natom) = NINT(werte(1:ianz))
            DO i=istart,ifinish
               iatoms(1) = i
               CALL molecularize(natom,iatoms, cen_ino, con_nnums,   &
                                 con_ino, con_nnames, con_name, itype, biso, clin, cqua)
            ENDDO
            DEALLOCATE(iatoms)
         ENDIF
      ENDIF
   ELSEIF(str_comp(cpara (1) , 'range', 1, lpara (1) , 5)) THEN
      CALL del_params (1, ianz, cpara, lpara, maxw)
      CALL ber_params (ianz, cpara, lpara, werte, maxw)
      IF(ier_num == 0) THEN
         IF(ianz == 4) THEN
         istart  = NINT(werte(1))
         ifinish = NINT(werte(2))
         itype   = NINT(werte(3))
         biso    =      werte(4)
         clin    =      werte(5)
         cqua    =      werte(6)
         CALL molecularize_numbers(istart, ifinish, itype, biso, clin, cqua)
         ELSE
            ier_num = -6
            ier_typ = ER_FORT
            ier_msg(1) = 'Molecularize command needs parameters:'
            ier_msg(2) = ' range, <from> , <to>, <type>, <biso>'
         ENDIF
      ENDIF
   ELSE
      ier_num = -6
      ier_typ = ER_FORT
      ier_msg(1) = 'Molecularize command needs parameters:'
      ier_msg(2) = 'conn, <central> [, <exclude>...]'
      ier_msg(3) = 'range, <from> , <to>, <type>, <biso>'
   ENDIF
ELSE
   ier_num = -6
   ier_typ = ER_FORT
   ier_msg(1) = 'Molecularize command needs parameters:'
   ier_msg(2) = 'conn, <central> [, <exclude>...]'
   ier_msg(3) = 'range, <from> , <to>, <type>, <biso>'
ENDIF
!
if(allocated(con_ino)) deallocate(con_ino)
if(allocated(con_name)) deallocate(con_name)
!
END SUBROUTINE do_molecularize
!
!
!*******************************************************************************
!
SUBROUTINE molecularize(natom,iatoms, cen_ino, con_nnums, con_ino,   &
                        con_nnames, con_name, itype, biso, clin, cqua)
!
! Groups atoms starting with the first into a new molecule with new molecule type
! Atoms no. 2 etc and their connectivity are excluded.
!
USE discus_allocate_appl_mod
USE crystal_mod
USE prop_para_mod
USE errlist_mod
IMPLICIT NONE
!
INTEGER,                                   INTENT(IN)  :: natom
INTEGER,          DIMENSION(1:natom),      INTENT(IN)  :: iatoms
logical,          DIMENSION(0:MAXSCAT  ),  INTENT(IN)  :: cen_ino
INTEGER,                                   INTENT(IN)  :: con_nnums
INTEGER,          DIMENSION(1:con_nnums),  INTENT(IN)  :: con_ino
INTEGER,                                   INTENT(IN)  :: con_nnames
character(len=*), DIMENSION(1:con_nnames), INTENT(IN)  :: con_name
integer                                  , intent(in)  :: itype
REAL(kind=PREC_DP)                                     , intent(in)  :: biso
REAL(kind=PREC_DP)                       , intent(in)  :: clin
REAL(kind=PREC_DP)                       , intent(in)  :: cqua
!
INTEGER   :: jatom
!
INTEGER   :: n_new,n_atom, n_type, n_mole
INTEGER   :: i, j
INTEGER   :: imole
logical   :: liscat
!
liscat = .FALSE.
loop_iscat: do i=1, MAXSCAT
   if(cen_ino(cr_iscat(1,iatoms(1)))) then
      liscat = .TRUE.
      exit loop_iscat
   endif
enddo loop_iscat
!
if(.NOT. liscat) return
!
jatom   = iatoms(1)                        ! Set central atom number
IF(.NOT. btest(cr_prop(jatom),1)) THEN      ! Atom is not yet inside a molecule
!
!  ktype = cr_iscat(1,iatoms(1))
!  ino_max = get_connectivity_numbers(ktype)  ! Get number of connectivity definitions
!
   ALLOCATE(t_list(1:cr_natoms))              ! Make temporary list of atoms already in the molecule
   t_list(:)     = .FALSE.
   t_list(jatom) = .TRUE.
!
!  Accumulate all neighbors of the current atom into temp_list
!
   IF(natom==1) THEN
      CALL mol_add_conn(jatom, con_nnums, con_ino)
   ELSE
      CALL mol_add_excl(jatom, natom, iatoms)
   ENDIF
   n_new = 0
   DO i=1, UBOUND(t_list,1)
      IF(t_list(i)) THEN
         n_new = n_new + 1
      ENDIF
   ENDDO
!
   IF(n_new == 1) THEN
      DEALLOCATE(t_list)
      RETURN
   ENDIF
!
   n_mole = mole_num_mole + 1   ! Reserve new molecule
   n_type = max(itype, mole_num_type)      ! If necessary Reserve new molecule type
   n_atom = mole_off(mole_num_mole) + mole_len(mole_num_mole) + n_new

   IF ( n_mole > MOLE_MAX_MOLE  .or.  &
        n_type > MOLE_MAX_TYPE  .or.  &
        n_atom > MOLE_MAX_ATOM      ) THEN
      n_mole = MAX(n_mole,MOLE_MAX_MOLE)
      n_type = MAX(n_type,MOLE_MAX_TYPE)
      n_atom = MAX(n_atom,MOLE_MAX_ATOM)
      CALL alloc_molecule(1, 1, n_mole, n_type, n_atom)
      IF ( ier_num /= 0) THEN
         RETURN
      ENDIF
   ENDIF
   mole_num_mole = mole_num_mole + 1
   if(itype >  mole_num_type) then             ! Molecule type does not exist
      mole_num_type = mole_num_type + 1
   endif
   mole_len(mole_num_mole)  = n_new
   mole_off(mole_num_mole)  = mole_off(mole_num_mole-1)+mole_len(mole_num_mole-1)
   mole_type(mole_num_mole) = itype            ! User defined type
   mole_char(mole_num_mole) = MOLE_ATOM
   mole_file(mole_num_mole) = ' '
   mole_biso(itype)         = biso
   mole_clin(itype)         = clin
   mole_cqua(itype)         = cqua
   mole_num_atom = mole_off(mole_num_mole) + mole_len(mole_num_mole)
!
   j = 1
   mole_cont(mole_off(mole_num_mole)+j) = jatom        ! Put central atom into first molecule entry
   cr_prop(jatom) = ibset(cr_prop(jatom),PROP_MOLECULE)
   cr_mole(jatom) = mole_num_mole
   t_list(jatom) = .false.                             ! Clear entry as atom is already in molecule
   DO i=1,UBOUND(t_list,1)
      IF(t_list(i)) THEN
         j = j + 1
         mole_cont(mole_off(mole_num_mole)+j) = i
         cr_prop(i) = ibset(cr_prop(i),PROP_MOLECULE)
         cr_mole(i) = mole_num_mole
      ENDIF
   ENDDO
!
   DEALLOCATE(t_list)
!
ELSE                                          ! Atom is already inside a molecule
!
   ALLOCATE(t_list(1:cr_natoms))              ! Make temporary list of atoms already in the molecule
   t_list(:)     = .FALSE.
   t_list(jatom) = .TRUE.
!
!  Accumulate all neighbors of the current atom into temp_list
!
   IF(natom==1) THEN
      CALL mol_add_conn(jatom, con_nnums, con_ino)
   ELSE
      CALL mol_add_excl(jatom, natom, iatoms)
   ENDIF
   imole = cr_mole(jatom)                    ! Atom is in this molecule
   n_new = 0
   DO i=1, UBOUND(t_list,1)
      IF(t_list(i) .AND. cr_mole(i)/=imole) THEN
         n_new = n_new + 1
      ENDIF
   ENDDO
!
   CALL molecule_shift(n_new, imole+1)
   j = mole_len(imole)
   DO i=1,UBOUND(t_list,1)
      IF(t_list(i).AND.cr_mole(i)/=imole) THEN
         j = j + 1
         mole_cont(mole_off(imole)+j) = i
         cr_prop(i) = ibset(cr_prop(i),PROP_MOLECULE)
         cr_mole(i) = imole
      ENDIF
   ENDDO
   mole_len(imole) = mole_len(imole) + n_new
!
   DEALLOCATE(t_list)
!
ENDIF
!
END SUBROUTINE molecularize
!
!*******************************************************************************
!
SUBROUTINE molecularize_sub(jatom,natom,iatoms,max_sub,n_sub, sub_list)
!
! Groups atoms starting with the first into a new molecule with new molecule type
! Atoms no. 2 etc and their connectivity are excluded.
!
USE discus_allocate_appl_mod
USE crystal_mod
USE errlist_mod
IMPLICIT NONE
!
INTEGER,                                  INTENT(IN)  :: jatom
INTEGER,                                  INTENT(IN)  :: natom
INTEGER, DIMENSION(1:natom),              INTENT(IN)  :: iatoms
INTEGER,                                  INTENT(IN)  :: max_sub
INTEGER,                                  INTENT(OUT) :: n_sub
INTEGER, DIMENSION(1:max_sub),            INTENT(OUT) :: sub_list
!
INTEGER   :: i
!
integer,               parameter :: loc_nnums = 1       ! All further recursions
integer, dimension(1), parameter :: loc_ino   = (/-1/)  ! use all connectivities
!
!ALLOCATE(t_list(1:max_sub))              ! Make temporary list of atoms already in the molecule
ALLOCATE(t_list(1:cr_natoms))              ! Make temporary list of atoms already in the molecule
t_list(:)     = .FALSE.
t_list(jatom) = .TRUE.
!
!  Accumulate all neighbors of the current atom into temp_list
!
IF(natom==0) THEN
   CALL mol_add_conn(jatom, loc_nnums, loc_ino)
ELSE
   CALL mol_add_excl(jatom, natom, iatoms)
ENDIF
n_sub = 0
DO i=1, UBOUND(t_list,1)
   IF(t_list(i)) THEN
      n_sub = n_sub + 1
      sub_list(n_sub) = i
   ENDIF
ENDDO
!write(*,*) ' Centr ', jatom
!write(*,*) ' Excl  ', natom,' : ', iatoms(:)
!write(*,*) ' GROUP ', n_sub
!write(*,*) ' GROUP ', sub_list(:)
!
!write(*,*) ' STATUS ', t_list
DEALLOCATE(t_list)
!
END SUBROUTINE molecularize_sub
!
!*******************************************************************************
!
RECURSIVE SUBROUTINE mol_add_conn(katom, con_nnums, con_ino)
!
USE crystal_mod
IMPLICIT NONE
!
INTEGER,                          INTENT(IN) :: katom
INTEGER,                          INTENT(IN) :: con_nnums
INTEGER, DIMENSION(1:con_nnums),  INTENT(IN) :: con_ino
!
INTEGER   :: ktype
INTEGER   :: ino_max  ! Number of connectivities around central atom
INTEGER   :: ino, inoo
INTEGER                              :: c_natoms  ! Number of atoms connected
INTEGER, DIMENSION(:  ), ALLOCATABLE :: c_list    ! List of atoms connected to current
INTEGER, DIMENSION(:,:), ALLOCATABLE :: c_offs    ! Offsets for atoms connected to current
!
INTEGER :: latom
INTEGER :: i
!
integer,               parameter :: loc_nnums = 1       ! All further recursions
integer, dimension(1), parameter :: loc_ino   = (/-1/)  ! use all connectivities
!
!  Accumulate all neighbors of the current atom into temp_neig
!
ktype = cr_iscat(1,katom)
ino_max = get_connectivity_numbers(ktype)  ! Get number of connectivity definitions
loop_ino: DO ino = 1, ino_max                        ! Loop over all connectivity types
   inoo = -1
   if(con_ino(1)==-1) then                 ! All connectivities allowed
      inoo = ino
   else                                    ! Restricted connectivities, search
      loop_search: do i=1,con_nnums
         if(ino==con_ino(i)) then          ! Found matching connectivity
            inoo = ino
            exit loop_search
         endif
      enddo loop_search
   endif
   if(inoo==-1) cycle loop_ino             ! This connectivity does not match
!
   CALL get_connectivity_list(katom, ktype, inoo, c_list, c_offs, c_natoms)
   IF(c_natoms>0) THEN
   DO i=1,c_natoms                         ! Loop over all neighbors
      IF(.NOT.t_list(c_list(i))) THEN      ! Not yet in the molecule
         t_list(c_list(i)) = .TRUE.        ! Add to temporary list
         latom = c_list(i)                 ! Now search current neighbor for further ones
!        Apply periodic boundary shifts to keep molecule contiguous
         IF(c_offs(1,i)/=0 .OR. c_offs(2,i)/=0 .OR. c_offs(3,i)/=0 ) THEN
            cr_pos(:,latom) = cr_pos(:,latom) + c_offs(:,i)
            CALL conn_update(latom, REAL(c_offs(:,i),kind=PREC_DP))
         ENDIF
         CALL mol_add_conn(latom, loc_nnums, loc_ino)
      ENDIF
   ENDDO
   DEALLOCATE(c_list)                      ! Remove temporary list of atoms
   DEALLOCATE(c_offs)                      ! Remove temporary offsets
   ENDIF
ENDDO loop_ino
!
END SUBROUTINE mol_add_conn
!
!*******************************************************************************
!
RECURSIVE SUBROUTINE mol_add_excl(katom, n_excl, excl)
!
USE crystal_mod
IMPLICIT NONE
!
INTEGER , INTENT(IN) :: katom
INTEGER , INTENT(IN) :: n_excl
INTEGER , DIMENSION(1:n_excl), INTENT(IN) :: excl
!
INTEGER   :: ktype
INTEGER   :: ino_max  ! Number of connectivities around central atom
INTEGER   :: ino, inoo
INTEGER                              :: c_natoms  ! Number of atoms connected
INTEGER, DIMENSION(:  ), ALLOCATABLE :: c_list    ! List of atoms connected to current
INTEGER, DIMENSION(:,:), ALLOCATABLE :: c_offs    ! Offsets for atoms connected to current
!
INTEGER :: latom
INTEGER :: i, j
!
!  Accumulate all neighbors of the current atom into temp_neig
!
ktype = cr_iscat(1,katom)
ino_max = get_connectivity_numbers(ktype)  ! Get number of connectivity definitions
DO ino = 1, ino_max                        ! Loop over all connectivity types
   inoo = ino
   CALL get_connectivity_list(katom, ktype, inoo, c_list, c_offs, c_natoms)
   search: DO i=1,c_natoms                         ! Loop over all neighbors
      DO j=1,n_excl
         IF(c_list(i)==excl(j)) CYCLE search
      ENDDO
      IF(.NOT.t_list(c_list(i))) THEN      ! Not yet in the molecule
         t_list(c_list(i)) = .TRUE.        ! Add to temporary list
         latom = c_list(i)                 ! Now search current neighbor for further ones
!        Apply periodic boundary shifts to keep molecule contiguous
         IF(c_offs(1,i)/=0 .OR. c_offs(2,i)/=0 .OR. c_offs(3,i)/=0 ) THEN
            cr_pos(:,latom) = cr_pos(:,latom) + real(c_offs(:,i), kind=PREC_DP)
            CALL conn_update(latom, REAL(c_offs(:,i),kind=PREC_DP))
         ENDIF
         CALL mol_add_excl(latom, n_excl, excl)
      ENDIF
   ENDDO search
   if(allocated(c_list)) DEALLOCATE(c_list)                      ! Remove temporary list of atoms
   if(allocated(c_offs)) DEALLOCATE(c_offs)                      ! Remove temporary offsets
ENDDO
END SUBROUTINE mol_add_excl
!
!*******************************************************************************
!
SUBROUTINE molecularize_numbers(istart,ifinish, new_type, biso, clin, cqua)
!
! Groups atoms istart to ifinish into a molecule
! with molecule type new_type. 
!
USE discus_allocate_appl_mod
USE crystal_mod
USE prop_para_mod
USE errlist_mod
use discus_show_menu
IMPLICIT NONE
!
INTEGER, INTENT(IN)  :: istart
INTEGER, INTENT(IN)  :: ifinish
INTEGER, INTENT(IN)  :: new_type
REAL(kind=PREC_DP)   , INTENT(IN)  :: biso
REAL(kind=PREC_DP)   , INTENT(IN)  :: clin
REAL(kind=PREC_DP)   , INTENT(IN)  :: cqua
!
INTEGER   :: n_new,n_atom, n_type, n_mole
INTEGER   :: i, j
!
IF(istart  < 1 .OR. istart  > cr_natoms .OR.  &
   ifinish < 1 .OR. ifinish > cr_natoms     ) THEN
   ier_num = -137
   ier_typ = ER_APPL
   WRITE(ier_msg(1),'(a4,i8,a13,2i8)') 'n[1]',cr_natoms,' lower/upper ',istart,ifinish
   RETURN
ENDIF
DO i=istart,ifinish
   IF(btest(cr_prop(i),1)) THEN      ! Atom is already inside a molecule
      ier_num = -138
      ier_typ = ER_APPL
      WRITE(ier_msg(1),'(a9,i8)') 'Atom no. ',i
      RETURN
   ENDIF
ENDDO
!
n_new = ifinish - istart + 1            ! no of atoms in new molecule
!
   n_mole = mole_num_mole + 1           ! Reserve new molecule
   n_type = MAX(mole_num_type,new_type) ! If necessary new molecule type number
                                        ! (offset+length) of previous + current length
   n_atom = mole_off(mole_num_mole) + mole_len(mole_num_mole) + n_new

   IF ( n_mole > MOLE_MAX_MOLE  .or.  &
        n_type > MOLE_MAX_TYPE  .or.  &
        n_atom > MOLE_MAX_ATOM      ) THEN
      n_mole = INT(MAX(n_mole,MOLE_MAX_MOLE)*1.25)        ! Add 25%
      n_type = MAX(n_type,MOLE_MAX_TYPE)
      n_atom = INT(MAX(n_atom+n_new,MOLE_MAX_ATOM)*1.25)  ! reserve space for two more molecules+25%
      CALL alloc_molecule(1, 1, n_mole, n_type, n_atom)
      IF ( ier_num /= 0) THEN
         RETURN
      ENDIF
   ENDIF
   mole_num_mole = mole_num_mole + 1
   mole_num_type = MAX(mole_num_type, new_type)
   mole_len(mole_num_mole)  = n_new
   mole_off(mole_num_mole)  = mole_off(mole_num_mole-1)+mole_len(mole_num_mole-1)
   mole_type(mole_num_mole) = new_type
   mole_char(mole_num_mole) = MOLE_ATOM
   mole_file(mole_num_mole) = ' '
   mole_biso(mole_type(mole_num_mole)) = biso
   mole_clin(mole_type(mole_num_mole)) = clin
   mole_cqua(mole_type(mole_num_mole)) = cqua
   mole_num_atom = mole_off(mole_num_mole) + n_new
   j = 0
   DO i=istart,ifinish
         j = j + 1
         mole_cont(mole_off(mole_num_mole)+j) = i
         cr_prop(i) = ibset(cr_prop(i),PROP_MOLECULE)
         cr_mole(i) = mole_num_mole
   ENDDO
!
END SUBROUTINE molecularize_numbers
!
!*******************************************************************************
!
SUBROUTINE molecule_shift(nshift, istart)
!
!   Shifts the atoms in molecule istart and up by nshift atoms further up
!
USE discus_allocate_appl_mod
USE molecule_mod
USE errlist_mod
IMPLICIT NONE
!
INTEGER, INTENT(IN) :: nshift
INTEGER, INTENT(IN) :: istart
!
INTEGER :: i,j,k
INTEGER :: n_mole
INTEGER :: n_type
INTEGER :: n_atom
!
n_atom = MAX(mole_off(mole_num_mole)+mole_len(mole_num_mole)+ nshift, MOLE_MAX_ATOM)
IF(n_atom > MOLE_MAX_ATOM) THEN
   n_mole = MOLE_MAX_MOLE
   n_type = MOLE_MAX_TYPE
   n_atom = INT(MAX(n_atom, MOLE_MAX_ATOM)*1.25)  ! reserve space for two more molecules+25%
   CALL alloc_molecule(1, 1, n_mole, n_type, n_atom)
   IF ( ier_num /= 0) THEN
      RETURN
   ENDIF
ENDIF
!
DO i=mole_num_mole, istart, -1
   DO j=mole_len(i), 1, -1
      k = mole_off(i) + j
      mole_cont(k+nshift) = mole_cont(k)
   ENDDO
   mole_off(i) = mole_off(i) + nshift
ENDDO
IF(istart <= mole_num_mole) THEN
   DO j=1,nshift 
      mole_cont(mole_off(istart-1)+mole_len(istart-1) + j) = 0
   ENDDO
ENDIF
mole_num_atom = mole_num_atom + nshift
!
END SUBROUTINE molecule_shift
!
!*******************************************************************************
!
end module molecule_func_mod
