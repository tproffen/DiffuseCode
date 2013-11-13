module discus_wrap
use iso_c_binding
use config_mod
use crystal_mod
use pdf_mod
use pdf_menu
use structur
use allocate_appl_mod
use stack_rese_mod
use update_cr_dim_mod
implicit none

contains
  
  subroutine get_cr_pos(cr_pos_out,n) bind(C)
    real(c_float), intent (out), dimension(3,n) :: cr_pos_out
    integer(c_int), intent(in), value :: n
    cr_pos_out=cr_pos
  end subroutine get_cr_pos
  
  subroutine set_cr_pos(cr_pos_in,n) bind(C)
    real(c_float), intent (in), dimension(3,n) :: cr_pos_in
    integer(c_int), intent(in), value :: n
    cr_pos=cr_pos_in
  end subroutine set_cr_pos
  
  subroutine get_cr_iscat(cr_iscat_out,n) bind(C)
    integer(c_int), intent(out), dimension(n) :: cr_iscat_out
    integer(c_int), value :: n
    !print*,size(cr_iscat),nmax
    !print*,cr_iscat
    cr_iscat_out = cr_iscat
  end subroutine get_cr_iscat
  
  subroutine get_pdf(out,nmi,n) bind(C)
    real(c_double), intent(out), dimension(n) :: out
    integer(c_int), intent(in), value :: nmi, n
    integer :: i,j
    j=0
    do i=nmi,nmi+n-1
       j=j+1
       out(j)=pdf_skal * pdf_calc(i)
    end do
  end subroutine get_pdf
  
  subroutine pdf_show_c(str_in,length) bind(C)
    character(kind=c_char), dimension(length), intent(in) :: str_in
    integer(c_int), value :: length
    character(kind=c_char,len=length) :: str_out
    str_out = transfer(str_in,str_out)
    call pdf_show(str_out)
  end subroutine pdf_show_c
  
  subroutine get_crystal_name(str) bind(C)
    character(kind=c_int), dimension(len_trim(cr_name)+1) :: str
    str = transfer(cr_name,str)
    str(len_trim(cr_name)+1)=c_null_char
  end subroutine get_crystal_name
  
  subroutine get_crystal_spcgr(str) bind(C)
    character(kind=c_int), dimension(len_trim(cr_name)+1) :: str
    str = transfer(cr_spcgr,str)
    str(len_trim(cr_spcgr)+1)=c_null_char    
  end subroutine get_crystal_spcgr
  
  subroutine get_atom_type(str,n) bind(C)
    character(kind=c_char), dimension(5) :: str
    integer(c_int), value :: n
    !print*,size(cr_at_lis),len(cr_at_lis),maxscat
    !print*,cr_at_lis(n)
    str=transfer(cr_at_lis(n),str)
    str(len_trim(cr_at_lis(n))+1)=c_null_char
  end subroutine get_atom_type
  
  subroutine read_cell(fname,fname_length,dim1,dim2,dim3) bind(C)
    character(kind=c_char), dimension(1024), intent(in) :: fname
    integer(c_int), value :: fname_length,dim1,dim2,dim3
    character(kind=c_char,len=1024) :: strucfile
    integer :: natoms,nscats,n_mole,n_type,n_atom,iatom,ce_natoms,ncells,i,j,k,l,n
    strucfile = transfer(fname,strucfile)
    strucfile(fname_length+1:) = ' '
    !print*,strucfile
    call rese_cr()
    cr_icc(1)=dim1
    cr_icc(2)=dim2
    cr_icc(3)=dim3
    call test_file(strucfile, natoms, nscats, n_mole, n_type, &
         n_atom, -1 , .false.)
    call readcell(strucfile)
    iatom = cr_icc (1) * cr_icc (2) * cr_icc (3) * cr_natoms
    if (iatom.gt.nmax) then
       call alloc_crystal ( maxscat, int(iatom * 1.1))
    endif
    ce_natoms = cr_natoms
    cr_ncatoms = cr_natoms
    cr_natoms = 0
    do k = 1, cr_icc (3)
       do j = 1, cr_icc (2)
          do i = 1, cr_icc (1)
             do n = 1, ce_natoms
                cr_natoms = cr_natoms + 1 
                cr_iscat (cr_natoms) = cr_iscat (n)
                cr_pos (1, cr_natoms) = cr_pos (1, n) + float ( i - 1)
                cr_pos (2, cr_natoms) = cr_pos (2, n) + float ( j - 1)
                cr_pos (3, cr_natoms) = cr_pos (3, n) + float ( k - 1)
                cr_prop (cr_natoms) = cr_prop (n)
             enddo
          enddo
       enddo
    enddo
    call update_cr_dim
    ncells = cr_icc (1) * cr_icc (2)* cr_icc (3)
    do l = 1, 3
       cr_dim0 (l, 1) = float (nint (cr_dim (l, 1) ) )
       cr_dim0 (l, 2) = float (nint (cr_dim (l, 2) ) )
    enddo
    call do_stack_rese    
  end subroutine read_cell
  
  subroutine alloc_pdf_f() BIND(C)
    CALL alloc_pdf( pdf_nscat, pdf_ndat, pdf_nbnd )
  end subroutine alloc_pdf_f
end module discus_wrap
