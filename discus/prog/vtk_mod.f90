module vtk_mod
  
  contains
  
  subroutine vtk_write (value, laver)
    use diffuse_mod
    use crystal_mod
    use output_mod
    implicit none
    integer, intent(in) :: value
    logical, intent(in) :: laver
    real                :: dx,dy,dz
    
    dx=sqrt(out_vi(1,1)**2+out_vi(2,1)**2+out_vi(3,1)**2)
    dy=sqrt(out_vi(1,2)**2+out_vi(2,2)**2+out_vi(3,2)**2)
    dz=sqrt(out_vi(1,3)**2+out_vi(2,3)**2+out_vi(3,3)**2)
    
    open(10,file=outfile,form='FORMATTED',&
         access='SEQUENTIAL',action='WRITE',status='REPLACE')
    write(10,'(A)')'# vtk DataFile Version 3.0'
    write(10,'(A)')'DISCUS OUTPUT'
    write(10,'(A)')'ASCII'
    write(10,'(A)')'DATASET STRUCTURED_POINTS'
    write(10,'(A,3(X,I0))')'DIMENSIONS',out_inc(1),out_inc(2),out_inc(3)
    write(10,'(A,3(X,F0.4))')'ORIGIN',out_eck(1,1),out_eck(2,1),out_eck(3,1)
    write(10,'(A,3(X,F0.4))')'SPACING',dx,dy,dz
    write(10,'(A,X,I0)'),'POINT_DATA',out_inc(1)*out_inc(2)*out_inc(3)
    write(10,'(A)')'SCALARS values float'
    write(10,'(A)')'LOOKUP_TABLE default'
    write(10,'(E13.6E2)') dsi
    
    close(10)
    
  end subroutine vtk_write

end module vtk_mod
