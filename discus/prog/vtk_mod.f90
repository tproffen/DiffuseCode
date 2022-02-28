module vtk_mod
  contains
  
  subroutine vtk_write ()
    use diffuse_mod
    use crystal_mod
    use output_mod
    use errlist_mod
use precision_mod
USE support_mod
    implicit none
    real(kind=PREC_DP)  :: dx,dy,dz
    
    dx=sqrt(out_vi(1,1)**2+out_vi(2,1)**2+out_vi(3,1)**2)
    dy=sqrt(out_vi(1,2)**2+out_vi(2,2)**2+out_vi(3,2)**2)
    dz=sqrt(out_vi(1,3)**2+out_vi(2,3)**2+out_vi(3,3)**2)
    
    call oeffne (2, outfile, 'unknown')
    if (ier_num.ne.0) return
    
    write(2,'(A)') '# vtk DataFile Version 3.0'
    write(2,'(A)') 'DISCUS OUTPUT'
    write(2,'(A)') 'ASCII'
    write(2,'(A)') 'DATASET STRUCTURED_POINTS'
    write(2,'(A,3(1X,I0))') 'DIMENSIONS',out_inc(1),out_inc(2),out_inc(3)
    write(2,'(A,3(1X,F0.4))') 'ORIGIN',out_eck(1,1),out_eck(2,1),out_eck(3,1)
    write(2,'(A,3(1X,F0.4))') 'SPACING',dx,dy,dz
    write(2,'(A,1X,I0)') 'POINT_DATA',out_inc(1)*out_inc(2)*out_inc(3)
    write(2,'(A)') 'SCALARS values float'
    write(2,'(A)') 'LOOKUP_TABLE default'
    write(2,'(E13.6E2)') dsi
    
    close(2)
    
  end subroutine vtk_write

end module vtk_mod
