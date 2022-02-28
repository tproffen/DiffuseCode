MODULE powder_tables_mod
!
USE errlist_mod
!
IMPLICIT NONE
!
PUBLIC
!
INTEGER, DIMENSION(:), ALLOCATABLE :: powder_istl   ! Look up for sin(theta)/lambda
!
CONTAINS
!*****7*****************************************************************
      SUBROUTINE powder_sinet 
!+                                                                      
!     This routine initialises the real sine lookup table and           
!     is called only at the first Powder run.                           
!-                                                                      
      USE discus_config_mod 
      USE debye_mod 
      USE diffuse_mod 
!                                                                       
      USE prompt_mod 
      USE precision_mod
      IMPLICIT none 
!                                                                       
      REAL(PREC_DP) twopi, xmult, xarg, xt 
      INTEGER i 
!                                                                       
      WRITE (output_io, 1000) 
!                                                                       
      xt = 1.0d0 / REAL(I2PI, KIND=KIND(0.0D0)) 
      twopi = 8.0d0 * atan(1.0d0) 
!                                                                       
      sinetab(0) = 0.0D0
      DO i = 1, MASK 
         xmult       = REAL(i, KIND=KIND(0.0D0)) * xt 
         xarg        = twopi * xmult 
         sinetab (i) = DBLE( dsin (xarg)        )
      ENDDO 
      ffour = .false. 
!                                                                       
 1000 FORMAT     (' Computing powder sine lookup table ...') 
      END SUBROUTINE powder_sinet                   
!*****7*****************************************************************
      SUBROUTINE powder_stltab (n_points, xstart, xdelta)
!+                                                                      
!     Sets up an integer array containing the corresponding integers    
!     to the formfactor table for each sin(theta)/lambda.               
!-                                                                      
      USE discus_config_mod 
      USE discus_allocate_appl_mod
      USE crystal_mod 
      USE debye_mod
      USE diffuse_mod 
!                                                                       
      USE prompt_mod 
      USE precision_mod 
      IMPLICIT none 
!
      INTEGER, INTENT(IN) :: n_points   ! number of points in reciprocal space
      REAL(PREC_DP), INTENT(IN) :: xstart     ! starting point in reciprocal space
      REAL(PREC_DP), INTENT(IN) :: xdelta     ! starting point in reciprocal space
!                                                                       
      REAL(kind=PREC_DP) :: q2 
      INTEGER i
      INTEGER     :: n_qxy   = 1 ! Maximum number of points in reciprocal space 
      INTEGER     :: n_nscat = 1 ! Maximum Number of atom type
      INTEGER     :: n_natom = 1 ! Maximum Number of atom type
      INTEGER     :: astatus ! Allocation status
!
!      n_qxy   = MAX(1,n_points)
      n_qxy   = MAX(  n_points      , MAXQXY,MAXDQXY)
!     n_nscat = 1
!     n_natom = 1
!                                                                       
!     IF (four_log) then 
         WRITE (output_io, 1000) 
!     ENDIF 
!                                                                       
      IF (n_points > MAX(MAXQXY, MAXDQXY)  .OR.          &
          n_points > ubound(istl,1)        .or.          &
          cr_nscat>DIF_MAXSCAT              ) THEN
         n_nscat = MAX(n_nscat,cr_nscat,DIF_MAXSCAT)
         n_natom = MAX(n_natom,cr_natoms,DIF_MAXAT)
         CALL alloc_diffuse (n_qxy,  n_nscat,  n_natom )
         IF (ier_num.ne.0) THEN
            RETURN
         ENDIF
      ENDIF
      IF(ALLOCATED(powder_istl)) DEALLOCATE(powder_istl)
      ALLOCATE(powder_istl(0:n_qxy), stat=astatus)
      powder_istl(:) = 0
!
      DO i = 0, n_points
         q2 = ( (xstart + xdelta  * REAL(i,KIND=KIND(0.0D0)) ) **2) / 4.0D0
         istl (i) = nint(sqrt(q2) * (1.0D0 / CFINC) ) 
!DBG                                                                    
!DBG      if(i.eq.150) then                                             
!DBG        write(*,*) 'q2  ', q2                                       
!DBG        write(*,*) 'q   ', sqrt(q2)                                 
!DBG        write(*,*) 'istl',istl(i)                                   
!DBG      endif                                                         
!                                                                       
         IF (istl (i) .gt.CFPKT) then 
            ier_num = - 3 
            ier_typ = ER_FOUR 
            RETURN 
         ENDIF 
!                                                                       
      ENDDO 
      powder_istl(:) = istl(:)
!                                                                       
 1000 FORMAT     (' Computing sin(theta)/lambda table ...') 
      END SUBROUTINE powder_stltab                  
!*****7*****************************************************************
END MODULE powder_tables_mod
