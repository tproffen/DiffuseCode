MODULE refine_allocate_appl
!-
!     Contains data and routines to allocate application dependent arrays
!
!     SUBROUTINE show_config        ! Shows the current configuration
!
!  The generic procedure alloc_appl acts as interface to:
!
!     SUBROUTINE do_allocate_appl   (zeile, lcomm)  ! Allocate user defined array sizes
!     SUBROUTINE do_deallocate_appl (zeile, lcomm, size) ! Goes back to the default array sizes
!     SUBROUTINE alloc_default                      ! Allocate default array sizes, mostly used at program start
!     SUBROUTINE alloc_params ( n_params)       ! Allocate the number of constraints
!+
USE allocate_generic
USE errlist_mod 
!
PRIVATE
!PUBLIC  :: refine_alloc_appl          ! Generic interface for all allocations
PUBLIC  :: refine_do_allocate_appl
PUBLIC  :: refine_do_deallocate_appl
PUBLIC  :: refine_alloc_default
PUBLIC  :: alloc_params
PUBLIC  :: alloc_params_fix
PUBLIC  :: refine_show_config
!
!
CONTAINS
!
SUBROUTINE refine_do_allocate_appl(zeile,lcomm)
!
USE ber_params_mod
USE get_params_mod
USE precision_mod
USE str_comp_mod
IMPLICIT NONE
!
!
CHARACTER (LEN=*), INTENT(IN)            :: zeile
INTEGER          , INTENT(INOUT)         :: lcomm
!
INTEGER , PARAMETER                      :: MAXW=10
CHARACTER (LEN=MAX(PREC_STRING,LEN(zeile))), DIMENSION(1:MAXW)  :: cpara
INTEGER             , DIMENSION(1:MAXW)  :: lpara
REAL(KIND=PREC_DP)  , DIMENSION(1:MAXW)  :: werte
INTEGER                                  :: ianz
!
!
CALL get_params (zeile, ianz, cpara, lpara, maxw, lcomm) 
IF (ier_num.eq.0) then 
   IF ( ianz == 0 ) THEN
      CALL refine_show_config
   ELSE IF (str_comp (cpara (1) , 'default'   , 1, lpara (1) , 7) )  THEN
      CALL refine_alloc_default
   ELSE IF (str_comp (cpara (1) , 'params', 1, lpara (1) ,10) )  THEN
      IF ( ianz == 2 ) THEN
         CALL del_params (1, ianz, cpara, lpara, maxw) 
         CALL ber_params (ianz, cpara, lpara, werte, maxw) 
         IF (ier_num.eq.0) then 
           CALL alloc_params ( NINT (werte(1)))
         ELSE 
            ier_num = - 6 
            ier_typ = ER_COMM 
         ENDIF 
      ELSE 
         ier_num = - 6 
         ier_typ = ER_COMM 
      ENDIF 
   ELSE IF (str_comp (cpara (1) , 'show', 1, lpara (1) , 4) )  THEN
      CALL refine_show_config
   ELSE 
      ier_num = - 6 
      ier_typ = ER_COMM 
   ENDIF 
ENDIF 
!
END SUBROUTINE refine_do_allocate_appl
!
SUBROUTINE refine_do_deallocate_appl(zeile,lcomm)
!
USE get_params_mod
USE precision_mod
USE str_comp_mod
IMPLICIT NONE
!
!
CHARACTER (LEN=*), INTENT(IN)            :: zeile     ! input command line
INTEGER          , INTENT(INOUT)         :: lcomm     ! command line length
!
INTEGER , PARAMETER                      :: MAXW=10
CHARACTER (LEN=MAX(PREC_STRING,LEN(zeile))), DIMENSION(1:MAXW)  :: cpara
INTEGER             , DIMENSION(1:MAXW)  :: lpara
INTEGER                                  :: ianz
!
!
CALL get_params (zeile, ianz, cpara, lpara, maxw, lcomm)
IF (ier_num.eq.0) then
   IF ( ianz == 0 ) THEN
      ier_num = - 6
      ier_typ = ER_COMM
   ELSE IF (str_comp (cpara (1) , 'all'       , 1, lpara (1) , 3) )  THEN
      CALL dealloc_params
   ELSE IF (str_comp (cpara (1) , 'params', 1, lpara (1) ,10) )  THEN
      CALL dealloc_params
   ELSE
      ier_num = - 6
      ier_typ = ER_COMM
   ENDIF
ENDIF
!
END SUBROUTINE refine_do_deallocate_appl
!
SUBROUTINE refine_show_config
!
!USE constraint
!
USE prompt_mod
!
IMPLICIT NONE
!
!      CONSTRAINTS
!
!IF ( ALLOCATED ( constr_line) .AND. ALLOCATED ( constr_length) ) THEN
!   write(output_io,3100) MAX_CONSTR
!   write(output_io,3110) constr_number
!ELSE
   write(output_io,3500) 
!END IF
!
1100 format(' Maximum number of members     ',i8)
1200 format(' Maximum number of parameters  ',i8)
1300 format(' Current number of members     ',i8)
1400 format(' Current number of parameters  ',i8)
1500 format(' POPULATION is  not allocated yet! REFINE will not run! ')
3100 format(' Maximum number of constraints ',i8)
3110 format(' Current number of constraints ',i8)
3500 format(' CONSTRAINTS are not allocated yet! REFINE will not run! ')
!
END SUBROUTINE refine_show_config
!
SUBROUTINE refine_alloc_default
!
IMPLICIT NONE
!
CALL alloc_params ( 1 )
CALL alloc_params_fix ( 1 )
!     CALL alloc_population ( 1,  1    )
!     CALL alloc_backup     ( 20)
!     CALL alloc_socket_nprogs ( 2, 1)
!
END SUBROUTINE refine_alloc_default
!
!*******************************************************************************
!
SUBROUTINE alloc_params ( n_params)
!-
!     Allocate the number of constraints avaliable to DIFFEV
!+
USE refine_params_mod
!
IMPLICIT NONE
!!
!!      
INTEGER, INTENT(IN)  :: n_params
!!
INTEGER              :: all_status
LOGICAL              :: lstat
!!
lstat = .TRUE.
!!
!    
CALL alloc_arr ( refine_params  ,1,n_params, all_status, ' ')
lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!!
CALL alloc_arr ( refine_fixed   ,1,n_params, all_status, ' ')
lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!!
CALL alloc_arr ( refine_range   ,1,n_params, 1, 2, all_status, 1.00000D0)
lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!!
CALL alloc_arr ( refine_p       ,1,n_params, all_status, 1.00000D0)
lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!!
CALL alloc_arr ( refine_f       ,1,n_params, all_status, 1.00000D0)
lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!!
CALL alloc_arr ( refine_dp      ,1,n_params, all_status, 0.00000D0)
lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!!
CALL alloc_arr ( refine_shift   ,1,n_params, all_status, 0.00500D0)
lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!!
CALL alloc_arr ( refine_nderiv  ,1,n_params, all_status, 2      )
lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!!
CALL alloc_arr ( refine_kderiv  ,1,n_params, all_status, 0      )
lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!!
IF( lstat ) THEN                        ! Success
   REF_MAXPARAM    = n_params
   ier_typ       = 0
   ier_num       = 0
!         IF ( all_status == 1 ) THEN
!            constr_number = 0
!            ier_typ       = 1
!            ier_num       = ER_COMM
!            ier_msg(1)    = 'Constraints'
!         ENDIF
ELSE                                    ! Failure
   REF_MAXPARAM     =  0
!         constr_number  =  0
   ier_num        = -3
   ier_typ        = ER_COMM
   ier_msg(1)     = 'Paramters'
END IF
!!
END SUBROUTINE alloc_params
!
!*******************************************************************************
!
SUBROUTINE alloc_params_fix ( n_params)
!-
!     Allocate the number of constraints avaliable to DIFFEV
!+
USE refine_params_mod
!
IMPLICIT NONE
!!
!!      
INTEGER, INTENT(IN)  :: n_params
!!
INTEGER              :: all_status
LOGICAL              :: lstat
!!
lstat = .TRUE.
!!
!    
CALL alloc_arr ( refine_fixed   ,1,n_params, all_status, ' ')
lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!!
CALL alloc_arr ( refine_f       ,1,n_params, all_status, 1.00000D0)
lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!!
IF( lstat ) THEN                        ! Success
   REF_MAXPARAM_FIX = n_params
   ier_typ       = 0
   ier_num       = 0
!         IF ( all_status == 1 ) THEN
!            constr_number = 0
!            ier_typ       = 1
!            ier_num       = ER_COMM
!            ier_msg(1)    = 'Constraints'
!         ENDIF
ELSE                                    ! Failure
   REF_MAXPARAM_FIX =  0
!         constr_number  =  0
   ier_num        = -3
   ier_typ        = ER_COMM
   ier_msg(1)     = 'Paramters'
END IF
!!
END SUBROUTINE alloc_params_fix
!
!
SUBROUTINE dealloc_params
!!-
!!     Deallocate the number of constraints avaliable to DIFFEV
!!     To avoid possible pitfals with old code, the arrays are simply
!!     reallocated to a size of 1.
!!+
!    USE constraint
IMPLICIT NONE
!!
!      CALL alloc_params ( 1)
!      constr_number =  1
!!
END SUBROUTINE dealloc_params
!
!
END MODULE refine_allocate_appl
