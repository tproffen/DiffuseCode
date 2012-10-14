MODULE allocate_appl
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
!     SUBROUTINE alloc_constraint ( n_constr)       ! Allocate the number of constraints
!     SUBROUTINE alloc_population ( n_pop, n_dimx)  ! Allocate the number of members, parameters
!+
USE allocate_generic
USE config
!
PRIVATE
PUBLIC  :: alloc_appl          ! Generic interface for all allocations
!PUBLIC  :: do_deallocate_appl
PUBLIC  :: show_config
!
INTERFACE  alloc_appl
   MODULE PROCEDURE do_allocate_appl, do_deallocate_appl, alloc_default, &
                    alloc_constraint, alloc_population
END INTERFACE alloc_appl
!
CONTAINS
!
    SUBROUTINE do_allocate_appl(zeile,lcomm)
!
    IMPLICIT NONE
!
    include 'errlist.inc'
!
    CHARACTER (LEN=*), INTENT(IN)            :: zeile
    INTEGER          , INTENT(IN)            :: lcomm
!
    INTEGER , PARAMETER                      :: MAXW=10
    CHARACTER (LEN=1024), DIMENSION(1:MAXW)  :: cpara
    INTEGER             , DIMENSION(1:MAXW)  :: lpara
    REAL                , DIMENSION(1:MAXW)  :: werte
    INTEGER                                  :: ianz
!
    LOGICAL  :: str_comp
!
    CALL get_params (zeile, ianz, cpara, lpara, maxw, lcomm) 
    IF (ier_num.eq.0) then 
       IF ( ianz == 0 ) THEN
          CALL show_config
       ELSE IF (str_comp (cpara (1) , 'default'   , 1, lpara (1) , 7) )  THEN
          CALL alloc_default
       ELSE IF (str_comp (cpara (1) , 'constraint', 1, lpara (1) ,10) )  THEN
          IF ( ianz == 2 ) THEN
             CALL del_params (1, ianz, cpara, lpara, maxw) 
             CALL ber_params (ianz, cpara, lpara, werte, maxw) 
             IF (ier_num.eq.0) then 
               CALL alloc_constraint ( NINT (werte(1)))
             ELSE 
                ier_num = - 6 
                ier_typ = ER_COMM 
             ENDIF 
          ELSE 
             ier_num = - 6 
             ier_typ = ER_COMM 
          ENDIF 
       ELSE IF (str_comp (cpara (1) , 'population', 1, lpara (1) ,10) )  THEN
          IF ( ianz == 3 ) THEN
             CALL del_params (1, ianz, cpara, lpara, maxw) 
             CALL ber_params (ianz, cpara, lpara, werte, maxw) 
             IF (ier_num.eq.0) then 
               CALL alloc_population ( NINT (werte(1)), NINT (werte(2)))
             ELSE 
                ier_num = - 6 
                ier_typ = ER_COMM 
             ENDIF 
          ELSE 
             ier_num = - 6 
             ier_typ = ER_COMM 
          ENDIF 
       ELSE IF (str_comp (cpara (1) , 'show', 1, lpara (1) , 4) )  THEN
          CALL show_config
       ELSE 
          ier_num = - 6 
          ier_typ = ER_COMM 
       ENDIF 
    ENDIF 
!
    END SUBROUTINE do_allocate_appl
!
    SUBROUTINE do_deallocate_appl(zeile,lcomm, flag)
!
       IMPLICIT NONE
!
       include 'errlist.inc'
!
       CHARACTER (LEN=*), INTENT(IN)            :: zeile     ! input command line
       INTEGER          , INTENT(IN)            :: lcomm     ! command line length
       CHARACTER (LEN=*), INTENT(IN)            :: flag      ! flag deallocation
!
       INTEGER , PARAMETER                      :: MAXW=10
       CHARACTER (LEN=1024), DIMENSION(1:MAXW)  :: cpara
       INTEGER             , DIMENSION(1:MAXW)  :: lpara
       REAL                , DIMENSION(1:MAXW)  :: werte
       INTEGER                                  :: ianz
!
       LOGICAL  :: str_comp
!
            CALL get_params (zeile, ianz, cpara, lpara, maxw, lcomm)
            IF (ier_num.eq.0) then
               IF ( ianz == 0 ) THEN
                  ier_num = - 6
                  ier_typ = ER_COMM
               ELSE IF (str_comp (cpara (1) , 'all'       , 1, lpara (1) , 3) )  THEN
                  CALL dealloc_constraint
                  CALL dealloc_population
               ELSE IF (str_comp (cpara (1) , 'constraint', 1, lpara (1) ,10) )  THEN
                  CALL dealloc_constraint
               ELSE IF (str_comp (cpara (1) , 'population', 1, lpara (1) ,10) )  THEN
                  CALL dealloc_population
               ELSE
                  ier_num = - 6
                  ier_typ = ER_COMM
               ENDIF
            ENDIF
!
    END SUBROUTINE do_deallocate_appl
!
    SUBROUTINE show_config
!
       USE constraint
       USE population
!
       IMPLICIT NONE
!
       include 'prompt.inc'
!
!      POPULATION
!
       IF ( ALLOCATED ( pop_x) ) THEN    ! There are many arrays, lets check just one
          write(output_io,1100) MAXPOP
          write(output_io,1200) MAXDIMX
          write(output_io,1300) pop_n
          write(output_io,1400) pop_dimx
       ELSE
          write(output_io,1500) 
       END IF
!
!      CONSTRAINTS
!
       IF ( ALLOCATED ( constr_line) .AND. ALLOCATED ( constr_length) ) THEN
          write(output_io,3100) MAX_CONSTR
          write(output_io,3110) constr_number
       ELSE
          write(output_io,3500) 
       END IF
!
1100 format(' Maximum number of members     ',i8)
1200 format(' Maximum number of parameters  ',i8)
1300 format(' Current number of members     ',i8)
1400 format(' Current number of parameters  ',i8)
1500 format(' POPULATION is  not allocated yet! DIFFEV will not run! ')
3100 format(' Maximum number of constraints ',i8)
3110 format(' Current number of constraints ',i8)
3500 format(' CONSTRAINTS are not allocated yet! DIFFEV will not run! ')
!
    END SUBROUTINE show_config
!
    SUBROUTINE alloc_default
!
      IMPLICIT NONE
!
      CALL alloc_constraint ( 1 )
      CALL alloc_population ( 1,  1    )
!
    END SUBROUTINE alloc_default
!
!
    SUBROUTINE alloc_constraint ( n_constr)
!-
!     Allocate the number of constraints avaliable to DIFFEV
!+
      USE constraint
!
      IMPLICIT NONE
!
      include 'errlist.inc'
!      
      INTEGER, INTENT(IN)  :: n_constr
!
      INTEGER              :: all_status
      LOGICAL              :: lstat
      INTEGER              :: size_of
!
      lstat = .TRUE.
      constr_size_of = 0
!
    
      CALL alloc_arr ( constr_line  ,1,n_constr, all_status, ' ', size_of)
      lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
      constr_size_of =  constr_size_of + size_of
!
      CALL alloc_arr ( constr_length,1,n_constr, all_status, 1, size_of)
      lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
      constr_size_of =  constr_size_of + size_of
!
      IF( lstat ) THEN                        ! Success
         MAX_CONSTR    = n_constr
         ier_typ       = 0
         ier_num       = 0
         IF ( all_status == 1 ) THEN
            constr_number = 0
            ier_typ       = 1
            ier_num       = ER_COMM
            ier_msg(1)    = 'Constraints'
         ENDIF
      ELSE                                    ! Failure
         MAX_CONSTR     =  0
         constr_number  =  0
         constr_size_of =  0
         ier_num        = -2
         ier_typ        = ER_COMM
         ier_msg(1)     = 'Constraints'
         RETURN
      END IF
!
      RETURN
    END SUBROUTINE alloc_constraint
!
    SUBROUTINE alloc_population ( n_pop, n_dimx)
!-
!     Allocate the number of members and parameters avaliable to DIFFEV
!+
      USE population
!
      IMPLICIT NONE
!
      include 'errlist.inc'
!      
      INTEGER, INTENT(IN)  :: n_pop
      INTEGER, INTENT(IN)  :: n_dimx
!
      INTEGER              :: all_status
      LOGICAL              :: lstat
      INTEGER              :: size_of
!
      lstat = .TRUE.
      pop_size_of = 0
!
      CALL alloc_arr ( pop_name      ,1,n_dimx, all_status, 'PARA0000', size_of)
      lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
      pop_size_of = pop_size_of + size_of
!
      CALL alloc_arr ( pop_type      ,1,n_dimx, all_status, POP_REAL, size_of)
      lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
      pop_size_of = pop_size_of + size_of
!
      CALL alloc_arr ( pop_lname     ,1,n_dimx, all_status, 8, size_of)
      lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
      pop_size_of = pop_size_of + size_of
!
      CALL alloc_arr ( pop_refine    ,1,n_dimx, all_status, .false., size_of)
      lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
      pop_size_of = pop_size_of + size_of
!
      CALL alloc_arr ( pop_ad_sigma  ,1,n_dimx, all_status, .true., size_of)
      lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
      pop_size_of = pop_size_of + size_of
!
      CALL alloc_arr ( pop_ad_lsigma ,1,n_dimx, all_status, .true., size_of)
      lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
      pop_size_of = pop_size_of + size_of
!
      CALL alloc_arr ( pop_para      ,1,n_dimx, all_status, 0.0, size_of)
      lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
      pop_size_of = pop_size_of + size_of
!
      CALL alloc_arr ( pop_xmin      ,1,n_dimx, all_status, 0.0, size_of)
      lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
      pop_size_of = pop_size_of + size_of
!
      CALL alloc_arr ( pop_xmax      ,1,n_dimx, all_status, 0.0, size_of)
      lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
      pop_size_of = pop_size_of + size_of
!
      CALL alloc_arr ( pop_pmin      ,1,n_dimx, all_status, 0.0, size_of)
      lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
      pop_size_of = pop_size_of + size_of
!
      CALL alloc_arr ( pop_pmax      ,1,n_dimx, all_status, 0.0, size_of)
      lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
      pop_size_of = pop_size_of + size_of
!
      CALL alloc_arr ( pop_smin      ,1,n_dimx, all_status, 0.0, size_of)
      lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
      pop_size_of = pop_size_of + size_of
!
      CALL alloc_arr ( pop_smax      ,1,n_dimx, all_status, 0.0, size_of)
      lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
      pop_size_of = pop_size_of + size_of
!
      CALL alloc_arr ( pop_sigma     ,1,n_dimx, all_status, 1.0, size_of)
      lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
      pop_size_of = pop_size_of + size_of
!
      CALL alloc_arr ( pop_lsig      ,1,n_dimx, all_status, 0.5, size_of)
      lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
      pop_size_of = pop_size_of + size_of
!
      CALL alloc_arr ( pop_sig_ad    ,1,n_dimx, all_status, 0.5, size_of)
      lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
      pop_size_of = pop_size_of + size_of
!
      CALL alloc_arr ( pop_lsig_ad   ,1,n_dimx, all_status, 0.05, size_of)
      lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
      pop_size_of = pop_size_of + size_of
!
!
      CALL alloc_arr ( child_val     ,1,n_pop , all_status, 0.0, size_of)
      lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
      pop_size_of = pop_size_of + size_of
!
      CALL alloc_arr ( trial_val     ,1,n_pop , all_status, 0.0, size_of)
      lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
      pop_size_of = pop_size_of + size_of
!
      CALL alloc_arr ( parent_val    ,1,n_pop , all_status, 0.0, size_of)
      lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
      pop_size_of = pop_size_of + size_of
!
!
      CALL alloc_arr ( pop_x ,1,n_dimx, 1,n_pop , all_status, 0.0, size_of)
      lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
      pop_size_of = pop_size_of + size_of
!
      CALL alloc_arr ( pop_t ,1,n_dimx, 1,n_pop , all_status, 0.0, size_of)
      lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
      pop_size_of = pop_size_of + size_of
!
      CALL alloc_arr ( child ,1,n_dimx, 1,n_pop , all_status, 0.0, size_of)
      lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
      pop_size_of = pop_size_of + size_of
!
      CALL alloc_arr ( trial ,1,n_dimx, 1,n_pop , all_status, 0.0, size_of)
      lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
      pop_size_of = pop_size_of + size_of
!
      IF( lstat ) THEN                        ! Success
         MAXDIMX       = n_dimx
         MAXPOP        = n_pop
         ier_typ       = 0
         ier_num       = 0
         IF ( all_status == 1 ) THEN
            pop_dimx      = 0
            pop_n         = 0
            pop_c         = 0
            ier_typ       = 1
            ier_num       = ER_COMM
            ier_msg(1)    = 'Population'
         ENDIF
      ELSE                                    ! Failure
         pop_dimx      =  0
         pop_n         =  0
         pop_c         =  0
         pop_size_of   =  0
         ier_num       = -2
         ier_typ       = ER_COMM
         ier_msg(1)    = 'Population'
         RETURN
      END IF
!
    END SUBROUTINE alloc_population
!
    SUBROUTINE dealloc_constraint
!-
!     Deallocate the number of constraints avaliable to DIFFEV
!     To avoid possible pitfals with old code, the arrays are simply
!     reallocated to a size of 1.
!+
    USE constraint
      IMPLICIT NONE
!
      CALL alloc_constraint ( 1)
      constr_number =  1
!
    END SUBROUTINE dealloc_constraint
!
    SUBROUTINE dealloc_population
!-
!     Deallocate the number of members and parameters avaliable to DIFFEV
!     To avoid possible pitfals with old code, the arrays are simply
!     reallocated to a size of 1.
!+
      USE blk_appl
!
      IMPLICIT NONE
!
      CALL alloc_population ( 1, 1)
      CALL init_population
!
    END SUBROUTINE dealloc_population
!
!
END MODULE allocate_appl
