MODULE diffev_allocate_appl
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
USE diffev_config
USE errlist_mod 
!
PRIVATE
!PUBLIC  :: diffev_alloc_appl          ! Generic interface for all allocations
PUBLIC  :: diffev_do_allocate_appl
PUBLIC  :: diffev_do_deallocate_appl
PUBLIC  :: diffev_alloc_default
PUBLIC  :: alloc_backup
PUBLIC  :: dealloc_backup
PUBLIC  :: alloc_socket_nprogs
PUBLIC  :: alloc_population
PUBLIC  :: alloc_constraint
PUBLIC  :: diffev_show_config
!
!
CONTAINS
!
    SUBROUTINE diffev_do_allocate_appl(zeile,lcomm)
!
USE ber_params_mod
    USE get_params_mod
USE precision_mod
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
    LOGICAL  :: str_comp
!
    CALL get_params (zeile, ianz, cpara, lpara, maxw, lcomm) 
    IF (ier_num.eq.0) then 
       IF ( ianz == 0 ) THEN
          CALL diffev_show_config
       ELSE IF (str_comp (cpara (1) , 'default'   , 1, lpara (1) , 7) )  THEN
          CALL diffev_alloc_default
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
          CALL diffev_show_config
       ELSE 
          ier_num = - 6 
          ier_typ = ER_COMM 
       ENDIF 
    ENDIF 
!
    END SUBROUTINE diffev_do_allocate_appl
!
    SUBROUTINE diffev_do_deallocate_appl(zeile,lcomm)
!
       USE get_params_mod
USE precision_mod
       IMPLICIT NONE
!
!
       CHARACTER (LEN=*), INTENT(IN)            :: zeile     ! input command line
       INTEGER          , INTENT(INOUT)         :: lcomm     ! command line length
!
       INTEGER , PARAMETER                      :: MAXW=10
       CHARACTER (LEN=MAX(PREC_STRING,LEN(zeile))), DIMENSION(1:MAXW)  :: cpara
       INTEGER             , DIMENSION(1:MAXW)  :: lpara
!      REAL                , DIMENSION(1:MAXW)  :: werte
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
    END SUBROUTINE diffev_do_deallocate_appl
!
    SUBROUTINE diffev_show_config
!
       USE constraint
       USE population
!
       USE prompt_mod
!
       IMPLICIT NONE
!
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
    END SUBROUTINE diffev_show_config
!
    SUBROUTINE diffev_alloc_default
!
      IMPLICIT NONE
!
      CALL alloc_constraint ( 1 )
      CALL alloc_population ( 1,  1    )
      CALL alloc_backup     ( 20)
      CALL alloc_socket_nprogs ( 2, 1)
!
    END SUBROUTINE diffev_alloc_default
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
         ier_num        = -3
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
      USE random_state_mod
      USE variable_mod
!
      IMPLICIT NONE
!
!      
      INTEGER, INTENT(IN)  :: n_pop
      INTEGER, INTENT(IN)  :: n_dimx
!
      INTEGER              :: all_status
      INTEGER              :: nseeds
      LOGICAL              :: lstat
      INTEGER              :: size_of
!
      lstat = .TRUE.
      pop_size_of = 0
!
      CALL alloc_arr ( pop_name      ,0,n_dimx, all_status, 'PARA0000', size_of)
      lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
      pop_name(0) = 'Rvalue'
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
      CALL alloc_arr ( child_val     ,1,n_pop , 0,15,all_status, 0.0, size_of)
      lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
      pop_size_of = pop_size_of + size_of
!
      CALL alloc_arr ( trial_val     ,1,n_pop , 0,15,all_status, 0.0, size_of)
      lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
      pop_size_of = pop_size_of + size_of
!
      CALL alloc_arr ( parent_val    ,1,n_pop , 0,15,all_status, 0.0, size_of)
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
      nseeds = random_nseeds()                ! too be improved for different compilers
!     nseeds = 12
      CALL alloc_arr ( pop_random ,0,nseeds, 1,n_pop , all_status, 0  , size_of)
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
            var_val(var_ref+1) = 0  ! Update global user variable
            var_val(var_ref+2) = 0  ! Update global user variable
            var_val(var_ref+3) = 0  ! Update global user variable
            ier_typ       = 1
            ier_num       = ER_COMM
            ier_msg(1)    = 'Population'
         ENDIF
      ELSE                                    ! Failure
         pop_dimx      =  0
         pop_n         =  0
         pop_c         =  0
         var_val(var_ref+1) = 0  ! Update global user variable
         var_val(var_ref+2) = 0  ! Update global user variable
         var_val(var_ref+3) = 0  ! Update global user variable
         pop_size_of   =  0
         ier_num       = -3
         ier_typ       = ER_COMM
         ier_msg(1)    = 'Population'
         RETURN
      END IF
!
    END SUBROUTINE alloc_population
!
!
    SUBROUTINE alloc_backup ( n_pop)
!-
!     Allocate the number of backups     avaliable to DIFFEV
!+
      USE population
!
      IMPLICIT NONE
!
!      
      INTEGER, INTENT(IN)  :: n_pop
!
      INTEGER              :: all_status
      LOGICAL              :: lstat
      INTEGER              :: size_of
!
      lstat = .TRUE.
!
!
      CALL alloc_arr ( pop_back_fil  ,1,n_pop, all_status, ' ', size_of)
      lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
      CALL alloc_arr ( pop_back_ext  ,1,n_pop, all_status, ' ', size_of)
      lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
      CALL alloc_arr ( pop_back_trg  ,1,n_pop, all_status, ' ', size_of)
      lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
      CALL alloc_arr ( pop_back_fil_l,1,n_pop, all_status, 1, size_of)
      lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
      CALL alloc_arr ( pop_back_ext_l,1,n_pop, all_status, 1, size_of)
      lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
      CALL alloc_arr ( pop_back_trg_l,1,n_pop, all_status, 1, size_of)
      lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
      IF( lstat ) THEN                        ! Success
         MAXBACK       = n_pop 
         ier_typ       = 0
         ier_num       = 0
         IF ( all_status == 1 ) THEN
            pop_back_number = 0
            ier_typ       = 1
            ier_num       = ER_COMM
            ier_msg(1)    = 'Backups'
         ENDIF
      ELSE                                    ! Failure
         MAXBACK        =  0
         pop_back_number=  0
         ier_num        = -3
         ier_typ        = ER_COMM
         ier_msg(1)     = 'Backups'
         RETURN
      END IF
!
      RETURN
    END SUBROUTINE alloc_backup
!
    SUBROUTINE alloc_socket_nprogs ( nprog, nproc)
!-
!     Allocate the number of programs that may be started via sockets
!+
      USE run_mpi_mod
!
      IMPLICIT NONE
!
!      
      INTEGER, INTENT(IN)  :: nprog
      INTEGER, INTENT(IN)  :: nproc
!
      INTEGER              :: all_status
      LOGICAL              :: lstat
      INTEGER              :: size_of
!
      lstat = .TRUE.
      progs_size_of = 0
!
      CALL alloc_arr ( prog_entry ,1,nprog, all_status, ' ', size_of)
      lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
      progs_size_of =   progs_size_of + size_of
!
      CALL alloc_arr ( socket_id ,1,nprog, all_status, 0  , size_of)
      lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
      progs_size_of =   progs_size_of + size_of
!
      CALL alloc_arr (   port_id ,1,nprog, 1, nproc, all_status, 0  , size_of)
      lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
      progs_size_of =   progs_size_of + size_of
!
      IF( lstat ) THEN                        ! Success
         RUN_MPI_MAXPROG = nprog
         ier_typ       = 0
         ier_num       = 0
         IF ( all_status == 1 ) THEN
            ier_typ       = 1
            ier_num       = ER_COMM
            ier_msg(1)    = 'Socket_programs'
         ENDIF
      ELSE                                    ! Failure
         RUN_MPI_MAXPROG =  0
         progs_size_of =   0
         ier_num         = -3
         ier_typ         = ER_COMM
         ier_msg(1)      = 'Socket_programs'
         RETURN
      END IF
!
      RETURN
    END SUBROUTINE alloc_socket_nprogs
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
      USE diffev_blk_appl
!
      IMPLICIT NONE
!
      CALL alloc_population ( 1, 1)
      CALL init_population
!
    END SUBROUTINE dealloc_population
!
    SUBROUTINE dealloc_backup
!-
!     Deallocate the number of constraints avaliable to DIFFEV
!     To avoid possible pitfals with old code, the arrays are simply
!     reallocated to a size of 1.
!+
    USE population
      IMPLICIT NONE
!
      CALL alloc_backup ( 1)
      pop_back_number =  0
!
    END SUBROUTINE dealloc_backup
!
!
END MODULE diffev_allocate_appl
