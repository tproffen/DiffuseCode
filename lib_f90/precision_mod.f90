MODULE precision_mod
!
use iso_fortran_env
!
INTEGER, PARAMETER:: PREC_INT_BYTE  = 1
INTEGER, PARAMETER:: PREC_INT_SHORT = 2
INTEGER, PARAMETER:: PREC_INT_WORD  = 4
INTEGER, PARAMETER:: PREC_INT_LONG  = 8
INTEGER, PARAMETER:: PREC_INT_LARGE=MAX(SELECTED_INT_KIND(8) , &
                                        SELECTED_INT_KIND(16) ) 
INTEGER, PARAMETER:: PREC_SP=SELECTED_REAL_KIND(p= 6      )  ! single precision
INTEGER, PARAMETER:: PREC_DP=SELECTED_REAL_KIND(p=15,r=307)  ! double precision
INTEGER, PARAMETER:: PREC_QP=SELECTED_REAL_KIND(p=30,r=307)  ! quad   precision
INTEGER, PARAMETER:: PREC_HP=SELECTED_REAL_KIND(p=30,r=607)  ! quad   precision
!
INTEGER, PARAMETER:: PREC_STRING  = 1024
INTEGER, PARAMETER:: PREC_LSTRING = 2048
!
integer(kind=PREC_INT_BYTE), parameter :: PREC_IBYTE  = 1
integer(kind=PREC_INT_BYTE), parameter :: PREC_ISHORT = 1
integer(kind=PREC_INT_BYTE), parameter :: PREC_IWORD  = 1
integer(kind=PREC_INT_BYTE), parameter :: PREC_ILONG  = 1
integer(kind=PREC_INT_BYTE), parameter :: PREC_ILARGE = 1
integer(kind=PREC_INT_BYTE), parameter :: PREC_HUGE_IBYTE  = huge(PREC_IBYTE )
integer(kind=PREC_INT_BYTE), parameter :: PREC_HUGE_ISHORT = huge(PREC_ISHORT)
integer(kind=PREC_INT_BYTE), parameter :: PREC_HUGE_IWORD  = huge(PREC_IWORD )
integer(kind=PREC_INT_BYTE), parameter :: PREC_HUGE_ILONG  = huge(PREC_ILONG )
integer(kind=PREC_INT_BYTE), parameter :: PREC_HUGE_ILARGE = huge(PREC_ILARGE)
!
END MODULE precision_mod
