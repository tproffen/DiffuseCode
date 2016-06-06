MODULE precision_mod
!
INTEGER, PARAMETER:: PREC_INT_LARGE=MAX(SELECTED_INT_KIND(8) , &
                                        SELECTED_INT_KIND(16) ) 
INTEGER, PARAMETER:: PREC_DP=KIND(0.D0)  ! double precision
END MODULE precision_mod
