MODULE arrays_mod
!
! Performs matrix operations on user defined variables
!
IMPLICIT NONE
!
CONTAINS
!
!*******************************************************************************
!
SUBROUTINE arr_matmul(line, length)
!
! Multiplies two user matrices, or a user matrix with a scalar
!
USE ber_params_mod
USE errlist_mod
USE get_params_mod
USE variable_mod
!
IMPLICIT NONE
!
CHARACTER(LEN=*), INTENT(INOUT) :: line
INTEGER         , INTENT(INOUT) :: length
!
INTEGER, PARAMETER :: MAXP = 3
CHARACTER(LEN=1024), DIMENSION(MAXP) :: cpara
INTEGER            , DIMENSION(MAXP) :: lpara
INTEGER            , DIMENSION(MAXP) :: look
REAL               , DIMENSION(MAXP) :: werte
INTEGER :: ianz, iianz
INTEGER :: i, j, k, ndel
REAL    :: scalar
!
CALL get_params (line, ianz, cpara, lpara, MAXP, length)
!
look(:) = -1
IF(ianz==2) THEN                       ! Single input operation
   cpara(3) = '1.0'
   lpara(3) = 3
   ianz     = 3
ENDIF
IF(ianz==3) THEN
   locate: DO i=1, ianz                ! Test is parameters are user variables
      search: DO j=var_sys+1,var_num
         IF(cpara(i)==var_name(j)) THEN
            look(i) = j
            CYCLE locate
         ENDIF
      ENDDO search
   ENDDO locate
   IF(look(1)==-1) THEN
      ier_num = -24
      ier_typ = ER_FORT
      ier_msg(1) = 'Offending variable '//cpara(i)(1:lpara(i))
      RETURN
   ENDIF
!
!  Try to calculate parameters 2 and 3 to detect single valued expressions
!  After the calculations scalar holds the product of any scalar valued 
!  expression in params 2 and 3
!  The options are reduced to:
!  scalar * scalar
!  scalar * matrix
!  matrix * scalar
!  matrix * matrix
!
   werte(:) = 0.0
   scalar = 1.0
   ndel   = 1
   CALL del_params (ndel, ianz, cpara, lpara, MAXP)
   iianz = 1
   CALL ber_params(iianz, cpara, lpara, werte, MAXP)
   IF(ier_num==0) THEN       ! Could calculate cpara (2)
      look(2) = -1           ! Cancel cpara(2)
      scalar = scalar*werte(1)
   ELSEIF(look(2)==-1) THEN  ! Error on cpara(2) but no matrix 
      RETURN
   ENDIF
   ndel   = 1
   CALL del_params (ndel, ianz, cpara, lpara, MAXP)
   iianz = 1
   CALL ber_params(iianz, cpara, lpara, werte, MAXP)
   IF(ier_num==0) THEN       ! Could calculate cpara (3)
      look(3) = -1           ! Cancel cpara(3)
      scalar = scalar*werte(1)
   ELSEIF(look(3)==-1) THEN  ! Error on cpara(3) but no matrix 
      RETURN
   ENDIF
   CALL no_error
!
   ier_num = -44
   ier_typ = ER_FORT
   IF(look(2)==-1 .AND. look(3)==-1) THEN                            ! scalar * scalar
      var_field(var_entry(look(1)))%var_value(:,:) = scalar
      ier_num = 0
      ier_typ = ER_NONE
!
   ELSEIF(look(2)==-1 .AND. look(3)>0) THEN
      IF(var_entry(look(3))>0) THEN   ! scalar * matrix
      IF(var_field(var_entry(look(1)))%var_shape(1)==var_field(var_entry(look(3)))%var_shape(1) .AND. &
         var_field(var_entry(look(1)))%var_shape(2)==var_field(var_entry(look(3)))%var_shape(2) ) THEN
      DO i=1, var_field(var_entry(look(1)))%var_shape(1)
         DO k=1, var_field(var_entry(look(1)))%var_shape(2)
            var_field(var_entry(look(1)))%var_value(i,k) =   &
            scalar           * var_field(var_entry(look(3)))%var_value(i,k)
         ENDDO
      ENDDO
      ier_num = 0
      ier_typ = ER_NONE
      ENDIF
      ENDIF
   ELSEIF(look(3)==-1 .AND. look(2)>0) THEN
      IF(var_entry(look(2))>0) THEN   ! matrix * scalar
      IF(var_field(var_entry(look(1)))%var_shape(1)==var_field(var_entry(look(2)))%var_shape(1) .AND. &
         var_field(var_entry(look(1)))%var_shape(2)==var_field(var_entry(look(2)))%var_shape(2) ) THEN
      DO i=1, var_field(var_entry(look(1)))%var_shape(1)
         DO k=1, var_field(var_entry(look(1)))%var_shape(2)
            var_field(var_entry(look(1)))%var_value(i,k) =   &
            scalar           * var_field(var_entry(look(2)))%var_value(i,k)
         ENDDO
      ENDDO
      ier_num = 0
      ier_typ = ER_NONE
      ENDIF
      ENDIF
!                                                                        ! matrix * matrix
   ELSEIF(look(2)>0 .AND. look(3)>0 .AND. var_entry(look(2))>0 .AND. var_entry(look(3))>0) THEN

   IF(var_field(var_entry(look(1)))%var_shape(1)==var_field(var_entry(look(2)))%var_shape(1) .AND. &
      var_field(var_entry(look(1)))%var_shape(2)==var_field(var_entry(look(3)))%var_shape(2) .AND. &
      var_field(var_entry(look(2)))%var_shape(2)==var_field(var_entry(look(3)))%var_shape(1) )   THEN
      DO i=1, var_field(var_entry(look(1)))%var_shape(1)
         DO k=1, var_field(var_entry(look(1)))%var_shape(2)
            var_field(var_entry(look(1)))%var_value(i,k) = 0
            DO j=1,var_field(var_entry(look(2)))%var_shape(2)
               var_field(var_entry(look(1)))%var_value(i,k) =    &
                 var_field(var_entry(look(1)))%var_value(i,k)    &
               + var_field(var_entry(look(2)))%var_value(i,j) *  &
                 var_field(var_entry(look(3)))%var_value(j,k)
            ENDDO
         ENDDO
      ENDDO
      ier_num = 0
      ier_typ = ER_NONE
   ENDIF
   ELSE
      ier_num = -44
      ier_typ = ER_FORT
      RETURN
   ENDIF
!
ELSE
   ier_num = -6
   ier_typ = ER_COMM
ENDIF
!
END SUBROUTINE arr_matmul
!
!*******************************************************************************
!
SUBROUTINE arr_matadd(line, length)
!
! Multiplies two user matrices, or a user matrix with a scalar
!
USE ber_params_mod
USE errlist_mod
USE get_params_mod
USE variable_mod
!
IMPLICIT NONE
!
CHARACTER(LEN=*), INTENT(INOUT) :: line
INTEGER         , INTENT(INOUT) :: length
!
INTEGER, PARAMETER :: MAXP = 4
INTEGER, PARAMETER :: MAXS = 1
CHARACTER(LEN=1024), DIMENSION(MAXP) :: cpara
INTEGER            , DIMENSION(MAXP) :: lpara
INTEGER            , DIMENSION(MAXP) :: look
CHARACTER(LEN=1024), DIMENSION(MAXS) :: ccpara
INTEGER            , DIMENSION(MAXS) :: llpara
REAL               , DIMENSION(MAXS) :: wwerte
INTEGER :: ianz, iianz
INTEGER :: i, j, k
REAL    :: scalar
!
CALL get_params (line, ianz, cpara, lpara, MAXP, length)
!
look(:) = -1
IF(ianz==3) THEN                       ! Single input operation
   cpara(4) = cpara(3)
   lpara(4) = lpara(3)
   cpara(3) = '1.0'
   lpara(3) = 3
   ianz     = 4
   iianz    = 0                        ! Signal that scalar is known
   scalar   = 1.0
ELSE
   iianz    = 1
   scalar   = 1.0
ENDIF
IF(ianz==4) THEN
   locate: DO i=1, ianz                ! Test is parameters are user variables
      IF(i==3) CYCLE locate            ! Para 3 is a scalar factor
      search: DO j=var_sys+1,var_num
         IF(cpara(i)==var_name(j)) THEN
            look(i) = j
            CYCLE locate
         ENDIF
      ENDDO search
      IF(look(i)==-1) THEN
         ier_num = -24
         ier_typ = ER_FORT
         ier_msg(1) = 'Offending variable '//cpara(i)(1:lpara(i))
         RETURN
      ENDIF
   ENDDO locate
!
   IF(iianz==1) THEN
!
!  Try to calculate parameter 3 to detect single valued expressions
!
      ccpara(1) = cpara(3)
      llpara(1) = lpara(3)
      wwerte(:) = 0.0
      scalar = 1.0
      iianz = 1
      CALL ber_params(iianz, ccpara, llpara, wwerte, MAXS)
      IF(ier_num==0) THEN       ! Could calculate cpara (2)
         look(3) = -1           ! Cancel cpara(2)
         scalar = wwerte(1)
      ELSE                      ! Error on cpara(3)
         RETURN
      ENDIF
      CALL no_error
   ENDIF
!
!  IF(look(1)>0 .AND. look(2)>0 .AND. look(4)>0) THEN     ! We have three variable
   IF(var_entry(look(1))>0 .AND. var_entry(look(2))>0 .AND. var_entry(look(4))>0) THEN ! three matrices
      IF(var_field(var_entry(look(1)))%var_shape(1)==var_field(var_entry(look(2)))%var_shape(1) .AND. &
         var_field(var_entry(look(1)))%var_shape(2)==var_field(var_entry(look(2)))%var_shape(2) .AND. &
         var_field(var_entry(look(1)))%var_shape(1)==var_field(var_entry(look(4)))%var_shape(1) .AND. &
         var_field(var_entry(look(1)))%var_shape(2)==var_field(var_entry(look(4)))%var_shape(2) ) THEN
         DO i=1, var_field(var_entry(look(1)))%var_shape(1)
            DO k=1, var_field(var_entry(look(1)))%var_shape(2)
               var_field(var_entry(look(1)))%var_value(i,k) =                  &
                                  var_field(var_entry(look(2)))%var_value(i,k) + &
               scalar           * var_field(var_entry(look(4)))%var_value(i,k)
            ENDDO
         ENDDO
      ELSE
         ier_num = -44
         ier_typ = ER_FORT
      ENDIF
   ELSE
      ier_num = -48
      ier_typ = ER_FORT
   ENDIF
!
ELSE
   ier_num = -6
   ier_typ = ER_COMM
ENDIF
!
END SUBROUTINE arr_matadd
!
!******************************************************************************
!
SUBROUTINE arr_detmat(line, length)
!
!  Calculate the determinant of a matrix, must be [n x n]
!
USE errlist_mod
USE get_params_mod
USE matrix_mod
USE param_mod
USE variable_mod
!
IMPLICIT NONE
!
CHARACTER(LEN=*), INTENT(INOUT) :: line
INTEGER         , INTENT(INOUT) :: length
!
INTEGER, PARAMETER :: MAXP = 1
CHARACTER(LEN=1024), DIMENSION(MAXP) :: cpara
INTEGER            , DIMENSION(MAXP) :: lpara
INTEGER            , DIMENSION(MAXP) :: look
INTEGER :: ianz
INTEGER :: i, j
!
CALL get_params (line, ianz, cpara, lpara, MAXP, length)
!
look(:) = -1
IF(ianz==1) THEN
   locate: DO i=1, ianz
      search: DO j=var_sys+1,var_num
         IF(cpara(i)==var_name(j)) THEN
            look(i) = j
            CYCLE locate
         ENDIF
      ENDDO search
         ier_num = -24
         ier_typ = ER_FORT
         ier_msg(1) = 'Offending variable '//cpara(i)(1:lpara(i))
         RETURN
   ENDDO locate
!
   IF(var_field(var_entry(look(1)))%var_shape(1)==var_field(var_entry(look(1)))%var_shape(2) ) THEN 
      IF(var_field(var_entry(look(1)))%var_shape(1)==1) THEN
         res_para(1) = var_field(var_entry(look(1)))%var_value(1,1)
         res_para(0) = 1
      ELSEIF(var_field(var_entry(look(1)))%var_shape(1)==2) THEN
         res_para(1) = det2(var_field(var_entry(look(1)))%var_value)
         res_para(0) = 1
      ELSEIF(var_field(var_entry(look(1)))%var_shape(1)==3) THEN
         res_para(1) = det3(var_field(var_entry(look(1)))%var_value)
         res_para(0) = 1
      ELSEIF(var_field(var_entry(look(1)))%var_shape(1)==4) THEN
         res_para(1) = det4(var_field(var_entry(look(1)))%var_value)
         res_para(0) = 1
      ENDIF
   ELSE
      ier_num = -44
      ier_typ = ER_FORT
      RETURN
   ENDIF
!
ELSE
   ier_num = -6
   ier_typ = ER_COMM
ENDIF
!
END SUBROUTINE arr_detmat
!
!******************************************************************************
!
SUBROUTINE arr_invmat(line, length)
!
!  Calculate the inverse matrix
!
USE errlist_mod
USE get_params_mod
USE matrix_mod
USE variable_mod
!
IMPLICIT NONE
!
CHARACTER(LEN=*), INTENT(INOUT) :: line
INTEGER         , INTENT(INOUT) :: length
!
INTEGER, PARAMETER :: MAXP = 2
CHARACTER(LEN=1024), DIMENSION(MAXP) :: cpara
INTEGER            , DIMENSION(MAXP) :: lpara
INTEGER            , DIMENSION(MAXP) :: look
INTEGER :: ianz
INTEGER :: i, j
!
CALL get_params (line, ianz, cpara, lpara, MAXP, length)
!
look(:) = -1
IF(ianz==2) THEN
   locate: DO i=1, ianz
      search: DO j=var_sys+1,var_num
         IF(cpara(i)==var_name(j)) THEN
            look(i) = j
            CYCLE locate
         ENDIF
      ENDDO search
      ier_num = -24
      ier_typ = ER_FORT
      ier_msg(1) = 'Offending variable '//cpara(i)(1:lpara(i))
      RETURN
   ENDDO locate
!
   IF(var_entry(look(1))>0 .AND. var_entry(look(2))>0 ) THEN
      IF(var_field(var_entry(look(1)))%var_shape(1)==var_field(var_entry(look(1)))%var_shape(2) .AND. &
         var_field(var_entry(look(2)))%var_shape(1)==var_field(var_entry(look(2)))%var_shape(2) .AND. &
         var_field(var_entry(look(1)))%var_shape(1)==var_field(var_entry(look(2)))%var_shape(1) )   THEN
         IF(var_field(var_entry(look(1)))%var_shape(1)==1) THEN
            IF(var_field(var_entry(look(1)))%var_value(1,1)/=0.0) THEN 
               var_field(var_entry(look(1)))%var_value(1,1) = &
             1./var_field(var_entry(look(2)))%var_value(1,1)
            ELSE
               ier_num = -45
               ier_typ = ER_FORT
            ENDIF
         ELSEIF(var_field(var_entry(look(1)))%var_shape(1)==2) THEN
            CALL matinv2(var_field(var_entry(look(2)))%var_value,  &
                         var_field(var_entry(look(1)))%var_value)
         ELSEIF(var_field(var_entry(look(1)))%var_shape(1)==3) THEN
            CALL matinv3(var_field(var_entry(look(2)))%var_value,  &
                         var_field(var_entry(look(1)))%var_value)
         ELSEIF(var_field(var_entry(look(1)))%var_shape(1)==4) THEN
            CALL matinv4(var_field(var_entry(look(2)))%var_value,  &
                         var_field(var_entry(look(1)))%var_value)
         ELSE
            ier_num = -47
            ier_typ = ER_FORT
         ENDIF
      ELSE
         ier_num = -44
         ier_typ = ER_FORT
      ENDIF
   ELSE
      ier_num = -46
      ier_typ = ER_FORT
   ENDIF
!
ELSE
   ier_num = -6
   ier_typ = ER_COMM
ENDIF
!
END SUBROUTINE arr_invmat
!
!******************************************************************************
!
SUBROUTINE arr_transpose(line, length)
!
!  Calculate the transpose matrix
!
USE errlist_mod
USE get_params_mod
USE matrix_mod
USE variable_mod
!
IMPLICIT NONE
!
CHARACTER(LEN=*), INTENT(INOUT) :: line
INTEGER         , INTENT(INOUT) :: length
!
INTEGER, PARAMETER :: MAXP = 2
CHARACTER(LEN=1024), DIMENSION(MAXP) :: cpara
INTEGER            , DIMENSION(MAXP) :: lpara
INTEGER            , DIMENSION(MAXP) :: look
INTEGER :: ianz
INTEGER :: i, j
!
CALL get_params (line, ianz, cpara, lpara, MAXP, length)
!
look(:) = -1
IF(ianz==2) THEN
   locate: DO i=1, ianz
      search: DO j=var_sys+1,var_num
         IF(cpara(i)==var_name(j)) THEN
            look(i) = j
            CYCLE locate
         ENDIF
      ENDDO search
      ier_num = -24
      ier_typ = ER_FORT
      ier_msg(1) = 'Offending variable '//cpara(i)(1:lpara(i))
      RETURN
   ENDDO locate
!
   IF(var_entry(look(1))>0 .AND. var_entry(look(2))>0 ) THEN
      IF(var_field(var_entry(look(1)))%var_shape(1)==var_field(var_entry(look(2)))%var_shape(2) .AND. &
         var_field(var_entry(look(1)))%var_shape(2)==var_field(var_entry(look(2)))%var_shape(1) ) THEN
         DO i=1,var_field(var_entry(look(1)))%var_shape(1)
            DO j=1,var_field(var_entry(look(1)))%var_shape(2)
               var_field(var_entry(look(1)))%var_value(i,j) = &
                   var_field(var_entry(look(2)))%var_value(j,i)
            ENDDO
         ENDDO
      ELSE
         ier_num = -44
         ier_typ = ER_FORT
      ENDIF
   ELSE
      ier_num = -46
      ier_typ = ER_FORT
   ENDIF
!
ELSE
   ier_num = -6
   ier_typ = ER_COMM
ENDIF
!
END SUBROUTINE arr_transpose
!
!******************************************************************************
!
END MODULE arrays_mod
