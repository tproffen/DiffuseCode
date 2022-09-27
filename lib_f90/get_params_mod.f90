MODULE get_params_mod
!
private
!
public get_params
public get_params_blank
public del_params
!
CONTAINS
!
!*****7***************************************************************  
!
SUBROUTINE get_params (string, ianz, cpara, lpara, nwerte, laenge) 
!-                                                                      
!     Reads parameters that have to be separated by ",". Each           
!       expression between "," is evaluated to give the value of        
!       the parameter. If the routine is called with a negative         
!     length, leading SPACE will NOT be removed !                       
!+                                                                      
      USE blanks_mod
      USE errlist_mod 
      IMPLICIT none 
!                                                                       
CHARACTER(LEN=*)                      , INTENT(IN)    :: string 
INTEGER                               , INTENT(OUT)   :: ianz
INTEGER                               , INTENT(IN)    :: nwerte 
CHARACTER(LEN=*   ), DIMENSION(nwerte), INTENT(OUT)   :: cpara
INTEGER            , DIMENSION(nwerte), INTENT(OUT)   :: lpara
INTEGER                               , INTENT(INOUT) :: laenge 
!                                                                       
      INTEGER i, lll 
      INTEGER ipos
      INTEGER level 
      LOGICAL quote 
      LOGICAL search 
      LOGICAL rmblk 
!                                                                       
      ipos    = 1 
      ianz    = 0 
      ier_num = 0 
      ier_typ = ER_NONE 
      DO i = 1, nwerte 
         cpara (i) = ' ' 
         lpara (i) = 1 
      ENDDO
!
      IF (laenge.eq.0) return 
!                                                                       
      IF (laenge.lt.0) then 
         laenge = - laenge 
         rmblk = .false. 
      ELSE 
         rmblk = .true. 
      ENDIF 
!
      IF ( LEN_TRIM(string).eq.0) THEN 
        laenge = 0
        RETURN
      ENDIF
!                                                                       
      search = .true. 
      level = 0 
      ianz = 1 
      lpara (ianz) = 0 
      DO i = 1, laenge 
      quote = string (i:i) .eq.'"'.or.string (i:i) .eq.'''' 
      IF (quote) then 
         search = .not.search 
      ENDIF 
      IF (search) then 
         IF (level.eq.0.and.string (i:i) .eq.',') then 
            IF (ianz.lt.nwerte) then 
               ianz = ianz + 1 
               lpara (ianz) = 0 
            ELSE 
               ier_num = - 17 
               ier_typ = ER_COMM 
               RETURN 
            ENDIF 
         ELSE 
            lpara (ianz) = lpara (ianz) + 1 
            cpara (ianz) (lpara (ianz) :lpara (ianz) ) = string (i:i) 
            IF (string (i:i) .eq.'('.or.string (i:i) .eq.'[') then 
               level = level + 1 
            ENDIF 
            IF (string (i:i) .eq.')'.or.string (i:i) .eq.']') then 
               level = level - 1 
            ENDIF 
         ENDIF 
      ELSE 
         lpara (ianz) = lpara (ianz) + 1 
         cpara (ianz) (lpara (ianz) :lpara (ianz) ) = string (i:i) 
         IF (string (i:i) .eq.'('.or.string (i:i) .eq.'[') then 
            level = level + 1 
         ENDIF 
         IF (string (i:i) .eq.')'.or.string (i:i) .eq.']') then 
            level = level - 1 
         ENDIF 
      ENDIF 
      ENDDO 
!
! RBN JULY 2018 ALWAYS remove leading and trailing banks !!!
!                                                                       
!     remove leading blanks if length was given positive                
!                                                                       
      IF (rmblk) then 
         DO i = 1, ianz 
         lll = lpara (i) 
         CALL rem_bl (cpara (i), lll) 
         lpara (i) = lll 
         ENDDO 
      ENDIF 
      DO i = 1, ianz 
         lll = lpara (i) 
         CALL rem_leading_bl (cpara (i), lll) 
         lpara (i) = lll 
      ENDDO 
!
!     remove trailing blanks if length was given positive
!
!     IF(rmblk) THEN
         DO i=1,ianz
            lpara(i) = LEN_TRIM(cpara(i))
         ENDDO
!     ENDIF
!                                                                       
!      Check parameter length to catch things like ...,,...             
!      ianz is decremented until ALL get_parameter calls do             
!      an error check...                                                
!                                                                       
      DO i = 1, ianz 
      IF (lpara (i) .eq.0) then 
         ier_num = - 2 
         ier_typ = ER_COMM 
         ianz = max (0, ianz - 1) 
         RETURN 
      ENDIF 
      ENDDO 
!                                                                       
!      Check if we are still in search mode, ie. missing ""s            
!      ianz is decremented until ALL get_parameter calls do             
!      an error check...                                                
!                                                                       
      IF (.not.search) then 
         ier_num = - 30 
         ier_typ = ER_FORT 
         ianz = max (0, ianz - 1) 
      ENDIF 
!                                                                       
      IF (level.ne.0) then 
         ier_num = - 9 
         ier_typ = ER_FORT 
         ianz = max (0, ianz - 1) 
      ENDIF 
!                                                                       
      RETURN 
!CCC                                                                    
!     OLD CODE, KEPT UNTIL ERROR CHECK IS COMPLETED!                    
!                                                                       
!DBG      iko=index(string(ipos:laenge),',')                            
!DBG      if(iko.ne.0) then                                             
!DBG        lll = laenge - ipos + 1                                     
!DBG        ikk=suche_nach(string(ipos:laenge),lll)+1                   
!DBG        if(ikk.eq.laenge) then                                      
!DBG          ikk=ikk+1                                                 
!DBG        endif                                                       
!DBG        do while(ikk.gt.1 .and. ianz.lt.nwerte-1)                   
!DBG          ianz       =ianz+1                                        
!DBG          cpara(ianz)=string(ipos:ipos+ikk-2)                       
!DBG          lpara(ianz)= (ipos+ikk-2) - ipos + 1                      
!DBG          ipos       =ipos+ikk                                      
!DBG          if(ipos.le.laenge) then                                   
!DBG            iko=INDEX(string(ipos:laenge),',')                      
!DBG            if(iko.ne.0) then                                       
!DBG              lll = laenge - ipos + 1                               
!DBG              ikk=suche_nach(string(ipos:laenge),lll)+1             
!DBG              if(ikk+ipos.gt.laenge) goto 10                        
!DBG            else                                                    
!DBG              ikk=0                                                 
!DBG            endif                                                   
!DBG          else                                                      
!DBG            goto 10                                                 
!DBG          endif                                                     
!DBG        ENDDO                                                       
!DBG      endif                                                         
!DBG10      continue                                                    
!DBG      if(ipos.le.laenge .and. string(ipos:laenge).ne. ' ') then     
!DBG        ianz       =ianz+1                                          
!DBG        cpara(ianz)=string(ipos:laenge)                             
!DBG        lpara(ianz)=laenge - ipos + 1                               
!DBG      endif                                                         
!DBG      ier_num = 0                                                   
!DBG      ier_typ = ER_NONE                                             
!DBGc                                                                   
!DBGc     remove leading blanks if length was given positive            
!DBGc                                                                   
!DBG      if (rmblk) then                                               
!DBG        do i=1,ianz                                                 
!DBG          lll      = lpara(i)                                       
!DBG          call rem_bl(cpara(i),lll)                                 
!DBG          lpara(i) = lll                                            
!DBG        ENDDO                                                       
!DBG      endif                                                         
      END SUBROUTINE get_params                     
!
!*****7***************************************************************  
!
SUBROUTINE get_params_blank (string, ianz, cpara, lpara, nwerte, laenge) 
!-                                                                      
!     Reads parameters that have to be separated by  a blank " ".
!+                                                                      
      USE errlist_mod 
      USE charact_mod 
      IMPLICIT none 
!                                                                       
!                                                                       
CHARACTER (LEN=* )                      ,INTENT(IN)  :: string 
INTEGER                                 ,INTENT(OUT) :: ianz
INTEGER                                 ,INTENT(IN)  :: nwerte 
CHARACTER(LEN=*   ),DIMENSION(1:nwerte), INTENT(OUT) :: cpara
INTEGER            ,DIMENSION(1:nwerte), INTENT(OUT) :: lpara
INTEGER                                 ,INTENT(IN ) :: laenge
!
      INTEGER :: i,j,k
      LOGICAL :: no_par
!
      cpara(:) = ' '
      j = 0
      k = 0
      no_par = .true.
main: DO i=1, laenge
         IF(no_par) THEN
            IF(string(i:i)==' ' .or. string==TAB) THEN  ! No param, cycle
               CYCLE main
            ELSE
               no_par = .false.                         !parameter ==> on
               j = j + 1                                !increment param no
               k = 1
               cpara(j)(k:k) = string(i:i)
               lpara(j)      = 1
            ENDIF
         ELSE
            IF(string(i:i)==' ' .or. string==TAB) THEN  ! No param, cycle
               no_par = .true.
            ELSE
               k = k + 1
               cpara(j)(k:k) = string(i:i)
               lpara(j)      = lpara(j) + 1
            ENDIF
         ENDIF
      ENDDO main
      ianz = j
!
END SUBROUTINE get_params_blank
!
!*****7***************************************************************  
!
SUBROUTINE del_params (ndel, ianz, cpara, lpara, nwerte) 
!-                                                                      
!     Deletes the first <ndel> parameters from the list and             
!     moves the rest to the beginning of the arrays.                    
!+                                                                      
USE errlist_mod 
!
IMPLICIT none 
!                                                                             
INTEGER                                 ,INTENT(IN )   :: ndel
INTEGER                                 ,INTENT(INOUT) :: ianz
INTEGER                                 ,INTENT(IN)    :: nwerte 
CHARACTER(LEN=*   ),DIMENSION(1:nwerte), INTENT(INOUT) :: cpara
INTEGER            ,DIMENSION(1:nwerte), INTENT(INOUT) :: lpara
!                                                                       
INTEGER i 
!                                                                       
IF (ndel.lt.ianz) then 
   DO i = ndel + 1, ianz 
     cpara (i - ndel) = cpara (i) 
     lpara (i - ndel) = lpara (i) 
   ENDDO 
   ianz = ianz - ndel 
ELSEIF(ndel==ianz) THEN
   cpara(:) = ' '
   lpara(:) = 1
   ianz = 0
ENDIF 
!                                                                       
END SUBROUTINE del_params                     
!
!*****7***************************************************************  
!
END MODULE get_params_mod
