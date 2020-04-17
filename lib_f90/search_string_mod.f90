MODULE search_string_mod
!
CONTAINS
!
!*****7**************************************************************** 
!
INTEGER FUNCTION suche_vor (string, i) 
!-                                                                      
!     searches for the last '(' or '.') in the string                   
!+                                                                      
IMPLICIT none 
!                                                                       
CHARACTER (LEN=*), INTENT(IN   ) :: string 
INTEGER          , INTENT(INOUT) :: i
!
      INTEGER level 
      LOGICAL l_hyp , l_quote
!                                                                       
      l_hyp = .false. 
      l_quote = .false. 
      level = 0 
      DO while (i.gt.0.and.level.ge.0) 
      IF (string (i:i) .eq.')'.or.string (i:i) .eq.']') then 
         level = level + 1 
      ELSEIF (string (i:i) .eq.''''.and..not.l_hyp) then 
         level = level + 1 
         l_hyp = .true. 
      ELSEIF (string (i:i) .eq.''''.and.l_hyp) then 
         level = level - 1 
         l_hyp = .false. 
      ELSEIF (string (i:i) .eq.'"'.and..not.l_quote) then 
         level = level + 1 
         l_quote = .true. 
      ELSEIF (string (i:i) .eq.'"'.and.l_quote) then 
         level = level - 1 
         l_quote = .false. 
      ELSEIF (string (i:i) .eq.'('.or.string (i:i) .eq.'[') then 
         level = level - 1 
      ELSEIF (level.eq.0.and.string (i:i) .eq.',') then 
         level = - 1 
      ELSEIF (i.gt.4) then 
         IF (string (i - 4:i) .eq.'.not.'.or.string (i - 4:i) .eq.'.and.') then
            level = - 1 
         ELSEIF (i.gt.3) then 
            IF (string (i - 3:i) .eq.'.lt.'.or. &
                string (i - 3:i) .eq.'.le.'.or. &
                string (i - 3:i) .eq.'.gt.'.or. &
                string (i - 3:i) .eq.'.ge.'.or. &
                string (i - 3:i) .eq.'.eq.'.or. &
                string (i - 3:i) .eq.'.ne.'.or. &
                string (i - 3:i) .eq.'.or.'.or. &
                string (i - 4:i) .eq.'.eqv.'.or. &
                string (i - 4:i) .eq.'.xor.'     ) then        
               level = - 1 
            ENDIF 
         ENDIF 
      ELSEIF (i.gt.3) then 
         IF (string (i - 3:i) .eq.'.lt.'.or. &
             string (i - 3:i) .eq.'.le.'.or. &
             string (i - 3:i) .eq.'.gt.'.or. &
             string (i - 3:i) .eq.'.ge.'.or. &
             string (i - 3:i) .eq.'.eq.'.or. &
             string (i - 3:i) .eq.'.ne.'.or. &
             string (i - 3:i) .eq.'.or.'      ) THEN
            level = - 1 
         ENDIF 
      ENDIF 
      i = i - 1 
      ENDDO 
      suche_vor = i + 2 
END FUNCTION suche_vor                        
!*****7**************************************************************** 
      INTEGER FUNCTION suche_nach (string, laenge) 
!-                                                                      
!     searches for the first '(' or '.') in the string                  
!+                                                                      
      IMPLICIT none 
!                                                                       
CHARACTER (LEN=*), INTENT(IN   ) :: string 
INTEGER          , INTENT(IN   ) :: laenge
!
      INTEGER i, level 
!                                                                       
      i = 1 
      level = 0 
      DO while (i.le.laenge.and.level.ge.0) 
      IF (string (i:i) .eq.')'.or.string (i:i) .eq.']') then 
         level = level - 1 
      ELSEIF (string (i:i) .eq.'('.or.string (i:i) .eq.'[') then 
         level = level + 1 
      ELSEIF (level.eq.0.and.string (i:i) .eq.',') then 
         level = - 1 
      ELSEIF (i.lt.laenge-4) then 
         IF (string (i:i + 4) .eq.'.not.'.or.string (i:i + 4) .eq.'.and.') then
            level = - 1 
         ELSEIF (i.lt.laenge-3) then 
            IF (string (i:i + 3) .eq.'.lt.'.or.   &
                string (i:i + 3) .eq.'.le.'.or.   &
                string (i:i + 3) .eq.'.gt.'.or.   &
                string (i:i + 3) .eq.'.ge.'.or.   &
                string (i:i + 3) .eq.'.eq.'.or.   &
                string (i:i + 3) .eq.'.ne.'.or.   &
                string (i:i + 3) .eq.'.or.'.or.   &
                string (i:i + 4) .eq.'.eqv.'.or.   &
                string (i:i + 4) .eq.'.xor.') then        
               level = - 1 
            ENDIF 
         ENDIF 
      ELSEIF (i.lt.laenge-3) then 
         IF (string (i:i + 3) .eq.'.lt.'.or.   &
             string (i:i + 3) .eq.'.le.'.or.   &
             string (i:i + 3) .eq.'.gt.'.or.   &
             string (i:i + 3) .eq.'.ge.'.or.   &
             string (i:i + 3) .eq.'.eq.'.or.   &
             string (i:i + 3) .eq.'.ne.'.or.   &
             string (i:i + 3) .eq.'.or.') THEN
            level = - 1 
         ENDIF 
      ENDIF 
      i = i + 1 
      ENDDO 
      suche_nach = i - 2 
      END FUNCTION suche_nach                       
!
!*****7**************************************************************** 
!
INTEGER FUNCTION suche_vor2 (string, i) 
!-                                                                      
!     Searches for the last '(' or '*' or '+' or '-'                    
!     in the string                                                     
!+                                                                      
      IMPLICIT none 
!                                                                       
CHARACTER (LEN=*), INTENT(IN   ) :: string 
INTEGER          , INTENT(INOUT) :: i
!
      LOGICAL lcont 
      INTEGER j 
!                                                                       
      suche_vor2 = 0 
      lcont = .true. 
      DO while (i.gt.0.and.lcont) 
      IF (string (i:i) .eq.'('.or.string (i:i) .eq.'*'.or.string (i:i)  .eq.'/') then
!     &                 string(i:i).eq.'/'               .or.           
!     &                (string(i:i).eq.'-' .and. i.eq.1) .or.           
!     &                (string(i:i).eq.'+' .and. i.eq.1)                
         lcont = .false. 
      ELSEIF (i.gt.1) then 
         IF ( (string (i:i) .eq.'-'.or.string (i:i) .eq.'+') .and. &
              (string (i - 1:i - 1) .ne.'E'.and.string (i - 1:i - 1) .ne.'e' .AND. &
               string (i - 1:i - 1) .ne.'D'.and.string (i - 1:i - 1) .ne.'d' &
              ) ) then                                                         
            lcont = .false. 
         ENDIF 
      ENDIF 
      IF (lcont) then 
         i = i - 1 
      ENDIF 
      ENDDO 
!DBG                                                                    
!DBG      if(j.eq.-1) then                                              
!DBG      do while(i.gt.0       .and. string(i:i).ne.'(' .and.          
!DBG     &   string(i:i).ne.'*'   .and.                                 
!DBG     &   string(i:i).ne.'/'   .and.                                 
!DBG     &   (string(i:i).ne.'-'   .or. (i.gt.1             .and.       
!DBG     &   string(i:i).eq.'-'   .and.                                 
!DBG     &   (string(i-1:i-1).eq.'E'.or.string(i-1:i-1).eq.'e') ) .or.  
!DBG     &   (string(i:i).eq.'-'   .and. i.eq.1)               ) .and.  
!DBG     &   (string(i:i).ne.'+'   .or. (i.gt.1             .and.       
!DBG     &   string(i:i).eq.'+'   .and.                                 
!DBG     &   (string(i-1:i-1).eq.'E'.or.string(i-1:i-1).eq.'e') ) .or.  
!DBG     &   (string(i:i).eq.'+'   .and. i.eq.1                  ) ))   
!DBG        i=i-1                                                       
!DBG      write(*,*) i                                                  
!DBG      ENDDO                                                         
!DBG      endif                                                         
      suche_vor2 = i + 1 
      IF (i.gt.0) then 
         IF (i.gt.1.and. (string (i:i) .eq.'-'.or.string (i:i) .eq.'+') ) then                                                         
            IF (i.eq.1) then 
               suche_vor2 = i 
            ELSE 
               j = i - 1 
               DO while (j.gt.1.and.string (j:j) .eq.' ') 
               j = j - 1 
               ENDDO 
               IF (string (j:j) .eq.'*'.or. &
                   string (j:j) .eq.'/'.or. &
                   string (j:j) .eq.'+'.or. &
                   string (j:j) .eq.'-'.or. &
                   string (j:j) .eq.'(') then
                  suche_vor2 = i 
               ENDIF 
            ENDIF 
         ENDIF 
      ENDIF 
      END FUNCTION suche_vor2                       
!
!*****7**************************************************************** 
!
INTEGER FUNCTION suche_nach2 (string, laenge) 
!-                                                                      
!     Searches for the first '(' or '*' or '/' or '+' or '-'            
!     in the string                                                     
!       im string                                                       
!+                                                                      
      IMPLICIT none 
!                                                                       
CHARACTER (LEN=*), INTENT(IN   ) :: string 
INTEGER          , INTENT(IN   ) :: laenge
!
      LOGICAL lcont 
      INTEGER i 
!                                                                       
      suche_nach2 = 0 
      i = 1 
      IF (laenge.eq.0) return 
      DO while (i.le.laenge.and.string (i:i) .eq.' ') 
      i = i + 1 
      ENDDO 
      IF (string (i:i) .eq.'-'.or.string (i:i) .eq.'+') then 
         i = i + 1 
      ENDIF 
      IF (laenge.gt.1) then 
         lcont = .true. 
         DO while (i.le.laenge.and.lcont) 
            IF (string (i:i) .eq.')'.or.string (i:i) .eq.'*'.or. &
                string (i:i) .eq.'/'.or.                         &
               (string (i:i) .eq.'-'.and.i.eq.1) .or.            &
               (string (i:i) .eq.'+'.and.i.eq.1) ) then                                
            lcont = .false. 
            ELSEIF (i.gt.1) then 
               IF ( (string (i:i) .eq.'-'.or.          &
                     string (i:i) .eq.'+'     )  .and. &
                    (                                  &
                    (string (i - 1:i - 1) .ne.'E'.and. &
                     string (i - 1:i - 1) .ne.'e'.AND. &
                     string (i - 1:i - 1) .ne.'D'.and. &
                     string (i - 1:i - 1) .ne.'d')     &
                    )                                  &
                    ) then                                             
               lcont = .false. 
               ENDIF 
            ENDIF 
            IF (lcont) then 
               i = i + 1 
            ENDIF 
         ENDDO 
!DBG                                                                    
         IF (i.eq. - 1111) then 
            DO while (i.le.laenge.and.                                &
               string (i:i) .ne.')' .and.  string (i:i) .ne.'*' .and. &
               string (i:i) .ne.'/' .and.                             &
              (string (i:i) .ne.'-' .or.                              &
               (i.gt.1             .and.                              &
                string (i:i) .eq.'-'.and.                             &
               (string (i - 1:i - 1) .eq.'E'.or.                      &
                string (i - 1:i - 1) .eq.'e' .OR.                     &
                string (i - 1:i - 1) .eq.'D'.or.                      &
                string (i - 1:i - 1) .eq.'d')                         &
               ) ) .and.                                              &
              (string (i:i) .ne.'+'.or.                               &
               (i.gt.1             .and.                              &
                string (i:i) .eq.'+' .and.                            &
                (string (i - 1:i - 1) .eq.'E'.or.                     &
                 string (i - 1: i - 1) .eq.'e' .OR.                   &
                 string (i - 1:i - 1) .eq.'D'.OR.                     &
                 string (i - 1: i - 1) .eq.'d'                        &
                ) ) ) )                                       
            i = i + 1 
            ENDDO 
         ENDIF 
         suche_nach2 = i - 1 
      ELSE 
         suche_nach2 = i 
      ENDIF 
      END FUNCTION suche_nach2                      
!
!*****7**************************************************************** 
!
INTEGER FUNCTION suche_vor_hoch (string, i) 
!-                                                                      
!     searches for the last '(' or '.') in the string                   
!+                                                                      
      IMPLICIT none 
!                                                                       
CHARACTER (LEN=*), INTENT(IN)    :: string 
INTEGER          , INTENT(INOUT) :: i
!
      INTEGER level 
      LOGICAL l_quote
!                                                                       
      l_quote = .false. 
      level = 0 
      DO while (i.gt.0.and.level.ge.0) 
      IF (string (i:i) .eq.'"'.and..not.l_quote) then 
         level = level + 1 
         l_quote = .true. 
      ELSEIF (string (i:i) .eq.'"'.and.l_quote) then 
         level = level - 1 
         l_quote = .false. 
      ENDIF
      i = i - 1 
      ENDDO 
      suche_vor_hoch = i + 2 
END FUNCTION suche_vor_hoch
!
!*****7**************************************************************** 
!
INTEGER FUNCTION suche_nach_hoch (string, laenge) 
!-                                                                      
!     searches for the last '(' or '.') in the string                   
!+                                                                      
      IMPLICIT none 
!                                                                       
CHARACTER (LEN=*), INTENT(IN)    :: string 
INTEGER, INTENT(IN) :: laenge
!
      INTEGER i, level 
      LOGICAL l_quote
!                                                                       
      i = 1
      l_quote = .false. 
      level = 0 
      DO while (i.le.laenge.and.level.ge.0) 
      IF (string (i:i) .eq.'"'.and..not.l_quote) then 
         level = level + 1 
         l_quote = .true. 
      ELSEIF (string (i:i) .eq.'"'.and.l_quote) then 
         level = level - 1 
         l_quote = .false. 
      ENDIF
      i = i + 1 
      ENDDO 
      suche_nach_hoch = i - 2 
END FUNCTION suche_nach_hoch
!
!****7***************************************************************** 
!
END MODULE search_string_mod
