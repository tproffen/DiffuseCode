#include <stdlib.h>
#include <unistd.h>
#include <stdio.h>
#include <string.h>
#include <readline/readline.h>
#include <readline/history.h>
#
/* -------------------------------------------------------------------------- */
FCreadline(int len, char *myline, char prompt[]){
/*
@(#)FCreadline.sh  return line from readline(3c) to Fortran. John S. Urban, 20100323

   Simple procedure that uses readline in "normal" (i.e. non-callback) mode. 

   len    -- number of characters in argument "myline"
   myline -- Fortran CHARACTER variable to recieve the line read by readline(3c) 
   prompt -- prompt string to preceed read

*/
   char *line;                            /* readline(3c) will return the read line to this pointer */
   int i;                                 /* counter for padding returned line with spaces */

   using_history();  
   line=readline(prompt);                 /* use readline(3c) to read a line of input in edit mode */
   if (line == NULL) {
     line = strdup("exit");
   }
   if(strlen(line) >0) {                  /* save non-zero strings only */
      add_history(line);
   }

   strncpy(myline,line,len);         /* copy line returned by readline(3c) to MYLINE up to length of MYLINE */

   for(i=strlen(line);i<len;i++){    /* starting with null, pad with spaces to end */
     myline[i]=' ';
   }

   /* free memory used to return line from readline(3c) */
   free(line);       
}
