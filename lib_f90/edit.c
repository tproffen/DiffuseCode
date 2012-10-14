/******************************************************/
/* New version linking to readline library            */
/******************************************************/

#ifdef VMS
#define cread_ CREAD
#define cinit_ CINIT
#endif

#ifdef hpux
#define cread_ cread
#define cinit_ cinit
#endif

#ifdef __linux__
#define _XOPEN_SOURCE
#include <sys/types.h>  /* for kill(2) */
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*----------------------------------------------------*/
void cinit_()

{
#ifdef READLINE
	printf (" Command line editing enabled ..\n\n");
#else
	printf (" Command line editing disabled ..\n\n");
#endif
}

/*----------------------------------------------------*/
void cread_( char *prom, int *lp,
             char *inp, int *li )
{
        char *cprom;
	char *cinp;
        size_t len;


#ifdef READLINE
        if( *lp <= 0 || *li <= 0 ) {  /* sanity check */
            *li = -1;
             return;
        }
        cprom = malloc( (size_t)*lp+1 );
        if( !cprom ) {
            *li = -1;
             return;
        }
        strncpy( cprom, (char *)prom, (size_t)*lp );
        cprom[*lp] = '\0';
        cinp = readline(cprom);

        if(cinp && *cinp) {
           len = strlen(cinp);
           *li = len < *li ? len : *li;
           memcpy( inp, cinp, (size_t)*li );
           add_history(cinp);
           free (cinp);
	}

        free (cprom);
#else
        *li = -1;
#endif
        return;
}

