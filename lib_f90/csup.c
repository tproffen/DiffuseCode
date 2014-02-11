/*******************************************************************/
/*                                                                 */
/* Various helper routines in C                                    */
/*                                                                 */
/*******************************************************************/

#ifdef WIN32
#include "win32-glob.h"
#else
#include <glob.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

glob_t globbuf;

/*******************************************************************/
/* Routines to glob files                                          */
/*   integer ifiles(mask,lm):        Find files matching "mask"    */
/*   subroutine getfile(fname,il,i): Get filename number "i"       */
/*   subroutine freefiles()        : Free space alloc. by ifiles   */
/*******************************************************************/
int ifiles_ (unsigned char *mask, int *l)
{
	char *cmask;

	cmask=malloc((int)*l+1);
	strncpy(cmask,(char *)mask,(int)*l);
	cmask[*l]='\0';

#ifdef WIN32
glob(cmask, 0, NULL, &globbuf);
#else
globbuf.gl_offs = 1;
glob(cmask, GLOB_DOOFFS, NULL, &globbuf);
#endif
	return globbuf.gl_pathc; 
}

void getfile_ (unsigned char *file, int *l, int *i)
{
	int	len;

	len=strlen(globbuf.gl_pathv[*i]);
	memcpy(file,globbuf.gl_pathv[*i],len);
}

void freefiles_()
{
	globfree(&globbuf);
}
