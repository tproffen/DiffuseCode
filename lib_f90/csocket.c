/********************************************************************/
/*                                                                 */
/* Interface to socket communication routines from FORTRAN          */
/*                                                                  */
/********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <string.h>

#ifdef WIN32
#include <winsock2.h>
#else
#include <sys/socket.h>
#include <netinet/in.h>
#include <netdb.h>
#endif

/********************************************************************/

int socket_init_ (int *sock, int *port, unsigned int *local)
{
	struct   sockaddr_in sin;
 
	/* create the socket */

	if((*sock=socket(AF_INET, SOCK_STREAM, 0))==-1) {
	  perror("socket create failed");
	}
	memset(&sin, 0, sizeof(sin));
	sin.sin_family = AF_INET;
	if (*local == 1) {
  	  sin.sin_addr.s_addr = inet_addr("127.0.0.1");
          printf(" Running in local mode (127.0.0.1) ..\n");
        } else {
	  sin.sin_addr.s_addr = INADDR_ANY;
          printf(" Running in network mode ..\n");
        }
	printf(" Listening to port %d ..\n",*port);
	sin.sin_port = htons(*port);

	/* socket must map into filespace */

	if (bind(*sock, (struct sockaddr *) &sin, sizeof(sin)) == -1) {
	  perror("socket bind problem");
        }

	/* prepare for listening */
	
	if (listen(*sock, -1) == -1) perror("socket listen problem");
	return (0);
}

/********************************************************************/

int socket_accept_ (int *sock, int *cid, unsigned char *host,
                     int *lhost)
{

	char     hostname[100];
	char	 iall[16];
	char	 ireq[16];
	int	 addrlen, test;
	struct   sockaddr_in pin,pout;
	struct   hostent *hp;

	memcpy(hostname,host,(size_t)*lhost);
	hostname[*lhost]='\0';
	if ((hp = gethostbyname(hostname)) == 0) {perror("gethostbyname");}
	pout.sin_addr.s_addr = ((struct in_addr *)(hp->h_addr))->s_addr;
        memcpy(iall,inet_ntoa(pout.sin_addr),16);
        printf(" Allowing connections from %s ..\n",iall);

	/* wait for a client to talk to us */

	do {
          addrlen = sizeof(pin); 
	  if ((*cid = accept(*sock,(struct sockaddr *)&pin,&addrlen)) == -1) {
	    perror("socket accept problem");
  	  }

          memcpy(ireq,inet_ntoa(pin.sin_addr),16);
	  test=((strcmp(ireq,iall)==0) || (strcmp(ireq,"127.0.0.1")==0));
	  if (test==0) {
            printf(" REJECTED connection from %s:%d ..\n",
                     inet_ntoa(pin.sin_addr),ntohs(pin.sin_port));
	    close(*cid);
	  }
	} while (test==0);

        printf(" Accepted connection from %s:%d ..\n",
                 inet_ntoa(pin.sin_addr),ntohs(pin.sin_port));
	return(0);
}

/********************************************************************/

int socket_close_ (int *sock)

{
	close(*sock);
	return(0);
} 

/********************************************************************/

int socket_get_ (int *sock, unsigned char *str, int *is)

{
	int  slen;
	char cstr[256];

	slen=recv(*sock, cstr, sizeof(cstr), 0);
	if (slen == -1) {perror("Recv error");}

	*is=slen-1;
	memcpy(str,cstr,(size_t)*is);

	return(0);
}
/********************************************************************/

int socket_connect_ (int *sock, unsigned char *host, int *lhost,
                      int *port)
{
        char   hostname[100];
	struct sockaddr_in sin;
	struct sockaddr_in pin;
	struct hostent *hp;

	*sock=0;

        memcpy(hostname,host,(size_t)*lhost);
        hostname[*lhost]='\0';

	/* go find out about the desired host machine */

	if ((hp = gethostbyname(hostname)) == 0) {perror("gethostbyname");}

	/* fill in the socket structure with host information */

	memset(&pin, 0, sizeof(pin));
	pin.sin_family = AF_INET;
	pin.sin_addr.s_addr = ((struct in_addr *)(hp->h_addr))->s_addr;
	pin.sin_port = htons(*port);

	/* grab an Internet domain socket */

	if ((*sock=socket(AF_INET, SOCK_STREAM, 0)) == -1) {perror("socket");}

	/* connect to PORT on HOST */

	printf (" Connecting to %s:%d ..\n",inet_ntoa(pin.sin_addr),
                ntohs(pin.sin_port));
	if (connect(*sock,(struct sockaddr *)  &pin, sizeof(pin)) == -1) {
		perror("connect");
	}
	return(0);
}

/********************************************************************/

int socket_send_ (int *sock, unsigned char *cmd, int *ic)

{
        char cstr[256];
        int  il;

	il=*ic + 1;
        memcpy(cstr,cmd,(size_t) il);
        cstr[il-1]='\n';
        cstr[il]='\0';
	if (send(*sock,cstr,(size_t) il,0) == -1) {perror("send problem");}
}

