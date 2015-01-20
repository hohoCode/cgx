#ifndef DISK_H
#define DISK_H

#include "ComTypes.h"

typedef struct disk_handler_s {

	/* reference related setup */
	FILE *reffile;
	int refmaxlen;
	int refisfirst;
	long refapplen; /* the approximate length of the whole reference */

	/* query related setup */
	FILE *qryfile;
	void *seq;
	int con;
	int qryapplen; /* the approximate length of a query (used to make our 
				   life easier when allocating buffers)*/

	/* reference target */
	FILE *reftargetfile;
	long tarapplen; /* the approximate length of the whole reference */
	
	FILE *align;
	char* wordpossibility;
	char* destinationDirectory;
} disk_handler_t;


#define DISK_INIT_HANDLER {0, 0, 0}

int diskQrySetLoad(disk_handler_t *handler,
	qry_set_t *qryset,
	int minmatch);

void diskQrySetFreeNamesArray(qry_set_t *qryset);

int diskRefLoad(disk_handler_t *handler,
	ref_t *ref,
	ref_t *prevref, char* temp);

void diskSetRefMaxSize(disk_handler_t *handler, 
	int refmaxlen);

int diskInit(disk_handler_t *handler, 
	char *reffile, 
	char *qryfile,
	int refmaxlen,  
	char* reftargetfile, 
	char* align,
    char* wordpossibility,
    char* destinationDirectory);

void diskFinalize(disk_handler_t *handler);

#endif /* DISK_H*/
