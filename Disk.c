
#include "Disk.h"
#include "ComTypes.h"

static void	setApproximateReferenceLength(disk_handler_t *handler) {
	/* seek to the end */
	fseek(handler->reffile, 0L, SEEK_END);
	handler->refapplen = ftell(handler->reffile);

	/* seek back to the beginning */
	fseek(handler->reffile, 0L, SEEK_SET);
}

static void	setApproximateQueryLength(disk_handler_t *handler) {
	/* seek to the end */
	fseek(handler->qryfile, 0L, SEEK_END);
	handler->qryapplen = ftell(handler->qryfile);

	/* seek back to the beginning */
	fseek(handler->qryfile, 0L, SEEK_SET);
}

static void	setApproximateTargetLength(disk_handler_t *handler) {
	/* seek to the end */
	fseek(handler->reftargetfile, 0L, SEEK_END);
	handler->tarapplen = ftell(handler->reftargetfile);

	/* seek back to the beginning */
	fseek(handler->reftargetfile, 0L, SEEK_SET);
}

int diskInit(disk_handler_t *handler, char *reffile, 
			char *qryfile, 
			int refmaxlen, 
			char* reftargetfile, 
			char* align, 
			char* wordscdec,
			char* destinationDirectory) {
	assert(handler);
	assert(reffile);
	assert(qryfile);
	assert(reftargetfile);
	assert(align);
	assert(!handler->reffile);
	assert(!handler->qryfile);
	assert(!handler->reftargetfile);
	assert(!handler->align);
	assert(!handler->wordpossibility);
	
	handler->reffile = fopen(reffile, "r");
	if (!handler->reffile) {
		fprintf(stderr, "Can not open reference file \"%s\"\n", reffile);
		return 0;
	}

	handler->qryfile = fopen(qryfile, "r");
	if (!handler->qryfile) {
		fprintf(stderr, "Can not open query file \"%s\"\n", qryfile);
		fclose(handler->reffile);
		handler->reffile = NULL;
		return 0;
	}

	handler->reftargetfile = fopen(reftargetfile, "r");
	if (!handler->reftargetfile) {
		fprintf(stderr, "Can not open reference file \"%s\"\n", reftargetfile);
		return 0;
	}

	handler->align = fopen(align, "r");
	if (!handler->align) {
		fprintf(stderr, "Can not open reference file \"%s\"\n", align);
		return 0;
	}

	handler->wordpossibility = wordscdec;
	handler->destinationDirectory = destinationDirectory;
	
	setApproximateReferenceLength(handler);
	setApproximateQueryLength(handler);
	setApproximateTargetLength(handler);

	handler->refmaxlen	= handler->refapplen;//refmaxlen > handler->refapplen ? handler->refapplen : refmaxlen;
	handler->refisfirst	= 1;
	//handler->seq = (void *)kseq_init(handler->qryfile);
	handler->con = 0;
	return 1;
}


void diskFinalize(disk_handler_t *handler) {

	assert(handler);
	assert(handler->reffile);
	assert(handler->qryfile);
	assert(handler->reftargetfile);	
	assert(handler->align);
	
	fclose(handler->reffile);
	fclose(handler->qryfile);
	fclose(handler->reftargetfile);
	fclose(handler->align);

	handler->reffile = NULL;
	handler->qryfile = NULL;
	handler->reftargetfile = NULL;
	handler->align = NULL;
}


void diskSetRefMaxSize(disk_handler_t *handler, int refmaxlen) {
	assert(handler);
	assert(refmaxlen > 0);

	handler->refmaxlen = refmaxlen;
}
