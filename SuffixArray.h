#ifndef SUFFIX_ARRAY_H
#define SUFFIX_ARRAY_H

#include "ComTypes.h"


void suffixArrayConstruct(ref_t *ref, int yoho, int* real);

void *suffixArraySearchInit(int refbufsize, int qrysbufsize, int size, int sasize);
void suffixArraySearchFinalize(void * handler);
void suffixArraySearch(void *handler,
	ref_set_t *refset,
	qry_set_t *qryset, 
	int minmatch,
	timing_t *timing, 
	int twoway);

int suffixArrayGetEquivalentMaxRefLen(int bufsize, int fingerlen);
int suffixArrayGetRequiredRefBufferSize(int reflen, int fingerlen);
int suffixArrayGetRequiredQrySetResultBufferSize(qry_set_t *qryset, int minmatchlen);
int recursion_lcp(int L, int R, int *lcp, int *lcpright, int *lcpleft);
void suffixArraySearchFinalize_One(void *handler);

#endif /* SUFFIX_ARRAY_H */
