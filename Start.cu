#include "ComTypes.h"
#include "Disk.h"
#include "SuffixArray.h"
#include "Timer.h"
#include "ExtractPair.h"
#include "uthash/uthash.h"
#include "PrintResults.h"

using namespace std;
#define REF_LOAD_BATCH 3
#define PRINTCACHELIMIT 1048576
#define PREALLOCATEDVOCSIZE 500000

static void finalizeRefSet(ref_set_t *refset) {
	int i;
	hashtbl *s, *tmp, *s2;
	hash_intchar *s1, *tmp1, *s3;

	for (i = 0; i < refset->count; i++) {
		cudaFreeHost(refset->refs[i].buf);
		cudaFreeHost(refset->refs[i].buf);

		HASH_ITER(hh, refset->refs[i].users, s, tmp) {
			HASH_DEL(refset->refs[i].users, s);
			//free(s->name);
			free(s);
		}
		HASH_ITER(hh, refset->refs[i].intchar, s1, tmp1) {
			HASH_DEL(refset->refs[i].intchar, s1);
			//free(s1->name);
			free(s1);
		}

		HASH_ITER(hh, refset->refs_target[i].users_target, s2, tmp) {
			HASH_DEL(refset->refs_target[i].users_target, s2);
			free(s2);
		}
		HASH_ITER(hh, refset->refs_target[i].intchar_target, s3, tmp1) {
			HASH_DEL(refset->refs_target[i].intchar_target, s3);
			free(s3);
		}
		//free(refset->refs[i].buf);
		//free(refset->refs[i].tokindex);
		//free(refset->refs[i].tokindex_reverse);
	}
	free(refset->refs);
	free(refset->refs_target);
}

void constructQryIndex(qry_set_t * qryset, disk_handler_t * handler, ref_set_t *refset){
	/////Read in queries line by line, and split them into tokens
	/////Record the tokens index, and the starting index inside each queries
	/////And record the number of tokens  

	ref_t * ref = &refset->refs[0];	
	char * line = NULL;
	size_t len = 0;
	int qcount = 0;    
	const char* delim = " ";    
	int count = 0;
	int read = 0;
	int max = -1;
	int percount = 0;

	int* toklen = (int*)malloc(2*sizeof(int)*(handler->qryapplen));
	int* tokindex_qryindex = (int*)malloc(sizeof(int)*(handler->qryapplen));
	if(handler->qryfile == NULL || toklen == NULL)
	{
		printf("not able to open file to construct tokindex\n") ;
		exit(0) ;
	}	
	assert(toklen !=NULL && ref->users != NULL);
	//fprintf(stderr, "Q: ");
	while(read = getline(&line, &len, handler->qryfile)!= -1){        
		qcount++;    
		if (line[len - 1] == '\n'){
			line[len - 1] = '\0';
		}
		char* token = strtok(line, delim);        
		toklen[qcount-1] = count;
		percount = 0;		
		//fprintf(stderr, "Query %d - Count%d\n", qcount-1, count);
		while(token != NULL && !isspace(*token)){            
			percount++;
			struct my_struct *sss = (struct my_struct*)malloc(sizeof(struct my_struct));   
			if(token[strlen(token)-1] == '\n'){
				token[strlen(token)-1] = '\0';
			}
			HASH_FIND_STR(ref->users, token, sss);            
			if(!sss){                                
				/*s = (hashtbl*)malloc(sizeof(hashtbl));                
				  s->id = HASH_COUNT(users) + 1;
				  last = s->id;
				  s->name = (char*)malloc(sizeof(char)*(strlen(token)+1));                
				  strcpy(s->name, token);                
				  HASH_ADD_KEYPTR( hh, users, s->name, strlen(s->name), s);*/
				toklen[handler->qryapplen+count] = -1;
				//fprintf(stderr, "-1 ");
			} else{
				toklen[handler->qryapplen+count] = sss->id;            
				//fprintf(stderr, "%d ", sss->id);
			}		
			tokindex_qryindex[count] = qcount-1;
			//printf("TOKEN: %s\n", token);           
			count++;		  	
			token = strtok(NULL, delim);           
		}
		if (percount > max){
			max = percount;
		}
	}
	fprintf(stderr, "\nMax length of queries is %d\n", max);
	qryset->totaltokens = count;
	qryset->resbufsize = count*sizeof(result_t_two);
	qryset->qryscount = qcount;
	qryset->qrysbufsize = (qryset->totaltokens+qryset->qryscount)*sizeof(int);

	//qryset->qrysbuf = (int*)malloc(qryset->qrysbufsize);
	//qryset->result_two = (result_t_two *)malloc(qryset->resbufsize);
	//qryset->connectoffset = (int*)malloc(qryset->totaltokens*sizeof(int));
	cudaMallocHost((void **)&(qryset->qrysbuf), qryset->qrysbufsize);
	cudaMallocHost((void **)&(qryset->result_two), qryset->resbufsize);
	cudaMallocHost((void **)&(qryset->connectoffset), qryset->totaltokens*sizeof(int));
	cudaMallocHost((void **)&(qryset->tokindex_qryindex), qryset->totaltokens*sizeof(int));

	memcpy(qryset->qrysbuf, toklen, qryset->qryscount*sizeof(int));
	memcpy(qryset->tokindex_qryindex, tokindex_qryindex, qryset->totaltokens*sizeof(int));
	memcpy(qryset->qrysbuf+qryset->qryscount, toklen+handler->qryapplen, qryset->totaltokens*sizeof(int));
	qryset->qrysoffsettok = qryset->qrysbuf;
	free(toklen);	
	free(tokindex_qryindex);
}

void add_intchar(int user_id, char *name, hash_intchar** all) {
	hash_intchar *s = (hash_intchar *)malloc(sizeof(hash_intchar));
	s->id = user_id;
	s->name = (char*)malloc(sizeof(char)*(strlen(name)+1)); 
	memcpy((void*)s->name, name, strlen(name)+1);
	HASH_ADD_INT(*all,id, s);  /* id: name of key field */
}

static void initRefTargetSet(ref_set_t *refset, disk_handler_t *handler, uint	*devicemem) {
	//fprintf(stderr, "\nLoading the target reference\n");
	refset->refs_target = (ref_t_target*)calloc(1, sizeof(ref_t_target));
	assert(refset->refs_target);
	ref_t_target * ref = &refset->refs_target[0];			 
	int* ref_temp = (int*)malloc(sizeof(int)*(handler->tarapplen+ 3));
	int* senidx = (int*)malloc(sizeof(int)*(handler->tarapplen));  ///Number of lines!!
	assert(ref_temp!=NULL);
	assert(senidx != NULL);

	hashtbl *s = NULL;    	
	ref->users_target = NULL;
	ref->intchar_target = NULL;

	int count = 0;    
	int qcount = 0;
	int read = 0;
	const char* delim = " ";    
	char* line = NULL;    
	size_t len = 0;
	int last = -1;
	senidx[0] = 0;
	///Direct Access to the string: from ID to char* pointers
	//ref->vocabulary.push_back(NULL);
	//ref->vocabulary.push_back(NULL);
	char** preVocabulary = (char**)malloc(sizeof(char*)*(PREALLOCATEDVOCSIZE));  //Not working
	while(read = getline(&line, &len, handler->reftargetfile)!= -1){        
		qcount++; 
		size_t newbuflen = strlen(line);
		if (line[newbuflen - 1] == '\n')
			line[newbuflen - 1] = '\0';
		char* token = strtok(line, delim);        
		while(token != NULL && !isspace(*token)){            
			struct my_struct *sss;
			if(token[strlen(token)-1] == '\n'){
				token[strlen(token)-1] = '\0';
			}
			HASH_FIND_STR(ref->users_target, token,sss);            
			if(!sss){                                
				s = (hashtbl*)malloc(sizeof(hashtbl));                
				s->id = HASH_COUNT(ref->users_target) + 2;
				last = s->id;
				s->name = (char*)malloc(sizeof(char)*(strlen(token)+1)); 
				memcpy((void*)s->name, token, strlen(token)+1);
				HASH_ADD_KEYPTR( hh, ref->users_target, s->name, strlen(s->name), s);                
				ref_temp[count] = s->id;
				add_intchar(s->id, token, &ref->intchar_target);
				//ref->vocabulary.(char*)(s->name);
				preVocabulary[s->id] = (char*)s->name;
			} else{
				ref_temp[count] = sss->id;            
			}
			//printf("TOKEN: %s - %d - %d \n", token, strlen(s->name), qcount);           
			count++;           				  	
			token = strtok(NULL, delim);           
		}
		ref_temp[count] = 1;
		count++;
		senidx[qcount]=count;
	}
	///Put this direct access char* array into ref
	if(PREALLOCATEDVOCSIZE<(HASH_COUNT(ref->users_target) + 2)){
		fprintf(stderr, "-->Error: Make sure it is bigger. Check.\n");
		exit(0);
	}
	ref->vocabulary = (char**)malloc(sizeof(char*)*(HASH_COUNT(ref->users_target) + 2));
	memcpy(ref->vocabulary, preVocabulary, sizeof(char*)*(HASH_COUNT(ref->users_target) + 2));
	free(preVocabulary);

	ref_temp[count] = 1;
	count++;
	last++;
	ref_temp[count] = last;
	count++;
	ref->toklen = count;
	fprintf(stderr, "Target Reference toklen number is %d HASH_COUNT TARGET %d INTCHARCOUN %d\n", count, HASH_COUNT(ref->users_target), HASH_COUNT(ref->intchar_target)); 
	long bufsize = count*sizeof(int)+4;
	long intoGPU = count*sizeof(int);
	//fprintf(stderr, "%f MB Memory for target reference will be allocated on GPU\n", (double)intoGPU/(1024*1024));
	//assert(intoGPU < *devicemem);
	//*devicemem -= intoGPU;
	ref->bufsize = bufsize;
	ref->buf = (int*)malloc(ref->bufsize);
	//cudaMallocHost((void **)(&ref->buf),ref->bufsize);

	ref->sentenceind = (int*)malloc(sizeof(int)*(qcount+1));  ///Number of lines!!
	assert(ref->buf);

	ref_temp[count] = ref_temp[count+1] = ref_temp[count+2] = 0;
	memcpy(ref->buf, ref_temp, ref->toklen*sizeof(int));
	memcpy(ref->sentenceind, senidx, (qcount+1)*sizeof(int));

	ref->str = ref->buf;
	free(ref_temp);
	free(senidx);
	//fprintf(stderr, "\nSub target reference count: %d\n",refset->count);	
}

static void initRefSet(ref_set_t *refset, disk_handler_t *handler, uint	*devicemem) {
	int iseof = 0;
	ref_t *prevref = NULL;
	mytimer_t __t;

	refset->refs = (ref_t *)calloc(REF_LOAD_BATCH, sizeof(ref_t));
	assert(refset->refs);
	refset->count = 0;
	timer_start(&__t);
	fprintf(stderr, "\nLoading the reference\n");
	ref_t * ref = &refset->refs[refset->count];
	int* ref_temp = (int*)malloc(sizeof(int)*(handler->refapplen + 3));
	uint8_t * Plong = (uint8_t*)malloc(sizeof(uint8_t)*(handler->refapplen + 3));
	int* senidx = (int*)malloc(sizeof(int)*(handler->refapplen));  ///Number of lines!!
	assert(senidx != NULL);
	assert(Plong != NULL);
	assert(ref_temp!=NULL);

	hashtbl *s = NULL;    
	ref->users = NULL;
	ref->intchar = NULL;

	int count = 0;    
	int qcount = 0;    
	int read = 0;
	const char* delim = " ";    
	char* line = NULL;    
	size_t len = 0;
	int last = -1;
	uint8_t localcount = 0;
	senidx[0] = 0;
	///Direct Access to the string: from ID to char* pointers
	char** preVocabulary = (char**)malloc(sizeof(char*)*(PREALLOCATEDVOCSIZE)); 
	while(read = getline(&line, &len, handler->reffile)!= -1){        
		qcount++;        
		size_t newbuflen = strlen(line);
		if (line[newbuflen - 1] == '\n')
			line[newbuflen - 1] = '\0';
		char* token = strtok(line, delim);        
		localcount = 0;
		while(token != NULL && !isspace(*token)){            
			struct my_struct *sss;    
			if(token[strlen(token)-1] == '\n'){
				token[strlen(token)-1] = '\0';
			}
			HASH_FIND_STR(ref->users, token,sss);            
			if(!sss){                                
				s = (hashtbl*)malloc(sizeof(hashtbl));                
				s->id = HASH_COUNT(ref->users) + 2;
				last = s->id;

				s->name = (char*)malloc(sizeof(char)*(strlen(token)+1)); 
				memcpy((void*)s->name, token, strlen(token)+1);
				HASH_ADD_KEYPTR( hh, ref->users, s->name, strlen(s->name), s);                
				ref_temp[count] = s->id;
				add_intchar(s->id, token, &ref->intchar);
				preVocabulary[s->id] = (char*)s->name;
			} else{
				ref_temp[count] = sss->id;            
			}
			Plong[count] = localcount;
			localcount++;
			//printf("TOKEN: %s - %d - %d\n", token, strlen(s->name), f);           
			count++;           				  	
			token = strtok(NULL, delim);           
		}
		ref_temp[count] = 1;
		Plong[count] = 0;
		count++;
		senidx[qcount]=count;
	}
	///Put this direct access char* array into ref
	//Security check
	if(PREALLOCATEDVOCSIZE<(HASH_COUNT(ref->users) + 2)){
		fprintf(stderr, "-->Error: Make sure it is bigger. Check.\n");
		exit(0);
	}
	ref->vocabulary = (char**)malloc(sizeof(char*)*(HASH_COUNT(ref->users) + 2));
	memcpy(ref->vocabulary, preVocabulary, sizeof(char*)*(HASH_COUNT(ref->users) + 2));
	free(preVocabulary);

	ref_temp[count] = 1;
	Plong[count] = 0;
	count++;
	last++;
	ref_temp[count] = last;
	count++;
	ref->toklen = count;
	fprintf(stderr, "Reference toklen number is %d HASH_COUNT %d INTCHARCOUN %d\n", count, HASH_COUNT(ref->users), HASH_COUNT(ref->intchar)); 

	size_t bufsize = count*4*sizeof(int)+4;
	size_t sa_size = count*sizeof(int);
	size_t Psize = count*sizeof(uint8_t);
	size_t intoGPU = count*4*sizeof(int);
	fprintf(stderr, "%f MB Memory for reference will be allocated on GPU\n", (double)bufsize*4/(6*1024*1024));
	assert(intoGPU < *devicemem);
	*devicemem -= intoGPU;
	ref->bufsize = bufsize;
	ref->distinctTokenCount = HASH_COUNT(ref->users) + 2;

	/////Memory Allocation buf and sa separation!
	cudaMallocHost((void **)(&ref->P), Psize);
	cudaMallocHost((void **)(&ref->buf), bufsize);
	cudaMallocHost((void **)(&ref->sa),sa_size);
	cudaMallocHost((void **)(&ref->str),sa_size);

	/*ref->P = (uint8_t*)malloc(Psize);	
	  ref->buf = (int*)malloc(bufsize);
	  ref->sa = (int*)malloc(sa_size);
	  ref->str = (int*)malloc(sa_size);*/

	ref->sentenceind = (int*)malloc(sizeof(int)*(qcount+1));  ///Number of lines!!

	assert(ref->buf);
	ref_temp[count] = ref_temp[count+1] = ref_temp[count+2] = 0;
	memcpy(ref->str, ref_temp, ref->toklen*sizeof(int));
	memcpy(ref->P, Plong, ref->toklen*sizeof(uint8_t));
	memcpy(ref->sentenceind, senidx, (qcount+1)*sizeof(int));

	suffixArrayConstruct(ref, last, ref_temp);
	ref->sentence_count = qcount;

	free(ref_temp);
	free(Plong);
	free(senidx);
	/* Increment the count and resize the ref set
	   array in case we ran out of space */
	refset->count++;
	if ((refset->count % REF_LOAD_BATCH) == 0) {
		refset->refs = (ref_t *)realloc(refset->refs, 
				sizeof(ref_t) * 
				(refset->count + REF_LOAD_BATCH));
		assert(refset->refs);
	}
	prevref = &refset->refs[refset->count - 1];

	timer_stop(&__t);
	/*fprintf(stderr, "\nSub reference count: %d\n"
	  "Loading the reference completed: %f Sec, HASHCOUNT REF %d INTCHAR COUNT %d\n", 
	  refset->count, timer_elapsed(&__t)/1000.0, HASH_COUNT(users), HASH_COUNT(intchar));    */
}

static void 
finalizeQrySetBuffers(qry_set_t *qryset, 
		int refsetcount) {	

	cudaFreeHost(qryset->result_two);
	cudaFreeHost(qryset->qrysbuf);
	cudaFreeHost(qryset->connectoffset);
	cudaFreeHost(qryset->result_connect);
}

void recordTime(options_t * options,
		ref_set_t *refset,
		timing_t * timing, 
		float overall) {

	int i = 0;

	double total = 
		timing->suffixarray +
		timing->qrysin +
		timing->refsin +
		timing->resout +
		timing->kernel +
		timing->printout +
		timing->kernel2 + 
		timing->extractin +
		timing->extractkernel;

	if (options->timefile) {
		FILE * fh = fopen(options->timefile, "a");
		fprintf(fh, 
				"total: %f , "
				"kernel 1: %f , "
				/*"kernel 2: %f , "*/
				"printmatches: %f , "
				"querystoGPU: %f , "
				"outputfromGPU: %f , "
				"referencetoGPU: %f , "
				"queriesfromdisk: %f , "
				"refpreprocessing: %f , "
				//"reflen: %d , "
				//"subreflen: %d , "
				"subrefcount: %d , "
				"subrefsize: %d , "
				"minmatchlen: %d\n",
				total,
				timing->kernel,
				//timing->kernel2,
				timing->printout,
				timing->qrysin,
				timing->resout,
				timing->refsin,
				timing->qrysload,
				timing->suffixarray,
				//reflen,
				/*refset->refs[0].len,*/
				refset->count,
				refset->refs[0].bufsize,
				options->minmatchlen);
		fclose(fh);
	}

	/* print the results again to stderr */
	fprintf(stderr, 
			"\n\ntotal: %f , "
			"kernel 1: %f , "
			/*"kernel 2: %f , "*/
			"printmatches: %f , "
			"querystoGPU: %f, "
			"outputfromGPU: %f, "
			"referencetoGPU: %f, "
			"queriesfromdisk: %f, "
			"refpreprocessing: %f, "
			"extraction transfer: %f, "
			"extraction kernel: %f\n",
			total,
			timing->kernel,
			//timing->kernel2,
			timing->printout,
			timing->qrysin,
			timing->resout,
			timing->refsin,
			timing->qrysload,
			timing->suffixarray,
			timing->extractin,
			timing->extractkernel);

}

int getDeviceMemory() {


	cudaDeviceProp props;
	int device;
	cudaError_t err;

	cudaGetDevice(&device);
	//assert(err == cudaSuccess);

	cudaGetDeviceProperties(&props, device);
	//assert(err == cudaSuccess);

	return int((double(props.totalGlobalMem) * 95.0/100.0));

}

void 
start(options_t * options) {

	disk_handler_t handler = DISK_INIT_HANDLER;
	ref_set_t refset;
	qry_set_t qryset;


	uint devicemem = 2000*1024*1024;//getDeviceMemory();
	//cudaDeviceSetLimit(cudaLimitPrintfFifoSize, PRINTCACHELIMIT);
	//fprintf(stderr, "Print cache limit is %f MB\n", (double)PRINTCACHELIMIT/(1024*1024));

	mytimer_t __t, __overall,__overall_sub;
	timing_t timing = {
		0,
		0,
		0,
		0,
		0,
		0,
		0
	};

	/* REF_DEV_BUF_RATIO is the ratio of device memory that 
	   will be used for the reference, the rest is for queries */
	int maxreflen = 
		suffixArrayGetEquivalentMaxRefLen(int(double(devicemem)*REF_DEV_BUF_RATIO),
				options->fingerlen);

	timer_start(&__overall);
	int rc = 
		diskInit(&handler,
				options->reffile,
				options->qryfile,
				maxreflen,
				options->reftargetfile,
				options->align, 
				options->wordscdec,
				options->destinationDirectory);

	if (!rc) {
		return;
	}

	fprintf(stderr, 
			"\nAvailable Device Memory: %.3fMB\n"
			"Approximate Reference Length: %.3fMB\n"
			"Approximate Query Length: %d\n"
			"Sub Reference Length: %.3fMB\n", 
			(double)devicemem / (1024.0*1024.0),
			(double)handler.refapplen / (1024.0*1024.0),
			handler.qryapplen,
			(double)handler.refmaxlen / (1024.0*1024.0));

	timer_start(&__t);

	///Read Reference Text In, Construct the SA Array
	initRefSet(&refset, &handler, &devicemem);
	initRefTargetSet(&refset, &handler, &devicemem);

	//map<string, categ> word_score = initWordPossibility(&handler, users, users_target);
	//Already cuda memory allocated.
	lexicalFileCuda* cudaLexFile = initWordPossibilityIntKey(&handler, 	&refset);

	timing.suffixarray = timer_stop(&__t);

	//fprintf(stderr, "\nLoading query set\n");		
	timer_start(&__overall_sub);
	timer_start(&__t);
	///Read Qryset In
	constructQryIndex(&qryset, &handler, &refset);
	timing.qrysload += timer_stop(&__t);
	//fprintf(stderr, "\nLoading query complete\n");

	//////////Reading the Alignment Files - This is an one-time offline operation, so that timings can be deducted
	timer_start(&__t);
	initAlignment(&refset, &handler, &devicemem);		
	double alignment_deduction = timer_stop(&__t); ////

	void * ctx = suffixArraySearchInit(
			refset.refs[0].toklen*sizeof(int)*2, 
			qryset.qrysbufsize, 
			qryset.totaltokens*sizeof(int), 
			refset.refs[0].toklen*sizeof(int));
	int iseof = 0;
	int flag = 0;
	int twoway = 1;
	/* Actual matching */
	suffixArraySearch(ctx, &refset, &qryset, 1, &timing, twoway);
	/* Print out the results */
	timer_start(&__t);
	//printResults_two_file(&refset, &qryset, options, users, intchar);
	timing.printout += timer_stop(&__t);
	/////////////////////////////////// Extract Pair
	fprintf(stderr, "Start Extract Pair\n");
	//////////////////////////////////////////////
	suffixArraySearchFinalize_One(ctx);////Array Cleanup - timings are included.
	ExtractPairs_Large_Data_Gappy(ctx, 
			&refset, 
			&qryset,
			&timing, 
			cudaLexFile,
			&handler);		//Changed to new function with gappy.
	///////////////////////////////////
	//fprintf(stderr, "Finish Extract Pair\n");

	timer_stop(&__overall_sub);
	timer_stop(&__overall);
	/*fprintf(stderr, 
			"Excluding Preprocessing Time - Qry Loading&Construction/Suffix Extraction/Consistent Pairs Extraction/Pair Scoring&Counting/All GPU Costs/All Rules PrintOut/Disk IO %f - TOTAL %f\n", 
			timer_elapsed(&__overall_sub)-alignment_deduction, timer_elapsed(&__overall)-alignment_deduction);*/

	print_query_GPU_Gappy(
			qryset.continousQueryWithID,		
			qryset.oneGapQueryWithID,
			qryset.twoGapQueryWithID,
			qryset.qryscount,
			qryset.globalOnPairsUpDownContinous,
			qryset.globalOnPairsUpDownGappy,
			qryset.globalOnPairsUpDownTwoGap,
			qryset.fast_speed,
			qryset.fast_speed_one_gap,
			qryset.fast_speed_two_gap,
			options->destinationDirectory,
			qryset.distinctContinousCount,//continous noGapsearch count
			qryset.distinctOneGapCount,//oneGapSearch count
			qryset.distinctTwoGapCount);//TwoGapSearch count
	/*print_query_GPU_Continous(
	  qryset.continousQueryWithID, 
	  qryset.qryscount,
	  qryset.globalOnPairsUpDownContinous,
	  qryset.fast_speed,
	  options->destinationDirectory);*/
	//recordTime(options, &refset, &timing, timer_elapsed(&__overall));

	// clean up!
	extractPairFinalize(ctx);
	suffixArraySearchFinalize(ctx);
	finalizeQrySetBuffers(&qryset, refset.count);
	finalizeRefSet(&refset);
	diskFinalize(&handler);
}



