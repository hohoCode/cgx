#ifndef COM_TYPES_H
#define COM_TYPES_H

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <assert.h>
#include <ctype.h>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <iterator>
#include <algorithm>
#include <map>
#include <vector>
#include "uthash/uthash.h"
#include "uthash/utarray.h"
#include "Timer.h"

/* this constant defines the ratio of device memory that 
will be used for the reference, the rest is for queries */
#define REF_DEV_BUF_RATIO (9.3/10.0)

/* the ratio of device memory reported by cuda runtime 
that is available for allocation. I estimated this 
ratio experimentally */
#define FREE_GLOBAL_MEM_RATIO 95.0/100.0

/* defines the maximum query and query-name lengths */
#define QRY_MAX_LENGTH (1024)
#define QRY_MAX_NAME_LENGTH (1024)

/* defines the state needed for each query: 1 for the null, and 2 ints 
for the offset and the length */
#define QRY_STATE_SIZE (1 + 2 * sizeof(int))

// calculates the space needed to store the results in the GPU
#define RESULT_DATA_SIZE(qrystotlen, qryscount, minmatchlen) \
	((qrystotlen - (qryscount * (minmatchlen - 1))) * sizeof(result_t))

#define MAX_rule_span 15//original 15
#define MAX_rule_span_pattern 15
#define MAX_rule_symbols 5

#define MAX_rule_span_double 28//Double the Max_rule_span for in kernel memory purpose
#define MIN_gap_size 1
#define MAX_nonterminals 2


#define MAXSCORE 99 //Feature Lexical: MaxEF

/***ALL OCCURNECES String Matching****/
#define ONEGAP_PREALLOCATION 60000000
#define PRECOMPUTECOUNT 100 //200
#define ONEGAP_PRECOMPUT_PREALLOCATION 60000000
#define TWOGAP_PREALLOCATION 70000000  
#define LINEARBLOCK 700
#define ONEGAP_ENU_PREALLOCATION 20000000
#define TWOGAP_ENU_PREALLOCATION 35000000

#define isSample true
#define SAMPLER 300
#define SAMPLER_ONEGAP 65
#define SAMPLER_TWOGAP 70

typedef struct options_s {

	char *reffile;
	char *qryfile;
	char *reftargetfile;
	char *align;
	char *timefile;
	char *wordscdec;
	int	 minmatchlen;
	int	 fingerlen;
	char* destinationDirectory;
} options_t;

typedef struct timing_s {
	double suffixarray;
	double qrysload;
	double qrysin;
	double refsin;
	double resout;
	double kernel;
	double printout;
	double kernel2;
	double extractin;
	double extractkernel;
} timing_t;

typedef struct result_s {
	int up;
	int down;
} result_t;

typedef struct result_s_two {
	int up;
	int down;
	int firstfindhit;
	int firstfindhitL;
    int firstfindhitR;
	int longestind;
	int longestmatch;
} result_t_two;

typedef struct qrySID1{
	unsigned int index;
	uint8_t type; 
	//1. precomp -> based on one gap seed, get one gap rule: 1->1/2 precomp
	//2. precomp -> based on two gaps seed, get twp gap rule: 2->2 precomp
	//3. no precomp -> based on one gap seed, get one gap rule: 1->1/2 no precomp
	//4. no precomp -> based on two gaps seed, get two gap rule: 2->2 no precomp
} qrySID;

typedef struct my_struct_lexicon {
	const char  *name;          /* key */
	int* alinged;  //Only those aligned
	//int all_suffix; ///All matches
	UT_hash_handle hh;         /* makes this structure hashable */
} hash_lexicon;

typedef struct my_struct {
	const char  *name;          /* key */
	int id;
	UT_hash_handle hh;         /* makes this structure hashable */
} hashtbl;

typedef struct my_struct2 {
	const char  *name;          
	int id;                   /* key */
	UT_hash_handle hh;         /* makes this structure hashable */
} hash_intchar;

typedef struct my_struct3 {   /// Only Used in CreatedLexicon() in Extractpairs.cu, to get the target string faster
	const char  *name;          /* key */
	int* id;
	UT_hash_handle hh;         /* makes this structure hashable */
} hashtbl_aux;

typedef struct __attribute__((__packed__)) gappytyp {
//typedef struct gappytyp { //ROLL BACK 
//If all four adds up, should minus 1 get the last position of aXb.
	int qrystart;
	uint8_t qrystart_len;
	uint8_t qryend_len;
	uint8_t gap;
} gappy;//gappy enum

typedef struct __attribute__((__packed__)) gappytyp2 {
//typedef struct gappytyp2 { //ROLL BACK
	unsigned int blockid; //for the one gappy search index id
	unsigned int gap2;
	uint8_t qryend_len;
} twogappy;//gappy enum

typedef struct __attribute__((__packed__)) gappytyp_se2 {
//typedef struct gappytyp_se2 { //ROLL BACK
	unsigned int blockid; //for the one gappy search index id
	unsigned int gap2;///Change it to unsigned int, no longer uint8_t
	uint8_t qryend_len;
	unsigned int position; //positions on the twogap array
	int start_on_salist;
	int end_on_salist;
} two_gappy_search;//gappy to be searched

typedef struct __attribute__((__packed__)) gappytyp_se {
//typedef struct gappytyp_se { //ROLL BACK
	int qrystart;
	uint8_t qrystart_len;
	uint8_t qryend_len;
	uint8_t gap;
	unsigned int position; //positions on the onegap array
	int start_on_salist;
	int end_on_salist;
} gappy_search;//gappy to be searched

typedef struct __attribute__((__packed__)) gappytyp_sa {
//typedef struct gappytyp_sa { //
	unsigned int position;// this is the block ID on onegapsearch array! Not that position
	unsigned int str_position;
	uint8_t length; //on str array
} oneGapOnSA;//one gappy results on str

typedef struct __attribute__((__packed__)) gappytyp_sa2 {
//typedef struct gappytyp_sa2 { //
	unsigned int position;// this is the block ID on twogapsearch array! Not that position
	unsigned int str_position;
	uint8_t length; //on str array - same as above - first end of the string aXb's b
	uint8_t length2; //on str array - the extra second real end of the stirng aXbXc's c
} twoGapOnSA;//one gappy results on str

typedef struct __attribute__((__packed__)) gapPattern1 {
//typedef struct gapPattern1{ //ROLL BACK
	int pattern[MAX_rule_symbols];	
	uint8_t number; // 1 - Max_rule_symbols
} gapPattern;

typedef struct __attribute__((__packed__)) gapPattern2 {
//typedef struct gapPattern2{ //ROLL BACK
	int pattern[MAX_rule_symbols-4];	
	uint8_t number; // 1 - Max_rule_symbols
	unsigned int blockid;//for sort purpose - one gap block id
} twoGapPattern;

typedef struct __attribute__((__packed__)) gappytyp2old {
//typedef struct gappytyp2old { //ROLL BACK
	unsigned int gappy_index;
	int qryend2;
	uint8_t qryend_len;
	uint8_t gap;
} two_gappy;

typedef struct __attribute__((__packed__)) gappytyp3 {
//typedef struct gappytyp3 { //ROLL BACK
	unsigned int gappy_index;
	unsigned int within_each_precomp;
	int qryend2;
	uint8_t qryend_len;
	uint8_t gap;
} two_gappy_precomp;

typedef struct __attribute__((__packed__)) gappytyp4 {
//typedef struct gappytyp4 { //ROLL BACK
	int gappy_index;
	unsigned int ref_str_start; //start
	uint8_t end;
	uint8_t gap1;//use gap1 and gap2 to determine the gap length, one of them could be zero.
	uint8_t gap1_1;	
} rule_onegap;

typedef struct __attribute__((__packed__)) gappytyp5 {
//typedef struct gappytyp5 { //ROLL BACK
	int twogappyindex; //detect whether it is bigger than two_gappy_precomp's count first
	unsigned int ref_str_start;
	uint8_t end;
	uint8_t gap1; //first gap length
	uint8_t gap1_1;
	uint8_t gap2;
	uint8_t gap2_1;	
} rule_twogap;

typedef struct red_dup_s { // CUDA output
	int blocknumber;	
	const char* lexical;
	int f; // all aligned source
	float all_suffix_fsample_score; //source with all suffix
	int all_suffix_fsample; //source with all suffix
	int* paircount; //fe pair count
	float bb; //source with all suffix
	float aa; //source with all suffix
	float MaxLexFgivenE;
	float MaxLexEgivenF;
} red_dup_t;

typedef struct qry_set_s {

	int *qrysbuf; 
	int qrysbufsize;

	int *qrys;

	int		 *connectoffset;
	int		 *qrysoffsettok;	/* list of queries offsets 
								within the qrys buffer */
	int		  qryscount;	/* number of queries */

	int		resbufsize;	/* size of the result buffer */

	int totaltokens;
	int totalconnect;
	result_t* result;
	result_t* result_connect;
	result_t_two* result_two;
	int* tokindex_qryindex;

	//Continous Data	
	////////////////////////////////////	
	std::vector<std::vector<unsigned int> > continousQueryWithID;
	unsigned int distinctContinousCount;

	//Gappy Data Structures
	////////////////////////////////////	
	unsigned int onegapcount_enu;//Done
	gappy* onegap;//Done Sorted
	gapPattern* onegapPattern;//Done sorted
	
	unsigned int distinctOneGapCount;//Done
	gappy_search* oneGapSearch;//Done

	oneGapOnSA* oneGapSA;//Done	sorted
	unsigned int countOneGapSA;//Done
	///////////////////////////////////////////
	unsigned int twogapcount_enu;//Done
	twogappy* twoGap;//Done
	gapPattern2* twogapPattern;//Done
	
	unsigned int distinctTwoGapCount;//Done
	two_gappy_search* twoGapSearch;//Done

	twoGapOnSA* twoGapSA;//Done
	unsigned int countTwoGapSA;//Done

	std::vector<std::vector<unsigned int> > twoGapQueryWithID;
	std::vector<std::vector<unsigned int> > oneGapQueryWithID;

	//two_gappy* twogaps;
	//unsigned int twogapscount;

	///////////////////////////////////////////Processed Up and Down with converted ids
	result_t* globalOnPairsUpDownContinous;//bnum - ab
	result_t* globalOnPairsUpDownGappy;//one gap - aXb, Xab, abX
	result_t* globalOnPairsUpDownTwoGap;//two gap - aXbXc, XaXb, aXbX

	////////////////////////////////Statistics all in here for contionous/one gap/two gap
	red_dup_t* fast_speed;
	red_dup_t* fast_speed_one_gap;
	red_dup_t* fast_speed_two_gap;
} qry_set_t;

typedef struct precomp_st_end1{
	unsigned int start;    
	unsigned int end;
} precomp_st_end;

typedef struct ref_s_tar {
	int	*buf;     
	long	 bufsize; /* size of the buffer */
	int *str;	/* reference string */
	uint8_t *L_tar;
	uint8_t *R_tar;
	unsigned int toklen;    /*the length of the sa, lcp, lcpcr, lcpcl and rank arrays */
	int* sentenceind;
	char** vocabulary;//from s->id to s->name on target lan side for faster speed	

	//Change global to local
	hashtbl *users_target;
	hash_intchar *intchar_target;
} ref_t_target;

typedef struct saind_s { // CUDA input
	int start;
	int end;	
	int matchlen;
	int string_start;
} saind_t;

typedef struct __attribute__((__packed__)) res_phrase_s { // CUDA output
	int tar_start;
	int blocknumber;
	uint8_t tar_end;
} res_phrase_t;

typedef struct workKeytyp {
     int ch; //first key using int as ids
     int eng; //second key
} wordKey;

typedef struct categ {
     float val1;
     float val2;
} category;

typedef struct taskRes_t { //search final results here
	 float MaxLexFgivenE;
	 float MaxLexEgivenF;
} lexTaskResults;

typedef struct lexicalFile { //search range and results here
     categ* values; 
     workKeytyp* keys; //each contain two parts
     unsigned int count;
} lexicalFileCuda;

typedef struct __attribute__((__packed__)) lexicalTaskt {
//typedef struct lexicalTaskt {
     unsigned int fastSpeedId;
	 ///Source side
	 int sourcePattern[MAX_rule_symbols];  
	 uint8_t sourcePatternCounter;
	 ///Target side
	 unsigned int targetStart;
	 uint8_t end;
	 uint8_t gap1; //first gap length
	 uint8_t gap1_1;
	 uint8_t gap2;
	 uint8_t gap2_1; 
} lexicalTask;

typedef struct targetCate{
     int target_start;
     const char* str;
} tarCate;

typedef struct precompute1 {
     int length;
     int token;
	 long start;
	 long end;
} precomp_tok;

typedef struct __attribute__((__packed__)) precompute_enu{
//typedef struct precompute_enu {
        int token_a;//aXb
        int token_b;
		unsigned int start;
		unsigned int length;
		bool reverse;
} precompute_enu;

typedef struct __attribute__((__packed__)) precompute_enu_21{
//typedef struct precompute_enu_2{ //roll back
        unsigned int index;//aXb
		unsigned int start;
		uint8_t length;
} precompute_enu_2;

typedef struct __attribute__((__packed__)) precompute_enu_31{
//typedef struct precompute_enu_31{ ////ROLL BACK
		unsigned int start;
		uint8_t length;//The last position of aXb's b. inclusive
} precompute_enu_3;

typedef struct ref_s {

	int	*buf;     /* This buffer contains all the state presented below 
				  contiguously in the following order: lcp, lcpleft, lcpright, rank
				  This would make it easier when moving data to the GPU */
	size_t bufsize; /* size of the buffer */

	int *str;	/* reference string, separate */

	int	*sa;	/* reference suffix array separate */
	int	*lcp;	/* longest common prefix array */
	int  *rank;
	int *lcpleft;  /*LCP Table Random Version */
	int *lcpright;
	unsigned int toklen;    /*the length of the sa, lcp, lcpcr, lcpcl and rank arrays */

	/* If this is a sub reference in a super reference then: */
	int begin;   /* is the beginning index of this 
				 sub reference in the super reference */
	int end;     /* is the end index of this 
				 sub reference in the super reference */
	uint8_t *P;
	unsigned int *RLP;
	int *sentenceind;
	int sentence_count;
	int distinctTokenCount; //HashCount(user)+2

	////Precomputation
	int* frequentList;
	precomp_st_end* precomp_index;//size: TOP*TOP
	precompute_enu_3* precomp_onegap;
	unsigned int precomp_count;
	int* featureMissingCount; //size: TOP*TOP
	
	char** vocabulary;//from s->id to s->name on source lan side for faster speed

	///Change global to local	
	hashtbl *users;
	hash_intchar *intchar;
} ref_t;

typedef struct ref_set_s {
        ref_t * refs;
        ref_t_target * refs_target;
        int count;
} ref_set_t;

typedef struct context_s {
        ref_t ref_d;
        ref_t_target ref_target_d;
        qry_set_t qryset_d;
        result_t *resbuf_d;
        saind_t* blocks_d;
        result_t * globalOnPairsUpDown_d;
        res_phrase_t* out_pair;		
} context_t;

#endif /* COM_TYPES_H */
