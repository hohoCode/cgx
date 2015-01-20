#ifndef EXtract_H
#define EXtract_H

#include "ComTypes.h"
#include "Disk.h"
#include <stdlib.h>
#include <string>
#include <stdio.h>
#include <algorithm>
#include <vector>
#include <iostream>

using namespace std;


void initAlignment(ref_set_t *refset, 
        disk_handler_t *handler, 
        uint *devicemem);

lexicalFileCuda* initWordPossibilityIntKey(disk_handler_t *handler, 
		ref_set_t* refset);

std::map<std::string, categ> initWordPossibility(disk_handler_t *handler,
		hashtbl* refstrid, 
		hashtbl* targetstrid);

void ExtractPairs(void *handler, 
        ref_set_t *refset, 
        qry_set_t *qryset, 
        timing_t *timing, 
        std::map<std::string, categ> worder);

void extractPairFinalize(void *handler);

red_dup_t* createLexicon( int count, 
        res_phrase_t* out_res, 
        int* source_str, 
        int* target_str, 
        hashtbl_aux** lexic, 
        hash_lexicon** target_count, 
        hash_lexicon** foreig_count, 
        int* lexicon_index, 
        saind_t * tmp_blocks,
        map<string, categ> word_score,
        vector<char*> sourceName,
        int global,
        red_dup_t* fast_speed,
        char** target_vocabulary);

red_dup_t* createLexiconGappy(
        int* source_str, 
        int* target_str, 
        int* lexicon_index, 
        saind_t * tmp_blocks,//gapSearch for continous phrases
        map<string, categ> word_score,
        vector<char*> sourceName,
        int global,
        rule_onegap* oneGapRule,
				rule_twogap* twoGapRule,
				unsigned int count_one_gap,//Total number of one gap results
				unsigned int count_two_gap,//Total number of two gap results
				int	seperatorOneGap,//Seperator between aXb; Xab, abX
				int	seperatorTwoGap_one,//XabX;aXbXc
				int seperatorTwoGap_two,//aXbXc;XaXb,aXbX
				unsigned int distinctOneGapCount,//Distinct pattern count
				unsigned int distinctTwoGapCount,//pattern count - gappy ID max
				gapPattern* onegapPattern,//source side pattern
				gapPattern2* twogapPattern,
				char** target_vocabulary,
				char** source_vocabulary,
				red_dup_t* fast_speed,
				gappy_search* oneGapSearch,
				two_gappy_search* twoGapSearch);

red_dup_t* createLexiconTwoGap(
        int* source_str, 
        int* target_str, 
        int* lexicon_index, 
        saind_t * tmp_blocks,//gapSearch for continous phrases
        map<string, categ> word_score,
        vector<char*> sourceName,
        int global,
        rule_onegap* oneGapRule,
				rule_twogap* twoGapRule,
				unsigned int count_one_gap,//Total number of one gap results
				unsigned int count_two_gap,//Total number of two gap results
				int	seperatorOneGap,//Seperator between aXb; Xab, abX
				int	seperatorTwoGap_one,//XabX;aXbXc
				int seperatorTwoGap_two,//aXbXc;XaXb,aXbX
				unsigned int distinctOneGapCount,//Distinct pattern count
				unsigned int distinctTwoGapCount,//pattern count - gappy ID max
				gapPattern* onegapPattern,//source side pattern
				gapPattern2* twogapPattern,
				char** target_vocabulary,
				char** source_vocabulary,
				red_dup_t* fast_speed,
				gappy_search* oneGapSearch,
				two_gappy_search* twoGapSearch);

red_dup_t* createLexiconFast( int count, 
        res_phrase_t* out_res, 
        int* source_str, 
        int* target_str, 
        hashtbl_aux** lexic, 
        hash_lexicon** target_count, 
        hash_lexicon** foreig_count, 
        int* lexicon_index, 
        saind_t * tmp_blocks,
        lexicalTask* maxEF,
        vector<char*> sourceName,
        int global,
        red_dup_t* fast_speed,
        char** target_vocabulary,
        unsigned int* lexicalTaskCounter);

red_dup_t* createLexiconGappyFast(
        int* source_str, 
        int* target_str, 
        int* lexicon_index, 
        saind_t * tmp_blocks,//gapSearch for continous phrases
		lexicalTask* maxEF,
        vector<char*> sourceName,
        int global,
        rule_onegap* oneGapRule,
				rule_twogap* twoGapRule,
				unsigned int count_one_gap,//Total number of one gap results
				unsigned int count_two_gap,//Total number of two gap results
				int	seperatorOneGap,//Seperator between aXb; Xab, abX
				int	seperatorTwoGap_one,//XabX;aXbXc
				int seperatorTwoGap_two,//aXbXc;XaXb,aXbX
				unsigned int distinctOneGapCount,//Distinct pattern count
				unsigned int distinctTwoGapCount,//pattern count - gappy ID max
				gapPattern* onegapPattern,//source side pattern
				gapPattern2* twogapPattern,
				char** target_vocabulary,
				char** source_vocabulary,
				red_dup_t* fast_speed,
				gappy_search* oneGapSearch,
				two_gappy_search* twoGapSearch,
				unsigned int* lexicalTaskCounter,
				oneGapOnSA* oneGapSA,
				ref_t* refSource);

red_dup_t* createLexiconTwoGapFast(
        int* source_str, 
        int* target_str, 
        int* lexicon_index, 
        saind_t * tmp_blocks,//gapSearch for continous phrases
		lexicalTask* maxEF,
        vector<char*> sourceName,
        int global,
        rule_onegap* oneGapRule,
				rule_twogap* twoGapRule,
				unsigned int count_one_gap,//Total number of one gap results
				unsigned int count_two_gap,//Total number of two gap results
				int	seperatorOneGap,//Seperator between aXb; Xab, abX
				int	seperatorTwoGap_one,//XabX;aXbXc
				int seperatorTwoGap_two,//aXbXc;XaXb,aXbX
				unsigned int distinctOneGapCount,//Distinct pattern count
				unsigned int distinctTwoGapCount,//pattern count - gappy ID max
				gapPattern* onegapPattern,//source side pattern
				gapPattern2* twogapPattern,
				char** target_vocabulary,
				char** source_vocabulary,
				red_dup_t* fast_speed,
				gappy_search* oneGapSearch,
				two_gappy_search* twoGapSearch,
				unsigned int* lexicalTaskCounter,
				oneGapOnSA* oneGapSA,
				ref_t* refSource);	

void ExtractPairs_Large_Data(void *handler,
        ref_set_t *refset,
        qry_set_t *qryset,
        timing_t *timing,
        std::map<std::string, categ> worder,
        disk_handler_t *diskInfo);

void ExtractPairs_Large_Data_Gappy(void *handler,
        ref_set_t *refset,
        qry_set_t *qryset,
        timing_t *timing,
        lexicalFileCuda* cudaLexFile,//already cuda memory allocated
        //std::map<std::string, categ> worder,
        disk_handler_t *diskInfo);

#endif /* EXtract_H END */

