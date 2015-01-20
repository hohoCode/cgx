#ifndef PRINT_RESULTS_H
#define PRINT_RESULTS_H

#include "ComTypes.h"

using namespace std;

void printResults_two(
		ref_set_t * refset, 
		qry_set_t *qryset, 
		options_t * options, 
		hashtbl* uu, 
		hash_intchar* gg);
		
void printResults_two_file(
		ref_set_t * refset, 
		qry_set_t *qryset, 
		options_t * options, 
		hashtbl* uu, 
		hash_intchar* gg);

void print_query_GPU_Continous( 		
		std::vector<std::vector<unsigned int> >& queryglobal, 
		int qrycounter,
		result_t* globalOnPairsUpDown,
		red_dup_t* fast_speed,
		char* destDir);

void print_query_GPU_Gappy( 		
		vector<vector<unsigned int> > &queryglobal,
		vector<vector<unsigned int> > &queryIdGappy, 
		vector<vector<unsigned int> > &queryIdTwoGap,	
		int qrycounter,
		result_t* globalOnPairsUpDownContinous,
		result_t* globalOnPairsUpDownGappy,
		result_t* globalOnPairsUpDownTwoGap,
		red_dup_t* fast_speed,
		red_dup_t* fast_speed_one,
		red_dup_t* fast_speed_two,
		char* destDir,
		unsigned int global,
		unsigned int distinctOneGapCount,
		unsigned int distinctTwoGapCount);
		

#endif /* PRINT_H */
