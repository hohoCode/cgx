#ifndef GAPPY_LOOK_H
#define GAPPY_LOOK_H

#include "ComTypes.h"

	__global__ void precomp(
			precompute_enu* oneGapPrecomp,
			int* refstr, 
			int *refsa,
			precompute_enu_2* onegap_precomp,
			unsigned int* counter,
			unsigned int* RLP,
			uint8_t* L_tar,
			uint8_t* R_tar,
			int* featureMissingCount_d);

__global__ void oneGapLookUpSA(	
		int* refstr, 
		int* refsa,
		unsigned int toklen,
		int *qrysoffsettok,
		int qryscount,    
		int tokenscount,
		result_t_two* qryresult,
		int* connectoffset,
		int totalconnect,
		result_t* result_connect,
		int* tokindex_qryindex,
		oneGapOnSA* oneGapSA,	
		unsigned int precomp_onegap_count,
		int* frequentList,
		precomp_st_end* precomp_index,
		precompute_enu_3* precomp_onegap,
		unsigned int* counter,
		gappy_search* oneGapSearch,
		int distinctOneGapCount,
		unsigned int* RLP,
		uint8_t* L_tar,
		uint8_t* R_tar);
		
__global__ void twoGapLookUpSA(	
			int* refstr, 
			//	int* refsa,
			unsigned int toklen,
			int *qrysoffsettok,
			int qryscount,    
			int tokenscount,
			result_t_two* qryresult,
			int* connectoffset,
			int totalconnect,
			result_t* result_connect,
			int* tokindex_qryindex,
			twoGapOnSA* twoGapSA,	
			unsigned int precomp_onegap_count,
			int* frequentList,
			precomp_st_end* precomp_index,
			precompute_enu_3* precomp_onegap,
			unsigned int* counter,
			two_gappy_search* twoGapSearch,
			int distinctTwoGapCount,
			int distinctOneGapCount,
			gappy_search* oneGapSearch,
			oneGapOnSA* oneGapSA,
			unsigned int oneGapSACount,
			unsigned int* RLP,
			uint8_t* L_tar,
			uint8_t* R_tar);
			
#endif /* GAPPY_LOOK_H */
