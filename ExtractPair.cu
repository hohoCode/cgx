#include "ComTypes.h"
#include "ExtractPair.h"
#include <thrust/device_vector.h>
#include <thrust/sort.h>
#include <thrust/version.h>

using namespace std;

#define THREADS_PER_BLOCK 512
#define THREADS_PER_BLOCK2 128

#define THREADS_PER_BLOCK_GLOBALPAIRS 64
#define TARGETALLOCATION 10000000 
#define ROUND(X) (int)(X+0.5)
#define LONGEST 15
#define LONGESTCHSOURCE 5
#define N 1 
#define M 1///Iterations for Larger Data inisde Kernel
#define TEMPSET 5 //Physical memory limitations

#define ALLOCATIONS 50000000 //64507758 - fbisUN no sample - continous out_pair
#define PREALLOCATION_TWOGAP_RULE 50000000 //OneGapRule pre
#define PREALLOCATION_ONEGAP_RULE 50000000 //twoGapRule pre

#define LEXFILE_LINE_PRE 5000000 //Read Lex File, pre
#define PREALLOCATION_LEX_TASK 50000000 //Lexical Task, sum of all lexicon counts

struct lexFileCompare
{
	__host__ __device__
	bool operator()(const workKeytyp& a, const workKeytyp& b) const
	{
		return (a.ch < b.ch) || ((a.ch== b.ch) && (a.eng < b.eng) );
	}
};

int compareUser(const void *v1, const void *v2)
{
	const red_dup_t* u1 = (red_dup_t*) v1;
	const red_dup_t* u2 = (red_dup_t*) v2;
	return u1->blocknumber > u2->blocknumber;
}

int compareUserTotal(const void *v1, const void *v2)
{
	const res_phrase_s* u1 = (res_phrase_s*) v1;
	const res_phrase_s* u2 = (res_phrase_s*) v2;
	return (u1->blocknumber > u2->blocknumber) || 
		((u1->blocknumber == u2->blocknumber) && (u1->tar_start > u2->tar_start));
}

__host__ __device__ saind_t binary_search_hit(red_dup_t* a, int low, int high, int target) {
	saind_t output;
	output.end = output.matchlen = output.start = -1;
	while (low <= high) {
		int middle = low + (high - low)/2;
		if (target < a[middle].blocknumber)
			high = middle - 1;
		else if (target > a[middle].blocknumber)
			low = middle + 1;
		else {
			output.start = low;
			output.end = high;
			output.matchlen = middle;
			return output;
		}
	}
	return output;
}

__host__ __device__ int binary_search_up(red_dup_t* a, int low, int high, int target) {
	int temp = low;
	while (low <= high) {
		int middle = low + (high - low)/2;
		if (target < a[middle].blocknumber)
			high = middle - 1;
		else if (target > a[middle].blocknumber)
			low = middle + 1;
		else {
			temp = middle;
			low = middle + 1;
		}
	}
	return temp;
}

__host__ __device__ int binary_search_down(red_dup_t* a, int low, int high, int target) {
	int temp = high;
	while (low <= high) {
		int middle = low + (high - low)/2;
		if (target < a[middle].blocknumber)
			high = middle - 1;
		else if (target > a[middle].blocknumber)
			low = middle + 1;
		else {
			temp = middle;
			high = middle - 1;
		}
	}
	return temp;
}

__host__ __device__ bool consistent(int start, int end, uint8_t* L_target, uint8_t* R_target, int start_chk, int end_chk, int startpos_source){
	int k;
	unsigned char L;
	unsigned char R;
	bool ok = true;
	unsigned char min_L = 255;
	unsigned char max_R = 0;
	for(k = start; k <= end && ok; k++){
		L = L_target[k];
		R = R_target[k];
		
		if (L==255 || R == 255){
			ok = true; //Change fixes original false to remove non-recog symbols.
		} else if (k == start){
			min_L = L;
			max_R = R;
		} else {
			if (min_L > L) {
				min_L = L;
			}
			if (max_R < R) {
				max_R = R;
			}
		} 
	}
	
	if(startpos_source+min_L != start_chk || startpos_source+max_R != end_chk){
		ok = false;
	}
	return ok;
}

__device__ bool checkBoundaryFast(//Need to return extra info, but no target checking needed
		unsigned int start, 
		unsigned int ender,
		uint8_t* L_tar,
		uint8_t* R_tar,
		unsigned int* RLP,
		unsigned char* min_LL,
		unsigned char* max_RR,
		int* sen_target_begin,
		int* tempind){
	unsigned char L = 0;
	unsigned char R = 0;

	(*sen_target_begin) = -1;
	unsigned char min_L = 255;
	unsigned char max_R = 0;
	(*tempind) = 0;
	unsigned int temp;	

	for(int k = start; k <= ender; k++){
		temp = RLP[k];
		L = (temp >> 24) & 0xFF;
		R = (temp >> 16) & 0xFF;
		if ( (L == 255 || R == 255) && (k == start || k == ender) ){
			k = ender + 1;
			return false;			
		} else if ( (L == 255 || R == 255) ) {
			L = 255;     
		} else if (k == start){		
			(*tempind) = k - ((temp >> 8) & 0xFF) - 1 ;	
			
			if ((*tempind) == -1){
				(*sen_target_begin) = 0;
			} else {
				(*sen_target_begin) = RLP[*tempind];
			}
			min_L = L;
			max_R = R;
		} else {
			if (min_L > L) {
				min_L = L;
			}
			if (max_R < R) {
				max_R = R;
			}
		} 
	}

	if (min_L <= max_R && max_R - min_L < MAX_rule_span){
		(*tempind)++;
		*min_LL = min_L;
		*max_RR = max_R;

		//This consistent checking can be turned off.
		//printf("target start %d - tar end %d - sour start %d - sour end %d\n", ss, tt, current_str, t);
		return true;
	}
	
	return false;
}

__device__ bool checkBoundaryFast2(//Normal interface, but no target checking needed
		unsigned int start, 
		unsigned int ender,
		uint8_t* L_tar,
		uint8_t* R_tar,
		unsigned int* RLP,
		unsigned int* target_start,
		unsigned int* target_end){
	unsigned char L = 0;
	unsigned char R = 0;
	int sen_target_begin = -1;
	unsigned char min_L = 255;
	unsigned char max_R = 0;
	int tempind = 0;
	unsigned int temp;	

	for(int k = start; k <= ender; k++){
		temp = RLP[k];
		L = (temp >> 24) & 0xFF;
		R = (temp >> 16) & 0xFF;
		if ( (L == 255 || R == 255) && (k == start || k == ender) ){
			k = ender + 1; 
			return false;
		} else if ( (L == 255 || R == 255) ) {
			L = 255;     
		} else if (k == start){		
			tempind = k - ((temp >> 8) & 0xFF) - 1 ;	
			//printf("L %u R %u - BK x %d y %d - bnum %d - tempind %d - k %d - Pi %u\n", L, R, blockIdx.x, blockIdx.y, bnum, tempind, k, ((RLP[k] >> 8) & 0xFF));
			if (tempind == -1){
				sen_target_begin = 0;
			} else {
				sen_target_begin = RLP[tempind];
			}
			min_L = L;
			max_R = R;
		} else {
			if (min_L > L) {
				min_L = L;
			}
			if (max_R < R) {
				max_R = R;
			}
		} 
	}

	(*target_start) = min_L + sen_target_begin;
	(*target_end) = max_R + sen_target_begin;

	if (min_L <= max_R && max_R - min_L < MAX_rule_span){
		tempind++;
		return true;//consistent( (*target_start),  (*target_end), L_tar, R_tar, start, ender, tempind);		
	}			

	return false;
}

__device__ uint8_t checkBoundary( //Target checking needed, and need to tell error codes difference
		unsigned int start, 
		unsigned int ender,
		uint8_t* L_tar,
		uint8_t* R_tar,
		unsigned int* RLP,
		unsigned int* target_start,
		unsigned int* target_end){
/*Return error codes:
0: normal false - aXb false only.
1: true and ok
2: front error - aXbX false, aXb false
3: end error - XaXb false, aXb false
4: both front and end error - XaXb false, aXb false, aXbX false
*/
	unsigned char L = 0;
	unsigned char R = 0;
	int sen_target_begin = -1;
	unsigned char min_L = 255;
	unsigned char max_R = 0;
	int tempind = 0;
	unsigned int temp;	
	
	uint8_t frontOrEndWrong = 0; //fix the corner case bug; aXb first a is not good or last b is not good.

	for(int k = start; k <= ender; k++){
		temp = RLP[k];
		L = (temp >> 24) & 0xFF;
		R = (temp >> 16) & 0xFF;
		if ( (L == 255 || R == 255) && (k == start || k == ender) ){
			if(start==ender&&frontOrEndWrong==0){
				frontOrEndWrong =4; //all XaXb false, aXb false, aXbX false.
			}else if(frontOrEndWrong==0&&k==start){
				frontOrEndWrong =2; //aXbX false, aXb false.	
			} else if (frontOrEndWrong==0&&k==ender){
				frontOrEndWrong =3;	//XaXb false, aXb false.
			} else if (frontOrEndWrong!=0){
				frontOrEndWrong =4;	//all XaXb false, aXb false, aXbX false.
			} else {
				printf("Not possible inside checkBoundary function\n");
				return 4;
			}			
			//Cannot do this to prestop the iteration, still need min_L and max_R

			if (k == start){		
				tempind = k - ((temp >> 8) & 0xFF) - 1 ;				
				if (tempind == -1){
					sen_target_begin = 0;
				} else {
					sen_target_begin = RLP[tempind];
				}
			}
		} else if ( (L == 255 || R == 255) ) {
			L = 255;     
		} else if (k == start){		
			tempind = k - ((temp >> 8) & 0xFF) - 1 ;	
			//printf("L %u R %u - BK x %d y %d - bnum %d - tempind %d - k %d - Pi %u\n", L, R, blockIdx.x, blockIdx.y, bnum, tempind, k, ((RLP[k] >> 8) & 0xFF));
			if (tempind == -1){
				sen_target_begin = 0;
			} else {
				sen_target_begin = RLP[tempind];
			}
			min_L = L;
			max_R = R;
		} else {
			if (min_L > L) {
				min_L = L;
			}
			if (max_R < R) {
				max_R = R;
			}
		} 
	}

	(*target_start) = min_L + sen_target_begin;
	(*target_end) = max_R + sen_target_begin;//Keep these guys

	if(frontOrEndWrong!=0){ //Bug Fix aXb
		return frontOrEndWrong;
	} else if (min_L <= max_R && max_R - min_L < MAX_rule_span){
		tempind++;

		if(consistent( (*target_start),  (*target_end), L_tar, R_tar, start, ender, 
		tempind)){
			return 1; //true
		}
	}			

	return 0;//normal false.

}

__device__ bool testUnAligned(int posA, uint8_t* L_target, uint8_t* R_target, int begin){
	if (posA >= begin && L_target[posA] == 255 && R_target[posA] == 255){
		return true;
	} 
	return false;
}

__global__ void extractConsistentPairs_OneGap(
		int reflen,
		uint8_t* L_target,
		uint8_t* R_target,
		unsigned int* RLP,
		unsigned int distinctOneGapSearchCount, 	
		int* counter_1gap,
		int* counter_2gap,
		bool looseH,
		int* refstr,
		rule_twogap* twoGapRule_ab_d,
		rule_onegap* oneGapRule_ab_d,
		unsigned int oneGapSACount,
		oneGapOnSA* oneGapSA,
		gappy_search* oneGapSearch,
		precomp_st_end* precomp_index,
		precompute_enu_3* precomp_onegap,
		bool sample ) {

	int oneBlockId = blockIdx.y + gridDim.y * blockIdx.x;
	if(oneBlockId >= distinctOneGapSearchCount){
		return;
	}

	int startSA = oneGapSearch[oneBlockId].start_on_salist;
	int endSA = oneGapSearch[oneBlockId].end_on_salist;

	if(startSA == -1 && endSA == -1){
		return;
	}

	if(startSA == -1 || endSA == -1){
		printf("How can this be possible only one minus 1?\n");
		return;
	}
	if (endSA >= oneGapSACount || endSA < startSA){
		printf("not possible1 in kernel one gap phrase extraction kernel\n");
		return;
	}

	int dis = 1 + endSA - startSA;

	uint8_t startLen = oneGapSearch[oneBlockId].qrystart_len;
	uint8_t endLen = oneGapSearch[oneBlockId].qryend_len;

	unsigned int gap1_start = 0;
	unsigned int gap1_end = 0;
	unsigned int target_start = 0;
	unsigned int target_end = 0;
	unsigned int gap2_start = 0;
	unsigned int gap2_end = 0;
	bool next = true;
	bool firstGap = true;
	int nextpos = 0;

	unsigned int current_str = 0;
	unsigned int ender = 0;
	unsigned char firstEnd = 0;
	uint8_t i = 1;
	bool left = true;
	bool right = true;
	bool preCompIndicator = false;
	unsigned int precomputationIndex = 0;
	uint8_t reNext = 1;
	
	//Precomputation
	if(dis == 1 && oneGapSA[startSA].length==0){
		preCompIndicator = true;
		precomputationIndex = oneGapSA[startSA].str_position;
		startSA = precomp_index[precomputationIndex].start;
		endSA = precomp_index[precomputationIndex].end;
		//Precomputation case
		dis = 1 + endSA - startSA;
		//If dis is <=0, then it must be filtered already in precomp with gap checking!!!!
		
		//debug purpose
		if(startLen!=1
				|| endLen!=1){
			printf("Good to know! One gap precomputation extraction kernel Wrong!\n");
			return;
		}
	} else if (oneGapSA[startSA].length==0) {
		printf("Thsi is not possible inside one gap extraction kernel!!\n");
		return;
	} 

	int threadx = threadIdx.x;

	unsigned char min_L;
	unsigned char max_R;
	int sen_target_begin;
	int tempind;	

	//Sampling variables
	bool pass = false;	
	if (!sample){
		pass = true;
	} else if (dis <= SAMPLER_ONEGAP) {
		pass = true;
	}
	
	bool rightThread = false;
	int desci = 0;
	float stepsize = (float)(dis) /(float)SAMPLER_ONEGAP;
	bool flager = true;
	int togo = -1;

	while(threadx < dis){		
		rightThread = false;
		desci = 0;
		flager = true;
		while (!pass && desci < SAMPLER_ONEGAP && flager) {
			togo = ROUND( desci* stepsize);    	
			if (togo == threadx) {
				rightThread = true;
				flager = false;
			} else if (togo > threadx){
				flager = false;
			}	
			desci++;
		}

		///Pass is for <Sampler, rightThread is for uniform selection, only right 
		// thread can come in and do matching.
		if(pass || rightThread){
			gap1_start = 0;
			gap1_end = 0;
			target_start = 0;
			target_end = 0;
			gap2_start = 0;
			gap2_end = 0;
			next = true;	
			firstGap = true;
			left = true;
			right = true;
				
			sen_target_begin = -1;
			tempind = -1;
			min_L = 255;
			max_R = 0;
			
			//Precomputation
			if(preCompIndicator){
				current_str = precomp_onegap[startSA+threadx].start;
				//last position, not the length; should be length - 1.
				firstEnd = precomp_onegap[startSA+threadx].length;
			} else{
				//debug purpose can be removed
				if(oneGapSA[startSA+threadx].position!=oneBlockId){
					printf("WONRG inside one gap phrase extraction kernel oneBlockId %d| position %d\n",oneBlockId,oneGapSA[startSA+threadx].position );
					return;
				}
				/////////////////////////////
				current_str = oneGapSA[startSA+threadx].str_position;
				firstEnd = oneGapSA[startSA+threadx].length;
			}
			//Check first gap
			//This checkBoundary can be changed to Sure version.
			//The first gap must be good.
			if(current_str + firstEnd - endLen > reflen){
				printf("Wrong! Precomp %d |Index %d |start %u|firstEnd %u|startLen %u|endLen %u|reflen %d|precomputationIndex %u|oneblock Id %d|threadx %d|startSA %d|endSA %d|dis %d\n",
					preCompIndicator, startSA+threadx, current_str, firstEnd, startLen, endLen, 
					reflen, precomputationIndex, oneBlockId,
					threadx, startSA, endSA, dis);
				return;
			}

			ender = current_str + firstEnd;
			
			firstGap = checkBoundaryFast(current_str + startLen, 
					ender - endLen,
					L_target, 
					R_target, 
					RLP, 
					&min_L, 
					&max_R,
					&sen_target_begin,
					&tempind);//main
					
			///Debug purpose
			if(!firstGap){
				printf("First Gap Wrong!! Precomp %d |Index %d |start %u|firstEnd %u|startLen %u|endLen %u|reflen %d|precomputationIndex %u|oneblock Id %d|threadx %d|startSA %d|endSA %d|dis %d\n",
					preCompIndicator, startSA+threadx, current_str, firstEnd, startLen, endLen, 
					reflen, precomputationIndex, oneBlockId,
					threadx, startSA, endSA, dis);
				return;
			}
			if(tempind == -1 || sen_target_begin == -1 || min_L > max_R){
				printf("Not initilized - first gap! tempind %d|sen_target_begin %d|min_L %d|max_R %d\n",
					tempind, sen_target_begin, min_L, max_R);
				return;
			}
			///Debug end

			gap1_start = min_L + sen_target_begin;
			gap1_end = max_R + sen_target_begin;
			
			////////////////Check aXb
			//Check the whole one gappy thing aXb
			if(firstGap){			
				/*Return error codes:
				0: normal false - aXb false only.
				1: true and ok
				2: front error - aXbX false, aXb false
				3: end error - XaXb false, aXb false
				4: both front and end error - XaXb false, aXb false, aXbX false
				*/
				reNext = checkBoundary(current_str, 
						ender,
						L_target, 
						R_target, 
						RLP, 
						&target_start, 
						&target_end);//main 		
				
				min_L = target_start - sen_target_begin;
				max_R = target_end - sen_target_begin;

				/*if(oneBlockId == 548){
					printf("aXb check| i %d, threadx %d, startSA %d, target_start %d, target_end %d, sen_target_begin %d|reNext %d|min_L %d max_R %d|gap1_start %d gap1_end %d||current_str %d firstEnd %d\n",
						i, threadx, startSA, target_start, target_end, sen_target_begin, reNext, min_L, max_R, gap1_start, gap1_end, current_str, firstEnd);
				}*/
				
				if(reNext == 0){
					next = false;
				} else if (reNext == 1){
					next = true;
				} else if (reNext == 2){
					next = false;
					right = false;
				} else if (reNext == 3){
					next = false;
					left = false;
				} else if (reNext == 4){
					next = false;
					left = false;
					right = false;
				}

				//Debug
				if((target_start == 0&& target_end==0) ||(min_L > max_R) || gap1_start < target_start || gap1_end > target_end){
					printf("Not initilized aXb? One gap kernel| blockID %d|target_start %d target_end %d|gap1_start %d gap1_end %d\n", 
						oneBlockId, target_start, target_end, gap1_start, gap1_end);
					return;
				}
			}
			
			if(next&&firstGap){
				nextpos = atomicAdd(counter_1gap, 1);
				oneGapRule_ab_d[nextpos].ref_str_start = target_start;
				oneGapRule_ab_d[nextpos].end = target_end -  target_start;
				oneGapRule_ab_d[nextpos].gap1 = gap1_start - target_start;
				oneGapRule_ab_d[nextpos].gap1_1= gap1_end - target_start;
				oneGapRule_ab_d[nextpos].gappy_index = oneBlockId;
				if(oneBlockId >= distinctOneGapSearchCount || gap1_end > target_end 
					|| gap1_start < target_start){
					printf("Thsi is not possible inside ExtractOneGap\n");
					return;
				}
			}

			////////////////Check two gap phrase: XaXb, aXbX
			unsigned int originalGapStart;
			unsigned int originalGapEnd;
			unsigned int temp;
			
			unsigned char L = 0;
			unsigned char R = 0;

			unsigned char min_XaXb = 255;
			unsigned char max_XaXb = 0;

			unsigned char min_aXbX = 255;
			unsigned char max_aXbX = 0;

			if(firstGap&&startLen + endLen + 1 + 1<= MAX_rule_symbols){
				target_start = 0;
				target_end = 0;
				i = 1;			

				originalGapStart = gap1_start;
				originalGapEnd = gap1_end;
				gap1_start = 0;
				gap1_end = 0;

				while(firstEnd + 1 + i <= MAX_rule_span
						&& (left || right)){
					//check XaXb
					if (left&&(int)(current_str - i) >=0 && refstr[current_str - i] >= 2 ){
						target_start = 0;
						target_end = 0;
						gap1_start = 0;
						gap1_end = 0;
						next = true;		

						//Check left gap X
						temp = RLP[current_str - i];
						L = (temp >> 24) & 0xFF;
						R = (temp >> 16) & 0xFF;
						
						if (L == 255 || R == 255){
							next = false;
							if (i == 1){
								left = false;
							}
						} else {
							if (min_XaXb > L) {
								min_XaXb = L;
							}
							if (max_XaXb < R) {
								max_XaXb = R;
							}
						} 
						
						///debugging purpose
						if (next && (min_XaXb > max_XaXb)){
							printf("Not possible - debuging one gap extraction \n");
							return;
						}
						//end debugging
						
						if (max_XaXb - min_XaXb >= MAX_rule_span){
							next = false;
							left = false;
						}
						
						//Check first gap
						if (next){
							gap1_start = sen_target_begin + min_XaXb;
							gap1_end = sen_target_begin + max_XaXb;
											
							next = consistent(gap1_start, 
											gap1_end, 
											L_target, 
											R_target, 
											current_str - i, 
											current_str - 1, 
											tempind);
						}			
						
						//Check the whole gappy XaXb thing
						if(next){ 
							if((min_XaXb < min_L)){
								target_start = sen_target_begin + min_XaXb;
							} else {
								target_start = sen_target_begin + min_L;	
							}
											
							if((max_XaXb < max_R)){
								target_end = sen_target_begin + max_R;
							} else {
								target_end = sen_target_begin + max_XaXb;	
							}
											
							///Debugging
							if (next&&(target_start > target_end)){
								printf("Not possible! Xab target_start %d | end %d | sen_target_begin %d | (min_L_Xab %d | min_L %d) | (max_R_Xab %d | max_R %d)\n",
										target_start, target_end, sen_target_begin, min_XaXb, min_L, max_XaXb, max_R);
								return;
							}
							//End of debugging
						
							if (target_end- target_start >= MAX_rule_span){
								next = false;
								left = false;
							}
							
							/*if(oneBlockId == 79){
								printf("aXbX||i %d -> next %d | right %d | target_start %d target_end %d | startSA %d EndSA %d | threadx %d dis %d|sen_tagret_begin %d | max_R %d, min_L %d|max_XaXb %d, min_XaXb %d\n", 
									i, next, right, target_start, target_end, startSA, endSA, threadx, dis,
									sen_target_begin, max_R, min_L, max_XaXb, min_XaXb);
							}*/
							
							if (next){
								next = consistent(target_start, 
													target_end, 
													L_target, 
													R_target, 
													current_str - i, 
													ender, 
													tempind);	
							}
						}
						
						/*if(oneBlockId == 79){
							printf("XaXb||i %d -> next %d | right %d | target_start %d target_end %d | originalGapStart %d originalGapEnd %d | gap2_start %d gap2_end %d | startSA %d EndSA %d | current_str %d ender %d | threadx %d dis %d\n", 
								i, next, right,
								target_start, target_end,
								originalGapStart, originalGapEnd,
								gap1_start, gap1_end,
								startSA, endSA,
								current_str, ender, threadx, dis);
						}*/

						if(next){
							nextpos = atomicAdd(counter_2gap, 1);
							twoGapRule_ab_d[nextpos].ref_str_start = target_start;
							twoGapRule_ab_d[nextpos].end = target_end -  target_start;//to the end
							twoGapRule_ab_d[nextpos].gap1 = gap1_start - target_start;
							twoGapRule_ab_d[nextpos].gap1_1= gap1_end - target_start;
							twoGapRule_ab_d[nextpos].gap2 = originalGapStart - target_start;
							twoGapRule_ab_d[nextpos].gap2_1 = originalGapEnd - target_start;
							twoGapRule_ab_d[nextpos].twogappyindex = 
								(oneBlockId);
							//Note the same here!
							left = false;
						}
					}else{
						left = false;
					}

					/////check aXbX
					if (right&&refstr[ender + i] >= 2 ){					
						target_start = 0;
						target_end = 0;
						next = true;
						gap2_start = 0;
						gap2_end = 0;


						//Check right gap X
						temp = RLP[ender + i];
						L = (temp >> 24) & 0xFF;
						R = (temp >> 16) & 0xFF;
						
						if (L == 255 || R == 255){
							next = false;
							if (i == 1){
								right = false;
							}
						} else {
							if (min_aXbX > L) {
								min_aXbX = L;
							}
							if (max_aXbX < R) {
								max_aXbX = R;
							}
						} 
						
						///debugging purpose
						if (next && (min_aXbX > max_aXbX)){
							printf("Not possible - debuging in aXbX \n");
							return;
						}
						//end debugging
						
						if (max_aXbX - min_aXbX >= MAX_rule_span){
							next = false;
							right = false;
						}
						
						//Check second gap aXbX
						if (next){
							gap2_start = sen_target_begin + min_aXbX;
							gap2_end = sen_target_begin + max_aXbX;

							///Debug
							if ((gap2_start > gap2_end)){
								printf("Not possible!\n");
								return;
							}

							next = consistent(gap2_start, 
											gap2_end, 
											L_target, 
											R_target, 
											ender + 1, 
											ender + i, 
											tempind);
						}			
						
						//Check the whole gappy thing
						if(next){ 
							if((min_aXbX < min_L)){
								target_start = sen_target_begin + min_aXbX;
							} else {
								target_start = sen_target_begin + min_L;	
							}
												
							if((max_aXbX < max_R)){
								target_end = sen_target_begin + max_R;
							} else {
								target_end = sen_target_begin + max_aXbX;	
							}

							///Debugging
							if (next&&(target_start > target_end)){
								printf("Not possible! abX target_start %d | end %d | sen_target_begin %d | (min_L_Xab %d | min_L %d) | (max_R_Xab %d | max_R %d)\n",
											target_start, target_end, sen_target_begin, min_aXbX, min_L, max_aXbX, max_R);
								return;
							}
							//End of debugging
							
							if (target_end - target_start >= MAX_rule_span){
								next = false;
								right = false;
							}
						
							/*if(oneBlockId == 121){
								printf("aXbX||i %d -> next %d | right %d | target_start %d target_end %d | startSA %d EndSA %d | threadx %d dis %d|sen_tagret_begin %d | max_R %d, min_L %d|max_aXbX %d, min_aXbX %d\n", 
									i, next, right, target_start, target_end, startSA, endSA, threadx, dis,
									sen_target_begin, max_R, min_L, max_aXbX, min_aXbX);
							}*/
							
							if (next){
								next = consistent(target_start, 
									target_end, 
									L_target, 
									R_target, 
									current_str, 
									ender+i, 
									tempind);	
							}
						}

						if(next){
							nextpos = atomicAdd(counter_2gap, 1);
							twoGapRule_ab_d[nextpos].ref_str_start = target_start;
							twoGapRule_ab_d[nextpos].end = target_end -  target_start;//to the end
							twoGapRule_ab_d[nextpos].gap1 = originalGapStart - target_start;
							twoGapRule_ab_d[nextpos].gap1_1= originalGapEnd - target_start;
							twoGapRule_ab_d[nextpos].gap2 = gap2_start - target_start;
							twoGapRule_ab_d[nextpos].gap2_1 = gap2_end - target_start;
							twoGapRule_ab_d[nextpos].twogappyindex = distinctOneGapSearchCount
								+ oneBlockId;
							right = false;
						}		
					}else{
						right = false;
					}
					i++;
				}
				///End
			}
		}		
		threadx += blockDim.x;
	}

}

__global__ void extractConsistentPairs_TwoGap(
		int reflen,
		uint8_t* L_target,
		uint8_t* R_target,
		unsigned int* RLP,
		unsigned int distinctTwoGapSearchCount, 	
		int* counter,
		bool looseH,
		int* refstr,
		rule_twogap* twoGapRule_ab_d,
		two_gappy_search* twoGapSearch,
		unsigned int twoGapSACount,
		twoGapOnSA* twoGapSA,
		gappy_search* oneGapSearch,
		unsigned int distinctOneGapSearchCount,
		bool sample) {

	int twoBlockId = blockIdx.y + gridDim.y * blockIdx.x;
	if(twoBlockId >= distinctTwoGapSearchCount){
		return;
	}

	int startSA = twoGapSearch[twoBlockId].start_on_salist;
	int endSA = twoGapSearch[twoBlockId].end_on_salist;
	if(startSA == -1 && endSA == -1){
		return;
	}
	if (endSA >= twoGapSACount || endSA < startSA || startSA == -1 || endSA == -1){
		printf("not possible1 in kernel two gap phrase extraction kernel|endSa %d | twoGapSACount %d | startSA %d|twoBlockId %d\n",
				endSA, twoGapSACount, startSA, twoBlockId);
		return;
	}

	int threadx = threadIdx.x;
	int dis = endSA - startSA+1;
	int oneBlockId = twoGapSearch[twoBlockId].blockid;

	uint8_t startLen = oneGapSearch[oneBlockId].qrystart_len;
	uint8_t endLen = oneGapSearch[oneBlockId].qryend_len;

	unsigned int gap1_start = 0;
	unsigned int gap1_end = 0;
	unsigned int target_start = 0;
	unsigned int target_end = 0;
	unsigned int gap2_start = 0;
	unsigned int gap2_end = 0;
	bool next = true;
	int nextpos = 0;

	unsigned int current_str = 0;
	unsigned int firstEnd = 0;
	unsigned int secondEnd = 0;
	uint8_t reNext = 1;

	//Sampling variables
	bool pass = false;	
	if (!sample){
		pass = true;
	} else if (dis <= SAMPLER_TWOGAP) {
		pass = true;
	}
	
	bool rightThread = false;
	int desci = 0;
	float stepsize = (float)(dis) /(float)SAMPLER_TWOGAP;
	bool flager = true;
	int togo = -1;

	while(threadx < dis){
		rightThread = false;
		desci = 0;
		flager = true;
		while (!pass && desci < SAMPLER_TWOGAP && flager) {
			togo = ROUND( desci* stepsize);    	
			if (togo == threadx) {
				rightThread = true;
				flager = false;
			} else if (togo > threadx){
				flager = false;
			}	
			desci++;
		}

		///Pass is for <Sampler, rightThread is for uniform selection, only right 
		// thread can come in and do matching.
		if(pass || rightThread){
			//debug purpose can be removed
			if(twoGapSA[startSA+threadx].position!=twoBlockId){
				printf("WONRG inside two gap phrase extraction kernel\n");
				return;
			}
			/////////////////////////////
			gap1_start = 0;
			gap1_end = 0;
			target_start = 0;
			target_end = 0;
			gap2_start = 0;
			gap2_end = 0;
			next = true;		
			current_str = twoGapSA[startSA+threadx].str_position;
			firstEnd = twoGapSA[startSA+threadx].length;
			secondEnd = twoGapSA[startSA+threadx].length2;

			//Check left gap
			next = checkBoundaryFast2(current_str + startLen, 
					current_str + firstEnd - endLen,
					L_target, 
					R_target, 
					RLP, 
					&gap1_start, 
					&gap1_end);//main

			//Check right gaps
			if(next){
				next = checkBoundaryFast2(current_str+firstEnd+1, 
						current_str+secondEnd-twoGapSearch[twoBlockId].qryend_len,
						L_target, 
						R_target, 
						RLP, 
						&gap2_start, 
						&gap2_end);//main
			}
			///Debug purpose
			if(!next){
				printf("This is not possible inside two gaps extraction kernel\n");
				return;
			}
			///Debug end

			//Check the whole gappy thing
			if(next){					
				reNext = checkBoundary(current_str, 
						current_str+secondEnd,
						L_target, 
						R_target, 
						RLP, 
						&target_start, 
						&target_end);//main 						
				if(reNext==1){
					next = true;
				} else {
					next = false;
				}
			}
			
			/*printf("2GapExt| OK %d | OneID %d | current_str %d | firstEnd %d | secondEnd %d | tar_start %d | tar_end %d | gap1_start %d | gap1_end %d | gap2_start %d | gap2_end %d | startSA %d | endSA %d\n",
				next, oneBlockId, current_str, firstEnd, secondEnd, target_start, target_end, 
				gap1_start, gap1_end, gap2_start, gap2_end, startSA, endSA);
			*/
			if(next){
				nextpos = atomicAdd(counter, 1);
				twoGapRule_ab_d[nextpos].ref_str_start = target_start;
				twoGapRule_ab_d[nextpos].end = target_end -  target_start;//to the end
				twoGapRule_ab_d[nextpos].gap1 = gap1_start - target_start;
				twoGapRule_ab_d[nextpos].gap1_1= gap1_end - target_start;
				twoGapRule_ab_d[nextpos].gap2 = gap2_start - target_start;
				twoGapRule_ab_d[nextpos].gap2_1 = gap2_end - target_start;
				twoGapRule_ab_d[nextpos].twogappyindex = twoBlockId;
			}	
		}
		threadx += blockDim.x;
	}
}

__global__ void extractConsistentPairs_Gappy(
		int*refsa,
		int reflen,
		uint8_t* L_target,
		uint8_t* R_target,
		unsigned int* RLP,
		saind_t* block,
		int globalc, 	
		int start_block,
		res_phrase_t* out_pair,
		int* counter,
		bool looseH,
		int* refstr,
		rule_onegap* oneGapRule_ab_d,
		rule_twogap* twoGapRule_ab_d,
		int* count_1gap_d,
		int* count_2gap_d,
		bool sample) {

	int bnum = blockIdx.y + gridDim.y * blockIdx.x +  start_block;

	if(bnum >= globalc){
		return;
	}

	int start = block[bnum].start;
	int end = block[bnum].end;
	if( threadIdx.x+start > end ){
		return;
	}
	int longestmatch = block[bnum].matchlen;
	if (longestmatch <1){
		return;
	}

	int current = threadIdx.x+start ;
	int k = 0;
	int current_str = 0;
	unsigned char L = 0;
	unsigned char R = 0;
	int sen_target_begin = -1;
	unsigned char min_L = 255;
	unsigned char max_R = 0;
	int tempind = 0;
	unsigned int temp;
	int nextpos = 0;

	//Used inside
	unsigned char i = 1;		
	unsigned int gap1_start = 0;
	unsigned int gap1_end = 0;	
	unsigned int gap2_start = 0;
	unsigned int gap2_end = 0;

	unsigned int target_start = 0;
	unsigned int target_end = 0;		
	bool next = true;
	int ender = -1;

	bool abX = true;
	bool Xab = true;
	bool XabX = true;
	bool ab = true;
	bool XabNoSuccess = true;//not equal to Xab, good so nolonger need checks 
	bool abXNoSuccess = true;
	
	uint8_t XabCount = 0;
	uint8_t abXCount = 0;
	
	unsigned char min_L_Xab = 255;
	unsigned char max_R_Xab = 0;		
	unsigned char min_L_abX = 255;
	unsigned char max_R_abX = 0;		

	unsigned char min_L_XabX = 255;
	unsigned char max_R_XabX = 0;

	///Sampling variable	
	bool pass = false;

	if (!sample) {
		pass = true;
	} else if (1 + end - start <= SAMPLER) {
		pass = true;
	}
	
	bool rightThread = false;
	int desci = 0;
	float stepsize = (float)(1 + end - start) /(float)SAMPLER;
	bool flager = true;
	int togo = -1;

	while(current <= end){      
		rightThread = false;
		desci = 0;
		flager = true;
		while (!pass && desci < SAMPLER && flager) {
			togo = ROUND( desci* stepsize);    	
			if (togo == current - start ) {
				rightThread = true;
				flager = false;
			} else if (togo > current - start ){
				flager = false;
			}	
			desci++;
		}	    
		///Pass is for <Sampler, rightThread is for uniform selection, only right 
		// thread can come in and do matching.
		if (pass || rightThread ) {				
			tempind = 0;
			min_L = 255;
			max_R = 0;
			current_str = refsa[current];

			abX = true;
			Xab = true;
			XabX = true;		
			ab = true;
			XabNoSuccess = true;
			abXNoSuccess = true;
			XabCount = 0;
			abXCount = 0;

			for(k = current_str; k < current_str + longestmatch; k++){
				temp = RLP[k];
				L = (temp >> 24) & 0xFF;
				R = (temp >> 16) & 0xFF;

				if(k == current_str){
					tempind = k - ((temp >> 8) & 0xFF) - 1 ;	
					//printf("L %u R %u - BK x %d y %d - bnum %d - tempind %d - k %d - Pi %u\n", L, R, blockIdx.x, blockIdx.y, bnum, tempind, k, ((RLP[k] >> 8) & 0xFF));
					if (tempind == -1){
						sen_target_begin = 0;
					} else {
						sen_target_begin = RLP[tempind];
					}
				}

				if ( (L == 255 || R == 255) 
						&& (k == current_str || k == current_str + longestmatch-1 ) ){
					//k = current_str + longestmatch + 1; //continue, donot stop
					ab = false;
					if(k==current_str){
						abXNoSuccess = false;
					} else {
						XabNoSuccess = false;
					}
				} else if ( (L == 255 || R == 255) ) {
					L = 255;     
				} else {
					if (min_L > L) {
						min_L = L;
					}
					if (max_R < R) {
						max_R = R;
					}
				} 
			}

			/*if(bnum==10){
						printf("XabX right before ab process 1|abX %d\n", abX);																						
			}*/
					
			if (min_L > max_R || max_R - min_L >= MAX_rule_span){
				abX = false;
				Xab = false;
				XabX = false;		
				ab = false;
				//continue;
			}
			/*
			if(bnum==10){
						printf("XabX before ab process 2|abX %d|min_L %d max_R %d\n", abX, min_L, max_R);																						
			}*/
					
			tempind++;
			ender = current_str+longestmatch - 1;
			if (ab){			
				/*int ss = min_L + sen_target_begin;
				  int tt = max_R + sen_target_begin;
				  int t = current_str + longestmatch - 1;*/
				nextpos = 0;
				//printf("target start %d - tar end %d - sour start %d - sour end %d\n", ss, tt, current_str, t);
				if (consistent(
							min_L + sen_target_begin, 
							max_R + sen_target_begin, 
							L_target, 
							R_target, 
							current_str, 
							ender, 
							tempind)){

					//Store Results!! ss tt current_str t	
					nextpos = atomicAdd(counter, 1);
					//out_pair[nextpos].start = current_str;
					out_pair[nextpos].tar_start = min_L + sen_target_begin;
					out_pair[nextpos].tar_end = max_R - min_L;
					out_pair[nextpos].blocknumber = bnum;				
				}				
			}		

			///ONE AND TWO GAP
			///////////////////////////////////////////////
			///Starting adding places from left and right!
			/////////////////////////////////////////////// 
			if (longestmatch + 1 > MAX_rule_symbols){
				abX = false;
				Xab = false;
			}
			if(longestmatch + 2 > MAX_rule_symbols){
				XabX = false;
			}

			i = 1;		
			min_L_Xab = 255;
			max_R_Xab = 0;		
			min_L_abX = 255;
			max_R_abX = 0;		

			min_L_XabX = 255;
			max_R_XabX = 0;
			/*if(bnum==10){
						printf("XabX before left/right process|abX %d\n", abX);																						
			}*/
					
			while(longestmatch + i <= MAX_rule_span 				
					&& (abXNoSuccess || XabNoSuccess || XabX)){
				if (Xab && current_str - i >=0 
						&& refstr[current_str - i] >= 2 ){
					next = true;				
					//Check left gap X
					temp = RLP[current_str - i];
					L = (temp >> 24) & 0xFF;
					R = (temp >> 16) & 0xFF;

					if (L == 255 || R == 255){
						next = false;
						if (i == 1){
							Xab = false;
							XabX = false;
						}
					} else {
						if (min_L_Xab > L) {
							min_L_Xab = L;
						}
						if (max_R_Xab < R) {
							max_R_Xab = R;
						}
					} 

					///debugging purpose
					if (next && (min_L_Xab > max_R_Xab)){
						printf("Not possible - debuging \n");
						//|| max_R - min_L >= LONGEST)){
						/*abX = false;
						  Xab = false;
						  XabX = false;		
						  ab = false;*/
						return;
					}
					//end debugging

					if (max_R_Xab - min_L_Xab >= MAX_rule_span){
						next = false;
						Xab = false;
					}
					
					//Check first gap
					if (next){
						gap1_start = sen_target_begin + min_L_Xab;
						gap1_end = sen_target_begin + max_R_Xab;
						if ((gap1_start > gap1_end)){
								printf("Not possible!\n");
								return;
						}
						
						next = consistent(gap1_start, 
								gap1_end, 
								L_target, 
								R_target, 
								current_str - i, 
								current_str - 1, 
								tempind);

						//Update the Xab count, the best shot on X side.
						if(next){
							XabCount = i;
						}
					}			

					//Check the whole gappy thing
					if(XabNoSuccess&&next){	
						if((min_L_Xab < min_L)){
							target_start = sen_target_begin + min_L_Xab;
						} else {
							target_start = sen_target_begin + min_L;	
						}
						
						if((max_R_Xab < max_R)){
							target_end = sen_target_begin + max_R;
						} else {
							target_end = sen_target_begin + max_R_Xab;	
						}
						
						///Debugging
						if (next&&(target_start > target_end)){
							printf("Not possible! Xab target_start %d | end %d | sen_target_begin %d | (min_L_Xab %d | min_L %d) | (max_R_Xab %d | max_R %d)\n",
									target_start, target_end, sen_target_begin, min_L_Xab, min_L, max_R_Xab, max_R);
							return;
						}
						//End of debugging

						if (target_end- target_start >= MAX_rule_span){
							next = false;
							Xab = false;
						}

						if (next){
							next = consistent(target_start, 
									target_end, 
									L_target, 
									R_target, 
									current_str - i, 
									ender, 
									tempind);	
						}
					}

					if(XabNoSuccess&&next){
						nextpos = atomicAdd(count_1gap_d, 1);						
						if(nextpos > PREALLOCATION_ONEGAP_RULE - 1000){
							printf("Error: Prealocated space for one gap is not enough\n");
							return;
						}
						oneGapRule_ab_d[nextpos].ref_str_start = target_start;
						oneGapRule_ab_d[nextpos].end = target_end -  target_start;
						oneGapRule_ab_d[nextpos].gap1 = gap1_start - target_start;
						oneGapRule_ab_d[nextpos].gap1_1= gap1_end - target_start;
						oneGapRule_ab_d[nextpos].gappy_index = bnum;									
						XabNoSuccess = false;
					}
					}else{
						Xab = false;
					}

					////////////////////////////
					///////////Adding from right
					////////////////////////////
					if (abX&&refstr[ender + i] >= 2){
						next = true;				
						//Check left gap X
						temp = RLP[ender + i];
						L = (temp >> 24) & 0xFF;
						R = (temp >> 16) & 0xFF;

						if (L == 255 || R == 255){
							next = false;
							if (i == 1){
								abX = false;
								XabX = false;
							}
						} else {
							if (min_L_abX > L) {
								min_L_abX = L;
							}
							if (max_R_abX < R) {
								max_R_abX = R;
							}
						} 

						///debugging purpose
						if (next && (min_L_abX > max_R_abX)){
							printf("Not possible - debuging in abX \n");
							/*abX = false;
							  Xab = false;
							  XabX = false;		
							  ab = false;*/
							return;
						}
						//end debugging
						if (max_R_abX - min_L_abX >= MAX_rule_span){
							next = false;
							abX = false;
						}
						//Check second gap abX
						if (next){
							gap1_start = sen_target_begin + min_L_abX;
							gap1_end = sen_target_begin + max_R_abX;
							if ((gap1_start > gap1_end)){
								printf("Not possible!\n");
								return;
							}
							next = consistent(gap1_start, 
									gap1_end, 
									L_target, 
									R_target, 
									ender + 1, 
									ender + i, 
									tempind);
							//Update the abX count, the best shot on X side.
							if(next){
								abXCount = i;
							}
						}			

						//Check the whole gappy thing
						if(abXNoSuccess&&next){	
							if((min_L_abX < min_L)){
								target_start = sen_target_begin + min_L_abX;
							} else {
								target_start = sen_target_begin + min_L;	
							}
							
							if((max_R_abX < max_R)){
								target_end = sen_target_begin + max_R;
							} else {
								target_end = sen_target_begin + max_R_abX;	
							}

							///Debugging
							if (next&&(target_start > target_end)){
								printf("Not possible! abX target_start %d | end %d | sen_target_begin %d | (min_L_Xab %d | min_L %d) | (max_R_Xab %d | max_R %d)\n",
									target_start, target_end, sen_target_begin, min_L_Xab, min_L, max_R_Xab, max_R);
								return;
							}
							//End of debugging

							if (target_end - target_start >= MAX_rule_span){
								next = false;
								abX = false;
							}

							if (next){
								next = consistent(target_start, 
										target_end, 
										L_target, 
										R_target, 
										current_str, 
										ender+i, 
										tempind);	
							}
						}

						if(abXNoSuccess&&next){
							nextpos = atomicAdd(count_1gap_d, 1);
							oneGapRule_ab_d[nextpos].ref_str_start = target_start;
							oneGapRule_ab_d[nextpos].end = target_end -  target_start;
							oneGapRule_ab_d[nextpos].gap1 = gap1_start - target_start;
							oneGapRule_ab_d[nextpos].gap1_1= gap1_end - target_start;
							oneGapRule_ab_d[nextpos].gappy_index = globalc + bnum;
							abXNoSuccess = false;
						}
					}else{
						abX = false;
					}

					///////////////////////////////////////////////
					///////////////////////////////////////////////
					///////////////XabX Two Gap////////////////////
					if(XabX&&(abX||Xab)){
					
						if(XabCount == i){
							min_L_XabX = 255;
							max_R_XabX = 0;
							for(uint8_t icount = 1; XabX&&icount <=abXCount; icount++){
								next = true;
								if(icount + XabCount + longestmatch <= MAX_rule_span){
									temp = RLP[ender + icount];
									L = (temp >> 24) & 0xFF;
									R = (temp >> 16) & 0xFF;
									
									if (L == 255 || R == 255){
										next = false;
										if (i == 1){
											printf("No possible here!!XabX should be checked 1|XabCount %d|abXcount %d\n", 
											XabCount, abXCount);									
											return;
										}
									} else {
										if (min_L_XabX > L) {
											min_L_XabX = L;
										}
										if (max_R_XabX < R) {
											max_R_XabX = R;
										}
									}
								}else {
									next = false;
									icount = abXCount+1;
								}
								
								if (next&&max_R_XabX - min_L_XabX >= MAX_rule_span){
									//printf("Not possible to reach this 1!XabX Should be checked already max_R_XabX %d | min_L_XabX %d\n", 
									//max_R_XabX, min_L_XabX);								
									next = false;
									icount = abXCount+1;
								}
								
								//Check XabX's second X; first Gap is good already.
								if(next){
									gap2_start = sen_target_begin + min_L_XabX;
									gap2_end = sen_target_begin + max_R_XabX;
									if ((min_L_XabX > max_R_XabX)){
										printf("Not possible! XabX min>max left X fix. Min %d max %d\n", min_L_XabX, max_R_XabX);
										return;
									}
									
									next = consistent(gap2_start, 
										gap2_end, 
										L_target, 
										R_target, 
										ender + 1, 
										ender + icount, 
										tempind);
								}
								//Debug
								/*if(bnum==10){
										printf("XabX check right gap: %d, min_L_XabX %d, max_R_XabX %d, sen_target_begin %d|next %d|XabCount %d\n",
											i, min_L_XabX, max_R_XabX, sen_target_begin, next, XabCount);																						
								}*/
									
								//Check the whole XabX thing
								if (next){			
									//Compare X1, ab, X2 min and max.
									if(min_L_XabX < min_L_Xab){
										temp = min_L_XabX;
									} else {
										temp = min_L_Xab;
									}
									if(temp > min_L){
										temp = min_L;
									} 
									target_start = sen_target_begin+temp;

									if(max_R_XabX < max_R_Xab){
										temp = max_R_Xab;
									} else {
										temp = max_R_XabX;
									}
									if(temp < max_R){
										temp = max_R;
									} 
									target_end = sen_target_begin+temp;
									if(target_start > target_end){
										printf("How can thsi be possible Xab\n");
										return;
									}
									if(target_end - target_start >= MAX_rule_span){
										next = false;
										icount = abXCount+1;
									}

									if(next){
										next = consistent(target_start, 
											target_end, 
											L_target, 
											R_target, 
											current_str - XabCount, 
											ender + icount, 
											tempind);
									}
									
									
									//Debug
									/*if(bnum==10){
										printf("XabX check right gap final one: %d, min_L_Xab %d, max_R_Xab %d, sen_target_begin %d, gap1_start %d, gap1_end %d, gap2_start %d, gap2_end %d|next %d|XabCount %d\n",
											i, min_L_Xab, max_R_Xab, sen_target_begin, gap1_start, gap1_end, gap2_start, gap2_end, next, XabCount);																						
									}*/
									
									if(next){
										gap1_start = sen_target_begin + min_L_Xab;
										gap1_end = sen_target_begin + max_R_Xab;

										nextpos = atomicAdd(count_2gap_d, 1);
										if(nextpos > PREALLOCATION_TWOGAP_RULE - 1000){
											printf("Error: Prealocated space for two gap is not enough\n");
											return;
										}
										twoGapRule_ab_d[nextpos].ref_str_start = target_start;
										twoGapRule_ab_d[nextpos].end = target_end -  target_start;
										twoGapRule_ab_d[nextpos].gap1 = gap1_start - target_start;
										twoGapRule_ab_d[nextpos].gap1_1= gap1_end - target_start;
										twoGapRule_ab_d[nextpos].gap2 = gap2_start - target_start;
										twoGapRule_ab_d[nextpos].gap2_1 = gap2_end - target_start;
										twoGapRule_ab_d[nextpos].twogappyindex = bnum;
										XabX = false;
									}
								}
							}
						}

						//Check Xab left side gap X.
						if (XabX&&abXCount == i){
							min_L_XabX = 255;
							max_R_XabX = 0;
							for(uint8_t icount = 1; XabX&&icount <=XabCount; icount++){
								next = true;
								if(icount + abXCount + longestmatch <= MAX_rule_span){
									//Check left gap X
									temp = RLP[current_str - icount];
									L = (temp >> 24) & 0xFF;
									R = (temp >> 16) & 0xFF;

									if (L == 255 || R == 255){
										next = false;
										if (i == 1){
											printf("No possible here!!XabX should be checked 2|XabCount %d|abXcount %d\n", 
											XabCount, abXCount);
											return;
										}
									} else {
										if (min_L_XabX > L) {
											min_L_XabX = L;
										}
										if (max_R_XabX < R) {
											max_R_XabX = R;
										}
									} 
								}else {
									icount = XabCount+1;
									next = false;			
								}

								if (next&&max_R_XabX - min_L_XabX >= MAX_rule_span){
									//Stop this for loop
									//printf("Not possible to reach this 2!XabX Should be checked already max_R_XabX %d | min_L_XabX %d\n", 
									//	max_R_XabX, min_L_XabX);								
									icount = XabCount+1;
									next = false;
								}
								
								//Check XabX's first X; right Gap is good already.
								if(next){
									gap1_start = sen_target_begin + min_L_XabX;
									gap1_end = sen_target_begin + max_R_XabX;
									if ((min_L_XabX > max_R_XabX)){
										printf("Not possible! XabX min>max right X fix. Min %d max %d\n", min_L_XabX, max_R_XabX);
										//printf("Not possible! XabX min>max\n");
										return;
									}
									
									next = consistent(gap1_start, 
										gap1_end, 
										L_target, 
										R_target, 
										current_str -icount,
										current_str -1, 
										tempind);
								}
								
								//Debug
								/*if(bnum==10){
										printf("XabX check left gap: %d, min_L_XabX %d, max_R_XabX %d, sen_target_begin %d|next %d|abXCount %d|icount %d|XabCount %d\n",
											i, min_L_XabX, max_R_XabX, sen_target_begin, next, abXCount, icount, XabCount);																						
								}*/
								//Check the whole XabX thing
								if (next){			
									//Compare X1, ab, X2 min and max.
									if(min_L_XabX < min_L_abX){
										temp = min_L_XabX;
									} else {
										temp = min_L_abX;
									}
									if(temp > min_L){
										temp = min_L;
									} 
									target_start = sen_target_begin+temp;

									if(max_R_XabX < max_R_abX){
										temp = max_R_abX;
									} else {
										temp = max_R_XabX;
									}
									if(temp < max_R){
										temp = max_R;
									} 
									target_end = sen_target_begin+temp;
									if(target_start > target_end){
										printf("How can thsi be possible Xab\n");
										return;
									}
									
									if(target_end - target_start >= MAX_rule_span){
										next = false;
										icount = XabCount+1;
									}
									if(next){
										next = consistent(target_start, 
											target_end, 
											L_target, 
											R_target, 
											current_str - icount, 
											ender + abXCount, 
											tempind);
									}
									
									//Debug
									/*if(bnum==10){
										printf("XabX check left gap final one: %d, min_L_abX %d, max_R_abX %d, sen_target_begin %d, gap1_start %d, gap1_end %d|next %d|abXCount %d|target_start %d target_end %d|min_L %d, max_R %d|icount %d\n",
											i, min_L_abX, max_R_abX, sen_target_begin, gap1_start, gap1_end, next, abXCount,
											target_start, target_end, min_L, max_R, icount);		
									}*/
									
									if(next){
										gap2_start = sen_target_begin + min_L_abX;
										gap2_end = sen_target_begin + max_R_abX;
										
										nextpos = atomicAdd(count_2gap_d, 1);
										twoGapRule_ab_d[nextpos].ref_str_start = target_start;
										twoGapRule_ab_d[nextpos].end = target_end -  target_start;
										twoGapRule_ab_d[nextpos].gap1 = gap1_start - target_start;
										twoGapRule_ab_d[nextpos].gap1_1= gap1_end - target_start;
										twoGapRule_ab_d[nextpos].gap2 = gap2_start - target_start;
										twoGapRule_ab_d[nextpos].gap2_1 = gap2_end - target_start;
										twoGapRule_ab_d[nextpos].twogappyindex = bnum;
										XabX = false;
									}
								}
							}						
						}					
					} else {
						XabX = false;		
					}

					
					//////////////////////////////////////////////
					//Sync to prevent spin
					if(XabX == false){
						if(Xab == false && XabNoSuccess == true){
							XabNoSuccess = false;
						}
						if(abX == false && abXNoSuccess == true){
							abXNoSuccess = false;
						}
					}
					i++;
				}
		}
		current += blockDim.x;
		}   		
	}

	__global__ void extractConsistentPairs(
			int*refsa,
			int reflen,
			uint8_t* L_target,
			uint8_t* R_target,
			unsigned int* RLP,
			saind_t* block,
			int globalc, 	
			int start_block,
			res_phrase_t* out_pair,
			int* counter,
			bool looseH) {

		int bnum = blockIdx.y + gridDim.y * blockIdx.x +  start_block;

		if(bnum >= globalc){
			return;
		}

		int start = block[bnum].start;
		int end = block[bnum].end;
		if( threadIdx.x+start > end ){
			return;
		}
		int longestmatch = block[bnum].matchlen;
		if (longestmatch <1){
			return;
		}

		int current = threadIdx.x+start ;
		int k = 0;
		int current_str = 0;
		unsigned char L = 0;
		unsigned char R = 0;
		int sen_target_begin = -1;
		unsigned char min_L = 255;
		unsigned char max_R = 0;
		int tempind = 0;
		unsigned int temp;
		int nextpos = 0;
		int ok = 1;


		while(current <= end){                   	  
			ok = 1;
			tempind = 0;
			min_L = 255;
			max_R = 0;
			current_str = refsa[current];
			for(k = current_str; k < current_str + longestmatch; k++){
				temp = RLP[k];
				L = (temp >> 24) & 0xFF;
				R = (temp >> 16) & 0xFF;
				if ( (L == 255 || R == 255) && (k == current_str || k == current_str + longestmatch-1 ) ){
					k = current_str + longestmatch + 1;
					ok = 0;
				} else if ( (L == 255 || R == 255) ) {
					L = 255;     
				} else if (k == current_str){		
					tempind = k - ((temp >> 8) & 0xFF) - 1 ;	
					//printf("L %u R %u - BK x %d y %d - bnum %d - tempind %d - k %d - Pi %u\n", L, R, blockIdx.x, blockIdx.y, bnum, tempind, k, ((RLP[k] >> 8) & 0xFF));
					if (tempind == -1){
						sen_target_begin = 0;
					} else {
						sen_target_begin = RLP[tempind];
					}
					min_L = L;
					max_R = R;
				} else {
					if (min_L > L) {
						min_L = L;
					}
					if (max_R < R) {
						max_R = R;
					}
				} 
			}

			if (ok == 1 && min_L <= max_R && max_R - min_L < LONGEST ){
				tempind++;
				int ss = min_L + sen_target_begin;
				int tt = max_R + sen_target_begin;
				int t = current_str + longestmatch - 1;
				nextpos = 0;
				//printf("target start %d - tar end %d - sour start %d - sour end %d\n", ss, tt, current_str, t);
				if (consistent(ss, tt, L_target, R_target, current_str, t, tempind)){
					//Store Results!! ss tt current_str t	
					/*if (*counter > ALLOCATIONS - 100){
					  return;
					  }*/
					nextpos = atomicAdd(counter, 1);
					//printf("C %d |", nextpos );
					/*if(bnum == 15){
					  printf("BK x %d y %d - bnum %d - ThreadX %d - sour s%d - sour t%d - tar s%d - tar len%d - count %d - npos %d\n", blockIdx.x, blockIdx.y, bnum, threadIdx.x, current_str, current_str+longestmatch, ss, max_R - min_L, *counter, nextpos);
					  }*/
					//out_pair[nextpos].start = current_str;
					out_pair[nextpos].tar_start = ss;
					out_pair[nextpos].tar_end = max_R - min_L;
					out_pair[nextpos].blocknumber = bnum;

					////////LOOSE HEURISTIC
					if(looseH){
						//LEFT CHECK BOUNDARY
						int left_boundary = 1;
						while(testUnAligned(ss - left_boundary, L_target, R_target, sen_target_begin)){
							left_boundary++ ;
						}
						left_boundary--;
						///Right CHECK BOUNDAY
						int right_boundary = 1;
						while(testUnAligned(tt + right_boundary, L_target, R_target,sen_target_begin)){
							right_boundary++ ;
						}
						right_boundary--;

						for(int i = left_boundary; i >= 0; i --){
							for(int j = right_boundary; j >= 0; j --){
								if(i != 0 || j != 0 ){
									nextpos = atomicAdd(counter, 1);
									out_pair[nextpos].tar_start = ss - left_boundary;
									out_pair[nextpos].tar_end = max_R - min_L + right_boundary;
									out_pair[nextpos].blocknumber = bnum;
								}
							}
						}
					}
				}				
			}			
			current += blockDim.x;
		}   		
	}

	__global__ void extractConsistentPairsSample(
			int*refsa,
			int reflen,
			uint8_t* L_target,
			uint8_t* R_target,
			unsigned int* RLP,
			saind_t* block,
			int globalc, 	
			res_phrase_t* out_pair,
			int* counter,
			int looseH) {

		int bnum = blockIdx.y + gridDim.y * blockIdx.x;

		if(bnum >= globalc){
			return;
		}

		int start = block[bnum].start;
		int end = block[bnum].end;
		if( threadIdx.x+start > end ){
			return;
		}
		int longestmatch = block[bnum].matchlen;
		if (longestmatch <1){
			return;
		}

		int current = threadIdx.x+start ;
		int k = 0;
		int current_str = 0;
		unsigned int L = 0;
		unsigned int R = 0;
		int sen_target_begin = -1;
		unsigned int min_L = 255;
		unsigned int max_R = 0;
		int tempind = 0;
		unsigned int temp;
		int nextpos = 0;
		int ok = 1;
		bool pass = false;

		if (1 + end - start <= SAMPLER) {
			pass = true;
		}

		bool rightThread = false;
		int desci = 0;
		float stepsize = (float)(1 + end - start) /(float)SAMPLER;
		bool flager = true;
		int togo = -1;

		while(current <= end){       
			rightThread = false;
			desci = 0;
			flager = true;
			while (desci < SAMPLER && flager && !pass) {
				togo = ROUND( desci* stepsize);    	
				if (togo == current - start ) {
					rightThread = true;
					flager = false;
				} else if (togo > current - start ){
					flager = false;
				}	
				desci++;
			}	    
			///Pass is for <Sampler, rightThread is for uniform selection, only right 
			// thread can come in and do matching.
			if (pass || rightThread ) {
				ok = 1;
				tempind = 0;
				min_L = 255;
				max_R = 0;
				current_str = refsa[current];
				for(k = current_str; k < current_str + longestmatch; k++){
					temp = RLP[k];
					L = (temp >> 24) & 0xFF;
					R = (temp >> 16) & 0xFF;
					if ( (L == 255 || R == 255) && (k == current_str || k == current_str + longestmatch-1 ) ){
						k = current_str + longestmatch + 1;
						ok = 0;
					} else if ( (L == 255 || R == 255) ) {
						L = 255;
					} else if (k == current_str){		
						tempind = k - ((temp >> 8) & 0xFF) - 1 ;	
						//printf("L %u R %u - BK x %d y %d - bnum %d - tempind %d - k %d - Pi %u\n", L, R, blockIdx.x, blockIdx.y, bnum, tempind, k, ((RLP[k] >> 8) & 0xFF));
						if (tempind == -1){
							sen_target_begin = 0;
						} else {
							sen_target_begin = RLP[tempind];
						}
						min_L = L;
						max_R = R;
					} else {
						if (min_L > L) {
							min_L = L;
						}
						if (max_R < R) {
							max_R = R;
						}
					} 
				}

				if (ok == 1 && min_L <= max_R  && max_R - min_L < LONGEST ){
					tempind++;
					int ss = min_L + sen_target_begin;
					int tt = max_R + sen_target_begin;
					int t = current_str + longestmatch - 1;
					nextpos = 0;
					//printf("target start %d - tar end %d - sour start %d - sour end %d\n", ss, tt, current_str, t);
					if (consistent(ss, tt, L_target, R_target, current_str, t, tempind)){
						//Store Results!! ss tt current_str t	
						nextpos = atomicAdd(counter, 1);
						//printf("BK x %d y %d - bnum %d - ThreadX %d - sour s%d - sour t%d - tar s%d - tar t%d - count %d\n", blockIdx.x, blockIdx.y, bnum, threadIdx.x, current_str, current_str+longestmatch, ss, ss+max_R - min_L, counter);
						//out_pair[nextpos].start = current_str;
						//out_pair[nextpos].end = longestmatch - 1;
						out_pair[nextpos].tar_start = ss;
						out_pair[nextpos].tar_end = max_R - min_L;
						out_pair[nextpos].blocknumber = bnum;

						////////LOOSE HEURISTIC
						if(looseH){
							//LEFT CHECK BOUNDARY
							int left_boundary = 1;
							while(testUnAligned(ss - left_boundary, L_target, R_target, sen_target_begin)){
								left_boundary++ ;
							}
							left_boundary--;
							///Right CHECK BOUNDAY
							int right_boundary = 1;
							while(testUnAligned(tt + right_boundary, L_target, R_target,sen_target_begin)){
								right_boundary++ ;
							}
							right_boundary--;

							for(int i = left_boundary; i >= 0; i --){
								for(int j = right_boundary; j >= 0; j --){
									if(i != 0 || j != 0 ){
										nextpos = atomicAdd(counter, 1);
										out_pair[nextpos].tar_start = ss - left_boundary;
										out_pair[nextpos].tar_end = max_R - min_L + right_boundary;
										out_pair[nextpos].blocknumber = bnum;
									}
								}
							}
						}
					}				
				}			
			}
			current += blockDim.x;
		}   		
	}

	__global__ void extractGlobalPairsUpDown(
			int global,
			red_dup_t* fast_speed,
			int count,
			result_t* globalOnPairsUpDown_d) {
		int globalindex = threadIdx.x + blockDim.x * blockIdx.x;
		if(globalindex >= global){
			return;
		}
		bool goUp = blockIdx.y%2;

		saind_t output = binary_search_hit(fast_speed, 0, count-1, globalindex); //-1?
		if (output.matchlen == -1 ){
			globalOnPairsUpDown_d[globalindex].up = -1;
			globalOnPairsUpDown_d[globalindex].down = -1;
			return;
		}
		if (goUp){
			int up = binary_search_up(fast_speed, output.matchlen, output.end, globalindex);
			globalOnPairsUpDown_d[globalindex].up = up;
		} else {
			int down = binary_search_down(fast_speed, output.start, output.matchlen, globalindex);
			globalOnPairsUpDown_d[globalindex].down = down;
		}	
	}

	__device__ float searchLexFile(
			int sourceId, 
			int targetId, 
			wordKey* lexFileKeys,
			categ* lexFileValues,
			unsigned int lexFileCount,
			bool oneOrTwo){
	
		unsigned int low = 0;
		unsigned int high = lexFileCount;
		unsigned int middle;
		float output = 0.0;
		bool stop = false;
		while (low <= high && !stop) {
			middle = low + (high - low)/2;
			if (sourceId < lexFileKeys[middle].ch){
				high = middle - 1;		
			} else if (sourceId > lexFileKeys[middle].ch){				
				low = middle + 1;		
			} else if (targetId < lexFileKeys[middle].eng) {
				high = middle - 1;
			} else if (targetId > lexFileKeys[middle].eng) {				
				low = middle + 1;				
			} else {
				if(oneOrTwo){
					output = lexFileValues[middle].val1;
				} else {
					output = lexFileValues[middle].val2;
				}
				stop = true;
			}
		}
	
		return output;
	}	
	
	__global__ void lexicalTaskMaxEF(
			int oneGapBound, //one gap lexicon boundary
			int twoGapBound, //two gap lexicon boundary
			int lexicalTaskCounter, //Total number of lexicon
			lexicalTask* maxEF_d, //all tasks are here: MaxEF
			int* targetStr, //target string
			unsigned int targetToklen,//target string boundary checking
			wordKey* lexFileKeys,//search among those keys
			categ* lexFileValues,//retrieve the values and get max
			unsigned int lexFileCount,//search range
			lexTaskResults* resultLex_d){
	
		unsigned int blockId = gridDim.x* blockIdx.y + blockIdx.x;
		unsigned int lexId = blockId*blockDim.x+threadIdx.x;
	
		if (lexId >= lexicalTaskCounter){
			return;
		}
	
		int j = 0;
		float fgivene = 0;
		float egivenf = 0;
		int jj = 0;
		float max_scorer_1 = 0.0;
	
		int sourcePatternCounter = maxEF_d[lexId].sourcePatternCounter;
		bool firstIn = true;
		float resTemp = 0;
		if (lexId < oneGapBound){ //one gap
			int targetEnd = maxEF_d[lexId].targetStart+maxEF_d[lexId].end;
			int gap1Start = maxEF_d[lexId].targetStart+maxEF_d[lexId].gap1;
			int gap1End = maxEF_d[lexId].targetStart+maxEF_d[lexId].gap1_1;
			
			///MaxLexFGivenE
			for(j = 0; j < sourcePatternCounter; j++){
				max_scorer_1 = 0;
				firstIn = true;
				for(jj = maxEF_d[lexId].targetStart; jj <= targetEnd; jj++){
					if (jj < gap1Start || jj > gap1End){
						if (firstIn) {
							resTemp = searchLexFile(
								maxEF_d[lexId].sourcePattern[j], 
								-1, 
								lexFileKeys,
								lexFileValues,
								lexFileCount,
								false); 
							if(resTemp > max_scorer_1){
								max_scorer_1 = resTemp;
							}
							firstIn = false;
						}
						
						resTemp = searchLexFile(
								maxEF_d[lexId].sourcePattern[j], 
								targetStr[jj], 
								lexFileKeys,
								lexFileValues,
								lexFileCount,
								false); 
						if(resTemp > max_scorer_1){
							max_scorer_1 = resTemp;
						}					
					}
				}
				if (max_scorer_1 > 0){
					fgivene += -log10(max_scorer_1);
				} else {
					fgivene += MAXSCORE;
				}
			}
	
			///MaxLexEGivenF
			for(jj = maxEF_d[lexId].targetStart; jj <= targetEnd; jj++){
				if (jj < gap1Start || jj > gap1End){								
					max_scorer_1 = 0;
					firstIn = true;
					for(j = 0; j < sourcePatternCounter; j++){
						if (firstIn) {
							resTemp = searchLexFile(
								-1,
								targetStr[jj], 
								lexFileKeys,
								lexFileValues,
								lexFileCount,
								true);	
							if(resTemp > max_scorer_1){
								max_scorer_1 = resTemp;
							}
							firstIn = false;
						}
						
						resTemp = searchLexFile(
								maxEF_d[lexId].sourcePattern[j], 
								targetStr[jj], 
								lexFileKeys,
								lexFileValues,
								lexFileCount,
								true);	
						if(resTemp > max_scorer_1){
							max_scorer_1 = resTemp;
						}
					}
					if (max_scorer_1 > 0){
						egivenf += -log10(max_scorer_1);
					} else {
						egivenf += MAXSCORE;
					}
				}
			}
			
		} else if (lexId < twoGapBound){ //two gap
			int targetEnd = maxEF_d[lexId].targetStart+maxEF_d[lexId].end;
			int gap1Start = maxEF_d[lexId].targetStart+maxEF_d[lexId].gap1;
			int gap1End = maxEF_d[lexId].targetStart+maxEF_d[lexId].gap1_1;
			int gap2Start = maxEF_d[lexId].targetStart+maxEF_d[lexId].gap2;
			int gap2End = maxEF_d[lexId].targetStart+maxEF_d[lexId].gap2_1;
			
			///MaxLexFGivenE
			for(j = 0; j < sourcePatternCounter; j++){
				max_scorer_1 = 0;
				firstIn = true;
				for(jj = maxEF_d[lexId].targetStart; jj <= targetEnd; jj++){
					if ((jj < gap1Start || jj > gap1End) 
						&& (jj < gap2Start || jj > gap2End) ){
						if (firstIn) {
							resTemp = searchLexFile(
								maxEF_d[lexId].sourcePattern[j], 
								-1, 
								lexFileKeys,
								lexFileValues,
								lexFileCount,
								false); 
							if(resTemp > max_scorer_1){
								max_scorer_1 = resTemp;
							}
							firstIn = false;
						}
						
						resTemp = searchLexFile(
								maxEF_d[lexId].sourcePattern[j], 
								targetStr[jj], 
								lexFileKeys,
								lexFileValues,
								lexFileCount,
								false); 
						if(resTemp > max_scorer_1){
							max_scorer_1 = resTemp;
						}					
					}
				}
				if (max_scorer_1 > 0){
					fgivene += -log10(max_scorer_1);
				} else {
					fgivene += MAXSCORE;
				}
			}
	
			///MaxLexEGivenF
			for(jj = maxEF_d[lexId].targetStart; jj <= targetEnd; jj++){
				if ((jj < gap1Start || jj > gap1End) 
					&& (jj < gap2Start || jj > gap2End) ){
					max_scorer_1 = 0;
					firstIn = true;
					for(j = 0; j < sourcePatternCounter; j++){
						if (firstIn) {
							resTemp = searchLexFile(
								-1,
								targetStr[jj], 
								lexFileKeys,
								lexFileValues,
								lexFileCount,
								true);	
							if(resTemp > max_scorer_1){
								max_scorer_1 = resTemp;
							}
							firstIn = false;
						}
						
						resTemp = searchLexFile(
								maxEF_d[lexId].sourcePattern[j], 
								targetStr[jj], 
								lexFileKeys,
								lexFileValues,
								lexFileCount,
								true);	
						if(resTemp > max_scorer_1){
							max_scorer_1 = resTemp;
						}
					}
					if (max_scorer_1 > 0){
						egivenf += -log10(max_scorer_1);
					} else {
						egivenf += MAXSCORE;
					}
				}
			}
		} else { //continous		
			int targetEnd = maxEF_d[lexId].targetStart+maxEF_d[lexId].end;
			
			///MaxLexFGivenE
			for(j = 0; j < sourcePatternCounter; j++){
				max_scorer_1 = 0;
				firstIn = true;
				for(jj = maxEF_d[lexId].targetStart; jj <= targetEnd; jj++){
					if (firstIn) {
						resTemp = searchLexFile(
								maxEF_d[lexId].sourcePattern[j], 
								-1, 
								lexFileKeys,
								lexFileValues,
								lexFileCount,
								false); 
						if(resTemp > max_scorer_1){
							max_scorer_1 = resTemp;
						}
						/*printf("lexId %d | -1 | Score: %f, MaxScore %f\n", 
							lexId,
							resTemp,
							max_scorer_1);*/
						firstIn = false;
					}
						
					resTemp = searchLexFile(
								maxEF_d[lexId].sourcePattern[j], 
								targetStr[jj], 
								lexFileKeys,
								lexFileValues,
								lexFileCount,
								false); 
					if(resTemp > max_scorer_1){
						max_scorer_1 = resTemp;
					}					
					/*printf("lexId %d | Score: %f, MaxScore %f\n", 
							lexId,
							resTemp,
							max_scorer_1);*/
				}
				if (max_scorer_1 > 0){
					fgivene += -log10(max_scorer_1);
				} else {
					fgivene += MAXSCORE;
				}
				/*printf("lexId %d | fgivene %f\n",
					lexId,
					fgivene);*/
			}
	
			///MaxLexEGivenF
			for(jj = maxEF_d[lexId].targetStart; jj <= targetEnd; jj++){
					max_scorer_1 = 0;
					firstIn = true;
					for(j = 0; j < sourcePatternCounter; j++){
						if (firstIn) {
							resTemp = searchLexFile(
								-1,
								targetStr[jj], 
								lexFileKeys,
								lexFileValues,
								lexFileCount,
								true);	
							if(resTemp > max_scorer_1){
								max_scorer_1 = resTemp;
							}
							firstIn = false;
						}
						
						resTemp = searchLexFile(
								maxEF_d[lexId].sourcePattern[j], 
								targetStr[jj], 
								lexFileKeys,
								lexFileValues,
								lexFileCount,
								true);	
						if(resTemp > max_scorer_1){
							max_scorer_1 = resTemp;
						}
					}
					if (max_scorer_1 > 0){
						egivenf += -log10(max_scorer_1);
					} else {
						egivenf += MAXSCORE;
					}
			}
		}
		
		resultLex_d[lexId].MaxLexFgivenE = fgivene;
		resultLex_d[lexId].MaxLexEgivenF = egivenf; 
	}

	void extractPairFinalize(void *handler) {
		context_t * ctx = (context_t *)handler;
		cudaFree(ctx->out_pair);
		cudaFree(ctx->blocks_d);
		cudaFree(ctx->ref_target_d.L_tar);
		cudaFree(ctx->ref_target_d.R_tar);
	}

	lexicalFileCuda* initWordPossibilityIntKey(disk_handler_t *handler, ref_set_t* refset){		
		mytimer_t t;
		timer_start(&t);

		hashtbl* refstrid = refset->refs[0].users;
		hashtbl* targetstrid = refset->refs_target[0].users_target;
		
		fprintf(stderr, "Start Word Possibility as INT KEY, loading CHHASH %d ENHASH %d\n", HASH_COUNT(refstrid), HASH_COUNT(targetstrid));
		lexicalFileCuda* returnVal = (lexicalFileCuda*)malloc(sizeof(lexicalFileCuda));

		//make sure hard code
		category* valuesOnCPU = (category*)malloc(LEXFILE_LINE_PRE*sizeof(category));		
		wordKey* keysOnCPU = (wordKey*)malloc(LEXFILE_LINE_PRE*sizeof(wordKey));
		
		int counter = 0;
		ifstream file(handler->wordpossibility);
		if (!file.is_open()){
			fprintf(stderr, "The Word Possibility File is not Found!\n");
			exit(0);
		}

		while (file.good()) {
			string chinese, english;
			categ val;
			file>> chinese;
			file>> english;
			file>> val.val1;
			file>> val.val2;

			//cerr << "file get "<< chinese<< " "<<english<<" "<<val1<<" "<<val2<< endl;
			struct my_struct *sss;
			HASH_FIND_STR(refstrid, chinese.c_str(), sss);            
			if(!sss){                                				
				if (chinese != "NULL"){					
					cout<< "Ch Not Available!!! " << chinese <<endl; 
					continue;
				}
			} 
			//cerr<<"chinese finished ";
			struct my_struct *sss1;		   
			HASH_FIND_STR(targetstrid, english.c_str(), sss1); 		   
			if(!sss1){ 						
				if (english != "NULL"){					
					cout<< "En Not Available!!! " << english <<endl; 
					continue;
				}
			} 
			//cerr<<"English finished "<<endl;
			wordKey key;
			if (!sss && !sss1){//For 'NULL' case that appear inside bitexts, deal problems.
				key.ch = -1;
				key.eng = -1;
			} else if (sss && !sss1) {
				key.ch = sss->id;
				key.eng = -1;
				/*char buffer[200];
				sprintf (buffer, "%d|-1", sss->id);
				key = string(buffer);*/
			} else if (!sss && sss1) {
				key.ch = -1;
				key.eng = sss1->id;
				/*char buffer[200];
				sprintf (buffer, "-1|%d", sss1->id);
				key = string(buffer);*/
			} else {
				key.ch = sss->id;
				key.eng = sss1->id;
				/*char buffer[200];
				sprintf(buffer, "%d|%d", sss->id, sss1->id);
				key = string(buffer);*/
			}
			//cerr<<"get "<<key<<endl;
			
			//Push
			valuesOnCPU[counter]= val;
			keysOnCPU[counter]=key;
			counter++;
		}	

		/*map<string, categ>::iterator itt;
		  for(itt = output.begin(); itt != output.end(); ++itt){
		  cerr << (*itt).first << ": " << (*itt).second.val1<< " "<< (*itt).second.val2  << endl;
		  }*/		  
		file.close();
		fprintf(stderr, "Lex File Word Possibility COUNTER: %d\n", counter);
		//////////////Now deal with CUDA side and sort them
		cudaMalloc((void**)&(returnVal->values), counter*sizeof(category));
		cudaMalloc((void**)&(returnVal->keys), counter*sizeof(wordKey));
		cudaMemcpy(returnVal->values, valuesOnCPU, counter*sizeof(category), cudaMemcpyHostToDevice);
		cudaMemcpy(returnVal->keys, keysOnCPU, counter*sizeof(wordKey), cudaMemcpyHostToDevice);

		thrust::device_ptr<categ> valuesThrust = thrust::device_pointer_cast(returnVal->values);
		thrust::device_ptr<workKeytyp> keysThrust = thrust::device_pointer_cast(returnVal->keys);
		thrust::device_ptr<workKeytyp> keysThrustEnd = thrust::device_pointer_cast(returnVal->keys+counter);
		//sort by key! with thrust
		thrust::sort_by_key(keysThrust, keysThrustEnd, valuesThrust, lexFileCompare());

		///Debug		
		/*workKeytyp* debugKeys = (workKeytyp*)malloc(counter*sizeof(workKeytyp));
		cudaMemcpy(debugKeys, returnVal->keys, counter*sizeof(wordKey), cudaMemcpyDeviceToHost);
		for(int i =0; i< counter; i++){
			printf("CUDALEX|| i %d -> ch %d | eng %d\n", i, debugKeys[i].ch,debugKeys[i].eng);
		}*/
		
		/////////////
		timer_stop(&t);
		fprintf(stderr, "Loading Lex File time: %f\n", timer_elapsed(&t)/1000);
		free(valuesOnCPU);
		free(keysOnCPU);
		returnVal->count = counter;
		
		return returnVal;
	}

	map<string, categ> initWordPossibility(disk_handler_t *handler,
			hashtbl* refstrid, 
			hashtbl* targetstrid){
		fprintf(stderr, "Start word possibility loading CHHASH %d ENHASH %d\n", 
			HASH_COUNT(refstrid), HASH_COUNT(targetstrid));
		map<string, categ> output;
		int counter = 0;
		ifstream file(handler->wordpossibility);
		if (!file.is_open()){
			fprintf(stderr, "The Word Possibility File is not Found!\n");
			exit(0);
		}

		while (file.good()) {
			string chinese, english;
			float val1, val2;
			file>> chinese;
			file>> english;
			file>> val1;
			file>> val2;

			//cerr << "file get "<< chinese<< " "<<english<<" "<<val1<<" "<<val2<< endl;
			char* chC = new char[chinese.size()+1];
			char* enC = new char[english.size()+1];
			strcpy(chC, chinese.c_str());
			strcpy(enC, english.c_str());

			struct my_struct *sss;
			HASH_FIND_STR(refstrid, chC, sss);            
			if(!sss){                                				
				if (chinese != "NULL"){					
					cout<< "Ch Not Available!!! " << chinese << " al " <<chC <<endl; 
					continue;
				}
			} 
			//cerr<<"chinese finished ";
			struct my_struct *sss1;		   
			HASH_FIND_STR(targetstrid, enC, sss1); 		   
			if(!sss1){ 						
				if (english != "NULL"){					
					cout<< "En Not Available!!! " << english << " al "<<enC <<endl; 
					continue;
				}
			} 
			//cerr<<"English finished "<<endl;
			string key = "";
			if (!sss && !sss1){
				key = "-1|-1";
			} else if (sss && !sss1) {
				char buffer[200];
				sprintf (buffer, "%d|-1", sss->id);
				key = string(buffer);
			} else if (!sss && sss1) {            
				char buffer[200];
				sprintf (buffer, "-1|%d", sss1->id);
				key = string(buffer);
			} else {
				char buffer[200];
				sprintf(buffer, "%d|%d", sss->id, sss1->id);
				key = string(buffer);
			}
			//cerr<<"get "<<key<<endl;
			categ val;
			val.val1 = val1;
			val.val2 = val2;
			output[key] = val;
			counter++;
		}	

		if (output.size() != counter){
			fprintf(stderr, "Word Possibility File COUNTER != Map Size\n");		
			exit(0);
		}

		/*map<string, categ>::iterator itt;
		  for(itt = output.begin(); itt != output.end(); ++itt){
		  cerr << (*itt).first << ": " << (*itt).second.val1<< " "<< (*itt).second.val2  << endl;
		  }*/
		file.close();
		fprintf(stderr, "Word Possibility File COUNTER: %d\n", counter);
		return output;
	}

	void initAlignment(ref_set_t *refset, disk_handler_t *handler, uint *devicemem){
		//fprintf(stderr, "\nLoading the Alignment Reference\n");
		assert(refset->refs_target);
		assert(refset->refs);
		ref_t_target * ref_target = &refset->refs_target[0];			 
		ref_t * ref_source = &refset->refs[0];			 
		cudaMallocHost((void **)(&ref_target->L_tar), ref_target->toklen*sizeof(uint8_t));

		cudaMallocHost((void **)(&ref_target->R_tar), ref_target->toklen*sizeof(uint8_t));
		uint8_t* L_source = (uint8_t*) malloc(ref_source->toklen* sizeof(uint8_t)); 
		uint8_t* R_source = (uint8_t*) malloc(ref_source->toklen* sizeof(uint8_t)); 	
		memset(ref_target->L_tar, 255, ref_target->toklen*sizeof(uint8_t));
		memset(ref_target->R_tar, 255, ref_target->toklen*sizeof(uint8_t));
		memset(L_source, 255, ref_source->toklen*sizeof(uint8_t));
		memset(R_source, 255, ref_source->toklen*sizeof(uint8_t));

		int qcount = -1;    
		int read = 0;
		const char* delim = " -";    
		char* line = NULL;    
		size_t len = 0;

		char* token = NULL;

		while(read = getline(&line, &len, handler->align)!= -1){        
			//fprintf(stderr, "%d Line %s\n", qcount+1, line);
			qcount++;
			token = NULL;
			//components = NULL;
			size_t newbuflen = strlen(line);
			if (line[newbuflen - 1] == '\n') 
				line[newbuflen - 1] = '\0';
			token = strtok(line, delim);        
			while(token != NULL && !isspace(*token)){            
				int sour_no = atoi(token);			
				//fprintf(stderr, "%d -> ", sour_no);
				token = strtok(NULL, delim);
				if (token==NULL){
					printf("Not possible!\n");
					exit(0);
				}
				int tar_no = atoi(token);
				//fprintf(stderr, "%d | ", tar_no);
				//printf("Source %d Target %d\n", sour_no, tar_no);
				if(sour_no >= 255 || tar_no >= 255 || sour_no < 0 || tar_no <0){
					printf("Not possible, too long sentence\n");
					exit(1);
				}
				int inde = ref_source->sentenceind[qcount]+sour_no;
				//printf("Index %d\n", inde);
				if (L_source[inde] ==255 || R_source[inde] ==255){
					L_source[inde] = tar_no; //Min
					R_source[inde] = tar_no; //Max
				} else if (tar_no > R_source[inde] ){
					R_source[inde] = tar_no;
				} else if (tar_no < L_source[inde] ){
					L_source[inde] = tar_no;
				}

				if (ref_target->L_tar[ref_target->sentenceind[qcount] + tar_no] ==255 || ref_target->R_tar[ref_target->sentenceind[qcount] + tar_no] ==255){
					ref_target->L_tar[ref_target->sentenceind[qcount] + tar_no] = sour_no; //Min
					ref_target->R_tar[ref_target->sentenceind[qcount] + tar_no] = sour_no; //Max
				} else if (sour_no > ref_target->R_tar[ref_target->sentenceind[qcount] + tar_no] ){
					ref_target->R_tar[ref_target->sentenceind[qcount] + tar_no] = sour_no;
				} else if (sour_no < ref_target->L_tar[ref_target->sentenceind[qcount] + tar_no] ){
					ref_target->L_tar[ref_target->sentenceind[qcount] + tar_no] = sour_no;
				}
				if(L_source[inde] > R_source[inde]){
					fprintf(stderr, "Not possible!! L > R\n");
					exit(1);
				}
				token = strtok(NULL, delim); 			
			}
			//fprintf(stderr,"\n");

		}
		//fprintf(stderr, "Start RLP generation...\n");
		/////////////////////////////////////////RLP Generation
		cudaMallocHost((void **)(&ref_source->RLP), ref_source->toklen*sizeof(unsigned int));
		int i = 0;

		qcount = 1;
		for(; i < ref_source->toklen-1; i++){
			if(i == ref_source->sentenceind[qcount]-1){
				ref_source->RLP[i] = ref_target->sentenceind[qcount];
				//printf("i %d - qcount %d - RLP sentence index %d\n", i, qcount, ref_target->sentenceind[qcount]);
				qcount++;
			} else{
				unsigned char c = 0;
				ref_source->RLP[i] = ((L_source[i] << 24) | ( R_source[i] << 16) | (ref_source->P[i] << 8) | c);
				//printf("i %d - qcount %d - L %u - R %u - P %u\n", i, qcount, (L_source[i]) , ( R_source[i]) , (ref_source->P[i]));
			}
		}
		//fprintf(stderr, "Finish RLP generation... Generated qcount %d\n", qcount);
		free(L_source);
		free(R_source);
		cudaFreeHost(ref_source->P);
		free(ref_target->sentenceind);
		free(ref_source->sentenceind);
		////////////////////////////////////////	
	}

	//vector<vector<unsigned int> > GenerateBlocks(
	void GenerateBlocks(
			void *handler, 
			ref_set_t * refset, 
			qry_set_t * qryset, 
			saind_t * tmp_blocks, 
			unsigned int* global,   
			vector<char*> *sourceName,
			vector<vector<unsigned int> > &  qryglobal,
			char** vocabulary) {			
		mytimer_t t;
		//hashtbl * cutdup = NULL;
		map<string, int> checkD;

		int refindex, j=0;
		context_t * ctx = (context_t *) handler;

		fprintf(stderr, "Generating Blocks of Results\n");	
		timer_start(&t);

		int qryindex;
		fprintf(stderr, "QRYCOUNT is %d\n", qryset->qryscount);
		//vector<vector<unsigned int> >  qryglobal(qryset->qryscount);
		char str1[300];
		char str2[100];
		char str3[100];
		char sname[200];
		char str_source[450];
		long totalcheck = 0;
		map<int, int> removalDUP;	

		for (qryindex = 0; qryindex < qryset->qryscount; qryindex++) {
			vector<unsigned int> go;	  
			for (refindex = 0; refindex < refset->count; refindex++) {
				ref_t * ref = &refset->refs[refindex];

				result_t_two* qryresult = qryset->result_two;
				int end = -1;
				if (qryindex != qryset->qryscount-1){
					end = qryset->qrysoffsettok[qryindex+1];
				} else {
					end = qryset->totaltokens;
				}
				//fprintf(stderr, "Generating Blocks of Results %d - S%d E%d\n", qryindex,qryset->qrysoffsettok[qryindex], end);

				for (j = qryset->qrysoffsettok[qryindex]; j < end; j++) {
					//fprintf(stderr, "Generating %d - longest %d\n", j, qryresult[j].longestmatch);
					map<string, int>::iterator iter;
					if (qryresult[j].longestmatch > 0) { /// length is 1
						sprintf(str1,"%d",qryresult[j].up);
						sprintf(str2,"|%d|1",qryresult[j].down);
						strcat(str1,str2);
						string checking(str1);
						iter = checkD.find(checking);///////HERE!
						if(iter == checkD.end()){
							tmp_blocks[*global].start = qryresult[j].up;
							tmp_blocks[*global].end = qryresult[j].down;
							tmp_blocks[*global].string_start = ref->sa[qryresult[j].up];
							tmp_blocks[*global].matchlen = 1;
							////////////////////
							strcpy(str_source,"");
							for(int spane = qryresult[j].up; spane <  1+qryresult[j].up; spane++){
								if (spane == qryresult[j].up){
									sprintf(sname,"%s", vocabulary[ref->str[ref->sa[spane]]]);
								} else {
									sprintf(sname," %s", vocabulary[ref->str[ref->sa[spane]]]);
								}
								//fprintf(stderr, "\t\tPre Sname: %s\n", sname);
								strcat(str_source, sname);			
							}	

							char* str_ss = /*new char[strlen(str_source)+1]; */ (char*)malloc(sizeof(char)*(strlen(str_source)+1));
							copy(str_source, str_source+strlen(str_source)+1 , str_ss);
							(*sourceName)[*global] = str_ss;
							/////////////////
							go.push_back(*global);
							(*global)++;
							totalcheck += qryresult[j].down - qryresult[j].up;
							checkD.insert(make_pair(checking, (*global)-1));
							pair<map<int, int>::iterator, bool> res = removalDUP.insert(make_pair((*global)-1, -1));
						} else {
							pair<map<int, int>::iterator, bool> res = removalDUP.insert(make_pair(iter->second, -1));
							if ( res.second ) {//Not Inserted already
								go.push_back(iter->second);
								//utarray_push_back(go, &(iter->second));
							}                        
						}
					}
					if (qryresult[j].longestmatch > 1) {
						int cc=qryset->connectoffset[j];
						int ct;
						for(ct = 2; ct <= qryresult[j].longestmatch && ct <= LONGESTCHSOURCE; ct++, cc++){
							//fprintf(stderr, " Inside second longest %d - up %d - low %d", ct, qryset->result_connect[cc].up, qryset->result_connect[cc].down);

							sprintf(str1,"%d",qryset->result_connect[cc].up);
							sprintf(str2,"|%d",qryset->result_connect[cc].down);
							sprintf(str3,"|%d",ct);
							strcat(str1,str2);
							strcat(str1,str3);
							string checking(str1);
							iter = checkD.find(checking);
							//HASH_FIND_STR(cutdup, str1, sss);    


							if(iter == checkD.end()){
								//fprintf(stderr, "INSIDE second block 1\n");

								tmp_blocks[*global].start = qryset->result_connect[cc].up;
								tmp_blocks[*global].end = qryset->result_connect[cc].down;
								tmp_blocks[*global].matchlen = ct;                            
								tmp_blocks[*global].string_start = ref->sa[qryset->result_connect[cc].up];
								//fprintf(stderr, "INSIDE second block 2\n");
								///////////////////////////
								strcpy(str_source,"");
								int gogogo = ref->sa[qryset->result_connect[cc].up];
								//fprintf(stderr,"Yoyo shi 1.5999 Sure?\n");
								for(int spane = gogogo; spane < ct + gogogo ; spane++){

									if (spane == gogogo){
										sprintf(sname,"%s",vocabulary[ref->str[spane]]);
									} else {
										sprintf(sname," %s",vocabulary[ref->str[spane]]);
									}
									//fprintf(stderr, "\t\tPre Sname: %s\n", sname);
									strcat(str_source, sname);			
								} 

								char* str_ss = /*new char[strlen(str_source)+1]*/ (char*)malloc(sizeof(char)*(strlen(str_source)+1));
								copy(str_source, str_source+strlen(str_source)+1 , str_ss); //memcpy((void*)str_ss, str_source, strlen(str_source)+1);
								(*sourceName)[*global] = str_ss;
								///////////////////////////								
								//utarray_push_back(go, global);///HERE
								go.push_back(*global);
								(*global)++;
								totalcheck += qryset->result_connect[cc].down - qryset->result_connect[cc].up;
								checkD.insert(make_pair(checking, (*global)-1));
								pair<map<int, int>::iterator, bool> res = removalDUP.insert(make_pair((*global) - 1, -1));
							} else {
								pair<map<int, int>::iterator, bool> res = removalDUP.insert(make_pair(iter->second, -1));
								//fprintf(stderr, "INSIDE second block else part\n");
								if ( res.second ) {//Not Inserted already
									//fprintf(stderr, "INSIDE second block sels s\n");
									go.push_back(iter->second);
									//utarray_push_back(go, &(iter->second));
								}                        
							}
							//PRINT("||%d %d||%d\n", upper, lower, lower - upper + 1);
						}
					}
				}
			}
			qryglobal[qryindex] = go;
			removalDUP.clear();
		}
		if(TEMPSET*(qryset->totaltokens) < *global){
			fprintf(stderr, "ERROR! Memory Space Exceeded Inside GenerateBlocks! Fixable Problem.\n");
			exit(0);
		}
		checkD.clear();
		timer_stop(&t);	
		//fprintf(stderr, "generating blocks time: %f - total check %d\n", timer_elapsed(&t)/1000, totalcheck);
		//return qryglobal;
	}

	void ExtractPairs(void *handler,
			ref_set_t *refset, 
			qry_set_t *qryset, 
			timing_t *timing, 
			map<string, categ> word_score) {

		saind_t* tmp_blocks = (saind_t*)malloc(1.5*(qryset->totaltokens)*sizeof(saind_t));
		char** sourceNameBlock = 
			(char**)malloc(1.5*(qryset->totaltokens)*sizeof(char*));
		int global = 0;
		//UT_array ** qryglobal = GenerateBlocks(handler, refset, qryset, tmp_blocks, 
		//        &global, sourceNameBlock, intchar);
		mytimer_t __t;
		context_t * ctx = (context_t *) handler;
		ref_t_target * ref_target_d = &ctx->ref_target_d;			 
		ref_t * ref_source_d = &ctx->ref_d;			 
		qry_set_t *qryset_d = &ctx->qryset_d;
		ref_t_target * ref_target = &refset->refs_target[0];			 
		ref_t * ref_source = &refset->refs[0];			 
		int count = 0, *count_d;

		cudaMalloc((void**)&count_d, sizeof(int));
		cudaMemcpy(count_d, &count, sizeof(int), cudaMemcpyHostToDevice);	
		///L and R in target side - CUDA mem location and copy
		///RLP in source side - CUDA Mem and copy	
		fprintf(stderr, "ref_source_d->RLP %f MB, ref_target_d->L_tar %f MB, ref_target_d->R_tar %f MB, ctx->out_pair %f MB, Ref SA %f MB, Global %f MB\n", 
				(double)ref_source->toklen*sizeof(int)/(1024*1024), 
				(double)ref_target->toklen*sizeof(uint8_t)/(1024*1024), 
				(double)ref_target->toklen*sizeof(uint8_t)/(1024*1024), 
				(double)ALLOCATIONS*sizeof(res_phrase_t)/(1024*1024),
				(double)ref_source_d ->toklen*sizeof(int)/(1024*1024),
				(double)global*sizeof(saind_t)/(1024*1024));

		timer_start(&__t);
		cudaMalloc((void**)&(ref_source_d->RLP), ref_source->toklen*sizeof(int));
		cudaMalloc((void**)&(ref_target_d->L_tar), ref_target->toklen*sizeof(uint8_t));
		cudaMalloc((void**)&(ref_target_d->R_tar), ref_target->toklen*sizeof(uint8_t));
		cudaMalloc((void**)&(ctx->out_pair), ALLOCATIONS*sizeof(res_phrase_t));
		timing->extractin += timer_stop(&__t);

		fprintf(stderr, "Start Extraction Data Transfer\n");
		timer_start(&__t);
		cudaMemcpy(ref_target_d->R_tar, ref_target->R_tar, ref_target->toklen*sizeof(uint8_t), cudaMemcpyHostToDevice);
		cudaMemcpy(ref_target_d->L_tar, ref_target->L_tar, ref_target->toklen*sizeof(uint8_t), cudaMemcpyHostToDevice);
		cudaMemcpy(ref_source_d->RLP, ref_source->RLP, ref_source->toklen*sizeof(int), cudaMemcpyHostToDevice);
		timing->extractin += timer_stop(&__t);
		timer_start(&__t);

		dim3  block(THREADS_PER_BLOCK, 1);
		dim3  grid(LINEARBLOCK, (global + LINEARBLOCK - 1)/LINEARBLOCK);

		fprintf(stderr, "Start Kernel on Extration of Pairs - Global %d - Grid X %d | Y %d\n", global, grid.x, grid.y);
		assert(ref_source_d ->sa!=NULL
				&&ref_target_d->L_tar!=NULL
				&&ref_target_d->R_tar!=NULL
				&&ref_source_d->RLP!=NULL
				&&ctx->blocks_d!=NULL
				&&ctx->out_pair!=NULL);

		//bool isSample = false;//true;
		bool looseH = false;//Loose Heuristic Switch
		if(!isSample){
			extractConsistentPairs<<<grid, block>>>(
					ref_source_d ->sa,
					ref_source_d ->toklen,
					ref_target_d->L_tar,
					ref_target_d->R_tar,
					ref_source_d->RLP,
					ctx->blocks_d,
					global,
					0,
					ctx->out_pair,
					count_d,
					looseH);
		} else {
			extractConsistentPairsSample<<<grid, block>>>(
					ref_source_d ->sa,
					ref_source_d ->toklen,
					ref_target_d->L_tar,
					ref_target_d->R_tar,
					ref_source_d->RLP,
					ctx->blocks_d,
					global,
					ctx->out_pair,
					count_d,
					looseH);
		}
		cudaThreadSynchronize();
		timing->extractkernel += timer_stop(&__t);
		cudaError_t error = cudaGetLastError();
		fprintf(stderr, "Kernel Extraction time: %f\n", timer_elapsed(&__t)/1000);
		if(error != cudaSuccess)
		{
			// print the CUDA error message and exit if any
			fprintf(stderr, "CUDA error second After kernel 1 - old: %s \n", cudaGetErrorString(error));
			exit(-1);
		}	
		cudaMemcpy(&count, count_d, sizeof(int), cudaMemcpyDeviceToHost);
		fprintf(stderr, "Found %d pairs!!\n", count);

		cudaFreeHost(ref_target->L_tar);
		cudaFreeHost(ref_target->R_tar);
		cudaFreeHost(ref_source->RLP);
		cudaFreeHost(ref_source->sa);
		//cudaFreeHost(ref_source->rank);
		//cudaFreeHost(ref_target->sentenceind);
		//cudaFreeHost(ref_source->sentenceind);
		//return;///Here!
		res_phrase_t* out_res;
		cudaMallocHost((void **)&(out_res), count*sizeof(res_phrase_t));
		assert(out_res != NULL);
		///out_res and fast_speed could be large.
		cudaMemcpy(out_res, ctx->out_pair, count*sizeof(res_phrase_t), cudaMemcpyDeviceToHost);   

		hashtbl_aux* lexic = NULL;
		hash_lexicon* target_count = NULL;
		hash_lexicon* foreig_count = NULL;
		//hash_intchar* index_lexi_table = NULL;
		fprintf(stderr, "Start Creating Lexincon\n");

		int lexicon_count = 0;
		red_dup_t* fast_speed;/* = createLexicon(
					 count, 
					 out_res, 
					 ref_source ->str, 
					 ref_target->str, 
					 intchar, 
					 intchar_target,
					 &lexic, 
					 &target_count,
					 &foreig_count, 
					 &lexicon_count,
					 tmp_blocks,
					 word_score,
					 sourceNameBlock,
					 global);
				       */
		assert(lexicon_count != 0 && fast_speed != NULL);
		fprintf(stderr, "Lexicon count is %d\n", lexicon_count);
		qsort(fast_speed, lexicon_count, sizeof(red_dup_t), compareUser);

		fprintf(stderr, "Create Lexincon DONE\n");

		bool cpu = 0;
		if (cpu){
			fprintf(stderr, "Start Print Grammar\n");
			/*print_query_CPU(count, 
			  out_res, 
			  ref_source ->str, 
			  ref_target->str, 
			  intchar, 
			  intchar_target, 
			  lexic, 
			  target_count, 
			  foreig_count,
			  qryglobal, 
			  qryset->qryscount);*/
		} else {
			fprintf(stderr, "Start Extract Global Pairs - GPU Mode\n");
			result_t* globalOnPairsUpDown = (result_t*)malloc(global*sizeof(result_t));
			red_dup_t* fast_speed_d = NULL;

			cudaFree(ref_source_d ->sa);
			cudaFree(ref_source_d ->RLP);
			cudaFree(ref_target_d->L_tar);
			cudaFree(ref_target_d->R_tar);
			cudaFree(ctx->out_pair);
			cudaFreeHost(out_res);
			//cudaFree(ctx->blocks_d);

			cudaMalloc((void**)&(ctx->globalOnPairsUpDown_d), global*sizeof(result_t));		
			cudaMalloc((void**)&(fast_speed_d), lexicon_count*sizeof(red_dup_t));
			cudaMemcpy(fast_speed_d, fast_speed, lexicon_count*sizeof(red_dup_t), cudaMemcpyHostToDevice);  
			assert(ctx->globalOnPairsUpDown_d!=NULL && fast_speed_d!=NULL);

			///For each global, get its up and lower boundary
			timer_start(&__t);
			dim3  block1(THREADS_PER_BLOCK_GLOBALPAIRS, 1);
			dim3  grid1((global + THREADS_PER_BLOCK_GLOBALPAIRS - 1)/THREADS_PER_BLOCK_GLOBALPAIRS, 2);

			extractGlobalPairsUpDown<<<grid1, block1>>>(
					global,
					fast_speed_d,
					lexicon_count,
					ctx->globalOnPairsUpDown_d);

			cudaThreadSynchronize();
			cudaMemcpy(globalOnPairsUpDown, ctx->globalOnPairsUpDown_d, global*sizeof(result_t), cudaMemcpyDeviceToHost);
			timing->extractkernel += timer_stop(&__t);

			cudaError_t error1 = cudaGetLastError();
			if(error1 != cudaSuccess)
			{
				// print the CUDA error message and exit if any 
				fprintf(stderr, "CUDA error second After kernel 1 - old: %s \n", cudaGetErrorString(error1));
				exit(-1);
			}

			fprintf(stderr, "Start Print Grammar - GPU Mode\n");         

		}	
	}

	__global__ void extractConsistentPairsGappy_TwoGap_Precomp(
			int toklen,
			uint8_t* L_tar,
			uint8_t* R_tar,
			unsigned int* RLP,
			two_gappy_precomp* twogaps_precomp,
			unsigned int twogapscount_precomp,	
			precompute_enu_3* precomp_onegap,
			unsigned int precomp_count,		
			unsigned int* count_two_gap_d,
			rule_twogap* twoGapRule_d,
			bool loose){

		int index = blockDim.x*blockIdx.x+threadIdx.x;

		if(index >= twogapscount_precomp){
			return;
		}

		unsigned int within_each_precomp = twogaps_precomp[index].within_each_precomp;
		unsigned int current_str = precomp_onegap[within_each_precomp].start;
		unsigned int ender =  current_str + 
			precomp_onegap[within_each_precomp].length + 
			twogaps_precomp[index].gap;

		unsigned int target_start = 0;
		unsigned int target_end = 0;
		if(!checkBoundary(current_str, ender, L_tar, R_tar, RLP, &target_start, &target_end)){ //main
			return;
		}

		if(target_start==0&&target_end == 0 || target_start > target_end){
			printf("MUST BE WRON[G!!\n");
			return;
		}

		unsigned int gap1_start = 0;
		unsigned int gap1_end = 0;

		if(!checkBoundary(current_str+1, current_str + precomp_onegap[within_each_precomp].length-1, L_tar, R_tar, RLP, &gap1_start, &gap1_end)){ //gap 1
			return;
		}

		if(gap1_start==0&&gap1_end == 0 || gap1_start > gap1_end){
			printf("MUST BE WRON[G!! -Gap 1\n");
			return;
		}

		unsigned int gap2_start = 0;
		unsigned int gap2_end = 0;

		if(!checkBoundary(current_str + precomp_onegap[within_each_precomp].length+1, 	
					ender - twogaps_precomp[index].qryend_len, 
					L_tar, 
					R_tar, 
					RLP, 
					&gap2_start, 
					&gap2_end)){ //gap 2
			return;
		}

		if(gap2_start==0&&gap2_end == 0 || gap2_start > gap2_end){
			printf("MUST BE WRON[G!! -Gap 2\n");
			return;
		}

		int nextpos = atomicAdd(count_two_gap_d, 1);

		twoGapRule_d[nextpos].ref_str_start = target_start;
		twoGapRule_d[nextpos].twogappyindex = (-1)*index; //if negative then from precomp
		twoGapRule_d[nextpos].end = target_end - target_start;
		twoGapRule_d[nextpos].gap1 = gap1_start - target_start; // relative position
		twoGapRule_d[nextpos].gap1_1 = gap1_end - target_start; // relative position 	
		twoGapRule_d[nextpos].gap2 = gap2_start - target_start; // relative position	
		twoGapRule_d[nextpos].gap2_1 = gap2_end - target_start; // relative position
	}
	
struct twoGapResCompare
{
	__host__ __device__
		bool operator()(const rule_twogap& a, const rule_twogap& b) const
		{
			return (a.twogappyindex < b.twogappyindex);
		}

};

struct oneGapResCompare
{
	__host__ __device__
		bool operator()(const rule_onegap& a, const rule_onegap& b) const
		{
			return (a.gappy_index < b.gappy_index);
		}

};

struct continousResCompare
{
	__host__ __device__
		bool operator()(const res_phrase_s& a, const res_phrase_s& b) const
		{
			return (a.blocknumber < b.blocknumber);
		}

};

void ExtractPairs_Large_Data_Gappy(void *handler,
		ref_set_t *refset, 
		qry_set_t *qryset, 
		timing_t *timing, 		
        lexicalFileCuda* cudaLexFile,
		disk_handler_t *diskInfo) {

	mytimer_t __t;
	char* destDir = diskInfo->destinationDirectory;	
	context_t * ctx = (context_t *) handler;
	ref_t_target * ref_target_d = &ctx->ref_target_d;			 
	ref_t * ref_source_d = &ctx->ref_d; 		 
	qry_set_t *qryset_d = &ctx->qryset_d;
	ref_t_target * ref_target = &refset->refs_target[0];			 
	ref_t * ref_source = &refset->refs[0];	

	size_t freeMem = 0;
	size_t totalMem = 0;
	size_t allocMem = 0;

	cudaMemGetInfo(&freeMem, &totalMem);  
	fprintf(stderr, "Memory avaliable: Free: %lu, Total: %lu - Start of Extract.cu\n",freeMem/(1024*1024), totalMem/(1024*1024));

	saind_t* tmp_blocks = (saind_t*)malloc(sizeof(saind_t)*TEMPSET*(qryset->totaltokens));
	unsigned int global = 0;    
	vector<char*> sourceNameBlock(TEMPSET*(qryset->totaltokens));// = (char**)malloc(2*(qryset->totaltokens)*sizeof(char*));
	
	vector<vector<unsigned int> > qryglobal(qryset->qryscount);
	GenerateBlocks(handler, 
					refset, 
					qryset, 
					tmp_blocks, 
					&global, 
					&sourceNameBlock, 
					qryglobal,
					ref_source->vocabulary);
	if(qryglobal.size()==0 || global > TEMPSET*qryset->totaltokens || global <= 0){
		fprintf(stderr, "Function generateblocks() Failed Due to Physical Memory Limitations\n");
		return;
	}
	
	cudaMalloc((void**)&(ctx->blocks_d), (global)*sizeof(saind_t));
	cudaMemcpy(ctx->blocks_d, tmp_blocks, (global)*sizeof(saind_t), cudaMemcpyHostToDevice);
	
	///L and R in target side - CUDA mem location and copy
	///RLP in source side - CUDA Mem and copy	
	fprintf(stderr, "ref_source_d->RLP %f MB, ref_target_d->L_tar %f MB, ref_target_d->R_tar %f MB, ctx->out_pair %f MB, Ref SA %f MB, Global %f MB\n", 
			(double)ref_source->toklen*sizeof(int)/(1024*1024), 
			(double)ref_target->toklen*sizeof(uint8_t)/(1024*1024), 
			(double)ref_target->toklen*sizeof(uint8_t)/(1024*1024), 
			(double)ALLOCATIONS*sizeof(res_phrase_t)/(1024*1024),
			(double)ref_source_d ->toklen*sizeof(int)/(1024*1024),
			(double)global*sizeof(saind_t)/(1024*1024));

	timer_start(&__t);
	cudaMalloc((void**)&(ref_source_d->sa), ref_source->toklen*sizeof(int));
	cudaMalloc((void**)&(ref_source_d->RLP), ref_source->toklen*sizeof(int));
	cudaMalloc((void**)&(ref_target_d->L_tar), ref_target->toklen*sizeof(uint8_t));
	cudaMalloc((void**)&(ref_target_d->R_tar), ref_target->toklen*sizeof(uint8_t));
	cudaMalloc((void**)&(ctx->out_pair), ALLOCATIONS*sizeof(res_phrase_t));
	timing->extractin += timer_stop(&__t);

	//fprintf(stderr, "Start Extraction Data Transfer\n");
	timer_start(&__t);
	cudaMemcpy(ref_source_d->sa, (void *)ref_source->sa, ref_source_d->toklen*sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(ref_target_d->R_tar, ref_target->R_tar, ref_target->toklen*sizeof(uint8_t), cudaMemcpyHostToDevice);
	cudaMemcpy(ref_target_d->L_tar, ref_target->L_tar, ref_target->toklen*sizeof(uint8_t), cudaMemcpyHostToDevice);
	cudaMemcpy(ref_source_d->RLP, ref_source->RLP, ref_source->toklen*sizeof(int), cudaMemcpyHostToDevice);
	timing->extractin += timer_stop(&__t);

	dim3  block(THREADS_PER_BLOCK, 1);

	assert(ref_source_d ->sa!=NULL
			&&ref_target_d->L_tar!=NULL
			&&ref_target_d->R_tar!=NULL
			&&ref_source_d->RLP!=NULL
			&&ctx->blocks_d!=NULL
			&&ctx->out_pair!=NULL);

	res_phrase_t* out_res;
	cudaMallocHost((void **)&(out_res), N*ALLOCATIONS*sizeof(res_phrase_t));
	assert(out_res != NULL);

	rule_twogap* twoGapRule;
	cudaMallocHost((void **)&twoGapRule, 2*PREALLOCATION_TWOGAP_RULE*sizeof(rule_twogap));
	assert(twoGapRule != NULL);

	unsigned int count_two_gap = 0;
	unsigned int count_one_gap = 0;

	rule_onegap* oneGapRule;
	cudaMallocHost((void **)&(oneGapRule), 2*PREALLOCATION_TWOGAP_RULE*sizeof(rule_onegap));
	assert(oneGapRule != NULL);

	int prev_cout = 0;
	int prev_block = 0;
	int start_block = 0;
	int end_block = 0;
	bool looseH = false;//false;//Loose Heuristic Switch - large data
	int seperatorOneGap;
	int seperatorTwoGap[2];

	int count = 0, *count_d;
	cudaMalloc((void**)&count_d, sizeof(int));
	cudaMemcpy(count_d, &count, sizeof(int), cudaMemcpyHostToDevice);	

	cudaMemGetInfo(&freeMem, &totalMem);  
	fprintf(stderr, "Memory avaliable: Free: %lu, Total: %lu - Before Allocation for one gap results holder\n",freeMem/(1024*1024), totalMem/(1024*1024));

	rule_onegap* oneGapRule_ab_d;
	cudaMalloc((void**)&(oneGapRule_ab_d), PREALLOCATION_ONEGAP_RULE*sizeof(rule_onegap));

	cudaMemGetInfo(&freeMem, &totalMem);  
	fprintf(stderr, "Memory avaliable: Free: %lu, Total: %lu - After Allocation for one gap holders\n",freeMem/(1024*1024), totalMem/(1024*1024));

	rule_twogap* twoGapRule_ab_d;
	cudaMalloc((void**)&(twoGapRule_ab_d), PREALLOCATION_TWOGAP_RULE*sizeof(rule_twogap));

	cudaMemGetInfo(&freeMem, &totalMem);  
	fprintf(stderr, "Memory avaliable: Free: %lu, Total: %lu - After Allocation for twoGaps holders\n",freeMem/(1024*1024), totalMem/(1024*1024));

	for(int i = 0; i < M; i++){		
		int count_1gap = 0, *count_1gap_d;
		cudaMalloc((void**)&count_1gap_d, sizeof(int));
		cudaMemcpy(count_1gap_d, &count_1gap, sizeof(int), cudaMemcpyHostToDevice);	

		int count_2gap = 0, *count_2gap_d;
		cudaMalloc((void**)&count_2gap_d, sizeof(int));
		cudaMemcpy(count_2gap_d, &count_2gap, sizeof(int), cudaMemcpyHostToDevice);	

		start_block = prev_block;
		if (i == M - 1){
			end_block = global;
		} else {
			end_block = ROUND((float)global*(i+1)/(float)M);
		}	    
		prev_block = end_block;

		dim3  grid(LINEARBLOCK, (end_block - start_block + 1 + LINEARBLOCK - 1)/LINEARBLOCK);
		//fprintf(stderr, "Start Kernel Round %d on Pair Extration - Global %d - Grid X %d | Y %d| Start Block %d | End Block %d\n", i, global, grid.x, grid.y, start_block, end_block);
		timer_start(&__t);

		assert(ref_source_d ->str!=NULL
				&&oneGapRule_ab_d!=NULL
				&&twoGapRule_ab_d!=NULL);

		extractConsistentPairs_Gappy<<<grid, block>>>(
					ref_source_d->sa,
					ref_source->toklen,
					ref_target_d->L_tar,
					ref_target_d->R_tar,
					ref_source_d->RLP,
					ctx->blocks_d,
					end_block,
					start_block,
					ctx->out_pair,
					count_d,
					looseH,
					ref_source_d->str,
					oneGapRule_ab_d,
					twoGapRule_ab_d,
					count_1gap_d,
					count_2gap_d,
					isSample);
			/*extractConsistentPairsSample<<<grid, block>>>(
					ref_source_d->sa,
					ref_source->toklen,
					ref_target_d->L_tar,
					ref_target_d->R_tar,
					ref_source_d->RLP,
					ctx->blocks_d,
					global,
					ctx->out_pair,
					count_d,
					looseH);*/
		cudaThreadSynchronize();
		timing->extractkernel += timer_stop(&__t);
		cudaError_t error = cudaGetLastError();
		fprintf(stderr, "-> Kernel extractConsistentPairs_Gappy ab, abX, Xab time: %f\n", timer_elapsed(&__t)/1000);
		if(error != cudaSuccess)
		{
			// print the CUDA error message and exit if any
			fprintf(stderr, "CUDA error second After kernel extractConsistentPairs_Gappy: %s \n", cudaGetErrorString(error));
			exit(-1);
		}	
		cudaMemcpy(&count, count_d, sizeof(int), cudaMemcpyDeviceToHost);
		fprintf(stderr, "Found %d cotinous pairs!!\n", count);

		cudaMemcpy(&count_1gap, count_1gap_d, sizeof(int), cudaMemcpyDeviceToHost);
		fprintf(stderr, "Found %d 1gap gappy pairs!!\n", count_1gap);

		cudaMemcpy(&count_2gap, count_2gap_d, sizeof(int), cudaMemcpyDeviceToHost);
		fprintf(stderr, "Found %d 2gap gappy pairs!!\n", count_2gap);

		///out_res and fast_speed could be large.
		if (count + prev_cout > N*ALLOCATIONS){
			fprintf(stderr, "Out Of Result Range! %d\n", count+prev_cout);
			return;
		}
		
		///Sort Continouse Phrase results
		thrust::device_ptr<res_phrase_s> out_res_thrust = thrust::device_pointer_cast(ctx->out_pair);
		thrust::sort(out_res_thrust, out_res_thrust+count, continousResCompare());

		cudaMemcpy(out_res+prev_cout, ctx->out_pair, count*sizeof(res_phrase_t), cudaMemcpyDeviceToHost);   
		prev_cout += count;
		cudaFree(ctx->out_pair);
		cudaFree(ref_source_d ->sa);
		cudaFree(ctx->blocks_d);

		cudaMemGetInfo(&freeMem, &totalMem);  
		fprintf(stderr, "Memory avaliable: Free: %lu, Total: %lu - After del continous pairs \n",freeMem/(1024*1024), totalMem/(1024*1024));

		///Sort Two Gap results
		timer_start(&__t);	
		thrust::device_ptr<gappytyp5> twoGapRes_gappy_thrust = thrust::device_pointer_cast(twoGapRule_ab_d);
		thrust::sort(twoGapRes_gappy_thrust, twoGapRes_gappy_thrust+count_2gap, twoGapResCompare());
		cudaMemcpy(twoGapRule, twoGapRule_ab_d, count_2gap*sizeof(rule_twogap), cudaMemcpyDeviceToHost);   
		count_two_gap += count_2gap;
		cudaFree(twoGapRule_ab_d);	

		cudaMemGetInfo(&freeMem, &totalMem);  
		fprintf(stderr, "Memory avaliable: Free: %lu, Total: %lu - After del two gap pairs \n",freeMem/(1024*1024), totalMem/(1024*1024));

		///Sort One Gap results
		thrust::device_ptr<gappytyp4> oneGapRes_gappy_thrust = thrust::device_pointer_cast(oneGapRule_ab_d);
		thrust::sort(oneGapRes_gappy_thrust, oneGapRes_gappy_thrust+count_1gap, oneGapResCompare());

		timer_stop(&__t);
		cudaMemcpy(oneGapRule, oneGapRule_ab_d, count_1gap*sizeof(rule_onegap), cudaMemcpyDeviceToHost);   
		count_one_gap += count_1gap;
		cudaFree(oneGapRule_ab_d);

		//Seperator Book keeping
		seperatorOneGap = count_1gap;
		seperatorTwoGap[0] = count_2gap;
		
		fprintf(stderr, "-> Gappy Results Sorting Thrust time: %f\n", timer_elapsed(&__t)/1000);	  
	}

	/////////////////////////////////////////
	//////////////GAPPY Phrase EXTRACTION////
	/////////////////////////////////////////

	/////////////////////////////////////////
	////BASED ON TWO GAP SEED!! 
	/////////////////////////////////////////

	//Two gap search
	cudaMalloc((void**)&(qryset_d->twoGapSearch), 
			(qryset->distinctTwoGapCount)*sizeof(two_gappy_search));
	cudaMemcpy(qryset_d->twoGapSearch, qryset->twoGapSearch, (qryset->distinctTwoGapCount)*sizeof(two_gappy_search), cudaMemcpyHostToDevice);	

	//two gap SA
	cudaMalloc((void**)&(qryset_d->twoGapSA), 
			(qryset->countTwoGapSA)*sizeof(twoGapOnSA));
	cudaMemcpy(qryset_d->twoGapSA, qryset->twoGapSA, (qryset->countTwoGapSA)*sizeof(twoGapOnSA), cudaMemcpyHostToDevice);	

	//One gap search
	cudaMalloc((void**)&(qryset_d->oneGapSearch), 
			(qryset->distinctOneGapCount)*sizeof(gappy_search));
	cudaMemcpy(qryset_d->oneGapSearch, qryset->oneGapSearch, (qryset->distinctOneGapCount)*sizeof(gappy_search), cudaMemcpyHostToDevice);	

	cudaMalloc((void**)&(twoGapRule_ab_d), (PREALLOCATION_TWOGAP_RULE)*sizeof(rule_twogap));

	cudaMemGetInfo(&freeMem, &totalMem);  
	fprintf(stderr, "Memory avaliable: Free: %lu, Total: %lu - After Allocation for twoGaps \n",freeMem/(1024*1024), totalMem/(1024*1024));

	count = 0;
	cudaMemcpy(count_d, &count, sizeof(int), cudaMemcpyHostToDevice);	

	dim3  gridTwoGap(LINEARBLOCK, (qryset->distinctTwoGapCount + LINEARBLOCK - 1)/LINEARBLOCK);
	dim3  blockTwoGap(THREADS_PER_BLOCK, 1);

	timer_start(&__t);

	fprintf(stderr, "Grid x - %d; y - %d\n", gridTwoGap.x, gridTwoGap.y);
	extractConsistentPairs_TwoGap<<<gridTwoGap, blockTwoGap>>>(
				ref_source->toklen,
				ref_target_d->L_tar,
				ref_target_d->R_tar,
				ref_source_d->RLP,
				qryset->distinctTwoGapCount,
				count_d,
				looseH,
				ref_source_d->str,
				twoGapRule_ab_d,
				qryset_d->twoGapSearch,
				qryset->countTwoGapSA,
				qryset_d->twoGapSA,
				qryset_d->oneGapSearch,
				qryset->distinctOneGapCount,
				isSample);

	cudaThreadSynchronize();
	timer_stop(&__t);
	cudaError_t error = cudaGetLastError();

	if(error != cudaSuccess)
	{
		// print the CUDA error message and exit if any
		fprintf(stderr, "CUDA error second After kernel extractConsistentPairs Two Gap Seeds : %s \n", cudaGetErrorString(error));
		exit(-1);
	}	

	fprintf(stderr, "-> Kernel extractConsistentPairs Two Gap Seeds aXbXc time: %f\n", timer_elapsed(&__t)/1000);
	cudaMemcpy(&count, count_d, sizeof(int), cudaMemcpyDeviceToHost);
	fprintf(stderr, "Found %d two gap phrase extraction pairs!!\n", count);

	///Sort Two Gap results
	timer_start(&__t);	
	thrust::device_ptr<gappytyp5> twoGapRes_1gap_thrust = thrust::device_pointer_cast(twoGapRule_ab_d);
	//sort with thrust based on gappy index attribute - ID
	thrust::sort(twoGapRes_1gap_thrust, twoGapRes_1gap_thrust+count, twoGapResCompare());
	timer_stop(&__t);
	fprintf(stderr, "-> Gappy Results Sorting Thrust time: %f\n", timer_elapsed(&__t)/1000);	  

	cudaMemcpy(twoGapRule+seperatorTwoGap[0], twoGapRule_ab_d, count*sizeof(rule_twogap), cudaMemcpyDeviceToHost);   
	count_two_gap += count;
	seperatorTwoGap[1] = count_two_gap;

	cudaFree(qryset_d->twoGapSA);
	cudaFree(qryset_d->twoGapSearch);
	//cudaFree(twoGapRule_ab_d);

	/////////////////////////////////////////
	////BASED ON One GAP SEED!! 
	/////////////////////////////////////////

	cudaMemGetInfo(&freeMem, &totalMem);  
	fprintf(stderr, "Memory avaliable: Free: %lu, Total: %lu - Start One Gap Seed Phrase Extraction Kernel Continous RULE EXTRACTION \n",freeMem/(1024*1024), totalMem/(1024*1024));

	//One gap SA
	cudaMalloc((void**)&(qryset_d->oneGapSA), 
			(qryset->countOneGapSA)*sizeof(oneGapOnSA));
	cudaMemcpy(qryset_d->oneGapSA, qryset->oneGapSA, qryset->countOneGapSA*sizeof(oneGapOnSA), cudaMemcpyHostToDevice);	

	cudaMalloc((void**)&(oneGapRule_ab_d), 
			(PREALLOCATION_ONEGAP_RULE)*sizeof(rule_onegap));

	cudaMemGetInfo(&freeMem, &totalMem);  
	fprintf(stderr, "Memory avaliable: Free: %lu, Total: %lu - After Allocation for twoGaps \n",freeMem/(1024*1024), totalMem/(1024*1024));

	int count_1 = 0, *count_1_d;
	cudaMalloc((void**)&count_1_d, sizeof(int));
	cudaMemcpy(count_1_d, &count_1, sizeof(int), cudaMemcpyHostToDevice);	

	int count_2 = 0, *count_2_d;
	cudaMalloc((void**)&count_2_d, sizeof(int));
	cudaMemcpy(count_2_d, &count_2, sizeof(int), cudaMemcpyHostToDevice);	

	cudaMalloc((void**)&(ref_source_d->precomp_index), PRECOMPUTECOUNT*PRECOMPUTECOUNT*sizeof(precomp_st_end));
	cudaMemcpy(ref_source_d->precomp_index, ref_source->precomp_index, PRECOMPUTECOUNT*PRECOMPUTECOUNT*sizeof(precomp_st_end), cudaMemcpyHostToDevice);

	cudaMalloc((void**)&(ref_source_d->precomp_onegap), ref_source->precomp_count*sizeof(precompute_enu_3));
	cudaMemcpy(ref_source_d->precomp_onegap, ref_source->precomp_onegap, ref_source->precomp_count*sizeof(precompute_enu_3), cudaMemcpyHostToDevice);

	cudaMemGetInfo(&freeMem, &totalMem);  
	fprintf(stderr, "Memory avaliable: Free: %lu, Total: %lu - Right before One Gap Kernel\n",freeMem/(1024*1024), totalMem/(1024*1024));

	error = cudaGetLastError();

	if(error != cudaSuccess)
	{
		// print the CUDA error message and exit if any
		fprintf(stderr, "CUDA error right before extractConsistentPairs One Gap: %s \n", cudaGetErrorString(error));
		exit(-1);
	}	

	assert(ref_source_d->precomp_onegap
		&& ref_source_d->precomp_index
		&& ref_source_d->RLP
		&& ref_target_d->L_tar
		&& ref_target_d->R_tar
		&& qryset_d->oneGapSA
		&& qryset_d->oneGapSearch);
	
	dim3  gridOneGap(LINEARBLOCK, (qryset->distinctOneGapCount + LINEARBLOCK - 1)/LINEARBLOCK);
	dim3  blockOneGap(THREADS_PER_BLOCK, 1);

	timer_start(&__t);

	fprintf(stderr, "ref_source_d->RLP %f MB, ref_target_d->L_tar %f MB, ref_target_d->R_tar %f MB, PrecompONEGAP %f MB, qryOneGAPonSA results %f MB\n",   (double)ref_source->toklen*sizeof(int)/(1024*1024),
                                (double)ref_target->toklen*sizeof(uint8_t)/(1024*1024),
                                (double)ref_target->toklen*sizeof(uint8_t)/(1024*1024),
                                (double)ref_source->precomp_count*sizeof(precompute_enu_3)/(1024*1024),
                                (double)qryset->countOneGapSA*sizeof(oneGapOnSA)/(1024*1024));

	extractConsistentPairs_OneGap<<<gridOneGap, blockOneGap>>>(
				ref_source->toklen,
				ref_target_d->L_tar,
				ref_target_d->R_tar,
				ref_source_d->RLP,
				qryset->distinctOneGapCount,
				count_1_d,
				count_2_d,
				looseH,
				ref_source_d->str,
				twoGapRule_ab_d,
				oneGapRule_ab_d,
				qryset->countOneGapSA,
				qryset_d->oneGapSA,
				qryset_d->oneGapSearch,
				ref_source_d->precomp_index,
				ref_source_d->precomp_onegap,
				isSample);

	cudaThreadSynchronize();
	timer_stop(&__t);
	error = cudaGetLastError();

	if(error != cudaSuccess)
	{
		// print the CUDA error message and exit if any
		fprintf(stderr, "CUDA error second After kernel extractConsistentPairs One Gap Seeds : %s \n", cudaGetErrorString(error));
		exit(-1);
	}	

	fprintf(stderr, "-> Kernel extractConsistentPairs One Gap Seeds aXb time: %f\n", timer_elapsed(&__t)/1000);
	cudaMemcpy(&count_1, count_1_d, sizeof(int), cudaMemcpyDeviceToHost);
	fprintf(stderr, "Found %d one gap phrase extraction pairs!!\n", count_1);
	cudaMemcpy(&count_2, count_2_d, sizeof(int), cudaMemcpyDeviceToHost);
	fprintf(stderr, "Found %d two gap phrase extraction pairs!!\n", count_2);

	cudaFree(qryset_d->oneGapSA);
	cudaFree(qryset_d->oneGapSearch);
	cudaFree(ref_source_d->precomp_index);
	cudaFree(ref_source_d->precomp_onegap);
	cudaFree(ref_source_d->str);
	cudaFree(ref_source_d->RLP);
	cudaFree(ref_target_d->L_tar);
	cudaFree(ref_target_d->R_tar);

	///Sort Two Gap results
	timer_start(&__t);	
	thrust::device_ptr<gappytyp5> twoGapRes_thrust = thrust::device_pointer_cast(twoGapRule_ab_d);
	//sort with thrust based on gappy index attribute - ID
	thrust::sort(twoGapRes_thrust, twoGapRes_thrust+count_2, twoGapResCompare());
	//Move the resulting two gap phrases out. Update the seperator mark
	cudaMemcpy(twoGapRule+seperatorTwoGap[1], twoGapRule_ab_d, count_2*sizeof(rule_twogap), cudaMemcpyDeviceToHost);   
	count_two_gap += count_2;
	//seperatorTwoGap[1] = count_two_gap;
	cudaFree(twoGapRule_ab_d);
	
	///Sort One Gap results
	thrust::device_ptr<gappytyp4> oneGapRes_thrust = thrust::device_pointer_cast(oneGapRule_ab_d);
	//sort with thrust based on gappyindex attribute - ID
	thrust::sort(oneGapRes_thrust, oneGapRes_thrust+count_1, oneGapResCompare());
	//Move the resulting two gap phrases out. Update the seperator mark
	cudaMemcpy(oneGapRule+count_one_gap, oneGapRule_ab_d, count_1*sizeof(rule_onegap), cudaMemcpyDeviceToHost);   
	seperatorOneGap = count_one_gap;
	count_one_gap += count_1;
	cudaFree(oneGapRule_ab_d);

	timer_stop(&__t);
	fprintf(stderr, "-> Gappy Results Sorting Thrust time: %f\n", timer_elapsed(&__t)/1000);	  

	printf("In total, %d pairs of One Gap phrases and %d pairs of Two Gap Phrases\n",count_one_gap, count_two_gap);	

	cudaMemGetInfo(&freeMem, &totalMem);  
	fprintf(stderr, "Memory avaliable: Free: %lu, Total: %lu - After two gap extraction Done!\n",freeMem/(1024*1024), totalMem/(1024*1024));
	///////////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////

	cudaMemGetInfo(&freeMem, &totalMem);  
	fprintf(stderr, "Memory avaliable: Free: %lu, Total: %lu - After Everything Done!\n",freeMem/(1024*1024), totalMem/(1024*1024));

	cudaFreeHost(ref_target->L_tar);
	cudaFreeHost(ref_target->R_tar);
	cudaFreeHost(ref_source->RLP);
	cudaFreeHost(ref_source->sa);
	cudaFreeHost(ref_source->buf);

	///////////////////////////////////
	///////////////////////////////////
	///////////////////////////////////
	///////////////////////////////////
	///////////////////////////////////
	////Feature Generation Starts!!!///	
	fprintf(stderr, "Start Creating Lexincon\n");

	lexicalTask* maxEF;
	cudaMallocHost((void **)&(maxEF), 
			(PREALLOCATION_LEX_TASK)*sizeof(lexicalTask));//Hard code. Notice!
	unsigned int lexicalTaskCounter = 0;
	////////////////////////////////
	
	int lexicon_count_gappy = 0;
	red_dup_t* fast_speed_one_gap1 = NULL; //(red_dup_t*) malloc(count*sizeof(red_dup_t));
	cudaMallocHost((void **)&(fast_speed_one_gap1), 
			((int)(count_one_gap))*sizeof(red_dup_t));//Hard code. Notice!
	assert(fast_speed_one_gap1!= NULL);

	red_dup_t* fast_speed_one_gap = createLexiconGappyFast(
			ref_source->str, 
			ref_target->str, 
			&lexicon_count_gappy,
			tmp_blocks,
			maxEF,
			sourceNameBlock,
			global,
			oneGapRule,
			twoGapRule,
			count_one_gap,//Total number of one gap results
			count_two_gap,//Total number of two gap results
			seperatorOneGap,//Seperator between Xab, abX; aXb
			seperatorTwoGap[0],//XabX;aXbXc
			seperatorTwoGap[1],//aXbXc;XaXb,aXbX
			qryset->distinctOneGapCount,//Distinct pattern count
			qryset->distinctTwoGapCount,//pattern count - gappy ID max
			qryset->onegapPattern,//source side pattern
			qryset->twogapPattern,
			ref_target->vocabulary,
			ref_source->vocabulary,
			fast_speed_one_gap1,
			qryset->oneGapSearch,
			qryset->twoGapSearch,
			&lexicalTaskCounter,
			qryset->oneGapSA,
			ref_source);//access pattern id, searchGap, position
	assert(lexicon_count_gappy!= 0 
		&& fast_speed_one_gap!= NULL 
		&& lexicalTaskCounter == lexicon_count_gappy
		&& lexicalTaskCounter < PREALLOCATION_LEX_TASK);
	fprintf(stderr, "Lexicon count for aXb, Xab, abX is %d\n", lexicon_count_gappy);

	//No need sort and cannot do sort
	//qsort(fast_speed_one_gap, lexicon_count_gappy, sizeof(red_dup_t), compareUser);
	int totalGappy = 2*global+qryset->distinctOneGapCount;

	qryset->globalOnPairsUpDownGappy= (result_t*)malloc(totalGappy*sizeof(result_t));
	for(int icc = 0; icc< totalGappy; icc++){
		qryset->globalOnPairsUpDownGappy[icc].down = -1;
		qryset->globalOnPairsUpDownGappy[icc].up = -1;					
	}
	for(unsigned int icc = 0; icc < lexicon_count_gappy; icc++){
		if(icc == 0|| fast_speed_one_gap[icc].blocknumber!= 
				fast_speed_one_gap[icc-1].blocknumber){
			qryset->globalOnPairsUpDownGappy[fast_speed_one_gap[icc].blocknumber].down = icc;
		}
		qryset->globalOnPairsUpDownGappy[fast_speed_one_gap[icc].blocknumber].up = icc;
	}

	////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////
	int lexicon_count_two_gap = 0;
	red_dup_t* fast_speed_two_gap1 = NULL; //(red_dup_t*) malloc(count*sizeof(red_dup_t));
	cudaMallocHost((void **)&(fast_speed_two_gap1), 
			((int)(count_two_gap))*sizeof(red_dup_t));//Hard code. Notice!
	assert(fast_speed_two_gap1!= NULL);

	red_dup_t* fast_speed_two_gap = createLexiconTwoGapFast(
			ref_source->str, 
			ref_target->str, 
			&lexicon_count_two_gap,
			tmp_blocks,
			maxEF,
			sourceNameBlock,
			global,
			oneGapRule,
			twoGapRule,
			count_one_gap,//Total number of one gap results
			count_two_gap,//Total number of two gap results
			seperatorOneGap,//Seperator between Xab, abX; aXb
			seperatorTwoGap[0],//XabX;aXbXc
			seperatorTwoGap[1],//aXbXc;XaXb,aXbX
			qryset->distinctOneGapCount,//Distinct pattern count
			qryset->distinctTwoGapCount,//pattern count - gappy ID max
			qryset->onegapPattern,//source side pattern
			qryset->twogapPattern,
			ref_target->vocabulary,
			ref_source->vocabulary,
			fast_speed_two_gap1,
			qryset->oneGapSearch,
			qryset->twoGapSearch,
			&lexicalTaskCounter,
			qryset->oneGapSA,
			ref_source);//access pattern id, searchGap, position

	assert(lexicon_count_two_gap!= 0 
		&& fast_speed_two_gap!= NULL
		&& lexicalTaskCounter == (lexicon_count_gappy+lexicon_count_two_gap)
		&& lexicalTaskCounter < PREALLOCATION_LEX_TASK);
	fprintf(stderr, "Lexicon count for aXbXc, XabX, XaXb, aXbX is %d\n", lexicon_count_two_gap);

	//No sort please
	//qsort(fast_speed_two_gap, lexicon_count_two_gap, sizeof(red_dup_t), compareUser);
	int totalTwoGap = global+
		2*qryset->distinctOneGapCount+
		qryset->distinctTwoGapCount;
	qryset->globalOnPairsUpDownTwoGap= (result_t*)malloc(totalTwoGap*sizeof(result_t));
	for(int icc = 0; icc< totalTwoGap; icc++){
		qryset->globalOnPairsUpDownTwoGap[icc].down = -1;
		qryset->globalOnPairsUpDownTwoGap[icc].up = -1;					
	}
	for(unsigned int icc = 0; icc < lexicon_count_two_gap; icc++){
		if(icc == 0|| fast_speed_two_gap[icc].blocknumber!= 
				fast_speed_two_gap[icc-1].blocknumber){
			qryset->globalOnPairsUpDownTwoGap[fast_speed_two_gap[icc].blocknumber].down = icc;
		}
		qryset->globalOnPairsUpDownTwoGap[fast_speed_two_gap[icc].blocknumber].up = icc;
	}
	
	//Debug
	cout<<"Total Two Gap ID possible range:"<<totalTwoGap<<"  real two gap lexicon range:"<<lexicon_count_two_gap<<endl;
	/*
	for(unsigned int icc =0; icc < totalTwoGap; icc++){
		if(qryset->globalOnPairsUpDownTwoGap[icc].down !=-1 &&
				qryset->globalOnPairsUpDownTwoGap[icc].up !=-1 ){
			if(qryset->globalOnPairsUpDownTwoGap[icc].down >= lexicon_count_two_gap || qryset->globalOnPairsUpDownTwoGap[icc].up >= lexicon_count_two_gap ){
				printf("Exceed!!!%d -> Down %d | Up %d\n", icc, qryset->globalOnPairsUpDownTwoGap[icc].down, qryset->globalOnPairsUpDownTwoGap[icc].up);		
			} else {
				printf("%d -> Down %d | Up %d\n", icc, qryset->globalOnPairsUpDownTwoGap[icc].down, qryset->globalOnPairsUpDownTwoGap[icc].up);		
			}
		} else if (qryset->globalOnPairsUpDownTwoGap[icc].down ==-1 &&
				qryset->globalOnPairsUpDownTwoGap[icc].up ==-1 ){
			continue;		
		} else {
			printf("WRONG!!! %d -> Down %d | Up %d\n", icc, qryset->globalOnPairsUpDownTwoGap[icc].down, qryset->globalOnPairsUpDownTwoGap[icc].up);		
		}
	}*/
	
	//////////////////////////////////////////////////
	//////////////////////////////////////////////////
	/////////////Start processing nongappy features///
	//////////////////////////////////////////////////
	//////////////////////////////////////////////////

	hashtbl_aux* lexic = NULL;
	hash_lexicon* target_count = NULL;
	hash_lexicon* foreig_count = NULL;

	int lexicon_count = 0;
	red_dup_t* fast_speed1 = NULL; //(red_dup_t*) malloc(count*sizeof(red_dup_t));
	cudaMallocHost((void **)&(fast_speed1), 
			((int)(0.55*prev_cout))*sizeof(red_dup_t));//Hard code. Notice!
	assert(fast_speed1!= NULL);

	red_dup_t* fast_speed = createLexiconFast(
			prev_cout, 
			out_res, 
			ref_source->str, 
			ref_target->str, 
			&lexic, 
			&target_count,
			&foreig_count, 
			&lexicon_count,
			tmp_blocks,
			maxEF,
			sourceNameBlock,
			global,
			fast_speed1,
			ref_target->vocabulary,
			&lexicalTaskCounter);

	assert(lexicon_count != 0 
		&& fast_speed != NULL
		&& lexicalTaskCounter == (lexicon_count_gappy+lexicon_count_two_gap+lexicon_count)
		&& lexicalTaskCounter < PREALLOCATION_LEX_TASK);
	fprintf(stderr, "Lexicon count for continous ab is %d\n", lexicon_count);

	//To be changed to thrust sorting! here
	//qsort(fast_speed, lexicon_count, sizeof(red_dup_t), compareUser);
	free(tmp_blocks);

	//fprintf(stderr, "Start Extract Global Pairs - GPU Mode\n");
	qryset->globalOnPairsUpDownContinous = (result_t*)malloc(global*sizeof(result_t));
	red_dup_t* fast_speed_d = NULL;

	cudaFreeHost(out_res);
	//cudaFree(ctx->blocks_d);

	cudaMalloc((void**)&(ctx->globalOnPairsUpDown_d), global*sizeof(result_t));		
	cudaMalloc((void**)&(fast_speed_d), lexicon_count*sizeof(red_dup_t));
	cudaMemcpy(fast_speed_d, fast_speed, lexicon_count*sizeof(red_dup_t), cudaMemcpyHostToDevice);  
	assert(ctx->globalOnPairsUpDown_d!=NULL && fast_speed_d!=NULL);

	///For each global, get its up and lower boundary
	timer_start(&__t);
	dim3  block1(THREADS_PER_BLOCK_GLOBALPAIRS, 1);
	dim3  grid1((global + THREADS_PER_BLOCK_GLOBALPAIRS - 1)/THREADS_PER_BLOCK_GLOBALPAIRS, 2);

	extractGlobalPairsUpDown<<<grid1, block1>>>(
			global,
			fast_speed_d,
			lexicon_count,
			ctx->globalOnPairsUpDown_d);

	cudaThreadSynchronize();
	cudaMemcpy(qryset->globalOnPairsUpDownContinous, ctx->globalOnPairsUpDown_d, global*sizeof(result_t), cudaMemcpyDeviceToHost);
	timing->extractkernel += timer_stop(&__t);

	cudaError_t error1 = cudaGetLastError();
	if(error1 != cudaSuccess)
	{
		// print the CUDA error message and exit if any 
		fprintf(stderr, "CUDA error second After kernel galobal pairs up down: %s \n", cudaGetErrorString(error1));
		exit(-1);
	}
	
	cudaFree(fast_speed_d);
	cudaFree(ctx->globalOnPairsUpDown_d);

	///////////////////////////////////////////////////
	///LexicalTask: MaxEgivenF and MaxFGivenE on GPU///
	///////////////////////////////////////////////////	
	timer_start(&__t);	

	lexicalTask* maxEF_d = NULL;
	cudaMalloc((void**)&(maxEF_d), lexicalTaskCounter*sizeof(lexicalTask));
	cudaMemcpy(maxEF_d, maxEF, lexicalTaskCounter*sizeof(lexicalTask), cudaMemcpyHostToDevice);  

	cudaMalloc((void**)&(ref_target_d->str), ref_target->toklen*sizeof(int));
	cudaMemcpy(ref_target_d->str, ref_target->str, ref_target->toklen*sizeof(int), cudaMemcpyHostToDevice);

	lexTaskResults* resultLex_d = NULL; 
	lexTaskResults* resultLex = NULL; 
	
	cudaMallocHost((void **)&(resultLex), lexicalTaskCounter*sizeof(lexTaskResults));
	cudaMalloc((void**)&(resultLex_d), lexicalTaskCounter*sizeof(lexTaskResults));
	
	cudaMemGetInfo(&freeMem, &totalMem);  
	fprintf(stderr, "Memory avaliable: Free: %lu, Total: %lu - Before Lexical task kernel!\n",freeMem/(1024*1024), totalMem/(1024*1024));

	assert(maxEF && maxEF_d && resultLex && ref_target_d->str && resultLex_d);
	dim3  blockLex(THREADS_PER_BLOCK_GLOBALPAIRS, 1);
	dim3  gridLex(100, (lexicalTaskCounter+blockLex.x*100 - 1)/(blockLex.x*100));

	lexicalTaskMaxEF<<<gridLex, blockLex>>>(
		lexicon_count_gappy, //one gap lexicon boundary
		lexicon_count_gappy+lexicon_count_two_gap, //two gap lexicon boundary
		lexicalTaskCounter, //Total number of lexicon
		maxEF_d, //all tasks are here: MaxEF
		ref_target_d->str,
		ref_target->toklen,//target string boundary checking
		cudaLexFile->keys,//search among those keys
		cudaLexFile->values,//retrieve the values and get max
		cudaLexFile->count,//search range
		resultLex_d);
	cudaThreadSynchronize();
	
	error1 = cudaGetLastError();
	if(error1 != cudaSuccess)
	{
		// print the CUDA error message and exit if any 
		fprintf(stderr, "CUDA error second After kernel lexical task: %s \n", cudaGetErrorString(error1));
		exit(-1);
	}	
	cudaMemcpy(resultLex, resultLex_d, lexicalTaskCounter*sizeof(lexTaskResults), cudaMemcpyDeviceToHost);

	for(int iccc = 0; iccc < lexicalTaskCounter; iccc++){
		if(iccc < lexicon_count_gappy){ //one gap
			fast_speed_one_gap[maxEF[iccc].fastSpeedId].MaxLexEgivenF = 
						resultLex[iccc].MaxLexEgivenF;
			fast_speed_one_gap[maxEF[iccc].fastSpeedId].MaxLexFgivenE = 
						resultLex[iccc].MaxLexFgivenE;
		} else if (iccc < lexicon_count_gappy+lexicon_count_two_gap) { // two gap
			fast_speed_two_gap[maxEF[iccc].fastSpeedId].MaxLexEgivenF = 
						resultLex[iccc].MaxLexEgivenF;
			fast_speed_two_gap[maxEF[iccc].fastSpeedId].MaxLexFgivenE = 
						resultLex[iccc].MaxLexFgivenE;
		} else { //continous
			fast_speed[maxEF[iccc].fastSpeedId].MaxLexEgivenF = 
						resultLex[iccc].MaxLexEgivenF;
			fast_speed[maxEF[iccc].fastSpeedId].MaxLexFgivenE = 
						resultLex[iccc].MaxLexFgivenE;
		}
	}
	
	cudaFree(maxEF_d);
	cudaFree(ref_target_d->str);
	cudaFree(resultLex_d);
	cudaFreeHost(maxEF);
	cudaFreeHost(resultLex);
	
	timer_stop(&__t);	
	fprintf(stderr, "-> Lexical Task CPU/GPU time: %f\n", timer_elapsed(&__t)/1000);	
	///////////////////////////////////////////////////
	//////////////////////////PRINT////////////////////
	///////////////////////////////////////////////////
	qryset->continousQueryWithID = qryglobal;
	qryset->distinctContinousCount  = global;
	qryset->fast_speed = fast_speed;
	qryset->fast_speed_one_gap = fast_speed_one_gap;
	qryset->fast_speed_two_gap = fast_speed_two_gap;

}

void ExtractPairs_Large_Data(void *handler,
		ref_set_t *refset, 
		qry_set_t *qryset, 
		timing_t *timing, 
		map<string, categ> word_score,
		disk_handler_t *diskInfo) {
	mytimer_t __t;	
	char* destDir = diskInfo->destinationDirectory;	
	context_t * ctx = (context_t *) handler;
	ref_t_target * ref_target_d = &ctx->ref_target_d;			 
	ref_t * ref_source_d = &ctx->ref_d; 		 
	qry_set_t *qryset_d = &ctx->qryset_d;
	ref_t_target * ref_target = &refset->refs_target[0];			 
	ref_t * ref_source = &refset->refs[0];		

	saind_t* tmp_blocks = (saind_t*)malloc(sizeof(saind_t)*TEMPSET*(qryset->totaltokens));
	unsigned int global = 0;    
	vector<char*> sourceNameBlock(TEMPSET*(qryset->totaltokens));// = (char**)malloc(2*(qryset->totaltokens)*sizeof(char*));
	
	//qryglobal = GenerateBlocks(handler, refset, qryset, tmp_blocks, &global, &sourceNameBlock, intchar, qryglobal);	
	vector<vector<unsigned int> >  qryglobal(qryset->qryscount);
	GenerateBlocks(handler, 
					refset, 
					qryset, 
					tmp_blocks, 
					&global, 
					&sourceNameBlock, 
					qryglobal,
					ref_source->vocabulary);
	if(qryglobal.size()==0 || global > TEMPSET*qryset->totaltokens){
		fprintf(stderr, "Function generateblocks() Failed Due to Physical Memory Limitations\n");
		return;
	}	
	
	cudaMalloc((void**)&(ctx->blocks_d), (global)*sizeof(saind_t));
	cudaMemcpy(ctx->blocks_d, tmp_blocks, (global)*sizeof(saind_t), cudaMemcpyHostToDevice);

	///L and R in target side - CUDA mem location and copy
	///RLP in source side - CUDA Mem and copy	
	fprintf(stderr, "ref_source_d->RLP %f MB, ref_target_d->L_tar %f MB, ref_target_d->R_tar %f MB, ctx->out_pair %f MB, Ref SA %f MB, Global %f MB\n", 
			(double)ref_source->toklen*sizeof(int)/(1024*1024), 
			(double)ref_target->toklen*sizeof(uint8_t)/(1024*1024), 
			(double)ref_target->toklen*sizeof(uint8_t)/(1024*1024), 
			(double)ALLOCATIONS*sizeof(res_phrase_t)/(1024*1024),
			(double)ref_source_d ->toklen*sizeof(int)/(1024*1024),
			(double)global*sizeof(saind_t)/(1024*1024));

	timer_start(&__t);
	cudaMalloc((void**)&(ref_source_d->RLP), ref_source->toklen*sizeof(int));
	cudaMalloc((void**)&(ref_target_d->L_tar), ref_target->toklen*sizeof(uint8_t));
	cudaMalloc((void**)&(ref_target_d->R_tar), ref_target->toklen*sizeof(uint8_t));
	cudaMalloc((void**)&(ctx->out_pair), ALLOCATIONS*sizeof(res_phrase_t));
	timing->extractin += timer_stop(&__t);

	//fprintf(stderr, "Start Extraction Data Transfer\n");
	timer_start(&__t);
	cudaMemcpy(ref_target_d->R_tar, ref_target->R_tar, ref_target->toklen*sizeof(uint8_t), cudaMemcpyHostToDevice);
	cudaMemcpy(ref_target_d->L_tar, ref_target->L_tar, ref_target->toklen*sizeof(uint8_t), cudaMemcpyHostToDevice);
	cudaMemcpy(ref_source_d->RLP, ref_source->RLP, ref_source->toklen*sizeof(int), cudaMemcpyHostToDevice);
	timing->extractin += timer_stop(&__t);

	dim3  block(THREADS_PER_BLOCK, 1);

	assert(ref_source_d->sa!=NULL
			&&ref_target_d->L_tar!=NULL
			&&ref_target_d->R_tar!=NULL
			&&ref_source_d->RLP!=NULL
			&&ctx->blocks_d!=NULL
			&&ctx->out_pair!=NULL);

	res_phrase_t* out_res;
	cudaMallocHost((void **)&(out_res), N*ALLOCATIONS*sizeof(res_phrase_t));
	assert(out_res != NULL);

	int prev_cout = 0;
	int prev_block = 0;
	int start_block = 0;
	int end_block = 0;
	//bool isSample = false;//true;//rue;    
	bool looseH = false;//false;//Loose Heuristic Switch - large data
	for(int i = 0; i < M; i++){		
		int count = 0, *count_d;
		cudaMalloc((void**)&count_d, sizeof(int));
		cudaMemcpy(count_d, &count, sizeof(int), cudaMemcpyHostToDevice);	

		start_block = prev_block;
		if (i == M - 1){
			end_block = global;
		} else {
			end_block = ROUND((float)global*(i+1)/(float)M);
		}	    
		prev_block = end_block;

		dim3  grid(LINEARBLOCK, (end_block - start_block + 1 + LINEARBLOCK - 1)/LINEARBLOCK);
		//fprintf(stderr, "Start Kernel Round %d on Pair Extration - Global %d - Grid X %d | Y %d| Start Block %d | End Block %d\n", i, global, grid.x, grid.y, start_block, end_block);
		timer_start(&__t);

		if(!isSample){
			extractConsistentPairs<<<grid, block>>>(
					ref_source_d->sa,
					ref_source->toklen,
					ref_target_d->L_tar,
					ref_target_d->R_tar,
					ref_source_d->RLP,
					ctx->blocks_d,
					end_block,
					start_block,
					ctx->out_pair,
					count_d,
					looseH);
		} else {
			extractConsistentPairsSample<<<grid, block>>>(
					ref_source_d->sa,
					ref_source->toklen,
					ref_target_d->L_tar,
					ref_target_d->R_tar,
					ref_source_d->RLP,
					ctx->blocks_d,
					global,
					ctx->out_pair,
					count_d,
					looseH);
		}
		cudaThreadSynchronize();
		timing->extractkernel += timer_stop(&__t);
		cudaError_t error = cudaGetLastError();
		fprintf(stderr, "Kernel Extraction time: %f\n", timer_elapsed(&__t)/1000);
		if(error != cudaSuccess)
		{
			// print the CUDA error message and exit if any
			fprintf(stderr, "CUDA error second After kernel extractConsistentPairs: %s \n", cudaGetErrorString(error));
			exit(-1);
		}	
		cudaMemcpy(&count, count_d, sizeof(int), cudaMemcpyDeviceToHost);
		fprintf(stderr, "Found %d pairs!!\n", count);

		///out_res and fast_speed could be large.
		if (count + prev_cout > N*ALLOCATIONS){
			fprintf(stderr, "Out Of Result Range! %d\n", count+prev_cout);
			return;
		}
		cudaMemcpy(out_res+prev_cout, ctx->out_pair, count*sizeof(res_phrase_t), cudaMemcpyDeviceToHost);   
		prev_cout += count;
	}
	//cudaFree(ctx->out_pair);
	//cudaFree(ref_source_d ->sa);
	//cudaFree(ctx->blocks_d);

	cudaFreeHost(ref_target->L_tar);
	cudaFreeHost(ref_target->R_tar);
	cudaFreeHost(ref_source->RLP);
	cudaFreeHost(ref_source->sa);
	cudaFreeHost(ref_source->buf);

	hashtbl_aux* lexic = NULL;
	hash_lexicon* target_count = NULL;
	hash_lexicon* foreig_count = NULL;
	fprintf(stderr, "Total pairs %d\nStart Creating Lexincon\n", prev_cout);

	int lexicon_count = 0;
	red_dup_t* fast_speed1 = NULL; //(red_dup_t*) malloc(count*sizeof(red_dup_t));
	cudaMallocHost((void **)&(fast_speed1), ((int)(0.53*prev_cout))*sizeof(red_dup_t));
	assert(fast_speed1!= NULL);

	red_dup_t* fast_speed = createLexicon(
			prev_cout, 
			out_res, 
			ref_source ->str, 
			ref_target->str, 
			&lexic, 
			&target_count,
			&foreig_count, 
			&lexicon_count,
			tmp_blocks,
			word_score,
			sourceNameBlock,
			global,
			fast_speed1,
			ref_target->vocabulary);

	assert(lexicon_count != 0 && fast_speed != NULL);
	fprintf(stderr, "Lexicon count is %d\n", lexicon_count);
	qsort(fast_speed, lexicon_count, sizeof(red_dup_t), compareUser);

	free(tmp_blocks);

	bool cpu = 0;
	if (cpu){
		fprintf(stderr, "Start Print Grammar\n");
	} else {
		//fprintf(stderr, "Start Extract Global Pairs - GPU Mode\n");
		result_t* globalOnPairsUpDown = (result_t*)malloc(global*sizeof(result_t));
		red_dup_t* fast_speed_d = NULL;

		cudaFree(ref_source_d ->RLP);
		cudaFree(ref_target_d->L_tar);
		cudaFree(ref_target_d->R_tar);
		cudaFreeHost(out_res);
		//cudaFree(ctx->blocks_d);

		cudaMalloc((void**)&(ctx->globalOnPairsUpDown_d), global*sizeof(result_t));		
		cudaMalloc((void**)&(fast_speed_d), lexicon_count*sizeof(red_dup_t));
		cudaMemcpy(fast_speed_d, fast_speed, lexicon_count*sizeof(red_dup_t), cudaMemcpyHostToDevice);  
		assert(ctx->globalOnPairsUpDown_d!=NULL && fast_speed_d!=NULL);

		///For each global, get its up and lower boundary
		timer_start(&__t);
		dim3  block1(THREADS_PER_BLOCK_GLOBALPAIRS, 1);
		dim3  grid1((global + THREADS_PER_BLOCK_GLOBALPAIRS - 1)/THREADS_PER_BLOCK_GLOBALPAIRS, 2);

		extractGlobalPairsUpDown<<<grid1, block1>>>(
				global,
				fast_speed_d,
				lexicon_count,
				ctx->globalOnPairsUpDown_d);

		cudaThreadSynchronize();
		cudaMemcpy(globalOnPairsUpDown, ctx->globalOnPairsUpDown_d, global*sizeof(result_t), cudaMemcpyDeviceToHost);
		timing->extractkernel += timer_stop(&__t);

		cudaError_t error1 = cudaGetLastError();
		if(error1 != cudaSuccess)
		{
			// print the CUDA error message and exit if any 
			fprintf(stderr, "CUDA error second After kernel extractGlobalPairsUpDown: %s -old \n", cudaGetErrorString(error1));
			exit(-1);
		}
	}	
}
