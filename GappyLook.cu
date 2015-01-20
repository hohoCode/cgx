#include "ComTypes.h"
#include "GappyLook.h"


__device__ int existPrecomputation(int* topPrecomputationList, int tokenA, int tokenB){
	int start = 0;
	int end = PRECOMPUTECOUNT-1;
	int middle;
	bool flag_a = true;
	while (end - start >= 0 && flag_a) {
		middle = (start + end) >> 1;
		if (topPrecomputationList[middle] > tokenA){
			end = middle - 1;						
		} else if (topPrecomputationList[middle] < tokenA){
			start = middle + 1;
		} else {
			flag_a = false;			
		}
	}

	start = 0;
	end = PRECOMPUTECOUNT-1;
	bool flag_b = true;
	int middle2;
	while (end - start >= 0 && flag_b) {
		middle2 = (start + end) >> 1;
		if (topPrecomputationList[middle2] > tokenB){
			end = middle2 - 1;						
		} else if (topPrecomputationList[middle2] < tokenB){
			start = middle2 + 1;
		} else {
			flag_b = false;			
		}
	}

	if(!flag_a && !flag_b){
		return middle*PRECOMPUTECOUNT+middle2;
	}
	return -1;
}


__device__ bool checkBoundaryGap(//Target checking required, not no interface changes.
		unsigned int start, 
		unsigned int ender,
		uint8_t* L_tar,
		uint8_t* R_tar,
		unsigned int* RLP){
	unsigned char L = 0;
	unsigned char R = 0;
	int sen_target_begin = -1;
	unsigned char min_L = 255;
	unsigned char max_R = 0;
	int tempind = 0;
	unsigned int temp;
	
	int target_start;
	int target_end;
	int k;
	bool returnVal;
	for(k = start; k <= ender; k++){
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

	if (min_L <= max_R && max_R - min_L < MAX_rule_span){
		tempind++;
		target_start = min_L + sen_target_begin;
		target_end = max_R + sen_target_begin;
		//printf("target start %d - tar end %d - sour start %d - sour end %d\n", ss, tt, current_str, t);
		//return consistent( (*target_start),  (*target_end), L_tar, R_tar, start, ender, tempind);
		//bool consistent(int start, int end, uint8_t* L_target, uint8_t* R_target, int start_chk, int end_chk, int startpos_source)
		returnVal = true;
		min_L = 255;
		max_R = 0;
		for(k = target_start; k <= target_end; k++){
			L = L_tar[k];
			R = R_tar[k];
			
			if (L==255 || R == 255){
				returnVal = true; //Change fixes original false to remove non-recog symbols.
			} else if (k == target_start){
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
			
		if(tempind+min_L != start || tempind+max_R != ender){
			returnVal = false;
		}
		return returnVal;
	}else {	
		return false;
	}
}

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
		uint8_t* R_tar) {

	int blockId = blockIdx.y*gridDim.x+blockIdx.x;
	
	/*if(blockId == 0){
		for(int i = 0; i<distinctOneGapCount;i++){
			printf("BlockId %d - qrystart %d - searchInd %d - gap %d - startlen %d\n", 
				i, 
				oneGapSearch[i].qrystart, 
				oneGapSearch[i].qrystart+oneGapSearch[i].gap + oneGapSearch[i].qrystart_len,
				oneGapSearch[i].gap,
				oneGapSearch[i].qrystart_len);
		}
	}*/
	if (blockId >= distinctOneGapCount){
		return;
	}
	
	int longest_len_start = oneGapSearch[blockId].qrystart_len;
	int longest_len_end = oneGapSearch[blockId].qryend_len;
	int range1_up = -1;
	int range1_down = -1;
	int range2_up = -1;
	int range2_down = -1;

	int precomputationIndex = -1;
	int tokindex = oneGapSearch[blockId].qrystart;//aXb's a, first letter of a
	int search_tokindex = tokindex+oneGapSearch[blockId].gap + 
		oneGapSearch[blockId].qrystart_len;//aXb's b, first letter of b.
	int nextpos = -1;
	
	//if(threadIdx.x==0){ //threadIdx.x
	/*printf("BlockId %d - qrystart %d - searchInd %d - gap %d - startlen %d\n", 
			blockId, 
			tokindex, 
			search_tokindex,
			oneGapSearch[blockId].gap,
			oneGapSearch[blockId].qrystart_len);*/
	//}
	if(oneGapSearch[blockId].gap==0||tokindex <0/*||tokindex > 11||search_tokindex>11*/){
		printf("Wrong!!\n");
		/*printf("Wrong!!! BlockId %d - qrystart %d - searchInd %d - gap %d - startlen %d\n", 
			blockId, 
			tokindex, 
			search_tokindex,
			oneGapSearch[blockId].gap,
			oneGapSearch[blockId].qrystart_len);*/
		return;
	}
	int tobesearch_start;
	int tobesearch_end;
	int dis, dis2; 
	int threadx = threadIdx.x;
	int cc;
	bool forwardOrBack = true;

	//Sanity check
	if (qryresult[search_tokindex].longestmatch < longest_len_end || 
			qryresult[tokindex].longestmatch < longest_len_start){
		printf("Not possible in side one gap SA]\n");
		return;
	}

	//precomputation possible
	precomputationIndex = existPrecomputation(frequentList, 
			qrysoffsettok[qryscount+tokindex+longest_len_start-1], 
			qrysoffsettok[qryscount+search_tokindex]);		
	
	if (precomputationIndex == -1){
		if(longest_len_start == 1){
			range1_up = qryresult[tokindex].up;
			range1_down = qryresult[tokindex].down;
			dis = range1_down - range1_up;
		} else {
			cc=connectoffset[tokindex] + longest_len_start - 2;
			range1_up = result_connect[cc].up;
			range1_down = result_connect[cc].down;
			dis = range1_down - range1_up;
		}

		if(longest_len_end == 1){
			range2_up = qryresult[search_tokindex].up;
			range2_down = qryresult[search_tokindex].down;
			dis2 = range2_down - range2_up;
		} else {
			cc=connectoffset[search_tokindex] + longest_len_end - 2;
			range2_up = result_connect[cc].up;
			range2_down = result_connect[cc].down;
			dis2 = range2_down - range2_up;
		}

		if(dis <= dis2){
			tobesearch_start = range1_up;
			tobesearch_end = range1_down;
			forwardOrBack = true;
		} else {
			dis = dis2;
			tobesearch_start = range2_up;
			tobesearch_end = range2_down;
			forwardOrBack = false;
		}

	}else {
		tobesearch_start = precomp_index[precomputationIndex].start;
		tobesearch_end = precomp_index[precomputationIndex].end;
		dis = tobesearch_end - tobesearch_start;		
	}

	if(precomputationIndex != -1 && longest_len_start == 1 && longest_len_end == 1
		 && dis >= 0){
		//Sanity check
		if (precomputationIndex < 0){
			printf("Not possible! Precomp less than zero\n");
			return;
		}
		if (threadx == 0){
			nextpos = atomicAdd(counter, 1);
			oneGapSA[nextpos].position  = blockId;
			oneGapSA[nextpos].str_position = precomputationIndex;
			oneGapSA[nextpos].length = 0;
		}
		return;
	}	
	
	int move = 0;
	bool flager = true;
	int precompstart;
	int precomplen;

	int backoff = 0;
	bool stop = false;
	int forward = -1;
	int gostart = 0;
	int temp =-1;
	int matchcount;

	while(threadx <= dis){
		move = 0;
		flager = true;
		if (precomputationIndex != -1){
			precompstart = precomp_onegap[tobesearch_start+threadx].start;
			precomplen = precomp_onegap[tobesearch_start+threadx].length;

			if (precomplen + 1 + longest_len_start - 1 + 
					longest_len_end - 1 > MAX_rule_span ){
				flager = false;
			}

			//Check the previous tokens are the same otherwise stop
			if (flager && longest_len_start > 1){
				backoff = 0;
				stop = false;
				while(flager && !stop){
					backoff++;
					if(precompstart-backoff<0||refstr[precompstart-backoff]!=
							qrysoffsettok[qryscount+tokindex+longest_len_start-1-backoff]){
						flager= false;
					}
					if(longest_len_start - backoff <= 1){
						stop = true;
					}
				}
			}

			//Check the B and afterwards
			if(flager&& longest_len_end > 1){
				forward = 1;
				while(forward < longest_len_end && flager){
					forward++;
					if (refstr[precompstart + precomplen + forward-1] 
							!= qrysoffsettok[qryscount+search_tokindex+forward-1]){
						flager = false;
					}
				}
			}

			//Everything is OK, record results
			if(flager){
				nextpos = atomicAdd(counter, 1);
				oneGapSA[nextpos].position	= blockId;
				oneGapSA[nextpos].str_position = precompstart - 
						longest_len_start + 1;
				oneGapSA[nextpos].length = precomplen + longest_len_start - 1 + 
						longest_len_end - 1;  // to the end aXb's b					
			}
		} else if (forwardOrBack){
			///Forward search, start from aXb's a, because a's range is smaller.
			gostart = refsa[threadx+tobesearch_start];
			move = 0;
			temp = -1;
			while(flager){
				if(move == 0){
					//Check gap, if the gap is a <2 between two sentences, 
					//then skip
					temp = refstr[gostart+longest_len_start];
					if(temp < 2){
						flager = false;
					}
				}
				temp = refstr[gostart+longest_len_start+MIN_gap_size+move];
				if (temp < 2 ){
					flager = false;
				} else if (flager && temp == qrysoffsettok[qryscount+search_tokindex]){
					//onegap[nextpos].gap = longest_start_iter+MIN_gap_size+move;
					matchcount = 1;
					stop = false;
					while(!stop&&matchcount < longest_len_end
							&& longest_len_start+MIN_gap_size+move + 1 + matchcount <= MAX_rule_span){ 
						////debugging
						if (qrysoffsettok[qryscount+search_tokindex+matchcount]< 2){
							printf("This is not possible - onegap kernel");
							return;
						}
						//debugging
						
						backoff = refstr[gostart+longest_len_start+MIN_gap_size+move+matchcount];
						if (backoff < 2){
							stop = true;
							flager = false;
						} else if(backoff == 
								qrysoffsettok[qryscount+search_tokindex+matchcount]){
							matchcount++;
						}	else {
							stop = true;
						}				
					}

					if(matchcount == longest_len_end && 
						 checkBoundaryGap(gostart+longest_len_start, 
							gostart+longest_len_start+MIN_gap_size+move+longest_len_end-1-longest_len_end, 
							L_tar, 
							R_tar, 
							RLP)){
						nextpos = atomicAdd(counter, 1);
						oneGapSA[nextpos].position	= blockId;
						oneGapSA[nextpos].str_position = gostart;
						oneGapSA[nextpos].length = 
							longest_len_start+MIN_gap_size+move+longest_len_end-1; 
						// to the end aXb's b	
					}
				} 

				move++;//+move is the correct pos of b's starting pos
				if(longest_len_start+MIN_gap_size+move + longest_len_end > MAX_rule_span){
					flager = false;
				}
			}						
		} else {
			///Backward search, start from aXb's b's first character, because b's range is smaller.
			gostart = refsa[threadx+tobesearch_start];
			move = 0;
			temp = -1;
			while(flager){
				if(move == 0){
					//Check gap, if the gap is a <2 between two sentences, 
					//then skip
					temp = refstr[gostart-1];
					if(temp < 2){
						flager = false;
					}
				}
				if(gostart-1-MIN_gap_size-move<0){
					temp = -1;
				} else {
					temp = refstr[gostart-1-MIN_gap_size-move];
				}
				
				if (temp < 2 ){
					flager = false;
				} else if (flager && temp == qrysoffsettok[qryscount+tokindex+longest_len_start-1]){
					//onegap[nextpos].gap = longest_start_iter+MIN_gap_size+move;
					matchcount = 1;
					stop = false;
					while(!stop&&matchcount < longest_len_start
							&& longest_len_end+MIN_gap_size+move + 1 + matchcount <= MAX_rule_span){
						///debuginig
						if (qrysoffsettok[qryscount+tokindex+longest_len_start-1-matchcount]
								< 2){
							printf("One gap lookup kernel- This is not possible\n");
							return;
						} 
						///Debugging
						if (gostart-1-MIN_gap_size-move-matchcount < 0){
							backoff = -1;							
						} else {
							backoff = 
							refstr[gostart-1-MIN_gap_size-move-matchcount];
						}
						
						if (backoff < 2){
							stop = true;//stop this iteration
							flager = false;//stop this ref position
						} else if(backoff == 
								qrysoffsettok[qryscount+tokindex+longest_len_start-1-matchcount]){
							matchcount++;
						} else {
							stop = true;
						}						
					}

					if(matchcount == longest_len_start
						&& 
						checkBoundaryGap(gostart-1-MIN_gap_size-move+1, 
							gostart-1, 
							L_tar, 
							R_tar, 
							RLP)){
						nextpos = atomicAdd(counter, 1);
						oneGapSA[nextpos].position	= blockId;
						oneGapSA[nextpos].str_position = gostart-1-MIN_gap_size-move-longest_len_start+1;
						oneGapSA[nextpos].length = longest_len_end+MIN_gap_size+move+longest_len_start-1;
						//to the end aXb's b	 
					}
				} 

				move++;
				if(longest_len_start+MIN_gap_size+move + longest_len_end > MAX_rule_span){
					flager = false;
				}
			}
		}
		threadx += blockDim.x;
	}

	}

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
			uint8_t* R_tar) {

				unsigned int twoBlockId = blockIdx.y*gridDim.x+blockIdx.x;
				if (twoBlockId >= distinctTwoGapCount){
					return;
				}
				unsigned int oneBlockId = twoGapSearch[twoBlockId].blockid;
				if (oneBlockId >= distinctOneGapCount){
					printf("Not possible inside twoGaponSa kernel! oneId %u | distinctOneGapCount %u\n", oneBlockId, distinctOneGapCount);
					return;
				}

				//No precomputation here
				int startSA = oneGapSearch[oneBlockId].start_on_salist;
				int endSA = oneGapSearch[oneBlockId].end_on_salist;

				if(startSA == -1 && endSA == -1){
					return;
				}

				if (endSA > oneGapSACount || endSA < startSA){
					printf("not possible1 in kernel two gap sa1| endSA %d | oneGapSACount %u | startSA %d\n", endSA, oneGapSACount, startSA);
					return;
				}	

				int threadx = threadIdx.x;
				int dis = endSA - startSA+1;
				int gostart = 0;
				int move = 0;
				bool flager = true;
				int temp = -1;

				//for aXbXc. Add them all minus 1. aXbXc's c starting position.
				//watch for c's longest_end_range
				/*int search_tokindex = oneGapSearch[oneBlockId].qrystart
				  +oneGapSearch[oneBlockId].qrystart_len
				  +oneGapSearch[oneBlockId].gap
				  +oneGapSearch[oneBlockId].qryend_len - 1
				  +twoGapSearch[twoBlockId].gap2;*/
				int search_tokindex = twoGapSearch[twoBlockId].gap2;
				int longest_len_end = twoGapSearch[twoBlockId].qryend_len;
				
				if(longest_len_end!=1){
					printf("Should be one!!\n");
					return;
				}
				
				int nextpos = -1;
				int precomputationIndex = -1;
				int tobesearch_start = 0;

				unsigned int precompstart;
				uint8_t precomplen;
				
				int matchcount;
				int backoff;
				bool stop=false;
				
				int preCache = qrysoffsettok[qryscount+search_tokindex];
				///for debug purpose				
				if (preCache < 2){
					printf("not possible in kernel two gap on SA! query checking wrong\n");
					return;
				}
				///end debug
				
				if(dis == 1 && oneGapSA[startSA].length==0){

					//Precomputation case
					precomputationIndex = oneGapSA[startSA].str_position;
					dis = precomp_index[precomputationIndex].end - 
						precomp_index[precomputationIndex].start + 1;
					tobesearch_start = precomp_index[precomputationIndex].start;

					//check
					if(oneGapSearch[oneBlockId].qrystart_len!=1
							|| oneGapSearch[oneBlockId].qryend_len!=1
							|| precomputationIndex < 0){
						printf("Good to know! Two gap on SA kernel\n");
						return;
					}					

					//precomputation case
					while(threadx < dis){
						precompstart =  precomp_onegap[tobesearch_start+threadx].start;
						precomplen = precomp_onegap[tobesearch_start+threadx].length;
						//last position, not real length
						gostart = precompstart + precomplen;					
						//gostart = oneGapSA[startSA+threadx].str_position+
						//	oneGapSA[startSA+threadx].length;
						move = 0;
						temp = -1;
						flager = true;
						while(flager){
							if(move == 0){
								//Check gap, if the gap is a <2 between two sentences, 
								//then skip
								temp = refstr[gostart+MIN_gap_size];
								if(temp < 2){
									flager = false;
								}
							}
							temp = refstr[gostart+1+MIN_gap_size+move];

							if (precomplen+
									1 + MIN_gap_size+move + 1>MAX_rule_span){
								flager = false;
							}
							
							if (temp < 2 ){
								flager = false;
							} else if (flager && temp == preCache){
								matchcount = 1;
								stop = false;
								while(!stop&&matchcount < longest_len_end
										&& precomplen+matchcount+MIN_gap_size+move + 
										1 + 1 <= MAX_rule_span ){
									////debugging
									if (qrysoffsettok[qryscount+search_tokindex+matchcount]< 2){
										printf("This is not possible - onegap kernel");
										return;
									}
									//debugging
									
									backoff = refstr[gostart+1+MIN_gap_size+move+matchcount];
									if (backoff < 2){
										stop = true;
										flager = false;
									} else if(backoff == 
											qrysoffsettok[qryscount+search_tokindex+matchcount]){
										matchcount++;
									} else {
										stop = true;
									}
								}

								if(matchcount == longest_len_end
									&& checkBoundaryGap(precompstart+precomplen+1, 
											precompstart+1+precomplen+MIN_gap_size+move-1, 
											L_tar, 
											R_tar, 
											RLP)){
									nextpos = atomicAdd(counter, 1);
									twoGapSA[nextpos].position	= twoBlockId;
									twoGapSA[nextpos].str_position = precompstart;
									twoGapSA[nextpos].length = precomplen;// to the end aXb's b	
									twoGapSA[nextpos].length2 = 
										precomplen+
										1+MIN_gap_size+move+longest_len_end-1; 
									//to the end of aXbXc's c, end of c if c has multiple characters.
								}
							} 			
							move++;
						}	
						threadx+=blockDim.x;
					}

				} else {/*if(dis!=1&& oneGapSA[startSA].length!=0){*/
					//None precomputation case
					while(threadx < dis){
						precompstart = oneGapSA[startSA+threadx].str_position;
						precomplen = oneGapSA[startSA+threadx].length;
						///check debug
						if(precomplen==0){
							printf("Not possible!!! inside kernel two gap on sa - length==0\n");
							return;							
						}
						/////				
						gostart = precompstart+precomplen;//last position, not real length
						move = 0;
						temp = -1;
						flager = true;
						while(flager){
							if(move == 0){
								//Check gap, if the gap is a <2 between two sentences, 
								//then skip
								temp = refstr[gostart+MIN_gap_size];
								if(temp < 2){
									flager = false;
								}
							}
							temp = refstr[gostart+1+MIN_gap_size+move];

							if (precomplen+
									1 + MIN_gap_size+move + 1>MAX_rule_span){
								flager = false;
							}
							
							if (temp < 2 ){
								flager = false;
							} else if (flager && temp == preCache){
								matchcount = 1;
								stop = false;
								while(!stop&&matchcount < longest_len_end
										&& precomplen+matchcount+MIN_gap_size+move + 
										1 + 1 <= MAX_rule_span ){
									////debugging
									if (qrysoffsettok[qryscount+search_tokindex+matchcount]< 2){
										printf("This is not possible - onegap kernel");
										return;
									}
									//debugging
									
									backoff = refstr[gostart+1+MIN_gap_size+move+matchcount];
									if (backoff < 2){
										stop = true;
										flager = false;
									} else if(backoff == qrysoffsettok[qryscount+search_tokindex+matchcount]){
										matchcount++;
									} else {
										stop = true;
									}
								}

								if(matchcount == longest_len_end
									&& checkBoundaryGap(precompstart+precomplen+1, 
											precompstart+1+precomplen+MIN_gap_size+move-1, 
											L_tar, 
											R_tar, 
											RLP)){
									nextpos = atomicAdd(counter, 1);
									twoGapSA[nextpos].position	= twoBlockId;
									twoGapSA[nextpos].str_position = precompstart;
									twoGapSA[nextpos].length = precomplen;// to the end aXb's b	
									twoGapSA[nextpos].length2 = 
										precomplen+
										1+MIN_gap_size+move+longest_len_end-1; 
									//to the end of aXbXc's c, end of c if c has multiple characters.
								}
							} 			
							move++;
						}	
						threadx+=blockDim.x;
					}
				}
			}


	__global__ void precomp(
			precompute_enu* oneGapPrecomp,
			int* refstr, 
			int *refsa,
			precompute_enu_2* onegap_precomp,
			unsigned int* counter,
			unsigned int* RLP,
			uint8_t* L_tar,
			uint8_t* R_tar,
			int* featureMissingCount_d) {
		//max_rule_span -> 15; Maximum rule span
		//max_rule_symbols -> 5; Maximum number of symbols (terminals + nontermals) in a rule
		//min_gap_size -> 1; Minimum gap size
		//max_nonterminals -> 2; Maximum number of nonterminals in a rule
		//max_phrase_len -> 4; Maximum frequent phrase length

		int index = blockIdx.x+gridDim.x*blockIdx.y;

		if(index < PRECOMPUTECOUNT*PRECOMPUTECOUNT){	
			featureMissingCount_d[index]=0;
			/*if(index == PRECOMPUTECOUNT*PRECOMPUTECOUNT-1){
			  printf("Got it!\n");
			  }*/
			int tid = threadIdx.x + oneGapPrecomp[index].start;
			//printf("TID %d\n", index);
			unsigned int end = 0;

			unsigned int nextpos;
			int temp;
			int move=0;
			bool flager = true;
			int gostart = -1;

			end = oneGapPrecomp[index].length + oneGapPrecomp[index].start;
			bool reverser = oneGapPrecomp[index].reverse;

			while(tid < end){
				move = 0;
				flager = true;

				if (reverser){
					gostart = refsa[tid];
					if(refstr[gostart]!=oneGapPrecomp[index].token_a  ||
							refstr[refsa[oneGapPrecomp[index].start]] != refstr[gostart]){
						printf("Not possible!!!!!! %d %d\n", refstr[gostart], oneGapPrecomp[index].token_a);
						return;
					}
					while(flager){
						if(move == 0){
							temp = refstr[gostart+MIN_gap_size];
							if(temp < 2){
								flager = false;
							}
						}
						temp = refstr[gostart+1+MIN_gap_size+move];
						if(temp < 2 ){
							flager = false;
						}
						else if (flager && temp == oneGapPrecomp[index].token_b){
							if(checkBoundaryGap(gostart+1, 
									gostart+move+1+MIN_gap_size-1, 
									L_tar, 
									R_tar, 
									RLP)){
								nextpos = atomicAdd(counter, 1);
								if(nextpos > ONEGAP_PRECOMPUT_PREALLOCATION - 1000){
									printf("Stop!!!\n");
									return;
								}
								onegap_precomp[nextpos].index = index;
								onegap_precomp[nextpos].start = gostart;
								//string array position
								//The last position of aXb's b. inclusive
								onegap_precomp[nextpos].length = move+1+MIN_gap_size;
							} else {
								atomicAdd(&(featureMissingCount_d[index]), 1);
							}
						}
						move++;
						if(1+MIN_gap_size+move + 1> MAX_rule_span){
							flager = false;
						}
					}					
				} else {
					gostart = refsa[tid];
					if(refstr[gostart]!=oneGapPrecomp[index].token_b ||
							refstr[refsa[oneGapPrecomp[index].start]] != refstr[gostart]){
						printf("Not possible!!!!!! %d %d\n", refstr[gostart], oneGapPrecomp[index].token_a);
					}
					while(flager){
						if(move == 0 && gostart-MIN_gap_size >= 0){
							temp = refstr[gostart-MIN_gap_size];
							if(temp < 2){
								flager = false;
							}
						}
						if (flager&&gostart-1-MIN_gap_size-move >= 0){						
							temp = refstr[gostart-1-MIN_gap_size-move];
							if(temp < 2 ){
								flager = false;
							}
							else if (flager && temp == oneGapPrecomp[index].token_a){
								if(checkBoundaryGap(gostart-1-MIN_gap_size-move+1, 
										gostart-1, 
										L_tar, 
										R_tar, 
										RLP)){
									nextpos = atomicAdd(counter, 1);
									onegap_precomp[nextpos].index = index;
									onegap_precomp[nextpos].start = gostart-1-MIN_gap_size-move;
									onegap_precomp[nextpos].length = move+1+MIN_gap_size;
								} else {
									atomicAdd(&(featureMissingCount_d[index]), 1);
								}
							}
						} else {
							flager = false;
						}

						move++;
						if(1+MIN_gap_size+move + 1> MAX_rule_span){
							flager = false;
						}
					}
				}
				tid+= blockDim.x;
			}		
			//index += gridDim.x;
			//printf("gridDim %d\n",gridDim.x);
		}	
	}
