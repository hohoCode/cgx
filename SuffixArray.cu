#include "ComTypes.h"
#include "Timer.h"
#include <thrust/device_vector.h>
#include <thrust/sort.h>
//#include <thrust/host_vector.h>
#include <thrust/version.h>

#include "GappyLook.h"

#define ISSPACE(x) ((x) == ' ' || (x) == '\t' || (x) == '\n')
#define cache 10
#define ISPUNCSPA(x) ((x) == ' ' || (x) == '\t' || (x) == '.' || (x) == '?'  || (x) == ','  || (x) == '!'  || (x) == ':'  || (x) == ';'  || (x) == '}'  || (x) == ']'  || (x) == ')'  || (x) == '*'  || (x) == '@' || (x) == '\n' || (x) == '"' || (x) == '\'')

#define COMP1(s1, s2, newb, mymatchcount, qrylimit, ok)		\
	mymatchcount = 0;\
ok = 0;\
do {									\
	int ref = *(newb+ s1 + mymatchcount);				\
	int qry = *(s2 + mymatchcount); 					\
	while (mymatchcount < qrylimit && (ref ==qry) && ref != 1 && qry != -1) {\
		mymatchcount++;\
		ref = *(newb+ s1 + mymatchcount);\
		qry = *(s2 + mymatchcount);\
	}\
	if (qry==-1 || mymatchcount == qrylimit){\
		ok = 1;\
	}\
	/*printf("str %c - q %c - myrc %d - count %d - saind %d\n",a, b, myrc, mymatchcount, s1);*/\
} while(0)

	struct twoGapEnumerationCompare
	{
		__host__ __device__
			bool operator()(const gapPattern2& a, const gapPattern2& b) const
			{
				if(a.blockid != b.blockid){
					return (a.blockid < b.blockid);
				} else if(a.number != b.number){
					return (a.number<b.number); 
				} else {
					for(int i=0; i<b.number;i++){
						if(a.pattern[i]!=b.pattern[i]){
							return a.pattern[i] < b.pattern[i];
						}
					}
				}
				return true;
			}
	};

	struct oneGapEnumerationCompare
	{
		__host__ __device__
			bool operator()(const gapPattern1& a, const gapPattern1& b) const
			{
				if(a.number != b.number){
					return (a.number<b.number); 
				} else {
					for(int i=0; i<b.number;i++){
						if(a.pattern[i]!=b.pattern[i]){
							return a.pattern[i] < b.pattern[i];
						}
					}
				}
				return true;
			}
	};

	struct oneGapSACompare
	{
		__host__ __device__
			bool operator()(const gappytyp_sa& a, const gappytyp_sa& b) const
			{
				return (a.position < b.position) ||((a.position == b.position) 
					&& a.str_position < b.str_position);
			}
	};

	struct twoGapSACompare
	{
		__host__ __device__
			bool operator()(const gappytyp_sa2& a, const gappytyp_sa2& b) const
			{
				return (a.position < b.position)||
					((a.position == b.position) && a.str_position < b.str_position);
			}
	};

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

__global__ void suffixArrayFindConnectionTwoWayTDI(
		int* refstr, 
		int *refsa,
		/*int *reflcp,*/
		int * reflcpleft,
		int * reflcpright,
		int reflen,
		int *qrysoffsettok,
		int qryscount,    
		int tokenscount,
		result_t_two* result,
		int* connectoffset,
		int totalconnect,
		result_t* result_connect) {

	int tokindex = threadIdx.x + blockDim.x * blockIdx.x;
	if(tokindex >= tokenscount){
		return;
	}
	int longestmatch = result[tokindex].longestmatch ;
	if (longestmatch <=1){
		return;
	}
	int start = connectoffset[tokindex];		
	if(start < 0 ){
		return;
	}
	int seq = blockIdx.y/2;
	int up = blockIdx.y%2;		

	int match = 2 + seq;
	start += seq;
	if (match > longestmatch){
		return;
	}
	int L, R, M, Rlcp = 0, skip = 0, Llcp = 0;
	int startREF = -1;
	int temp = -1;
	int longest = -1;
	int longlen = -1;
	int firstfindhit = -1;
	int firstfindhitL = -1;		
	int firstfindhitR = -1;		
	int firstfindhitlen = -1;
	int holdtemp = -1;		
	int foundexactlcp = 0;
	int a, b;
	int* query = qrysoffsettok + qryscount+ tokindex;
	int LL = result[tokindex].firstfindhitL;
	int MM = result[tokindex].firstfindhit;
	int RR = result[tokindex].firstfindhitR;

	while(match <= longestmatch){
		L = LL;
		R = RR;
		//////////////////////////////////////////////////////////////////////////////////////
		startREF = -1;
		temp = -1;
		holdtemp = -1;
		longlen = 0;
		foundexactlcp = 0;
		firstfindhit = -1;
		firstfindhitL = -1;
		firstfindhitR = -1;
		Llcp = 0;
		Rlcp = 0;
		while (R - L > 1) {
			longlen = 0;
			if(L == LL && R == RR){
				M = MM;
			} else {
				M = (L + R) >> 1;
			}
			//printf("#Tokindex %d - L %d - R %d - M%d\n", tokindex, L, R, M);
			if (Llcp >= Rlcp){
				longlen = Llcp;
				if (L == M-1){
					skip = reflcpleft[M]; 
				} else{
					holdtemp = (L + M)>>1;
					skip = fminf(reflcpleft[holdtemp], reflcpright[holdtemp]);
				}

				if(longlen < skip){
					L = M;
				} else if (longlen > skip) {
					R = M;
					Rlcp = skip;
				} else {					
					//printf("Comparison LlCP/SKIP %d\n", skip);
					startREF = refsa[M] + longlen;
					a = *(query + longlen);
					b = refstr[startREF];
					if (a == -1) {
						foundexactlcp = 1;
						printf("THIS IS NOT POSSIBLE!!!!\n");                                             	
						break;
					}
					//printf("TBD Llcp %d - REF %d - Tokindex %d - Qryindex %d - qry%c - str%c\n", longlen, startREF, tokindex, qryindex, a, b);							
					if(a!= -1 && b!=1 ){
						temp = a - b;
						while (a != -1 && b !=1 && temp == 0) {									
							longest = M;
							//printf("Middle Matched  Llcp %d - REF %d - Tokindex %d - Qryindex %d - qry%c - str%c\n", longlen, startREF, tokindex, qryindex, a, b);
							longlen++;
							/*if (longlen >= toklen - tokindex) {									
							  foundexactlcp = 1;
							  break;
							  }*/
							startREF++;

							if(firstfindhit == -1 && M >=0 && longlen >= match){
								firstfindhit = M;
								firstfindhitL = L;
								firstfindhitR = R;
								firstfindhitlen = longlen;
								foundexactlcp = 1;
								break;
							}
							a = *(query + longlen);
							b = refstr[startREF];
							if (a == -1) {									
								foundexactlcp = 1;
								printf("THIS IS NOT POSSIBLE!!!!\n");          
								break;
							}
							if(a != -1 && b !=1){
								temp = a - b;
							}
						}
						if (foundexactlcp == 1){
							//printf("Can you believe it? Tokindex%d - match%d\n", tokindex, match);
							break;
						}
					}

					if (a == -1){
						R = M;
						L = M;//break;
					} else if (b ==1) {
						L = M;
						Llcp = longlen;
					} else if (temp > 0){
						L = M;
						Llcp = longlen;
					} else {
						R = M;
						Rlcp = longlen;
					}
				}

			}else {
				longlen = Rlcp;
				if (R == M+1){
					skip = reflcpright[M]; 
				} else{
					holdtemp = (R + M)>>1;
					skip = fminf(reflcpleft[holdtemp], reflcpright[holdtemp]);
				}

				//printf("Rlcp Iteration - R %d - L %d - M %d - Llcp %d - Rlcp %d - M'M %d - TID %d - Qryindex %d\n", R, L, M, Llcp, Rlcp, skip, tokindex, qryindex);
				if(longlen < skip){
					R = M;                        
				} else if (longlen > skip) {
					L = M;
					Llcp = skip;
				} else {					
					//printf("Comparison LlCP/SKIP %d\n", skip);
					startREF = refsa[M] + longlen;	
					a = *(query + longlen);
					b = refstr[startREF];
					if (a == -1) {									
						foundexactlcp = 1;
						break;
					}
					if(a!= -1 && b!=1){
						temp = a - b;					
						while (a != -1 && b !=1 && temp == 0) {									
							longest = M;
							//printf("Middle Matched Ref = %c - Llcp %d - REF %d - tokindex %d - qryindex %d\n", *(refstr+startREF), longlen, startREF, tokindex, qryindex);
							longlen++;
							startREF++;

							if(firstfindhit == -1 && M >=0 && longlen >= match){
								firstfindhit = M;
								firstfindhitL = L;
								firstfindhitR = R;
								firstfindhitlen = longlen;
								foundexactlcp = 1;
								break;
							}
							a = *(query + longlen);
							b = refstr[startREF];
							if (a == -1) {									
								foundexactlcp = 1;
								printf("THIS IS NOT POSSIBLE!!!!\n"); 
								break;
							}	
							if(a != -1  && b !=1){
								temp = a - b;
							}
						}
						if (foundexactlcp == 1){
							//printf("Can you believe it? Tokindex%d - match%d\n", tokindex, match);
							break;
						}
					}

					if (a == -1){
						R = M;
						L = M;//break;
					} else if (b ==1) {
						L = M;
						Llcp = longlen;
					} else if (temp > 0){
						L = M;
						Llcp = longlen;
					} else {
						R = M;
						Rlcp = longlen;
					}
				}

			}
		}

		if (firstfindhit == -1){
			printf("Firstfindhit is -1 - THIS NOT POSSIBLE!!!!\n %d - not possible!!!! longestmatch %d - totalcon %d - match%d\n", tokindex, longestmatch, totalconnect, match); 
			break;
		}

		if (firstfindhit != -1 && longlen > 0 && foundexactlcp == 1){
			/// Here we have found something, this is the final step
			/// We are going to locate the upper bound, and lower bound - Final step
			///Now I have L R and longest.					
			///Deal with upper and lower bound!
			longest = firstfindhit;			
			//printf("QID %d - LongLEN %d - FirstFINDHIT %d\n", qryindex, longlen, longest); 
			if (up == 1){
				R = firstfindhit;
				L = firstfindhitL;
				while (R - L > 1){
					M = (L + R) >> 1;
					holdtemp = (R + M)>>1;
					if (R == M+1){
						skip = reflcpright[M]; 
					} else{
						skip = fminf(reflcpleft[holdtemp], reflcpright[holdtemp]); //LCP(M, M')
					}
					//printf("!! Matched= R%d - L%d - M%d - LCP(R,M)%d\n", R, L, M, skip);
					if (skip >= match){
						longest = M;
						R = M;
						firstfindhitlen = skip;
					} else {
						L = M;
					}
				}			   	
				////Write results
				result_connect[start].up = longest;
			} else {
				R = firstfindhitR;
				L = firstfindhit;
				while (R - L > 1){
					M = (L + R) >> 1;
					holdtemp = (L + M)>>1;

					if (L == M-1){
						skip = reflcpleft[M]; 
					} else{
						skip = fminf(reflcpleft[holdtemp], reflcpright[holdtemp]); //LCP(M, M')
					}

					//printf("!! Matched= R%d - L%d - M%d - LCP(L,M)%d\n", R, L, M, skip);
					if (skip >= match){
						longest = M;
						L = M;
						firstfindhitlen = skip;
					} else {
						R = M;
					}
				}
				////Write results
				result_connect[start].down = longest;
			}
		}

		//////////////////////////////////////////////////////////////////////////////////////
		match += gridDim.y/2;
		start += gridDim.y/2;
	}
}

__global__ void suffixArrayFindLwRwKernelTwoWayTDI(
		int* refstr, 
		int *refsa,
		/*int *reflcp,*/
		int * reflcpleft,
		int * reflcpright,
		int reflen,
		int *qrysoffsettok,
		int qryscount,    
		int tokenscount,
		result_t_two*result) {

	int qryindex = threadIdx.y + blockDim.y * blockIdx.y;
	int tokindex = threadIdx.x;

	if (qryindex >= qryscount ){
		return;
	}

	int queryoffset = qrysoffsettok[qryindex];
	int toklen = -1;
	if (qryindex == qryscount - 1){
		toklen = tokenscount - queryoffset;
	} else {
		toklen = qrysoffsettok[qryindex+1] - queryoffset;
	}
	int up = blockIdx.x % 2;

	//printf("NOeeD %d - toklen %d - QID %d\n", tokindex, toklen, qryindex);
	if (tokindex >= toklen){
		//printf("EXCEED %d - toklen %d - QID %d\n", tokindex, toklen, qryindex);
		return;
	}


	//char *query = &(qrys[qrysoffset[qryindex]]);
	int L, R, M, Rlcp = 0, skip = 0, Llcp = 0;
	int startREF = -1;
	int temp = -1;
	int longest = -1;
	int firstfindhit = -1;
	int firstfindhitL = -1;		
	int firstfindhitR = -1;		
	int firstfindhitlen = -1;
	int holdtemp = -1;
	int longlen = -1;
	int foundexactlcp = 0;
	int ok = 0;
	int a, b;
	result += queryoffset;		

	//while ( tokindex < toklen ){
	//__syncthreads();
	if(result[tokindex].longestmatch >0&& up == 1){
		result[tokindex].longestmatch = 0;
	}
	//printf("tokindex %d - toklen %d - iteration %d - qryindex %d - reflen %d\n", tokindex, toklen, iteration, qryindex, reflen);
	a = qrysoffsettok[qryscount+ queryoffset + tokindex];  
	if(a == -1){
		return;
		//goto END;
		/*iteration++;
		  tokindex += iteration*blockDim.x;
		  continue;*/
	}

	L = 0;
	R = reflen - 1;
	Llcp = 0;
	Rlcp = 0;

	/**First check the reference bundary see whether there is a match or not*/
	/**Starting from Right side */
	foundexactlcp = 0;			
	longest = -1;
	firstfindhit = -1;			
	firstfindhitL = -1;		
	firstfindhitR = -1;		
	firstfindhitlen = -1;

	longlen = -1;
	ok = 0;
	COMP1(refsa[R], qrysoffsettok +qryscount+ queryoffset + tokindex , refstr, Rlcp, toklen-tokindex, ok);//Ref -rc- Query
	if (Rlcp >0 && ok == 1){			
		foundexactlcp = 1; ///Found longest!
		longest = R;
		longlen = Rlcp;
		if (up != 1){
			result[tokindex].longestmatch = Rlcp;
			//result[tokindex].longestind = R;//First Longest
		}
	}
	//printf("RIGHT UP %d - up/down %d - Rlcp %d - dingdong %d - exactlcp %d - qryindex %d\n", up, R, Rlcp, dingdong, exactlcp, qryindex);

	if (Rlcp > 0){
		firstfindhit = R;
		firstfindhitL = L;
		firstfindhitR = R;
		firstfindhitlen = Rlcp;
		if (up != 1){
			result[tokindex].down = R;
			//result[tokindex].lcpdown = Rlcp;					
			//goto END;
			return;
			//iteration++;
			//tokindex += iteration*blockDim.x;
			//printf("UP %d - up/down %d - Rlcp %d\n", up, R, Rlcp);
			//continue;
		} else {
			result[tokindex].up = R;
			//result[tokindex].lcpup = Rlcp;
		}
	}

	/****Now move to the left side, check whether there is a match or not*/
	ok = 0;			

	/***If we have found the longest match in the bundary, we do not need to find another one in the middle*/
	/***If we have not found any match in the two bundary, we need cotinue to find a hit in the middle*/
	/***Once there is a find on the longest match, we break and move to the next step*/
	if(foundexactlcp == 0 ){
		// Look for Lw
		/** This is using binary seearch for finding a longest match in the middle*/
		/**Binary seach with non-standard LCP*/
		startREF = -1;
		temp = -1;
		holdtemp = -1;
		longlen = 0;
		int* query = qrysoffsettok + qryscount+ queryoffset + tokindex;  
		while (R - L > 1) {
			longlen = 0;
			M = (L + R) >> 1;
			//printf("#Tokindex %d - L %d - R %d - M%d\n", tokindex, L, R, M);
			if (Llcp >= Rlcp){
				longlen = Llcp;
				if (L == M-1){
					skip = reflcpleft[M]; 
				} else{
					holdtemp = (L + M)>>1;
					skip = fminf(reflcpleft[holdtemp], reflcpright[holdtemp]);
				}

				//printf("Llcp Iteration - R %d - L %d - M %d - Llcp %d - Rlcp %d - MM' %d - TID %d - qrydindex %d\n", R, L, M, Llcp, Rlcp, skip, tokindex, qryindex);
				if(longlen < skip){
					L = M;
				} else if (longlen > skip) {
					R = M;
					Rlcp = skip;
				} else {					
					//printf("Comparison LlCP/SKIP %d\n", skip);
					startREF = refsa[M] + longlen;
					a = *(query + longlen);
					b = refstr[startREF];
					if (longlen >= toklen - tokindex || a == -1) {
						foundexactlcp = 1;
						//printf("BREAK QID %d TID %d\n", qryindex, tokindex);                                    
						break;
					}
					//printf("TBD Llcp %d - REF %d - Tokindex %d - Qryindex %d - qry%c - str%c\n", longlen, startREF, tokindex, qryindex, a, b);							
					if(a!= -1 && b!=1 ){
						temp = a - b;
						while (a != -1 && b !=1 && temp == 0) {									
							longest = M;
							//printf("Middle Matched  Llcp %d - REF %d - Tokindex %d - Qryindex %d - qry%c - str%c\n", longlen, startREF, tokindex, qryindex, a, b);
							longlen++;
							startREF++;

							if(firstfindhit == -1 && M >=0){
								firstfindhit = M;
								firstfindhitL = L;
								firstfindhitR = R;
							}
							if (firstfindhit == M){
								firstfindhitlen = longlen;
							}
							//printf("Middle Matched  Llcp %d - REF %d - Tokindex %d - Qryindex %d - qry%c - str%c - 2nd - longlen %d - exactlcp %d\n", longlen, startREF, tokindex, qryindex, a, b, longlen, exactlcp);
							if (longlen >= toklen - tokindex) {									
								foundexactlcp = 1;
								//printf("BREAK QID %d TID %d\n", qryindex, tokindex);                                    
								break;
							}
							a = *(query + longlen);
							b = refstr[startREF];
							if (a == -1){
								foundexactlcp = 1;									
								break;
							}	
							if(a != -1 && b !=1){
								temp = a - b;
							}
						}
						if (foundexactlcp == 1){
							break;
						}
					}

					if (a == -1){
						R = M;
						L = M;//break;
					} else if (b ==1) {
						L = M;
						Llcp = longlen;
					} else if (temp > 0){
						L = M;
						Llcp = longlen;
					} else {
						R = M;
						Rlcp = longlen;
					}
				}

			}else {
				longlen = Rlcp;
				if (R == M+1){
					skip = reflcpright[M]; 
				} else{
					holdtemp = (R + M)>>1;
					skip = fminf(reflcpleft[holdtemp], reflcpright[holdtemp]);
				}

				//printf("Rlcp Iteration - R %d - L %d - M %d - Llcp %d - Rlcp %d - M'M %d - TID %d - Qryindex %d\n", R, L, M, Llcp, Rlcp, skip, tokindex, qryindex);
				if(longlen < skip){
					R = M;                        
				} else if (longlen > skip) {
					L = M;
					Llcp = skip;
				} else {					
					//printf("Comparison LlCP/SKIP %d\n", skip);
					startREF = refsa[M] + longlen;	
					a = *(query + longlen);
					b = refstr[startREF];
					if (longlen >= toklen - tokindex || a == -1) {									
						foundexactlcp = 1;
						//printf("BREAK QID %d TID %d\n", qryindex, tokindex);                                    
						break;
					}
					if(a!= -1 && b!=1){
						temp = a - b;					
						while (a != -1 && b !=1 && temp == 0) {									
							longest = M;
							//printf("Middle Matched Ref = %c - Llcp %d - REF %d - tokindex %d - qryindex %d\n", *(refstr+startREF), longlen, startREF, tokindex, qryindex);
							longlen++;
							startREF++;

							if(firstfindhit == -1 && M >=0){
								firstfindhit = M;
								firstfindhitL = L;
								firstfindhitR = R;
							}
							if (firstfindhit == M){
								firstfindhitlen = longlen;
							}
							if (longlen >= toklen - tokindex) {									
								foundexactlcp = 1;
								//printf("BREAK QID %d TID %d\n", qryindex, tokindex);                                    
								break;
							}				
							a = *(query + longlen);
							b = refstr[startREF];
							if (a == -1){
								foundexactlcp = 1;									
								break;
							}				
							if(a != -1  && b !=1){
								temp = a - b;
							}
						}
						if (foundexactlcp == 1){
							break;
						}
					}

					if (a == -1){
						R = M;
						L = M;//break;
					} else if (b ==1) {
						L = M;
						Llcp = longlen;
					} else if (temp > 0){
						L = M;
						Llcp = longlen;
					} else {
						R = M;
						Rlcp = longlen;
					}
				}

			}
		}
	}

	if(longlen >0 && up == 1){
		result[tokindex].longestmatch = longlen;
		//result[tokindex].longestind = longest;
	}
	if (firstfindhit == -1 && longlen >  0){
		printf("THIS FIRSTHIT NOT POSSIBLE!\n");
	}

	if (firstfindhit != -1 && longlen > 0){
		/// Here we have found something, this is the final step
		/// We are going to locate the upper bound, and lower bound - Final step
		///Now I have L R and longest.					
		///Deal with upper and lower bound!
		longest = firstfindhit;
		if(up == 1){
			result[tokindex].firstfindhit = firstfindhit;
			result[tokindex].firstfindhitL = firstfindhitL;
			result[tokindex].firstfindhitR = firstfindhitR;
		}

		//printf("QID %d - LongLEN %d - FirstFINDHIT %d\n", qryindex, longlen, longest); 
		if (up == 1){
			R = firstfindhit;
			L = firstfindhitL;
			while (R - L > 1){
				M = (L + R) >> 1;
				holdtemp = (R + M)>>1;
				if (R == M+1){
					skip = reflcpright[M]; 
				} else{
					skip = fminf(reflcpleft[holdtemp], reflcpright[holdtemp]); //LCP(M, M')
				}
				//printf("!! Matched= R%d - L%d - M%d - LCP(R,M)%d\n", R, L, M, skip);
				if (skip >= 1){
					longest = M;
					R = M;
					//firstfindhitlen = skip;
				} else {
					L = M;
				}
			}			   	
			////Write results
			result[tokindex].up = longest;
			//result[tokindex].lcpup = firstfindhitlen;				
		} else {
			R = firstfindhitR;
			L = firstfindhit;
			while (R - L > 1){
				M = (L + R) >> 1;
				holdtemp = (L + M)>>1;

				if (L == M-1){
					skip = reflcpleft[M]; 
				} else{
					skip = fminf(reflcpleft[holdtemp], reflcpright[holdtemp]); //LCP(M, M')
				}

				//printf("!! Matched= R%d - L%d - M%d - LCP(L,M)%d\n", R, L, M, skip);
				if (skip >= 1){
					longest = M;
					L = M;
					//firstfindhitlen = skip;
				} else {
					R = M;
				}
			}	

			////Write results
			result[tokindex].down = longest;
			//result[tokindex].lcpdown = firstfindhitlen;
		}

	}

}

void * suffixArraySearchInit(int refbufsize, int qrysbufsize, int connectsize, int sasize) {
	cudaError_t err;
	
	size_t freeMem = 0;
	size_t totalMem = 0;
	size_t allocMem = 0;

	cudaMemGetInfo(&freeMem, &totalMem);  
	fprintf(stderr, "Memory avaliable: Free: %lu, Total: %lu. Start of Suffix.cu\n",freeMem/(1024*1024), totalMem/(1024*1024));

	context_t * ctx = (context_t *)malloc(sizeof(context_t));

	/* allocate the queries buffer */   
	cudaMalloc((void**)&(ctx->qryset_d.qrysbuf), qrysbufsize);

	/* allocate space for the suffix array data (SA, LCP),
	   the reference itself will be allocated separately */
	cudaMalloc((void**)&(ctx->ref_d.buf), refbufsize);
	cudaMalloc((void**)&(ctx->ref_d.sa), sasize);
	cudaMalloc((void**)&(ctx->ref_d.str), sasize);
	cudaMalloc((void**)&(ctx->qryset_d.connectoffset), connectsize);
	return (void *)ctx;
}

void suffixArraySearchFinalize_One(void *handler) {

	context_t * ctx = (context_t *)handler;

	cudaFree(ctx->qryset_d.qrysbuf);
	cudaFree(ctx->qryset_d.result_two);	
	cudaFree(ctx->qryset_d.connectoffset);
	cudaFree(ctx->qryset_d.result_connect);
	cudaFree(ctx->qryset_d.onegap);
	//cudaFree(ctx->ref_d.buf);//lcp, rcp stuff
	//cudaFree(ctx->ref_d.str);
}


void suffixArraySearchFinalize(void *handler) {

	context_t * ctx = (context_t *)handler;	

	cudaFree(ctx->ref_d.sa);
	free(ctx);
}


__global__ void twoGapEnumeration(
		int* qrysoffsettok,
		int qryscount,    
		int tokenscount,
		result_t_two* qryresult,
		int totalconnect,
		int* tokindex_qryindex,
		gappy_search* searchGap,
		twogappy* twoGap,
		twoGapPattern* twogapPattern,
		unsigned int* counter,
		int distinctOneGapCount,
		gappy* oneGap,
		unsigned int oneGapEnuCount) {

	int searchPos = blockIdx.y*gridDim.x+blockIdx.x;
	if (searchPos >= distinctOneGapCount){
		return;
	}
	if(searchGap[searchPos].start_on_salist == -1||
			searchGap[searchPos].end_on_salist==-1){
		return;
	}

	int limit_symbol = MAX_rule_symbols - 1 - 1 - 
		searchGap[searchPos].qrystart_len - searchGap[searchPos].qryend_len;
	if(limit_symbol < -1){
		printf("Not possible symbol number in one gap kernel not ok! two Gap enu - limit_symbol %d; searchGap[searchPos].qrystart_len %d; searchGap[searchPos].qryend_len %d\n", limit_symbol, 
				searchGap[searchPos].qrystart_len,
				searchGap[searchPos].qryend_len);
		return;
	}		
	if(limit_symbol < 1){
		return;
	}

	unsigned int ender; //end positions on the onegap array
	if(searchPos == distinctOneGapCount - 1){
		ender = oneGapEnuCount;
		//Debug
		if(oneGapEnuCount <= searchGap[searchPos].position){
			printf("Not possible here Inside Two gap Enumeration\n");
			return;
		}
	} else {
		ender = searchGap[searchPos+1].position;
	}

	int threadx = threadIdx.x+searchGap[searchPos].position;
	int searchStart;
	int search_tokindex;

	int longest_len_end;
	int longest_end_iter = 1;
	unsigned int nextpos = 0;
	int i = 0;
	int qryindex;
	int end; //The first word of the next sentence

	while(threadx < ender){
		searchStart = oneGap[threadx].qrystart+	
			oneGap[threadx].qrystart_len+
			oneGap[threadx].gap+
			oneGap[threadx].qryend_len-1;
		search_tokindex =  searchStart+MIN_gap_size+1;

		//should minus the Min_gap_size
		if(searchStart <= tokenscount-1){
				//&& qrysoffsettok[qryscount+search_tokindex-1] != -1
				//&& qrysoffsettok[qryscount+search_tokindex] != -1){
			qryindex = tokindex_qryindex[searchStart];	
			if (qryindex != qryscount-1){
				end = qrysoffsettok[qryindex+1];
			} else {
				end = tokenscount;
			}

			while(search_tokindex < end ){//Note: && qrysoffsettok[qryscount+search_tokindex] != -1 check? - no!
				longest_len_end = qryresult[search_tokindex].longestmatch;
				longest_end_iter = 1;

				while(longest_end_iter<= limit_symbol
						&& longest_end_iter <= longest_len_end
						&& search_tokindex - oneGap[threadx].qrystart +  longest_end_iter - 1 <= MAX_rule_span_pattern){
					nextpos = atomicAdd(counter, 1);
					twoGap[nextpos].qryend_len = longest_end_iter;
					twoGap[nextpos].gap2 = search_tokindex;//Change it to unsigned int
					twoGap[nextpos].blockid =  searchPos;

					/*printf("Two Gap ENU|Pos %d - tokindex %d - qrystart_len %d - qryend_len1 %d - search_tok1 %d - qryend_len2 %d - search_tok2 %d|real %d|OneID %d\n", 		
						nextpos, oneGap[threadx].qrystart, oneGap[threadx].qrystart_len, oneGap[threadx].qryend_len, searchStart, longest_len_end, search_tokindex, qrysoffsettok[qryscount+search_tokindex], searchPos);*/
					///get the real aXb for sorting purpose.
					for(i=0; i < MAX_rule_symbols-4; i++){
						if(i < longest_end_iter){
							twogapPattern[nextpos].pattern[i] = 
								qrysoffsettok[qryscount+search_tokindex+i];							
						} else {
							twogapPattern[nextpos].pattern[i] = -2; 
						}
					}
					twogapPattern[nextpos].number = longest_end_iter;
					twogapPattern[nextpos].blockid = searchPos;
					longest_end_iter++; 			
				}			
				search_tokindex++; 
			}		
		}
		threadx += blockDim.x;
	}

}

__global__ void oneGapEnumeration(
		int* qrysoffsettok,
		int qryscount,    
		int tokenscount,
		result_t_two* qryresult,
		int totalconnect,
		int* tokindex_qryindex,
		gappy* onegap,
		gapPattern* onegapPattern,
		unsigned int* counter) {
	//max_rule_span -> 15; Maximum rule span
	//max_rule_symbols -> 5; Maximum number of symbols (terminals + nontermals) in a rule
	//min_gap_size -> 1; Minimum gap size
	//max_nonterminals -> 2; Maximum number of nonterminals in a rule
	//max_phrase_len -> 4; Maximum frequent phrase length

	int tokindex = threadIdx.x+blockIdx.y*blockDim.x;
	if(tokindex >= tokenscount-1){
		return;
	}
	int qryindex = tokindex_qryindex[tokindex];
	int end; //The first word of the next sentence

	if (qryindex != qryscount-1){
		end = qrysoffsettok[qryindex+1];
	} else {
		end = tokenscount;
	}

	if(qryindex >= qryscount 
			|| qryindex < 0 
			|| tokindex == end-1 
			|| tokindex == end-2){
		return;
	}
	int longest_len_start = qryresult[tokindex].longestmatch;
	int search_tokindex;
	int longest_len_end;
	int longest_end_iter = 1;
	int longest_start_iter = 1;
	unsigned int nextpos = 0;
	int i = 0;
	bool gappyYes = false; ///for print gappy pattern; bug fix; indicate whether it passes gaps or not
	uint8_t prePos = 0;

	/*if(qrysoffsettok[qryscount+tokindex]==-1){
		printf("Interesting MISTAKE inside one gap enu: longest match %d\n", longest_len_start);
		return;
	}*/
	
	while(longest_start_iter <= longest_len_start){
		search_tokindex = tokindex+longest_start_iter+MIN_gap_size;
		while(search_tokindex < end && search_tokindex - tokindex <= MAX_rule_span_pattern){
			if(qrysoffsettok[qryscount+search_tokindex]!=-1){
				longest_len_end = qryresult[search_tokindex].longestmatch;
				longest_end_iter = 1;

				while(longest_start_iter + 1 + longest_end_iter<= MAX_rule_symbols
						&& longest_end_iter <= longest_len_end
						&& search_tokindex - tokindex + longest_end_iter - 1<= MAX_rule_span_pattern){
					nextpos = atomicAdd(counter, 1);
					onegap[nextpos].qrystart = tokindex;
					onegap[nextpos].qrystart_len = longest_start_iter;
					onegap[nextpos].qryend_len = longest_end_iter;
					onegap[nextpos].gap = search_tokindex - tokindex - longest_start_iter;
					/*printf("One Gap ENU|Pos %d - tokindex %d - qrystart_len %d - qryend_len %d - search_to %d - gap %d\n", 
						nextpos, tokindex,longest_start_iter, longest_end_iter, search_tokindex, 
						onegap[nextpos].gap);*/
					///get the real aXb for sorting purpose.
					gappyYes = false;
					prePos = 0;
					for(i=0; i < MAX_rule_symbols; i++){
						if(i < longest_start_iter + 1 + longest_end_iter){
							if(!gappyYes){
								if(i == longest_start_iter + 1 - 1){
									onegapPattern[nextpos].pattern[i] = -1;
									gappyYes = true;
									prePos = i;
									//represents the gap
								} else {
									onegapPattern[nextpos].pattern[i] = 
										qrysoffsettok[qryscount+tokindex+i];
									/*if(qrysoffsettok[qryscount+tokindex+i]==-1){
									  printf("Wrong! Kernel one gap enu 1! -1|longest_len_end %d |longest_len_iter %d|\n", 
									  longest_len_end, longest_end_iter);
									  }*/
								}
							} else {
								onegapPattern[nextpos].pattern[i] = 
									qrysoffsettok[qryscount+search_tokindex+i-1-prePos];
								/*if(search_tokindex+i-1-prePos>search_tokindex){
									printf("Wrong!!! inside one gap enu pattern generation\n");
									return;
								}									
								if(qrysoffsettok[qryscount+search_tokindex+i-1-prePos]==-1){
								  printf("Wrong! Kernel one gap enu 2! -1|longest_len_end %d |longest_len_iter %d|\n", 
								  longest_len_end, longest_end_iter);
								  }*/
							}
						} else {
							onegapPattern[nextpos].pattern[i] = -2;
						}
					}
					onegapPattern[nextpos].number = longest_start_iter+1+longest_end_iter;
					longest_end_iter++;				
				}			
			}
			search_tokindex++; 
		}
		longest_start_iter++;
	}			
}

__global__ void zeroOneDiff(
		uint8_t* zeroOneDiffArray_d,
		gapPattern* onegapPattern,
		unsigned int counter) {
	
	int blocks = gridDim.x*blockIdx.y+blockIdx.x;
	unsigned int index = threadIdx.x+blocks*blockDim.x;
	
	bool stop = true;
	if(index < counter){
		//zeroOneDiffArray_d[index] = 0;
		if(index ==0){
			zeroOneDiffArray_d[index] = 1;
		} else {
			if (onegapPattern[index-1].number != onegapPattern[index].number){		
				zeroOneDiffArray_d[index] = 1;
			} else {
				for(int i =0; i < onegapPattern[index].number && stop; i++){
					if(onegapPattern[index-1].pattern[i] 
							!= onegapPattern[index].pattern[i]){
						zeroOneDiffArray_d[index] = 1;
						stop = false;
					}
				}		
			}
		}
	}
}

__global__ void zeroOneDiffTwoGap(
		uint8_t* zeroOneDiffArray_d,
		twoGapPattern* twoGapPattern,
		unsigned int counter) {
			
	int blocks = gridDim.x*blockIdx.y+blockIdx.x;
	unsigned int index = threadIdx.x+blocks*blockDim.x;
	
	//unsigned int index = threadIdx.x+blockIdx.y*blockDim.x;
	if(index >= counter){
		return;
	}

	if(index ==0){
		zeroOneDiffArray_d[index] = 1;
		return;
	}

	if (twoGapPattern[index-1].number != twoGapPattern[index].number){		
		zeroOneDiffArray_d[index] = 1;
		return;
	} else if (twoGapPattern[index-1].blockid != twoGapPattern[index].blockid){
		zeroOneDiffArray_d[index] = 1;
		return;
	}

	for(int i =0; i < twoGapPattern[index].number; i++){
		if(twoGapPattern[index-1].pattern[i] 
				!= twoGapPattern[index].pattern[i]){
			zeroOneDiffArray_d[index] = 1;
			return;
		}
	}

	zeroOneDiffArray_d[index] = 0;		
}


	int compareUserTotal1(const void *v1, const void *v2)
	{
		const precomp_tok* u1 = (precomp_tok*) v1;
		const precomp_tok* u2 = (precomp_tok*) v2;
		return (u1->length < u2->length);
	}


	int compareUserTotal2(const void *v1, const void *v2)
	{
		const precomp_tok* u1 = (precomp_tok*) v1;
		const precomp_tok* u2 = (precomp_tok*) v2;
		return (u1->token  > u2->token);
	}

	int compareUserTotal3(const void *v1, const void *v2)
	{
		const precompute_enu_2* u1 = (precompute_enu_2*) v1;
		const precompute_enu_2* u2 = (precompute_enu_2*) v2;
		return (u1->index  > u2->index) 
			|| ((u1->index == u2->index) && (u1->start > u2->start )) 
			|| ((u1->index == u2->index) && (u1->start == u2->start ) && (u1->length > u2->length));
	}

	void preComputation(
			int* sa_d, 
			int* str_d, 
			ref_t* ref,
			unsigned int* RLP_d,
			uint8_t* L_target_d,
			uint8_t* R_target_d){
		mytimer_t __t;
		timer_start(&__t);

		int i = 0;
		int prevtoken;
		size_t freeMem = 0;
		size_t totalMem = 0;
		size_t allocMem = 0;

		while(ref->str[ref->sa[i]] < 2){
			i++;
			prevtoken = ref->str[ref->sa[i]];
		}
		int prevTokCount = 0;
		int counter_toplist = 0;
		precomp_tok* toplist = (precomp_tok*)malloc(ref->distinctTokenCount*sizeof(precomp_tok));
		i = 0;
		for(; i< ref->toklen; i++){
			if (ref->str[ref->sa[i]] < 2){
				continue;
			}
			if(prevtoken != ref->str[ref->sa[i]]){		
				toplist[counter_toplist].end = i-1;
				toplist[counter_toplist].length = prevTokCount;
				toplist[counter_toplist].token = prevtoken;
				counter_toplist++;
				prevtoken =  ref->str[ref->sa[i]]; 	
				prevTokCount = 1;
			} else {
				prevTokCount++;
			}
		}	
		toplist[counter_toplist].end = i-1;
		toplist[counter_toplist].length = prevTokCount;
		toplist[counter_toplist].token = prevtoken;
		counter_toplist++;
		qsort(toplist, counter_toplist, sizeof(precomp_tok), compareUserTotal1);
		qsort(toplist, PRECOMPUTECOUNT, sizeof(precomp_tok), compareUserTotal2);

		//verification - top 1000 chinese words in Suffix array - looks good.
		/*for(int j=0; j< PRECOMPUTECOUNT; j++){
		  if (toplist[j].token < 2){			
		  fprintf(stderr, "NOT POSSIBLE! LESS Than 2\n");
		//continue;
		}
		struct my_struct *s;
		HASH_FIND_INT(intchar, &toplist[j].token, s );	
		if(!s){
		printf("CANNOT FIND this WORD! - IMPOSSIBLE!!  -->");
		}
		printf("%d -> LEN %d | Tok %d - %s| E %d\n", j, toplist[j].length, toplist[j].token, s->name, toplist[j].end);		
		}*/

		ref->frequentList = (int*)malloc(PRECOMPUTECOUNT*sizeof(int));
		for(int j=0; j< PRECOMPUTECOUNT; j++){
			if (toplist[j].token < 2){			
				fprintf(stderr, "NOT POSSIBLE! LESS Than 2\n");
			}
			ref->frequentList[j] = toplist[j].token;
		}

		precompute_enu* oneGapPrecomp = (precompute_enu*)malloc(PRECOMPUTECOUNT*PRECOMPUTECOUNT*sizeof(precompute_enu));
		for(int cc=0; cc<PRECOMPUTECOUNT ; cc++){		
			for(int jj=0; jj<PRECOMPUTECOUNT; jj++){
				if(toplist[jj].length >= toplist[cc].length){
					oneGapPrecomp[cc*PRECOMPUTECOUNT+jj].length = toplist[cc].length;
					oneGapPrecomp[cc*PRECOMPUTECOUNT+jj].token_a = toplist[cc].token;				
					oneGapPrecomp[cc*PRECOMPUTECOUNT+jj].token_b = toplist[jj].token;
					oneGapPrecomp[cc*PRECOMPUTECOUNT+jj].start = toplist[cc].end - toplist[cc].length + 1;
					oneGapPrecomp[cc*PRECOMPUTECOUNT+jj].reverse = true;
				} else {
					oneGapPrecomp[cc*PRECOMPUTECOUNT+jj].length = toplist[jj].length;
					oneGapPrecomp[cc*PRECOMPUTECOUNT+jj].token_a = toplist[cc].token;				
					oneGapPrecomp[cc*PRECOMPUTECOUNT+jj].token_b = toplist[jj].token;
					oneGapPrecomp[cc*PRECOMPUTECOUNT+jj].start = toplist[jj].end - toplist[jj].length + 1;
					oneGapPrecomp[cc*PRECOMPUTECOUNT+jj].reverse = false;
				}
			}
		}	

		//Verification on OneGapPrecomp
		/*for(int icao = 0; icao<PRECOMPUTECOUNT*PRECOMPUTECOUNT; icao++){
		  int cc = icao / PRECOMPUTECOUNT;
		  int jj = icao % PRECOMPUTECOUNT;

		  printf("%d - ID A %d; ID B %d - Len A %d; Len B %d - End len %d \n",
		  icao,
		  oneGapPrecomp[icao].token_a,
		  oneGapPrecomp[icao].token_b,
		  toplist[cc].length,
		  toplist[jj].length,
		  oneGapPrecomp[icao].length);
		  }*/

		////CUDA Starts
		int* featureMissingCount_d;
    cudaMalloc((void**)&(featureMissingCount_d), 
                  PRECOMPUTECOUNT*PRECOMPUTECOUNT*sizeof(int));
		precompute_enu* oneGapPrecomp_d;
		cudaMalloc((void**)&(oneGapPrecomp_d), PRECOMPUTECOUNT*PRECOMPUTECOUNT*sizeof(precompute_enu));
		cudaMemcpy(oneGapPrecomp_d, oneGapPrecomp, PRECOMPUTECOUNT*PRECOMPUTECOUNT*sizeof(precompute_enu), cudaMemcpyHostToDevice);

		cudaMemGetInfo(&freeMem, &totalMem);  
		fprintf(stderr, "Memory avaliable: Free: %lu, Total: %lu\n",freeMem/(1024*1024), totalMem/(1024*1024));

		unsigned int count = 0, *count_d;
		cudaMalloc((void**)&count_d, sizeof(unsigned int));
		cudaMemcpy(count_d, &count, sizeof(unsigned int), cudaMemcpyHostToDevice);

		precompute_enu_2* oneGapPrecomp_result_d;
		//precompute_enu_2* oneGapPrecomp_result;
		cudaMalloc((void**)&(oneGapPrecomp_result_d), sizeof(precompute_enu_2)*ONEGAP_PRECOMPUT_PREALLOCATION);

		cudaMemGetInfo(&freeMem, &totalMem);  
		fprintf(stderr, "Memory avaliable: ONEGAP - Free: %lu, Total: %lu\n",freeMem/(1024*1024), totalMem/(1024*1024));

#define THREADS_PER_BLOCK 512 
		dim3 block3(THREADS_PER_BLOCK, 1);
		//dim3 grid3(20000, 1);//Tune this and see the difference!
		dim3 grid3(1000, (PRECOMPUTECOUNT*PRECOMPUTECOUNT+1000-1)/1000);

		assert(str_d!=NULL && sa_d != NULL && oneGapPrecomp_result_d != NULL
			&& RLP_d!=NULL && L_target_d!=NULL && R_target_d!=NULL);
		precomp<<<grid3, block3>>>(
			oneGapPrecomp_d, 
			str_d, 
			sa_d, 
			oneGapPrecomp_result_d, 
			count_d,
			RLP_d,
			L_target_d,
			R_target_d,
			featureMissingCount_d);
		cudaThreadSynchronize();

		cudaError_t error = cudaGetLastError();
		if(error != cudaSuccess)
		{
			fprintf(stderr, "CUDA error second After Precomputation Kernel: %s \n", cudaGetErrorString(error));
			exit(-1);
		}	
		cudaMemGetInfo(&freeMem, &totalMem);  
		fprintf(stderr, "Memory avaliable: Free: %lu, Total: %lu\n",freeMem/(1024*1024), totalMem/(1024*1024));

		cudaMemcpy(&count, count_d, sizeof(unsigned int), cudaMemcpyDeviceToHost);
		fprintf(stderr, "Found %u pairs!!\n", count);

		precompute_enu_2* precomp_onegap_large;
		cudaMallocHost((void**)&(precomp_onegap_large), count*sizeof(precompute_enu_2));
		cudaMemcpy(precomp_onegap_large, oneGapPrecomp_result_d, count*sizeof(precompute_enu_2), cudaMemcpyDeviceToHost);

		ref->featureMissingCount = (int*)malloc(sizeof(int)*PRECOMPUTECOUNT*PRECOMPUTECOUNT);
    cudaMemcpy(ref->featureMissingCount, featureMissingCount_d, 
                        sizeof(int)*PRECOMPUTECOUNT*PRECOMPUTECOUNT, cudaMemcpyDeviceToHost);
		///Process results
		cudaFree(featureMissingCount_d);
		cudaFree(oneGapPrecomp_result_d);	
		cudaFree(oneGapPrecomp_d);	
		free(oneGapPrecomp);
		free(toplist);

		qsort(precomp_onegap_large, count, sizeof(precompute_enu_2), compareUserTotal3);
		cudaMallocHost((void **)&(ref->precomp_index), PRECOMPUTECOUNT*PRECOMPUTECOUNT*sizeof(precomp_st_end));
		unsigned int prevC=0;
		bool firsttime = true;
		for(unsigned int ic = 0; ic<PRECOMPUTECOUNT*PRECOMPUTECOUNT; ic++ ){
			firsttime = true;
			ref->precomp_index[ic].start = 1;
			ref->precomp_index[ic].end = 0;
			while(precomp_onegap_large[prevC].index == ic){
				if(firsttime){
					ref->precomp_index[ic].start = prevC;	
					firsttime = false;
				}
				ref->precomp_index[ic].end = prevC;
				prevC++;
			}
		}

		/*for(unsigned int ic = 0; ic<PRECOMPUTECOUNT*PRECOMPUTECOUNT; ic++ ){
		  if(ref->precomp_index[ic].start == -1){
		  continue;
		  }
		  printf("%d -> s %d | e%d\n", ic, ref->precomp_index[ic].start, ref->prec);
		  }*/

		cudaMallocHost((void**)&(ref->precomp_onegap), count*sizeof(precompute_enu_3));
		for(unsigned int yoshi = 0; yoshi < count; yoshi++){
			if(precomp_onegap_large[yoshi].start+precomp_onegap_large[yoshi].length > 
				ref->toklen || precomp_onegap_large[yoshi].length <= 0 
				|| precomp_onegap_large[yoshi].length+1 > MAX_rule_span){
				printf("Not possible inside precomputation!|index %d|end %d\n", 
				yoshi, precomp_onegap_large[yoshi].start+precomp_onegap_large[yoshi].length);
			}
			ref->precomp_onegap[yoshi].start = precomp_onegap_large[yoshi].start;
			ref->precomp_onegap[yoshi].length = precomp_onegap_large[yoshi].length;
		}
		cudaFreeHost(precomp_onegap_large);
		ref->precomp_count = count;
		timer_stop(&__t);	
		fprintf(stderr, "-> Precmputation time: %f\n", timer_elapsed(&__t)/1000);	
	}

void suffixArraySearch(
		void *handler,
		ref_set_t *refset, 
		qry_set_t *qryset, 
		int minmatch,
		timing_t *timing,
		int twoway) {

	mytimer_t __t;
	cudaError_t err;
	context_t * ctx = (context_t *)handler;

	ref_t *ref_d	= &ctx->ref_d;
	ref_t_target* ref_target_d = &ctx->ref_target_d;	

	qry_set_t *qryset_d = &ctx->qryset_d;
	size_t freeMem = 0;
	size_t totalMem = 0;
	size_t allocMem = 0;

	/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
	qryset_d->qrysoffsettok = (int *)qryset_d->qrysbuf; 
	qryset_d->qryscount  = qryset->qryscount;
	qryset_d->totaltokens =qryset->totaltokens; 

	timer_start(&__t);
	/* move the query set to the device */
	cudaMemcpy(qryset_d->qrysoffsettok, qryset->qrysbuf, qryset->qrysbufsize, cudaMemcpyHostToDevice);
	timing->qrysin += timer_stop(&__t);

	std::cout << "Thrust v" << THRUST_MAJOR_VERSION << "." << THRUST_MINOR_VERSION << std::endl;

#define THREADS_PER_BLOCK1 128 
#define Qry_Per_BLOCK 4 

	dim3  block2(THREADS_PER_BLOCK1, Qry_Per_BLOCK);//TBD!!!!
	dim3  grid2(2 ,(qryset->qryscount+block2.y-1)/block2.y);

	int refindex;
	for (refindex = 0;	refindex < refset->count;	refindex++) {
		ref_t * ref_h = &refset->refs[refindex];
		ref_t_target* ref_target_h = &refset->refs_target[refindex];	
		
		//ref_d->str 	  = (int *)ref_d->buf;
		//ref_d->sa      = (int *)ref_d->buf + (ref_h->sa - ref_h->str);
		ref_d->lcp	  = (int *)(ref_d->buf + (ref_h->lcp - ref_h->buf));
		ref_d->rank   = (int *)(ref_d->buf + (ref_h->rank - ref_h->buf));
		ref_d->lcpleft = (int *)(ref_d->buf + (ref_h->lcpleft - ref_h->buf));
		ref_d->lcpright = (int *)(ref_d->buf + (ref_h->lcpright - ref_h->buf));
		ref_d->toklen = ref_h->toklen;

		int datasize = ref_d->toklen * sizeof(int)*2;
		// Copy the sub reference to the device
		timer_start(&__t);
		cudaMemcpy(ref_d->sa, (void *)ref_h->sa, ref_h->toklen*sizeof(int), cudaMemcpyHostToDevice);
		cudaMemcpy(ref_d->str, (void *)ref_h->str, ref_h->toklen*sizeof(int), cudaMemcpyHostToDevice);		
		cudaMemcpy(ref_d->buf, (void *)ref_h->buf, datasize, cudaMemcpyHostToDevice);

		cudaMalloc((void**)&(qryset_d->result_two), qryset->resbufsize);
		timing->refsin += timer_stop(&__t);

		/////////////////////////////////////////////
		///Precompute First! - One time offline!/////
		/////////////////////////////////////////////
		cudaMalloc((void**)&(ref_d->RLP), ref_h->toklen*sizeof(unsigned int));
		cudaMalloc((void**)&(ref_target_d->L_tar), ref_target_h->toklen*sizeof(uint8_t));
		cudaMalloc((void**)&(ref_target_d->R_tar), ref_target_h->toklen*sizeof(uint8_t));
				
		cudaMemcpy(ref_target_d->R_tar, ref_target_h->R_tar, ref_target_h->toklen*sizeof(uint8_t), cudaMemcpyHostToDevice);
		cudaMemcpy(ref_target_d->L_tar, ref_target_h->L_tar, ref_target_h->toklen*sizeof(uint8_t), cudaMemcpyHostToDevice);
		cudaMemcpy(ref_d->RLP, ref_h->RLP, ref_h->toklen*sizeof(unsigned int), cudaMemcpyHostToDevice);
		
		preComputation(ref_d->sa, 
							ref_d->str, 
							ref_h, 
							ref_d->RLP,
							ref_target_d->L_tar,
							ref_target_d->R_tar);
		cudaFree(ref_d->RLP);
		cudaFree(ref_target_d->L_tar);
		cudaFree(ref_target_d->R_tar);
		
		////////////////////////////////////////////
		/////////////////Suffix Array lookUp////////
		////////////////////////////////////////////
		timer_start(&__t);
		
		cudaMemGetInfo(&freeMem, &totalMem);  
		fprintf(stderr, "Memory avaliable pass1: Free: %lu, Total: %lu\n",freeMem/(1024*1024), totalMem/(1024*1024));
		suffixArrayFindLwRwKernelTwoWayTDI<<<grid2, block2>>>(
				ref_d->str,
				ref_d->sa,
				ref_d->lcpleft,
				ref_d->lcpright,
				ref_d->toklen,
				qryset_d->qrysoffsettok,
				qryset->qryscount,
				qryset->totaltokens,
				qryset_d->result_two);

		cudaThreadSynchronize();

		timing->kernel += timer_stop(&__t);
		cudaError_t error = cudaGetLastError();
		fprintf(stderr, "Kernel 1 time: %f\n", timer_elapsed(&__t)/1000);
		if(error != cudaSuccess)
		{
			// print the CUDA error message and exit if any
			printf("CUDA error second After kernel 1: %s \n", cudaGetErrorString(error));
			exit(0);
		}
		// Copy the results from the device
		timer_start(&__t);
		cudaMemcpy(qryset->result_two,qryset_d->result_two,qryset->resbufsize, cudaMemcpyDeviceToHost);
		
	error = cudaGetLastError();
		if(error != cudaSuccess)
		{
			// print the CUDA error message and exit if any
			printf("CUDA error second IN Between Memcpy: %s \n", cudaGetErrorString(error));
			exit(0);
		}
		int ii = 0;
		qryset->totalconnect = 0;
		for(;ii< qryset->totaltokens;ii++){
			//fprintf(stderr, "%d -> qryresults longest match %d\n", ii, qryset->result_two[ii].longestmatch);
			if (qryset->result_two[ii].longestmatch-1 > 0){
				qryset->connectoffset[ii] = qryset->totalconnect;
				qryset->totalconnect += qryset->result_two[ii].longestmatch-1;
			}else {
				qryset->connectoffset[ii] = -1;
			}
		}

		timing->resout += timer_stop(&__t);
		timer_start(&__t);
#define THREADS_PER_BLOCK1 128
#define LONGEST_PER_BLOCK 6 
		dim3  block(THREADS_PER_BLOCK1, 1);
		dim3  grid((qryset->totaltokens+block.x -1)/block.x ,LONGEST_PER_BLOCK);
		//printf("qryset->totalconnect %d\n", qryset->totalconnect);
		cudaMalloc((void**)&(qryset_d->result_connect), qryset->totalconnect*sizeof(result_t));
		cudaMemcpy(qryset_d->connectoffset, qryset->connectoffset, qryset->totaltokens*sizeof(int), cudaMemcpyHostToDevice);
		timing->qrysin += timer_stop(&__t); 

		timer_start(&__t);
		
		cudaMemGetInfo(&freeMem, &totalMem);  
		fprintf(stderr, "Memory avaliable pass2: Free: %lu, Total: %lu\n",freeMem/(1024*1024), totalMem/(1024*1024));
		suffixArrayFindConnectionTwoWayTDI<<<grid, block>>>(
				ref_d->str,
				ref_d->sa,
				/*ref_d->lcp,*/						
				ref_d->lcpleft,
				ref_d->lcpright,
				ref_d->toklen,
				qryset_d->qrysoffsettok,
				qryset->qryscount,
				qryset->totaltokens,
				qryset_d->result_two,
				qryset_d->connectoffset,
				qryset->totalconnect,
				qryset_d->result_connect);

		cudaThreadSynchronize();					
		timing->kernel += timer_stop(&__t);
		timer_start(&__t);
		error = cudaGetLastError();
		if(error != cudaSuccess)
		{
			// print the CUDA error message and exit if any
			printf("CUDA error second Second Kernel: %s \n", cudaGetErrorString(error));
			//exit(-1);
		}

		cudaMallocHost((void **)&(qryset->result_connect), qryset->totalconnect*sizeof(result_t));
		cudaMemcpy(qryset->result_connect,qryset_d->result_connect, qryset->totalconnect*sizeof(result_t), cudaMemcpyDeviceToHost);

		timing->resout += timer_stop(&__t);
		fprintf(stderr, "Result out time: %f second\n", timer_elapsed(&__t)/1000);      
		fprintf(stderr, "Total Continous Phrase Kernel time: %f second\n", timing->kernel/1000);     

		//////////////////////////////////////////////////////////////////
		////////////////////////////////////////Gappy Phrase Lookup///////
		//////////////////////////////////////////////////////////////////
		cudaFree(ref_d->buf);

		//////////////////////////////////////////////////One gap enu!
#define THREADS_PER_BLOCK_GAP 32
		dim3  block_gap_enu(THREADS_PER_BLOCK_GAP, 1);
		dim3  grid_one_gap_enu(1, (qryset->totaltokens+THREADS_PER_BLOCK_GAP-1)/THREADS_PER_BLOCK_GAP);

		cudaMalloc((void**)&(qryset_d->tokindex_qryindex), qryset->totaltokens*sizeof(int));
		cudaMemcpy(qryset_d->tokindex_qryindex, (void *)qryset->tokindex_qryindex, qryset->totaltokens*sizeof(int), cudaMemcpyHostToDevice);

		cudaMemGetInfo(&freeMem, &totalMem);  
		fprintf(stderr, "Memory avaliable: Free: %lu, Total: %lu\n",freeMem/(1024*1024), totalMem/(1024*1024));

		gappy* onegap;
		cudaMalloc((void**)&(onegap), sizeof(gappy)*ONEGAP_ENU_PREALLOCATION);
		cudaMemGetInfo(&freeMem, &totalMem);  
		fprintf(stderr, "Memory avaliable for one gap without pattern: Free: %lu, Total: %lu\n",freeMem/(1024*1024), totalMem/(1024*1024));

		gapPattern* onegapPattern;
		cudaMalloc((void**)&(onegapPattern), sizeof(gapPattern)*ONEGAP_ENU_PREALLOCATION);
		cudaMemGetInfo(&freeMem, &totalMem);  
		fprintf(stderr, "Memory avaliable for one gap with pattern: Free: %lu, Total: %lu\n",freeMem/(1024*1024), totalMem/(1024*1024));

		unsigned int count = 0, *count_d;
		cudaMalloc((void**)&count_d, sizeof(unsigned int));
		cudaMemcpy(count_d, &count, sizeof(unsigned int), cudaMemcpyHostToDevice);

		error = cudaGetLastError();
		if(error != cudaSuccess)
		{
			fprintf(stderr, "CUDA error second Before kernel One Gap First - Phrase Lookup: %s \n", cudaGetErrorString(error));
			exit(-1);
		}	
		
		//Debug
		//fprintf(stderr, "Tokenscount %d | qryresult[tokenscount-1] %d\n", qryset->totaltokens, qryset->result_two[qryset->totaltokens-2].longestmatch);	
		timer_start(&__t);			
		assert(qryset_d->tokindex_qryindex!=NULL 
				&& onegap!=NULL
				&& qryset_d->qrysoffsettok != NULL);

		oneGapEnumeration<<<grid_one_gap_enu, block_gap_enu>>>(
				qryset_d->qrysoffsettok,
				qryset->qryscount,
				qryset->totaltokens,
				qryset_d->result_two,//qryresult
				qryset->totalconnect,
				qryset_d->tokindex_qryindex,
				onegap,
				onegapPattern,
				count_d);

		cudaThreadSynchronize();		
		error = cudaGetLastError();
		if(error != cudaSuccess)
		{
			fprintf(stderr, "CUDA error second After kernel One Gap Enumer: %s \n", cudaGetErrorString(error));
			exit(-1);
		}	
		cudaMemcpy(&count, count_d, sizeof(unsigned int), cudaMemcpyDeviceToHost);
		timer_stop(&__t);
		fprintf(stderr, "Found %u pairs for one gap enumeration!!\n", count);
		fprintf(stderr, "-> One Gap Enumeration time: %f\n", timer_elapsed(&__t)/1000);      

		/////////////Sorts these gap aXb based on gappy patterns
		timer_start(&__t);		
		thrust::device_ptr<gappytyp> onegap_thrust = thrust::device_pointer_cast(onegap);
		thrust::device_ptr<gapPattern1> onegapPattern_thrust = thrust::device_pointer_cast(onegapPattern);
		thrust::device_ptr<gapPattern1> onegapPattern_thrust_end = thrust::device_pointer_cast(onegapPattern+count);

		//sort by key! with thrust
		thrust::sort_by_key(onegapPattern_thrust, onegapPattern_thrust_end, onegap_thrust, oneGapEnumerationCompare());

		cudaMemGetInfo(&freeMem, &totalMem);  		
		fprintf(stderr, "Memory avaliable: Free: %lu, Total: %lu - Duplication Memory - After Sort\n",freeMem/(1024*1024), totalMem/(1024*1024));

		//Debugiing!! - Show unsorted onegap pattern
		/*gapPattern* debugUnsorterPattern;
		cudaMallocHost((void**)&(debugUnsorterPattern), count*sizeof(gapPattern));
		cudaMemcpy(debugUnsorterPattern, onegapPattern, count*sizeof(gapPattern), cudaMemcpyDeviceToHost);
		for(int k=0;k < count; k++ ){
		///Debuging see sorted one gap pattern 
		printf("Sorted one gap| No. %d --> ", k);
		for(int j=0; j< MAX_rule_symbols; j++){
		printf("%d ", debugUnsorterPattern[k].pattern[j]);					
		}
		printf(" |\n");
		///Debugging end				
		}*/
		///
		
		///Get 1-0 bit array to identify the differences for later sequential scan process
#define THREADS_PER_BLOCK_DIFF 256
		dim3 block_gap_diff(THREADS_PER_BLOCK_DIFF, 1);
		dim3 grid_one_gap_diff(10, (count+THREADS_PER_BLOCK_DIFF*10-1)/(10*THREADS_PER_BLOCK_DIFF));

		uint8_t* zeroOneDiffArray_d;
		cudaMalloc((void**)&(zeroOneDiffArray_d), sizeof(uint8_t)*(count+1));
		cudaMemset(zeroOneDiffArray_d, 0, sizeof(uint8_t)*(count+1));

		cudaMemGetInfo(&freeMem, &totalMem);  
		fprintf(stderr, "Memory avaliable after zero one array: Free: %lu, Total: %lu\n",freeMem/(1024*1024), totalMem/(1024*1024));
		
		printf("ZeroOoneDiff -> %d %d\n", grid_one_gap_diff.x, grid_one_gap_diff.y);
		zeroOneDiff<<<grid_one_gap_diff, block_gap_diff>>>(
				zeroOneDiffArray_d,
				onegapPattern,
				count);

		cudaThreadSynchronize();		
		error = cudaGetLastError();
		if(error != cudaSuccess)
		{
			fprintf(stderr, "CUDA error second After kernel - Zero One Diff: %s \n", cudaGetErrorString(error));
			exit(-1);
		}	

		uint8_t* zeroOneDiffArray;
		cudaMallocHost((void **)&(zeroOneDiffArray), count*sizeof(uint8_t));
		cudaMemcpy(zeroOneDiffArray,zeroOneDiffArray_d, count*sizeof(uint8_t), cudaMemcpyDeviceToHost);

		////////////////Now get sequential scan done			
		qryset->onegapcount_enu = count;
		cudaMallocHost((void**)&(qryset->onegap), qryset->onegapcount_enu*sizeof(gappy));
		cudaMemcpy(qryset->onegap, onegap, qryset->onegapcount_enu*sizeof(gappy), cudaMemcpyDeviceToHost);

		gappy_search* oneGapSearch, *oneGapSearch_d;
		cudaMallocHost((void **)&(oneGapSearch), count*sizeof(gappy_search));

		std::vector<std::vector<unsigned int> > oneGapQueryWithID(qryset->qryscount);
		std::map<unsigned int, bool> checkDup;
		std::map<unsigned int, bool>::iterator iter_check;

		/////////Get the sorted stuff out and see gappy phrase		
		cudaMallocHost((void**)&(qryset->onegapPattern), qryset->onegapcount_enu*sizeof(gapPattern));
		cudaMemcpy(qryset->onegapPattern, onegapPattern, qryset->onegapcount_enu*sizeof(gapPattern), cudaMemcpyDeviceToHost);
		timer_stop(&__t);
		fprintf(stderr, "-> One Gap Enumeration Sorting time: %f\n", timer_elapsed(&__t)/1000);      
		
		//scan
		timer_start(&__t);
		unsigned int distinctOneGapCount =0;
		int initialPos = 0;
		for(int i = 0; i< count; i++){
			///Debuging see sorted one gap pattern 
			/*printf("No. %d --> ", i);
			  for(int j=0; j< MAX_rule_symbols; j++){
			  printf("%d ", qryset->onegapPattern[i].pattern[j]);					
			  }
			  printf(" ||| Diff %d\n",zeroOneDiffArray[i]);*/
			///Debugging end

			if(zeroOneDiffArray[i]==1){
				checkDup.clear();
				initialPos = i;
				oneGapSearch[distinctOneGapCount].position = i;
				oneGapSearch[distinctOneGapCount].gap = qryset->onegap[i].gap;
				oneGapSearch[distinctOneGapCount].qryend_len = qryset->onegap[i].qryend_len;
				oneGapSearch[distinctOneGapCount].qrystart = qryset->onegap[i].qrystart;
				oneGapSearch[distinctOneGapCount].qrystart_len = qryset->onegap[i].qrystart_len;
				oneGapSearch[distinctOneGapCount].end_on_salist = -1;
				oneGapSearch[distinctOneGapCount].start_on_salist = -1;
				
				if(qryset->onegap[i].gap<=0 || qryset->onegap[i].qrystart <0){
					printf("No. %d --> qrystart %d - searchInd %d - gap %d - qrystartlen %d\n", distinctOneGapCount,
						oneGapSearch[distinctOneGapCount].qrystart, 
						oneGapSearch[distinctOneGapCount].qrystart+oneGapSearch[distinctOneGapCount].gap + oneGapSearch[distinctOneGapCount].qrystart_len,
						oneGapSearch[distinctOneGapCount].gap,
						oneGapSearch[distinctOneGapCount].qrystart_len);
				}
				distinctOneGapCount ++;
				///Debuging see sorted one gap pattern 
				/* printf("aXb -> REAL SORTED ID %d: ", distinctOneGapCount-1);
				  for(int j=0; j< MAX_rule_symbols; j++){
				  printf("%d ", qryset->onegapPattern[i].pattern[j]);					
				  }
				  printf("||| qrystart_len %d | qryend_len %d\n", qryset->onegap[i].qrystart_len, qryset->onegap[i].qryend_len);	
				*/  
				///Debugging end
			} 
			int qid = qryset->tokindex_qryindex[qryset->onegap[i].qrystart];
			if(qid > qryset->qryscount){
				printf("Not possible, Qry count exceeded!\n");
				return;
			}
			///remove duplications
			iter_check = checkDup.find(qid);
			if(iter_check == checkDup.end()){
				checkDup.insert(std::make_pair(qid, false));
				oneGapQueryWithID[qid].push_back(distinctOneGapCount-1); //all to the 
				// uniform id, innitial ID/pos
			} 
		}		
		fprintf(stderr, "Distinct one gap pattern: %d\n", distinctOneGapCount);      
		timer_stop(&__t);
                fprintf(stderr, "-> One Gap Enumeration CPU Processing time: %f\n", timer_elapsed(&__t)/1000);

		cudaFree(onegap);
		cudaFree(onegapPattern);
		cudaFree(zeroOneDiffArray_d);
		cudaFreeHost(zeroOneDiffArray);
		cudaMemGetInfo(&freeMem, &totalMem);  
		fprintf(stderr, "Memory avaliable after one gap enumeration finish: Free: %lu, Total: %lu - After One Gap\n",freeMem/(1024*1024), totalMem/(1024*1024));

		////////////Now Look For distinct one gappy phrases in Suffix Array
		count = 0;
		cudaMemcpy(count_d, &count, sizeof(unsigned int), cudaMemcpyHostToDevice);

		cudaMalloc((void**)&(oneGapSearch_d), sizeof(gappy_search)*distinctOneGapCount);
		cudaMemcpy(oneGapSearch_d, oneGapSearch, distinctOneGapCount*sizeof(gappy_search), 
				cudaMemcpyHostToDevice);

		oneGapOnSA* oneGapSA, *oneGapSA_d;		
		cudaMalloc((void**)&(oneGapSA_d), sizeof(oneGapOnSA)*(ONEGAP_PREALLOCATION));

		cudaMalloc((void**)&(ref_d->precomp_index), PRECOMPUTECOUNT*PRECOMPUTECOUNT*sizeof(precomp_st_end));
		cudaMemcpy(ref_d->precomp_index, ref_h->precomp_index, PRECOMPUTECOUNT*PRECOMPUTECOUNT*sizeof(precomp_st_end), cudaMemcpyHostToDevice);

		cudaMalloc((void**)&(ref_d->frequentList), PRECOMPUTECOUNT*sizeof(int));
		cudaMemcpy(ref_d->frequentList, ref_h->frequentList, PRECOMPUTECOUNT*sizeof(int), cudaMemcpyHostToDevice);

		cudaMalloc((void**)&(ref_d->precomp_onegap), ref_h->precomp_count*sizeof(precompute_enu_3));
		cudaMemcpy(ref_d->precomp_onegap, ref_h->precomp_onegap, ref_h->precomp_count*sizeof(precompute_enu_3), cudaMemcpyHostToDevice);

		/////RLP Array memory allocation
		cudaMalloc((void**)&(ref_d->RLP), ref_h->toklen*sizeof(unsigned int));
		cudaMalloc((void**)&(ref_target_d->L_tar), ref_target_h->toklen*sizeof(uint8_t));
		cudaMalloc((void**)&(ref_target_d->R_tar), ref_target_h->toklen*sizeof(uint8_t));
					
		cudaMemcpy(ref_target_d->R_tar, ref_target_h->R_tar, ref_target_h->toklen*sizeof(uint8_t), cudaMemcpyHostToDevice);
		cudaMemcpy(ref_target_d->L_tar, ref_target_h->L_tar, ref_target_h->toklen*sizeof(uint8_t), cudaMemcpyHostToDevice);
		cudaMemcpy(ref_d->RLP, ref_h->RLP, ref_h->toklen*sizeof(unsigned int), cudaMemcpyHostToDevice);
			
		cudaMemGetInfo(&freeMem, &totalMem);  
		fprintf(stderr, "Memory avaliable right before one gap lookup kernel: Free: %lu, Total: %lu - After One Gap\n",freeMem/(1024*1024), totalMem/(1024*1024));

		assert(ref_d->sa 
				&& ref_d->str
				&& ref_d->toklen
				&& qryset_d->qrysoffsettok
				&& qryset_d->result_two
				&& qryset_d->connectoffset
				&& qryset_d->result_connect
				&& qryset_d->tokindex_qryindex
				&& oneGapSA_d
				&& ref_d->frequentList
				&& ref_d->precomp_index
				&& ref_d->precomp_onegap
				&& oneGapSearch_d);

#define THREADS_PER_BLOCK_SA 128
		dim3 block_gap_onegap_sa(THREADS_PER_BLOCK_SA, 1);
		dim3 grid_one_gap_sa(500, (distinctOneGapCount+500-1)/500);

		timer_start(&__t);	
		oneGapLookUpSA<<<grid_one_gap_sa, block_gap_onegap_sa>>>(
				ref_d->str,
				ref_d->sa,				
				ref_d->toklen,
				qryset_d->qrysoffsettok,
				qryset->qryscount,
				qryset->totaltokens,
				qryset_d->result_two,//qryresult
				qryset_d->connectoffset,
				qryset->totalconnect,
				qryset_d->result_connect,
				qryset_d->tokindex_qryindex,
				oneGapSA_d,
				ref_h->precomp_count,
				ref_d->frequentList,
				ref_d->precomp_index,
				ref_d->precomp_onegap,
				count_d,
				oneGapSearch_d,
				distinctOneGapCount,
				ref_d->RLP,
				ref_target_d->L_tar,
				ref_target_d->R_tar);

		cudaThreadSynchronize();		
		error = cudaGetLastError();
		if(error != cudaSuccess)
		{
			fprintf(stderr, "CUDA error second After kernel - oneGapLookUpSA: %s \n", cudaGetErrorString(error));
			exit(-1);
		}	

		unsigned int countOneGapSA = 0;
		cudaMemcpy(&countOneGapSA, count_d, sizeof(unsigned int), cudaMemcpyDeviceToHost);
		
		fprintf(stderr, "Found %u for one gap on SA!\n", countOneGapSA);

		timer_stop(&__t);
                fprintf(stderr, "-> One Gap Look up on SA Kernel time: %f\n", timer_elapsed(&__t)/1000);
		
		cudaFree(ref_d->sa);
		cudaFree(oneGapSearch_d);
		cudaFree(ref_d->RLP);
		cudaFree(ref_target_d->L_tar);
		cudaFree(ref_target_d->R_tar);
		
		cudaMemGetInfo(&freeMem, &totalMem);  
		fprintf(stderr, "Memory avaliable before thrust one gap sort: Free: %lu, Total: %lu - After One Gap\n",freeMem/(1024*1024), totalMem/(1024*1024));

		timer_start(&__t);	
		thrust::device_ptr<gappytyp_sa> onegapSA_thrust = thrust::device_pointer_cast(oneGapSA_d);
		thrust::device_ptr<gappytyp_sa> onegapSA_thrust_end = thrust::device_pointer_cast(oneGapSA_d+countOneGapSA);
		assert(onegapSA_thrust+countOneGapSA == onegapSA_thrust_end);	
		//sort with thrust based on position attribute - ID
		thrust::sort(onegapSA_thrust, onegapSA_thrust_end, oneGapSACompare());

	  	timer_stop(&__t);
                fprintf(stderr, "-> One Gap on SA Sorting Thrust time: %f\n", timer_elapsed(&__t)/1000);
	
		cudaMallocHost((void**)&(oneGapSA), countOneGapSA*sizeof(oneGapOnSA));
		cudaMemcpy(oneGapSA, oneGapSA_d, countOneGapSA*sizeof(oneGapOnSA), cudaMemcpyDeviceToHost);

		///Debugging
		/*for(int cc=0;cc<countOneGapSA;cc++){
		  printf("%d->%u %u %d\n", cc, oneGapSA[cc].position, oneGapSA[cc].str_position, oneGapSA[cc].length); //on str array);
		  }			*/
		////Debugging end

		cudaFree(oneGapSA_d);			

	 	timer_start(&__t);	
		//Get the start_on_salist/end_on_salist done
		for(int i =0; i<countOneGapSA; i++ ){
			if (oneGapSA[i].position >= distinctOneGapCount){
				printf("Counts not POSSIBLE! Wrong in sequential - postion %d, distinctC %d\n", oneGapSA[i].position, distinctOneGapCount);
				return;
			}
			/*if (oneGapSA[i].position == 548){
				printf("oneGapSA[i].str %d length %d|pattern: ", oneGapSA[i].str_position, oneGapSA[i].length);
				for(int iccc = oneGapSA[i].str_position; iccc <= oneGapSA[i].str_position+oneGapSA[i].length; iccc++){
					printf("%d ", ref_h->str[iccc]);
				}
				printf("\n");
			}*/
			if (oneGapSearch[oneGapSA[i].position].start_on_salist == -1){
				oneGapSearch[oneGapSA[i].position].start_on_salist = i;
			} 

			if (oneGapSearch[oneGapSA[i].position].end_on_salist == -1){
				oneGapSearch[oneGapSA[i].position].end_on_salist = i;
			} else if (oneGapSearch[oneGapSA[i].position].end_on_salist < i){
				oneGapSearch[oneGapSA[i].position].end_on_salist = i;
			}				
		}
		timer_stop(&__t);
                fprintf(stderr, "-> One Gap on SA Processing time: %f\n", timer_elapsed(&__t)/1000);

		///Debugging
		
		/*for(int i =0; i<distinctOneGapCount;i++){
		   printf("ONEGAP|qrystart %d, gap %d||one gap ID %d||First POS %d||start %d, end %d||startlen %d, endlen %d\n", 
		   oneGapSearch[i].qrystart, oneGapSearch[i].gap, i, oneGapSearch[i].position,
		   oneGapSearch[i].start_on_salist, oneGapSearch[i].end_on_salist,
		   oneGapSearch[i].qrystart_len, oneGapSearch[i].qryend_len);
		}*/
		//Book keeping
		cudaMallocHost((void**)&(qryset->oneGapSearch), sizeof(gappy_search)*distinctOneGapCount);
		memcpy(qryset->oneGapSearch, oneGapSearch, sizeof(gappy_search)*distinctOneGapCount);

		//Book keeping
		cudaMallocHost((void**)&(qryset->oneGapSA), sizeof(oneGapOnSA)*countOneGapSA);
		memcpy(qryset->oneGapSA, oneGapSA, countOneGapSA*sizeof(oneGapOnSA));

		qryset->distinctOneGapCount = distinctOneGapCount;
		qryset->countOneGapSA = countOneGapSA;
		qryset->oneGapQueryWithID = oneGapQueryWithID;
		fprintf(stderr, "After one gap look up: distinctOneGapCount %d; countOneGapSA %d\n", distinctOneGapCount, countOneGapSA);
		/////////////////////////////////
		/////////////////////////////////
		/////////////////////////////////
		/////////////////////////////////
		////Two and more gaps!!//////////	

		cudaMalloc((void**)&(oneGapSearch_d), sizeof(gappy_search)*distinctOneGapCount);
		cudaMemcpy(oneGapSearch_d, oneGapSearch, distinctOneGapCount*sizeof(gappy_search), 
				cudaMemcpyHostToDevice);

#define THREADS_PER_BLOCK_GAP_TWO 64
		dim3  block_two_gap_enu(THREADS_PER_BLOCK_GAP_TWO, 1);
		dim3  grid_two_gap_enu(LINEARBLOCK, (distinctOneGapCount+LINEARBLOCK-1)/LINEARBLOCK);

		cudaMalloc((void**)&(qryset_d->onegap), sizeof(gappy)*qryset->onegapcount_enu);
		cudaMemcpy(qryset_d->onegap, qryset->onegap, qryset->onegapcount_enu*sizeof(gappy), 
				cudaMemcpyHostToDevice);

		twogappy* twoGap_d;
		cudaMalloc((void**)&(twoGap_d), sizeof(twogappy)*TWOGAP_ENU_PREALLOCATION);
		cudaMemGetInfo(&freeMem, &totalMem);  
		fprintf(stderr, "Memory avaliable for two gap without pattern: Free: %lu, Total: %lu\n",freeMem/(1024*1024), totalMem/(1024*1024));

		twoGapPattern* twogapPattern_d;
		cudaMalloc((void**)&(twogapPattern_d), sizeof(twoGapPattern)*TWOGAP_ENU_PREALLOCATION);
		cudaMemGetInfo(&freeMem, &totalMem);  
		fprintf(stderr, "Memory avaliable for two gap with pattern: Free: %lu, Total: %lu\n",freeMem/(1024*1024), totalMem/(1024*1024));

		count = 0;
		cudaMemcpy(count_d, &count, sizeof(unsigned int), cudaMemcpyHostToDevice);

		assert(qryset_d->tokindex_qryindex!=NULL 
				&& oneGapSearch_d!=NULL
				&& qryset_d->qrysoffsettok != NULL
				&& twoGap_d != NULL
				&& twogapPattern_d != NULL
				&& qryset_d->qrysoffsettok != NULL
				&& qryset_d->result_two != NULL
				&& qryset_d->onegap != NULL);
		
		timer_start(&__t);
		twoGapEnumeration<<<grid_two_gap_enu, block_two_gap_enu>>>(
				qryset_d->qrysoffsettok,
				qryset->qryscount,
				qryset->totaltokens,
				qryset_d->result_two,//qryresult
				qryset->totalconnect,
				qryset_d->tokindex_qryindex,
				oneGapSearch_d,
				twoGap_d,
				twogapPattern_d,
				count_d,
				distinctOneGapCount,
				qryset_d->onegap,
				qryset->onegapcount_enu);

		cudaThreadSynchronize();		
		error = cudaGetLastError();
		if(error != cudaSuccess)
		{
			fprintf(stderr, "CUDA error second After kernel Two Gap Enumer: %s \n", cudaGetErrorString(error));
			exit(-1);
		}	

		unsigned int countTwoGapEnu = 0;
		cudaMemcpy(&countTwoGapEnu, count_d, sizeof(unsigned int), cudaMemcpyDeviceToHost);

		timer_stop(&__t);
		fprintf(stderr, "Found %u pairs for two gap enumeration!!\n", countTwoGapEnu);
		fprintf(stderr, "-> Two Gap Enumeration time: %f\n", timer_elapsed(&__t)/1000);      

		//Debugiing!! - Show unsorted onegap pattern
		/*gapPattern2* debugUnsorterPattern2;
		cudaMallocHost((void**)&(debugUnsorterPattern2), countTwoGapEnu*sizeof(gapPattern2));
		cudaMemcpy(debugUnsorterPattern2, twogapPattern_d, countTwoGapEnu*sizeof(gapPattern2), cudaMemcpyDeviceToHost);
		for(int k=0;k < countTwoGapEnu; k++ ){
		///Debuging see sorted one gap pattern 
		printf("No. %d --> ", k);
		for(int j=0; j< MAX_rule_symbols-4; j++){
		printf("%d ", debugUnsorterPattern2[k].pattern[j]);					
		}
		printf(" |\n");
		///Debugging end				
		}*/
		///
		timer_start(&__t);
		thrust::device_ptr<gappytyp2> twogap_thrust = thrust::device_pointer_cast(twoGap_d);
		thrust::device_ptr<gapPattern2> twogapPattern_thrust = thrust::device_pointer_cast(twogapPattern_d);

		//sort by key! with thrust
		thrust::sort_by_key(twogapPattern_thrust, twogapPattern_thrust+countTwoGapEnu, twogap_thrust, twoGapEnumerationCompare());

		cudaMemGetInfo(&freeMem, &totalMem);  		
		fprintf(stderr, "Memory avaliable: Free: %lu, Total: %lu - Two Gap Enu Duplication Memory - After Sort\n",freeMem/(1024*1024), totalMem/(1024*1024));

		//Debugiing!! - Show unsorted onegap pattern		
		/*cudaMallocHost((void**)&(debugUnsorterPattern2), countTwoGapEnu*sizeof(gapPattern2));
		cudaMemcpy(debugUnsorterPattern2, twogapPattern_d, countTwoGapEnu*sizeof(gapPattern2), cudaMemcpyDeviceToHost);
		for(int k=0;k < countTwoGapEnu; k++ ){
		///Debuging see sorted one gap pattern 
		printf("Sorted No. %d --> ", k);
		for(int j=0; j< MAX_rule_symbols-4; j++){
		printf("ID %d | %d ", debugUnsorterPattern2[k].blockid, debugUnsorterPattern2[k].pattern[j]);					
		}
		printf(" |\n");
		///Debugging end				
		}*/
		///
		
		///Get 1-0 bit array to identify the differences for later sequential scan process
		dim3 grid_two_gap_diff(10, (countTwoGapEnu+10*THREADS_PER_BLOCK_DIFF-1)/(10*THREADS_PER_BLOCK_DIFF));

		uint8_t* zeroOneDiffArrayTwoGap_d;
		cudaMalloc((void**)&(zeroOneDiffArrayTwoGap_d), sizeof(uint8_t)*(countTwoGapEnu+1));

		cudaMemGetInfo(&freeMem, &totalMem);  
		fprintf(stderr, "Memory avaliable after zero one array two gap: Free: %lu, Total: %lu\n",freeMem/(1024*1024), totalMem/(1024*1024));

		zeroOneDiffTwoGap<<<grid_two_gap_diff, block_gap_diff>>>(
				zeroOneDiffArrayTwoGap_d,
				twogapPattern_d,
				countTwoGapEnu);
		timer_stop(&__t);
		fprintf(stderr, "-> Two Gap Enumeration sort time: %f\n", timer_elapsed(&__t)/1000);

		cudaThreadSynchronize();		
		error = cudaGetLastError();
		if(error != cudaSuccess)
		{
			fprintf(stderr, "CUDA error second After kernel - Zero One Diff - Two gap: %s \n", cudaGetErrorString(error));
			exit(-1);
		}	

		uint8_t* zeroOneDiffArrayTwoGap;
		cudaMallocHost((void **)&(zeroOneDiffArrayTwoGap), countTwoGapEnu*sizeof(uint8_t));
		cudaMemcpy(zeroOneDiffArrayTwoGap, zeroOneDiffArrayTwoGap_d, countTwoGapEnu*sizeof(uint8_t), cudaMemcpyDeviceToHost);

		//////////////////////////////////////
		///Sequential scan! Two Gaps
		////////////////Now get sequential scan done

		//Book keeping
		qryset->twogapcount_enu = countTwoGapEnu;
		cudaMallocHost((void**)&(qryset->twoGap), 
				qryset->twogapcount_enu*sizeof(twogappy));
		cudaMemcpy(qryset->twoGap, twoGap_d, 
				qryset->twogapcount_enu*sizeof(twogappy), cudaMemcpyDeviceToHost);

		two_gappy_search* twoGapSearch, *twoGapSearch_d;
		cudaMallocHost((void **)&(twoGapSearch), countTwoGapEnu*sizeof(two_gappy_search));// hard coded - can be tune

		//Book keeping
		cudaMallocHost((void**)&(qryset->twogapPattern), 
				qryset->twogapcount_enu*sizeof(twoGapPattern));
		cudaMemcpy(qryset->twogapPattern, twogapPattern_d, 
				qryset->twogapcount_enu*sizeof(twoGapPattern), cudaMemcpyDeviceToHost);

		timer_start(&__t);
		std::vector<std::vector<unsigned int> > twoGapQueryWithID(qryset->qryscount);
		checkDup.clear();
		//scan
		int distinctTwoGapCount =0;
		initialPos = 0;
		for(int i = 0; i< countTwoGapEnu; i++){
			///Debuging see sorted one gap pattern 
			/*	printf("No. %d --> oneGapId %d %d| number %d ", i, qryset->twogapPattern[i].blockid, 
				qryset->twoGap[i].blockid, qryset->twogapPattern[i].number);
				for(int j=0; j< MAX_rule_symbols-4; j++){
				printf("%d ", qryset->twogapPattern[i].pattern[j]);					
				}
				printf("||| qryend_len %d ||| distinct %d\n", qryset->twoGap[i].qryend_len, distinctTwoGapCount);	*/
			///Debugging end

			if(zeroOneDiffArrayTwoGap[i]==1){
				checkDup.clear();
				initialPos = i;
				twoGapSearch[distinctTwoGapCount].blockid = 
					qryset->twoGap[i].blockid;
				twoGapSearch[distinctTwoGapCount].position = i;
				twoGapSearch[distinctTwoGapCount].qryend_len = qryset->twoGap[i].qryend_len;
				twoGapSearch[distinctTwoGapCount].gap2 = 
					qryset->twoGap[i].gap2;
				twoGapSearch[distinctTwoGapCount].end_on_salist = -1;
				twoGapSearch[distinctTwoGapCount].start_on_salist = -1;
				distinctTwoGapCount ++;					
			} 
			int qid = qryset->tokindex_qryindex[qryset->twoGap[i].gap2];
			if(qid > qryset->qryscount){
				printf("Not possible, Qry count exceeded! in two Gap\n");
				return;
			}
			///remove duplications
			iter_check = checkDup.find(qid);
			if(iter_check == checkDup.end()){
				checkDup.insert(std::make_pair(qid, false));
				twoGapQueryWithID[qid].push_back(distinctTwoGapCount-1); //all to the 
				// uniform id, innitial ID/pos on the two gap enu array
			} 
		}

		/////////Get the sorted stuff out and see gappy phrase		

		timer_stop(&__t);
		fprintf(stderr, "-> Two Gap Enumeration CPU Processing time: %f\n", timer_elapsed(&__t)/1000);      
		fprintf(stderr, "Disticnt two gap pattern: %d\n", distinctTwoGapCount);

		cudaFree(twogapPattern_d);
		cudaFree(zeroOneDiffArrayTwoGap_d);
		cudaFree(twoGap_d);
		cudaFreeHost(zeroOneDiffArrayTwoGap);

		cudaMemGetInfo(&freeMem, &totalMem);  
		fprintf(stderr, "Memory avaliable after two gap enumeration finish: Free: %lu, Total: %lu - After One Gap\n",freeMem/(1024*1024), totalMem/(1024*1024));

		//////////Now do the real Suffix Array matching!			
		//////////Now Look For distinct one gappy phrases in Suffix Array
		count = 0;
		cudaMemcpy(count_d, &count, sizeof(unsigned int), cudaMemcpyHostToDevice);

		cudaMalloc((void**)&(twoGapSearch_d), 
				sizeof(two_gappy_search)*distinctTwoGapCount);
		cudaMemcpy(twoGapSearch_d, twoGapSearch, 
				sizeof(two_gappy_search)*distinctTwoGapCount, 
				cudaMemcpyHostToDevice);

		cudaMalloc((void**)&(oneGapSA_d), sizeof(oneGapOnSA)*countOneGapSA);
		cudaMemcpy(oneGapSA_d, oneGapSA, sizeof(oneGapOnSA)*countOneGapSA,
				cudaMemcpyHostToDevice);

		twoGapOnSA* twoGapSA, *twoGapSA_d;		
		cudaMalloc((void**)&(twoGapSA_d), sizeof(twoGapOnSA)*(ONEGAP_PRECOMPUT_PREALLOCATION));

		/////RLP Array memory allocation
		cudaMalloc((void**)&(ref_d->RLP), ref_h->toklen*sizeof(unsigned int));
		cudaMalloc((void**)&(ref_target_d->L_tar), ref_target_h->toklen*sizeof(uint8_t));
		cudaMalloc((void**)&(ref_target_d->R_tar), ref_target_h->toklen*sizeof(uint8_t));
							
		cudaMemcpy(ref_target_d->R_tar, ref_target_h->R_tar, ref_target_h->toklen*sizeof(uint8_t), cudaMemcpyHostToDevice);
		cudaMemcpy(ref_target_d->L_tar, ref_target_h->L_tar, ref_target_h->toklen*sizeof(uint8_t), cudaMemcpyHostToDevice);
		cudaMemcpy(ref_d->RLP, ref_h->RLP, ref_h->toklen*sizeof(unsigned int), cudaMemcpyHostToDevice);
					
		cudaMemGetInfo(&freeMem, &totalMem);  
		fprintf(stderr, "Memory avaliable right before two gap look up: Free: %lu, Total: %lu - After One Gap\n",freeMem/(1024*1024), totalMem/(1024*1024));

#define THREADS_PER_BLOCK_SA_TWO 128
		dim3 block_gap_twogap_sa(THREADS_PER_BLOCK_SA_TWO, 1);
		dim3 grid_two_gap_sa(500, (distinctTwoGapCount+500-1)/500);

		timer_start(&__t);
		twoGapLookUpSA<<<grid_two_gap_sa, block_gap_twogap_sa>>>(
				ref_d->str,
				ref_d->toklen,
				qryset_d->qrysoffsettok,
				qryset->qryscount,
				qryset->totaltokens,
				qryset_d->result_two,//qryresult
				qryset_d->connectoffset,
				qryset->totalconnect,
				qryset_d->result_connect,
				qryset_d->tokindex_qryindex,
				twoGapSA_d,
				ref_h->precomp_count,
				ref_d->frequentList,
				ref_d->precomp_index,
				ref_d->precomp_onegap,
				count_d,
				twoGapSearch_d,
				distinctTwoGapCount,
				distinctOneGapCount,
				oneGapSearch_d,
				oneGapSA_d,
				countOneGapSA,
				ref_d->RLP,
				ref_target_d->L_tar,
				ref_target_d->R_tar);

		cudaThreadSynchronize();		
		error = cudaGetLastError();
		if(error != cudaSuccess)
		{
			fprintf(stderr, "CUDA error second After kernel twoGapLookUpSA: %s \n", cudaGetErrorString(error));
			exit(-1);
		}	

		unsigned int countTwoGapSA = 0;
		cudaMemcpy(&countTwoGapSA, count_d, sizeof(unsigned int), cudaMemcpyDeviceToHost);
		fprintf(stderr, "Found %d two gap on SA\n", countTwoGapSA);

		timer_stop(&__t);	
		fprintf(stderr, "-> Two Gap Look Up on SA Kernel: %f\n", timer_elapsed(&__t)/1000);

		cudaMallocHost((void**)&(twoGapSA), countTwoGapSA*sizeof(twoGapOnSA));

		cudaFree(twoGapSearch_d);
		cudaFree(oneGapSA_d);
		cudaFree(oneGapSearch_d);
		cudaFree(ref_d->frequentList);
		cudaFree(ref_d->precomp_index);	
		cudaFree(ref_d->precomp_onegap);
		cudaFree(ref_d->RLP);
		cudaFree(ref_target_d->L_tar);
		cudaFree(ref_target_d->R_tar);
		
		timer_start(&__t);	
		thrust::device_ptr<gappytyp_sa2> twogapSA_thrust = thrust::device_pointer_cast(twoGapSA_d);
		//sort with thrust based on position attribute - ID
		thrust::sort(twogapSA_thrust, twogapSA_thrust+countTwoGapSA, twoGapSACompare());
		
		timer_stop(&__t);
		fprintf(stderr, "-> Two Gap on SA Sorting Kernel: %f\n", timer_elapsed(&__t)/1000);

		cudaMemcpy(twoGapSA, twoGapSA_d, countTwoGapSA*sizeof(twoGapOnSA), cudaMemcpyDeviceToHost);
		
		timer_start(&__t);
		//Get the start_on_salist/end_on_salist done
		for(int i =0; i<countTwoGapSA; i++ ){
			if (twoGapSA[i].position >= distinctTwoGapCount){
				printf("Counts not POSSIBLE! Wrong sequential scan\n");
				return;
			}
			///Debuging see sorted one gap pattern 
			/*printf("No. %d --> position %d|str_position %d\n", i, twoGapSA[i].position,
					twoGapSA[i].str_position);*/
			///Debugging end
			
			if (twoGapSearch[twoGapSA[i].position].start_on_salist == -1){
				twoGapSearch[twoGapSA[i].position].start_on_salist = i;
			} 

			if (twoGapSearch[twoGapSA[i].position].end_on_salist == -1){
				twoGapSearch[twoGapSA[i].position].end_on_salist = i;
			} else if (twoGapSearch[twoGapSA[i].position].end_on_salist < i){
				twoGapSearch[twoGapSA[i].position].end_on_salist = i;
			}				
		}		
		
		//debugging
		/*for(int i =0; i< distinctTwoGapCount;i++){
			printf("No. %d --> OneGapId %d | SearchTok2 %d | start %d|end %d\n", i, twoGapSearch[i].blockid, twoGapSearch[i].gap2, twoGapSearch[i].start_on_salist,
					twoGapSearch[i].end_on_salist);
			///Debugging end
		}*/
			
		timer_stop(&__t);	
		fprintf(stderr, "-> Two Gap on SA CPU processing time: %f\n", timer_elapsed(&__t)/1000);

		//Book keeping
		cudaMallocHost((void**)&(qryset->twoGapSearch), 
				sizeof(two_gappy_search)*distinctTwoGapCount);
		memcpy(qryset->twoGapSearch, twoGapSearch, sizeof(two_gappy_search)*distinctTwoGapCount);
		qryset->distinctTwoGapCount = distinctTwoGapCount;

		//Book keeping
		cudaMallocHost((void**)&(qryset->twoGapSA), 
				sizeof(twoGapOnSA)*countTwoGapSA);
		memcpy(qryset->twoGapSA, twoGapSA, countTwoGapSA*sizeof(twoGapOnSA));
		qryset->countTwoGapSA = countTwoGapSA;
		qryset->twoGapQueryWithID = twoGapQueryWithID;

		fprintf(stderr, "Suffix Array Kernel Clean Up\n");
		cudaFree(twoGapSA_d);			
		cudaFreeHost(oneGapSearch);
		cudaFreeHost(oneGapSA);
		cudaFreeHost(twoGapSA);
		cudaFreeHost(twoGapSearch);

		cudaMemGetInfo(&freeMem, &totalMem);  
		fprintf(stderr, "Memory avaliable: Free: %lu, Total: %lu - End of Suffix.cu\n",freeMem/(1024*1024), totalMem/(1024*1024));
	}

}



////////////////////////////////////////////////	
////Deduplication on SEED one gap precomutation
////////////////////////////////////////////////
/*printf("Start Deduplication - %d\n", qryset->onegapcount);
  timer_start(&__t);
  thrust::device_ptr<gappytyp> dev_ptr = thrust::device_pointer_cast(qryset_d->onegap);

  cudaMemGetInfo(&freeMem, &totalMem);  
  fprintf(stderr, "Memory avaliable: Free: %lu, Total: %lu - Duplication Memory - Before Sort\n",freeMem/(1024*1024), totalMem/(1024*1024));

  thrust::sort(dev_ptr, dev_ptr+qryset->onegapcount, oneGapCompare());
  cudaMemGetInfo(&freeMem, &totalMem);  

  fprintf(stderr, "Memory avaliable: Free: %lu, Total: %lu - Duplication Memory - After Sort\n",freeMem/(1024*1024), totalMem/(1024*1024));

  cudaMemcpy(qryset->onegap, qryset_d->onegap, sizeof(gappy)*qryset->onegapcount, cudaMemcpyDeviceToHost);
  cudaFree(qryset_d->onegap);*/

//qsort(qryset->onegap, qryset->onegapcount, sizeof(gappy), compareUserTotal_gappy);
//ref_t* ref_hh = &refset->refs[0];

/*for(int i = 0; i< qryset->onegapcount; i++){
  if(qryset->onegap[i].gap+1 > MAX_rule_span){
  printf("WRRRONG! 1 %d\n", i);
  }

  if(qryset->onegap[i].qryend - qryset->onegap[i].qrystart <= 1){
  printf("WRRRONG! 2 %d\n", i);
  }

  if(qryset->onegap[i].qrystart_len + qryset->onegap[i].qryend_len + 1 > MAX_rule_symbols){
  printf("WRRRONG! 3 %d\n", i);
  }	

  if(qryset->tokindex_qryindex[qryset->onegap[i].qrystart] 
  != qryset->tokindex_qryindex[qryset->onegap[i].qryend+qryset->onegap[i].qryend_len-1]){
  printf("WRRRONG! 4 %d\n", i);
  }	

  if(ref_hh->str[qryset->onegap[i].refstr_s] 
  != qryset->qrysoffsettok[qryset->qryscount+qryset->onegap[i].qrystart]
  || ref_hh->str [qryset->onegap[i].refstr_s + qryset->onegap[i].gap] 
  != qryset->qrysoffsettok[qryset->qryscount+qryset->onegap[i].qryend + qryset->onegap[i].qryend_len-1]){
  printf("WRRRONG! 5 %d\n", i);
  }				
  }

  timer_stop(&__t);	
  fprintf(stderr, "-> Deduplication qsort time: %f\n", 
  timer_elapsed(&__t)/1000);*/
////////////////////////////////////////Normal Operations


//thrust::host_vector<gappytyp> h_structures(qryset->onegap, qryset->onegap+qryset->onegapcount);//qryset->onegapcount);//( //(qryset->onegapcount);// = onegap_thrust_h;//(qryset->onegap, qryset->onegap+qryset->onegapcount);// = onegap_thrust_h;
/*
   for(int i = 0; i< qryset->onegapcount; i++){
   h_structures[i].qrystart = qryset->onegap[i].qrystart;
   h_structures[i].qrystart_len = qryset->onegap[i].qrystart_len;
   h_structures[i].qryend =qryset->onegap[i].qryend;
   h_structures[i].qryend_len = qryset->onegap[i].qryend_len;
   h_structures[i].gap = qryset->onegap[i].gap;
   h_structures[i].refstr_s = qryset->onegap[i].refstr_s;
   }*/
//cudaFree(raw_ptr);

