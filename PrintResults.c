#include "ComTypes.h"
#include "PrintResults.h"

using namespace std;

#define OUTPUT_BUF_SIZE  32*1024
#define MIN(a, b) (a > b ? b : a)

static char output_buf[OUTPUT_BUF_SIZE];
static int output_buf_write_index = 0;

#define PRINT(format, args...)						\
	do {									\
	int nc = snprintf(&output_buf[output_buf_write_index],		\
	OUTPUT_BUF_SIZE - output_buf_write_index,	\
	format, args);					\
	if (nc >= (OUTPUT_BUF_SIZE - output_buf_write_index)) {		\
	output_buf[output_buf_write_index] = '\0';			\
	printf("%s", output_buf);					\
	output_buf_write_index = 0;					\
	output_buf[0] = '\0';						\
	nc = snprintf(output_buf,					\
	OUTPUT_BUF_SIZE,					\
	format, args);					\
	assert(nc < OUTPUT_BUF_SIZE);					\
	assert(nc >= 0);						\
	}									\
	output_buf_write_index += nc;					\
	} while (0)

#define myassert(cond, msg, var1, var2, var3, var4)	\
	do {							\
	if(!(cond)) {					\
	fprintf(stderr, msg, var1, var2, var3, var4);	\
	assert(0);					\
	}							\
	} while(0)

void printResults_two_file(ref_set_t * refset, qry_set_t * qryset, options_t * options, hashtbl* users, hash_intchar* intchar) {

	mytimer_t t;
	int refindex, j=0;

	output_buf_write_index = 0;
	output_buf[output_buf_write_index] = '\0';

	fprintf(stderr, "\nPrinting results out\n");
	timer_start(&t);
	int upper, lower;
	int qryindex;
	fprintf(stderr, "QRYCOUNT is %d\n", qryset->qryscount);
	for (qryindex = 0; qryindex < qryset->qryscount; qryindex++) {
		PRINT(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> Query: %d\n", qryindex);

		for (refindex = 0; refindex < refset->count; refindex++) {
			ref_t * ref = &refset->refs[refindex];
			int reflen = ref->toklen;
			result_t_two* qryresult = qryset->result_two;
			int end = -1;
			if (qryindex != qryset->qryscount-1){
				end = qryset->qrysoffsettok[qryindex+1];
			} else {
				end = qryset->totaltokens;
			}

			for (j = qryset->qrysoffsettok[qryindex]; j < end; j++) {
				if (qryresult[j].longestmatch > 0) {///TBD
					upper = qryresult[j].up;
					lower = qryresult[j].down;
					int longestmatch = qryresult[j].longestmatch;

					//fprintf(stdout, "lowlen %d - low %d - matchlenup %d - upper %d -	name %s - qrylen %d - reflen %d\n", lowlen, lower, uplen, upper, qryset->qrysname[qryindex],qryset->qryslen[qryindex], reflen);
					if (lower > reflen){
						printf("Lower > reflen? Possible?\n");
						continue;
					}

					PRINT("@Longest Match - %d\n", longestmatch);

					int UPP = ref->str[ref->sa[upper]]; 				
					struct my_struct *s;
					HASH_FIND_INT( intchar, &UPP, s );	/* s: output pointer */
					if(!s){ 			  
						printf("CANNOT FIND this WORD! - IMPOSSIBLE!\n");
					}

					PRINT("%s", s->name);

					int LOW = ref->str[ref->sa[lower]];
					struct my_struct *s1;
					HASH_FIND_INT( intchar, &LOW, s1);	/* s: output pointer */
					if(!s1){			   
						printf("CANNOT FIND this WORD! - IMPOSSIBLE!\n");
					}

					PRINT("||%d %d||%d\n", upper, lower, lower - upper + 1);

					struct my_struct *s2;
					HASH_FIND_INT( intchar, &(qryset->qrysoffsettok[qryset->qryscount + j]), s2);  /* s: output pointer */
					if(!s2){
						printf("CANNOT FIND this WORD! - IMPOSSIBLE!\n");
					} /*else{
					  PRINT("This token %d - %s\n-----------------------------\n", qryset->qrysoffsettok[qryset->qryscount + j], s2->name);
					  }*/

					if(qryset->qrysoffsettok[qryset->qryscount + j] != LOW || qryset->qrysoffsettok[qryset->qryscount + j] != UPP){
						PRINT("XXXXXXXXXXXXXXXX%d - tok %d - start %d\n",qryindex, end-j, j-qryset->qrysoffsettok[qryindex]);
					}
				}
				if (qryresult[j].longestmatch > 1) {					
					int cc=qryset->connectoffset[j];
					int ct;
					for(ct = 2; ct <= qryresult[j].longestmatch; ct++, cc++){
						upper = qryset->result_connect[cc].up;
						lower = qryset->result_connect[cc].down;
						int iii = 0;
						for(;iii < ct; iii++){
							int LOW = ref->str[ref->sa[lower]+iii];
							struct my_struct *s1;
							HASH_FIND_INT( intchar, &LOW, s1);	/* s: output pointer */
							if(!s1){			   
								printf("CANNOT FIND this WORD! - IMPOSSIBLE!\n");
							}
							PRINT("%s ", s1->name);
						}
						PRINT("||%d %d||%d\n", upper, lower, lower - upper + 1);
					}
				}
				
			}			
		}
	}

	/* print out whatever we have left in the buffer */
	if (output_buf_write_index > 0)
		printf("%s", output_buf);

	timer_stop(&t);
	fprintf(stderr, "print out time: %f\n", timer_elapsed(&t)/1000);
}



void printResults_two(ref_set_t * refset, qry_set_t * qryset, options_t * options, hashtbl* users, hash_intchar* intchar) {

	mytimer_t t;
	int refindex, j=0;

	output_buf_write_index = 0;
	output_buf[output_buf_write_index] = '\0';

	fprintf(stderr, "\nPrinting results out\n");
	timer_start(&t);
	int upper, lower, uplen, lowlen;
	int qryindex;
	fprintf(stderr, "QRYCOUNT is %d\n", qryset->qryscount);
	for (qryindex = 0; qryindex < qryset->qryscount; qryindex++) {
		PRINT(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> Query: %d\n", qryindex);

		for (refindex = 0; refindex < refset->count; refindex++) {
			ref_t * ref = &refset->refs[refindex];
			int reflen = ref->toklen;
			result_t_two* qryresult = qryset->result_two;
			int end = -1;
			if (qryindex != qryset->qryscount-1){
				end = qryset->qrysoffsettok[qryindex+1];
			} else {
				end = qryset->totaltokens;
			}
			//fprintf(stderr, "Qind %d - start %d - end %d\n", qryindex, qryset->qrysoffsettok[qryindex], end);
			for (j = qryset->qrysoffsettok[qryindex]; j < end; j++) {
				if (qryresult[j].longestmatch > 0) {///TBD
					upper = qryresult[j].up;
					lower = qryresult[j].down;

					int longestind = qryresult[j].longestind;
					int longestmatch = qryresult[j].longestmatch;

					//fprintf(stdout, "lowlen %d - low %d - matchlenup %d - upper %d -  name %s - qrylen %d - reflen %d\n", lowlen, lower, uplen, upper, qryset->qrysname[qryindex],qryset->qryslen[qryindex], reflen);
					if (lower > reflen){
						printf("Lower > reflen? Possible?\n");
						continue;
					}

					//assert(lower <= reflen);				
					PRINT("Longest Match - %d\n", longestmatch);

					int UPP = ref->str[ref->sa[upper]];					
					struct my_struct *s;
					HASH_FIND_INT( intchar, &UPP, s );  /* s: output pointer */
					if(!s){               
						printf("CANNOT FIND this WORD! - IMPOSSIBLE!\n");
					}

					PRINT("UPP - SA%d\tREF%d\tTOK%d\t%s\n",
						upper,
						UPP,
						j,
						s->name);

					int LOW = ref->str[ref->sa[lower]];
					struct my_struct *s1;
					HASH_FIND_INT( intchar, &LOW, s1);  /* s: output pointer */
					if(!s1){               
						printf("CANNOT FIND this WORD! - IMPOSSIBLE!\n");
					}
					PRINT("LOW - SA%d\tREF%d\tTOK%d\t%s\n",
						lower,
						LOW,
						j,
						s1->name);
					PRINT("Matched Len Matchcount %d - SAUpper %d - SALow %d\n", lower - upper + 1, upper, lower);

					struct my_struct *s2;
					HASH_FIND_INT( intchar, &(qryset->qrysoffsettok[qryset->qryscount + j]), s2);  /* s: output pointer */
					if(!s2){
						printf("CANNOT FIND this WORD! - IMPOSSIBLE!\n");
					} else{
						PRINT("This token %d - %s\n-----------------------------\n", qryset->qrysoffsettok[qryset->qryscount + j], s2->name);
					}

					if(qryset->qrysoffsettok[qryset->qryscount + j] != LOW || qryset->qrysoffsettok[qryset->qryscount + j] != UPP){
						PRINT("XXXXXXXXXXXXXXXX%d - tok %d - start %d\n",qryindex, end-j, j-qryset->qrysoffsettok[qryindex]);
					}
				}
				if (qryresult[j].longestmatch > 1) {					
					int cc=qryset->connectoffset[j];
					int ct;
					for(ct = 2; ct <= qryresult[j].longestmatch; ct++, cc++){
						upper = qryset->result_connect[cc].up;
						lower = qryset->result_connect[cc].down;
						int iii = 0;
						for(;iii < ct; iii++){
							int LOW = ref->str[ref->sa[lower]+iii];
							struct my_struct *s1;
							HASH_FIND_INT( intchar, &LOW, s1);	/* s: output pointer */
							if(!s1){			   
								printf("CANNOT FIND this WORD! - IMPOSSIBLE!\n");
							}
							PRINT("%s ", s1->name);
						}
						PRINT("|| SAUpper %d - SALow %d || Matched Len Matchcount %d\n", upper, lower, lower - upper + 1);
					}
				}
			}			
		}
	}

	/* print out whatever we have left in the buffer */
	if (output_buf_write_index > 0)
		printf("%s", output_buf);

	timer_stop(&t);
	fprintf(stderr, "print out time: %f\n", timer_elapsed(&t)/1000);
}

void print_query_GPU_Continous( 		
		vector<vector<unsigned int> > queryglobal, 
		int qrycounter,
		result_t* globalOnPairsUpDown,
		red_dup_t* fast_speed,
		char* destDir) {		
	FILE *fp;
	int i = 0;
	int ii = 0;
	int n;	
	int p;
	mytimer_t t;
	timer_start(&t);

	char filename[100];
	cerr<<"Start Printing NonGappy Phrases..."<<endl;

	for(; ii < qrycounter; ii++){			
		if (isSample){			
			//sprintf(filename,"/scratch0/huah/gpugrammar_sample_loose_fbis/grammar.%d",ii);
			sprintf(filename,"%s/grammar.noGap.%d.s", destDir, ii);
		} else {
			//sprintf(filename,"/scratch0/huah/gpugrammar_nonsample_loose_fbis/grammar.%d",ii);//dev/test
			sprintf(filename,"%s/grammar.noGap.%d.n", destDir, ii);
		}

		/*for(; ii < qrycounter; ii++){ 		
		if (isSample){
			sprintf(filename,"/scratch0/huah/gpugrammar_sample_loose_fbis/grammar.%d",ii);
		} else {
			sprintf(filename,"/scratch0/huah/gpugrammar_nonsample_loose_fbis/grammar.%d",ii);//dev/test
		}*/
		//ofstream out(filename, std::ios::out | std::ios::trunc);

		fp = fopen(filename, "w");
		if (fp == NULL ){
			fprintf(stderr, "Please check your file directory address for grammar rule files output. It is not valid. Program Exits.\n");
			exit(0);
		}		 

		//fprintf(fp,"[X] ||| [X,1] [X,2] ||| [X,1] [X,2] ||| IsMonotone=1.0\n[X] ||| [X,1] [X,2] ||| [X,2] [X,1] ||| IsSwapped=1.0\n");
		//cerr<<ii<<endl;

		for(int ic=0; ic< queryglobal[ii].size(); ic++){ //p=(int*)utarray_front(queryglobal[ii]);p!=NULL;p=(int*)utarray_next(queryglobal[ii], p)) {
			p = queryglobal[ii][ic];
			//fprintf(stderr, "NEW BLOCk START! %d - qry c %d\n", *p, ii);
			if (globalOnPairsUpDown[p].down == -1 || globalOnPairsUpDown[p].up == -1){
				continue;
			}
			for(i = globalOnPairsUpDown[p].down; i<= globalOnPairsUpDown[p].up; i++) { // < or <= ?

				//fprintf(fp, "%s, Score %f, c(e,f) %d, c(e) %d, c(f) %d\n", getcounter_target3->name, -log10((double)getcounter_target3->id/(double)getcounter_target->id), getcounter_target3->id, getcounter_target->id, getcounter_target2->id);					
				//fprintf(fp, "%s|||Score %f, c(e,f) %d, c(f) %d\n", getcounter_target3->name, -log10((double)getcounter_target3->id/(double)getcounter_target2->id), getcounter_target3->id, getcounter_target2->id);
				fprintf(fp, 
						"[X] ||| %s ||| EgivenFCoherent=%f SampleCountF=%f CountEF=%f MaxLexFgivenE=%f MaxLexEgivenF=%f IsSingletonF=%d IsSingletonFE=%d\n",
						fast_speed[i].lexical, 
						fast_speed[i].aa,
						fast_speed[i].all_suffix_fsample_score,
						fast_speed[i].bb,
						fast_speed[i].MaxLexFgivenE,
						fast_speed[i].MaxLexEgivenF,
						(fast_speed[i].f == 1),
						(*fast_speed[i].paircount) == 1);				
				/*	 fprintf(fp, "[X] ||| %s ||| %f %f %f %f %f %d %d\n", 
					 fast_speed[i].lexical,
					 fast_speed[i].aa,	
					 fast_speed[i].all_suffix_fsample_score,
					 fast_speed[i].bb,
					 fast_speed[i].MaxLexFgivenE,
					 fast_speed[i].MaxLexEgivenF,
					 (fast_speed[i].f == 1),
					 (*fast_speed[i].paircount)== 1);*/

			}
		}
		fclose(fp); 		
		n = 0;
	}
	timer_stop(&t);
	fprintf(stderr, "Print out time: %f\n", timer_elapsed(&t)/1000);
}

void printGapMode(
		int gapMode,	
		unsigned int id, 
		result_t* globalOnPairsUpDownContinous,
		result_t* globalOnPairsUpDownGappy,
		result_t* globalOnPairsUpDownTwoGap,
		red_dup_t* fast_speed,
		red_dup_t* fast_speed_one,
		red_dup_t* fast_speed_two,
		FILE* fp){

	if(gapMode == 0){
		if (globalOnPairsUpDownContinous[id].down != -1&&globalOnPairsUpDownContinous[id].up != -1){		
			//Deal with real ab phrase	
			for(int i = globalOnPairsUpDownContinous[id].down; i<= 
					globalOnPairsUpDownContinous[id].up; i++) { // < or <= ?
				fprintf(fp, 
						"[X] ||| %s ||| EgivenFCoherent=%f SampleCountF=%f CountEF=%f MaxLexFgivenE=%f MaxLexEgivenF=%f IsSingletonF=%d IsSingletonFE=%d\n",
						fast_speed[i].lexical, 
						fast_speed[i].aa,
						fast_speed[i].all_suffix_fsample_score,
						fast_speed[i].bb,
						fast_speed[i].MaxLexFgivenE,
						fast_speed[i].MaxLexEgivenF,
						(fast_speed[i].f == 1),
						(*fast_speed[i].paircount) == 1);								
			}
		}
	} else if (gapMode == 1){
		if(globalOnPairsUpDownGappy[id].down!=-1&&globalOnPairsUpDownGappy[id].up	!=-1){
			for(int i = globalOnPairsUpDownGappy[id].down; i<= globalOnPairsUpDownGappy[id].up; i++) { 
				fprintf(fp, 
						"[X] ||| %s ||| EgivenFCoherent=%f SampleCountF=%f CountEF=%f MaxLexFgivenE=%f MaxLexEgivenF=%f IsSingletonF=%d IsSingletonFE=%d\n",
						//fast_speed_one[i].blocknumber, 
						fast_speed_one[i].lexical, 
						fast_speed_one[i].aa,
						fast_speed_one[i].all_suffix_fsample_score,
						fast_speed_one[i].bb,
						fast_speed_one[i].MaxLexFgivenE,
						fast_speed_one[i].MaxLexEgivenF,
						(fast_speed_one[i].f == 1),
						(*fast_speed_one[i].paircount) == 1);								
			}
		}
	} else if (gapMode == 2) {
		if(globalOnPairsUpDownTwoGap[id].down!=-1&&globalOnPairsUpDownTwoGap[id].up!=-1){			
			for(int i = globalOnPairsUpDownTwoGap[id].down; i<= globalOnPairsUpDownTwoGap[id].up; i++) { 
				fprintf(fp, 
						"[X] ||| %s ||| EgivenFCoherent=%f SampleCountF=%f CountEF=%f MaxLexFgivenE=%f MaxLexEgivenF=%f IsSingletonF=%d IsSingletonFE=%d\n",
						//fast_speed_two[i].blocknumber, 
						fast_speed_two[i].lexical, 
						fast_speed_two[i].aa,
						fast_speed_two[i].all_suffix_fsample_score,
						fast_speed_two[i].bb,
						fast_speed_two[i].MaxLexFgivenE,
						fast_speed_two[i].MaxLexEgivenF,
						(fast_speed_two[i].f == 1),
						(*fast_speed_two[i].paircount) == 1);								
			}
		}
	} else {
		fprintf(stderr, "No such gap mode!\n");
		exit(0);
	}


}

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
		unsigned int distinctTwoGapCount) { 	  			

	FILE *fp;
	int i = 0;
	int ii = 0;
	int n;	
	int p;
	mytimer_t t;
	timer_start(&t);

	char filename[150];
	cerr<<"Start Printing Gappy Phrases..."<<endl;

	for(; ii < qrycounter; ii++){			
		if (isSample){			
			//sprintf(filename,"/scratch0/huah/gpugrammar_sample_loose_fbis/grammar.%d",ii);
			sprintf(filename,"%s/grammar.%d.s", destDir, ii);
		} else {
			//sprintf(filename,"/scratch0/huah/gpugrammar_nonsample_loose_fbis/grammar.%d",ii);//dev/test
			sprintf(filename,"%s/grammar.%d.n", destDir, ii);
		}
		fp = fopen(filename, "w");
		if (fp == NULL ){
			fprintf(stderr, "Please check your file directory address for grammar rule files output. It is not valid. Program Exits.\n");
			exit(0);
		}
		//fprintf(fp,"[X] ||| [X,1] [X,2] ||| [X,1] [X,2] ||| IsMonotone=1.0\n[X] ||| [X,1] [X,2] ||| [X,2] [X,1] ||| IsSwapped=1.0\n");
		//cerr<<ii<<endl;

		////Continous One - bnum - ab
		for(int ic=0; ic< queryglobal[ii].size(); ic++){ 
		//for(p=(int*)utarray_front(queryglobal[ii]);p!=NULL;p=(int*)utarray_next(queryglobal[ii], p)) {			
			//fprintf(stderr, "NEW BLOCk START! %d - qry c %d\n", *p, ii);
			//Deal with gappy abX, Xab, and XabX phrases			
			p = queryglobal[ii][ic];
			//abX
			printGapMode(
					1, 
					p+global, 
					globalOnPairsUpDownContinous,
					globalOnPairsUpDownGappy,
					globalOnPairsUpDownTwoGap,
					fast_speed,
					fast_speed_one,
					fast_speed_two,
					fp);	

			//Xab
			printGapMode(
					1, 
					p, 
					globalOnPairsUpDownContinous,
					globalOnPairsUpDownGappy,
					globalOnPairsUpDownTwoGap,
					fast_speed,
					fast_speed_one,
					fast_speed_two,
					fp);	

			//XabX
			printGapMode(
					2, 
					p, 
					globalOnPairsUpDownContinous,
					globalOnPairsUpDownGappy,
					globalOnPairsUpDownTwoGap,
					fast_speed,
					fast_speed_one,
					fast_speed_two,
					fp);	

			//ab
			printGapMode(
					0, 
					p, 
					globalOnPairsUpDownContinous,
					globalOnPairsUpDownGappy,
					globalOnPairsUpDownTwoGap,
					fast_speed,
					fast_speed_one,
					fast_speed_two,
					fp);			
		}

		////One Gap Phrase - aXb
		int gapId = -1;
		for(int si = 0; si < queryIdGappy[ii].size(); si++){		
			gapId = 2*global+queryIdGappy[ii][si];
			//aXb
			printGapMode(
					1, 
					gapId, 
					globalOnPairsUpDownContinous,
					globalOnPairsUpDownGappy,
					globalOnPairsUpDownTwoGap,
					fast_speed,
					fast_speed_one,
					fast_speed_two,
					fp);	

			gapId = global + distinctTwoGapCount + queryIdGappy[ii][si];
			//XaXb
			printGapMode(
					2, 
					gapId, 
					globalOnPairsUpDownContinous,
					globalOnPairsUpDownGappy,
					globalOnPairsUpDownTwoGap,
					fast_speed,
					fast_speed_one,
					fast_speed_two,
					fp);	

			gapId = global + distinctTwoGapCount + distinctOneGapCount + queryIdGappy[ii][si];
			
			//Debug
			if(gapId >= global + distinctTwoGapCount + 2*distinctOneGapCount){
				fprintf(stderr, "Error: Make sure your gap Id is within the range!\n");
				return;
			}
			//aXbX
			printGapMode(
					2, 
					gapId, 
					globalOnPairsUpDownContinous,
					globalOnPairsUpDownGappy,
					globalOnPairsUpDownTwoGap,
					fast_speed,
					fast_speed_one,
					fast_speed_two,
					fp);
		}
		////Two Gap Phrase - aXbXc
		for(int si = 0; si < queryIdTwoGap[ii].size(); si++){		
			
			gapId = global+queryIdTwoGap[ii][si];			
			/*if(ii==0){
				fprintf(stderr, "Two gap %d | original %d |||", gapId, queryIdTwoGap[ii][si]);
			}*/
			printGapMode(
					2, 
					gapId, 
					globalOnPairsUpDownContinous,
					globalOnPairsUpDownGappy,
					globalOnPairsUpDownTwoGap,
					fast_speed,
					fast_speed_one,
					fast_speed_two,
					fp);	
		}

		fclose(fp);
		n = 0;
	}
	timer_stop(&t);
	fprintf(stderr, "Print out time: %f\n", timer_elapsed(&t)/1000);
}


