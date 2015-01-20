#include "ComTypes.h"
#include <iostream>
#include <map>
#include "Timer.h"
#include "Disk.h"
#include "uthash/uthash.h"
#include "uthash/utarray.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <assert.h>
#include <algorithm>
#include <vector>

using namespace std;

red_dup_t* createLexicon(
		int count, 
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
		char** target_vocabulary) {

	mytimer_t t;
	timer_start(&t);
	int i = 0;
	int j = 0;
	bool flag = false;
	char str1[250];
	char combine[500];
	char tar_second_str[450];
	map<string,categ>::iterator it;
	(*lexicon_index) = 0;
	int* fsample_arr = new int[global];    
	memset(fsample_arr, 0, sizeof(int)*global);    
	fprintf(stderr, "Start continous lexicon creation process\n");

	int source_starter = -1;
	int source_ender = -1;
	for(; i < count; i++){		
		fsample_arr[out_res[i].blocknumber]++;
	}

	i = 0;
	//fprintf(stderr, "%s - go\n",sourceName[2]);
	for(; i < count; i++){
		flag = false;
		strcpy(combine, ""); 
		//cerr << "Count "<<i<<endl;
		source_starter = tmp_blocks[out_res[i].blocknumber].string_start;
		source_ender = source_starter + tmp_blocks[out_res[i].blocknumber].matchlen -1;
		int b = out_res[i].tar_start + out_res[i].tar_end;
		//cerr<<"source start\n";
		strcpy(tar_second_str, "| "); 		
		int jj = out_res[i].tar_start;
		//fprintf(stderr, "target start i %d - jj %d - b %d\n", i, jj,b);
		for(; jj <= b; jj++){
			char* directAcc = target_vocabulary[target_str[jj]];									  
			
			//cerr<<"Inside"<<endl;
			//fprintf(stderr, "before cat s1->name %s | str1 %s | tar_second_str %s\n", s1->name, str1, tar_second_str);
			if (jj == out_res[i].tar_start){
				sprintf(str1,"%s",directAcc);
			} else {
				sprintf(str1," %s",directAcc);
			}			
			strcat(tar_second_str, str1);
		}
		if(flag){
			continue;
		}

		///Collect Phrase Count of c(e,f)
		//cerr<<"Yoshi"<<out_res[i].blocknumber<<endl;
		strcat(combine, sourceName[out_res[i].blocknumber]);
		//cerr<<"Yoshi2"<<endl;//gogo
		strcat(combine, " ||" );
		strcat(combine, tar_second_str);

		struct my_struct3 *sss;
		HASH_FIND_STR(*lexic, combine, sss);            
		if(!sss){
			hashtbl_aux* s = (hashtbl_aux*)malloc(sizeof(hashtbl_aux));   
			s->id = (int*)malloc(sizeof(int));
			(*s->id) = 1;
			//cerr<<"c(e, f) inside before malloc start\n";
			s->name = (char*)malloc(sizeof(char)*(strlen(combine)+1)); 
			memcpy((void*)s->name, combine, strlen(combine)+1);
			HASH_ADD_KEYPTR( hh, *lexic, s->name, strlen(s->name), s);  
			//cerr<<"c(e, f) inside after malloc start\n";

			//for f in fphrase.words:
			//		for e in ewords:
			//			score = self.ttable.get_score(f, e, 1)
			//		if score > max_score:
			//			max_score = score
			//mlfe += -log10(max_score) if max_score > 0 else MAXSCORE

			///////////MaxLexFgivenE 
			float fgivene = 0;
			float egivenf = 0;
			j = 0;
			jj = 0;
			float max_scorer_1 = 0.0;
			float max_scorer_2 = 0.0;
			stringstream oos1;            
			for(j =source_starter; j <=  source_ender; j++){///	fgivene
				max_scorer_1 = 0;
				for(jj = out_res[i].tar_start; jj <= b; jj++){				
					if (jj == out_res[i].tar_start) {
						oos1 << source_str[j] << "|-1";	
						it = word_score.find(oos1.str());
						if (it != word_score.end() && it->second.val2 > max_scorer_1){
							max_scorer_1 = it->second.val2;
						}
						oos1.str("");
					}
					oos1 << source_str[j] << "|"<< target_str[jj];
					it = word_score.find(oos1.str());
					if (it != word_score.end() && it->second.val2 > max_scorer_1){
						max_scorer_1 = it->second.val2;
					}
					oos1.str("");
				}
				if (max_scorer_1 > 0){
					fgivene += -log10(max_scorer_1);
				} else {
					fgivene += MAXSCORE;
				}
			}		

			///////////MaxLexEgivenF
			j = 0;
			jj = 0;
			oos1.str("");
			for(jj = out_res[i].tar_start; jj <= b; jj++){///            		
				max_scorer_2 = 0;
				for(j = source_starter; j <=  source_ender; j++){
					if (j == source_starter) {
						oos1 << "-1|" << target_str[jj];	
						it = word_score.find(oos1.str());
						if (it != word_score.end() && it->second.val1 > max_scorer_2){
							max_scorer_2 = it->second.val1;
						}
						oos1.str("");
					}
					oos1 << source_str[j] << "|"<< target_str[jj];
					it = word_score.find(oos1.str());
					if (it != word_score.end() && it->second.val1 > max_scorer_2){
						max_scorer_2 = it->second.val1;
					}
					oos1.str("");
				}
				if (max_scorer_2 > 0){
					egivenf += -log10(max_scorer_2);
				} else {
					egivenf += MAXSCORE;
				}
			}        

			////Speed UP!
			fast_speed[*lexicon_index].blocknumber = out_res[i].blocknumber;
			fast_speed[*lexicon_index].lexical = s->name;  
			fast_speed[*lexicon_index].all_suffix_fsample = 1 + tmp_blocks[out_res[i].blocknumber].end - tmp_blocks[out_res[i].blocknumber].start;
			fast_speed[*lexicon_index].all_suffix_fsample_score = log10(1 + (fast_speed[*lexicon_index].all_suffix_fsample));
			fast_speed[*lexicon_index].f = fsample_arr[out_res[i].blocknumber];
			fast_speed[*lexicon_index].paircount = s->id;
			fast_speed[*lexicon_index].MaxLexEgivenF = egivenf;
			fast_speed[*lexicon_index].MaxLexFgivenE = fgivene;
			(*lexicon_index)++;
		} else {
			(*sss->id)++;
		}
	}
	//fprintf(stderr, "OK\n");
	for(int i = 0; i< *lexicon_index; i++){
		fast_speed[i].aa = 
			-log10((float)(*fast_speed[i].paircount)/(float)fast_speed[i].all_suffix_fsample);
		fast_speed[i].bb = log10(1 + (*fast_speed[i].paircount));
	}
	timer_stop(&t);
	fprintf(stderr, "Create continous lexicon time: %f\n", timer_elapsed(&t)/1000);
	delete[] fsample_arr;	
	//free(tmp_blocks);
	return fast_speed;
}

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
		two_gappy_search* twoGapSearch) {

			fprintf(stderr, "Start one gap lexicon creation process\n");
			mytimer_t t;
			timer_start(&t);

			bool flag = false;
			char str1[250];
			char combine[500];
			char tar_second_str[450];
			char sourceStr[250];
			map<string,categ>::iterator it;
			(*lexicon_index) = 0;
			int source_starter = -1;
			int source_ender = -1;

			unsigned int position = -1;
			unsigned int jj;
			int j;

			map<string, int*> deDuplication;
			map<string, int*>::iterator iteratorDeDup;	
			char* directAcc = NULL;

			int* pattern;
			int number;
			vector<int> targetStringIds;
			vector<int> sourceStringIds;

			int* fsample_arr = new int[2*global+distinctOneGapCount];    
			memset(fsample_arr, 0, sizeof(int)*(2*global+distinctOneGapCount));    

			//fsample feature
			unsigned int convertedId = 0;
			for(unsigned int i =0; i < count_one_gap; i++){		
				if(i < seperatorOneGap){
					convertedId = oneGapRule[i].gappy_index;			
				} else {
					convertedId = 2*global + oneGapRule[i].gappy_index;
				}
				fsample_arr[convertedId]++;
			}

			for(unsigned int i = 0; i < count_one_gap; i++){
				//cerr << "Count "<<i<<endl;
				//cerr<<"Yoshi"<<out_res[i].blocknumber<<endl;
				flag = false;
				strcpy(combine, ""); 		
				directAcc = NULL;
				//Collect Source String
				if(i==0 || oneGapRule[i].gappy_index != oneGapRule[i-1].gappy_index){			
					strcpy(sourceStr, ""); 
					sourceStringIds.clear();
					if(i < seperatorOneGap){				
						convertedId = oneGapRule[i].gappy_index;
						if(oneGapRule[i].gappy_index<global){
							//Xab - need shift by one and -1
							strcat(sourceStr, "[X,1] ");		
							strcat(sourceStr, sourceName[convertedId]);	
							source_starter = tmp_blocks[convertedId].string_start;
							source_ender = source_starter + tmp_blocks[convertedId].matchlen -1;								
						} else {
							//abX - need shift by one
							strcat(sourceStr, sourceName[convertedId-global]);			
							strcat(sourceStr, " [X,1]");
							source_starter = tmp_blocks[convertedId-global].string_start;
							source_ender = source_starter + tmp_blocks[convertedId-global].matchlen -1;
						}
						
						//fprintf(stderr, "i->%d\n",i);
						for(int iccc = source_starter; iccc <= source_ender; iccc++){
							sourceStringIds.push_back(source_str[iccc]);
						}
						//fprintf(stderr, "\n");
					} else {
						//Debug
						if(oneGapRule[i].gappy_index<0){
							printf("Not possible for this oneGap rule|gappy_index %d, count %d\n", 
									oneGapRule[i].gappy_index, i);
							exit(0);
						}
						///End
						convertedId = 2*global + oneGapRule[i].gappy_index;

						//Deal with aXb case, where ID should be positive always. No shift needed
						number = onegapPattern[oneGapSearch[oneGapRule[i].gappy_index].position].number;
						pattern = onegapPattern[oneGapSearch[oneGapRule[i].gappy_index].position].pattern;

						for(jj = 0; jj < number; jj++){
							if(pattern[jj]>=0){
								sourceStringIds.push_back(pattern[jj]);
								directAcc = source_vocabulary[pattern[jj]];	
								//Debug purpose on source side
								//End debug 						
								if(directAcc == NULL){
									fprintf(stderr, "Please.. i - %d\n", i);
								}								
							} else {
								directAcc = "[X,1]";
							}

							if (jj == 0){
								sprintf(str1,"%s",directAcc);
							} else {
								sprintf(str1," %s",directAcc);
							}			
							strcat(sourceStr, str1);				
						}
					}

					//Clean the hash table
					deDuplication.clear();					
				}

				unsigned int targetStart = oneGapRule[i].ref_str_start;
				unsigned int targetEnd = oneGapRule[i].ref_str_start + oneGapRule[i].end;
				unsigned int gap1Start = oneGapRule[i].ref_str_start + oneGapRule[i].gap1;
				unsigned int gap1End = oneGapRule[i].ref_str_start + oneGapRule[i].gap1_1;

				///get Target String aXb or abX or Xab;
				strcpy(tar_second_str, " ||| "); 		        
				directAcc = NULL;
				//fprintf(stderr, "target start i %d - jj %d - b %d\n", i, jj,b);
				targetStringIds.clear();

				for(jj = targetStart; jj <= targetEnd; jj++){
					if (jj >= gap1Start && jj <= gap1End){
						directAcc = "[X,1]";				
					} else {
						directAcc = target_vocabulary[target_str[jj]];			
						targetStringIds.push_back(target_str[jj]);

						//Debug purpose
						if(directAcc == NULL){
							fprintf(stderr, "NULL direct access\n");
							exit(0);
						}
					}			

					if (jj == targetStart){
						sprintf(str1,"%s",directAcc);
					} else {
						sprintf(str1," %s",directAcc);
					}			
					strcat(tar_second_str, str1);

					if (jj >= gap1Start && jj <= gap1End){
						jj = gap1End;
					}
				}

				if(flag){
					continue;
				}

				///Collect Phrase Count of c(e,f)
				//cerr<<"Yoshi2"<<endl;//gogo
				strcat(combine, sourceStr);
				strcat(combine, tar_second_str);
				/*if(out_res[i].blocknumber == 15){
				  cerr<<"c(e, f) e f"<<combine<<"  count: "<<i<<" blocknumber: "<<out_res[i].blocknumber<<endl;
				  }*/
				std::string targetStr(tar_second_str);
				//struct my_struct3 *sss;
				//HASH_FIND_STR(*lexic, tar_second_str, sss);            
				iteratorDeDup = deDuplication.find(targetStr);

				if(iteratorDeDup == deDuplication.end()){
					//hashtbl_aux* s = (hashtbl_aux*)malloc(sizeof(hashtbl_aux));   
					int* counter = (int*)malloc(sizeof(int));
					(*counter) = 1;

					char* lexString = (char*)malloc(sizeof(char)*(strlen(combine)+1)); 
					memcpy((void*)lexString, combine, strlen(combine)+1);
					//cerr<<"c(e, f) inside after malloc start\n";

					//Insert to map			
					deDuplication.insert(pair<string,int*>(targetStr, counter));

					//for f in fphrase.words:
					//			for e in ewords:
					//			score = self.ttable.get_score(f, e, 1)
					//		if score > max_score:
					//		max_score = score
					//mlfe += -log10(max_score) if max_score > 0 else MAXSCORE

					///////////MaxLexFgivenE 
					float fgivene = 0;
					float egivenf = 0;
					j = 0;
					jj = 0;
					float max_scorer_1 = 0.0;
					float max_scorer_2 = 0.0;
					stringstream oos1;            
					for(j =0; j < sourceStringIds.size(); j++){///	fgivene
						max_scorer_1 = 0;
						for(jj = 0; jj < targetStringIds.size(); jj++){	
							if (jj == 0) {
								oos1 << sourceStringIds[j] << "|-1";	
								it = word_score.find(oos1.str());
								if (it != word_score.end() && it->second.val2 > max_scorer_1){
									max_scorer_1 = it->second.val2;
								}
								oos1.str("");
							}
							oos1 << sourceStringIds[j] << "|"<< targetStringIds[jj];
							it = word_score.find(oos1.str());
							if (it != word_score.end() && it->second.val2 > max_scorer_1){
								max_scorer_1 = it->second.val2;
							}
							oos1.str("");
						}
						if (max_scorer_1 > 0){
							fgivene += -log10(max_scorer_1);
						} else {
							fgivene += MAXSCORE;
						}
					}		

					///////////MaxLexEgivenF
					j = 0;
					jj = 0;
					oos1.str("");
					for(jj = 0; jj < targetStringIds.size(); jj++){				
						max_scorer_2 = 0;
						for(j = 0; j < sourceStringIds.size(); j++){
							if (j == 0) {
								oos1 << "-1|" << targetStringIds[jj];	
								it = word_score.find(oos1.str());
								if (it != word_score.end() && it->second.val1 > max_scorer_2){
									max_scorer_2 = it->second.val1;
								}
								oos1.str("");
							}
							oos1 << sourceStringIds[j] << "|"<< targetStringIds[jj];
							it = word_score.find(oos1.str());
							if (it != word_score.end() && it->second.val1 > max_scorer_2){
								max_scorer_2 = it->second.val1;
							}
							oos1.str("");
						}
						if (max_scorer_2 > 0){
							egivenf += -log10(max_scorer_2);
						} else {
							egivenf += MAXSCORE;
						}
					}        

					////Speed UP!
					fast_speed[*lexicon_index].blocknumber = convertedId;            
					if(i < seperatorOneGap){
						int realId = convertedId;
						if(convertedId>=global){
							realId = convertedId - global;
						} else if (convertedId >= 2*global){
							//Debug purpose
							printf("Not correct!! Converted Id in gappy lexicon\n");
							exit(0);
							//Debug end
						}
						
						fast_speed[*lexicon_index].all_suffix_fsample = 1 + 
							tmp_blocks[realId].end - 
							tmp_blocks[realId].start;						
					} else {
						fast_speed[*lexicon_index].all_suffix_fsample = 1 + 
							oneGapSearch[oneGapRule[i].gappy_index].end_on_salist - 
							oneGapSearch[oneGapRule[i].gappy_index].start_on_salist;
					}
					fast_speed[*lexicon_index].f = fsample_arr[convertedId];
					fast_speed[*lexicon_index].lexical = lexString;  
					fast_speed[*lexicon_index].all_suffix_fsample_score = log10(1 + (fast_speed[*lexicon_index].all_suffix_fsample));
					fast_speed[*lexicon_index].paircount = counter;
					fast_speed[*lexicon_index].MaxLexEgivenF = egivenf;
					fast_speed[*lexicon_index].MaxLexFgivenE = fgivene;

					(*lexicon_index)++;
				} else {
					(*(iteratorDeDup->second))++;
				}
			}
			//fprintf(stderr, "OK\n");
			for(int i = 0; i< *lexicon_index; i++){
				fast_speed[i].aa = 
					-log10((float)(*fast_speed[i].paircount)/(float)fast_speed[i].all_suffix_fsample);
				fast_speed[i].bb = log10(1 + (*fast_speed[i].paircount));
			}
			timer_stop(&t);
			fprintf(stderr, "Create One Gap Lexicon time: %f\n", timer_elapsed(&t)/1000);

			delete[] fsample_arr;	
			return fast_speed;
		}

		red_dup_t* createLexiconFast(
				int count, 
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
				unsigned int* lexicalTaskCounter) {
		
			mytimer_t t;
			timer_start(&t);
			int i = 0;
			int j = 0;
			bool flag = false;
			char str1[400];
			char combine[700];
			char tar_second_str[650];
			map<string,categ>::iterator it;
			(*lexicon_index) = 0;
			int* fsample_arr = new int[global];    
			memset(fsample_arr, 0, sizeof(int)*global);    
			fprintf(stderr, "Start continous lexicon creation process\n");
		
			int source_starter = -1;
			int source_ender = -1;
			for(; i < count; i++){		
				fsample_arr[out_res[i].blocknumber]++;
			}
		
			i = 0;
			//fprintf(stderr, "%s - go\n",sourceName[2]);
			for(; i < count; i++){
				flag = false;
				strcpy(combine, ""); 
				//cerr << "Count "<<i<<endl;
				source_starter = tmp_blocks[out_res[i].blocknumber].string_start;
				source_ender = source_starter + tmp_blocks[out_res[i].blocknumber].matchlen -1;
				int b = out_res[i].tar_start + out_res[i].tar_end;
				//cerr<<"source start\n";
				strcpy(tar_second_str, "| ");		
				int jj = out_res[i].tar_start;
				//fprintf(stderr, "target start i %d - jj %d - b %d\n", i, jj,b);
				for(; jj <= b; jj++){
					//Debug
					char* directAcc = target_vocabulary[target_str[jj]];									  
					//End of debug					
					if (jj == out_res[i].tar_start){						
						sprintf(str1,"%s",directAcc);
					} else {
						sprintf(str1," %s",directAcc);
					}								
					//cerr<<"Yoshi0"<<endl;
					strcat(tar_second_str, str1);
				}
				if(flag){
					continue;
				}
		
				///Collect Phrase Count of c(e,f)
				//cerr<<"Yoshi"<<out_res[i].blocknumber<<endl;
				strcat(combine, sourceName[out_res[i].blocknumber]);
				//cerr<<"Yoshi2"<<endl;//gogo
				strcat(combine, " ||" );
				//cerr<<"Yoshi2.1"<<endl;//gogo
				strcat(combine, tar_second_str);

				struct my_struct3 *sss;
				HASH_FIND_STR(*lexic, combine, sss);			
				if(!sss){
					hashtbl_aux* s = (hashtbl_aux*)malloc(sizeof(hashtbl_aux));   
					s->id = (int*)malloc(sizeof(int));
					(*s->id) = 1;
					//cerr<<"c(e, f) inside before malloc start\n";
					s->name = (char*)malloc(sizeof(char)*(strlen(combine)+1)); 
					memcpy((void*)s->name, combine, strlen(combine)+1);
					HASH_ADD_KEYPTR( hh, *lexic, s->name, strlen(s->name), s);	
					
					//for f in fphrase.words:
					//		for e in ewords:
					//			score = self.ttable.get_score(f, e, 1)
					//		if score > max_score:
					//			max_score = score
					//mlfe += -log10(max_score) if max_score > 0 else MAXSCORE

					///////////Offload this task to GPU
					maxEF[*lexicalTaskCounter].fastSpeedId = *lexicon_index;
					maxEF[*lexicalTaskCounter].sourcePatternCounter = tmp_blocks[out_res[i].blocknumber].matchlen;							
					for(j =source_starter; j <= source_ender; j++){///	
						maxEF[*lexicalTaskCounter].sourcePattern[j-source_starter] = source_str[j];
					}
					maxEF[*lexicalTaskCounter].targetStart = out_res[i].tar_start;
					maxEF[*lexicalTaskCounter].end = out_res[i].tar_end;

					//Debug
			/*		if(out_res[i].blocknumber == 0){
						printf("lexicalTaskCounter %d| count %d | start %d | end %d | %s\n", 
							*lexicalTaskCounter, 
							maxEF[*lexicalTaskCounter].sourcePatternCounter,
							source_starter,
							source_ender,
							combine);
						for(j =source_starter; j <= source_ender; j++){///	
							printf("Source pattern: %d\n", maxEF[*lexicalTaskCounter].sourcePattern[j-source_starter]);
						}			
						for(j =maxEF[*lexicalTaskCounter].targetStart; j <= maxEF[*lexicalTaskCounter].end; j++){///	
							printf("Target pattern: %d\n", target_str[j]);
						}			
					}*/										

					(*lexicalTaskCounter)++;
							
					////Speed UP!
					fast_speed[*lexicon_index].blocknumber = out_res[i].blocknumber;
					fast_speed[*lexicon_index].lexical = s->name;  
					fast_speed[*lexicon_index].all_suffix_fsample = 1 + tmp_blocks[out_res[i].blocknumber].end - tmp_blocks[out_res[i].blocknumber].start;
					if(isSample&&fast_speed[*lexicon_index].all_suffix_fsample > SAMPLER){
						fast_speed[*lexicon_index].all_suffix_fsample = SAMPLER;
					}
					fast_speed[*lexicon_index].all_suffix_fsample_score = log10(1 + (fast_speed[*lexicon_index].all_suffix_fsample));
					fast_speed[*lexicon_index].f = fsample_arr[out_res[i].blocknumber];
					fast_speed[*lexicon_index].paircount = s->id;
					//fast_speed[*lexicon_index].MaxLexEgivenF = egivenf;
					//fast_speed[*lexicon_index].MaxLexFgivenE = fgivene;
					(*lexicon_index)++;
				} else {
					(*sss->id)++;
				}
			}
			//fprintf(stderr, "OK\n");
			for(int i = 0; i< *lexicon_index; i++){
				fast_speed[i].aa = 
					-log10((float)(*fast_speed[i].paircount)/(float)fast_speed[i].all_suffix_fsample);
				fast_speed[i].bb = log10(1 + (*fast_speed[i].paircount));
			}
			timer_stop(&t);
			fprintf(stderr, "Create continous lexicon time: %f\n", timer_elapsed(&t)/1000);
			delete[] fsample_arr;	
			//free(tmp_blocks);
			return fast_speed;
		}
		
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
				int seperatorOneGap,//Seperator between aXb; Xab, abX
				int seperatorTwoGap_one,//XabX;aXbXc
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
				ref_t* refSource) {
		
					fprintf(stderr, "Start one gap lexicon creation process\n");
					mytimer_t t;
					timer_start(&t);
		
					bool flag = false;
					char str1[250];
					char combine[500];
					char tar_second_str[450];
					char sourceStr[250];
					map<string,categ>::iterator it;
					(*lexicon_index) = 0;
					int source_starter = -1;
					int source_ender = -1;
		
					unsigned int position = -1;
					unsigned int jj;
					int j;
		
					map<string, int*> deDuplication;
					map<string, int*>::iterator iteratorDeDup;	
					char* directAcc = NULL;
		
					int* pattern;
					int number;
					vector<int> sourceStringIds;
		
					int* fsample_arr = new int[2*global+distinctOneGapCount];	 
					memset(fsample_arr, 0, sizeof(int)*(2*global+distinctOneGapCount));    
		
					//fsample feature
					unsigned int convertedId = 0;
					for(unsigned int i =0; i < count_one_gap; i++){ 	
						if(i < seperatorOneGap){
							convertedId = oneGapRule[i].gappy_index;			
						} else {
							convertedId = 2*global + oneGapRule[i].gappy_index;
						}
						fsample_arr[convertedId]++;
					}
		
					for(unsigned int i = 0; i < count_one_gap; i++){
						//cerr << "Count "<<i<<endl;
						//cerr<<"Yoshi"<<out_res[i].blocknumber<<endl;
						flag = false;
						strcpy(combine, "");		
						directAcc = NULL;
						//Collect Source String
						if(i==0 || oneGapRule[i].gappy_index != oneGapRule[i-1].gappy_index
							|| i == seperatorOneGap){			
							strcpy(sourceStr, ""); 
							sourceStringIds.clear();
							if(i < seperatorOneGap){				
								convertedId = oneGapRule[i].gappy_index;
								if(oneGapRule[i].gappy_index<global){
									//Xab - need shift by one and -1
									strcat(sourceStr, "[X,1] ");		
									strcat(sourceStr, sourceName[convertedId]); 
									source_starter = tmp_blocks[convertedId].string_start;
									source_ender = source_starter + tmp_blocks[convertedId].matchlen -1;								
								} else {
									//abX - need shift by one
									strcat(sourceStr, sourceName[convertedId-global]);			
									strcat(sourceStr, " [X,1]");
									source_starter = tmp_blocks[convertedId-global].string_start;
									source_ender = source_starter + tmp_blocks[convertedId-global].matchlen -1;
								}
								
								//fprintf(stderr, "i->%d\n",i);
								for(int iccc = source_starter; iccc <= source_ender; iccc++){
									sourceStringIds.push_back(source_str[iccc]);
								}
								//fprintf(stderr, "\n");
							} else {
								//Debug
								if(oneGapRule[i].gappy_index<0){
									printf("Not possible for this oneGap rule|gappy_index %d, count %d\n", 
											oneGapRule[i].gappy_index, i);
									exit(0);
								}
								///End
								convertedId = 2*global + oneGapRule[i].gappy_index;
		
								//Deal with aXb case, where ID should be positive always. No shift needed
								number = onegapPattern[oneGapSearch[oneGapRule[i].gappy_index].position].number;
								pattern = onegapPattern[oneGapSearch[oneGapRule[i].gappy_index].position].pattern;
		
								for(jj = 0; jj < number; jj++){
									if(pattern[jj]>=0){
										sourceStringIds.push_back(pattern[jj]);
										directAcc = source_vocabulary[pattern[jj]]; 
										if(directAcc == NULL){
											fprintf(stderr, "Please.. i - %d\n", i);
										}								
									} else {
										directAcc = "[X,1]";
									}
		
									if (jj == 0){
										sprintf(str1,"%s",directAcc);
									} else {
										sprintf(str1," %s",directAcc);
									}			
									strcat(sourceStr, str1);				
								}
							}
		
							//Clean the hash table
							deDuplication.clear();					
						}
		
						unsigned int targetStart = oneGapRule[i].ref_str_start;
						unsigned int targetEnd = oneGapRule[i].ref_str_start + oneGapRule[i].end;
						unsigned int gap1Start = oneGapRule[i].ref_str_start + oneGapRule[i].gap1;
						unsigned int gap1End = oneGapRule[i].ref_str_start + oneGapRule[i].gap1_1;
		
						///get Target String aXb or abX or Xab;
						strcpy(tar_second_str, " ||| ");				
						directAcc = NULL;
						//fprintf(stderr, "target start i %d - jj %d - b %d\n", i, jj,b);
						//targetStringIds.clear();
		
						for(jj = targetStart; jj <= targetEnd; jj++){
							if (jj >= gap1Start && jj <= gap1End){
								directAcc = "[X,1]";				
							} else {
								directAcc = target_vocabulary[target_str[jj]];			
								//targetStringIds.push_back(target_str[jj]);
		
								//Debug purpose
								if(directAcc == NULL){
									fprintf(stderr, "NULL direct access\n");
									exit(0);
								}							
							}			
		
							if (jj == targetStart){
								sprintf(str1,"%s",directAcc);
							} else {
								sprintf(str1," %s",directAcc);
							}			
							strcat(tar_second_str, str1);
		
							if (jj >= gap1Start && jj <= gap1End){
								jj = gap1End;
							}
						}
		
						if(flag){
							continue;
						}
		
						///Collect Phrase Count of c(e,f)
						//cerr<<"Yoshi2"<<endl;//gogo
						strcat(combine, sourceStr);
						strcat(combine, tar_second_str);
						std::string targetStr(tar_second_str);
						iteratorDeDup = deDuplication.find(targetStr);
		
						if(iteratorDeDup == deDuplication.end()){
							int* counter = (int*)malloc(sizeof(int));
							(*counter) = 1;
		
							char* lexString = (char*)malloc(sizeof(char)*(strlen(combine)+1)); 
							memcpy((void*)lexString, combine, strlen(combine)+1);
							//Insert to map 		
							deDuplication.insert(pair<string,int*>(targetStr, counter));
		
							//for f in fphrase.words:
							//			for e in ewords:
							//			score = self.ttable.get_score(f, e, 1)
							//		if score > max_score:
							//		max_score = score
							//mlfe += -log10(max_score) if max_score > 0 else MAXSCORE
		
							///////////Offload this task to GPU
							maxEF[*lexicalTaskCounter].fastSpeedId = *lexicon_index;
							maxEF[*lexicalTaskCounter].sourcePatternCounter = sourceStringIds.size();							
							for(j =0; j < sourceStringIds.size(); j++){///	
								maxEF[*lexicalTaskCounter].sourcePattern[j] = sourceStringIds[j];
							}
							maxEF[*lexicalTaskCounter].targetStart = targetStart;
							maxEF[*lexicalTaskCounter].end = oneGapRule[i].end;
							maxEF[*lexicalTaskCounter].gap1 = oneGapRule[i].gap1;
							maxEF[*lexicalTaskCounter].gap1_1 = oneGapRule[i].gap1_1;
							(*lexicalTaskCounter)++;
							
							////Speed UP!
							fast_speed[*lexicon_index].blocknumber = convertedId;			 
							if(i < seperatorOneGap){//Xab, abX
								int realId = convertedId;
								if(convertedId>=global){
									realId = convertedId - global;
								} else if (convertedId >= 2*global){
									//Debug purpose
									printf("Not correct!! Converted Id in gappy lexicon\n");
									exit(0);
									//Debug end
								}
								
								fast_speed[*lexicon_index].all_suffix_fsample = 1 + 
									tmp_blocks[realId].end - 
									tmp_blocks[realId].start;						
							} else {//aXb							
								fast_speed[*lexicon_index].all_suffix_fsample = 1 + 
									oneGapSearch[oneGapRule[i].gappy_index].end_on_salist - 
									oneGapSearch[oneGapRule[i].gappy_index].start_on_salist;

								///Precomputation values
								if(fast_speed[*lexicon_index].all_suffix_fsample==1
									&& oneGapSA[oneGapSearch[oneGapRule[i].gappy_index].start_on_salist].length==0){
									unsigned int precompIndex = 
										oneGapSA[oneGapSearch[oneGapRule[i].gappy_index].start_on_salist].str_position;
									fast_speed[*lexicon_index].all_suffix_fsample = 1 
										- refSource->precomp_index[precompIndex].start 
										+ refSource->precomp_index[precompIndex].end
										+ refSource->featureMissingCount[precompIndex];
								}
							}
							if(isSample&&fast_speed[*lexicon_index].all_suffix_fsample > SAMPLER){
								fast_speed[*lexicon_index].all_suffix_fsample = SAMPLER;
							}
							fast_speed[*lexicon_index].f = fsample_arr[convertedId];
							fast_speed[*lexicon_index].lexical = lexString;  
							fast_speed[*lexicon_index].all_suffix_fsample_score = log10(1 + (fast_speed[*lexicon_index].all_suffix_fsample));
							fast_speed[*lexicon_index].paircount = counter;
							//fast_speed[*lexicon_index].MaxLexEgivenF = egivenf;
							//fast_speed[*lexicon_index].MaxLexFgivenE = fgivene;
		
							(*lexicon_index)++;
						} else {
							(*(iteratorDeDup->second))++;
						}
					}
					//fprintf(stderr, "OK\n");
					for(int i = 0; i< *lexicon_index; i++){
						fast_speed[i].aa = 
							-log10((float)(*fast_speed[i].paircount)/(float)fast_speed[i].all_suffix_fsample);
						fast_speed[i].bb = log10(1 + (*fast_speed[i].paircount));
					}
					timer_stop(&t);
					fprintf(stderr, "Create One Gap Lexicon time: %f\n", timer_elapsed(&t)/1000);
		
					delete[] fsample_arr;	
					return fast_speed;
				}
		
		
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
				int seperatorOneGap,//Seperator between aXb; Xab, abX
				int seperatorTwoGap_one,//XabX;aXbXc
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
				ref_t* refSource) {
		
					fprintf(stderr, "Start two gap lexicon creation process\n");
					mytimer_t t;
					timer_start(&t);
		
					bool flag = false;
					char str1[250];
					char combine[500];
					char tar_second_str[450];
					char sourceStr[250];
					map<string,categ>::iterator it;
					(*lexicon_index) = 0;
					int source_starter = -1;
					int source_ender = -1;
		
					unsigned int position = -1;
					int jj;
					int j;
		
					map<string, int*> deDuplication;
					map<string, int*>::iterator iteratorDeDup;	
					char* directAcc = NULL;
		
					int* pattern;
					int number;
					//vector<int> targetStringIds;
					vector<int> sourceStringIds;
		
					int* fsample_arr = new int[global+2*distinctOneGapCount+distinctTwoGapCount];
					memset(fsample_arr, 0, sizeof(int)*(global+2*distinctOneGapCount+distinctTwoGapCount));    
		
					//fsample feature
					unsigned int convertedId;
					for(unsigned int i =0; i < count_two_gap; i++){ 	
						if(i < seperatorTwoGap_one){
							convertedId = twoGapRule[i].twogappyindex;	//bnum - XabX		
						} else if (i < seperatorTwoGap_two){
							convertedId = global + twoGapRule[i].twogappyindex; //twoGapID - aXbXc
						} else {
							convertedId = global+distinctTwoGapCount+twoGapRule[i].twogappyindex; //2*oneBlockId - XaXb;aXbX
						}
						fsample_arr[convertedId]++;
					}
		
					for(int i = 0; i < count_two_gap; i++){
						//cerr << "Count "<<i<<endl;
						//cerr<<"Yoshi"<<out_res[i].blocknumber<<endl;
						flag = false;
						strcpy(combine, "");		
						directAcc = NULL;
						//Collect Source String
						if(i==0 || twoGapRule[i].twogappyindex != twoGapRule[i-1].twogappyindex
							|| i == seperatorTwoGap_one || i == seperatorTwoGap_two){			
							strcpy(sourceStr, ""); 
							sourceStringIds.clear();
							if(i < seperatorTwoGap_one){				
								convertedId = twoGapRule[i].twogappyindex;
								//XabX - need shift by one and -1
								strcat(sourceStr, "[X,1] ");		
								strcat(sourceStr, sourceName[convertedId]); 									
								strcat(sourceStr, " [X,2]");
		
								source_starter = tmp_blocks[convertedId].string_start;
								source_ender = source_starter + tmp_blocks[convertedId].matchlen -1;
								for(int iccc = source_starter; iccc <= source_ender; iccc++){
									sourceStringIds.push_back(source_str[iccc]);
								}
							} else if (i < seperatorTwoGap_two) {				
								convertedId = global + twoGapRule[i].twogappyindex; //twoGapID - aXbXc
		
								//Deal with aXb case, where ID should be positive always. No shift needed
								int oneGapBlockId = twoGapSearch[twoGapRule[i].twogappyindex].blockid;
		
								number = onegapPattern[oneGapSearch[oneGapBlockId].position].number;
								pattern = onegapPattern[oneGapSearch[oneGapBlockId].position].pattern;
		
								for(jj = 0; jj < number; jj++){
									if(pattern[jj]>=0){
										sourceStringIds.push_back(pattern[jj]);
										directAcc = source_vocabulary[pattern[jj]]; 										
									} else {
										directAcc = "[X,1]";
									}
		
									if (jj == 0){
										sprintf(str1,"%s",directAcc);
									} else {
										sprintf(str1," %s",directAcc);
									}			
									strcat(sourceStr, str1);				
								}
		
								strcat(sourceStr, " [X,2]");			
		
								//Deal with the rest part				
								number = twogapPattern[twoGapSearch[twoGapRule[i].twogappyindex].position].number;
								pattern = twogapPattern[twoGapSearch[twoGapRule[i].twogappyindex].position].pattern;
		
								for(jj = 0; jj < number; jj++){
									if(pattern[jj]>=0){
										sourceStringIds.push_back(pattern[jj]);
										directAcc = source_vocabulary[pattern[jj]]; 										
									} else {
										directAcc = "[X,N]";
										printf("Not possible here in two gap source string!!\n");
										exit(0);
									}
		
									sprintf(str1," %s",directAcc);
									strcat(sourceStr, str1);				
								}
		
							} else {
								//Debug
								if(twoGapRule[i].twogappyindex<0){
									printf("Not possible for this twoGap rule|gappy_index %d, count %d\n", 
											twoGapRule[i].twogappyindex, i);
									exit(0);
								}
								///End				
								convertedId = global+distinctTwoGapCount+twoGapRule[i].twogappyindex; //2*oneBlockId - XaXb;aXbX
								int oneGapBlockId;
								bool XaXb = true;
								if(convertedId >= global+distinctTwoGapCount + distinctOneGapCount){
									oneGapBlockId = twoGapRule[i].twogappyindex - distinctOneGapCount;
									XaXb = false;
								} else {
									oneGapBlockId = twoGapRule[i].twogappyindex;				
									strcat(sourceStr, "[X,1]"); 		
								}
		
								//Deal with aXb case, where ID should be positive always. No shift needed
								number = onegapPattern[oneGapSearch[oneGapBlockId].position].number;
								pattern = onegapPattern[oneGapSearch[oneGapBlockId].position].pattern;
		
								for(jj = 0; jj < number; jj++){
									if(pattern[jj]>=0){
										sourceStringIds.push_back(pattern[jj]);
										directAcc = source_vocabulary[pattern[jj]]; 
										
									} else if (XaXb){
										directAcc = "[X,2]";
									} else {
										directAcc = "[X,1]";
									}
		
									if (jj == 0&& !XaXb){
										sprintf(str1,"%s",directAcc);
									} else {
										sprintf(str1," %s",directAcc);
									}			
									strcat(sourceStr, str1);				
								}
		
								if(!XaXb){							
									strcat(sourceStr, " [X,2]");				
								}
							}
							//Clean the hash table
							deDuplication.clear();
						}
		
						unsigned int targetStart = twoGapRule[i].ref_str_start;
						unsigned int targetEnd = twoGapRule[i].ref_str_start + twoGapRule[i].end;
						unsigned int gap1Start = twoGapRule[i].ref_str_start + twoGapRule[i].gap1;
						unsigned int gap1End = twoGapRule[i].ref_str_start + twoGapRule[i].gap1_1;
						unsigned int gap2Start = twoGapRule[i].ref_str_start + twoGapRule[i].gap2;
						unsigned int gap2End = twoGapRule[i].ref_str_start + twoGapRule[i].gap2_1;
		
						///get Target String aXb or abX or Xab;
						strcpy(tar_second_str, " ||| ");				
						directAcc = NULL;
						//fprintf(stderr, "target start i %d - jj %d - b %d\n", i, jj,b);
						//targetStringIds.clear();
		
						for(jj = targetStart; jj <= targetEnd; jj++){
							if (jj >= gap1Start && jj <= gap1End){
								directAcc = "[X,1]";				
							} else if (jj >= gap2Start && jj <= gap2End) {
								directAcc = "[X,2]";
							} else {
								directAcc = target_vocabulary[target_str[jj]];			
								//targetStringIds.push_back(target_str[jj]);
							}			
		
							if (jj == targetStart){
								sprintf(str1,"%s",directAcc);
							} else {
								sprintf(str1," %s",directAcc);
							}			
							strcat(tar_second_str, str1);
		
							if (jj >= gap1Start && jj <= gap1End){
								jj = gap1End;
							} else if (jj >= gap2Start && jj <= gap2End){
								jj = gap2End;
							}
						}
		
						if(flag){
							continue;
						}
		
						///Collect Phrase Count of c(e,f)
						strcat(combine, sourceStr);
						strcat(combine, tar_second_str);
						std::string targetStr(tar_second_str);
						iteratorDeDup = deDuplication.find(targetStr);
		
						if(iteratorDeDup == deDuplication.end()){
							int* counter = (int*)malloc(sizeof(int));
							(*counter) = 1;
		
							char* lexString = (char*)malloc(sizeof(char)*(strlen(combine)+1)); 
							memcpy((void*)lexString, combine, strlen(combine)+1);
							//Insert to map 		
							deDuplication.insert(pair<string,int*>(targetStr, counter));
		
							//for f in fphrase.words:
							//			for e in ewords:
							//			score = self.ttable.get_score(f, e, 1)
							//		if score > max_score:
							//		max_score = score
							//mlfe += -log10(max_score) if max_score > 0 else MAXSCORE
		
							///////////Offload this task to GPU
							maxEF[*lexicalTaskCounter].fastSpeedId = *lexicon_index;
							maxEF[*lexicalTaskCounter].sourcePatternCounter = sourceStringIds.size();							
							for(j =0; j < sourceStringIds.size(); j++){///	
								maxEF[*lexicalTaskCounter].sourcePattern[j] = sourceStringIds[j];
							}
							maxEF[*lexicalTaskCounter].targetStart = targetStart;
							maxEF[*lexicalTaskCounter].end = twoGapRule[i].end;
							maxEF[*lexicalTaskCounter].gap1 = twoGapRule[i].gap1;
							maxEF[*lexicalTaskCounter].gap1_1 = twoGapRule[i].gap1_1;
							maxEF[*lexicalTaskCounter].gap2 = twoGapRule[i].gap2;
							maxEF[*lexicalTaskCounter].gap2_1 = twoGapRule[i].gap2_1;
							(*lexicalTaskCounter)++;
							
							////Speed UP!
							fast_speed[*lexicon_index].blocknumber = convertedId;		  
		
							int realId;
							if(i < seperatorTwoGap_one){
								realId = twoGapRule[i].twogappyindex;	//bnum - XabX		
								fast_speed[*lexicon_index].all_suffix_fsample = 1 + 
									tmp_blocks[realId].end - 
									tmp_blocks[realId].start;
							} else if (i < seperatorTwoGap_two){
								realId = twoGapRule[i].twogappyindex; //twoGapID - aXbXc
								fast_speed[*lexicon_index].all_suffix_fsample = 1 + 
									twoGapSearch[realId].end_on_salist - 
									twoGapSearch[realId].start_on_salist;
							} else {
								//2*oneBlockId - XaXb;aXbX
								if(convertedId >= global+distinctTwoGapCount + distinctOneGapCount){
									realId = twoGapRule[i].twogappyindex - distinctOneGapCount; 				
								} else {
									realId = twoGapRule[i].twogappyindex;
								}
		
								//Debug
								if(convertedId < global+distinctTwoGapCount){
									printf("Not possile global %d | distinctTwoGapCount %d | convertedId %d | i%d | twogapindex %d| realId %d\n", 
											global, distinctTwoGapCount, convertedId, i, twoGapRule[i].twogappyindex, realId);
									exit(0);
								}
								//Debug
		
								fast_speed[*lexicon_index].all_suffix_fsample = 1 + 
									oneGapSearch[realId].end_on_salist - 
									oneGapSearch[realId].start_on_salist;	

								///Precomputation values
								if(fast_speed[*lexicon_index].all_suffix_fsample==1
									&& oneGapSA[oneGapSearch[realId].start_on_salist].length==0){
									unsigned int precompIndex = oneGapSA[oneGapSearch[realId].start_on_salist].str_position;
									fast_speed[*lexicon_index].all_suffix_fsample = 1 
										- refSource->precomp_index[precompIndex].start 
										+ refSource->precomp_index[precompIndex].end
										+ refSource->featureMissingCount[precompIndex];
								}
							}			
							if(isSample&&fast_speed[*lexicon_index].all_suffix_fsample > SAMPLER){
								fast_speed[*lexicon_index].all_suffix_fsample = SAMPLER;
							}
							
							fast_speed[*lexicon_index].f = fsample_arr[convertedId];
							fast_speed[*lexicon_index].lexical = lexString;  
							fast_speed[*lexicon_index].all_suffix_fsample_score = log10(1 + (fast_speed[*lexicon_index].all_suffix_fsample));
							fast_speed[*lexicon_index].paircount = counter;
							//fast_speed[*lexicon_index].MaxLexEgivenF = egivenf;
							//fast_speed[*lexicon_index].MaxLexFgivenE = fgivene;
		
							(*lexicon_index)++;
						} else {
							(*(iteratorDeDup->second))++;
						}
					}
					//fprintf(stderr, "OK\n");
					for(int i = 0; i< *lexicon_index; i++){
						fast_speed[i].aa = 
							-log10((float)(*fast_speed[i].paircount)/(float)fast_speed[i].all_suffix_fsample);
						fast_speed[i].bb = log10(1 + (*fast_speed[i].paircount));
					}
					timer_stop(&t);
					fprintf(stderr, "Create Two Gap Lexicon time: %f\n", timer_elapsed(&t)/1000);
					delete[] fsample_arr;
		
					return fast_speed;
				}

