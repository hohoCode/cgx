#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "ComTypes.h"
#include "Timer.h"

#define PTR(a) ((char *) (&(a)))

#define MIN(a, b) (a > b ? b : a)

// lexic. order for pairs
#define leqPair(a1, a2, b1, b2) \
    ((a1 < b1) || ((a1 == b1) && (a2 <= b2)))

// and triples
#define leqTriple(a1, a2, a3, b1, b2, b3) \
    ((a1 < b1) || ((a1 == b1) && leqPair(a2, a3, b2, b3)))

// stably sort a[0..n-1] to b[0..n-1] with keys in 0..K from r
static void radixPassInt(int* a, int* b, 
        int* r, int n, int K) {
    // counter array
    int* c = (int *)calloc(K + 1, sizeof(int));
    int i = 0;

    // count occurences
    for (i = 0;  i < n;  i++) { 
        c[r[a[i]]]++;
    }

    // exclusive prefix sums
    int sum, t;
    for (i = 0, sum = 0;  i <= K;  i++) { 
        t = c[i];
        c[i] = sum;  
        sum += t;
    }

    // sort
    for (i = 0;  i < n;  i++) { 
        b[c[r[a[i]]]++] = a[i];
    }

    free(c);
}

// find the suffix array SA of s[0..n-1] in {1..K}^n
// require s[n]=s[n+1]=s[n+2]=0, n>=2
static void suffixArrayInt(int* s, int* SA, int n, int K) {

    int n0 = (n+2)/3, n1 = (n+1)/3, n2 = n/3, n02 = n0+n2; 
    int* s12  = (int *)malloc((n02 + 3) * sizeof(int));  s12[n02] = s12[n02+1] = s12[n02+2] = 0; 
    int* SA12 = (int *)malloc((n02 + 3) * sizeof(int)); SA12[n02] = SA12[n02+1] = SA12[n02+2] = 0;
    int* s0   = (int *)malloc(n0 * sizeof(int));
    int* SA0  = (int *)malloc(n0 * sizeof(int));

    // generate positions of mod 1 and mod  2 suffixes
    // the "+(n0-n1)" adds a dummy mod 1 suffix if n%3 == 1
    int i = 0;
    int j = 0;
    for (i=0, j=0; i < n+(n0-n1);  i++) {
        if (i%3 != 0) s12[j++] = i;
    }

    // lsb radix sort the mod 1 and mod 2 triples
    radixPassInt(s12 , SA12, s+2, n02, K);
    radixPassInt(SA12, s12 , s+1, n02, K);  
    radixPassInt(s12 , SA12, s  , n02, K);

    // find lexicographic names of triples
    int name = 0, c0 = -1, c1 = -1, c2 = -1;
    for (i = 0;  i < n02;  i++) {
        if (s[SA12[i]] != c0 || s[SA12[i]+1] != c1 || s[SA12[i]+2] != c2) { 
            name++;  c0 = s[SA12[i]];  c1 = s[SA12[i]+1];  c2 = s[SA12[i]+2];
        }
        if (SA12[i] % 3 == 1) {
            // left half
            s12[SA12[i]/3] = name; 
        } else { 
            // right half
            s12[SA12[i]/3 + n0] = name; 
        } 
    }

    // recurse if names are not yet unique
    if (name < n02) {
        suffixArrayInt(s12, SA12, n02, name);
        // store unique names in s12 using the suffix array 
        for (i = 0;  i < n02;  i++) s12[SA12[i]] = i + 1;
    } else // generate the suffix array of s12 directly
        for (i = 0;  i < n02;  i++) SA12[s12[i] - 1] = i; 

    // stably sort the mod 0 suffixes from SA12 by their first character
    for (i = 0, j = 0; i < n02; i++)  {
        if (SA12[i] < n0) s0[j++] = 3*SA12[i];
    }
    radixPassInt(s0, SA0, s, n0, K);

    // merge sorted SA0 suffixes and sorted SA12 suffixes
    int p, t, k;
    for (p=0,  t=n0-n1,  k=0;  k < n;  k++) {
#define GetI() (SA12[t] < n0 ? SA12[t] * 3 + 1 : (SA12[t] - n0) * 3 + 2)
        i = GetI(); // pos of current offset 12 suffix
        j = SA0[p]; // pos of current offset 0  suffix
        if (SA12[t] < n0 ? 
                leqPair(s[i], s12[SA12[t] + n0], s[j], s12[j/3]) :
                leqTriple(s[i], s[i+1], s12[SA12[t]-n0+1], s[j],s[j+1], s12[j/3+n0])) { 
            // suffix from SA12 is smaller
            SA[k] = i; t++;
            if (t == n02) { 
                // done --- only SA0 suffixes left
                for (k++;  p < n0;  p++, k++) SA[k] = SA0[p];
            }
        } else { 
            SA[k] = j;  p++; 
            if (p == n0) {
                // done --- only SA12 suffixes left
                for (k++;  t < n02;  t++, k++) SA[k] = GetI();
            }
        }
    }

    free(s12); 
    free(SA12); 
    free(SA0); 
    free(s0);
}

int recursion_lcp(int L, int R, int* lcp, int *lcpright, int * lcpleft){
    //fprintf(stdout, "Inside Recursion L %d - R %d\n", L, R);    
    if(L == R -1){
        return lcp[R];
    } 
    int M = (L+R)/2;

    lcpleft[M] = recursion_lcp(L,M,lcp,lcpright,lcpleft);
    lcpright[M] = recursion_lcp(M,R,lcp,lcpright,lcpleft);
    return MIN(lcpleft[M],lcpright[M]);
}

static void buildLCPTable(ref_t *ref) {

    ref->rank = (int *)(ref->buf + ref->toklen*3);
    assert((long)ref->rank % sizeof(int) == 0);

    ref->lcp = (int *)(ref->buf + ref->toklen*2);
    assert((long)ref->lcp % sizeof(int) == 0);

    unsigned i = 0;
    for (i = 0; i < ref->toklen; i++) {
        ref->rank[ref->sa[i]] = i;
    }

    fprintf(stderr,"Start Building LCP and Rank Array\n");
    int h = 0;
    int k = 0;
    for (i = 0; i < ref->toklen; i++) {
        if (ref->rank[i] > 0) {
            k = ref->sa[ref->rank[i]-1];
            while (ref->str[ref->sa[ref->rank[i]]+h] == ref->str[k+h]) {
                h++;
            }
            ref->lcp[ref->rank[i]] = h;
            h = 0;
        }
    }

    ref->lcpleft = (int *)(ref->buf );
    ref->lcpright = (int *)(ref->buf + ref->toklen*1);

    fprintf(stderr,"Start Buliding NonStandard LCP\n");
    ref->lcpleft[0] = 0;
    ref->lcpleft[ref->toklen -1] = 0;
    ref->lcpright[0] = 0;
    ref->lcpright[ref->toklen -1] = 0;

    recursion_lcp(0, ref->toklen -1, ref->lcp, ref->lcpright, ref->lcpleft);
    /*   --- DEBUG
         for(i = 0 ; i< ref->toklen; i++){
         if(ref->lcp[i]!=ref->lcpleft[i+1] && i!=ref->toklen-1){
         fprintf(stderr, "LCP - %d - %d\n",i, ref->lcp[i]);
         fprintf(stderr, "LCPLEFT - %d - %d\n",i, ref->lcpleft[i+1]);
         } else if (ref->lcp[i]!=ref->lcpright[i-1] && i!=0){
         fprintf(stderr, "LCP - %d - %d\n",i, ref->lcp[i]);
         fprintf(stderr, "LCPRIGHT - %d - %d\n",i, ref->lcpright[i-1]);
         }
         }
         printf("Number of Tokens %d - Testing Finished\n", ref->toklen);
         exit(0);
         */
}


void suffixArrayConstruct(ref_t *ref, int last, int* ref_temp) {

    //Buf and SA is now Sperated!!!
    //ref->sa = (int*)ref->buf + ref->toklen; 
    if(ref->toklen <= 0 || last <= 0) { 
        fprintf(stderr, "Wrong size limits exceeded. %d\n");    
        exit(0);
    }

    clock_t start1 = clock();
    
    using namespace std;
    bool outtofile = false; //true
    bool write = false; //false
    
    if (!outtofile){
        suffixArrayInt(ref_temp, ref->sa, ref->toklen, last);
    } else {
        ifstream infile("sa_precomp.txt");        
        for(unsigned int i = 0; i< ref->toklen; i++ ){
            infile >> ref->sa[i];
        }
        //infile.read((char*)&(ref->sa), ref->toklen);
        infile.close();
    }

    if (write){
        ofstream outfile("sa_precomp.txt");
        copy(ref->sa, ref->sa+ref->toklen, ostream_iterator<int>(outfile, " "));
        /*for(unsigned int i = 0; i< ref->toklen; i++ ){
          outfile << ref->sa[i];
          }*/
        //outfile.write((char*) &(ref->sa), ref->toklen);
        outfile.close();
    } 
	

    ///Limit the size
    int cc = 0;
    fprintf(stderr, "%d - SA%d\n", cc, ref->sa[cc]);

    clock_t finish = clock();  
    fprintf(stderr, "SA Construction %.4f sec\n", (double)(finish - start1) / (double)CLOCKS_PER_SEC);

    /* Build the LCP array, will be placed right after the Suffix array */
    buildLCPTable(ref);
}

int suffixArrayGetEquivalentMaxRefLen(int bufsize,	int fingerlen) {

    assert(bufsize > 0);

    /* space for the reference string itself, the suffix array, the LCP and 
       rank arrays. The arrays are integers, hence the "sizeof(int)" */
    int maxreflen = bufsize / ((5 * sizeof(int)) + 1);

    /* make sure it ends at sizeof(int) bytes boundary so that the arrays 
       (which will be placed after the string) start at integer boundary */
    maxreflen = maxreflen % sizeof(int) == 0? maxreflen: ((maxreflen/sizeof(int)) - 1) * sizeof(int);

    /* one for the null and two for suffix array construction */
    maxreflen -= 3;

    return maxreflen;
}

int suffixArrayGetRequiredRefBufferSize(int maxreflen, 	int fingerlen) {

    /* one for the null and two for suffix array construction */
    maxreflen += 3;

    /* make sure it ends at sizeof(int) bytes boundary so that the arrays 
       (which will be placed after the string) start at integer boundary */
    maxreflen = maxreflen % sizeof(int) == 0? maxreflen: ((maxreflen/sizeof(int)) + 1)*sizeof(int);

    int bufsize = maxreflen +   /* space for the reference string itself */
        (maxreflen - 3) * sizeof(int) + /* space for the suffix array */
        (maxreflen - 3) * sizeof(int) + /* space for the lcp array */
        (maxreflen - 3) * sizeof(int) + /* space for LCPCR array*/
        (maxreflen - 3) * sizeof(int) +
        (maxreflen - 3) * sizeof(int);  /* space for the rank array */


    return bufsize;
}

