
#include <unistd.h>
#include "ComTypes.h"

options_t options = {

	NULL,	/* OPT_reffile */
	NULL,	/* OPT_qryfile */
	NULL,   /* OPT_reftargetfile */
	NULL,
	NULL,	/* OPT_timefile */
    NULL,
	1,		/* OPT_minmatchlen */
	10,		/* OPT_fingerlen */
	NULL,
};

extern void start(options_t *options);

void printHelp(char * exename) {
	printf("\nGPU source codes for gappy extraction. Please check your input arguments.\n"
		"\n",
		exename
		);

	exit(0);
}

void parseCommandLine(int argc, char ** argv) {

        int errflg = 0;
		int ch;
		optarg = NULL;

		while(!errflg && ((ch = getopt (argc, argv, "hl:t:s:")) != EOF)) {
			switch  (ch) {
			case 'h': printHelp(argv[0]); break;
			case '?': fprintf(stderr, "Unknown option %c\n", optopt); 
				errflg = 1; break;
			case 'l': options.minmatchlen = atoi(optarg); break;
			case 't': options.fingerlen = atoi(optarg); break;
			case 's': options.timefile = optarg; break;
			default: errflg = 1; break;
			};
		}

		if ((optind != argc-6) || errflg) { ///change from 2 to 3
			printHelp(argv[0]); 
		}

		if (options.fingerlen > 10 || options.fingerlen <= 0) {
			fprintf(stderr, "finger length must be between 1 and 10\n"); 
			exit(0);
		}

		options.reffile = argv[optind++];
		options.qryfile = argv[optind++];
		options.reftargetfile = argv[optind++];
		options.align = argv[optind++];
		options.wordscdec = argv[optind++];
		options.destinationDirectory = argv[optind++];
}

void printOptions() {
	fprintf(stderr, 
		"reference file: %s\n"
		"query file: %s\n"
		"ref target file: %s\n"
		"align file: %s\n"
		"Minimum match length: %d\n",
		options.reffile,
		options.qryfile,
		options.reftargetfile,
		options.align,
		options.minmatchlen);

}

int main(int argc, char **argv) {

		parseCommandLine(argc, argv);
		printOptions();
		start(&options);

		return 0;
}
