# cgx

Fast CUDA Grammar eXtractor for Statistical Machine Translation

Introduction
------------

This is a GPU grammar extractor for statistical machine translation. This tool can extract hierarchical translation grammars on GPU efficiently.

For more details, please refer to the following two papers:
- ``Approximate Pattern Matching on GPUs for On-Demand Extraction of Hierarchical Translation Grammars``.
Hua He, Jimmy Lin, and Adam Lopez. Transactions of the Association for Computational Linguistics (TACL). 2015 January. In press.
- ``Massively Parallel Suffix Array Queries and On-demand Phrase Extraction for Statistical Machine Translation using
GPUs``. Hua He, Jimmy Lin, and Adam Lopez. North American Chapter of the Association for Computational Linguistics (NAACL). 2013.

We will further update our codes significantly in our next release to make this tool easier to use... Stay tuned.

You are very welcome to share your usage experiences with us. Thank you.


Installation
------------

- Install CUDA library and CUDA driver on your machine. Please follow the instruction on Nvidia website:
https://developer.nvidia.com/cuda-zone

- This program requires the GPU device to have at least 4GB GPU memory. The codes work with Kephler/Fermi/pre-Fermi architecture GPUs.

- Please download thrust GPU library, simply git clone its repo from here: https://github.com/thrust/thrust . It should normally work well, just in case to avoid any future version incompatibility since thrust library is constantly updated, you can also switch to its git branch: 1.8.0 version.

- In provided `Makefile`, we have:
```
NVCC =nvcc -arch=compute_35 -code=sm_35
CUDA_INSTALL_PATH= /opt/common/cuda/cuda-5.5.22
OPT = -O3 -I./uthash/ -I/scratch0/huah/thrust/
```

- The above three variables in `MAKEFILE` need to be updated according to your runnning enviroment. For example, the CUDA library install path `CUDA_INSTALL_PATH` needs to be set to the corresponding path on your GPU machine; `OPT` needs thrust library directory's path; so is the computing version of your GPUs (For example, Tesla K20 is 3.5 therfore it is `-arch=compute_35 -code=sm_35`); etc. 

- In the main directory please compile the codes as below. Probably you will see lots of warnings please just ignore those as long as there are no errors. If errors that could probably be CUDA library related issues, please update your CUDA driver and library to the latest version.
```
make
```

- If you can get this executable file under bin/ directory, that means compilation is good:
```
bin/strmatchcuda
```

- Though our Cgx grammar extractor is technically independent from the popular SMT system `cdec`, currently you may still need its SMT decoder to use our generated hierarchical grammars and also its intermidiate output (lexical file as shown below), therefore please install `cdec` following its instructions and tutorials here:
http://www.cdec-decoder.org/guide/


Running
------------

- Command to run (one toy example is provided: hansards French-English parallel data under `./toy` directory):
```
./bin/strmatchcuda 
./toy/hansards.f
./toy/query.f
./toy/hansards.e
./toy/hansards.a
./toy/lex.bin
gpugrammar_temp
```
- Please provide the following as input arguments:
    1. First argument is the address of source side of parallel corpus (`toy/hansards.f` file).
    2. Second: query file (each line is a query sentence, `toy/query.f` file)
    3. Third: target side of your parallel data (`toy/hansards.e` file)
    4. Four: alignment file (same format as used in cdec)
    5. Five: lexical file for lexical feature calculation. It is directly from cdec's precomputation step and is in text format (not binary).
    6. Six: the output directory address. (this is the directory address to hold the output grammar files. Each query will output one grammar file. So all grammars will be in the directory). If this directory does not exist yet please create this directory first.
    7. All input files are having the same format as cdec's grammmar extractor. Please check the provided toy data or http://www.cdec-decoder.org/guide/tutorial.html for details.

- In the end if you see output log is like the below, this means everything has been done and it is in the last printing step (`IO step` as in the paper):
```
Start Printing Gappy Phrases...
```

- Once done just go check the gpugrammar_temp directory and there should be bunch of grammar files for queries (one file for each query).

- Short description of running process: first of all read in parallel corpus and queries, then construct suffix array and auxiliary data structures for indexing, then do GPU-based gappy phrase lookup on suffix array and hierarchical grammar extraction, finally calculate features and output. Please refer to papers for more details on each pass. 

- Currently the running of this program includes one-time costs (e.g. suffix array construction/precomputation) and real hierarchical grammar extraction costs. Therefore depending on your parallel corpus typically this can take some very acceptable extra time before real GPU grammar extraction starts. This program will get updates and such one-time costs will be separated in future release.

- Troubleshooting. Common installation problems are mainly from CUDA memory allocation/access side, causing weido memory problems during running. We encourage you to use the latest GPU driver and CUDA library.

Thanks.

