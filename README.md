# cgx
UltraFast GPU Grammar Extractor for Statistical Machine Translation

INSTALLATION Instruction
------------------------------------------------------------------
- please download thrust GPU library, which actually just needs git clone and that is it. Address is here: https://github.com/thrust/thrust

- Install CUDA library on your server (soemtimes this step people just stuck here I see many of issues in this step for example on stackoverflow).

- Would require the GPU device to have at least 3GB GPU memory. and latest GPU (Kephler architecture)  the best. But these codes should work with Kephler/Fermi/pre-Fermi architecture GPUs.

- In provided MAKEFILE file, we have:
```
NVCC =nvcc -arch=compute_35 -code=sm_35
CUDA_INSTALL_PATH= /opt/common/cuda/cuda-5.5.22
OPT = -O3 -I./uthash/ -I/scratch0/huah/thrust/
NVCCFLAGS = $(OPT) -use_fast_math -I. -I$(CUDA_INSTALL_PATH)/include 
```
- The above four variables in MAKEFILE need to be updated accoring to your runnning enviroment. For example, the CUDA library install path ($CUDA_INSTALL_PATH) needs to be set to the corresponding path on your GPU server; $OPT needs thrust library directory's path; The computing version of your GPU device should be set (For example, Tesla K20 is 3.5 therfore it is -arch=compute_35), etc. 

- In the main directory please compile the codes as below. Probably will see lots of warnings please just ignore those as long as there are no errors. If errors that could probably be CUDA library related issues, Please update your cuda driver and library to the latest version.
```
make
```

- If you can get this executable file under bin/ directory, that means compilation is good:
```
bin/strmatchcuda
```

- Command to run (one example on FBIS parallel data):
```
./bin/strmatchcuda 
fbis.zh 
allqueries.txt.2000 
fbis.en
fbis.aligned 
lex.bin 
gpugrammar_temp
```
- Please provide the following as input arguments:
    1. First argument is the address of chinese side of parallel corpus.
    2. Second: query file (each line is a query)
    3. Third: english side of FBIS data
    4. Four: Alignment file
    5. Five: lexical bin (which is directly from cdec precomputation step on FBIS data, just reuse its output for lexical features)
    6. Six: the output directory address. (this is the directory address to hold the output grammar files. Each query will output one grammar file. So all grammars will be in the directory). If this directory does not exist yet please create this directory first.
    7. Note the above example just assume these files are on the same directory, while actually they can be anywhere just specify the correct address. The above files are having the same file format as before.

- So in the end if you see output log is like the below, which means everything has been done and it is in the last printing step ('IO step' as in the paper):
```
Start Printing Gappy Phrases...
```
- Once done just go check the  gpugrammar_temp directory and there should be bunch of grammar files for queries (one file for each query).
- Troubleshooting. Common installation problems are mainly from CUDA memory allocation/access side, causing weido memory problems during running. We encourage you to use the latest GPU driver and CUDA library.
- We will further update our codes significantly in our next release to make this code easier to use.. Stay tuned.

Thanks.

