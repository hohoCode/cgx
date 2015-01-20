# cgx
UltraFast GPU Grammar Extractor for Statistical Machine Translation

INSTALLATION Instruction
------------------------------------------------------------------
1. please download thrust GPU library, which actually just needs git clone and that is it. Address is here: https://github.com/thrust/thrust

2. Install CUDA library on your server (soemtimes this step people just stuck here I see many of issues in this step for example on stackoverflow).

3. Would require the GPU device to have at least 3GB GPU memory. and latest GPU (Kephler architecture)  the best. But these codes should work with Kephler/Fermi/pre-Fermi architecture GPUs.

4. In provided MAKEFILE file, we have:
```
NVCC =nvcc -arch=compute_35 -code=sm_35
CUDA_INSTALL_PATH= /opt/common/cuda/cuda-5.5.22
OPT = -O3 -I./uthash/ -I/scratch0/huah/thrust/
NVCCFLAGS = $(OPT) -use_fast_math -I. -I$(CUDA_INSTALL_PATH)/include 
```
The above four variables need to be updated. For example, the CUDA library install path ($CUDA_INSTALL_PATH) needs to be set to the corresponding path on your GPU server; $OPT needs thrust library directory's path; The computing version of your GPU device (For example, Tesla K20 is 3.5 therfore it is -arch=compute_35), etc. Please update your enviroment variables accordingly.

5. In the main directory please compile the codes with this:
```
make
```
Probably will see lots of warnings but just ignore those as long as there are no errors. If errors probably that could be enviroment settings/configuration poblems. Can let me know if need help maybe I have seen some of weirdo before... Usually CUDA library related issues if compilation errors.

6. If you can get this executable file under bin/ directory, that means compilation is good:
```
bin/strmatchcuda
```

7. Command to run (one example on FBIS parallel data):
```
./bin/strmatchcuda 
fbis.zh 
allqueries.txt.2000 
fbis.en
fbis.aligned 
lex.bin 
gpugrammar_temp
```
Please provide the following addresses as input arguments:

- First argument is the address of chinese side of parallel corpus.
- Second: query file (each line is a query)
- Third: english side of FBIS data
- Four: Alignment file
- Five: lexical bin (which is directory from cdec precomputation step on FBIS data, I just reuse its output for lexical features)
- Six: the output directory address. (this is the directory address to hold the output grammar files. Each query will output one grammar file. So all grammars will be in thie directory if the program finishes successfully). If this directory does not exist yet please create this directory first.

Note the above example just assume these files are on the same directory, while actually they can be anywhere just specify the correct address. The above files are having the same file format as before.

8. So in the end if you see output log is like this:
```
Start Printing Gappy Phrases...
```
That usually means everything has been done and it is in the last printing step ('IO step' as in the paper).
9. Once done just go check the  gpugrammar_temp directory and there should be bunch of grammar files for queries (one file for each query).
10. Troubleshooting. Common installation problems are mainly from CUDA memory allocation/access side, causing weido memory problems during running. We encourage you to use the latest GPU driver and CUDA library.

Thanks.

