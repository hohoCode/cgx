OBJDIR = obj
BINDIR = bin
CC     = g++ # -g # -std=c++0x#--std=gnu++0x
NVCC =nvcc -arch=compute_35 -code=sm_35
CUDA_INSTALL_PATH= /opt/common/cuda/cuda-5.5.22
OPT = -O3 -I./uthash/ -I/scratch0/huah/thrust/#OMPOPT = -fopenmp -lgomp
NVCCFLAGS = $(OPT) -use_fast_math -I. -I$(CUDA_INSTALL_PATH)/include # -G -g #-I./uthash/ -I /scratch0/huah/thrust/ -I$(CUDA_INSTALL_PATH)/include
CFLAGS = $(OPT) $(OMPOPT) -Wall -Wno-format -msse4.2
LFLAGS = $(OPT) $(OMPOPT) -fPIC -lm  -L$(CUDA_INSTALL_PATH)/lib64 -lcudart -lcuda
DEPS   =  $(wildcard *.h)
SOURCES = $(wildcard *.c)
CUSOURCES = $(wildcard *.cu)
OBJS    = $(SOURCES:%.c=$(OBJDIR)/%.o)
CUOBJS  = $(CUSOURCES:%.cu=$(OBJDIR)/%.cu_o)

#-------------------------------------------------------------------------------
# etags command, used to generate tags file
#-------------------------------------------------------------------------------
ETAGSCMD = rm -f TAGS; find . -name '*.c' -o -name '*.h' | xargs etags
TARGET = $(BINDIR)/strmatchcuda

.PHONY: $(TARGET) clean clean-all


all : $(OBJDIR) $(BINDIR) $(TARGET)


$(TARGET): $(OBJS) $(CUOBJS)
	@echo
	@echo Linking ...
	$(CC) $(LFLAGS) -o $@ $(OBJS) $(CUOBJS)
	@rm -f *.linkinfo
	@$(ETAGSCMD)


$(OBJS): $(OBJDIR)/%.o: %.c $(DEPS)
	$(CC) $(CFLAGS) -c $< -o $@

$(CUOBJS): $(OBJDIR)/%.cu_o: %.cu $(DEPS)
	$(NVCC) $(NVCCFLAGS) -c $< -o $@

$(OBJDIR):
	@if test ! -d $(OBJDIR); then mkdir $(OBJDIR); fi

$(BINDIR):
	@if test ! -d $(BINDIR); then mkdir $(BINDIR); fi

objects: $(OBJDIR) $(OBJS) $(CUOBJS)

clean:
	@rm -f $(OBJDIR)/*.o $(OBJDIR)/*.cu_o *~ *.linkinfo $(TARGET)

clean-all:
	@rm -rf $(OBJDIR) $(BINDIR) *~ *.linkinfo TAGS

tags:
	@$(ETAGSCMD)
