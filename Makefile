# Makefile

BINDIR := bin
LIBDIR := lib
SRCDIR := src

CCFLAGS := -pedantic -g

# Libraries

# PC
#CC := h5pcc
#MISC_CLINKFLAGS = -L/usr/local/lib -lgsl -lgslcblas -lm -L/home/diogo/fftMPI_source/fftmpi-master/src  -lfft3dmpi -lfft2dmpi
#MISC_CCOMPILEFLAGS = -I/usr/local/include -I/home/diogo/fftMPI_source/fftmpi-master/src

# Mac
CC := h5pcc
MISC_CLINKFLAGS = -L/usr/lib64 -lgsl -lgslcblas -lm -L/home/dsimoes/fftmpi-master/src  -lfft3dmpi -lfft2dmpi
MISC_CCOMPILEFLAGS = -I/usr/include -I/home/dsimoes/fftmpi-master/src

# Accelerates
#CC := /data/dsimoes/HDF5_install/bin/h5pcc
#MISC_CLINKFLAGS = -L/data/dsimoes/gsl_install/lib -lgsl -lgslcblas -lm -L/data/dsimoes/fftMPI_source/fftmpi-master/src -lfft3dmpi -lfft2dmpi
#MISC_CCOMPILEFLAGS = -I/data/dsimoes/gsl_install/include -I/data/dsimoes/fftMPI_source/fftmpi-master/src

VPATH = main:src

SRC := $(wildcard $(SRCDIR)/*.c)
OBJ := $(patsubst %.c, $(BINDIR)/%.o, $(notdir $(SRC)))
INC := $(wildcard $(SRCDIR)/*.h)

lib: $(LIBDIR)/libquetzal.a

$(LIBDIR)/libquetzal.a: $(OBJ) 
	@echo make lib...
	ar ruv $@ $^
	ranlib $@

%.exe: $(BINDIR)/%.o $(LIBDIR)/libquetzal.a
	@echo LINKING!...
	$(CC) $(CCFLAGS) -I $(SRCDIR) $(MISC_CCOMPILEFLAGS) $< -o $(BINDIR)/$@ -L lib -lquetzal $(MISC_CLINKFLAGS)

$(BINDIR)/%.o: %.c | $(INC)
	@echo COMPILING!... $< $(MPI_ROOT) $(H5_ROOT)
	$(CC) $(CCFLAGS) -I $(SRCDIR) $(MISC_CCOMPILEFLAGS) -c $< -o $@ $(MISC_CLINKFLAGS)

######### clean

tilde := $(wildcard */*~) $(wildcard *~)
exe := $(wildcard */*.exe) $(wildcard *.exe)
obj := $(wildcard */*.o) $(wildcard *.o) $(wildcard */*.so) $(wildcard */*.pcm) $(wildcard */*.d)

clean:
	@echo cleaning dir...
	rm -f $(exe) $(obj) $(tilde) $(LIBDIR)/*
