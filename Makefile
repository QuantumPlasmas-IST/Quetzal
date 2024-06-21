# Makefile

BINDIR := bin
LIBDIR := lib
SRCDIR := src

CCFLAGS := -pedantic -g

CC := gcc

# Libraries

# PC
MPI_ROOT = /home/diogo/MPI_install
MPI_CCOMPILEFLAGS = -I$(MPI_ROOT)/include -I$(MPI_ROOT)/openmpi
MPI_CLINKFLAGS    = -L$(MPI_ROOT)/lib -lmpi

H5_ROOT = /home/diogo/hdf5_install
H5_CCOMPILEFLAGS = -I$(H5_ROOT)/include
H5_CLINKFLAGS    = -Wl,-rpath -Wl,$(H5_ROOT)/lib -L$(H5_ROOT)/lib -lhdf5 -lm
HDF5_ENABLE_PARALLEL = ON

# Accelerates
#MPI_CCOMPILEFLAGS = -I$(MPI_INCLUDE)
#MPI_CLINKFLAGS    = -Wl,-rpath -Wl,$(MPI_LIB) -L$(MPI_LIB) \
                    -lmpifort -lmpi -Wl,--enable-new-dtags

#H5_ROOT = /data/dsimoes/hdf5_install
#H5_CCOMPILEFLAGS = -I$(H5_ROOT)/include 
#H5_CLINKFLAGS    = -Wl,-rpath -Wl,$(H5_ROOT)/lib -L$(H5_ROOT)/lib -lhdf5 -lm
#HDF5_ENABLE_PARALLEL = ON

VPATH = main:src

SRC := $(wildcard $(SRCDIR)/*.c)
OBJ := $(patsubst %.c, $(BINDIR)/%.o, $(notdir $(SRC)))
INC := $(wildcard $(SRCDIR)/*.h)

lib: $(LIBDIR)/libFC.a

$(LIBDIR)/libFC.a: $(OBJ) 
	@echo make lib...
	ar ruv $@ $^
	ranlib $@
    
%.exe: $(BINDIR)/%.o $(LIBDIR)/libFC.a
	@echo LINKING!...
	$(CC) $(CCFLAGS) -I $(SRCDIR) $(MPI_CCOMPILEFLAGS) $(H5_CCOMPILEFLAGS) $< -o $(BINDIR)/$@ -L lib -lFC $(MPI_CLINKFLAGS) $(H5_CLINKFLAGS)

$(BINDIR)/%.o: %.c | $(INC)
	@echo COMPILING!... $< $(MPI_ROOT) $(H5_ROOT)
	$(CC) $(CCFLAGS) -I $(SRCDIR) $(MPI_CCOMPILEFLAGS) $(H5_CCOMPILEFLAGS) -c $< -o $@ $(MPI_CLINKFLAGS) $(H5_CLINKFLAGS)

######### clean

tilde := $(wildcard */*~) $(wildcard *~)
exe := $(wildcard */*.exe) $(wildcard *.exe)
obj := $(wildcard */*.o) $(wildcard *.o) $(wildcard */*.so) $(wildcard */*.pcm) $(wildcard */*.d)

clean:
	@echo cleaning dir...
	rm -f $(exe) $(obj) $(tilde) $(LIBDIR)/*