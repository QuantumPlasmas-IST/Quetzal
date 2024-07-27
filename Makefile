# Makefile

BINDIR := bin
LIBDIR := lib
SRCDIR := src

CCFLAGS := -pedantic -g

CC := h5pcc

# Libraries

# Laptop-PC
#MISC_CLINKFLAGS = -L/usr/local/lib -lgsl -lgslcblas -lm
#MISC_CCOMPILEFLAGS = -I/usr/local/include

# Desktop-Mac
MISC_CLINKFLAGS = -L/usr/lib64 -lgsl -lgslcblas -lm
MISC_CCOMPILEFLAGS = -I/usr/include

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
	$(CC) $(CCFLAGS) -I $(SRCDIR) $(MISC_CCOMPILEFLAGS) $< -o $(BINDIR)/$@ -L lib -lFC $(MISC_CLINKFLAGS)

$(BINDIR)/%.o: %.c | $(INC)
	@echo COMPILING!... $<
	$(CC) $(CCFLAGS) -I $(SRCDIR) $(MISC_CCOMPILEFLAGS) -c $< -o $@ $(MISC_CLINKFLAGS)

######### clean

tilde := $(wildcard */*~) $(wildcard *~)
exe := $(wildcard */*.exe) $(wildcard *.exe)
obj := $(wildcard */*.o) $(wildcard *.o) $(wildcard */*.so) $(wildcard */*.pcm) $(wildcard */*.d)

clean:
	@echo cleaning dir...
	rm -f $(exe) $(obj) $(tilde) $(LIBDIR)/*