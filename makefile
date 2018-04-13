######################################################
# Makefile for sparsenc
######################################################

TOP = .
SRCDIR := src
OBJDIR := src
INCLUDEDIR = include src
INC_PARMS = $(INCLUDEDIR:%=-I%)

UNAME := $(shell uname)
CC := gcc
ifeq ($(UNAME), Darwin)
	SED = gsed
	CC  = gcc-7
	#CC  = clang
	HAS_SSSE3 := $(shell sysctl -a | grep supplementalsse3)
	HAS_AVX2  := $(shell sysctl -a | grep avx2)
endif
ifeq ($(UNAME), Linux)
	SED = sed
	CC  = gcc
	HAS_NEON32  := $(shell grep -i neon /proc/cpuinfo)
	HAS_NEON64  := $(shell uname -a | grep -i aarch64)
	HAS_SSSE3 := $(shell grep -i ssse3 /proc/cpuinfo)
	HAS_AVX2  := $(shell grep -i avx2 /proc/cpuinfo)
endif

GNCENC  := $(OBJDIR)/common.o $(OBJDIR)/bipartite.o $(OBJDIR)/sncEncoder.o $(OBJDIR)/galois.o $(OBJDIR)/gaussian.o $(OBJDIR)/mt19937ar.o

CFLAGS0 = -Winline -std=c99 -lm -g -DNDEBUG $(INC_PARMS)
ifneq ($(HAS_NEON32),)
	CFLAGS1 = -DARM_NEON32 -mfloat-abi=hard -mfpu=neon -O3 -std=c99
	GNCENC  := $(OBJDIR)/common.o $(OBJDIR)/bipartite.o $(OBJDIR)/sncEncoder.o $(OBJDIR)/galois_neon.o $(OBJDIR)/gaussian.o $(OBJDIR)/mt19937ar.o
endif
ifneq ($(HAS_NEON64),)
	CFLAGS1 = -DARM_NEON64 -mfloat-abi-hard -mfpu=neon -O3 -std=c99
	GNCENC  := $(OBJDIR)/common.o $(OBJDIR)/bipartite.o $(OBJDIR)/sncEncoder.o $(OBJDIR)/galois_neon.o $(OBJDIR)/gaussian.o $(OBJDIR)/mt19937ar.o
endif
ifneq ($(HAS_SSSE3),)
	CFLAGS1 = -mssse3 -DINTEL_SSSE3
endif
ifneq ($(HAS_AVX2),)
	CFLAGS1 += -mavx2 -DINTEL_AVX2
endif
# Additional compile options
# CFLAGS2 = 

vpath %.h src include
vpath %.c src examples

DEFS    := sparsenc.h common.h galois.h decoderGG.h decoderOA.h decoderBD.h decoderCBD.h decoderPP.h
RECODER := $(OBJDIR)/sncRecoder.o $(OBJDIR)/sncRecoderBATS.o 
DECODER := $(OBJDIR)/sncDecoder.o
GGDEC   := $(OBJDIR)/decoderGG.o 
OADEC   := $(OBJDIR)/decoderOA.o $(OBJDIR)/pivoting.o
BDDEC   := $(OBJDIR)/decoderBD.o $(OBJDIR)/pivoting.o
CBDDEC  := $(OBJDIR)/decoderCBD.o
PPDEC   := $(OBJDIR)/decoderPP.o

.PHONY: all
all: sncDecoder sncDecoderFile sncRecoder2Hop sncRestore

libsparsenc.so: $(GNCENC) $(GGDEC) $(OADEC) $(BDDEC) $(CBDDEC) $(PPDEC) $(RECODER) $(DECODER)
	$(CC) -shared -o libsparsenc.so $^
	
sncRLNC: $(GNCENC) $(GGDEC) $(OADEC) $(BDDEC) $(CBDDEC) $(PPDEC) $(RECODER) $(DECODER) test.RLNC.c
	$(CC) -o $@ $(CFLAGS0) $(CFLAGS1) $^

#Test snc decoder
sncDecoders: libsparsenc.so test.decoders.c
	$(CC) -L. -lsparsenc -o $@ $(CFLAGS0) $(CFLAGS1) $^
#Test snc decoder linked statically
sncDecoderST: $(GNCENC) $(GGDEC) $(OADEC) $(BDDEC) $(CBDDEC) $(PPDEC) $(RECODER) $(DECODER) test.decoders.c
	$(CC) -o $@ $(CFLAGS0) $(CFLAGS1) $^
#Test snc store/restore decoder
sncRestore: libsparsenc.so test.restore.c
	$(CC) -L. -lsparsenc -o $@ $(CFLAGS0) $(CFLAGS1) $^
#Test decoder for files
sncDecodersFile: libsparsenc.so test.file.decoders.c
	$(CC) -L. -lsparsenc -o $@ $(CFLAGS0) $(CFLAGS1) $^
#Test recoder
sncRecoder2Hop: libsparsenc.so test.2hopRecoder.c
	$(CC) -L. -lsparsenc -o $@ $(CFLAGS0) $(CFLAGS1) $^
#Test recoder
sncRecoder-n-Hop: libsparsenc.so test.nhopRecoder.c
	$(CC) -L. -lsparsenc -o $@ $(CFLAGS0) $(CFLAGS1) $^
#Test recoder, statically linked
sncRecoder-n-Hop-ST: $(GNCENC) $(GGDEC) $(OADEC) $(BDDEC) $(CBDDEC) $(PPDEC) $(RECODER) $(DECODER) test.nhopRecoder.c
	$(CC) -o $@ $(CFLAGS0) $(CFLAGS1) $^
#Test recoder
sncRecoderFly: libsparsenc.so test.butterfly.c
	$(CC) -L. -lsparsenc -o $@ $(CFLAGS0) $(CFLAGS1) $^

sncHAPmulticast: libsparsenc.so test.HAPmulticast.c
	$(CC) -L. -lsparsenc -o $@ $(CFLAGS0) $(CFLAGS1) $^

sncD2Dmulticast: libsparsenc.so test.D2Dmulticast.c
	$(CC) -L. -lsparsenc -o $@ $(CFLAGS0) $(CFLAGS1) $^

sncRecoderNhopBATS: libsparsenc.so test.batsRecoder.c
	$(CC) -L. -lsparsenc -o $@ $(CFLAGS0) $(CFLAGS1) $^

$(OBJDIR)/%.o: $(OBJDIR)/%.c $(DEFS)
	$(CC) -c -fpic -o $@ $< $(CFLAGS0) $(CFLAGS1) $(CFLAGS2)

.PHONY: clean
clean:
	rm -f *.o $(OBJDIR)/*.o libsparsenc.so sncDecoders sncDecoderST sncDecodersFile sncRecoder2Hop sncRecoder-n-Hop sncRecoder-n-Hop-ST sncRecoderFly sncRestore sncRLNC sncHAPmulticast sncD2Dmulticast sncRecoderNhopBATS

install: libsparsenc.so
	cp include/sparsenc.h /usr/include/
	if [[ `uname -a | grep -o x86_64` == "x86_64" ]]; then \
		cp libsparsenc.so /usr/lib64/; \
	else \
		cp libsparsenc.so /usr/lib/; \
	fi

.PHONY: uninstall
uninstall:
	rm -f /usr/include/sparsenc.h
	if [[ `uname -a | grep -o x86_64` == "X86_64" ]]; then \
		rm -f /usr/lib64/libsparsenc.so; \
	else \
		rm -f /usr/lib/libsparsenc.so; \
	fi
