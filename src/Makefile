# IsoUMI Makefile
CC ?= cc
UNAME_S := $(shell uname -s 2>/dev/null)
ifeq ($(UNAME_S),Darwin)
OPENMP_CFLAGS ?=
OPENMP_LDFLAGS ?=
else
OPENMP_CFLAGS ?= -fopenmp
OPENMP_LDFLAGS ?= -fopenmp
endif

NATIVE ?= 0
ifeq ($(NATIVE),1)
NATIVE_CFLAGS ?= -march=native
else
NATIVE_CFLAGS ?=
endif

CFLAGS ?= -O3 $(NATIVE_CFLAGS) -Wall -Wextra -std=c11 $(OPENMP_CFLAGS)
LDFLAGS ?= $(OPENMP_LDFLAGS)
EXTRA_CFLAGS ?=
EXTRA_LDFLAGS ?=
COMMON_LIB_DIRS := $(foreach d,/usr/local/lib /opt/homebrew/lib,$(if $(wildcard $(d)),-L$(d)))
RHTSLIB_PREFIX ?= /Library/Frameworks/R.framework/Versions/4.2-arm64/Resources/library/Rhtslib
HTS_EXTRA_LIBS ?= -lz -lbz2 -llzma -lcurl
CRYPTO_LIBS ?=

PKGCONFIG ?= $(shell command -v pkg-config 2>/dev/null)
HAVE_HTSLIB_PC := $(shell if [ -n "$(PKGCONFIG)" ] && $(PKGCONFIG) --exists htslib 2>/dev/null; then echo yes; fi)
ifeq ($(origin HTSLIB_INC), undefined)
ifeq ($(HAVE_HTSLIB_PC),yes)
HTS_CFLAGS := $(shell $(PKGCONFIG) --cflags htslib)
HTS_LDFLAGS := $(shell $(PKGCONFIG) --libs htslib)
else
ifneq ($(wildcard $(RHTSLIB_PREFIX)/include/htslib/sam.h),)
HTSLIB_INC ?= $(RHTSLIB_PREFIX)/include
HTSLIB_LIB ?= $(RHTSLIB_PREFIX)/usrlib
else
HTSLIB_INC ?= /usr/local/include
HTSLIB_LIB ?= /usr/local/lib
endif
HTS_CFLAGS := -I$(HTSLIB_INC)
HTS_LDFLAGS := -L$(HTSLIB_LIB) $(COMMON_LIB_DIRS) -lhts $(HTS_EXTRA_LIBS) $(CRYPTO_LIBS)
endif
else
HTS_CFLAGS := -I$(HTSLIB_INC)
HTS_LDFLAGS := -L$(HTSLIB_LIB) $(COMMON_LIB_DIRS) -lhts $(HTS_EXTRA_LIBS) $(CRYPTO_LIBS)
endif

INC := $(HTS_CFLAGS)
LIBS := $(HTS_LDFLAGS)

OBJS = main.o cli.o io.o sj.o key.o umi.o vector.o key_build.o pipeline.o
BIN = isoumi

all: $(BIN)

%.o: %.c
	$(CC) $(CFLAGS) $(INC) $(EXTRA_CFLAGS) -c $< -o $@

$(BIN): $(OBJS)
	$(CC) $(CFLAGS) $(INC) -o $@ $(OBJS) $(LDFLAGS) $(LIBS) $(EXTRA_LDFLAGS)

.PHONY: all test clean

test: $(BIN)
	sh ../tests/smoke.sh ./$(BIN)
	sh ../tests/synthetic_truth.sh ./$(BIN)
	sh ../tests/bam_tags.sh ./$(BIN)

clean:
	rm -f $(OBJS) $(BIN)
