# LongUMI Makefile
CC ?= cc
UNAME_S := $(shell uname -s 2>/dev/null)
ifeq ($(UNAME_S),Darwin)
OPENMP_CFLAGS ?=
OPENMP_LDFLAGS ?=
else
OPENMP_CFLAGS ?= -fopenmp
OPENMP_LDFLAGS ?= -fopenmp
endif

CFLAGS ?= -O3 -march=native -Wall -Wextra -std=c11 $(OPENMP_CFLAGS)
LDFLAGS ?= $(OPENMP_LDFLAGS)
EXTRA_CFLAGS ?=
EXTRA_LDFLAGS ?=
COMMON_LIB_DIRS := -L/usr/local/lib -L/opt/homebrew/lib

PKGCONFIG ?= $(shell command -v pkg-config 2>/dev/null)
ifeq ($(origin HTSLIB_INC), undefined)
ifeq ($(PKGCONFIG),)
RHTSLIB_PREFIX := /Library/Frameworks/R.framework/Versions/4.2-arm64/Resources/library/Rhtslib
ifneq ($(wildcard $(RHTSLIB_PREFIX)/include/htslib/sam.h),)
HTSLIB_INC ?= $(RHTSLIB_PREFIX)/include
HTSLIB_LIB ?= $(RHTSLIB_PREFIX)/usrlib
else
HTSLIB_INC ?= /usr/local/include
HTSLIB_LIB ?= /usr/local/lib
endif
HTS_CFLAGS := -I$(HTSLIB_INC)
HTS_LDFLAGS := -L$(HTSLIB_LIB) $(COMMON_LIB_DIRS) -lhts -lz -lbz2 -llzma -lcurl -lcrypto
else
HTS_CFLAGS := $(shell pkg-config --cflags htslib)
HTS_LDFLAGS := $(shell pkg-config --libs htslib)
endif
else
HTS_CFLAGS := -I$(HTSLIB_INC)
HTS_LDFLAGS := -L$(HTSLIB_LIB) $(COMMON_LIB_DIRS) -lhts -lz -lbz2 -llzma -lcurl -lcrypto
endif

INC := $(HTS_CFLAGS)
LIBS := $(HTS_LDFLAGS)

OBJS = main.o cli.o io.o sj.o key.o umi.o vector.o key_build.o pipeline.o
BIN = longumi

all: $(BIN)

%.o: %.c
	$(CC) $(CFLAGS) $(INC) $(EXTRA_CFLAGS) -c $< -o $@

$(BIN): $(OBJS)
	$(CC) $(CFLAGS) $(INC) -o $@ $(OBJS) $(LDFLAGS) $(LIBS) $(EXTRA_LDFLAGS)

clean:
	rm -f $(OBJS) $(BIN)
