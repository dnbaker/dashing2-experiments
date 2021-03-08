.PHONY=all tests clean obj update
CXX=g++
CC=gcc

MAKE?=make

GMATCH=$(findstring g++,$(CXX))
GIT_VERSION := $(shell git describe --abbrev=4 --always)

WARNINGS=-Wall -Wextra -Wno-char-subscripts \
		 -Wpointer-arith -Wwrite-strings -Wdisabled-optimization \
		 -Wformat -Wcast-align -Wno-unused-function -Wno-unused-parameter \
		 -pedantic -DUSE_PDQSORT -Wunused-variable -Wno-attributes -Wno-cast-align \
        -Wno-ignored-attributes -Wno-missing-braces
EXTRA?=
INCLUDE+=-I. -Isketch/include -Isketch -Isketch/vec/blaze
EXTRA_LD?=
DBG?=-DNDEBUG
OS:=$(shell uname)
FLAGS=

OPT_MINUS_OPENMP= -O3 -funroll-loops\
	  -pipe -fno-strict-aliasing -march=native -mpclmul $(FLAGS) $(EXTRA)
OPT=$(OPT_MINUS_OPENMP) -fopenmp
XXFLAGS=-fno-rtti
CXXFLAGS=$(OPT) $(XXFLAGS) -std=c++14 $(WARNINGS) -DBONSAI_VERSION=\"$(GIT_VERSION)\" $(INCLUDE)
CXXFLAGS_MINUS_OPENMP=$(OPT_MINUS_OPENMP) $(XXFLAGS) -std=c++1z $(WARNINGS) -Wno-cast-align -Wno-gnu-zero-variadic-macro-arguments -DBONSAI_VERSION=\"$(GIT_VERSION)\"
CCFLAGS=$(OPT) $(CFLAGS) -std=c11 $(WARNINGS) -DBONSAI_VERSION=\"$(GIT_VERSION)\"
LIB=-lz
LD=-L. $(EXTRA_LD)

ifneq (,$(findstring g++,$(CXX)))
	ifeq ($(shell uname),Darwin)
		ifeq (,$(findstring clang,$(CXX)))
			POPCNT_CXX:=clang
		else
			POPCNT_CXX:=$(CXX)
		endif
	else
		POPCNT_CXX:=$(CXX)
	endif
endif

all: errexp


obj: $(OBJS) $(DOBJS) $(ZOBJS) $(ZW_OBJS)

libzstd.a:
	+cd zstd && $(MAKE) lib && cp lib/libzstd.a .. && cd ..

clhash.o: bin/clhash.c
	$(CC) -std=c99 -c $< -o $@ -O3 -march=native -Iinclude/

%.o: %.c
	$(CC) $(CFLAGS) $(INCLUDE) -DNDEBUG -c $< -o $@ $(LIB)

%.do: %.c
	$(CC) $(CFLAGS) $(INCLUDE) -g -c $< -o $@ $(LIB)

%.zo: %.c
	$(CC) $(CFLAGS) $(INCLUDE) -c $< -o $@ $(LIB) $(ZCOMPILE_FLAGS)

%.o: %.cpp $(wildcard include/bonsai/*.h) $(wildcard hll/include/sketch/*.h)
	$(CXX) $(CXXFLAGS) $(DBG) $(INCLUDE) $(LD) -DNDEBUG -c $< -o $@ $(LIB)

%: bin/%.o clhash.o
	$(CXX) $< clhash.o -o $@ $(CXXFLAGS)

clean:
	rm *.o
