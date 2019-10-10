CC=g++
OPT=-O3
STD=c++11
EXTRA_FLAGS=-g -save-temps -fverbose-asm -Wall -Wextra -Werror -pedantic -fno-omit-frame-pointer
CDEPS=-I/home/ranul/Downloads/benchmark/include
GARBAGE=*.o *.s *.i *.ii *.bc *.opt.yaml *.optrpt *objdump* *perf.data*

COURTEMANCHE_BM=CRM98_chunked_interchanged_multithread_bm CRM98_chunked_interchanged_bm CRM98_original_withflags_bm CRM98_original_bm

ifeq "$(CC)" "icpc"
LDEPS=-lbenchmark -lpthread -mkl
else
LDEPS=-lbenchmark -lpthread -lmvec
endif

ifeq "$(CC)" "g++"
VEC=-march=sandybridge -ftree-vectorize -fopenmp -ffast-math
endif
ifeq "$(CC)" "clang++"
VEC=-march=sandybridge -ftree-vectorize -fopenmp -ffast-math -fsave-optimization-record -Rpass=loop-vectorize -Rpass-missed=loop-vectorize -fno-math-errno
endif
ifeq "$(CC)" "icpc"
VEC=-xHost -qopenmp -qopt-report5 -no-prec-div
endif
# VEC=-march=sandybridge -fno-tree-vectorize -fopenmp -ffast-math
# VEC=-xHost -qopenmp -qopt-report5 -no-prec-div -no-vec
# VEC=

all: $(COURTEMANCHE_BM)

CRM98_chunked_interchanged_multithread_bm : Benchmark.cpp ChunkedArrayInterchangeMultithreadedTest_CRM98.cpp Vmath.o
	$(CC) -std=$(STD) $(OPT) $(EXTRA_FLAGS) $(CDEPS) $(VEC) $^ -o $@ $(LDEPS)

CRM98_chunked_interchanged_bm : Benchmark.cpp ChunkedArrayInterchangeTest_CRM98.cpp Vmath.o
	$(CC) -std=$(STD) $(OPT) $(EXTRA_FLAGS) $(CDEPS) $(VEC) $^ -o $@ $(LDEPS)

CRM98_original_withflags_bm : Benchmark.cpp OriginalTest_CRM98.cpp Vmath.o
	$(CC) -std=$(STD) $(OPT) $(EXTRA_FLAGS) $(CDEPS) $(VEC) $^ -o $@ $(LDEPS)

CRM98_original_bm : Benchmark.cpp OriginalTest_CRM98.cpp Vmath.cpp
	$(CC) -std=$(STD) $(OPT) $(EXTRA_FLAGS) $(CDEPS) $^ -o $@ -lbenchmark -lpthread

Vmath.o : Vmath.cpp
	$(CC) -std=$(STD) $(OPT) $(EXTRA_FLAGS) $(VEC) -c $^ -o $@

.PHONY : clean all
clean :
	rm $(COURTEMANCHE_BM) $(GARBAGE)


