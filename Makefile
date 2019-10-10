CC=g++
OPT=-O3
STD=c++11
EXTRA_FLAGS=-g -save-temps -fverbose-asm -Wall -Wextra -Werror -pedantic -fno-omit-frame-pointer
CDEPS=-I/home/ranul/Downloads/benchmark/include
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

all: chunked_interchanged_multithread_bm chunked_interchanged_bm original_bm original_bm_unmod

chunked_interchanged_multithread_bm : Benchmark.cpp ChunkedArrayInterchangeMultithreadedTest.cpp Vmath.o
	$(CC) -std=$(STD) $(OPT) $(EXTRA_FLAGS) $(CDEPS) $(VEC) $^ -o $@ $(LDEPS)

chunked_interchanged_bm : Benchmark.cpp ChunkedArrayInterchangeTest.cpp Vmath.o
	$(CC) -std=$(STD) $(OPT) $(EXTRA_FLAGS) $(CDEPS) $(VEC) $^ -o $@ $(LDEPS)

original_bm : Benchmark.cpp OriginalTest.cpp Vmath.o
	$(CC) -std=$(STD) $(OPT) $(EXTRA_FLAGS) $(CDEPS) $(VEC) $^ -o $@ $(LDEPS)

original_bm_unmod : Benchmark.cpp OriginalTest.cpp Vmath.cpp
	$(CC) -std=$(STD) $(OPT) $(EXTRA_FLAGS) $(CDEPS) $^ -o $@ -lbenchmark -lpthread

Vmath.o : Vmath.cpp
	$(CC) -std=$(STD) $(OPT) $(EXTRA_FLAGS) $(VEC) -c $^ -o $@ 

.PHONY : clean all
clean :
	rm chunked_interchanged_multithread_bm chunked_interchanged_bm chunked_bm original_bm original_bm_unmod Benchmark.o Benchmark.s Benchmark.i Benchmark.ii Benchmark.bc Benchmark.opt.yaml ChunkedArrayInterchangeTest.o ChunkedArrayInterchangeTest.s ChunkedArrayInterchangeTest.i ChunkedArrayInterchangeTest.ii ChunkedArrayInterchangeTest.bc ChunkedArrayInterchangeTest.opt.yaml ChunkedArrayInterchangeMultithreadedTest.o ChunkedArrayInterchangeMultithreadedTest.s ChunkedArrayInterchangeMultithreadedTest.i ChunkedArrayInterchangeMultithreadedTest.ii ChunkedArrayInterchangeMultithreadedTest.bc ChunkedArrayInterchangeMultithreadedTest.opt.yaml OriginalTest.o OriginalTest.s OriginalTest.i OriginalTest.ii OriginalTest.bc OriginalTest.opt.yaml Vmath.o Vmath.s Vmath.i Vmath.ii Vmath.bc Vmath.opt.yaml *.optrpt *objdump* perf.data*

