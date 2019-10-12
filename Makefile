CC=icpc
OPT=-O3
STD=c++11
EXTRA_FLAGS=-g -Wall -Wextra -Werror -pedantic -fno-omit-frame-pointer -save-temps
CDEPS=-I/home/ranul/Downloads/benchmark/include
GARBAGE=*.o *.s *.i *.ii *.bc *.opt.yaml *.optrpt *objdump* *perf.data*

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
VEC=-xHost -qopenmp -qopt-report5 -no-prec-div -no-prec-sqrt
endif
# VEC=-march=sandybridge -fno-tree-vectorize -fopenmp -ffast-math
# VEC=-xHost -qopenmp -qopt-report5 -no-vec -no-prec-div
# VEC=

COURTEMANCHE_BM=CRN98_chunked_interchanged_splitloops_multithread_bm \
		CRN98_chunked_interchanged_splitloops_bm \
		CRN98_chunked_interchanged_multithread_bm \
		CRN98_chunked_interchanged_bm \
		CRN98_original_withvecflags_bm \
		CRN98_original_bm
COURTEMANCHE_AP=CRN98_chunked_interchanged_splitloops_multithread_ap \
		CRN98_chunked_interchanged_splitloops_ap \
		CRN98_chunked_interchanged_multithread_ap \
		CRN98_chunked_interchanged_ap \
		CRN98_original_withvecflags_ap \
		CRN98_original_ap
FENTONKARMA_BM=FK_chunked_interchanged_multithreaded_bm \
	       FK_chunked_interchanged_bm \
	       FK_original_withvecflags_bm \
	       FK_original_bm
FENTONKARMA_AP=FK_chunked_interchanged_multithreaded_ap \
	       FK_chunked_interchanged_ap \
	       FK_original_withvecflags_ap \
	       FK_original_ap


FK_DEFS=-DNUM_IN_OUT_VARS=3 -DNUM_GATE_VARS=2
CRN98_DEFS=-DNUM_IN_OUT_VARS=21 -DNUM_GATE_VARS=15 -DHAS_CONCENTRATIONS

all: $(COURTEMANCHE_BM) $(FENTONKARMA_BM) $(COURTEMANCHE_AP) $(FENTONKARMA_AP)

courtemanche: $(COURTEMANCHE_BM) $(COURTEMANCHE_AP)

fenton: $(FENTONKARMA_BM) $(FENTONKARMA_AP)

actionpotential: $(COURTEMANCHE_AP) $(FENTONKARMA_AP)

CRN98_chunked_interchanged_multithread_bm : Benchmark.o Globals_CRN98.o ChunkedTimeIntegrate.cpp ChunkedArrayInterchangeTest_CRN98.cpp Vmath.o
	$(CC) -std=$(STD) $(CRN98_DEFS) $(OPT) $(EXTRA_FLAGS) -DMULTITHREAD $(CDEPS) $(VEC) $^ -o $@ $(LDEPS)

CRN98_chunked_interchanged_splitloops_multithread_bm : Benchmark.o Globals_CRN98.o ChunkedTimeIntegrate.cpp ChunkedArrayInterchangeSplitLoopsTest_CRN98.cpp Vmath.o
	$(CC) -std=$(STD) $(CRN98_DEFS) $(OPT) $(EXTRA_FLAGS) -DMULTITHREAD $(CDEPS) $(VEC) $^ -o $@ $(LDEPS)

CRN98_chunked_interchanged_splitloops_bm : Benchmark.o Globals_CRN98.o ChunkedTimeIntegrate.cpp ChunkedArrayInterchangeSplitLoopsTest_CRN98.cpp Vmath.o
	$(CC) -std=$(STD) $(CRN98_DEFS) $(OPT) $(EXTRA_FLAGS) $(CDEPS) $(VEC) $^ -o $@ $(LDEPS)

CRN98_chunked_interchanged_bm : Benchmark.o Globals_CRN98.o ChunkedTimeIntegrate.cpp ChunkedArrayInterchangeTest_CRN98.cpp Vmath.o
	$(CC) -std=$(STD) $(CRN98_DEFS) $(OPT) $(EXTRA_FLAGS) $(CDEPS) $(VEC) $^ -o $@ $(LDEPS)

CRN98_original_withvecflags_bm : Benchmark.o Globals_CRN98.o OriginalTest_CRN98.cpp Vmath.o
	$(CC) -std=$(STD) $(OPT) $(EXTRA_FLAGS) $(CDEPS) $(VEC) $^ -o $@ $(LDEPS)

CRN98_original_bm : Benchmark.cpp Globals_CRN98.cpp OriginalTest_CRN98.cpp Vmath.cpp
	$(CC) -std=$(STD) $(OPT) $(EXTRA_FLAGS) $(CDEPS) $^ -o $@ -lbenchmark -lpthread


CRN98_chunked_interchanged_multithread_ap : ActionPotentialTest.o Globals_CRN98.o ChunkedTimeIntegrate.cpp ChunkedArrayInterchangeTest_CRN98.cpp Vmath.o
	$(CC) -std=$(STD) $(CRN98_DEFS) $(OPT) $(EXTRA_FLAGS) -DMULTITHREAD $(CDEPS) $(VEC) $^ -o $@ $(LDEPS)

CRN98_chunked_interchanged_splitloops_multithread_ap : ActionPotentialTest.o Globals_CRN98.o ChunkedTimeIntegrate.cpp ChunkedArrayInterchangeSplitLoopsTest_CRN98.cpp Vmath.o
	$(CC) -std=$(STD) $(CRN98_DEFS) $(OPT) $(EXTRA_FLAGS) -DMULTITHREAD $(CDEPS) $(VEC) $^ -o $@ $(LDEPS)

CRN98_chunked_interchanged_splitloops_ap : ActionPotentialTest.o Globals_CRN98.o ChunkedTimeIntegrate.cpp ChunkedArrayInterchangeSplitLoopsTest_CRN98.cpp Vmath.o
	$(CC) -std=$(STD) $(CRN98_DEFS) $(OPT) $(EXTRA_FLAGS) $(CDEPS) $(VEC) $^ -o $@ $(LDEPS)

CRN98_chunked_interchanged_ap : ActionPotentialTest.o Globals_CRN98.o ChunkedTimeIntegrate.cpp ChunkedArrayInterchangeTest_CRN98.cpp Vmath.o
	$(CC) -std=$(STD) $(CRN98_DEFS) $(OPT) $(EXTRA_FLAGS) $(CDEPS) $(VEC) $^ -o $@ $(LDEPS)

CRN98_original_withvecflags_ap : ActionPotentialTest.o Globals_CRN98.o OriginalTest_CRN98.cpp Vmath.o
	$(CC) -std=$(STD) $(OPT) $(EXTRA_FLAGS) $(CDEPS) $(VEC) $^ -o $@ $(LDEPS)

CRN98_original_ap : ActionPotentialTest.cpp Globals_CRN98.cpp OriginalTest_CRN98.cpp Vmath.cpp
	$(CC) -std=$(STD) $(OPT) $(EXTRA_FLAGS) $(CDEPS) $^ -o $@ -lbenchmark -lpthread


FK_chunked_interchanged_multithreaded_bm : Benchmark.o Globals_FK.o ChunkedTimeIntegrate.cpp ChunkedArrayInterchangeTest_FK.cpp Vmath.o
	$(CC) -std=$(STD) $(FK_DEFS) $(OPT) $(EXTRA_FLAGS) -DMULTITHREAD $(CDEPS) $(VEC) $^ -o $@ $(LDEPS)

FK_chunked_interchanged_bm : Benchmark.o Globals_FK.o ChunkedTimeIntegrate.cpp ChunkedArrayInterchangeTest_FK.cpp Vmath.o
	$(CC) -std=$(STD) $(FK_DEFS) $(OPT) $(EXTRA_FLAGS) $(CDEPS) $(VEC) $^ -o $@ $(LDEPS)

FK_original_withvecflags_bm : Benchmark.o Globals_FK.o OriginalTest_FK.cpp Vmath.o
	$(CC) -std=$(STD) $(OPT) $(EXTRA_FLAGS) $(CDEPS) $(VEC) $^ -o $@ $(LDEPS)

FK_original_bm : Benchmark.cpp Globals_FK.cpp OriginalTest_FK.cpp Vmath.cpp
	$(CC) -std=$(STD) $(OPT) $(EXTRA_FLAGS) $(CDEPS) $^ -o $@ -lbenchmark -lpthread


FK_chunked_interchanged_multithreaded_ap : ActionPotentialTest.o Globals_FK.o ChunkedTimeIntegrate.cpp ChunkedArrayInterchangeTest_FK.cpp Vmath.o
	$(CC) -std=$(STD) $(FK_DEFS) $(OPT) $(EXTRA_FLAGS) -DMULTITHREAD $(CDEPS) $(VEC) $^ -o $@ $(LDEPS)

FK_chunked_interchanged_ap : ActionPotentialTest.o Globals_FK.o ChunkedTimeIntegrate.cpp ChunkedArrayInterchangeTest_FK.cpp Vmath.o
	$(CC) -std=$(STD) $(FK_DEFS) $(OPT) $(EXTRA_FLAGS) $(CDEPS) $(VEC) $^ -o $@ $(LDEPS)

FK_original_withvecflags_ap : ActionPotentialTest.o Globals_FK.o OriginalTest_FK.cpp Vmath.o
	$(CC) -std=$(STD) $(OPT) $(EXTRA_FLAGS) $(CDEPS) $(VEC) $^ -o $@ $(LDEPS)

FK_original_ap : ActionPotentialTest.cpp Globals_FK.cpp OriginalTest_FK.cpp Vmath.cpp
	$(CC) -std=$(STD) $(OPT) $(EXTRA_FLAGS) $(CDEPS) $^ -o $@ -lbenchmark -lpthread


Benchmark.o : Benchmark.cpp
	$(CC) -std=$(STD) $(OPT) $(EXTRA_FLAGS) $(CDEPS) $(VEC) -c $^ -o $@

ActionPotentialTest.o : ActionPotentialTest.cpp
	$(CC) -std=$(STD) $(OPT) $(EXTRA_FLAGS) $(VEC) -c $^ -o $@

Vmath.o : Vmath.cpp
	$(CC) -std=$(STD) $(OPT) $(EXTRA_FLAGS) $(VEC) -c $^ -o $@

Globals_FK.o : Globals_FK.cpp
	$(CC) -std=$(STD) $(OPT) $(EXTRA_FLAGS) $(VEC) -c $^ -o $@

Globals_CRN98.o : Globals_CRN98.cpp
	$(CC) -std=$(STD) $(OPT) $(EXTRA_FLAGS) $(VEC) -c $^ -o $@

.PHONY : clean all courtemanche fenton actionpotential
clean :
	rm $(COURTEMANCHE_BM) $(COURTEMANCHE_AP) $(FENTONKARMA_BM) $(FENTONKARMA_AP) $(GARBAGE)


