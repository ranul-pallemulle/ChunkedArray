CC=g++
OPT=-O3
STD=c++11
EXTRA_FLAGS=-g -Wall -Wextra -Werror -pedantic # -fno-omit-frame-pointer -save-temps
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
VEC=-march=sandybridge -ftree-vectorize -fopenmp -ffast-math -Rpass=loop-vectorize -Rpass-missed=loop-vectorize -fno-math-errno
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
ALIEVPANFILOV_BM=AP_chunked_interchanged_multithreaded_bm \
		 AP_chunked_interchanged_bm \
		 AP_original_withvecflags_bm \
		 AP_original_bm
ALIEVPANFILOV_AP=AP_chunked_interchanged_multithreaded_ap \
		 AP_chunked_interchanged_ap \
		 AP_original_withvecflags_ap \
		 AP_original_ap
FITZHUGHNAGUMO_BM=FN_chunked_interchanged_multithreaded_bm \
		  FN_chunked_interchanged_bm \
		  FN_original_withvecflags_bm \
		  FN_original_bm
FITZHUGHNAGUMO_AP=FN_chunked_interchanged_multithreaded_ap \
		  FN_chunked_interchanged_ap \
		  FN_original_withvecflags_ap \
		  FN_original_ap
FOX_BM=F02_original_withvecflags_bm \
       F02_original_bm

FOX_AP=F02_original_withvecflags_ap \
       F02_original_ap

LUO_BM=LR91_chunked_interchanged_multithreaded_bm \
	LR91_chunked_interchanged_bm \
	LR91_original_withvecflags_bm \
	LR91_original_bm

LUO_AP=LR91_chunked_interchanged_multithreaded_ap \
	LR91_chunked_interchanged_ap \
	LR91_original_withvecflags_ap \
	LR91_original_ap

PANDIT_BM=PGD03_original_withvecflags_bm \
	 PGD03_original_bm

PANDIT_AP=PGD03_original_withvecflags_ap \
	 PGD03_original_ap

FK_DEFS=-DNUM_IN_OUT_VARS=3 -DHAS_GATES -DNUM_GATE_VARS=2
CRN98_DEFS=-DNUM_IN_OUT_VARS=21 -DHAS_GATES -DNUM_GATE_VARS=15 -DHAS_CONCENTRATIONS
AP_DEFS=-DNUM_IN_OUT_VARS=2 -DHAS_CONCENTRATIONS
FN_DEFS=-DNUM_IN_OUT_VARS=2 -DHAS_CONCENTRATIONS
F02_DEFS=-DNUM_IN_OUT_VARS=13 -DHAS_GATES -DNUM_GATE_VARS=10 -DHAS_CONCENTRATIONS
LR91_DEFS=-DNUM_IN_OUT_VARS=8 -DHAS_GATES -DNUM_GATE_VARS=6 -DHAS_CONCENTRATIONS
PGD03_DEFS=-DNUM_IN_OUT_VARS=27

all: $(COURTEMANCHE_BM) $(FENTONKARMA_BM) $(COURTEMANCHE_AP) $(FENTONKARMA_AP) $(ALIEVPANFILOV_BM) $(ALIEVPANFILOV_AP) $(FITZHUGHNAGUMO_BM) $(FITZHUGHNAGUMO_AP) $(FOX_BM) $(FOX_AP) $(LUO_BM) $(LUO_AP) $(PANDIT_BM) $(PANDIT_AP)

courtemanche: $(COURTEMANCHE_BM) $(COURTEMANCHE_AP)

fenton: $(FENTONKARMA_BM) $(FENTONKARMA_AP)

aliev: $(ALIEVPANFILOV_BM) $(ALIEVPANFILOV_AP)

fitz: $(FITZHUGHNAGUMO_BM) $(FITZHUGHNAGUMO_AP)

fox: $(FOX_BM) $(FOX_AP)

luo: $(LUO_BM) $(LUO_AP)

pandit: $(PANDIT_BM) $(PANDIT_AP)

actionpotential: $(COURTEMANCHE_AP) $(FENTONKARMA_AP) $(ALIEVPANFILOV_AP) $(FITZHUGHNAGUMO_AP) $(FOX_AP) $(LUO_AP) $(PANDIT_AP)

benchmark: $(COURTEMANCHE_BM) $(FENTONKARMA_BM) $(ALIEVPANFILOV_BM) $(FITZHUGHNAGUMO_BM) $(FOX_BM) $(LUO_BM) $(PANDIT_BM)


CRN98_chunked_interchanged_multithread_bm : Benchmark.o Globals_CRN98.o ChunkedTimeIntegrate.cpp ChunkedArrayInterchangeTest_CRN98.cpp Vmath.o
	$(CC) -std=$(STD) $(CRN98_DEFS) $(OPT) $(EXTRA_FLAGS) -DMULTITHREAD $(CDEPS) $(VEC) $^ -o $@ $(LDEPS)

CRN98_chunked_interchanged_splitloops_multithread_bm : Benchmark.o Globals_CRN98.o ChunkedTimeIntegrate.cpp ChunkedArrayInterchangeSplitLoopsTest_CRN98.cpp Vmath.o
	$(CC) -std=$(STD) $(CRN98_DEFS) $(OPT) $(EXTRA_FLAGS) -DMULTITHREAD $(CDEPS) $(VEC) $^ -o $@ $(LDEPS)

CRN98_chunked_interchanged_splitloops_bm : Benchmark.o Globals_CRN98.o ChunkedTimeIntegrate.cpp ChunkedArrayInterchangeSplitLoopsTest_CRN98.cpp Vmath.o
	$(CC) -std=$(STD) $(CRN98_DEFS) $(OPT) $(EXTRA_FLAGS) $(CDEPS) $(VEC) $^ -o $@ $(LDEPS)

CRN98_chunked_interchanged_bm : Benchmark.o Globals_CRN98.o ChunkedTimeIntegrate.cpp ChunkedArrayInterchangeTest_CRN98.cpp Vmath.o
	$(CC) -std=$(STD) $(CRN98_DEFS) $(OPT) $(EXTRA_FLAGS) $(CDEPS) $(VEC) $^ -o $@ $(LDEPS)

CRN98_original_withvecflags_bm : Benchmark.o Globals_CRN98.o OriginalTimeIntegrate.cpp OriginalTest_CRN98.cpp Vmath.o
	$(CC) -std=$(STD) $(OPT) $(EXTRA_FLAGS) $(CDEPS) $(VEC) $^ -o $@ $(LDEPS)

CRN98_original_bm : Benchmark.cpp Globals_CRN98.cpp OriginalTimeIntegrate.cpp OriginalTest_CRN98.cpp Vmath.cpp
	$(CC) -std=$(STD) $(OPT) $(EXTRA_FLAGS) $(CDEPS) $^ -o $@ -lbenchmark -lpthread


CRN98_chunked_interchanged_multithread_ap : ActionPotentialTest.o Globals_CRN98.o ChunkedTimeIntegrate.cpp ChunkedArrayInterchangeTest_CRN98.cpp Vmath.o
	$(CC) -std=$(STD) $(CRN98_DEFS) $(OPT) $(EXTRA_FLAGS) -DMULTITHREAD $(CDEPS) $(VEC) $^ -o $@ $(LDEPS)

CRN98_chunked_interchanged_splitloops_multithread_ap : ActionPotentialTest.o Globals_CRN98.o ChunkedTimeIntegrate.cpp ChunkedArrayInterchangeSplitLoopsTest_CRN98.cpp Vmath.o
	$(CC) -std=$(STD) $(CRN98_DEFS) $(OPT) $(EXTRA_FLAGS) -DMULTITHREAD $(CDEPS) $(VEC) $^ -o $@ $(LDEPS)

CRN98_chunked_interchanged_splitloops_ap : ActionPotentialTest.o Globals_CRN98.o ChunkedTimeIntegrate.cpp ChunkedArrayInterchangeSplitLoopsTest_CRN98.cpp Vmath.o
	$(CC) -std=$(STD) $(CRN98_DEFS) $(OPT) $(EXTRA_FLAGS) $(CDEPS) $(VEC) $^ -o $@ $(LDEPS)

CRN98_chunked_interchanged_ap : ActionPotentialTest.o Globals_CRN98.o ChunkedTimeIntegrate.cpp ChunkedArrayInterchangeTest_CRN98.cpp Vmath.o
	$(CC) -std=$(STD) $(CRN98_DEFS) $(OPT) $(EXTRA_FLAGS) $(CDEPS) $(VEC) $^ -o $@ $(LDEPS)

CRN98_original_withvecflags_ap : ActionPotentialTest.o Globals_CRN98.o OriginalTimeIntegrate.cpp OriginalTest_CRN98.cpp Vmath.o
	$(CC) -std=$(STD) $(OPT) $(EXTRA_FLAGS) $(CDEPS) $(VEC) $^ -o $@ $(LDEPS)

CRN98_original_ap : ActionPotentialTest.cpp Globals_CRN98.cpp OriginalTimeIntegrate.cpp OriginalTest_CRN98.cpp Vmath.cpp
	$(CC) -std=$(STD) $(OPT) $(EXTRA_FLAGS) $(CDEPS) $^ -o $@ -lbenchmark -lpthread


FK_chunked_interchanged_multithreaded_bm : Benchmark.o Globals_FK.o ChunkedTimeIntegrate.cpp ChunkedArrayInterchangeTest_FK.cpp Vmath.o
	$(CC) -std=$(STD) $(FK_DEFS) $(OPT) $(EXTRA_FLAGS) -DMULTITHREAD $(CDEPS) $(VEC) $^ -o $@ $(LDEPS)

FK_chunked_interchanged_bm : Benchmark.o Globals_FK.o ChunkedTimeIntegrate.cpp ChunkedArrayInterchangeTest_FK.cpp Vmath.o
	$(CC) -std=$(STD) $(FK_DEFS) $(OPT) $(EXTRA_FLAGS) $(CDEPS) $(VEC) $^ -o $@ $(LDEPS)

FK_original_withvecflags_bm : Benchmark.o Globals_FK.o OriginalTimeIntegrate.cpp OriginalTest_FK.cpp Vmath.o
	$(CC) -std=$(STD) $(OPT) $(EXTRA_FLAGS) $(CDEPS) $(VEC) $^ -o $@ $(LDEPS)

FK_original_bm : Benchmark.cpp Globals_FK.cpp OriginalTimeIntegrate.cpp OriginalTest_FK.cpp Vmath.cpp
	$(CC) -std=$(STD) $(OPT) $(EXTRA_FLAGS) $(CDEPS) $^ -o $@ -lbenchmark -lpthread


FK_chunked_interchanged_multithreaded_ap : ActionPotentialTest.o Globals_FK.o ChunkedTimeIntegrate.cpp ChunkedArrayInterchangeTest_FK.cpp Vmath.o
	$(CC) -std=$(STD) $(FK_DEFS) $(OPT) $(EXTRA_FLAGS) -DMULTITHREAD $(CDEPS) $(VEC) $^ -o $@ $(LDEPS)

FK_chunked_interchanged_ap : ActionPotentialTest.o Globals_FK.o ChunkedTimeIntegrate.cpp ChunkedArrayInterchangeTest_FK.cpp Vmath.o
	$(CC) -std=$(STD) $(FK_DEFS) $(OPT) $(EXTRA_FLAGS) $(CDEPS) $(VEC) $^ -o $@ $(LDEPS)

FK_original_withvecflags_ap : ActionPotentialTest.o Globals_FK.o OriginalTimeIntegrate.cpp OriginalTest_FK.cpp Vmath.o
	$(CC) -std=$(STD) $(OPT) $(EXTRA_FLAGS) $(CDEPS) $(VEC) $^ -o $@ $(LDEPS)

FK_original_ap : ActionPotentialTest.cpp Globals_FK.cpp OriginalTimeIntegrate.cpp OriginalTest_FK.cpp Vmath.cpp
	$(CC) -std=$(STD) $(OPT) $(EXTRA_FLAGS) $(CDEPS) $^ -o $@ -lbenchmark -lpthread



AP_chunked_interchanged_multithreaded_bm : Benchmark.o Globals_AP.o ChunkedTimeIntegrate.cpp ChunkedArrayInterchangeTest_AP.cpp Vmath.o
	$(CC) -std=$(STD) $(AP_DEFS) $(OPT) $(EXTRA_FLAGS) -DMULTITHREAD $(CDEPS) $(VEC) $^ -o $@ $(LDEPS)

AP_chunked_interchanged_bm : Benchmark.o Globals_AP.o ChunkedTimeIntegrate.cpp ChunkedArrayInterchangeTest_AP.cpp Vmath.o
	$(CC) -std=$(STD) $(AP_DEFS) $(OPT) $(EXTRA_FLAGS) $(CDEPS) $(VEC) $^ -o $@ $(LDEPS)

AP_original_withvecflags_bm : Benchmark.o Globals_AP.o OriginalTimeIntegrate.cpp OriginalTest_AP.cpp Vmath.o
	$(CC) -std=$(STD) $(OPT) $(EXTRA_FLAGS) $(CDEPS) $(VEC) $^ -o $@ $(LDEPS)

AP_original_bm : Benchmark.cpp Globals_AP.cpp OriginalTimeIntegrate.cpp OriginalTest_AP.cpp Vmath.cpp
	$(CC) -std=$(STD) $(OPT) $(EXTRA_FLAGS) $(CDEPS) $^ -o $@ -lbenchmark -lpthread



AP_chunked_interchanged_multithreaded_ap : ActionPotentialTest.o Globals_AP.o ChunkedTimeIntegrate.cpp ChunkedArrayInterchangeTest_AP.cpp Vmath.o
	$(CC) -std=$(STD) $(AP_DEFS) $(OPT) $(EXTRA_FLAGS) -DMULTITHREAD $(CDEPS) $(VEC) $^ -o $@ $(LDEPS)

AP_chunked_interchanged_ap : ActionPotentialTest.o Globals_AP.o ChunkedTimeIntegrate.cpp ChunkedArrayInterchangeTest_AP.cpp Vmath.o
	$(CC) -std=$(STD) $(AP_DEFS) $(OPT) $(EXTRA_FLAGS) $(CDEPS) $(VEC) $^ -o $@ $(LDEPS)

AP_original_withvecflags_ap : ActionPotentialTest.o Globals_AP.o OriginalTimeIntegrate.cpp OriginalTest_AP.cpp Vmath.o
	$(CC) -std=$(STD) $(OPT) $(EXTRA_FLAGS) $(CDEPS) $(VEC) $^ -o $@ $(LDEPS)

AP_original_ap : ActionPotentialTest.cpp Globals_AP.cpp OriginalTimeIntegrate.cpp OriginalTest_AP.cpp Vmath.cpp
	$(CC) -std=$(STD) $(OPT) $(EXTRA_FLAGS) $(CDEPS) $^ -o $@ -lbenchmark -lpthread


FN_chunked_interchanged_multithreaded_bm : Benchmark.o Globals_FN.o ChunkedTimeIntegrate.cpp ChunkedArrayInterchangeTest_FN.cpp Vmath.o
	$(CC) -std=$(STD) $(FN_DEFS) $(OPT) $(EXTRA_FLAGS) -DMULTITHREAD $(CDEPS) $(VEC) $^ -o $@ $(LDEPS)

FN_chunked_interchanged_bm : Benchmark.o Globals_FN.o ChunkedTimeIntegrate.cpp ChunkedArrayInterchangeTest_FN.cpp Vmath.o
	$(CC) -std=$(STD) $(FN_DEFS) $(OPT) $(EXTRA_FLAGS) $(CDEPS) $(VEC) $^ -o $@ $(LDEPS)

FN_original_withvecflags_bm : Benchmark.o Globals_FN.o OriginalTimeIntegrate.cpp OriginalTest_FN.cpp Vmath.o
	$(CC) -std=$(STD) $(OPT) $(EXTRA_FLAGS) $(CDEPS) $(VEC) $^ -o $@ $(LDEPS)

FN_original_bm : Benchmark.cpp Globals_FN.cpp OriginalTimeIntegrate.cpp OriginalTest_FN.cpp Vmath.cpp
	$(CC) -std=$(STD) $(OPT) $(EXTRA_FLAGS) $(CDEPS) $^ -o $@ -lbenchmark -lpthread


FN_chunked_interchanged_multithreaded_ap : ActionPotentialTest.o Globals_FN.o ChunkedTimeIntegrate.cpp ChunkedArrayInterchangeTest_FN.cpp Vmath.o
	$(CC) -std=$(STD) $(FN_DEFS) $(OPT) $(EXTRA_FLAGS) -DMULTITHREAD $(CDEPS) $(VEC) $^ -o $@ $(LDEPS)

FN_chunked_interchanged_ap : ActionPotentialTest.o Globals_FN.o ChunkedTimeIntegrate.cpp ChunkedArrayInterchangeTest_FN.cpp Vmath.o
	$(CC) -std=$(STD) $(FN_DEFS) $(OPT) $(EXTRA_FLAGS) $(CDEPS) $(VEC) $^ -o $@ $(LDEPS)

FN_original_withvecflags_ap : ActionPotentialTest.o Globals_FN.o OriginalTimeIntegrate.cpp OriginalTest_FN.cpp Vmath.o
	$(CC) -std=$(STD) $(OPT) $(EXTRA_FLAGS) $(CDEPS) $(VEC) $^ -o $@ $(LDEPS)

FN_original_ap : ActionPotentialTest.cpp Globals_FN.cpp OriginalTimeIntegrate.cpp OriginalTest_FN.cpp Vmath.cpp
	$(CC) -std=$(STD) $(OPT) $(EXTRA_FLAGS) $(CDEPS) $^ -o $@ -lbenchmark -lpthread


F02_original_withvecflags_bm : Benchmark.o Globals_F02.o OriginalTimeIntegrate.cpp OriginalTest_F02.cpp Vmath.o
	$(CC) -std=$(STD) $(OPT) $(EXTRA_FLAGS) $(CDEPS) $(VEC) $^ -o $@ $(LDEPS)

F02_original_bm : Benchmark.cpp Globals_F02.cpp OriginalTimeIntegrate.cpp OriginalTest_F02.cpp Vmath.cpp
	$(CC) -std=$(STD) $(OPT) $(EXTRA_FLAGS) $(CDEPS) $^ -o $@ -lbenchmark -lpthread


F02_original_withvecflags_ap : ActionPotentialTest.o Globals_F02.o OriginalTimeIntegrate.cpp OriginalTest_F02.cpp Vmath.o
	$(CC) -std=$(STD) $(OPT) $(EXTRA_FLAGS) $(CDEPS) $(VEC) $^ -o $@ $(LDEPS)

F02_original_ap : ActionPotentialTest.cpp Globals_F02.cpp OriginalTimeIntegrate.cpp OriginalTest_F02.cpp Vmath.cpp
	$(CC) -std=$(STD) $(OPT) $(EXTRA_FLAGS) $(CDEPS) $^ -o $@ -lbenchmark -lpthread


LR91_chunked_interchanged_multithreaded_bm : Benchmark.o Globals_LR91.o ChunkedTimeIntegrate.cpp ChunkedArrayInterchangeTest_LR91.cpp Vmath.o
	$(CC) -std=$(STD) $(LR91_DEFS) $(OPT) $(EXTRA_FLAGS) -DMULTITHREAD $(CDEPS) $(VEC) $^ -o $@ $(LDEPS)

LR91_chunked_interchanged_bm : Benchmark.o Globals_LR91.o ChunkedTimeIntegrate.cpp ChunkedArrayInterchangeTest_LR91.cpp Vmath.o
	$(CC) -std=$(STD) $(LR91_DEFS) $(OPT) $(EXTRA_FLAGS) $(CDEPS) $(VEC) $^ -o $@ $(LDEPS)

LR91_original_withvecflags_bm : Benchmark.o Globals_LR91.o OriginalTimeIntegrate.cpp OriginalTest_LR91.cpp Vmath.o
	$(CC) -std=$(STD) $(OPT) $(EXTRA_FLAGS) $(CDEPS) $(VEC) $^ -o $@ $(LDEPS)

LR91_original_bm : Benchmark.cpp Globals_LR91.cpp OriginalTimeIntegrate.cpp OriginalTest_LR91.cpp Vmath.cpp
	$(CC) -std=$(STD) $(OPT) $(EXTRA_FLAGS) $(CDEPS) $^ -o $@ -lbenchmark -lpthread


LR91_chunked_interchanged_multithreaded_ap : ActionPotentialTest.o Globals_LR91.o ChunkedTimeIntegrate.cpp ChunkedArrayInterchangeTest_LR91.cpp Vmath.o
	$(CC) -std=$(STD) $(LR91_DEFS) $(OPT) $(EXTRA_FLAGS) -DMULTITHREAD $(CDEPS) $(VEC) $^ -o $@ $(LDEPS)

LR91_chunked_interchanged_ap : ActionPotentialTest.o Globals_LR91.o ChunkedTimeIntegrate.cpp ChunkedArrayInterchangeTest_LR91.cpp Vmath.o
	$(CC) -std=$(STD) $(LR91_DEFS) $(OPT) $(EXTRA_FLAGS) $(CDEPS) $(VEC) $^ -o $@ $(LDEPS)

LR91_original_withvecflags_ap : ActionPotentialTest.o Globals_LR91.o OriginalTimeIntegrate.cpp OriginalTest_LR91.cpp Vmath.o
	$(CC) -std=$(STD) $(OPT) $(EXTRA_FLAGS) $(CDEPS) $(VEC) $^ -o $@ $(LDEPS)

LR91_original_ap : ActionPotentialTest.cpp Globals_LR91.cpp OriginalTimeIntegrate.cpp OriginalTest_LR91.cpp Vmath.cpp
	$(CC) -std=$(STD) $(OPT) $(EXTRA_FLAGS) $(CDEPS) $^ -o $@ -lbenchmark -lpthread


PGD03_original_withvecflags_bm : Benchmark.o Globals_PGD03.o OriginalTimeIntegrate.cpp OriginalTest_PGD03.cpp Vmath.o
	$(CC) -std=$(STD) $(OPT) $(EXTRA_FLAGS) $(CDEPS) $(VEC) $^ -o $@ $(LDEPS)

PGD03_original_bm : Benchmark.cpp Globals_PGD03.cpp OriginalTimeIntegrate.cpp OriginalTest_PGD03.cpp Vmath.cpp
	$(CC) -std=$(STD) $(OPT) $(EXTRA_FLAGS) $(CDEPS) $^ -o $@ -lbenchmark -lpthread

PGD03_original_withvecflags_ap : ActionPotentialTest.o Globals_PGD03.o OriginalTimeIntegrate.cpp OriginalTest_PGD03.cpp Vmath.o
	$(CC) -std=$(STD) $(OPT) $(EXTRA_FLAGS) $(CDEPS) $(VEC) $^ -o $@ $(LDEPS)

PGD03_original_ap : ActionPotentialTest.cpp Globals_PGD03.cpp OriginalTimeIntegrate.cpp OriginalTest_PGD03.cpp Vmath.cpp
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

Globals_AP.o : Globals_AP.cpp
	$(CC) -std=$(STD) $(OPT) $(EXTRA_FLAGS) $(VEC) -c $^ -o $@

Globals_FN.o : Globals_FN.cpp
	$(CC) -std=$(STD) $(OPT) $(EXTRA_FLAGS) $(VEC) -c $^ -o $@

Globals_F02.o : Globals_F02.cpp
	$(CC) -std=$(STD) $(OPT) $(EXTRA_FLAGS) $(VEC) -c $^ -o $@

Globals_LR91.o : Globals_LR91.cpp
	$(CC) -std=$(STD) $(OPT) $(EXTRA_FLAGS) $(VEC) -c $^ -o $@

Globals_PGD03.o : Globals_PGD03.cpp
	$(CC) -std=$(STD) $(OPT) $(EXTRA_FLAGS) $(VEC) -c $^ -o $@

.PHONY : clean all courtemanche fenton aliev fitz fox luo pandit actionpotential benchmark
clean :
	rm $(COURTEMANCHE_BM) $(COURTEMANCHE_AP) $(FENTONKARMA_BM) $(FENTONKARMA_AP) $(ALIEVPANFILOV_BM) $(ALIEVPANFILOV_AP) $(FITZHUGHNAGUMO_BM) $(FITZHUGHNAGUMO_AP) $(FOX_BM) $(FOX_AP) $(LUO_BM) $(LUO_AP) $(PANDIT_BM) $(PANDIT_AP) $(GARBAGE)


