#include <iostream>
#include "SharedArray.hpp"
#include "benchmark/benchmark.h"

#ifndef NekDouble
#define NekDouble double
#endif

extern void init_globals(int);
extern void init_test(int);
extern void TimeIntegrate(const Array<OneD, const Array<OneD, NekDouble> >& inarray,
			        Array<OneD,       Array<OneD, NekDouble> >& outarray,
			  const NekDouble time);

extern Array<OneD, Array<OneD, NekDouble> > vSol;
extern Array<OneD, Array<OneD, NekDouble> > vWsp;
extern unsigned int substeps;
extern NekDouble timeStep;
NekDouble lastTime = 0;
NekDouble vTime = 0.0;

void BM_TimeIntegrate(benchmark::State& state)
{
    init_globals(state.range(0));
    init_test(state.range(0));
    for (auto _ : state)
    {
	TimeIntegrate(vSol, vWsp, vTime);
    }
}

BENCHMARK(BM_TimeIntegrate)->RangeMultiplier(2)->Range(1 << 7, 1 << 20)->UseRealTime();
BENCHMARK_MAIN();
