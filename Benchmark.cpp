#include <iostream>
#include "SharedArray.hpp"
#include "VmathArray.hpp"
#include "benchmark/benchmark.h"

#ifndef NekDouble
#define NekDouble double
#endif

extern void init_test(int);

extern void TimeIntegrate(const Array<OneD, const Array<OneD, NekDouble> >& inarray,
			        Array<OneD,       Array<OneD, NekDouble> >& outarray,
			  const NekDouble time);

Array<OneD, Array<OneD, NekDouble> > vSol;
Array<OneD, Array<OneD, NekDouble> > vWsp;

NekDouble lastTime = 0;
unsigned int substeps = 5;
NekDouble timeStep = 0.02;
NekDouble finTime = 30;
NekDouble vTime = 0.0;
unsigned int nSteps = finTime/timeStep;

void init(int n)
{
    vSol = Array<OneD, Array<OneD, NekDouble> > (1);
    vWsp = Array<OneD, Array<OneD, NekDouble> > (1);
    vSol[0] = Array<OneD, NekDouble> (n, -81.19);
    vWsp[0] = Array<OneD, NekDouble> (n, 0.0);
}


void BM_TimeIntegrate(benchmark::State& state)
{
    init(state.range(0));
    init_test(state.range(0));
    for (auto _ : state)
    {
	TimeIntegrate(vSol, vWsp, vTime);
    }
}

BENCHMARK(BM_TimeIntegrate)->RangeMultiplier(2)->Range(1 << 7, 1 << 20);
BENCHMARK_MAIN();
