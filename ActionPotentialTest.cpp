#include <iostream>
#include <string>
#include <chrono>

#include <cfenv>

#include "SharedArray.hpp"
#include "VmathArray.hpp"

#ifndef NekDouble
#define NekDouble double
#endif

using std::chrono::time_point;
using hr_clock = std::chrono::high_resolution_clock;
using std::chrono::duration_cast;
using micros = std::chrono::microseconds;

extern void init_globals(int);
extern void init_test(int);
extern void TimeIntegrate(const Array<OneD, const Array<OneD, NekDouble> >& inarray,
			        Array<OneD,       Array<OneD, NekDouble> >& outarray,
			  const NekDouble time);

extern Array<OneD, Array<OneD, NekDouble> > vSol;
extern Array<OneD, Array<OneD, NekDouble> > vWsp;
extern unsigned int substeps;
extern NekDouble timeStep;
extern NekDouble finTime;
extern NekDouble stimStrength;
extern NekDouble stimSteps;
NekDouble lastTime = 0;
NekDouble vTime = 0.0;
unsigned int nSteps = finTime/timeStep;

int main(int argc, char** argv)
{
    feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
    int num_points;
    if (argc > 1)
    {
	num_points = std::stoi(argv[1]);
    }
    else
    {
	std::cerr << "Usage: action_potential_[cellmodel] <num_phys/coef_points>" << std::endl;
	std::exit(1);
    }
    
    init_globals(num_points);
    init_test(num_points);

    for (unsigned int i = 0; i < nSteps; ++i)
    {
	auto start = hr_clock::now();
	TimeIntegrate(vSol, vWsp, vTime);
	if (i < stimSteps)
	{
	    vWsp[0][0] = stimStrength;
	}

	Vmath::Svtvp(num_points, timeStep, vWsp[0], 1, vSol[0], 1, vSol[0], 1);
	vTime += timeStep;
	auto end = hr_clock::now();
	auto elapsed = end - start;

	std::cout << "Time: " << vTime << "    ";
	std::cout << "Voltage: " << vSol[0][0] << "    ";
	std::cout << "CPU time: " << duration_cast<micros>(elapsed).count() << "us";
	std::cout << "\n";
    }
    std::cout << std::flush;
}
