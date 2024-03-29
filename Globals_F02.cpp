#include "SharedArray.hpp"

#ifndef NekDouble
#define NekDouble double
#endif

Array<OneD, Array<OneD, NekDouble> > vSol;
Array<OneD, Array<OneD, NekDouble> > vWsp;

unsigned int substeps = 2;
NekDouble timeStep = 0.02;
NekDouble finTime = 100;

NekDouble stimStrength = 100;
NekDouble stimSteps = 50;

void init_globals(int n)
{
    vSol = Array<OneD, Array<OneD, NekDouble> >(1);
    vWsp = Array<OneD, Array<OneD, NekDouble> >(1);
    vSol[0] = Array<OneD, NekDouble> (n, 0.0);
    vWsp[0] = Array<OneD, NekDouble> (n, 0.0);
}
