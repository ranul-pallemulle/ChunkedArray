#include "SharedArray.hpp"

#ifndef NekDouble
#define NekDouble double
#endif

Array<OneD, Array<OneD, NekDouble> > vSol;
Array<OneD, Array<OneD, NekDouble> > vWsp;

unsigned int substeps = 10;
NekDouble timeStep = 0.001;
NekDouble finTime = 5;

NekDouble stimStrength = 50;
NekDouble stimSteps = 100;

void init_globals(int n)
{
    vSol = Array<OneD, Array<OneD, NekDouble> >(1);
    vWsp = Array<OneD, Array<OneD, NekDouble> >(1);
    vSol[0] = Array<OneD, NekDouble> (n, -84.3801107371);
    vWsp[0] = Array<OneD, NekDouble> (n, 0.0);
}
