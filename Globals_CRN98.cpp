#include "SharedArray.hpp"

#ifndef NekDouble
#define NekDouble double
#endif

Array<OneD, Array<OneD, NekDouble> > vSol;
Array<OneD, Array<OneD, NekDouble> > vWsp;

unsigned int substeps = 5;
NekDouble timeStep = 0.02;
NekDouble finTime = 30;

void init_globals(int n)
{
    vSol = Array<OneD, Array<OneD, NekDouble> > (1);
    vWsp = Array<OneD, Array<OneD, NekDouble> > (1);
    vSol[0] = Array<OneD, NekDouble> (n, -81.19);
    vWsp[0] = Array<OneD, NekDouble> (n, 0.0);
}
