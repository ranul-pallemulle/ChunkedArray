#include <cmath>
#include <vector>
#include "SharedArray.hpp"
#include "VmathArray.hpp"

#ifndef NekDouble
#define NekDouble double
#endif

extern NekDouble lastTime;
extern unsigned int substeps;
extern int nq;
extern unsigned int n_var;
extern std::vector<int> gates;
extern std::vector<int> concentrations;
extern Array<OneD, Array<OneD, NekDouble> > cellSol;
extern Array<OneD, Array<OneD, NekDouble> > wsp;
extern Array<OneD, Array<OneD, NekDouble> > gates_tau;

extern void v_Update(const Array<OneD, const Array<OneD, NekDouble> >& inarray,
		           Array<OneD,       Array<OneD, NekDouble> >& outarray,
		     const NekDouble time);

void TimeIntegrate(const Array<OneD, const Array<OneD, NekDouble> >& inarray,
		         Array<OneD,       Array<OneD, NekDouble> >& outarray,
		   const NekDouble time)
{
    Vmath::Vcopy(nq, inarray[0], 1, cellSol[0], 1);
    NekDouble delta_t = (time - lastTime)/substeps;

    for (unsigned int i = 0; i < substeps - 1; ++i)
    {
	v_Update(cellSol, wsp, time);
	Vmath::Svtvp(nq, delta_t, wsp[0], 1, cellSol[0], 1, cellSol[0], 1);
	for (unsigned int j = 0; j < concentrations.size(); ++j)
	{
	    Vmath::Svtvp(nq, delta_t, wsp[concentrations[j]], 1,
			 cellSol[concentrations[j]], 1, cellSol[concentrations[j]], 1);
	}

	for (unsigned int j = 0; j < gates.size(); ++j)
	{
	    Vmath::Sdiv(nq, -delta_t, gates_tau[j], 1, gates_tau[j], 1);
	    Vmath::Vexp(nq, gates_tau[j], 1, gates_tau[j], 1);
	    Vmath::Vsub(nq, cellSol[gates[j]], 1, wsp[gates[j]], 1, cellSol[gates[j]], 1);
	    Vmath::Vvtvp(nq, cellSol[gates[j]], 1, gates_tau[j], 1,
			 wsp[gates[j]], 1, cellSol[gates[j]], 1);
	}
    }
    v_Update(cellSol, wsp, time);

    Vmath::Vcopy(nq, wsp[0], 1, outarray[0], 1);

    for (unsigned int j = 0; j < concentrations.size(); ++j)
    {
	Vmath::Svtvp(nq, delta_t, wsp[concentrations[j]], 1,
		     cellSol[concentrations[j]], 1, cellSol[concentrations[j]], 1);
    }


    for (unsigned int j = 0; j < gates.size(); ++j)
    {
	Vmath::Sdiv(nq, -delta_t, gates_tau[j], 1, gates_tau[j], 1);
	Vmath::Vexp(nq, gates_tau[j], 1, gates_tau[j], 1);
	Vmath::Vsub(nq, cellSol[gates[j]], 1, wsp[gates[j]], 1, cellSol[gates[j]], 1);
	Vmath::Vvtvp(nq, cellSol[gates[j]], 1, gates_tau[j], 1,
		     wsp[gates[j]], 1, cellSol[gates[j]], 1);
    }
    lastTime = time;
}
