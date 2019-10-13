#include <vector>
#include "SharedArray.hpp"
#include "VmathArray.hpp"
#include "Assertions.hpp"
#include "CellModelConstants_FN.hpp"

int nq;
unsigned int n_var;
std::vector<int> gates;
std::vector<int> concentrations;
Array<OneD, Array<OneD, NekDouble> > cellSol;
Array<OneD, Array<OneD, NekDouble> > wsp;
Array<OneD, Array<OneD, NekDouble> > gates_tau;

Array<OneD, NekDouble> m_uuu;

extern NekDouble lastTime;
extern unsigned int substeps;

void init_test(int n)
{
    gates = std::vector<int>();
    concentrations = std::vector<int>();
    concentrations.push_back(1);
    n_var = 2;
    nq = n;

    cellSol = Array<OneD, Array<OneD, NekDouble> > (n_var);
    wsp = Array<OneD, Array<OneD, NekDouble> > (n_var);
    gates_tau = Array<OneD, Array<OneD, NekDouble> > (gates.size());

    for (unsigned int i = 0; i < n_var; ++i) {
	cellSol[i] = Array<OneD, NekDouble> (nq);
	wsp[i] = Array<OneD, NekDouble> (nq);
    }
    for (unsigned int i = 0; i < gates.size(); ++i) {
	gates_tau[i] = Array<OneD, NekDouble> (nq);
    }

    m_uuu = Array<OneD, NekDouble> (nq, 0.0);

    Vmath::Fill(nq, 0.0, cellSol[0], 1);
    Vmath::Fill(nq, 0.0, cellSol[1], 1);
}

void v_Update(
	      const Array<OneD, const  Array<OneD, NekDouble> >&inarray,
	            Array<OneD,        Array<OneD, NekDouble> >&outarray,
	      const NekDouble)
{
    NekDouble m_gamma = 0.5;

    // compute u^2: m_u = u*u
    Vmath::Vmul(nq, &inarray[0][0], 1, &inarray[0][0], 1, &m_uuu[0], 1);

    // compute u^3: m_u = u*u*u
    Vmath::Vmul(nq, &inarray[0][0], 1, &m_uuu[0], 1, &m_uuu[0], 1);

    // For u: (1/m_epsilon)*( u*-u*u*u/3 - v )
    // physfield = u - (1.0/3.0)*u*u*u
    Vmath::Svtvp(nq, (-1.0/3.0), &m_uuu[0], 1, &inarray[0][0], 1, &outarray[0][0], 1);

    Vmath::Vsub(nq, &inarray[1][0], 1, &outarray[0][0], 1, &outarray[0][0], 1);
    Vmath::Smul(nq, -1.0/m_epsilon, &outarray[0][0], 1, &outarray[0][0], 1);

    // For v: m_epsilon*( u + m_beta - m_gamma*v )
    Vmath::Svtvp(nq, -1.0*m_gamma, &inarray[1][0], 1, &inarray[0][0], 1, &outarray[1][0], 1);
    Vmath::Sadd(nq, m_beta, &outarray[1][0], 1, &outarray[1][0], 1);
    Vmath::Smul(nq, m_epsilon, &outarray[1][0], 1, &outarray[1][0], 1);
}
