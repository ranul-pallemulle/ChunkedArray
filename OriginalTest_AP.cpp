#include <vector>
#include "SharedArray.hpp"
#include "VmathArray.hpp"
#include "Assertions.hpp"
#include "CellModelConstants_AP.hpp"

int nq;
unsigned int n_var;
std::vector<int> gates;
std::vector<int> concentrations;
Array<OneD, Array<OneD, NekDouble> > cellSol;
Array<OneD, Array<OneD, NekDouble> > wsp;
Array<OneD, Array<OneD, NekDouble> > gates_tau;

Array<OneD, NekDouble> uu;
Array<OneD, NekDouble> uuu;
Array<OneD, NekDouble> tmp1;
Array<OneD, NekDouble> tmp2;

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

    uu = Array<OneD, NekDouble> (nq, 0.0);
    uuu = Array<OneD, NekDouble> (nq, 0.0);
    tmp1 = Array<OneD, NekDouble> (nq, 0.0);
    tmp2 = Array<OneD, NekDouble> (nq, 0.0);

    Vmath::Fill(nq, 0.0, cellSol[0], 1);
    Vmath::Fill(nq, 0.0, cellSol[1], 1);
}

void v_Update(const Array<OneD, const Array<OneD, NekDouble> >& inarray,
	            Array<OneD,       Array<OneD, NekDouble> >& outarray,
	      const NekDouble)
{
    // inarray[0] holds initial physical u values throughout
    // inarray[1] holds initial physical v values throughout

    // compute u^2: m_u = u*u
    Vmath::Vmul(nq, &inarray[0][0], 1, &inarray[0][0], 1, &uu[0], 1);

    // compute u^3: m_u = u*u*u
    Vmath::Vmul(nq, &inarray[0][0], 1, &uu[0], 1, &uuu[0], 1);

    // --------------------------------------
    // Compute reaction term f(u,v)
    // --------------------------------------
    // Ru = au
    Vmath::Smul(nq, m_a, &inarray[0][0], 1, &tmp1[0], 1);
    // Ru = (-1-a)u*u + au
    Vmath::Svtvp(nq, (-1.0-m_a), &uu[0], 1, &tmp1[0], 1,
		 &tmp1[0], 1);
    // Ru = u*u*u - (1+a)u*u + au
    Vmath::Vadd(nq, &uuu[0], 1, &tmp1[0], 1, &tmp1[0], 1);
    // Ru = k(u*u*u - (1+a)u*u + au)
    Vmath::Smul(nq, m_k, &tmp1[0], 1, &tmp1[0], 1);
    // Ru = k(u*u*u - (1+a)u*u + au) + I_stim
    Vmath::Vadd(nq, &outarray[0][0], 1, &tmp1[0], 1, &outarray[0][0], 1);
    // Ru = k(u*u*u - (1+a)u*u + au) + uv + I_stim
    Vmath::Vvtvp(nq, &inarray[0][0], 1, &inarray[1][0], 1, &tmp1[0], 1,
		 &outarray[0][0], 1);
    // Ru = -k(u*u*u - (1+a)u*u + au) - uv - I_stim
    Vmath::Neg(nq, &outarray[0][0], 1);


    // --------------------------------------
    // Compute reaction term g(u,v)
    // --------------------------------------
    // tmp2 = mu2 + u
    Vmath::Sadd(nq, m_mu2, &inarray[0][0], 1, &tmp2[0], 1);

    // tmp2 = v/(mu2 + u)
    Vmath::Vdiv(nq, &inarray[1][0], 1, &tmp2[0], 1, &tmp2[0], 1);

    // tmp2 = mu1*v/(mu2 + u)
    Vmath::Smul(nq, m_mu1, &tmp2[0], 1, &tmp2[0], 1);

    // tmp2 = Eps + mu1*v/(mu2+u)
    Vmath::Sadd(nq, m_eps, &tmp2[0], 1, &tmp2[0], 1);

    // tmp1 = (-a-1) + u
    Vmath::Sadd(nq, (-m_a-1), &inarray[0][0], 1, &tmp1[0], 1);

    // tmp1 = k(u-a-1)
    Vmath::Smul(nq, m_k, &tmp1[0], 1, &tmp1[0], 1);

    // tmp1 = ku(u-a-1) + v
    Vmath::Vvtvp(nq, &inarray[0][0], 1, &tmp1[0], 1, &inarray[1][0], 1,
		 &tmp1[0], 1);

    // tmp1 = -ku(u-a-1)-v
    Vmath::Neg(nq, &tmp1[0], 1);

    // outarray = [Eps + mu1*v/(mu2+u)] * [-ku(u-a-1)-v]
    Vmath::Vmul(nq, &tmp1[0], 1, &tmp2[0], 1, &outarray[1][0], 1);
}
