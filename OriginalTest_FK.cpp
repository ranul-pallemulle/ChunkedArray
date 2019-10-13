#include <vector>
#include "SharedArray.hpp"
#include "VmathArray.hpp"
#include "Assertions.hpp"
#include "CellModelConstants_FK.hpp"

int nq;
unsigned int n_var;
std::vector<int> gates;
std::vector<int> concentrations;
Array<OneD, Array<OneD, NekDouble> > cellSol;
Array<OneD, Array<OneD, NekDouble> > wsp;
Array<OneD, Array<OneD, NekDouble> > gates_tau;

extern NekDouble lastTime;
extern unsigned int substeps;

void init_test(int n)
{
    gates = std::vector<int>();
    concentrations = std::vector<int>();
    gates.push_back(1);
    gates.push_back(2);
    n_var = 3;
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

    Vmath::Fill(nq, (NekDouble)0.0, cellSol[0], 1);
    Vmath::Fill(nq, (NekDouble)1.0, cellSol[1], 1);
    Vmath::Fill(nq, (NekDouble)1.0, cellSol[2], 1);
}

void v_Update(const Array<OneD, const Array<OneD, NekDouble> >& inarray,
	            Array<OneD,       Array<OneD, NekDouble> >& outarray,
	      const NekDouble)
{
    ASSERTL0(inarray.get() != outarray.get(),
	     "Must have different arrays for input and output.");

    int n = nq;
    int i = 0;

    const NekDouble *u = &inarray[0][0];
    const NekDouble *v = &inarray[1][0];
    const NekDouble *w = &inarray[2][0];
    NekDouble *u_new = &outarray[0][0];
    NekDouble *v_new = &outarray[1][0];
    NekDouble *w_new = &outarray[2][0];
    NekDouble *v_tau = &gates_tau[0][0];
    NekDouble *w_tau = &gates_tau[1][0];

    NekDouble J_fi, J_so, J_si, h1, h2, h3, alpha, beta;
    NekDouble V;

    for (i = 0; i < n; ++i) {

	V = (*u - V_0)/(V_fi - V_0);

	h1 = (V < u_c) ? 0.0 : 1.0;
	h2 = (V < u_v) ? 0.0 : 1.0;
	h3 = (V < u_r) ? 0.0 : 1.0;

	alpha = (1 - h1)/tau_w_minus;
	beta = h1/tau_w_plus;
	*w_tau = 1.0/(alpha + beta);
	*w_new = alpha * (*w_tau);

	alpha = (1 - h1)/(h2*tau_v1_minus + (1-h2)*tau_v2_minus);
	beta = h1/tau_v_plus;
	*v_tau = 1.0/(alpha + beta);
	*v_new = alpha * (*v_tau);

	J_fi = -(*v)*h1*(1 - V)*(V - u_c)/tau_d;

	J_so = V*(1-h3)*(1-k2*(*v))/tau_0 + h3/tau_r;

	J_si = -(*w)*(1 + tanh(k1*(V - u_csi)))/(2.0*tau_si);

	*u_new = -J_fi - J_so - J_si;
	*u_new *= C_m*(V_fi - V_0);

	++u, ++v, ++w, ++u_new, ++v_new, ++w_new, ++v_tau, ++w_tau;
    }
}
