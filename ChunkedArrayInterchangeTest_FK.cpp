#include <vector>
#include "SharedArray.hpp"
#include "VmathArray.hpp"
#include "Assertions.hpp"
#include "CellModelConstants_FK.hpp"
#include "ChunkedArray.hpp"

Array<OneD, Array<OneD, NekDouble> > cellSol;
NekChunkArray m_data;

extern NekDouble lastTime;
extern unsigned int substeps;

void init_test(int n)
{
    unsigned int n_var = 3;
    int nq = n;

    cellSol = Array<OneD, Array<OneD, NekDouble>> (n_var);
    for (unsigned int i = 0; i < n_var; ++i)
    {
	cellSol[i] = Array<OneD, NekDouble> (nq);
    }

    Vmath::Fill(nq, 0.0, cellSol[0], 1);
    Vmath::Fill(nq, 1.0, cellSol[1], 1);
    Vmath::Fill(nq, 1.0, cellSol[2], 1);

    m_data = NekChunkArray{cellSol};
}


void v_Update(NekChunkArray::ChunkUnit& chunk)
{
    PRAGMA_VECTORIZE_IVDEP
    PRAGMA_VECTORIZE_ALIGNED
    for (unsigned int i = 0; i < chunk.size; ++i)
    {
	const NekDouble& u = chunk.in0[i];
	const NekDouble& v = chunk.in1[i];
	const NekDouble& w = chunk.in2[i];
	NekDouble& u_new = chunk.out0[i];
	NekDouble& v_new = chunk.out1[i];
	NekDouble& w_new = chunk.out2[i];
	NekDouble& v_tau = chunk.gate0[i];
	NekDouble& w_tau = chunk.gate1[i];

	NekDouble J_fi, J_so, J_si, h1, h2, h3, alpha, beta;

	// Heavyside functions
	h1 = (u < u_c) ? 0.0 : 1.0;
	h2 = (u < u_v) ? 0.0 : 1.0;
	h3 = (u < u_r) ? 0.0 : 1.0;

	// w-gate
	alpha = (1-h1)/tau_w_minus;
	beta = h1/tau_w_plus;
	w_tau = 1.0 / (alpha + beta);
	w_new = alpha * w_tau;

	// v-gate
	alpha = (1-h1)/(h2*tau_v1_minus + (1-h2)*tau_v2_minus);
	beta = h1/tau_v_plus;
	v_tau = 1.0/(alpha + beta);
	v_new = alpha * v_tau;

	// J_fi
	J_fi = -v*h1*(1-u)*(u-u_fi)/tau_d;

	// J_so
	J_so = u*(1-h3)*(1-k2*v)/tau_0 + h3/tau_r;

	// J_si
	J_si = -w*(1 + tanh(k1*(u - u_csi)))/(2.0*tau_si);

	// u
	u_new = -J_fi - J_so - J_si;
    }	
}
