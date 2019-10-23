#include <vector>
#include "VmathArray.hpp"
#include "Assertions.hpp"
#include "CellModelConstants_FK.hpp"
#include "ChunkedArray.hpp"

NekChunkArray m_data;

extern NekDouble lastTime;
extern unsigned int substeps;

void init_test(int n)
{
    int nq = n;
    m_data = NekChunkArray(nq, {0.0, 1.0, 1.0});
}


void v_Update(NekChunkArray::ChunkUnit& chunk)
{
    NekDouble J_fi, J_so, J_si, h1, h2, h3, alpha, beta, V;
    
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
	
	V = (u - V_0)/(V_fi - V_0);

	// Heavyside functions
	h1 = (V < u_c) ? 0.0 : 1.0;
	h2 = (V < u_v) ? 0.0 : 1.0;
	h3 = (V < u_r) ? 0.0 : 1.0;

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
	J_fi = -v*h1*(1-V)*(V-u_c)/tau_d;

	// J_so
	J_so = V*(1-h3)*(1-k2*v)/tau_0 + h3/tau_r;

	// J_si
	J_si = -w*(1 + tanh(k1*(V - u_csi)))/(2.0*tau_si);

	// u
	u_new = -J_fi - J_so - J_si;
	u_new *= C_m*(V_fi - V_0);
    }	
}
