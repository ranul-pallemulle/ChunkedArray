#include <vector>
#include "SharedArray.hpp"
#include "VmathArray.hpp"
#include "Assertions.hpp"
#include "ChunkedArray.hpp"
#define NekDouble double

Array<OneD, Array<OneD, NekDouble> > cellSol;
NekChunkArray m_data;

extern NekDouble lastTime;
extern unsigned int substeps;

void init_test(int n)
{
    unsigned int n_var = 8;
    int nq = n;

    cellSol = Array<OneD, Array<OneD, NekDouble> > (n_var);
    for (unsigned int i = 0; i < n_var; ++i) {
	cellSol[i] = Array<OneD, NekDouble> (nq);
    }

    Vmath::Fill(nq, -84.3801107371,       cellSol[0],  1);
    Vmath::Fill(nq, 0.00171338077730188,  cellSol[1],  1);
    Vmath::Fill(nq, 0.982660523699656,    cellSol[2],  1);
    Vmath::Fill(nq, 0.989108212766685,    cellSol[3],  1);
    Vmath::Fill(nq, 0.00017948816388306,  cellSol[4],  1);
    Vmath::Fill(nq, 0.00302126301779861,  cellSol[5],  1);
    Vmath::Fill(nq, 0.999967936476325,    cellSol[6],  1);
    Vmath::Fill(nq, 0.0417603108167287,   cellSol[7],  1);

    m_data = NekChunkArray{cellSol};
}

void v_Update(NekChunkArray::ChunkUnit& chunk)
{
    const NekDouble& V_old = chunk.in0[i];
    const NekDouble& m_old = chunk.in1[i];
    const NekDouble& h_old = chunk.in2[i];
    const NekDouble& j_old = chunk.in3[i];
    const NekDouble& d_old = chunk.in4[i];
    const NekDouble& f_old = chunk.in5[i];
    const NekDouble& X_old = chunk.in6[i];
    const NekDouble& Cai_old = chunk.in7[i];

    const NekDouble E_si = 7.7 - (13.0287 * log(Cai_old / 1.0));
    const NekDouble i_si =  0.09 * d_old * f_old * (V_old - E_si);
    const NekDouble alpha_m = (0.32 * (V_old + 47.13)) /
	(1.0 - exp((-0.1) * (V_old + 47.13)));
    const NekDouble beta_m = 0.08 * exp((-V_old) / 11.0);
    const NekDouble m_gate_d_m_d_env_time = (alpha_m * (1.0 - m_old)) -
	(beta_m * m_old);
    const NekDouble beta_h = (V_old < (-40.0)) ?
	((3.56 * exp(0.079 * V_old)) + (310000.0 * exp(0.35 * V_old))) :
	(1.0 / (0.13 * (1.0 + exp((V_old + 10.66) / (-11.1)))));
    const NekDouble alpha_h = (V_old < (-40.0)) ?
	(0.135 * exp((80.0 + V_old) / (-6.8))) : 0.0;
    const NekDouble h_gate_d_h_d_env_time = (alpha_h * (1.0 - h_old)) - (beta_h * h_old);
    const NekDouble alpha_j = (V_old < (-40.0)) ?
	(((((-127140.0) * exp(0.2444 * V_old)) -
	   (3.474e-05 * exp((-0.04391) * V_old))) * (V_old + 37.78)) /
	 (1.0 + exp(0.311 * (V_old + 79.23)))) : 0.0;
    const NekDouble beta_j = (V_old < (-40.0)) ?
	((0.1212 * exp((-0.01052) * V_old)) /
	 (1.0 + exp((-0.1378) * (V_old + 40.14)))) :
	((0.3 * exp((-2.535e-07) * V_old)) / (1.0 + exp((-0.1) * (V_old + 32.0))));
    const NekDouble j_gate_d_j_d_env_time = (alpha_j * (1.0 - j_old)) - (beta_j * j_old);
    const NekDouble alpha_d = (0.095 * exp((-0.01) * (V_old - 5.0))) / (1.0 + exp((-0.072) * (V_old - 5.0)));
    const NekDouble beta_d = (0.07 * exp((-0.017) * (V_old + 44.0))) / (1.0 + exp(0.05 * (V_old + 44.0)));
    
}
