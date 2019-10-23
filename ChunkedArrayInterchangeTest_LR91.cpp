#include <vector>
#include "VmathArray.hpp"
#include "Assertions.hpp"
#include "ChunkedArray.hpp"
#define NekDouble double

NekChunkArray m_data;

extern NekDouble lastTime;
extern unsigned int substeps;

void init_test(int n)
{
    int nq = n;
    m_data = NekChunkArray(nq, {-84.3801107371,
				0.00171338077730188,
				0.982660523699656,
				0.989108212766685,
				0.00017948816388306,
				0.00302126301779861,
				0.999967936476325,
				0.0417603108167287});
}

void v_Update(NekChunkArray::ChunkUnit& chunk)
{
    PRAGMA_VECTORIZE_ALIGNED
    PRAGMA_VECTORIZE_IVDEP
    for (unsigned int i = 0; i < chunk.size; ++i)
    {
	const NekDouble& V_old = chunk.in0[i];
	const NekDouble& m_old = chunk.in1[i];
	const NekDouble& h_old = chunk.in2[i];
	const NekDouble& j_old = chunk.in3[i];
	const NekDouble& d_old = chunk.in4[i];
	const NekDouble& f_old = chunk.in5[i];
	const NekDouble& X_old = chunk.in6[i];
	const NekDouble& Cai_old = chunk.in7[i];

	NekDouble& V_new = chunk.out0[i];
	NekDouble& m_new = chunk.out1[i];
	NekDouble& h_new = chunk.out2[i];
	NekDouble& j_new = chunk.out3[i];
	NekDouble& d_new = chunk.out4[i];
	NekDouble& f_new = chunk.out5[i];
	NekDouble& X_new = chunk.out6[i];
	NekDouble& Cai_new = chunk.out7[i];

	NekDouble& m_gate = chunk.gate0[i];
	NekDouble& h_gate = chunk.gate1[i];
	NekDouble& j_gate = chunk.gate2[i];
	NekDouble& d_gate = chunk.gate3[i];
	NekDouble& f_gate = chunk.gate4[i];
	NekDouble& X_gate = chunk.gate5[i];

	const NekDouble E_si = 7.7 - (13.0287 * log(Cai_old / 1.0));
	const NekDouble i_si =  0.09 * d_old * f_old * (V_old - E_si);
	const NekDouble alpha_m = (0.32 * (V_old + 47.13)) /
	    (1.0 - exp((-0.1) * (V_old + 47.13)));
	const NekDouble beta_m = 0.08 * exp((-V_old) / 11.0);
	// const NekDouble m_gate_d_m_d_env_time = (alpha_m * (1.0 - m_old)) -
	(beta_m * m_old);
	const NekDouble beta_h = (V_old < (-40.0)) ?
	    ((3.56 * exp(0.079 * V_old)) + (310000.0 * exp(0.35 * V_old))) :
	    (1.0 / (0.13 * (1.0 + exp((V_old + 10.66) / (-11.1)))));
	const NekDouble alpha_h = (V_old < (-40.0)) ?
	    (0.135 * exp((80.0 + V_old) / (-6.8))) : 0.0;
	// const NekDouble h_gate_d_h_d_env_time = (alpha_h * (1.0 - h_old)) - (beta_h * h_old);
	const NekDouble alpha_j = (V_old < (-40.0)) ?
	    (((((-127140.0) * exp(0.2444 * V_old)) -
	       (3.474e-05 * exp((-0.04391) * V_old))) * (V_old + 37.78)) /
	     (1.0 + exp(0.311 * (V_old + 79.23)))) : 0.0;
	const NekDouble beta_j = (V_old < (-40.0)) ?
	    ((0.1212 * exp((-0.01052) * V_old)) /
	     (1.0 + exp((-0.1378) * (V_old + 40.14)))) :
	    ((0.3 * exp((-2.535e-07) * V_old)) / (1.0 + exp((-0.1) * (V_old + 32.0))));
	// const NekDouble j_gate_d_j_d_env_time = (alpha_j * (1.0 - j_old)) - (beta_j * j_old);
	const NekDouble alpha_d = (0.095 * exp((-0.01) * (V_old - 5.0))) / (1.0 + exp((-0.072) * (V_old - 5.0)));
	const NekDouble beta_d = (0.07 * exp((-0.017) * (V_old + 44.0))) / (1.0 + exp(0.05 * (V_old + 44.0)));
	// const NekDouble d_gate_d_d_d_env_time = (alpha_d * (1.0 - d_old)) - (beta_d * d_old);
	const NekDouble alpha_f = (0.012 * exp((-0.008) * (V_old + 28.0))) / (1.0 + exp(0.15 * (V_old + 28.0)));
	const NekDouble beta_f = (0.0065 * exp((-0.02) * (V_old + 30.0))) / (1.0 + exp((-0.2) * (V_old + 30.0)));
	// const NekDouble f_gate_d_f_d_env_time = (alpha_f * (1.0 - f_old)) - (beta_f * f_old);
	const NekDouble beta_X = (0.0013 * exp((-0.06) * (V_old + 20.0))) / (1.0 + exp((-0.04) * (V_old + 20.0)));
	const NekDouble alpha_X = (0.0005 * exp(0.083 * (V_old + 50.0))) / (1.0 + exp(0.057 * (V_old + 50.0)));
	// const NekDouble X_gate_d_X_d_env_time = (alpha_X * (1.0 - X_old)) - (beta_X * X_old);
	const NekDouble Cai_d_Cai_d_env_time = (((-0.0001) / 1.0) * i_si) + (0.07 * (0.0001 - Cai_old));

	const NekDouble membrane_R = 8314.0;
	const NekDouble membrane_T = 310.0;
	const NekDouble membrane_F = 96484.6;
	const NekDouble membrane_C = 1.0;
	const NekDouble membrane_I_stim = 0.0;
	const NekDouble g_Na = 23.0;
	const NekDouble Nao = 140.0;
	const NekDouble Nai = 18.0;
	const NekDouble E_Na = ((membrane_R * membrane_T) / membrane_F) * log(Nao / Nai);
	const NekDouble i_Na = g_Na * pow(m_old, 3.0) * h_old * j_old * (V_old - E_Na);
	const NekDouble Xi = (V_old > (-100.0)) ?
	    ((2.837 * (exp(0.04 * (V_old + 77.0)) - 1.0)) /
	     ((V_old + 77.0) * exp(0.04 * (V_old + 35.0)))) : 1.0;
	const NekDouble Ko = 5.4;
	const NekDouble g_K = 0.282 * sqrt(Ko / 5.4);
	const NekDouble PR_NaK = 0.01833;
	const NekDouble Ki = 145.0;
	const NekDouble E_K = ((membrane_R * membrane_T) / membrane_F) * log((Ko + (PR_NaK * Nao)) / (Ki + (PR_NaK * Nai)));
	const NekDouble i_K = g_K * X_old * Xi * (V_old - E_K);
	const NekDouble E_K1 = ((membrane_R * membrane_T) / membrane_F) * log(Ko / Ki);
	const NekDouble beta_K1 = ((0.49124 * exp(0.08032 * ((V_old + 5.476) - E_K1))) + (1.0 * exp(0.06175 * (V_old - (E_K1 + 594.31))))) / (1.0 + exp((-0.5143) * ((V_old - E_K1) + 4.753)));
	const NekDouble alpha_K1 = 1.02 / (1.0 + exp(0.2385 * ((V_old - E_K1) - 59.215)));
	const NekDouble K1_inf = alpha_K1 / (alpha_K1 + beta_K1);
	const NekDouble g_K1 = 0.6047 * sqrt(Ko / 5.4);
	const NekDouble i_K1 = g_K1 * K1_inf * (V_old - E_K1);
	const NekDouble g_Kp = 0.0183;
	const NekDouble Kp = 1.0 / (1.0 + exp((7.488 - V_old) / 5.98));
	const NekDouble E_Kp = E_K1;
	const NekDouble i_Kp = g_Kp * Kp * (V_old - E_Kp);
	const NekDouble E_b = -59.87;
	const NekDouble g_b = 0.03921;
	const NekDouble i_b = g_b * (V_old - E_b);
	const NekDouble d_V_d_env_time = ((-1.0) / membrane_C) * (membrane_I_stim + i_Na + i_si + i_K + i_K1 + i_Kp + i_b);
	const NekDouble m_inf = alpha_m/(alpha_m + beta_m);
	const NekDouble m_tau = 1.0/(alpha_m + beta_m);
	const NekDouble h_inf = alpha_h/(alpha_h + beta_h);
	const NekDouble h_tau = 1.0/(alpha_h + beta_h);
	const NekDouble j_inf = alpha_j/(alpha_j + beta_j);
	const NekDouble j_tau = 1.0/(alpha_j + beta_j);
	const NekDouble d_inf = alpha_d/(alpha_d + beta_d);
	const NekDouble d_tau = 1.0/(alpha_d + beta_d);
	const NekDouble f_inf = alpha_f/(alpha_f + beta_f);
	const NekDouble f_tau = 1.0/(alpha_f + beta_f);
	const NekDouble X_inf = alpha_X/(alpha_X + beta_X);
	const NekDouble X_tau = 1.0/(alpha_X + beta_X);

	V_new = d_V_d_env_time;
	m_new = m_inf;
	m_gate = m_tau;
	h_new = h_inf;
	h_gate = h_tau;
	j_new = j_inf;
	j_gate = j_tau;
	d_new = d_inf;
	d_gate = d_tau;
	f_new = f_inf;
	f_gate = f_tau;
	X_new = X_inf;
	X_gate = X_tau;
	Cai_new = Cai_d_Cai_d_env_time;	
    }
}
