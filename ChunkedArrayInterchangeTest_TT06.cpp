#include <vector>
#include "VmathArray.hpp"
#include "Assertions.hpp"
#include "CellModelConstants_TT06.hpp"
#include "ChunkedArray.hpp"

NekChunkArray m_data;

extern NekDouble lastTime;
extern unsigned int substeps;

void init_test(int n)
{
    int nq = n;
    m_data = NekChunkArray(nq, {-86.709,
				0.00448,
				0.476,
				0.0087,
				0.00155,
				0.7573,
				0.7225,
				3.164e-5,
				0.8009,
				0.9778,
				0.9953,
				0.3212,
				2.235e-8,
				0.00013,
				3.715,
				0.00036,
				0.9068,
				10.355,
				138.4});
}

void v_Update(NekChunkArray::ChunkUnit& chunk)
{
    PRAGMA_VECTORIZE_ALIGNED
    PRAGMA_VECTORIZE_IVDEP
    for (unsigned int i = 0; i < chunk.size; ++i)
    {
	const NekDouble& V_old = chunk.in0[i];
	const NekDouble& Xr1_old = chunk.in1[i];
	const NekDouble& Xr2_old = chunk.in2[i];
	const NekDouble& Xs_old = chunk.in3[i];
	const NekDouble& m_old = chunk.in4[i];
	const NekDouble& h_old = chunk.in5[i];
	const NekDouble& j_old = chunk.in6[i];
	const NekDouble& d_old = chunk.in7[i];
	const NekDouble& f_old = chunk.in8[i];
	const NekDouble& f2_old = chunk.in9[i];
	const NekDouble& fCass_old = chunk.in10[i];
	const NekDouble& s_old = chunk.in11[i];
	const NekDouble& r_old = chunk.in12[i];
	const NekDouble& Ca_i_old = chunk.in13[i];
	const NekDouble& Ca_SR_old = chunk.in14[i];
	const NekDouble& Ca_ss_old = chunk.in15[i];
	const NekDouble& R_prime_old = chunk.in16[i];
	const NekDouble& Na_i_old = chunk.in17[i];
	const NekDouble& K_i_old = chunk.in18[i];

	NekDouble& V_new = chunk.out0[i];
	NekDouble& Xr1_new = chunk.out1[i];
	NekDouble& Xr2_new = chunk.out2[i];
	NekDouble& Xs_new = chunk.out3[i];
	NekDouble& m_new = chunk.out4[i];
	NekDouble& h_new = chunk.out5[i];
	NekDouble& j_new = chunk.out6[i];
	NekDouble& d_new = chunk.out7[i];
	NekDouble& f_new = chunk.out8[i];
	NekDouble& f2_new = chunk.out9[i];
	NekDouble& fCass_new = chunk.out10[i];
	NekDouble& s_new = chunk.out11[i];
	NekDouble& r_new = chunk.out12[i];
	NekDouble& Ca_i_new = chunk.out13[i];
	NekDouble& Ca_SR_new = chunk.out14[i];
	NekDouble& Ca_ss_new = chunk.out15[i];
	NekDouble& R_prime_new = chunk.out16[i];
	NekDouble& Na_i_new = chunk.out17[i];
	NekDouble& K_i_new = chunk.out18[i];

	NekDouble& Xr1_gate = chunk.gate0[i];
	NekDouble& Xr2_gate = chunk.gate1[i];
	NekDouble& Xs_gate = chunk.gate2[i];
	NekDouble& m_gate = chunk.gate3[i];
	NekDouble& h_gate = chunk.gate4[i];
	NekDouble& j_gate = chunk.gate5[i];
	NekDouble& d_gate = chunk.gate6[i];
	NekDouble& f_gate = chunk.gate7[i];
	NekDouble& f2_gate = chunk.gate8[i];
	NekDouble& fCass_gate = chunk.gate9[i];
	NekDouble& s_gate = chunk.gate10[i];
	NekDouble& r_gate = chunk.gate11[i];

	const NekDouble membrane_R = 8314.472;
	const NekDouble membrane_T = 310.0;
	const NekDouble membrane_F = 96485.3415;
	const NekDouble membrane_Cm = 0.185;
	const NekDouble membrane_V_c = 0.016404;
	const NekDouble& K_o = k_0;
	const NekDouble E_K = ((membrane_R * membrane_T) / membrane_F) * log(K_o / K_i_old);
	const NekDouble beta_K1 = ((3.0 * exp(0.0002 * ((V_old - E_K) + 100.0))) + exp(0.1 * ((V_old - E_K) - 10.0))) / (1.0 + exp((-0.5) * (V_old - E_K)));
	const NekDouble alpha_K1 = 0.1 / (1.0 + exp(0.06 * ((V_old - E_K) - 200.0)));
	const NekDouble xK1_inf = alpha_K1 / (alpha_K1 + beta_K1);
	const NekDouble g_K1 = 5.405;
	const NekDouble i_K1 = g_K1 * xK1_inf * sqrt(K_o / 5.4) * (V_old - E_K);
	const NekDouble i_to = g_to * r_old * s_old * (V_old - E_K);
	const NekDouble g_Kr = 0.153;
	const NekDouble i_Kr = g_Kr * sqrt(K_o / 5.4) * Xr1_old * Xr2_old * (V_old - E_K);
	const NekDouble Na_o = 140.0;
	const NekDouble P_kna = 0.03;
	const NekDouble E_Ks = ((membrane_R * membrane_T) / membrane_F) * log((K_o + (P_kna * Na_o)) / (K_i_old + (P_kna * Na_i_old)));
	const NekDouble i_Ks = g_Ks * pow(Xs_old, 2.0) * (V_old - E_Ks);
	const NekDouble g_CaL = 3.98e-05;
	const NekDouble Ca_o = 2.0;
	const NekDouble i_CaL = (((g_CaL * d_old * f_old * f2_old * fCass_old * 4.0 * (V_old - 15.0) * pow(membrane_F, 2.0)) / (membrane_R * membrane_T)) * ((0.25 * Ca_ss_old * exp((2.0 * (V_old - 15.0) * membrane_F) / (membrane_R * membrane_T))) - Ca_o)) / (exp((2.0 * (V_old - 15.0) * membrane_F) / (membrane_R * membrane_T)) - 1.0);
	const NekDouble K_mk = 1.0;
	const NekDouble P_NaK = 2.724;
	const NekDouble K_mNa = 40.0;
	const NekDouble i_NaK = ((((P_NaK * K_o) / (K_o + K_mk)) * Na_i_old) / (Na_i_old + K_mNa)) / (1.0 + (0.1245 * exp(((-0.1) * V_old * membrane_F) / (membrane_R * membrane_T))) + (0.0353 * exp(((-V_old) * membrane_F) / (membrane_R * membrane_T))));
	const NekDouble g_Na = 14.838;
	const NekDouble E_Na = ((membrane_R * membrane_T) / membrane_F) * log(Na_o / Na_i_old);
	const NekDouble i_Na = g_Na * pow(m_old, 3.0) * h_old * j_old * (V_old - E_Na);
	const NekDouble g_bna = 0.00029;
	const NekDouble i_b_Na = g_bna * (V_old - E_Na);
	const NekDouble alpha = 2.5;
	const NekDouble gamma_c = 0.35;
	const NekDouble K_sat = 0.1;
	const NekDouble Km_Ca = 1.38;
	const NekDouble K_NaCa = 1000.0;
	const NekDouble Km_Nai = 87.5;
	const NekDouble i_NaCa = (K_NaCa * ((exp((gamma_c * V_old * membrane_F) / (membrane_R * membrane_T)) * pow(Na_i_old, 3.0) * Ca_o) - (exp(((gamma_c - 1.0) * V_old * membrane_F) / (membrane_R * membrane_T)) * pow(Na_o, 3.0) * Ca_i_old * alpha))) / ((pow(Km_Nai, 3.0) + pow(Na_o, 3.0)) * (Km_Ca + Ca_o) * (1.0 + (K_sat * exp(((gamma_c - 1.0) * V_old * membrane_F) / (membrane_R * membrane_T)))));
	const NekDouble E_Ca = ((0.5 * membrane_R * membrane_T) / membrane_F) * log(Ca_o / Ca_i_old);
	const NekDouble g_bca = 0.000592;
	const NekDouble i_b_Ca = g_bca * (V_old - E_Ca);
	const NekDouble g_pK = 0.0146;
	const NekDouble i_p_K = (g_pK * (V_old - E_K)) / (1.0 + exp((25.0 - V_old) / 5.98));
	const NekDouble K_pCa = 0.0005;
	const NekDouble g_pCa = 0.1238;
	const NekDouble i_p_Ca = (g_pCa * Ca_i_old) / (Ca_i_old + K_pCa);
	const NekDouble i_Stim_converter = 0.0;
	const NekDouble alpha_xr1 = 450.0 / (1.0 + exp(((-45.0) - V_old) / 10.0));
	const NekDouble beta_xr1 = 6.0 / (1.0 + exp((V_old + 30.0) / 11.5));
	const NekDouble tau_xr1 = 1.0 * alpha_xr1 * beta_xr1;
	const NekDouble xr1_inf = 1.0 / (1.0 + exp(((-26.0) - V_old) / 7.0));
	const NekDouble alpha_xr2 = 3.0 / (1.0 + exp(((-60.0) - V_old) / 20.0));
	const NekDouble beta_xr2 = 1.12 / (1.0 + exp((V_old - 60.0) / 20.0));
	const NekDouble tau_xr2 = 1.0 * alpha_xr2 * beta_xr2;
	const NekDouble xr2_inf = 1.0 / (1.0 + exp((V_old + 88.0) / 24.0));
	const NekDouble beta_xs = 1.0 / (1.0 + exp((V_old - 35.0) / 15.0));
	const NekDouble alpha_xs = 1400.0 / sqrt(1.0 + exp((5.0 - V_old) / 6.0));
	const NekDouble tau_xs = (1.0 * alpha_xs * beta_xs) + 80.0;
	const NekDouble xs_inf = 1.0 / (1.0 + exp(((-5.0) - V_old) / 14.0));
	const NekDouble alpha_m = 1.0 / (1.0 + exp(((-60.0) - V_old) / 5.0));
	const NekDouble beta_m = (0.1 / (1.0 + exp((V_old + 35.0) / 5.0))) + (0.1 / (1.0 + exp((V_old - 50.0) / 200.0)));
	const NekDouble tau_m = 1.0 * alpha_m * beta_m;
	const NekDouble m_inf = 1.0 / pow(1.0 + exp(((-56.86) - V_old) / 9.03), 2.0);
	const NekDouble h_inf = 1.0 / pow(1.0 + exp((V_old + 71.55) / 7.43), 2.0);
	const NekDouble beta_h = (V_old < (-40.0)) ? ((2.7 * exp(0.079 * V_old)) + (310000.0 * exp(0.3485 * V_old))) : (0.77 / (0.13 * (1.0 + exp((V_old + 10.66) / (-11.1)))));
	const NekDouble alpha_h = (V_old < (-40.0)) ? (0.057 * exp((-(V_old + 80.0)) / 6.8)) : 0.0;
	const NekDouble tau_h = 1.0 / (alpha_h + beta_h);
	const NekDouble j_inf = 1.0 / pow(1.0 + exp((V_old + 71.55) / 7.43), 2.0);
	const NekDouble alpha_j = (V_old < (-40.0)) ? ((((((-25428.0) * exp(0.2444 * V_old)) - (6.948e-06 * exp((-0.04391) * V_old))) * (V_old + 37.78)) / 1.0) / (1.0 + exp(0.311 * (V_old + 79.23)))) : 0.0;
	const NekDouble beta_j = (V_old < (-40.0)) ? ((0.02424 * exp((-0.01052) * V_old)) / (1.0 + exp((-0.1378) * (V_old + 40.14)))) : ((0.6 * exp(0.057 * V_old)) / (1.0 + exp((-0.1) * (V_old + 32.0))));
	const NekDouble tau_j = 1.0 / (alpha_j + beta_j);
	const NekDouble alpha_d = (1.4 / (1.0 + exp(((-35.0) - V_old) / 13.0))) + 0.25;
	const NekDouble gamma_d =  1.0 / (1.0 + exp((50.0 - V_old) / 20.0));
	const NekDouble beta_d = 1.4 / (1.0 + exp((V_old + 5.0) / 5.0));
	const NekDouble tau_d = (1.0 * alpha_d * beta_d) + gamma_d;
	const NekDouble d_inf = 1.0 / (1.0 + exp(((-8.0) - V_old) / 7.5));
	const NekDouble tau_f = (1102.5 * exp((-pow(V_old + 27.0, 2.0)) / 225.0)) + (200.0 / (1.0 + exp((13.0 - V_old) / 10.0))) + (180.0 / (1.0 + exp((V_old + 30.0) / 10.0))) + 20.0;
	const NekDouble f_inf = 1.0 / (1.0 + exp((V_old + 20.0) / 7.0));
	const NekDouble f2_inf = (0.67 / (1.0 + exp((V_old + 35.0) / 7.0))) + 0.33;
	const NekDouble tau_f2 = (562.0 * exp((-pow(V_old + 27.0, 2.0)) / 240.0)) + (31.0 / (1.0 + exp((25.0 - V_old) / 10.0))) + (80.0 / (1.0 + exp((V_old + 30.0) / 10.0)));
	const NekDouble tau_fCass = (80.0 / (1.0 + pow(Ca_ss_old / 0.05, 2.0))) + 2.0;
	const NekDouble fCass_inf = (0.6 / (1.0 + pow(Ca_ss_old / 0.05, 2.0))) + 0.4;
	const NekDouble s_inf = 1.0 / (1.0 + exp((V_old + s_inf_factor) / 5.0));
	const NekDouble tau_s = (s_tau_f1 * exp((-pow(V_old + s_tau_f2, 2.0)) / s_tau_f3)) + s_tau_f4 + s_tau_f5*((5.0 / (1.0 + exp((V_old - 20.0) / 5.0))) + 3.0);
	const NekDouble r_inf = 1.0 / (1.0 + exp((20.0 - V_old) / 6.0));
	const NekDouble tau_r = (9.5 * exp((-pow(V_old + 40.0, 2.0)) / 1800.0)) + 0.8;
	const NekDouble V_rel = 0.102;
	const NekDouble k1_prime = 0.15;
	const NekDouble max_sr = 2.5;
	const NekDouble EC = 1.5;
	const NekDouble min_sr = 1.0;
	const NekDouble kcasr = max_sr - ((max_sr - min_sr) / (1.0 + pow(EC / Ca_SR_old, 2.0)));
	const NekDouble k1 = k1_prime / kcasr;
	const NekDouble k3 = 0.06;
	const NekDouble O = (k1 * pow(Ca_ss_old, 2.0) * R_prime_old) / (k3 + (k1 * pow(Ca_ss_old, 2.0)));
	const NekDouble i_rel = V_rel * O * (Ca_SR_old - Ca_ss_old);
	const NekDouble Vmax_up = 0.006375;
	const NekDouble K_up = 0.00025;
	const NekDouble i_up = Vmax_up / (1.0 + (pow(K_up, 2.0) / pow(Ca_i_old, 2.0)));
	const NekDouble V_leak = 0.00036;
	const NekDouble i_leak = V_leak * (Ca_SR_old - Ca_i_old);
	const NekDouble V_xfer = 0.0038;
	const NekDouble i_xfer = V_xfer * (Ca_ss_old - Ca_i_old);
	const NekDouble k2_prime = 0.045;
	const NekDouble k2 = k2_prime * kcasr;
	const NekDouble k4 = 0.005;
	const NekDouble Buf_c = 0.2;
	const NekDouble K_buf_c = 0.001;
	const NekDouble Ca_i_bufc = 1.0 / (1.0 + ((Buf_c * K_buf_c) / pow(Ca_i_old + K_buf_c, 2.0)));
	const NekDouble K_buf_sr = 0.3;
	const NekDouble Buf_sr = 10.0;
	const NekDouble Ca_sr_bufsr = 1.0 / (1.0 + ((Buf_sr * K_buf_sr) / pow(Ca_SR_old + K_buf_sr, 2.0)));
	const NekDouble Buf_ss = 0.4;
	const NekDouble K_buf_ss = 0.00025;
	const NekDouble Ca_ss_bufss = 1.0 / (1.0 + ((Buf_ss * K_buf_ss) / pow(Ca_ss_old + K_buf_ss, 2.0)));
	const NekDouble V_sr = 0.001094;
	const NekDouble V_ss = 5.468e-05;
	const NekDouble d_Ca_i_d_env_time = Ca_i_bufc * (((((i_leak - i_up) * V_sr) / membrane_V_c) + i_xfer) - ((1.0 * ((i_b_Ca + i_p_Ca) - (2.0 * i_NaCa)) * membrane_Cm) / (2.0 * 1.0 * membrane_V_c * membrane_F)));
	const NekDouble d_Ca_SR_d_env_time = Ca_sr_bufsr * (i_up - (i_rel + i_leak));
	const NekDouble d_Ca_ss_d_env_time = Ca_ss_bufss * (((((-1.0) * i_CaL * membrane_Cm) / (2.0 * 1.0 * V_ss * membrane_F)) + ((i_rel * V_sr) / V_ss)) - ((i_xfer * membrane_V_c) / V_ss));
	const NekDouble d_R_prime_d_env_time = ((-k2) * Ca_ss_old * R_prime_old) + (k4 * (1.0 - R_prime_old));
	const NekDouble d_Na_i_d_env_time = (((-1.0) * (i_Na + i_b_Na + (3.0 * i_NaK) + (3.0 * i_NaCa))) / (1.0 * membrane_V_c * membrane_F)) * membrane_Cm;
	const NekDouble membrane_capacitance = 1.0;
	const NekDouble i_Stim = i_Stim_converter / membrane_capacitance;
	const NekDouble d_K_i_d_env_time = (((-1.0) * ((i_K1 + i_to + i_Kr + i_Ks + i_p_K + i_Stim) - (2.0 * i_NaK))) / (1.0 * membrane_V_c * membrane_F)) * membrane_Cm;
	const NekDouble d_V_d_env_time = ((-1.0) / 1.0) * (i_K1 + i_to + i_Kr + i_Ks + i_CaL + i_NaK + i_Na + i_b_Na + i_NaCa + i_b_Ca + i_p_K + i_p_Ca + i_Stim);

	V_new = d_V_d_env_time;
	Xr1_new = xr1_inf;
	Xr1_gate = tau_xr1;
	Xr2_new = xr2_inf;
	Xr2_gate = tau_xr2;
	Xs_new = xs_inf;
	Xs_gate = tau_xs;
	m_new = m_inf;
	m_gate = tau_m;
	h_new = h_inf;
	h_gate = tau_h;
	j_new = j_inf;
	j_gate = tau_j;
	d_new = d_inf;
	d_gate= tau_d;
	f_new = f_inf;
	f_gate=  tau_f;
	f2_new = f2_inf;
	f2_gate = tau_f2;
	fCass_new = fCass_inf;
	fCass_gate = tau_fCass;
	s_new = s_inf;
	s_gate = tau_s;
	r_new = r_inf;
	r_gate = tau_r;
	Ca_i_new = d_Ca_i_d_env_time;
	Ca_SR_new = d_Ca_SR_d_env_time;
	Ca_ss_new = d_Ca_ss_d_env_time;
	R_prime_new = d_R_prime_d_env_time;
	Na_i_new = d_Na_i_d_env_time;
	K_i_new = d_K_i_d_env_time;
    }
}
