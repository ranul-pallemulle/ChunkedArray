#include <vector>
#include "SharedArray.hpp"
#include "VmathArray.hpp"
#include "Assertions.hpp"
#include "CellModelConstants_CRN98.hpp"
#include "ChunkedArray.hpp"

Array<OneD, Array<OneD, NekDouble> > cellSol;
NekChunkArray m_data;

extern NekDouble lastTime;
extern unsigned int substeps;

void init_test(int n)
{
    unsigned int n_var = 21;
    int nq = n;
    
    cellSol = Array<OneD, Array<OneD, NekDouble>> (n_var);
    for (unsigned int i = 0; i < n_var; ++i) {
	cellSol[i] = Array<OneD, NekDouble> (nq);
    }

    Vmath::Fill(nq, -81.0,      cellSol[0],  1);
    Vmath::Fill(nq, 2.908e-03,  cellSol[1],  1);
    Vmath::Fill(nq, 9.649e-01,  cellSol[2],  1);
    Vmath::Fill(nq, 9.775e-01,  cellSol[3],  1);
    Vmath::Fill(nq, 3.043e-02,  cellSol[4],  1);
    Vmath::Fill(nq, 9.992e-01,  cellSol[5],  1);
    Vmath::Fill(nq, 4.966e-03,  cellSol[6],  1);
    Vmath::Fill(nq, 9.986e-01,  cellSol[7],  1);
    Vmath::Fill(nq, 3.296e-05,  cellSol[8],  1);
    Vmath::Fill(nq, 1.869e-02,  cellSol[9],  1);
    Vmath::Fill(nq, 1.367e-04,  cellSol[10], 1);
    Vmath::Fill(nq, 9.996e-01,  cellSol[11], 1);
    Vmath::Fill(nq, 7.755e-01,  cellSol[12], 1);
    Vmath::Fill(nq, 2.35e-112,  cellSol[13], 1);
    Vmath::Fill(nq, 1.0,        cellSol[14], 1);
    Vmath::Fill(nq, 0.9992,     cellSol[15], 1);
    Vmath::Fill(nq, 1.117e+01,  cellSol[16], 1);
    Vmath::Fill(nq, 1.013e-04,  cellSol[17], 1);
    Vmath::Fill(nq, 1.39e+02,   cellSol[18], 1);
    Vmath::Fill(nq, 1.488,      cellSol[19], 1);
    Vmath::Fill(nq, 1.488,      cellSol[20], 1);

    m_data = NekChunkArray{cellSol};
}

void v_Update(NekChunkArray::ChunkUnit& chunk)
{
    // Variables
    //  0   V    membrane potential
    //  1   m    fast sodium current m gate
    //  2   h    fast sodium current h gate
    //  3   j    fast sodium current j gate
    //  4   o_a  transient outward potassium o_a gate
    //  5   o_i  transient outward potassium o_i gate
    //  6   u_a  ultra-rapid delayed rectifier K current gate
    //  7   u_i  ultra-rapid delayed rectifier K current gate
    //  8   x_r  rapid delayed rectifier K current gate
    //  9   x_s  slow delayed rectifier K current gate
    //  10  d    L_type calcium gate
    //  11  f    L-type calcium gate
    //  12  f_Ca L-type calcium gate
    //  13  u    Ca release u gate
    //  14  v    Ca release v gate
    //  15  w    Ca release w gate
    //  16  Na_i Sodium
    //  17  Ca_i Calcium
    //  18  K_i  Potassium
    //  19  Ca_rel Calcium Rel
    //  20  Ca_up  Calcium up

    PRAGMA_VECTORIZE_IVDEP
    PRAGMA_VECTORIZE_ALIGNED
    for (unsigned int i = 0; i < chunk.size; ++i)
    {
	const NekDouble& Na_i_old = chunk.in16[i];
	const NekDouble& V_old = chunk.in0[i];
	const NekDouble& m_old = chunk.in1[i];
        const NekDouble& h_old = chunk.in2[i];
        const NekDouble& j_old = chunk.in3[i];
	NekDouble& V_new = chunk.out0[i];
	NekDouble& Na_i_new = chunk.out16[i];
	NekDouble& tmp_E_Na = chunk.out14[i];
	NekDouble& tmp_I_Na = chunk.out15[i];
	NekDouble& tmp_I_b_Na = chunk.out15[i];
	
	tmp_E_Na = R*T*log(Na_o/Na_i_old)/F;

	// Sodium I_Na, chunk.out15[i] == tmp_I_Na
	tmp_I_Na = C_m * g_Na * h_old * j_old * m_old * m_old * m_old *
	    (V_old - tmp_E_Na);
	V_new = -tmp_I_Na; // chunk.out0[i] == 0 at start
	Na_i_new = -tmp_I_Na;

	// Background current, sodium, chunk.out15[i] == tmp_I_b_Na
	tmp_I_b_Na = C_m * g_b_Na * (V_old - tmp_E_Na);
	V_new = V_new - tmp_I_b_Na;
	Na_i_new = Na_i_new - tmp_I_b_Na;
    }
    PRAGMA_VECTORIZE_IVDEP
    PRAGMA_VECTORIZE_ALIGNED
    for (unsigned int i = 0; i < chunk.size; ++i)
    {
	const NekDouble& V_old = chunk.in0[i];
	const NekDouble& K_i_old = chunk.in18[i];
	const NekDouble& oa_old = chunk.in4[i];
        const NekDouble& oi_old = chunk.in5[i];
	const NekDouble& ua_old = chunk.in6[i];
        const NekDouble& ui_old = chunk.in7[i];
	const NekDouble& xr_old = chunk.in8[i];
	const NekDouble& xs_old = chunk.in9[i];
	NekDouble& V_new = chunk.out0[i];
	NekDouble& K_i_new = chunk.out18[i];
	NekDouble& tmp_V_E_k = chunk.out14[i];
	NekDouble& tmp_I_K1 = chunk.out15[i];
	NekDouble& tmp_I_to = chunk.out15[i];
	NekDouble& tmp_I_kur = chunk.out15[i];
	NekDouble& tmp_I_Kr = chunk.out15[i];
	NekDouble& tmp_I_Ks = chunk.out15[i];
	
	// V - E_K, chunk.out14[i] == tmp_V_E_k
	tmp_V_E_k = V_old - R*T*log(K_o / K_i_old)/F;

	// Potassium I_K1, chunk.out15[i] == tmp_I_K1
	tmp_I_K1 = C_m * g_K1 * tmp_V_E_k / (1.0 + exp(0.07*(80.0 + V_old)));
	V_new = V_new - tmp_I_K1;
	K_i_new = -tmp_I_K1;

	// Transient Outward K+ current, chunk.out15[i] == tmp_I_to
	tmp_I_to = C_m * g_to * oa_old * oa_old * oa_old * oi_old * tmp_V_E_k;
        V_new = V_new - tmp_I_to;
        K_i_new = K_i_new - tmp_I_to;

	// Ultrarapid Delayed rectifier K+ current, chunk.out15[i] ==
	tmp_I_kur = C_m * g_Kur_scaling * ui_old * ua_old * ua_old * ua_old *
	    tmp_V_E_k *
	    (0.005 + (0.05 / (1.0 + exp(-1.0*(-15.0 + V_old)/13.0))));
	V_new = V_new - tmp_I_kur;
        K_i_new = K_i_new - tmp_I_kur;

	// Rapid delayed outward rectifier K+ current, chunk.out15[i]
	tmp_I_Kr = C_m * g_Kr * xr_old * tmp_V_E_k /
	    (1.0 + exp((15.0 + V_old) / 22.4));
        V_new = V_new - tmp_I_Kr;
	K_i_new = K_i_new - tmp_I_Kr;

	// Slow delayed outward rectifier K+ current, chunk.out15[i]
	tmp_I_Ks = C_m * g_Ks * xs_old * xs_old * tmp_V_E_k;
	V_new = V_new - tmp_I_Ks;
	K_i_new = K_i_new - tmp_I_Ks;
    }
    PRAGMA_VECTORIZE_IVDEP
    PRAGMA_VECTORIZE_ALIGNED
    for (unsigned int i = 0; i < chunk.size; ++i)
    {
	const NekDouble& V_old = chunk.in0[i];
	const NekDouble& Ca_i_old = chunk.in17[i];
        const NekDouble& d_old = chunk.in10[i];
        const NekDouble& f_old = chunk.in11[i];
        const NekDouble& f_Ca_old = chunk.in12[i];
	NekDouble& V_new = chunk.out0[i];
	NekDouble& tmp_I_b_Ca = chunk.out1[i];
	NekDouble& tmp_I_Ca_L = chunk.out2[i];
	// Background current, calcium, chunk.out1[i] == tmp_I_b_Ca
	tmp_I_b_Ca = C_m * g_b_Ca *
	    (V_old - (0.5 * R * T * log(Ca_o / Ca_i_old) / F));
	V_new = V_new - tmp_I_b_Ca;

	// L-Type Ca2+ current, outarray[2] == tmp_I_Ca_L
	tmp_I_Ca_L = C_m * g_Ca_L * f_old * f_Ca_old * d_old * (-65.0 + V_old);
	V_new = V_new - tmp_I_Ca_L;
    }
    PRAGMA_VECTORIZE_IVDEP
    PRAGMA_VECTORIZE_ALIGNED
    for (unsigned int i = 0; i < chunk.size; ++i)
    {
	const NekDouble& V_old = chunk.in0[i];
	const NekDouble& Na_i_old = chunk.in16[i];
	const NekDouble& Ca_i_old = chunk.in17[i];
	NekDouble& V_new = chunk.out0[i];
	NekDouble& Na_i_new = chunk.out16[i];
	NekDouble& K_i_new = chunk.out18[i];
	NekDouble& tmp_f_Na_k = chunk.out14[i];
	NekDouble& tmp = chunk.out11[i];
	NekDouble& tmp_I_Na_K = chunk.out15[i];
	NekDouble& tmp_I_Na_Ca = chunk.out3[i];
	NekDouble& tmp2 = chunk.out12[i];
	NekDouble& tmp_I_p_Ca = chunk.out4[i];
	// Na-K Pump Current, outarray[14] == tmp_f_Na_k
	tmp_f_Na_k = -F/R/T * V_old;
	tmp = 0.0365 * sigma * exp(tmp_f_Na_k);
	tmp_f_Na_k = 1.0 + tmp + 0.1245 * exp(-0.1 * F/R/T * V_old);
	tmp_I_Na_K = (C_m * I_Na_K_max * K_o / (K_o + K_i)) /
	    (tmp_f_Na_k * (1.0 + pow(K_m_Na_i / Na_i_old, 1.5)));
        V_new = V_new - tmp_I_Na_K;
	Na_i_new = -3.0 * tmp_I_Na_K + Na_i_new;
	K_i_new = 2.0 * tmp_I_Na_K + K_i_new;

	// Na-Ca exchanger current, chunk.out3[i] == tmp_I_Na_Ca,
	tmp = exp((gamma_d - 1.0) * F/R/T * V_old);
	tmp_I_Na_Ca = (K_m_Na * K_m_Na * K_m_Na + Na_o * Na_o * Na_o) *
	    (K_m_Ca + Ca_o) * (1.0 + K_sat * tmp);
	tmp2 = Ca_i_old * Na_o * Na_o * Na_o * tmp;
	tmp = C_m * I_NaCa_max * (Ca_o * Na_i_old * Na_i_old * Na_i_old *
				  exp(gamma_d * F/R/T * V_old) - tmp2);
	tmp_I_Na_Ca = tmp / tmp_I_Na_Ca;
	V_new = V_new - tmp_I_Na_Ca;
	Na_i_new = -3.0 * tmp_I_Na_Ca + Na_i_new;

	// Calcium pump current, chunk.out4[i] == tmp_I_p_Ca
	tmp_I_p_Ca = C_m * I_p_Ca_max * Ca_i_old / (0.0005 + Ca_i_old);
	V_new = V_new - tmp_I_p_Ca;

	// Scale currents by capacitance
	V_new = 1.0/C_m * V_new;

	// Scale sodium and potassium by FV_i
	Na_i_new = 1.0/F/V_i * Na_i_new;
	K_i_new = 1.0/F/V_i * K_i_new;
    }
    PRAGMA_VECTORIZE_IVDEP
    PRAGMA_VECTORIZE_ALIGNED
    for (unsigned int i = 0; i < chunk.size; ++i)
    {
	const NekDouble& Ca_up_old = chunk.in20[i];
	const NekDouble& Ca_rel_old = chunk.in19[i];
	const NekDouble& Ca_i_old = chunk.in17[i];
	const NekDouble& u_old = chunk.in13[i];
        const NekDouble& v_old = chunk.in14[i];
        const NekDouble& w_old = chunk.in15[i];
	NekDouble& tmp_I_tr = chunk.out5[i];
	NekDouble& tmp_I_up_leak = chunk.out6[i];
	NekDouble& tmp_I_up = chunk.out7[i];
	NekDouble& tmp_I_rel = chunk.out8[i];
	
	// I_tr, chunk.out5[i] == tmp_I_tr
	tmp_I_tr = (Ca_up_old - Ca_rel_old) / tau_tr;

	// I_up_leak, chunk.out6[i] == tmp_I_up_leak
	tmp_I_up_leak = NSR_I_up_max/NSR_I_Ca_max * Ca_up_old;
	
	// I_up, chunk.out7[i] == tmp_I_up
        tmp_I_up = NSR_I_up_max / (1.0 + (NSR_K_up / Ca_i_old));

	// I_rel, chunk.out8[i] == tmp_I_rel
	tmp_I_rel = JSR_K_rel * v_old * w_old * u_old * u_old *
	    (Ca_rel_old - Ca_i_old);
    }
    PRAGMA_VECTORIZE_IVDEP
    PRAGMA_VECTORIZE_ALIGNED
    for (unsigned int i = 0; i < chunk.size; ++i)
    {
	const NekDouble& Ca_i_old = chunk.in17[i];
	const NekDouble& Ca_rel_old = chunk.in19[i];
	NekDouble& Ca_i_new = chunk.out17[i];
	NekDouble& Ca_up_new = chunk.out20[i];
	NekDouble& Ca_rel_new = chunk.out19[i];
	NekDouble& tmp_I_up_leak = chunk.out6[i];
	NekDouble& tmp_I_up = chunk.out7[i];
	NekDouble& tmp_I_rel = chunk.out8[i];
	NekDouble& tmp_B1 = chunk.out9[i];
	NekDouble& tmp_B2 = chunk.out10[i];
	NekDouble& tmp_I_Na_Ca = chunk.out3[i];
	NekDouble& tmp_I_p_Ca = chunk.out4[i];
	NekDouble& tmp_I_Ca_L = chunk.out2[i];
	NekDouble& tmp_I_b_Ca = chunk.out1[i];
	NekDouble& tmp_I_tr = chunk.out5[i];
	NekDouble& tmp = chunk.out11[i];
	// B1, chunk.out9[i] == tmp_B1
        tmp_B1 = (JSR_V_rel * tmp_I_rel - JSR_V_up * tmp_I_up + 
		  JSR_V_up * tmp_I_up_leak + 
		  0.5 * (2.0 * tmp_I_Na_Ca - tmp_I_p_Ca - 
			 tmp_I_Ca_L - tmp_I_b_Ca) / F) / V_i;

	// B2, outarray[10] == tmp_B2
	tmp_B2 = Cmdn_max * Km_Cmdn /
	    ((Km_Cmdn + Ca_i_old) * (Km_Cmdn + Ca_i_old));
	tmp = Trpn_max * Km_Trpn /
	    ((Km_Trpn + Ca_i_old) * (Km_Trpn + Ca_i_old));
	tmp_B2 = 1.0 + tmp_B2 + tmp;
	
	// Calcium concentration (18)
        Ca_i_new = tmp_B1 / tmp_B2;
	
	// Calcium up (21)
	Ca_up_new = -JSR_V_rel/JSR_V_up * tmp_I_tr + tmp_I_up - tmp_I_up_leak;

	// Calcium rel (20)
        tmp = tmp_I_tr - tmp_I_rel;
        Ca_rel_new = tmp / (1.0 + Csqn_max * Km_Csqn /
	     ((Km_Csqn + Ca_rel_old) * (Km_Csqn + Ca_rel_old)));
    }
    PRAGMA_VECTORIZE_IVDEP
    PRAGMA_VECTORIZE_ALIGNED
    for (unsigned int i = 0; i < chunk.size; ++i)
    {
	const NekDouble& V_old = chunk.in0[i];
	NekDouble& m_new = chunk.out1[i];
	NekDouble& m_gate = chunk.gate0[i];
	NekDouble& h_new = chunk.out2[i];
	NekDouble& h_gate = chunk.gate1[i];
	NekDouble& alpha = chunk.out11[i];
	NekDouble& beta = chunk.out12[i];
	// m
	// chunk.in0[i] == v, chunk.in1[i] == x,
	// chunk.out1[i] == x_new, chunk.gate0[i] == x_tau
	alpha = (V_old == (-47.13)) ?
	    3.2 : (0.32*(V_old + 47.13))/(1.0-exp((-0.1)*(V_old + 47.13)));
	beta = 0.08*exp(-V_old/11.0);
	m_gate = 1.0/(alpha + beta);
        m_new = alpha*m_gate;

	// h
	// chunk.in0[i] == v, chunk.in2[i] == x,
	// chunk.out2[i] == x_new, chunk.gate1[i] == x_tau
	alpha = (V_old >= -40.0) ?
	    0.0 : 0.135*exp(-(V_old + 80.0)/6.8);
	beta = (V_old >= -40.0) ?
	    1.0/(0.13*(1.0+exp(-(V_old + 10.66)/11.1)))
	    : 3.56*exp(0.079*V_old)+ 310000.0*exp(0.35*V_old);
	h_gate = 1.0/(alpha + beta);
	h_new = alpha*h_gate;
    }
    PRAGMA_VECTORIZE_IVDEP
    PRAGMA_VECTORIZE_ALIGNED
    for (unsigned int i = 0; i < chunk.size; ++i)
    {
	const NekDouble& V_old = chunk.in0[i];
	NekDouble& j_new = chunk.out3[i];
	NekDouble& j_gate = chunk.gate2[i];
	NekDouble& oa_new = chunk.out4[i];
	NekDouble& oa_gate = chunk.gate3[i];
	NekDouble& alpha = chunk.out11[i];
	NekDouble& beta = chunk.out12[i];
	// j
	// chunk.in0[i] == v, chunk.in3[i] == x,
	// chunk.out3[i] == x_new, chunk.gate2[i] == x_tau
	alpha = (V_old >= -40.0) ?
	    0.0 :
	    (-127140.0*exp(0.2444*V_old)-3.474e-05*exp(-0.04391*V_old))*((V_old+37.78)/(1.0+exp(0.311*(V_old+79.23))));
	beta  = (V_old >= -40.0) ?
	    (0.3*exp(-2.535e-07*V_old)/(1.0+exp(-0.1*(V_old+32.0))))
	    : 0.1212*exp(-0.01052*V_old)/(1.0+exp(-0.1378*(V_old+40.14)));
	j_gate = 1.0/(alpha + beta);
	j_new = alpha*(j_gate);

	// oa
	// chunk.in0[i] == v, chunk.in4[i] == x,
	// chunk.out4[i] == x_new, chunk.gate3[i] == x_tau
	alpha = 0.65/(exp(-(V_old+10.0)/8.5) + exp(-(V_old-30.0)/59.0));
	beta  = 0.65/(2.5 + exp((V_old+82.0)/17.0));
	oa_gate = 1.0/K_Q10/(alpha + beta);
	oa_new = (1.0/(1.0+exp(-(V_old+20.47)/17.54)));
    }
    PRAGMA_VECTORIZE_IVDEP
    PRAGMA_VECTORIZE_ALIGNED
    for (unsigned int i = 0; i < chunk.size; ++i)
    {
	const NekDouble& V_old = chunk.in0[i];
	NekDouble& alpha = chunk.out11[i];
	NekDouble& beta = chunk.out12[i];
	NekDouble& oi_new = chunk.out5[i];
	NekDouble& oi_gate = chunk.gate4[i];
	NekDouble& ua_new = chunk.out6[i];
	NekDouble& ua_gate = chunk.gate5[i];
	NekDouble& ui_new = chunk.out7[i];
	NekDouble& ui_gate = chunk.gate6[i];
	NekDouble& xr_new = chunk.out8[i];
	NekDouble& xr_gate = chunk.gate7[i];
	
	// oi
	// chunk.in0[i] == v, chunk.in5[i] == x,
	// chunk.out5[i] == x_new, chunk.gate4[i] == x_tau
	alpha = 1.0/(18.53 + exp((V_old+113.7)/10.95));
	beta  = 1.0/(35.56 + exp(-(V_old+1.26)/7.44));
	oi_gate = 1.0/K_Q10/(alpha + beta);
	oi_new = (1.0/(1.0+exp((V_old+43.1)/5.3)));

	// ua
	// chunk.in0[i] == v, chunk.in6[i] == x,
	// chunk.out6[i] == x_new, chunk.gate5[i] == x_tau
	alpha = 0.65/(exp(-(V_old+10.0)/8.5)+exp(-(V_old-30.0)/59.0));
	beta  = 0.65/(2.5+exp((V_old+82.0)/17.0));
	ua_gate = 1.0/K_Q10/(alpha + beta);
	ua_new = 1.0/(1.0+exp(-(V_old+30.3)/9.6));

	// ui
	// chunk.in0[i] == v, chunk.in7[i] == x,
	// chunk.out7[i] == x_new, chunk.gate6[i] == x_tau
	alpha = 1.0/(21.0 + exp(-(V_old-185.0)/28.0));
	beta  = exp((V_old-158.0)/16.0);
	ui_gate = 1.0/K_Q10/(alpha + beta);
	ui_new = 1.0/(1.0+exp((V_old-99.45)/27.48));

	// xr
	// chunk.in0[i] == v, chunk.in8[i] == x,
	// chunk.out8[i] == x_new, chunk.gate7[i] == x_tau
	alpha = 0.0003*(V_old+14.1)/(1-exp(-(V_old+14.1)/5.0));
	beta  = 7.3898e-5*(V_old-3.3328)/(exp((V_old-3.3328)/5.1237)-1.0);
	xr_gate = 1.0/(alpha + beta);
	xr_new = 1.0/(1+exp(-(V_old+14.1)/6.5));
    }
    PRAGMA_VECTORIZE_IVDEP
    PRAGMA_VECTORIZE_ALIGNED
    for (unsigned int i = 0; i < chunk.size; ++i)
    {
	const NekDouble& V_old = chunk.in0[i];
	const NekDouble& Ca_i_old = chunk.in17[i];
	NekDouble& xs_new = chunk.out9[i];
	NekDouble& xs_gate = chunk.gate8[i];
	NekDouble& d_new = chunk.out10[i];
	NekDouble& d_gate = chunk.gate9[i];
	NekDouble& f_new = chunk.out11[i];
	NekDouble& f_gate = chunk.gate10[i];
	NekDouble& f_Ca_new = chunk.out12[i];
	NekDouble& f_Ca_gate = chunk.gate11[i];
	NekDouble& alpha = chunk.out11[i];
	NekDouble& beta = chunk.out12[i];
	// xs
	// chunk.in0[i] == v, chunk.in9[i] == x,
	// chunk.out9[i] == x_new, chunk.gate8[i] == x_tau
	alpha = 4e-5*(V_old-19.9)/(1.0-exp(-(V_old-19.9)/17.0));
	beta  = 3.5e-5*(V_old-19.9)/(exp((V_old-19.9)/9.0)-1.0);
	xs_gate = 0.5/(alpha + beta);
	xs_new = 1.0/sqrt(1.0+exp(-(V_old-19.9)/12.7));

	// d
	// chunk.in0[i] == v, chunk.in10[i] == x,
	// chunk.out10[i] == x_new, chunk.gate9[i] == x_tau
	d_gate = (1-exp(-(V_old+10.0)/6.24))/(0.035*(V_old+10.0)*(1+exp(-(V_old+10.0)/6.24)));
	d_new = 1.0/(1.0 + exp(-(V_old+10)/8.0));

	// f
	// chunk.in0[i] == v, chunk.in11[i] == x,
	// chunk.out11[i] == x_new, chunk.gate10[i] == x_tau
	f_gate = 9.0/(0.0197*exp(-0.0337*0.0337*(V_old+10.0)*(V_old+10.0))+0.02);
	f_new = exp((-(V_old + 28.0)) / 6.9) / (1.0 + exp((-(V_old + 28.0)) / 6.9));

	// f_Ca
	// chunk.in0[i] == v, chunk.in12[i] == x,
	// chunk.out12[i] == x_new, chunk.gate11[i] == x_tau
        f_Ca_gate = 2.0;
        f_Ca_new = 1.0/(1.0+Ca_i_old/0.00035);
    }
    PRAGMA_VECTORIZE_IVDEP
    PRAGMA_VECTORIZE_ALIGNED
    for (unsigned int i = 0; i < chunk.size; ++i)
    {
	const NekDouble& V_old = chunk.in0[i];
	NekDouble& tmp_Fn = chunk.out15[i];
	NekDouble& tmp_I_Ca_L = chunk.out2[i];
	NekDouble& tmp_I_Na_Ca = chunk.out3[i];
	NekDouble& tmp_I_rel = chunk.out8[i];
	NekDouble& u_new = chunk.out13[i];
	NekDouble& u_gate = chunk.gate12[i];
	NekDouble& v_new = chunk.out14[i];
	NekDouble& v_gate = chunk.gate13[i];
	NekDouble& w_new = chunk.out15[i];
	NekDouble& w_gate = chunk.gate14[i];
	
	// outarray[15] == tmp_Fn
	tmp_Fn = 0.5*5e-13*tmp_I_Ca_L/F + (-0.2*5e-13)*tmp_I_Na_Ca/F;
	tmp_Fn = 1e-12*JSR_V_rel*tmp_I_rel - tmp_Fn;

	// u
	// chunk.out15[i] == v, chunk.in13[i] == x,
	// chunk.out13[i] == x_new, chunk.gate12[i] == x_tau
	u_gate  = 8.0;
        u_new = 1.0/(1.0 + exp(-(tmp_Fn - 3.4175e-13)/1.367e-15));

	// v
	// chunk.out15[i] == v, chunk.in14[i] == x,
	// chunk.out14[i] == x_new, chunk.gate13[i] == x_tau
        v_gate  = 1.91 + 2.09/(1.0+exp(-(tmp_Fn - 3.4175e-13)/13.67e-16));
        v_new = 1.0 - 1.0/(1.0 + exp(-(tmp_Fn - 6.835e-14)/13.67e-16));

	// w
	// chunk.in0[i] == v, chunk.in15[i] == x,
	// chunk.out15[i] == x_new, chunk.gate14[i] == x_tau
        w_gate = 6.0*(1.0-exp(-(V_old-7.9)/5.0))/(1.0+0.3*exp(-(V_old-7.9)/5.0))/(V_old-7.9);
        w_new = 1.0 - 1.0/(1.0 + exp(-(V_old - 40.0)/17.0));
	
    }
}
