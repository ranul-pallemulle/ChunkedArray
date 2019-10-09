#include <vector>
#include "SharedArray.hpp"
#include "VmathArray.hpp"
#include "Assertions.hpp"
#include "CellModelConstants.hpp"

#define NUM_IN_OUT_VARS 21
#define NUM_GATE_VARS 15
#define HAS_CONCENTRATIONS

#include "ChunkedArray.hpp"

int nq;
unsigned int n_var;
Array<OneD, Array<OneD, NekDouble> > cellSol;
Array<OneD, Array<OneD, NekDouble> > wsp;
Array<OneD, Array<OneD, NekDouble> > gates_tau;

using NekChunkArray = Nektar::ChunkArray<double, 512, 8>;
NekChunkArray m_data;

extern NekDouble lastTime;
extern unsigned int substeps;
extern NekDouble finTime;

void init_test(int n)
{
    n_var = 21;
    nq = n;
    
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
    NekDouble alpha, beta;
	
    PRAGMA_VECTORIZE_IVDEP
    for (unsigned int i = 0; i < chunk.size; ++i)
    {
	// E_Na, chunk.out14[i] == tmp_E_Na
	chunk.out14[i] = R*T*log(Na_o/chunk.in16[i])/F;
    }
    PRAGMA_VECTORIZE_IVDEP
    for (unsigned int i = 0; i < chunk.size; ++i)
    {
	// Sodium I_Na, chunk.out15[i] == tmp_I_Na
	chunk.out15[i] = C_m * g_Na *
	    chunk.in2[i] * chunk.in3[i] *
	    chunk.in1[i] * chunk.in1[i] * chunk.in1[i] *
	    (chunk.in0[i] - chunk.out14[i]);
	chunk.out0[i] = -chunk.out15[i]; // chunk.out0[i] == 0 at start
	chunk.out16[i] = -chunk.out15[i];
    }
    PRAGMA_VECTORIZE_IVDEP
    for (unsigned int i = 0; i < chunk.size; ++i)
    {
	// Background current, sodium, chunk.out15[i] == tmp_I_b_Na
	chunk.out15[i] = C_m * g_b_Na *
	    (chunk.in0[i] - chunk.out14[i]);
	chunk.out0[i] = chunk.out0[i] - chunk.out15[i];
	chunk.out16[i] = chunk.out16[i] - chunk.out15[i];
    }
    PRAGMA_VECTORIZE_IVDEP
    for (unsigned int i = 0; i < chunk.size; ++i)
    {
	// V - E_K, chunk.out14[i] == tmp_V_E_k
	chunk.out14[i] = chunk.in0[i] -
	    R*T*log(K_o / chunk.in18[i])/F;

	// Potassium I_K1, chunk.out15[i] == tmp_I_K1
	chunk.out15[i] = C_m * g_K1 *
	    chunk.out14[i] / (1.0 + exp(0.07*(80.0 + chunk.in0[i])));
	chunk.out0[i] = chunk.out0[i] - chunk.out15[i];
	chunk.out18[i] = -chunk.out15[i];
    }
    PRAGMA_VECTORIZE_IVDEP
    for (unsigned int i = 0; i < chunk.size; ++i)
    {
	// Transient Outward K+ current, chunk.out15[i] == tmp_I_to
	chunk.out15[i] = C_m * g_to *
	    chunk.in4[i] * chunk.in4[i] * chunk.in4[i] *
	    chunk.in5[i] * chunk.out14[i];
	chunk.out0[i] = chunk.out0[i] - chunk.out15[i];
	chunk.out18[i] = chunk.out18[i] - chunk.out15[i];
    }
    PRAGMA_VECTORIZE_IVDEP
    for (unsigned int i = 0; i < chunk.size; ++i)
    {
	// Ultrarapid Delayed rectifier K+ current, chunk.out15[i] ==
	chunk.out15[i] = C_m * g_Kur_scaling * chunk.in7[i] *
	    chunk.in6[i] * chunk.in6[i] * chunk.in6[i] *
	    chunk.out14[i] * (0.005 +
			      (0.05 / (1.0 + exp(-1.0*(-15.0 + chunk.in0[i])/13.0))));
	chunk.out0[i] = chunk.out0[i] - chunk.out15[i];
	chunk.out18[i] = chunk.out18[i] - chunk.out15[i];
    }
    PRAGMA_VECTORIZE_IVDEP
    for (unsigned int i = 0; i < chunk.size; ++i)
    {
	// Rapid delayed outtward rectifier K+ current, chunk.out15[i]
	chunk.out15[i] = C_m * g_Kr * chunk.in8[i] *
	    chunk.out14[i] / (1.0 + exp((15.0 + chunk.in0[i]) / 22.4));
	chunk.out0[i] = chunk.out0[i] - chunk.out15[i];
	chunk.out18[i] = chunk.out18[i] - chunk.out15[i];
    }
    PRAGMA_VECTORIZE_IVDEP
    for (unsigned int i = 0; i < chunk.size; ++i)
    {
	// Slow delayed outtward rectifier K+ current, chunk.out15[i]
	chunk.out15[i] = C_m * g_Ks *
	    chunk.in9[i] * chunk.in9[i] * chunk.out14[i];
	chunk.out0[i] = chunk.out0[i] - chunk.out15[i];
	chunk.out18[i] = chunk.out18[i] - chunk.out15[i];
    }
    PRAGMA_VECTORIZE_IVDEP
    for (unsigned int i = 0; i < chunk.size; ++i)
    {
	// Background current, calcium, chunk.out1[i] == tmp_I_b_Ca
	chunk.out1[i] = C_m * g_b_Ca * (chunk.in0[i] - 
					(0.5 * R * T * log(Ca_o / chunk.in17[i]) / F));
	chunk.out0[i] = chunk.out0[i] - chunk.out1[i];
    }
    PRAGMA_VECTORIZE_IVDEP
    for (unsigned int i = 0; i < chunk.size; ++i)
    {
	// L-Type Ca2+ current, outtarray[2] == tmp_I_Ca_L
	chunk.out2[i] = C_m * g_Ca_L * chunk.in11[i] * chunk.in12[i] *
	    chunk.in10[i] * (-65.0 + chunk.in0[i]);
	chunk.out0[i] = chunk.out0[i] - chunk.out2[i];
    }
    PRAGMA_VECTORIZE_IVDEP
    for (unsigned int i = 0; i < chunk.size; ++i)
    {
	// Na-K Pump Current, outtarray[14] == tmp_f_Na_k
	chunk.out14[i] = -F/R/T * chunk.in0[i];
	chunk.out11[i] = 0.0365 * sigma * exp(chunk.out14[i]);
	chunk.out14[i] = 1.0 + chunk.out11[i] + 
	    0.1245 * exp(-0.1 * F/R/T * chunk.in0[i]);
	chunk.out15[i] = (C_m * I_Na_K_max * K_o / (K_o + K_i)) /
	    (chunk.out14[i] * (1.0 + pow(K_m_Na_i / chunk.in16[i], 1.5)));
	chunk.out0[i] = chunk.out0[i] - chunk.out15[i];
	chunk.out16[i] = -3.0 * chunk.out15[i] + chunk.out16[i];
	chunk.out18[i] = 2.0 * chunk.out15[i] + chunk.out18[i];
    }
    PRAGMA_VECTORIZE_IVDEP
    for (unsigned int i = 0; i < chunk.size; ++i)
    {
	// Na-Ca exchanger current, chunk.out3[i] == tmp_I_Na_Ca,
	chunk.out11[i] = exp((gamma_d - 1.0) * F/R/T * chunk.in0[i]);
	chunk.out3[i] = (K_m_Na * K_m_Na * K_m_Na + Na_o * Na_o * Na_o) *
	    (K_m_Ca + Ca_o) * 
	    (1.0 + K_sat * chunk.out11[i]);
	chunk.out12[i] = chunk.in17[i] * Na_o * Na_o * Na_o * chunk.out11[i];
	chunk.out11[i] = C_m * I_NaCa_max *
	    (Ca_o * chunk.in16[i] * chunk.in16[i] *
	     chunk.in16[i] * exp(gamma_d * F/R/T * chunk.in0[i])
	     - chunk.out12[i]);
	chunk.out3[i] = chunk.out11[i] / chunk.out3[i];
	chunk.out0[i] = chunk.out0[i] - chunk.out3[i];
	chunk.out16[i] = -3.0 * chunk.out3[i] + chunk.out16[i];
    }
    PRAGMA_VECTORIZE_IVDEP
    for (unsigned int i = 0; i < chunk.size; ++i)
    {
	// Calcium pump current, chunk.out4[i] == tmp_I_p_Ca
	chunk.out4[i] = C_m * I_p_Ca_max *
	    chunk.in17[i] / (0.0005 + chunk.in17[i]);
	chunk.out0[i] = chunk.out0[i] - chunk.out4[i];

	// Scale currents by capacitance
	chunk.out0[i] = 1.0/C_m * chunk.out0[i];

	// Scale sodium and potassium by FV_i
	chunk.out16[i] = 1.0/F/V_i * chunk.out16[i];
	chunk.out18[i] = 1.0/F/V_i * chunk.out18[i];

	// I_tr, chunk.out5[i] == tmp_I_tr
	chunk.out5[i] = (chunk.in20[i] - chunk.in19[i]) / tau_tr;

	// I_up_leak, chunk.out6[i] == tmp_I_up_leak
	chunk.out6[i] = NSR_I_up_max/NSR_I_Ca_max * chunk.in20[i];
    }
    PRAGMA_VECTORIZE_IVDEP
    for (unsigned int i = 0; i < chunk.size; ++i)
    {
	// I_up, chunk.out7[i] == tmp_I_up
	chunk.out7[i] = NSR_I_up_max / 
	    (1.0 + (NSR_K_up / chunk.in17[i]));

	// I_rel, chunk.out8[i] == tmp_I_rel
	chunk.out8[i] = JSR_K_rel * chunk.in14[i] * chunk.in15[i] *
	    chunk.in13[i] * chunk.in13[i] *
	    (chunk.in19[i] - chunk.in17[i]);

	// B1, chunk.out9[i] == tmp_B1
	chunk.out9[i] = (JSR_V_rel * chunk.out8[i] -
			 JSR_V_up * chunk.out7[i] + 
			 JSR_V_up * chunk.out6[i] + 
			 0.5 * (2.0 * chunk.out3[i] - chunk.out4[i] - 
				chunk.out2[i] - chunk.out1[i]) / F) / V_i;
    }
    PRAGMA_VECTORIZE_IVDEP
    for (unsigned int i = 0; i < chunk.size; ++i)
    {
	// B2, outarray[10] == tmp_B2
	chunk.out10[i] = Cmdn_max * Km_Cmdn /
	    ((Km_Cmdn + chunk.in17[i]) * (Km_Cmdn + chunk.in17[i]));
	chunk.out11[i] = Trpn_max * Km_Trpn /
	    ((Km_Trpn + chunk.in17[i]) * (Km_Trpn + chunk.in17[i]));
	chunk.out10[i] = 1.0 + chunk.out10[i] + chunk.out11[i];
	
	// Calcium concentration (18)
	chunk.out17[i] = chunk.out9[i] / chunk.out10[i];
	
	// Calcium up (21)
	chunk.out20[i] = -JSR_V_rel/JSR_V_up * chunk.out5[i] +
	    chunk.out7[i] - chunk.out6[i];

	// Calcium rel (20)
	chunk.out11[i] = chunk.out5[i] - chunk.out8[i];
	chunk.out19[i] = chunk.out11[i] / 
	    (1.0 + Csqn_max * Km_Csqn /
	     ((Km_Csqn + chunk.in19[i]) * (Km_Csqn + chunk.in19[i])));
    }
#pragma omp simd private(alpha,beta)
    for (unsigned int i = 0; i < chunk.size; ++i)
    {
	// m
	// chunk.in0[i] == v, chunk.in1[i] == x,
	// chunk.out1[i] == x_new, chunk.gate0[i] == x_tau
	alpha = (chunk.in0[i] == (-47.13)) ?
	    3.2 : (0.32*(chunk.in0[i] + 47.13))/(1.0-exp((-0.1)*(chunk.in0[i] + 47.13)));
	beta = 0.08*exp(-chunk.in0[i]/11.0);
	chunk.gate0[i] = 1.0/(alpha + beta);
	chunk.out1[i] = alpha*chunk.gate0[i];
    }
#pragma omp simd private(alpha,beta)
    for (unsigned int i = 0; i < chunk.size; ++i)
    {
	// h
	// chunk.in0[i] == v, chunk.in2[i] == x,
	// chunk.out2[i] == x_new, chunk.gate1[i] == x_tau
	alpha = (chunk.in0[i] >= -40.0) ?
	    0.0 : 0.135*exp(-(chunk.in0[i] + 80.0)/6.8);
	beta = (chunk.in0[i] >= -40.0) ?
	    1.0/(0.13*(1.0+exp(-(chunk.in0[i] + 10.66)/11.1)))
	    : 3.56*exp(0.079*chunk.in0[i])+ 310000.0*exp(0.35*chunk.in0[i]);
	chunk.gate1[i] = 1.0/(alpha + beta);
	chunk.out2[i] = alpha*chunk.gate1[i];
    }
#pragma omp simd private(alpha,beta)
    for (unsigned int i = 0; i < chunk.size; ++i)
    {
	// j
	// chunk.in0[i] == v, chunk.in3[i] == x,
	// chunk.out3[i] == x_new, chunk.gate2[i] == x_tau
	alpha = (chunk.in0[i] >= -40.0) ?
	    0.0 :
	    (-127140.0*exp(0.2444*chunk.in0[i])-3.474e-05*exp(-0.04391*chunk.in0[i]))*((chunk.in0[i]+37.78)/(1.0+exp(0.311*(chunk.in0[i]+79.23))));
	beta  = (chunk.in0[i] >= -40.0) ?
	    (0.3*exp(-2.535e-07*chunk.in0[i])/(1.0+exp(-0.1*(chunk.in0[i]+32.0))))
	    : 0.1212*exp(-0.01052*chunk.in0[i])/(1.0+exp(-0.1378*(chunk.in0[i]+40.14)));
	chunk.gate2[i] = 1.0/(alpha + beta);
	chunk.out3[i] = alpha*(chunk.gate2[i]);
    }
#pragma omp simd private(alpha,beta)
    for (unsigned int i = 0; i < chunk.size; ++i)
    {
	// oa
	// chunk.in0[i] == v, chunk.in4[i] == x,
	// chunk.out4[i] == x_new, chunk.gate3[i] == x_tau
	alpha = 0.65/(exp(-(chunk.in0[i]+10.0)/8.5) + exp(-(chunk.in0[i]-30.0)/59.0));
	beta  = 0.65/(2.5 + exp((chunk.in0[i]+82.0)/17.0));
	chunk.gate3[i] = 1.0/K_Q10/(alpha + beta);
	chunk.out4[i] = (1.0/(1.0+exp(-(chunk.in0[i]+20.47)/17.54)));
    }
#pragma omp simd private(alpha,beta)
    for (unsigned int i = 0; i < chunk.size; ++i)
    {
	// oi
	// chunk.in0[i] == v, chunk.in5[i] == x,
	// chunk.out5[i] == x_new, chunk.gate4[i] == x_tau
	alpha = 1.0/(18.53 + exp((chunk.in0[i]+113.7)/10.95));
	beta  = 1.0/(35.56 + exp(-(chunk.in0[i]+1.26)/7.44));
	chunk.gate4[i] = 1.0/K_Q10/(alpha + beta);
	chunk.out5[i] = (1.0/(1.0+exp((chunk.in0[i]+43.1)/5.3)));
    }
#pragma omp simd private(alpha,beta)
    for (unsigned int i = 0; i < chunk.size; ++i)
    {
	// ua
	// chunk.in0[i] == v, chunk.in6[i] == x,
	// chunk.out6[i] == x_new, chunk.gate5[i] == x_tau
	alpha = 0.65/(exp(-(chunk.in0[i]+10.0)/8.5)+exp(-(chunk.in0[i]-30.0)/59.0));
	beta  = 0.65/(2.5+exp((chunk.in0[i]+82.0)/17.0));
	chunk.gate5[i] = 1.0/K_Q10/(alpha + beta);
	chunk.out6[i] = 1.0/(1.0+exp(-(chunk.in0[i]+30.3)/9.6));
    }
#pragma omp simd private(alpha,beta)
    for (unsigned int i = 0; i < chunk.size; ++i)
    {
	// ui
	// chunk.in0[i] == v, chunk.in7[i] == x,
	// chunk.out7[i] == x_new, chunk.gate6[i] == x_tau
	alpha = 1.0/(21.0 + exp(-(chunk.in0[i]-185.0)/28.0));
	beta  = exp((chunk.in0[i]-158.0)/16.0);
	chunk.gate6[i] = 1.0/K_Q10/(alpha + beta);
	chunk.out7[i] = 1.0/(1.0+exp((chunk.in0[i]-99.45)/27.48));
    }
#pragma omp simd private(alpha,beta)
    for (unsigned int i = 0; i < chunk.size; ++i)
    {
	// xr
	// chunk.in0[i] == v, chunk.in8[i] == x,
	// chunk.out8[i] == x_new, chunk.gate7[i] == x_tau
	alpha = 0.0003*(chunk.in0[i]+14.1)/(1-exp(-(chunk.in0[i]+14.1)/5.0));
	beta  = 7.3898e-5*(chunk.in0[i]-3.3328)/(exp((chunk.in0[i]-3.3328)/5.1237)-1.0);
	chunk.gate7[i] = 1.0/(alpha + beta);
	chunk.out8[i] = 1.0/(1+exp(-(chunk.in0[i]+14.1)/6.5));
    }
#pragma omp simd private(alpha,beta)
    for (unsigned int i = 0; i < chunk.size; ++i)
    {
	// xs
	// chunk.in0[i] == v, chunk.in9[i] == x,
	// chunk.out9[i] == x_new, chunk.gate8[i] == x_tau
	alpha = 4e-5*(chunk.in0[i]-19.9)/(1.0-exp(-(chunk.in0[i]-19.9)/17.0));
	beta  = 3.5e-5*(chunk.in0[i]-19.9)/(exp((chunk.in0[i]-19.9)/9.0)-1.0);
	chunk.gate8[i] = 0.5/(alpha + beta);
	chunk.out9[i] = 1.0/sqrt(1.0+exp(-(chunk.in0[i]-19.9)/12.7));
    }
    PRAGMA_VECTORIZE_IVDEP
    for (unsigned int i = 0; i < chunk.size; ++i)
    {
	// d
	// chunk.in0[i] == v, chunk.in10[i] == x,
	// chunk.out10[i] == x_new, chunk.gate9[i] == x_tau
	chunk.gate9[i] = (1-exp(-(chunk.in0[i]+10.0)/6.24))/(0.035*(chunk.in0[i]+10.0)*(1+exp(-(chunk.in0[i]+10.0)/6.24)));
	chunk.out10[i] = 1.0/(1.0 + exp(-(chunk.in0[i]+10)/8.0));
    }
    PRAGMA_VECTORIZE_IVDEP
    for (unsigned int i = 0; i < chunk.size; ++i)
    {
	// f
	// chunk.in0[i] == v, chunk.in11[i] == x,
	// chunk.out11[i] == x_new, chunk.gate10[i] == x_tau
	chunk.gate10[i] = 9.0/(0.0197*exp(-0.0337*0.0337*(chunk.in0[i]+10.0)*(chunk.in0[i]+10.0))+0.02);
	chunk.out11[i] = exp((-(chunk.in0[i] + 28.0)) / 6.9) / (1.0 + exp((-(chunk.in0[i] + 28.0)) / 6.9));
    }
    PRAGMA_VECTORIZE_IVDEP
    for (unsigned int i = 0; i < chunk.size; ++i)
    {
	// f_Ca
	// chunk.in0[i] == v, chunk.in12[i] == x,
	// chunk.out12[i] == x_new, chunk.gate11[i] == x_tau
	chunk.gate11[i]  = 2.0;
	chunk.out12[i] = 1.0/(1.0+chunk.in17[i]/0.00035);
    }
    PRAGMA_VECTORIZE_IVDEP
    for (unsigned int i = 0; i < chunk.size; ++i)
    {
	// outarray[15] == tmp_Fn
	chunk.out15[i] = 0.5*5e-13*chunk.out2[i]/F + (-0.2*5e-13)*chunk.out3[i]/F;
	chunk.out15[i] = 1e-12*JSR_V_rel*chunk.out8[i] - chunk.out15[i];
    }
    PRAGMA_VECTORIZE_IVDEP
    for (unsigned int i = 0; i < chunk.size; ++i)
    {
	// u
	// chunk.out15[i] == v, chunk.in13[i] == x,
	// chunk.out13[i] == x_new, chunk.gate12[i] == x_tau
	chunk.gate12[i]  = 8.0;
	chunk.out13[i] = 1.0/(1.0 + exp(-(chunk.out15[i] - 3.4175e-13)/1.367e-15));
    }
    PRAGMA_VECTORIZE_IVDEP
    for (unsigned int i = 0; i < chunk.size; ++i)
    {
	// v
	// chunk.out15[i] == v, chunk.in14[i] == x,
	// chunk.out14[i] == x_new, chunk.gate13[i] == x_tau
	chunk.gate13[i]  = 1.91 + 2.09/(1.0+exp(-(chunk.out15[i] - 3.4175e-13)/13.67e-16));
	chunk.out14[i] = 1.0 - 1.0/(1.0 + exp(-(chunk.out15[i] - 6.835e-14)/13.67e-16));
    }
    PRAGMA_VECTORIZE_IVDEP
    for (unsigned int i = 0; i < chunk.size; ++i)
    {
	// w
	// chunk.in0[i] == v, chunk.in15[i] == x,
	// chunk.out15[i] == x_new, chunk.gate14[i] == x_tau
	chunk.gate14[i] = 6.0*(1.0-exp(-(chunk.in0[i]-7.9)/5.0))/(1.0+0.3*exp(-(chunk.in0[i]-7.9)/5.0))/(chunk.in0[i]-7.9);
	chunk.out15[i] = 1.0 - 1.0/(1.0 + exp(-(chunk.in0[i] - 40.0)/17.0));
    }
}


void TimeIntegrate(const Array<OneD, const Array<OneD, NekDouble> >& inarray,
		         Array<OneD,       Array<OneD, NekDouble> >& outarray,
		   const NekDouble time)
{
    // Copy in voltage
    m_data.fromInArray(inarray[0]);
    NekDouble delta_t = (time - lastTime)/substeps;

    #pragma omp parallel for
    for (unsigned int k = 0; k < m_data.num_elements(); ++k)
    {
	NekChunkArray::ChunkUnit& chunk = m_data[k];
	for (unsigned int i = 0; i < substeps - 1; ++i)
	{
	    v_Update(chunk);

	    PRAGMA_VECTORIZE_IVDEP
	    for (unsigned int j = 0; j < chunk.size; ++j)
	    {
		// Forward Euler on voltage
		chunk.in0[j] = delta_t * chunk.out0[j] + chunk.in0[j];

		// Forward Euler on concentrations
#define BOOST_PP_LOCAL_MACRO(n)						\
		chunk.in##n[j] = delta_t * chunk.out##n[j] + chunk.in##n[j];
#define BOOST_PP_LOCAL_LIMITS (CONC_START_IDX, CONC_END_IDX)
#include BOOST_PP_LOCAL_ITERATE()

		// Rush-Larsen integration on gating variables
#define BOOST_PP_LOCAL_MACRO(n)					\
		chunk.gate##n[j] = -delta_t / chunk.gate##n[j]; \
								\
		chunk.gate##n[j] = exp(chunk.gate##n[j]);	\
								\
		APPEND(chunk.in, BOOST_PP_ADD(n,1)[j]) =	\
		    APPEND(chunk.in, BOOST_PP_ADD(n,1)[j]) -	\
		    APPEND(chunk.out, BOOST_PP_ADD(n,1)[j]);	\
								\
		APPEND(chunk.in, BOOST_PP_ADD(n,1)[j]) =	\
		    APPEND(chunk.in, BOOST_PP_ADD(n,1)[j]) *	\
		    APPEND(chunk.gate, n[j]) +			\
		    APPEND(chunk.out, BOOST_PP_ADD(n,1)[j]);
		
#define BOOST_PP_LOCAL_LIMITS (0, NUM_GATE_VARS - 1)
#include BOOST_PP_LOCAL_ITERATE()
	    }
	}
	
	v_Update(chunk);

	PRAGMA_VECTORIZE_IVDEP
	for (unsigned int j = 0; j < chunk.size; ++j)
	{
	    // Forward Euler on concentrations
#define BOOST_PP_LOCAL_MACRO(n)						\
	    chunk.in##n[j] = delta_t * chunk.out##n[j] + chunk.in##n[j];
#define BOOST_PP_LOCAL_LIMITS (CONC_START_IDX, CONC_END_IDX)
#include BOOST_PP_LOCAL_ITERATE()

	    // Rush-Larsen integration on gating variables
#define BOOST_PP_LOCAL_MACRO(n)					\
	    chunk.gate##n[j] = -delta_t / chunk.gate##n[j];	\
								\
	    chunk.gate##n[j] = exp(chunk.gate##n[j]);		\
								\
	    APPEND(chunk.in, BOOST_PP_ADD(n,1)[j]) =		\
		APPEND(chunk.in, BOOST_PP_ADD(n,1)[j]) -	\
		APPEND(chunk.out, BOOST_PP_ADD(n,1)[j]);	\
								\
	    APPEND(chunk.in, BOOST_PP_ADD(n,1)[j]) =		\
		APPEND(chunk.in, BOOST_PP_ADD(n,1)[j]) *	\
		APPEND(chunk.gate, n[j]) +			\
		APPEND(chunk.out, BOOST_PP_ADD(n,1)[j]);
		
#define BOOST_PP_LOCAL_LIMITS (0, NUM_GATE_VARS - 1)
#include BOOST_PP_LOCAL_ITERATE()
	}
    }
    // Copy out new voltage
    m_data.toOutArray(outarray[0]);
    lastTime = time;
}
