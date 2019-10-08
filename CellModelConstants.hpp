#ifndef CELLMODELCONSTANTS_H
#define CELLMODELCONSTANTS_H

#include <cmath>

#define NekDouble double

NekDouble C_m = 100;      // picoF
NekDouble g_Na = 7.8;     // nanoS_per_picoF
NekDouble g_K1 = 0.09;    // nanoS_per_picoF
NekDouble g_Kr = 0.029411765;
NekDouble g_Ks = 0.12941176;
NekDouble g_b_Na = 0.0006744375;
NekDouble g_b_Ca = 0.001131;
NekDouble R = 8.3143;
NekDouble T = 310.0;
NekDouble F = 96.4867;
NekDouble Na_o = 140.0;   // millimolar
NekDouble K_o = 5.4;      // millimolar
NekDouble sigma = 1.0/7.0*(exp(Na_o/67.3)-1);
NekDouble K_i = 1.5;
NekDouble K_m_Na_i = 10.0;
NekDouble I_Na_K_max = 0.59933874;
NekDouble I_NaCa_max = 1600.0;
NekDouble gamma_d = 0.35;
NekDouble Ca_o = 1.8;
NekDouble K_m_Na = 87.5;
NekDouble K_m_Ca = 1.38;
NekDouble K_sat = 0.1;
NekDouble I_p_Ca_max = 0.275;
NekDouble Trpn_max = 0.07;
NekDouble Km_Trpn = 0.0005;
NekDouble Cmdn_max = 0.05;
NekDouble Csqn_max = 10.0;
NekDouble Km_Cmdn = 0.00238;
NekDouble Km_Csqn = 0.8;
NekDouble NSR_I_up_max = 0.005;
NekDouble NSR_I_Ca_max = 15.0;
NekDouble NSR_K_up = 0.00092;
NekDouble JSR_K_rel = 30.0;
NekDouble JSR_V_cell = 20100.0;
NekDouble JSR_V_rel = 0.0048 * JSR_V_cell;
NekDouble JSR_V_up = 0.0552 * JSR_V_cell;
NekDouble tau_tr = 180.0;
NekDouble K_Q10 = 3.0;
NekDouble V_i = 0.68*JSR_V_cell;
// eOriginal
NekDouble g_to = 0.1652;  // nanoS_per_picoF
NekDouble g_Kur_scaling = 1.0;
NekDouble g_Ca_L = 0.12375;

#endif /* CELLMODELCONSTANTS_H */
