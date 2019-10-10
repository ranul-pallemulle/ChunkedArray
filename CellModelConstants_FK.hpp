#ifndef CELLMODELCONSTANTS_FK_H
#define CELLMODELCONSTANTS_FK_H

#define NekDouble double

NekDouble C_m = 1;
NekDouble V_0 = -85;
NekDouble g_fi_max = 4;
NekDouble tau_r = 33.33;
NekDouble tau_si = 29;
NekDouble tau_0 = 12.5;
NekDouble tau_v_plus = 3.33;
NekDouble tau_v1_minus = 1250;
NekDouble tau_v2_minus = 19.6;
NekDouble tau_w_plus = 870;
NekDouble tau_w_minus = 41;
NekDouble u_c = 0.13;
NekDouble u_v = 0.04;
NekDouble u_r = 0.13;
NekDouble u_fi = 0.13;
NekDouble u_csi = 0.85;
NekDouble k1 = 10;
NekDouble k2 = 0.0;
NekDouble tau_d = C_m / g_fi_max;

#endif /* CELLMODELCONSTANTS_FK_H */
