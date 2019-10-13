#include <vector>
#include "SharedArray.hpp"
#include "VmathArray.hpp"
#include "Assertions.hpp"
#define NekDouble double

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
    gates.push_back(3);
    gates.push_back(4);
    gates.push_back(5);
    gates.push_back(6);
    concentrations.push_back(7);
    n_var = 8;
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

    Vmath::Fill(nq, -84.3801107371,       cellSol[0],  1);
    Vmath::Fill(nq, 0.00171338077730188,  cellSol[1],  1);
    Vmath::Fill(nq, 0.982660523699656,    cellSol[2],  1);
    Vmath::Fill(nq, 0.989108212766685,    cellSol[3],  1);
    Vmath::Fill(nq, 0.00017948816388306,  cellSol[4],  1);
    Vmath::Fill(nq, 0.00302126301779861,  cellSol[5],  1);
    Vmath::Fill(nq, 0.999967936476325,    cellSol[6],  1);
    Vmath::Fill(nq, 0.0417603108167287,   cellSol[7],  1);
}

void v_Update(const Array<OneD, const Array<OneD, NekDouble> >& inarray,
	            Array<OneD,       Array<OneD, NekDouble> >& outarray,
	      const NekDouble)
{
    for (int i = 0; i < nq; ++i)
    {
	// Inputs:
	// Time units: millisecond
	NekDouble var_chaste_interface__membrane__V = inarray[0][i];
	// Units: millivolt; Initial value: -84.3801107371
	NekDouble var_chaste_interface__fast_sodium_current_m_gate__m = inarray[1][i];
	// Units: dimensionless; Initial value: 0.00171338077730188
	NekDouble var_chaste_interface__fast_sodium_current_h_gate__h = inarray[2][i];
	// Units: dimensionless; Initial value: 0.982660523699656
	NekDouble var_chaste_interface__fast_sodium_current_j_gate__j = inarray[3][i];
	// Units: dimensionless; Initial value: 0.989108212766685
	NekDouble var_chaste_interface__slow_inward_current_d_gate__d = inarray[4][i];
	// Units: dimensionless; Initial value: 0.00302126301779861
	NekDouble var_chaste_interface__slow_inward_current_f_gate__f = inarray[5][i];
	// Units: dimensionless; Initial value: 0.999967936476325
	NekDouble var_chaste_interface__time_dependent_potassium_current_X_gate__X = inarray[6][i];
	// Units: dimensionless; Initial value: 0.0417603108167287
	NekDouble var_chaste_interface__intracellular_calcium_concentration__Cai = inarray[7][i];
	// Units: millimolar; Initial value: 0.00017948816388306

            
	// Mathematics
	NekDouble d_dt_chaste_interface__membrane__V;
	const NekDouble var_fast_sodium_current__j = var_chaste_interface__fast_sodium_current_j_gate__j; // dimensionless
	const NekDouble var_fast_sodium_current__h = var_chaste_interface__fast_sodium_current_h_gate__h; // dimensionless
	const NekDouble var_fast_sodium_current__m = var_chaste_interface__fast_sodium_current_m_gate__m; // dimensionless
	const NekDouble var_fast_sodium_current__V = var_chaste_interface__membrane__V; // millivolt
	const NekDouble var_slow_inward_current__d = var_chaste_interface__slow_inward_current_d_gate__d; // dimensionless
	const NekDouble var_slow_inward_current__f = var_chaste_interface__slow_inward_current_f_gate__f; // dimensionless
	const NekDouble var_slow_inward_current__V = var_chaste_interface__membrane__V; // millivolt
	const NekDouble var_slow_inward_current__Cai = var_chaste_interface__intracellular_calcium_concentration__Cai; // millimolar
	const NekDouble var_slow_inward_current__E_si = 7.7 - (13.0287 * log(var_slow_inward_current__Cai / 1.0)); // millivolt
	const NekDouble var_slow_inward_current__i_si = 0.09 * var_slow_inward_current__d * var_slow_inward_current__f * (var_slow_inward_current__V - var_slow_inward_current__E_si); // microA_per_cm2
	const NekDouble var_time_dependent_potassium_current__V = var_chaste_interface__membrane__V; // millivolt
	const NekDouble var_time_dependent_potassium_current__X = var_chaste_interface__time_dependent_potassium_current_X_gate__X; // dimensionless
#if 0 
	const NekDouble var_fast_sodium_current_m_gate__m = var_fast_sodium_current__m; // dimensionless
#endif
	const NekDouble var_fast_sodium_current_m_gate__V = var_fast_sodium_current__V; // millivolt
	const NekDouble var_fast_sodium_current_m_gate__alpha_m = (0.32 * (var_fast_sodium_current_m_gate__V + 47.13)) / (1.0 - exp((-0.1) * (var_fast_sodium_current_m_gate__V + 47.13))); // per_millisecond
	const NekDouble var_fast_sodium_current_m_gate__beta_m = 0.08 * exp((-var_fast_sodium_current_m_gate__V) / 11.0); // per_millisecond
#if 0 
	const NekDouble var_fast_sodium_current_m_gate__d_m_d_environment__time = (var_fast_sodium_current_m_gate__alpha_m * (1.0 - var_fast_sodium_current_m_gate__m)) - (var_fast_sodium_current_m_gate__beta_m * var_fast_sodium_current_m_gate__m); // per_millisecond
	const NekDouble var_fast_sodium_current__fast_sodium_current_m_gate__d_m_d_environment__time = var_fast_sodium_current_m_gate__d_m_d_environment__time; // per_millisecond
#endif
	const NekDouble var_fast_sodium_current_h_gate__V = var_fast_sodium_current__V; // millivolt
	const NekDouble var_fast_sodium_current_h_gate__beta_h = (var_fast_sodium_current_h_gate__V < (-40.0)) ? ((3.56 * exp(0.079 * var_fast_sodium_current_h_gate__V)) + (310000.0 * exp(0.35 * var_fast_sodium_current_h_gate__V))) : (1.0 / (0.13 * (1.0 + exp((var_fast_sodium_current_h_gate__V + 10.66) / (-11.1))))); // per_millisecond
	const NekDouble var_fast_sodium_current_h_gate__alpha_h = (var_fast_sodium_current_h_gate__V < (-40.0)) ? (0.135 * exp((80.0 + var_fast_sodium_current_h_gate__V) / (-6.8))) : 0.0; // per_millisecond
#if 0 
	const NekDouble var_fast_sodium_current_h_gate__h = var_fast_sodium_current__h; // dimensionless
	const NekDouble var_fast_sodium_current_h_gate__d_h_d_environment__time = (var_fast_sodium_current_h_gate__alpha_h * (1.0 - var_fast_sodium_current_h_gate__h)) - (var_fast_sodium_current_h_gate__beta_h * var_fast_sodium_current_h_gate__h); // per_millisecond
	const NekDouble var_fast_sodium_current__fast_sodium_current_h_gate__d_h_d_environment__time = var_fast_sodium_current_h_gate__d_h_d_environment__time; // per_millisecond
#endif
	const NekDouble var_fast_sodium_current_j_gate__V = var_fast_sodium_current__V; // millivolt
	const NekDouble var_fast_sodium_current_j_gate__alpha_j = (var_fast_sodium_current_j_gate__V < (-40.0)) ? (((((-127140.0) * exp(0.2444 * var_fast_sodium_current_j_gate__V)) - (3.474e-05 * exp((-0.04391) * var_fast_sodium_current_j_gate__V))) * (var_fast_sodium_current_j_gate__V + 37.78)) / (1.0 + exp(0.311 * (var_fast_sodium_current_j_gate__V + 79.23)))) : 0.0; // per_millisecond
	const NekDouble var_fast_sodium_current_j_gate__beta_j = (var_fast_sodium_current_j_gate__V < (-40.0)) ? ((0.1212 * exp((-0.01052) * var_fast_sodium_current_j_gate__V)) / (1.0 + exp((-0.1378) * (var_fast_sodium_current_j_gate__V + 40.14)))) : ((0.3 * exp((-2.535e-07) * var_fast_sodium_current_j_gate__V)) / (1.0 + exp((-0.1) * (var_fast_sodium_current_j_gate__V + 32.0)))); // per_millisecond
#if 0 
	const NekDouble var_fast_sodium_current_j_gate__j = var_fast_sodium_current__j; // dimensionless
	const NekDouble var_fast_sodium_current_j_gate__d_j_d_environment__time = (var_fast_sodium_current_j_gate__alpha_j * (1.0 - var_fast_sodium_current_j_gate__j)) - (var_fast_sodium_current_j_gate__beta_j * var_fast_sodium_current_j_gate__j); // per_millisecond
	const NekDouble var_fast_sodium_current__fast_sodium_current_j_gate__d_j_d_environment__time = var_fast_sodium_current_j_gate__d_j_d_environment__time; // per_millisecond
#endif
	const NekDouble var_slow_inward_current_d_gate__V = var_slow_inward_current__V; // millivolt
	const NekDouble var_slow_inward_current_d_gate__alpha_d = (0.095 * exp((-0.01) * (var_slow_inward_current_d_gate__V - 5.0))) / (1.0 + exp((-0.072) * (var_slow_inward_current_d_gate__V - 5.0))); // per_millisecond
#if 0 
	const NekDouble var_slow_inward_current_d_gate__d = var_slow_inward_current__d; // dimensionless
#endif
	const NekDouble var_slow_inward_current_d_gate__beta_d = (0.07 * exp((-0.017) * (var_slow_inward_current_d_gate__V + 44.0))) / (1.0 + exp(0.05 * (var_slow_inward_current_d_gate__V + 44.0))); // per_millisecond
#if 0 
	const NekDouble var_slow_inward_current_d_gate__d_d_d_environment__time = (var_slow_inward_current_d_gate__alpha_d * (1.0 - var_slow_inward_current_d_gate__d)) - (var_slow_inward_current_d_gate__beta_d * var_slow_inward_current_d_gate__d); // per_millisecond
	const NekDouble var_slow_inward_current__slow_inward_current_d_gate__d_d_d_environment__time = var_slow_inward_current_d_gate__d_d_d_environment__time; // per_millisecond
	const NekDouble var_slow_inward_current_f_gate__f = var_slow_inward_current__f; // dimensionless
#endif
	const NekDouble var_slow_inward_current_f_gate__V = var_slow_inward_current__V; // millivolt
	const NekDouble var_slow_inward_current_f_gate__alpha_f = (0.012 * exp((-0.008) * (var_slow_inward_current_f_gate__V + 28.0))) / (1.0 + exp(0.15 * (var_slow_inward_current_f_gate__V + 28.0))); // per_millisecond
	const NekDouble var_slow_inward_current_f_gate__beta_f = (0.0065 * exp((-0.02) * (var_slow_inward_current_f_gate__V + 30.0))) / (1.0 + exp((-0.2) * (var_slow_inward_current_f_gate__V + 30.0))); // per_millisecond
#if 0 
	const NekDouble var_slow_inward_current_f_gate__d_f_d_environment__time = (var_slow_inward_current_f_gate__alpha_f * (1.0 - var_slow_inward_current_f_gate__f)) - (var_slow_inward_current_f_gate__beta_f * var_slow_inward_current_f_gate__f); // per_millisecond
	const NekDouble var_slow_inward_current__slow_inward_current_f_gate__d_f_d_environment__time = var_slow_inward_current_f_gate__d_f_d_environment__time; // per_millisecond
	const NekDouble var_time_dependent_potassium_current_X_gate__X = var_time_dependent_potassium_current__X; // dimensionless
#endif
	const NekDouble var_time_dependent_potassium_current_X_gate__V = var_time_dependent_potassium_current__V; // millivolt
	const NekDouble var_time_dependent_potassium_current_X_gate__beta_X = (0.0013 * exp((-0.06) * (var_time_dependent_potassium_current_X_gate__V + 20.0))) / (1.0 + exp((-0.04) * (var_time_dependent_potassium_current_X_gate__V + 20.0))); // per_millisecond
	const NekDouble var_time_dependent_potassium_current_X_gate__alpha_X = (0.0005 * exp(0.083 * (var_time_dependent_potassium_current_X_gate__V + 50.0))) / (1.0 + exp(0.057 * (var_time_dependent_potassium_current_X_gate__V + 50.0))); // per_millisecond
#if 0 
	const NekDouble var_time_dependent_potassium_current_X_gate__d_X_d_environment__time = (var_time_dependent_potassium_current_X_gate__alpha_X * (1.0 - var_time_dependent_potassium_current_X_gate__X)) - (var_time_dependent_potassium_current_X_gate__beta_X * var_time_dependent_potassium_current_X_gate__X); // per_millisecond
	const NekDouble var_time_dependent_potassium_current__time_dependent_potassium_current_X_gate__d_X_d_environment__time = var_time_dependent_potassium_current_X_gate__d_X_d_environment__time; // per_millisecond
#endif
	const NekDouble var_intracellular_calcium_concentration__Cai = var_chaste_interface__intracellular_calcium_concentration__Cai; // millimolar
	const NekDouble var_intracellular_calcium_concentration__i_si = var_slow_inward_current__i_si; // microA_per_cm2
	const NekDouble var_intracellular_calcium_concentration__d_Cai_d_environment__time = (((-0.0001) / 1.0) * var_intracellular_calcium_concentration__i_si) + (0.07 * (0.0001 - var_intracellular_calcium_concentration__Cai)); // 'millimole per litre per millisecond'
#if 0
	const NekDouble var_chaste_interface__fast_sodium_current_m_gate__d_m_d_environment__time = var_fast_sodium_current__fast_sodium_current_m_gate__d_m_d_environment__time; // per_millisecond
	const NekDouble var_chaste_interface__fast_sodium_current_h_gate__d_h_d_environment__time = var_fast_sodium_current__fast_sodium_current_h_gate__d_h_d_environment__time; // per_millisecond
	const NekDouble var_chaste_interface__fast_sodium_current_j_gate__d_j_d_environment__time = var_fast_sodium_current__fast_sodium_current_j_gate__d_j_d_environment__time; // per_millisecond
	const NekDouble var_chaste_interface__slow_inward_current_d_gate__d_d_d_environment__time = var_slow_inward_current__slow_inward_current_d_gate__d_d_d_environment__time; // per_millisecond
	const NekDouble var_chaste_interface__slow_inward_current_f_gate__d_f_d_environment__time = var_slow_inward_current__slow_inward_current_f_gate__d_f_d_environment__time; // per_millisecond
	const NekDouble var_chaste_interface__time_dependent_potassium_current_X_gate__d_X_d_environment__time = var_time_dependent_potassium_current__time_dependent_potassium_current_X_gate__d_X_d_environment__time; // per_millisecond
#endif
	const NekDouble var_chaste_interface__intracellular_calcium_concentration__d_Cai_d_environment__time = var_intracellular_calcium_concentration__d_Cai_d_environment__time; // ___units_16
#if 0
	const NekDouble d_dt_chaste_interface__fast_sodium_current_m_gate__m = var_chaste_interface__fast_sodium_current_m_gate__d_m_d_environment__time; // per_millisecond
	const NekDouble d_dt_chaste_interface__fast_sodium_current_h_gate__h = var_chaste_interface__fast_sodium_current_h_gate__d_h_d_environment__time; // per_millisecond
	const NekDouble d_dt_chaste_interface__fast_sodium_current_j_gate__j = var_chaste_interface__fast_sodium_current_j_gate__d_j_d_environment__time; // per_millisecond
	const NekDouble d_dt_chaste_interface__slow_inward_current_d_gate__d = var_chaste_interface__slow_inward_current_d_gate__d_d_d_environment__time; // per_millisecond
	const NekDouble d_dt_chaste_interface__slow_inward_current_f_gate__f = var_chaste_interface__slow_inward_current_f_gate__d_f_d_environment__time; // per_millisecond
	const NekDouble d_dt_chaste_interface__time_dependent_potassium_current_X_gate__X = var_chaste_interface__time_dependent_potassium_current_X_gate__d_X_d_environment__time; // per_millisecond
#endif
	const NekDouble d_dt_chaste_interface__intracellular_calcium_concentration__Cai = var_chaste_interface__intracellular_calcium_concentration__d_Cai_d_environment__time; // 'millimole per litre per millisecond'
            
	const NekDouble var_membrane__R = 8314.0; // joule_per_kilomole_kelvin
	const NekDouble var_membrane__T = 310.0; // kelvin
	const NekDouble var_membrane__F = 96484.6; // coulomb_per_mole
	const NekDouble var_membrane__C = 1.0; // microF_per_cm2
	const NekDouble var_chaste_interface__membrane__I_stim = 0.0;
	const NekDouble var_membrane__I_stim = var_chaste_interface__membrane__I_stim; // microA_per_cm2
	const NekDouble var_fast_sodium_current__g_Na = 23.0; // milliS_per_cm2
	const NekDouble var_fast_sodium_current__R = var_membrane__R; // joule_per_kilomole_kelvin
	const NekDouble var_fast_sodium_current__F = var_membrane__F; // coulomb_per_mole
	const NekDouble var_ionic_concentrations__Nao = 140.0; // millimolar
	const NekDouble var_fast_sodium_current__Nao = var_ionic_concentrations__Nao; // millimolar
	const NekDouble var_ionic_concentrations__Nai = 18.0; // millimolar
	const NekDouble var_fast_sodium_current__Nai = var_ionic_concentrations__Nai; // millimolar
	const NekDouble var_fast_sodium_current__T = var_membrane__T; // kelvin
	const NekDouble var_fast_sodium_current__E_Na = ((var_fast_sodium_current__R * var_fast_sodium_current__T) / var_fast_sodium_current__F) * log(var_fast_sodium_current__Nao / var_fast_sodium_current__Nai); // millivolt
	const NekDouble var_fast_sodium_current__i_Na = var_fast_sodium_current__g_Na * pow(var_fast_sodium_current__m, 3.0) * var_fast_sodium_current__h * var_fast_sodium_current__j * (var_fast_sodium_current__V - var_fast_sodium_current__E_Na); // microA_per_cm2
	const NekDouble var_membrane__i_Na = var_fast_sodium_current__i_Na; // microA_per_cm2
	const NekDouble var_membrane__i_si = var_slow_inward_current__i_si; // microA_per_cm2
	const NekDouble var_time_dependent_potassium_current_Xi_gate__V = var_time_dependent_potassium_current__V; // millivolt
	const NekDouble var_time_dependent_potassium_current_Xi_gate__Xi = (var_time_dependent_potassium_current_Xi_gate__V > (-100.0)) ? ((2.837 * (exp(0.04 * (var_time_dependent_potassium_current_Xi_gate__V + 77.0)) - 1.0)) / ((var_time_dependent_potassium_current_Xi_gate__V + 77.0) * exp(0.04 * (var_time_dependent_potassium_current_Xi_gate__V + 35.0)))) : 1.0; // dimensionless
	const NekDouble var_time_dependent_potassium_current__Xi = var_time_dependent_potassium_current_Xi_gate__Xi; // dimensionless
	const NekDouble var_ionic_concentrations__Ko = 5.4; // millimolar
	const NekDouble var_time_dependent_potassium_current__Ko = var_ionic_concentrations__Ko; // millimolar
	const NekDouble var_time_dependent_potassium_current__g_K = 0.282 * sqrt(var_time_dependent_potassium_current__Ko / 5.4); // milliS_per_cm2
	const NekDouble var_time_dependent_potassium_current__PR_NaK = 0.01833; // dimensionless
	const NekDouble var_time_dependent_potassium_current__F = var_membrane__F; // coulomb_per_mole
	const NekDouble var_time_dependent_potassium_current__Nao = var_ionic_concentrations__Nao; // millimolar
	const NekDouble var_ionic_concentrations__Ki = 145.0; // millimolar
	const NekDouble var_time_dependent_potassium_current__Ki = var_ionic_concentrations__Ki; // millimolar
	const NekDouble var_time_dependent_potassium_current__Nai = var_ionic_concentrations__Nai; // millimolar
	const NekDouble var_time_dependent_potassium_current__T = var_membrane__T; // kelvin
	const NekDouble var_time_dependent_potassium_current__R = var_membrane__R; // joule_per_kilomole_kelvin
	const NekDouble var_time_dependent_potassium_current__E_K = ((var_time_dependent_potassium_current__R * var_time_dependent_potassium_current__T) / var_time_dependent_potassium_current__F) * log((var_time_dependent_potassium_current__Ko + (var_time_dependent_potassium_current__PR_NaK * var_time_dependent_potassium_current__Nao)) / (var_time_dependent_potassium_current__Ki + (var_time_dependent_potassium_current__PR_NaK * var_time_dependent_potassium_current__Nai))); // millivolt
	const NekDouble var_time_dependent_potassium_current__i_K = var_time_dependent_potassium_current__g_K * var_time_dependent_potassium_current__X * var_time_dependent_potassium_current__Xi * (var_time_dependent_potassium_current__V - var_time_dependent_potassium_current__E_K); // microA_per_cm2
	const NekDouble var_membrane__i_K = var_time_dependent_potassium_current__i_K; // microA_per_cm2
	const NekDouble var_time_independent_potassium_current__V = var_chaste_interface__membrane__V; // millivolt
	const NekDouble var_time_independent_potassium_current_K1_gate__V = var_time_independent_potassium_current__V; // millivolt
	const NekDouble var_time_independent_potassium_current__Ki = var_ionic_concentrations__Ki; // millimolar
	const NekDouble var_time_independent_potassium_current__R = var_membrane__R; // joule_per_kilomole_kelvin
	const NekDouble var_time_independent_potassium_current__F = var_membrane__F; // coulomb_per_mole
	const NekDouble var_time_independent_potassium_current__Ko = var_ionic_concentrations__Ko; // millimolar
	const NekDouble var_time_independent_potassium_current__T = var_membrane__T; // kelvin
	const NekDouble var_time_independent_potassium_current__E_K1 = ((var_time_independent_potassium_current__R * var_time_independent_potassium_current__T) / var_time_independent_potassium_current__F) * log(var_time_independent_potassium_current__Ko / var_time_independent_potassium_current__Ki); // millivolt
	const NekDouble var_time_independent_potassium_current_K1_gate__E_K1 = var_time_independent_potassium_current__E_K1; // millivolt
	const NekDouble var_time_independent_potassium_current_K1_gate__beta_K1 = ((0.49124 * exp(0.08032 * ((var_time_independent_potassium_current_K1_gate__V + 5.476) - var_time_independent_potassium_current_K1_gate__E_K1))) + (1.0 * exp(0.06175 * (var_time_independent_potassium_current_K1_gate__V - (var_time_independent_potassium_current_K1_gate__E_K1 + 594.31))))) / (1.0 + exp((-0.5143) * ((var_time_independent_potassium_current_K1_gate__V - var_time_independent_potassium_current_K1_gate__E_K1) + 4.753))); // per_millisecond
	const NekDouble var_time_independent_potassium_current_K1_gate__alpha_K1 = 1.02 / (1.0 + exp(0.2385 * ((var_time_independent_potassium_current_K1_gate__V - var_time_independent_potassium_current_K1_gate__E_K1) - 59.215))); // per_millisecond
	const NekDouble var_time_independent_potassium_current_K1_gate__K1_infinity = var_time_independent_potassium_current_K1_gate__alpha_K1 / (var_time_independent_potassium_current_K1_gate__alpha_K1 + var_time_independent_potassium_current_K1_gate__beta_K1); // dimensionless
	const NekDouble var_time_independent_potassium_current__K1_infinity = var_time_independent_potassium_current_K1_gate__K1_infinity; // dimensionless
	const NekDouble var_time_independent_potassium_current__g_K1 = 0.6047 * sqrt(var_time_independent_potassium_current__Ko / 5.4); // milliS_per_cm2
	const NekDouble var_time_independent_potassium_current__i_K1 = var_time_independent_potassium_current__g_K1 * var_time_independent_potassium_current__K1_infinity * (var_time_independent_potassium_current__V - var_time_independent_potassium_current__E_K1); // microA_per_cm2
	const NekDouble var_membrane__i_K1 = var_time_independent_potassium_current__i_K1; // microA_per_cm2
	const NekDouble var_plateau_potassium_current__g_Kp = 0.0183; // milliS_per_cm2
	const NekDouble var_plateau_potassium_current__V = var_chaste_interface__membrane__V; // millivolt
	const NekDouble var_plateau_potassium_current__Kp = 1.0 / (1.0 + exp((7.488 - var_plateau_potassium_current__V) / 5.98)); // dimensionless
	const NekDouble var_plateau_potassium_current__E_K1 = var_time_independent_potassium_current__E_K1; // millivolt
	const NekDouble var_plateau_potassium_current__E_Kp = var_plateau_potassium_current__E_K1; // millivolt
	const NekDouble var_plateau_potassium_current__i_Kp = var_plateau_potassium_current__g_Kp * var_plateau_potassium_current__Kp * (var_plateau_potassium_current__V - var_plateau_potassium_current__E_Kp); // microA_per_cm2
	const NekDouble var_membrane__i_Kp = var_plateau_potassium_current__i_Kp; // microA_per_cm2
	const NekDouble var_background_current__E_b =  -59.87; // millivolt
	const NekDouble var_background_current__g_b = 0.03921; // milliS_per_cm2
	const NekDouble var_background_current__V = var_chaste_interface__membrane__V; // millivolt
	const NekDouble var_background_current__i_b = var_background_current__g_b * (var_background_current__V - var_background_current__E_b); // microA_per_cm2
	const NekDouble var_membrane__i_b = var_background_current__i_b; // microA_per_cm2
	const NekDouble var_membrane__d_V_d_environment__time = ((-1.0) / var_membrane__C) * (var_membrane__I_stim + var_membrane__i_Na + var_membrane__i_si + var_membrane__i_K + var_membrane__i_K1 + var_membrane__i_Kp + var_membrane__i_b); // 'millivolt per millisecond'
	const NekDouble var_chaste_interface__membrane__d_V_d_environment__time = var_membrane__d_V_d_environment__time; // ___units_1
	d_dt_chaste_interface__membrane__V = var_chaste_interface__membrane__d_V_d_environment__time; // 'millivolt per millisecond'
	const NekDouble m_inf = var_fast_sodium_current_m_gate__alpha_m/(var_fast_sodium_current_m_gate__alpha_m + var_fast_sodium_current_m_gate__beta_m);
	const NekDouble m_tau = 1.0/(var_fast_sodium_current_m_gate__alpha_m + var_fast_sodium_current_m_gate__beta_m);
	const NekDouble h_inf = var_fast_sodium_current_h_gate__alpha_h/(var_fast_sodium_current_h_gate__alpha_h + var_fast_sodium_current_h_gate__beta_h);
	const NekDouble h_tau = 1.0/(var_fast_sodium_current_h_gate__alpha_h + var_fast_sodium_current_h_gate__beta_h);
	const NekDouble j_inf = var_fast_sodium_current_j_gate__alpha_j/(var_fast_sodium_current_j_gate__alpha_j + var_fast_sodium_current_j_gate__beta_j);
	const NekDouble j_tau = 1.0/(var_fast_sodium_current_j_gate__alpha_j + var_fast_sodium_current_j_gate__beta_j);
	const NekDouble d_inf = var_slow_inward_current_d_gate__alpha_d/(var_slow_inward_current_d_gate__alpha_d + var_slow_inward_current_d_gate__beta_d);
	const NekDouble d_tau = 1.0/(var_slow_inward_current_d_gate__alpha_d + var_slow_inward_current_d_gate__beta_d);
	const NekDouble f_inf = var_slow_inward_current_f_gate__alpha_f/(var_slow_inward_current_f_gate__alpha_f + var_slow_inward_current_f_gate__beta_f);
	const NekDouble f_tau = 1.0/(var_slow_inward_current_f_gate__alpha_f + var_slow_inward_current_f_gate__beta_f);
	const NekDouble X_inf = var_time_dependent_potassium_current_X_gate__alpha_X/(var_time_dependent_potassium_current_X_gate__alpha_X + var_time_dependent_potassium_current_X_gate__beta_X);
	const NekDouble X_tau = 1.0/(var_time_dependent_potassium_current_X_gate__alpha_X + var_time_dependent_potassium_current_X_gate__beta_X);

	outarray[0][i] = d_dt_chaste_interface__membrane__V;
	outarray[1][i] = m_inf;
	gates_tau[0][i] = m_tau;
	outarray[2][i] = h_inf;
	gates_tau[1][i] = h_tau;
	outarray[3][i] = j_inf;
	gates_tau[2][i] = j_tau;
	outarray[4][i] = d_inf;
	gates_tau[3][i] = d_tau;
	outarray[5][i] = f_inf;
	gates_tau[4][i] = f_tau;
	outarray[6][i] = X_inf;
	gates_tau[5][i] = X_tau;
	outarray[7][i] = d_dt_chaste_interface__intracellular_calcium_concentration__Cai;
    }
}
