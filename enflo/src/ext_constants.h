#ifndef __EXT_CONSTANTS_H__
#define __EXT_CONSTANTS_H__

#define   NVAR   5

// Coefficients for RK-schemes
static const double arks[] = {0.0, 3.0/4.0, 1.0/3.0};
static const double brks[] = {1.0, 1.0/4.0, 2.0/3.0};
static const double jameson_rks[] = {1.0/4.0,1.0/3.0,1.0/2.0, 1.0};
static const double rk4_rks[] = {1.0,1.0/2.0,1.0/2.0,1.0};

//sp_weno --> original one proposed for scalar con laws
//sp_wenoc --> modified SP-WENO for systems of conservation laws (fourth order perturbation)
//sp_weno_corr2 --> experimental modification to SP-WENOc
enum ReconstructionScheme {first, second,tvd_minmod, eno, weno,sp_weno, sp_wenoc, sp_weno_dl, sp_weno_corr2, minabs_lim};

enum FluxScheme {roe, kep, kep4, kepec, kepec4, roe_ec, roe_ec4, keps, keps4,
                 kep_tecno_roe, kepes_tecno_roe, kepes_tecno_rusanov, roe_tecno_roe, 
                 kepes_tecno4_roe, kepes_tecno4_rusanov, roe_tecno4_roe,};
enum TimeIntegrationScheme {rk1, ssprk3, jameson_rk4, rk4};
enum FlowModel {euler, ns_kep, ns_es};
enum MuModel {mu_constant, mu_sutherland, mu_power};
enum BCType {periodic, dirichlet, neumann, wall};

#endif
