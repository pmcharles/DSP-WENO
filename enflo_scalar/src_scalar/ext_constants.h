#ifndef __EXT_CONSTANTS_H__
#define __EXT_CONSTANTS_H__

#define   NVAR   1

// Coefficients for RK-schemes
static const double arks[] = {0.0, 3.0/4.0, 1.0/3.0};
static const double brks[] = {1.0, 1.0/4.0, 2.0/3.0};
static const double jameson_rks[] = {1.0/4.0,1.0/3.0,1.0/2.0, 1.0};
static const double rk4_rks[] = {1.0,1.0/2.0,1.0/2.0,1.0};

//sp_weno --> original one proposed for scalar con laws
//sp_wenoc --> modified SP-WENO for systems of conservation laws (fourth order perturbation)
//sp_weno_corr2 --> experimental modification to SP-WENOc
enum ReconstructionScheme {first, second, tvd_minmod, eno, weno, sp_weno, sp_wenoc, sp_weno_dl, sp_weno_corr2, minabs_lim};

enum FluxScheme {lin_ec, lin_ec4, lin_tecno, lin_tecno4, bur_ec, bur_ec4, bur_tecno, bur_tecno4 };
enum TimeIntegrationScheme {rk1, ssprk3, jameson_rk4, rk4};
enum FlowModel {linadv, burgers};
enum BCType {periodic, dirichlet, neumann, wall};

#endif
