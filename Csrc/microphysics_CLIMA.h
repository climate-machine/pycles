#pragma once
#include <math.h>
#include "parameters.h"
#include "thermodynamic_functions.h"
#include "advection_interpolation.h"
#include "entropies.h"

// CLIMA microphysics parameters
// see CLIMA Microphysics docs for info
#define rho_cloud_liq 1e3
#define rho_cloud_ice 916.7

#define r_0_rai 1e-3
#define m_0_rai 4./3 * pi * rho_cloud_liq * pow(r_0_rai, 3)
#define m_e_rai 3
#define a_0_rai pi * pow(r_0_rai, 2)
#define a_e_rai 2
#define v_e_rai 0.5

#define r_0_ice 1e-5
#define m_0_ice 4./3 * pi * rho_cloud_ice * pow(r_0_ice, 3)
#define m_e_ice 3

#define r_0_sno 1e-3
#define m_0_sno 0.1 * pow(r_0_sno, 2)
#define m_e_sno 2
#define a_0_sno 0.3 * pi * pow(r_0_sno, 2)
#define a_e_sno 2
#define v_0_sno pow(2, 9./4) * pow(r_0_sno, 1./4)
#define v_e_sno 0.25

// the above parameters are hardcoded when used inside gamma functions
#define gamma_3 2.
#define gamma_4 6.
#define gamma_5 24.
#define gamma_6 120.
#define gamma_7_2 3.3233509704478426
#define gamma_9_2 11.631728396567448
#define gamma_11_4 1.6083594219855457
#define gamma_13_4 2.549256966718529281826
#define gamma_13_2 287.8852778150443609963
#define gamma_21_8 1.456933205091971725255

#define n_0_rai 16 * 1e6
#define n_0_ice 2 * 1e7
#define mu_sno 4.36 * 1e9
#define nu_sno 0.63

#define C_drag 0.55

#define tau_cond_evap 10
#define tau_dep_sub 10

#define q_liq_threshold 5e-4
#define tau_acnv 1e3

#define r_ice_sno 62.5 * 1e-6

#define T_freeze  273.15

#define E_col_liq_rai 0.8
#define E_col_liq_sno 0.1
#define E_col_ice_rai 1.0
#define E_col_ice_sno 0.1
#define E_col_rai_sno 1.0

#define nu_air 1.6e-5
#define K_therm 2.4e-2
#define D_vapor 2.26e-5

//TODO - make sure they are consistent with generated parameters in PyCLES
#define L_vapor 2.501e6
#define L_subli 2.8334e6
#define L_fusion L_subli - L_vapor

#define a_rai_vent 1.5
#define b_rai_vent 0.53
#define a_sno_vent 0.65
#define b_sno_vent 0.44

// additional parameters for pycles implementation
#define max_iter 10
#define microph_eps 1e-3

double CLIMA_v_0_rai(double rho){

    double v_0_rai = sqrt(8./3 / C_drag * (rho_cloud_liq / rho - 1.) * g * r_0_rai);
    return v_0_rai;
}

double CLIMA_n_0_sno(double q_sno, double rho){

    double n_0_sno = mu_sno * pow(rho * fmax(0, q_sno), nu_sno);
    return n_0_sno;
}

double CLIMA_lambda_ice(double q_ice, double rho){

    double lambda_ice = 0.;
    if(q_ice > 0){
        lambda_ice = pow(
            m_0_ice * n_0_ice * gamma_4 / rho / q_ice / pow(r_0_ice, m_e_ice),
            1.0 / (m_e_ice + 1)
        );
    }
    return lambda_ice;
}

double CLIMA_lambda_rai(double q_rai, double rho){

    double lambda_rai = 0.;
    if(q_rai > 0){
        lambda_rai = pow(
            m_0_rai * n_0_rai * gamma_4 / rho / q_rai / pow(r_0_rai, m_e_rai),
            1.0 / (m_e_rai + 1)
        );
    }
    return lambda_rai;
}

double CLIMA_lambda_sno(double q_sno, double rho){

    double n_0_sno = CLIMA_n_0_sno(q_sno, rho);

    double lambda_sno = 0.;
    if(q_sno > 0){
        lambda_sno = pow(
            m_0_sno * n_0_sno * gamma_3 / rho / q_sno / pow(r_0_sno, m_e_sno),
            1.0 / (m_e_sno + 1)
        );
    }
    return lambda_sno;
}

double CLIMA_terminal_velocity_rai(double q_rai, double rho){

    double v_0_rai = CLIMA_v_0_rai(rho);
    double lambda_rai = CLIMA_lambda_rai(q_rai, rho);

    double term_vel = 0.;
    if(q_rai > 0){
        term_vel = v_0_rai * pow(lambda_rai * r_0_rai, -v_e_rai) * gamma_9_2 / gamma_4;
    }
    return term_vel;
}

double CLIMA_terminal_velocity_sno(double q_sno, double rho){

    double lambda_sno = CLIMA_lambda_sno(q_sno, rho);

    double term_vel = 0.;
    if(q_sno > 0){
        term_vel = v_0_sno * pow(lambda_sno * r_0_sno, -v_e_sno) * gamma_13_4 / gamma_3;
    }
    return term_vel;
}

void CLIMA_conv_q_liq_to_q_rai(double q_liq, double* qr_tendency_aut){

    *qr_tendency_aut = fmax(0., q_liq - q_liq_threshold) / tau_acnv;
    return;
}

void CLIMA_conv_q_ice_to_q_sno(double q_tot, double q_liq, double q_ice,
                               double T, double p, double rho,
                               struct LookupStruct *LT,
                               double (*lam_fp)(double),
                               double (*L_fp)(double, double),
                               double* qs_tendency_aut){
    double lam = lam_fp(T);
    double L = L_fp(T, lam);
    double pv_s = lookup(LT, T);

    double qv_sat = qv_star_c(p, q_tot, pv_s);
    double q_v = q_tot - q_liq - q_ice;
    double S = q_v/qv_sat - 1;

    double acnv_rate = 0.;
    if (q_ice > 0. && S > 0.){

        double G_param = 1. / (L_subli / K_therm / T * (L_subli / Rv / T - 1.) +
                           Rv * T / D_vapor / pv_s
                          );
        double lambda_ice = CLIMA_lambda_ice(q_ice, rho);

        acnv_rate =
            4. * pi * S * G_param * n_0_ice / rho *
            exp(-lambda_ice * r_ice_sno) *
            (
                pow(r_ice_sno, 2) / (m_e_ice) +
                (r_ice_sno * lambda_ice + 1) / pow(lambda_ice,2)
            );
    }

    *qs_tendency_aut = fmax(0., acnv_rate);
    return;
}

void CLIMA_q_liq_q_rai_accr(double q_liq, double q_rai, double rho,
                            double* qr_tendency_acc_q_liq){

    double v_0_rai = CLIMA_v_0_rai(rho);
    double lambda_rai = CLIMA_lambda_rai(q_rai, rho);

    double accr_rate = 0;
    if (q_liq > 0. && q_rai > 0.){
        accr_rate  =
            q_liq * E_col_liq_rai * n_0_rai * a_0_rai * v_0_rai / lambda_rai *
            gamma_7_2 /
            pow(lambda_rai * r_0_rai, a_e_rai + v_e_rai);
    }

    *qr_tendency_acc_q_liq = fmax(0., accr_rate);
    return;
}

void CLIMA_q_ice_q_sno_accr(double q_ice, double q_sno, double rho,
                            double* qs_tendency_acc_q_ice){

    double lambda_sno = CLIMA_lambda_sno(q_sno, rho);
    double n_0_sno = CLIMA_n_0_sno(q_sno, rho);

    double accr_rate = 0;
    if (q_ice > 0. && q_sno > 0.){
        accr_rate =
            q_ice * E_col_ice_sno * n_0_sno * a_0_sno * v_0_sno / lambda_sno *
            gamma_13_4 /
            pow(lambda_sno * r_0_sno, a_e_sno + v_e_sno);
    }

    *qs_tendency_acc_q_ice = accr_rate;
    return;
}

void CLIMA_q_liq_q_sno_accr(double q_liq, double q_sno, double rho, double T,
                            double* ql_tendency_acc_q_liq_q_sno,
                            double* qs_tendency_acc_q_liq_q_sno,
                            double* qr_tendency_acc_q_liq_q_sno){

    double lambda_sno = CLIMA_lambda_sno(q_sno, rho);
    double n_0_sno = CLIMA_n_0_sno(q_sno, rho);

    double accr_rate = 0;
    if (q_liq > 0. && q_sno > 0.){
        accr_rate =
            q_liq * E_col_liq_sno * n_0_sno * a_0_sno * v_0_sno / lambda_sno *
            gamma_13_4 /
            pow(lambda_sno * r_0_sno, a_e_sno + v_e_sno);
    }
    if (T < T_freeze){
        *qs_tendency_acc_q_liq_q_sno = accr_rate;
        *ql_tendency_acc_q_liq_q_sno = -accr_rate;
        *qr_tendency_acc_q_liq_q_sno = 0.;
    }
    else{
        *qs_tendency_acc_q_liq_q_sno = -accr_rate * cl / L_fusion * (T  - T_freeze);
        *ql_tendency_acc_q_liq_q_sno = -accr_rate;
        *qr_tendency_acc_q_liq_q_sno = accr_rate * (1.0 + cl / L_fusion * (T  - T_freeze));
    }
    return;
}

void CLIMA_q_ice_q_rai_accr(double q_ice, double q_rai, double rho,
                            double* qi_sink_acc_q_ice_q_rai){

    double v_0_rai = CLIMA_v_0_rai(rho);
    double lambda_rai = CLIMA_lambda_rai(q_rai, rho);

    double accr_rate = 0;
    if (q_ice > 0. && q_rai > 0.){
        accr_rate =
            q_ice * E_col_ice_rai * n_0_rai * a_0_rai * v_0_rai / lambda_rai *
            gamma_7_2 /
            pow(lambda_rai * r_0_rai, a_e_rai + v_e_rai);
    }

    *qi_sink_acc_q_ice_q_rai = fmax(0., accr_rate);
    return;
}

void CLIMA_accretion_rain_sink(double q_ice, double q_rai, double rho,
                               double* qr_sink_acc_q_ice_q_rai){

    double lambda_rai = CLIMA_lambda_rai(q_rai, rho);
    double lambda_ice = CLIMA_lambda_ice(q_ice, rho);

    double v_0_rai = CLIMA_v_0_rai(rho);

    double accr_rate = 0;
    if (q_ice > 0. && q_rai > 0.){
        accr_rate =
            E_col_ice_rai / rho * n_0_rai * n_0_ice * m_0_rai * a_0_rai * v_0_rai /
            lambda_ice / lambda_rai *
            gamma_13_2 /
            pow(r_0_rai * lambda_rai, m_e_rai + a_e_rai + v_e_rai);
    }

    *qr_sink_acc_q_ice_q_rai = accr_rate;
    return;
}

void CLIMA_accretion_snow_rain(double q_rai, double q_sno, double rho, double T,
                               double* qr_qs_accretion_rain, double* qr_qs_accretion_snow){

    double lambda_rai = CLIMA_lambda_rai(q_rai, rho);
    double lambda_sno = CLIMA_lambda_sno(q_sno, rho);

    double n_0_sno = CLIMA_n_0_sno(q_sno, rho);

    double v_rai = CLIMA_terminal_velocity_rai(q_rai, rho);
    double v_sno = CLIMA_terminal_velocity_sno(q_sno, rho);

    double accr_rate_snow = 0.;
    double accr_rate_rain = 0.;
    double accr_rate = 0.;
    if (q_rai > 0. && q_sno > 0){
        if (T > T_freeze){
             accr_rate = pi / rho * n_0_rai * n_0_sno * m_0_sno * E_col_rai_sno *
                fabs(v_rai - v_sno) / pow(r_0_sno, m_e_sno) * (
                    2. * gamma_3 / pow(lambda_rai, 3) /
                    pow(lambda_sno, m_e_sno + 1) +
                    2. * gamma_4 / pow(lambda_rai, 2) /
                    pow(lambda_sno, m_e_sno + 2) +
                    gamma_5 / lambda_rai /
                    pow(lambda_sno, m_e_sno + 3)
                );
        accr_rate_rain = accr_rate;
        accr_rate_snow = -accr_rate;
        }
        else{
            accr_rate = pi / rho * n_0_sno * n_0_rai * m_0_rai * E_col_rai_sno *
                fabs(v_sno - v_rai) / pow(r_0_rai, m_e_rai) * (
                    2. * gamma_4 / pow(lambda_sno, 3) /
                    pow(lambda_rai, m_e_rai + 1) +
                    2. * gamma_5 / pow(lambda_sno, 2) /
                    pow(lambda_rai, m_e_rai + 2) +
                    gamma_6 / lambda_sno /
                    pow(lambda_rai, m_e_rai + 3)
                );
        accr_rate_snow = accr_rate;
        accr_rate_rain = -accr_rate;
        }
    }
    *qr_qs_accretion_rain = accr_rate_rain;
    *qr_qs_accretion_snow = accr_rate_snow;
    return;
}

void CLIMA_rain_evaporation(
    double q_rai, double q_tot, double q_liq, double q_ice, double T, double p,
    double rho, struct LookupStruct *LT, double (*lam_fp)(double),
    double (*L_fp)(double, double), double* qr_tendency_evp){

    double v_0_rai = CLIMA_v_0_rai(rho);
    double lambda_rai = CLIMA_lambda_rai(q_rai, rho);

    double N_Sc = nu_air / D_vapor;

    double lam = lam_fp(T);
    double L = L_fp(T, lam);
    double pv_s = lookup(LT, T);

    double qv_sat = qv_star_c(p, q_tot, pv_s);
    double q_v = q_tot - q_liq - q_ice;
    double S = q_v/qv_sat - 1;

    double G_param = 1./ (L_vapor / K_therm / T * (L_vapor / Rv / T - 1.) +
                     Rv * T / D_vapor / pv_s);

    double evap_rate = 0.;
    if (q_rai > 0. && S < 0.){
        evap_rate =
            4. * pi * n_0_rai / rho * S * G_param / pow(lambda_rai, 2) * (
                a_rai_vent +
                b_rai_vent * pow(nu_air / D_vapor, 1./3) /
                pow(r_0_rai * lambda_rai, v_e_rai / 2) *
                pow(2 * v_0_rai / nu_air / lambda_rai, 1./2) *
                gamma_11_4
            );
    }
    *qr_tendency_evp = fmin(0., evap_rate);
    return;
}

void CLIMA_snow_deposition_sublimation(
    double q_sno, double q_tot, double q_liq, double q_ice, double T, double p,
    double rho, struct LookupStruct *LT, double (*lam_fp)(double),
    double (*L_fp)(double, double), double* qs_tendency_dep_sub){

    double n_0_sno = CLIMA_n_0_sno(q_sno, rho);
    double lambda_sno = CLIMA_lambda_sno(q_sno, rho);

    double N_Sc = nu_air / D_vapor;

    double lam = lam_fp(T);
    double L = L_fp(T, lam);
    double pv_s = lookup(LT, T);

    double qv_sat = qv_star_c(p, q_tot, pv_s);
    double q_v = q_tot - q_liq - q_ice;
    double S = q_v/qv_sat - 1;

    double G_param = 1./ (L_subli / K_therm / T * (L_subli / Rv / T - 1.) +
                     Rv * T / D_vapor / pv_s);

    double dep_sub_rate = 0.;
    if (q_sno > 0.){
        dep_sub_rate =
            4. * pi * n_0_sno / rho * S * G_param / pow(lambda_sno, 2) * (
                a_sno_vent +
                b_sno_vent * pow(nu_air / D_vapor, 1./3) /
                pow(r_0_sno * lambda_sno, v_e_sno / 2.) *
                pow(2 * v_0_sno / nu_air / lambda_sno, 1./2) *
                gamma_21_8
            );
    }
    *qs_tendency_dep_sub = dep_sub_rate;
    return;
}

void CLIMA_snow_melt(double q_sno, double rho, double T, double* qs_tendency_melt){

    double lambda_sno = CLIMA_lambda_sno(q_sno, rho);
    double n_0_sno = CLIMA_n_0_sno(q_sno, rho);

    double snow_melt_rate = 0.;
    if (q_sno > 0. && T > T_freeze){

        snow_melt_rate =
            4. * pi * n_0_sno / rho * K_therm / L_fusion * (T - T_freeze) / pow(lambda_sno, 2) * (
                a_sno_vent +
                b_sno_vent * pow(nu_air / D_vapor, 1./3) /
                pow(r_0_sno * lambda_sno, v_e_sno / 2.) *
                pow(2. * v_0_sno / nu_air / lambda_sno, 1./2) *
                gamma_21_8
            );
    }
    *qs_tendency_melt = snow_melt_rate;
    return;
}


void CLIMA_sedimentation_velocities(const struct DimStruct *dims,
                                    double* restrict density,
                                    double* restrict qr,
                                    double* restrict qs,
                                    double* restrict qr_velocity,
                                    double* restrict qs_velocity){

    const ssize_t istride = dims->nlg[1] * dims->nlg[2];
    const ssize_t jstride = dims->nlg[2];
    const ssize_t imin = 0;
    const ssize_t jmin = 0;
    const ssize_t kmin = 0;
    const ssize_t imax = dims->nlg[0];
    const ssize_t jmax = dims->nlg[1];
    const ssize_t kmax = dims->nlg[2];

    for(ssize_t i=imin; i<imax; i++){
        const ssize_t ishift = i * istride;
        for(ssize_t j=jmin; j<jmax; j++){
            const ssize_t jshift = j * jstride;
            for(ssize_t k=kmin-1; k<kmax+1; k++){
                const ssize_t ijk = ishift + jshift + k;

                double qr_tmp = fmax(qr[ijk], 0.0);
                double qs_tmp = fmax(qs[ijk], 0.0);

                qr_velocity[ijk] = -fmin(
                    CLIMA_terminal_velocity_rai(qr_tmp, density[k]),
                    10.0
                );
                qs_velocity[ijk] = -fmin(
                    CLIMA_terminal_velocity_sno(qs_tmp, density[k]),
                    10.0
                );
            }
        }
    }

    for(ssize_t i=imin; i<imax; i++){
        const ssize_t ishift = i * istride;
        for(ssize_t j=jmin; j<jmax; j++){
            const ssize_t jshift = j * jstride;
            for(ssize_t k=kmin; k<kmax-1 ; k++){
                const ssize_t ijk = ishift + jshift + k;

                qr_velocity[ijk] = interp_2(qr_velocity[ijk], qr_velocity[ijk+1]);
                qs_velocity[ijk] = interp_2(qs_velocity[ijk], qs_velocity[ijk+1]);
            }
        }
    }
    return;
}

void CLIMA_precipitation_sources(const struct DimStruct *dims, struct LookupStruct *LT,
                                 double (*lam_fp)(double), double (*L_fp)(double, double),
                                 double* restrict density, double* restrict p0,
                                 double* restrict temperature, double* restrict qt,
                                 double* restrict ql, double* restrict qi,
                                 double* restrict qr, double* restrict qs, double dt,
                                 double* restrict qr_tendency_micro, double* restrict qr_tendency,
                                 double* restrict qs_tendency_micro, double* restrict qs_tendency,
                                 double* restrict coll_rate, double* restrict evap_rate, double* restrict melt_rate){

    // helpers for storing tendencies for rain, snow, cloud liquid and ice and water vapor
    double qr_tendency_tmp, qs_tendency_tmp, ql_tendency_tmp, qi_tendency_tmp, qv_tendency_tmp;

    // helpers for tendencies from individual microphysics processes
    double qr_tendency_aut, qs_tendency_aut;
    double qr_tendency_acc_q_liq, qs_tendency_acc_q_ice;
    double ql_tendency_acc_q_liq_q_sno, qs_tendency_acc_q_liq_q_sno, qr_tendency_acc_q_liq_q_sno;
    double qi_sink_acc_q_ice_q_rai, qr_sink_acc_q_ice_q_rai;
    double qr_qs_accretion_rain, qr_qs_accretion_snow;
    double qr_tendency_evp, qs_tendency_dep_sub, qs_tendency_melt;

    // helpers for computing entropy sources
    double coll_rate_tmp, evap_rate_tmp, melt_rate_tmp;

    const ssize_t istride = dims->nlg[1] * dims->nlg[2];
    const ssize_t jstride = dims->nlg[2];
    const ssize_t imin = dims->gw;
    const ssize_t jmin = dims->gw;
    const ssize_t kmin = dims->gw;
    const ssize_t imax = dims->nlg[0]-dims->gw;
    const ssize_t jmax = dims->nlg[1]-dims->gw;
    const ssize_t kmax = dims->nlg[2]-dims->gw;

    for(ssize_t i=imin; i<imax; i++){
        const ssize_t ishift = i * istride;
        for(ssize_t j=jmin; j<jmax; j++){
            const ssize_t jshift = j * jstride;
            for(ssize_t k=kmin; k<kmax; k++){
                const ssize_t ijk = ishift + jshift + k;

                double qr_tmp = fmax(qr[ijk], 0.0);
                double qs_tmp = fmax(qs[ijk], 0.0);

                double qt_tmp = fmax(qt[ijk], 0.0);
                double ql_tmp = fmax(ql[ijk], 0.0);
                double qi_tmp = fmax(qi[ijk], 0.0);
                double qv_tmp = fmax(qt_tmp - ql_tmp - qi_tmp, 0.0);

                coll_rate_tmp = 0.0;
                evap_rate_tmp = 0.0;
                melt_rate_tmp = 0.0;

                // helpers for timestep limiting
                double time_added = 0.0, dt_;
                double rate_factor = 1.05;
                double rate_qv, rate_ql, rate_qi, rate_qr, rate_qs, rate_max;
                ssize_t iter_count = 0;
                do{
                    iter_count += 1;
                    qv_tendency_tmp = 0.0;
                    ql_tendency_tmp = 0.0;
                    qi_tendency_tmp = 0.0;
                    qr_tendency_tmp = 0.0;
                    qs_tendency_tmp = 0.0;

                    // compute the source terms from ...

                    //... autoconversion ...
                    // ql -> qr
                    qr_tendency_aut = 0.0;
                    CLIMA_conv_q_liq_to_q_rai(ql_tmp, &qr_tendency_aut);
                    qr_tendency_tmp += qr_tendency_aut;
                    ql_tendency_tmp -= qr_tendency_aut;
                    // qi -> qs
                    qs_tendency_aut = 0.0;
                    CLIMA_conv_q_ice_to_q_sno(qt_tmp, ql_tmp, qi_tmp,
                                              temperature[ijk], p0[k], density[k],
                                              LT, lam_fp, L_fp, &qs_tendency_aut);
                    qs_tendency_tmp += qs_tendency_aut;
                    qi_tendency_tmp -= qs_tendency_aut;

                    //... accretion ...
                    // qr + ql -> qr
                    qr_tendency_acc_q_liq = 0.0;
                    CLIMA_q_liq_q_rai_accr(
                        ql_tmp, qr_tmp, density[k], &qr_tendency_acc_q_liq
                    );
                    qr_tendency_tmp += qr_tendency_acc_q_liq;
                    ql_tendency_tmp -= qr_tendency_acc_q_liq;
                    // qs + qi -> qs
                    qs_tendency_acc_q_ice = 0.0;
                    CLIMA_q_ice_q_sno_accr(
                        qi_tmp, qs_tmp, density[k], &qs_tendency_acc_q_ice
                    );
                    qs_tendency_tmp += qs_tendency_acc_q_ice;
                    qi_tendency_tmp -= qs_tendency_acc_q_ice;
                    // qs + ql -> (T < T_freeze) ? qs : qr
                    qs_tendency_acc_q_liq_q_sno = 0.0;
                    ql_tendency_acc_q_liq_q_sno = 0.0;
                    qr_tendency_acc_q_liq_q_sno = 0.0;
                    CLIMA_q_liq_q_sno_accr(
                        ql_tmp, qs_tmp, density[k], temperature[ijk],
                        &ql_tendency_acc_q_liq_q_sno,
                        &qs_tendency_acc_q_liq_q_sno,
                        &qr_tendency_acc_q_liq_q_sno
                    );
                    ql_tendency_tmp += ql_tendency_acc_q_liq_q_sno;
                    qs_tendency_tmp += qs_tendency_acc_q_liq_q_sno;
                    qr_tendency_tmp += qr_tendency_acc_q_liq_q_sno;
                    // qr + qi -> qs
                    qi_sink_acc_q_ice_q_rai = 0.0;
                    qr_sink_acc_q_ice_q_rai = 0.0;
                    CLIMA_q_ice_q_rai_accr(
                        qi_tmp, qr_tmp, density[k], &qi_sink_acc_q_ice_q_rai
                    );

                    CLIMA_accretion_rain_sink(
                        qi_tmp, qr_tmp, density[k], &qr_sink_acc_q_ice_q_rai
                    );
                    qi_tendency_tmp -= qi_sink_acc_q_ice_q_rai;
                    qr_tendency_tmp -= qr_sink_acc_q_ice_q_rai;
                    qs_tendency_tmp += qi_sink_acc_q_ice_q_rai + qr_sink_acc_q_ice_q_rai;
                    // qr + qs -> (T < T_freeze) ? qs : qr
                    qr_qs_accretion_rain = 0.0;
                    qr_qs_accretion_snow = 0.0;
                    CLIMA_accretion_snow_rain(
                        qr_tmp, qs_tmp, density[k], temperature[ijk],
                        &qr_qs_accretion_rain, &qr_qs_accretion_snow
                    );
                    qs_tendency_tmp += qr_qs_accretion_rain;
                    qr_tendency_tmp += qr_qs_accretion_snow;

                    //... rain evaporation, snow deposition, sublimation, melting
                    // qr -> qv
                    qr_tendency_evp = 0.0;
                    CLIMA_rain_evaporation(
                        qr_tmp, qt_tmp, ql_tmp, qi_tmp, temperature[ijk], p0[k],
                        density[k], LT, lam_fp, L_fp, &qr_tendency_evp
                    );
                    qr_tendency_tmp += qr_tendency_evp;
                    qv_tendency_tmp -= qr_tendency_evp;

                    // qs -> qv (or reverse)
                    qs_tendency_dep_sub = 0.0;
                    CLIMA_snow_deposition_sublimation(
                        qs_tmp, qt_tmp, ql_tmp, qi_tmp, temperature[ijk], p0[k],
                        density[k], LT, lam_fp, L_fp, &qs_tendency_dep_sub
                    );
                    qs_tendency_tmp += qs_tendency_dep_sub;
                    qv_tendency_tmp -= qs_tendency_dep_sub;

                    // qs -> qr
                    qs_tendency_melt = 0.0;
                    CLIMA_snow_melt(qs_tmp, density[k], temperature[ijk], &qs_tendency_melt);
                    qs_tendency_tmp -= qs_tendency_melt;
                    qr_tendency_tmp += qs_tendency_melt;

                    // find the maximum substep time...
                    dt_ = dt - time_added;
                    //... check the source term magnitudes and
                    //    adjust the rates if necessary (factor of 1.05 is ad-hoc)
                    rate_qv = rate_factor * qv_tendency_tmp * dt_ / (-fmax(qv_tmp, microph_eps));
                    rate_ql = rate_factor * ql_tendency_tmp * dt_ / (-fmax(ql_tmp, microph_eps));
                    rate_qi = rate_factor * qi_tendency_tmp * dt_ / (-fmax(qi_tmp, microph_eps));
                    rate_qr = rate_factor * qr_tendency_tmp * dt_ / (-fmax(qr_tmp, microph_eps));
                    rate_qs = rate_factor * qs_tendency_tmp * dt_ / (-fmax(qs_tmp, microph_eps));
                    rate_max = fmax(rate_qv, fmax(rate_ql, fmax(rate_qi, fmax(rate_qr, rate_qs))));
                    if(rate_max > 1.0 && iter_count < max_iter){
                        //Limit the timestep, but don't allow it to become vanishingly small
                        //Don't adjust if we have reached the maximum iteration number
                        dt_ = fmax(dt_/rate_max, 1.0e-3);
                    }

                    // augument entropy source helpers
                    coll_rate_tmp += (
                        - qr_tendency_aut - qs_tendency_aut
                        - qr_tendency_acc_q_liq - qs_tendency_acc_q_ice
                        + ql_tendency_acc_q_liq_q_sno
                        - qi_sink_acc_q_ice_q_rai
                    ) * dt_;
                    evap_rate_tmp += (qr_tendency_evp + qs_tendency_dep_sub) * dt_;
                    melt_rate_tmp -= qs_tendency_melt * dt_;

                    // integrate forward in time
                    qv_tmp += qv_tendency_tmp * dt_;
                    ql_tmp += ql_tendency_tmp * dt_;
                    qi_tmp += qi_tendency_tmp * dt_;
                    qr_tmp += qr_tendency_tmp * dt_;
                    qs_tmp += qs_tendency_tmp * dt_;

                    // fend off negative values
                    qv_tmp = fmax(qv_tmp, 0.0);
                    ql_tmp = fmax(ql_tmp, 0.0);
                    qi_tmp = fmax(qi_tmp, 0.0);
                    qr_tmp = fmax(qr_tmp, 0.0);
                    qs_tmp = fmax(qs_tmp, 0.0);
                    qt_tmp = qv_tmp + ql_tmp + qi_tmp;

                    time_added += dt_ ;

                }while(time_added < dt);
                // augment precipitation tendencies for rain ...
                qr_tendency_micro[ijk] = (qr_tmp - qr[ijk])/dt;
                qr_tendency[ijk] += qr_tendency_micro[ijk];
                // ... and snow
                qs_tendency_micro[ijk] = (qs_tmp - qs[ijk])/dt;
                qs_tendency[ijk] += qs_tendency_micro[ijk];
                // update entropy source helpers
                coll_rate[ijk] = coll_rate_tmp / dt;
                evap_rate[ijk] = evap_rate_tmp / dt;
                melt_rate[ijk] = melt_rate_tmp / dt;
            }
        }
    }
    return;
}

void CLIMA_qt_source_formation(const struct DimStruct *dims,
                               double* restrict qr_tendency,
                               double* restrict qs_tendency,
                               double* restrict qt_tendency){

    const ssize_t istride = dims->nlg[1] * dims->nlg[2];
    const ssize_t jstride = dims->nlg[2];
    const ssize_t imin = dims->gw;
    const ssize_t jmin = dims->gw;
    const ssize_t kmin = dims->gw;
    const ssize_t imax = dims->nlg[0]-dims->gw;
    const ssize_t jmax = dims->nlg[1]-dims->gw;
    const ssize_t kmax = dims->nlg[2]-dims->gw;

    for(ssize_t i=imin; i<imax; i++){
        const ssize_t ishift = i * istride;
        for(ssize_t j=jmin; j<jmax; j++){
            const ssize_t jshift = j * jstride;
            for(ssize_t k=kmin; k<kmax; k++){
                const ssize_t ijk = ishift + jshift + k;

                qt_tendency[ijk] -= qr_tendency[ijk] + qs_tendency[ijk];
            }
        }
    }
    return;
}

// Source terms of entropy related to microphysics.
// See Pressel et al. 2015, Eq. 49-54
void CLIMA_entropy_source_formation(const struct DimStruct *dims, struct LookupStruct *LT,
                                    double (*lam_fp)(double), double (*L_fp)(double, double),
                                    double* restrict p0,
                                    double* restrict T, double* restrict Twet,
                                    double* restrict qt, double* restrict qv,
                                    double* restrict coll_rate, double* restrict evap_rate,
                                    double* restrict entropy_tendency){

    const ssize_t istride = dims->nlg[1] * dims->nlg[2];
    const ssize_t jstride = dims->nlg[2];
    const ssize_t imin = dims->gw;
    const ssize_t jmin = dims->gw;
    const ssize_t kmin = dims->gw;
    const ssize_t imax = dims->nlg[0]-dims->gw;
    const ssize_t jmax = dims->nlg[1]-dims->gw;
    const ssize_t kmax = dims->nlg[2]-dims->gw;

    // entropy tendencies from formation, evaporation,
    // deposition/sublimation of precipitation
    for(ssize_t i=imin; i<imax; i++){
        const ssize_t ishift = i * istride;
        for(ssize_t j=jmin; j<jmax; j++){
            const ssize_t jshift = j * jstride;
            for(ssize_t k=kmin; k<kmax; k++){
                const ssize_t ijk = ishift + jshift + k;

                const double lam_T = lam_fp(T[ijk]);
                const double L_fp_T = L_fp(T[ijk],lam_T);

                const double lam_Tw = lam_fp(Twet[ijk]);
                const double L_fp_Tw = L_fp(Twet[ijk],lam_Tw);

                const double pv_star_T = lookup(LT, T[ijk]);
                const double pv_star_Tw = lookup(LT,Twet[ijk]);

                const double pv = pv_c(p0[k], qt[ijk], qv[ijk]);
                const double pd = p0[k] - pv;

                const double sd_T = sd_c(pd, T[ijk]);
                const double sv_star_T = sv_c(pv_star_T,T[ijk] );
                const double sv_star_Tw = sv_c(pv_star_Tw, Twet[ijk]);

                const double S_P = sd_T - sv_star_T + L_fp_T/T[ijk];
                const double S_E = sv_star_Tw - L_fp_Tw/Twet[ijk] - sd_T;
                const double S_D = -Rv * log(pv/pv_star_T) + cpv * log(T[ijk]/Twet[ijk]);

                entropy_tendency[ijk] +=  S_P        * coll_rate[ijk];
                entropy_tendency[ijk] -= (S_E + S_D) * evap_rate[ijk];
            }
        }
    }
    return;
}

void CLIMA_entropy_source_heating(const struct DimStruct *dims,
                               double* restrict T, double* restrict Twet,
                               double* restrict qr, double* restrict qs,
                               double* restrict w_qr, double* restrict w_qs, double* restrict w,
                               double* restrict melt_rate,
                               double* restrict entropy_tendency){

    const ssize_t istride = dims->nlg[1] * dims->nlg[2];
    const ssize_t jstride = dims->nlg[2];
    const ssize_t imin = dims->gw;
    const ssize_t jmin = dims->gw;
    const ssize_t kmin = dims->gw;
    const ssize_t imax = dims->nlg[0]-dims->gw;
    const ssize_t jmax = dims->nlg[1]-dims->gw;
    const ssize_t kmax = dims->nlg[2]-dims->gw;
    const double dzi = 1.0/dims->dx[2];

    for(ssize_t i=imin; i<imax; i++){
        const ssize_t ishift = i * istride;
        for(ssize_t j=jmin; j<jmax; j++){
            const ssize_t jshift = j * jstride;
            for(ssize_t k=kmin; k<kmax; k++){
                const ssize_t ijk = ishift + jshift + k;

                const double w_adv_rai = w_qr[ijk] - w[ijk];
                const double w_adv_sno = w_qs[ijk] - w[ijk];

                if (qr[ijk] > 0.0){
                    // heat transfer to falling rain
                    if (w_adv_rai < 0.0){
                        entropy_tendency[ijk] += qr[ijk] * fabs(w_adv_rai) * cl * (Twet[ijk+1] - Twet[ijk]) * dzi / T[ijk];
                    }
                    else
                    {
                        entropy_tendency[ijk] += qr[ijk] * w_adv_rai * cl * (Twet[ijk] - Twet[ijk-1]) * dzi / T[ijk];
                    }
                }
                if (qs[ijk] > 0.0){
                    // heat transfer to falling snow
                    if (w_adv_sno < 0.0){
                        entropy_tendency[ijk] += qs[ijk] * fabs(w_adv_sno) * ci * (Twet[ijk+1] - Twet[ijk]) * dzi / T[ijk];
                    }
                    else
                    {
                        entropy_tendency[ijk] += qs[ijk] * w_adv_sno * ci * (Twet[ijk] - Twet[ijk-1]) * dzi / T[ijk];
                    }
                    // heat transfer to melting snow
                    entropy_tendency[ijk] += melt_rate[ijk] * L_fusion / T[ijk];
                }
            }
        }
    }
    return;
}

void CLIMA_entropy_source_drag(const struct DimStruct *dims,
                               double* restrict T,
                               double* restrict q_precip, double* restrict w_precip,
                               double* restrict entropy_tendency){

    const ssize_t istride = dims->nlg[1] * dims->nlg[2];
    const ssize_t jstride = dims->nlg[2];
    const ssize_t imin = dims->gw;
    const ssize_t jmin = dims->gw;
    const ssize_t kmin = dims->gw;
    const ssize_t imax = dims->nlg[0]-dims->gw;
    const ssize_t jmax = dims->nlg[1]-dims->gw;
    const ssize_t kmax = dims->nlg[2]-dims->gw;

    for(ssize_t i=imin; i<imax; i++){
        const ssize_t ishift = i * istride;
        for(ssize_t j=jmin; j<jmax; j++){
            const ssize_t jshift = j * jstride;
            for(ssize_t k=kmin; k<kmax; k++){
                const ssize_t ijk = ishift + jshift + k;

                entropy_tendency[ijk] += g * q_precip[ijk] * fabs(w_precip[ijk]) / T[ijk];
            }                                                     //TODO - shouldnt that be velocity relative to air velocity?
        }
    }
    return;
}

// diagnostics for output
void CLIMA_autoconversion_rain_wrapper(const struct DimStruct *dims,
                                       double* restrict ql,
                                       double* restrict qr_tendency){

    const ssize_t istride = dims->nlg[1] * dims->nlg[2];
    const ssize_t jstride = dims->nlg[2];
    const ssize_t imin = dims->gw;
    const ssize_t jmin = dims->gw;
    const ssize_t kmin = dims->gw;
    const ssize_t imax = dims->nlg[0]-dims->gw;
    const ssize_t jmax = dims->nlg[1]-dims->gw;
    const ssize_t kmax = dims->nlg[2]-dims->gw;

    for(ssize_t i=imin; i<imax; i++){
        const ssize_t ishift = i * istride;
        for(ssize_t j=jmin; j<jmax; j++){
            const ssize_t jshift = j * jstride;
            for(ssize_t k=kmin; k<kmax; k++){
                const ssize_t ijk = ishift + jshift + k;

                //compute the source terms
                double ql_tmp = fmax(ql[ijk], 0.0);

                CLIMA_conv_q_liq_to_q_rai(ql_tmp, &qr_tendency[ijk]);
            }
        }
    }
    return;
}
void CLIMA_conv_q_ice_to_q_sno_wrapper(const struct DimStruct *dims,
                                       double* restrict qt,
                                       double* restrict ql,
                                       double* restrict qi,
                                       double* restrict temperature,
                                       double* restrict p0,
                                       double* restrict density,
                                       struct LookupStruct *LT,
                                       double (*lam_fp)(double),
                                       double (*L_fp)(double, double),
                                       double* restrict qs_tendency){

    const ssize_t istride = dims->nlg[1] * dims->nlg[2];
    const ssize_t jstride = dims->nlg[2];
    const ssize_t imin = dims->gw;
    const ssize_t jmin = dims->gw;
    const ssize_t kmin = dims->gw;
    const ssize_t imax = dims->nlg[0]-dims->gw;
    const ssize_t jmax = dims->nlg[1]-dims->gw;
    const ssize_t kmax = dims->nlg[2]-dims->gw;

    for(ssize_t i=imin; i<imax; i++){
        const ssize_t ishift = i * istride;
        for(ssize_t j=jmin; j<jmax; j++){
            const ssize_t jshift = j * jstride;
            for(ssize_t k=kmin; k<kmax; k++){
                const ssize_t ijk = ishift + jshift + k;

                const double qt_tmp = fmax(qt[ijk], 0.0);
                const double ql_tmp = fmax(ql[ijk], 0.0);
                const double qi_tmp = fmax(qi[ijk], 0.0);

                CLIMA_conv_q_ice_to_q_sno(qt_tmp, ql_tmp, qi_tmp,
                                          temperature[ijk], p0[k], density[k],
                                          LT, lam_fp, L_fp, &qs_tendency);
            }
        }
    }
    return;
}
void CLIMA_q_liq_q_rai_accr_wrapper(const struct DimStruct *dims,
                                    double* restrict ql,
                                    double* restrict qr,
                                    double* restrict density,
                                    double* restrict qr_tendency){

    const ssize_t istride = dims->nlg[1] * dims->nlg[2];
    const ssize_t jstride = dims->nlg[2];
    const ssize_t imin = dims->gw;
    const ssize_t jmin = dims->gw;
    const ssize_t kmin = dims->gw;
    const ssize_t imax = dims->nlg[0]-dims->gw;
    const ssize_t jmax = dims->nlg[1]-dims->gw;
    const ssize_t kmax = dims->nlg[2]-dims->gw;

    for(ssize_t i=imin; i<imax; i++){
        const ssize_t ishift = i * istride;
        for(ssize_t j=jmin; j<jmax; j++){
            const ssize_t jshift = j * jstride;
            for(ssize_t k=kmin; k<kmax; k++){
                const ssize_t ijk = ishift + jshift + k;

                const double ql_tmp = fmax(ql[ijk], 0.0);
                const double qr_tmp = fmax(qr[ijk], 0.0);

                CLIMA_q_liq_q_rai_accr(
                    ql_tmp, qr_tmp, density[k], &qr_tendency[ijk]
                );
            }
        }
    }
    return;
}
void CLIMA_q_ice_q_sno_accr_wrapper(const struct DimStruct *dims,
                                    double* restrict qi,
                                    double* restrict qs,
                                    double* restrict density,
                                    double* restrict qs_tendency){

    const ssize_t istride = dims->nlg[1] * dims->nlg[2];
    const ssize_t jstride = dims->nlg[2];
    const ssize_t imin = dims->gw;
    const ssize_t jmin = dims->gw;
    const ssize_t kmin = dims->gw;
    const ssize_t imax = dims->nlg[0]-dims->gw;
    const ssize_t jmax = dims->nlg[1]-dims->gw;
    const ssize_t kmax = dims->nlg[2]-dims->gw;

    for(ssize_t i=imin; i<imax; i++){
        const ssize_t ishift = i * istride;
        for(ssize_t j=jmin; j<jmax; j++){
            const ssize_t jshift = j * jstride;
            for(ssize_t k=kmin; k<kmax; k++){
                const ssize_t ijk = ishift + jshift + k;

                const double qi_tmp = fmax(qi[ijk], 0.0);
                const double qs_tmp = fmax(qs[ijk], 0.0);

                CLIMA_q_ice_q_sno_accr(
                    qi_tmp, qs_tmp, density[k], &qs_tendency[ijk]
                );
            }
        }
    }
    return;
}
void CLIMA_q_liq_q_sno_accr_wrapper(const struct DimStruct *dims,
                                    double* restrict ql,
                                    double* restrict qs,
                                    double* restrict density,
                                    double* restrict temperature,
                                    double* restrict ql_tendency,
                                    double* restrict qs_tendency,
                                    double* restrict qr_tendency){

    const ssize_t istride = dims->nlg[1] * dims->nlg[2];
    const ssize_t jstride = dims->nlg[2];
    const ssize_t imin = dims->gw;
    const ssize_t jmin = dims->gw;
    const ssize_t kmin = dims->gw;
    const ssize_t imax = dims->nlg[0]-dims->gw;
    const ssize_t jmax = dims->nlg[1]-dims->gw;
    const ssize_t kmax = dims->nlg[2]-dims->gw;

    for(ssize_t i=imin; i<imax; i++){
        const ssize_t ishift = i * istride;
        for(ssize_t j=jmin; j<jmax; j++){
            const ssize_t jshift = j * jstride;
            for(ssize_t k=kmin; k<kmax; k++){
                const ssize_t ijk = ishift + jshift + k;

                const double ql_tmp = fmax(ql[ijk], 0.0);
                const double qs_tmp = fmax(qs[ijk], 0.0);

                CLIMA_q_liq_q_sno_accr(
                    ql_tmp, qs_tmp, density[k], temperature[ijk],
                    &ql_tendency[ijk], &qs_tendency[ijk], &qr_tendency[ijk]
                );
            }
        }
    }
    return;
}
void CLIMA_q_ice_q_rai_accr_wrapper(const struct DimStruct *dims,
                                    double* restrict qi,
                                    double* restrict qr,
                                    double* restrict density,
                                    double* restrict qi_tendency){

    const ssize_t istride = dims->nlg[1] * dims->nlg[2];
    const ssize_t jstride = dims->nlg[2];
    const ssize_t imin = dims->gw;
    const ssize_t jmin = dims->gw;
    const ssize_t kmin = dims->gw;
    const ssize_t imax = dims->nlg[0]-dims->gw;
    const ssize_t jmax = dims->nlg[1]-dims->gw;
    const ssize_t kmax = dims->nlg[2]-dims->gw;

    for(ssize_t i=imin; i<imax; i++){
        const ssize_t ishift = i * istride;
        for(ssize_t j=jmin; j<jmax; j++){
            const ssize_t jshift = j * jstride;
            for(ssize_t k=kmin; k<kmax; k++){
                const ssize_t ijk = ishift + jshift + k;

                const double qi_tmp = fmax(qi[ijk], 0.0);
                const double qr_tmp = fmax(qr[ijk], 0.0);
                CLIMA_q_ice_q_rai_accr(
                    qi_tmp, qr_tmp, density[k], &qi_tendency[ijk]
                );
            }
        }
    }
    return;
}
void CLIMA_accretion_rain_sink_wrapper(const struct DimStruct *dims,
                                       double* restrict qi,
                                       double* restrict qr,
                                       double* restrict density,
                                       double* restrict qr_tendency){

    const ssize_t istride = dims->nlg[1] * dims->nlg[2];
    const ssize_t jstride = dims->nlg[2];
    const ssize_t imin = dims->gw;
    const ssize_t jmin = dims->gw;
    const ssize_t kmin = dims->gw;
    const ssize_t imax = dims->nlg[0]-dims->gw;
    const ssize_t jmax = dims->nlg[1]-dims->gw;
    const ssize_t kmax = dims->nlg[2]-dims->gw;

    for(ssize_t i=imin; i<imax; i++){
        const ssize_t ishift = i * istride;
        for(ssize_t j=jmin; j<jmax; j++){
            const ssize_t jshift = j * jstride;
            for(ssize_t k=kmin; k<kmax; k++){
                const ssize_t ijk = ishift + jshift + k;

                const double qi_tmp = fmax(qi[ijk], 0.0);
                const double qr_tmp = fmax(qr[ijk], 0.0);
                CLIMA_accretion_rain_sink(
                    qi_tmp, qr_tmp, density[k], &qr_tendency[ijk]
                );
            }
        }
    }
    return;
}
void CLIMA_accretion_snow_rain_wrapper(const struct DimStruct *dims,
                                       double* restrict qr,
                                       double* restrict qs,
                                       double* restrict density,
                                       double* restrict temperature,
                                       double* restrict qs_tendency,
                                       double* restrict qr_tendency){

    const ssize_t istride = dims->nlg[1] * dims->nlg[2];
    const ssize_t jstride = dims->nlg[2];
    const ssize_t imin = dims->gw;
    const ssize_t jmin = dims->gw;
    const ssize_t kmin = dims->gw;
    const ssize_t imax = dims->nlg[0]-dims->gw;
    const ssize_t jmax = dims->nlg[1]-dims->gw;
    const ssize_t kmax = dims->nlg[2]-dims->gw;

    for(ssize_t i=imin; i<imax; i++){
        const ssize_t ishift = i * istride;
        for(ssize_t j=jmin; j<jmax; j++){
            const ssize_t jshift = j * jstride;
            for(ssize_t k=kmin; k<kmax; k++){
                const ssize_t ijk = ishift + jshift + k;

                const double qr_tmp = fmax(qr[ijk], 0.0);
                const double qs_tmp = fmax(qs[ijk], 0.0);

                CLIMA_accretion_snow_rain(
                    qr_tmp, qs_tmp, density[k], temperature[ijk],
                    &qs_tendency[ijk], &qr_tendency[ijk]
                );
            }
        }
    }
    return;
}
void CLIMA_rain_evaporation_wrapper(const struct DimStruct *dims,
                                    double* restrict qr,
                                    double* restrict qt,
                                    double* restrict ql,
                                    double* restrict qi,
                                    double* restrict temperature,
                                    double* restrict p0,
                                    double* restrict density,
                                    struct LookupStruct *LT,
                                    double (*lam_fp)(double),
                                    double (*L_fp)(double, double),
                                    double* restrict qr_tendency){

    const ssize_t istride = dims->nlg[1] * dims->nlg[2];
    const ssize_t jstride = dims->nlg[2];
    const ssize_t imin = dims->gw;
    const ssize_t jmin = dims->gw;
    const ssize_t kmin = dims->gw;
    const ssize_t imax = dims->nlg[0]-dims->gw;
    const ssize_t jmax = dims->nlg[1]-dims->gw;
    const ssize_t kmax = dims->nlg[2]-dims->gw;

    for(ssize_t i=imin; i<imax; i++){
        const ssize_t ishift = i * istride;
        for(ssize_t j=jmin; j<jmax; j++){
            const ssize_t jshift = j * jstride;
            for(ssize_t k=kmin; k<kmax; k++){
                const ssize_t ijk = ishift + jshift + k;

                const double qr_tmp = fmax(qr[ijk], 0.0);
                const double qt_tmp = fmax(qt[ijk], 0.0);
                const double ql_tmp = fmax(ql[ijk], 0.0);
                const double qi_tmp = fmax(qi[ijk], 0.0);

                CLIMA_rain_evaporation(
                    qr_tmp, qt_tmp, ql_tmp, qi_tmp, temperature[ijk], p0[k],
                    density[k], LT, lam_fp, L_fp, &qr_tendency[ijk]);
            }
        }
    }
    return;
}
void CLIMA_snow_deposition_sublimation_wrapper(
    const struct DimStruct *dims,
    double* restrict qs,
    double* restrict qt,
    double* restrict ql,
    double* restrict qi,
    double* restrict temperature,
    double* restrict p0,
    double* restrict density,
    struct LookupStruct *LT,
    double (*lam_fp)(double),
    double (*L_fp)(double, double),
    double* restrict qs_tendency){

    const ssize_t istride = dims->nlg[1] * dims->nlg[2];
    const ssize_t jstride = dims->nlg[2];
    const ssize_t imin = dims->gw;
    const ssize_t jmin = dims->gw;
    const ssize_t kmin = dims->gw;
    const ssize_t imax = dims->nlg[0]-dims->gw;
    const ssize_t jmax = dims->nlg[1]-dims->gw;
    const ssize_t kmax = dims->nlg[2]-dims->gw;

    for(ssize_t i=imin; i<imax; i++){
        const ssize_t ishift = i * istride;
        for(ssize_t j=jmin; j<jmax; j++){
            const ssize_t jshift = j * jstride;
            for(ssize_t k=kmin; k<kmax; k++){
                const ssize_t ijk = ishift + jshift + k;

                const double qs_tmp = fmax(qs[ijk], 0.0);
                const double qt_tmp = fmax(qt[ijk], 0.0);
                const double ql_tmp = fmax(ql[ijk], 0.0);
                const double qi_tmp = fmax(qi[ijk], 0.0);

                CLIMA_snow_deposition_sublimation(
                    qs_tmp, qt_tmp, ql_tmp, qi_tmp, temperature[ijk], p0[k],
                    density[k], LT, lam_fp, L_fp, &qs_tendency[ijk]
                );
            }
        }
    }
    return;
}
void CLIMA_snow_melt_wrapper(
    const struct DimStruct *dims,
    double* restrict qs,
    double* restrict density,
    double* restrict temperature,
    double* restrict qs_tendency){

    const ssize_t istride = dims->nlg[1] * dims->nlg[2];
    const ssize_t jstride = dims->nlg[2];
    const ssize_t imin = dims->gw;
    const ssize_t jmin = dims->gw;
    const ssize_t kmin = dims->gw;
    const ssize_t imax = dims->nlg[0]-dims->gw;
    const ssize_t jmax = dims->nlg[1]-dims->gw;
    const ssize_t kmax = dims->nlg[2]-dims->gw;

    for(ssize_t i=imin; i<imax; i++){
        const ssize_t ishift = i * istride;
        for(ssize_t j=jmin; j<jmax; j++){
            const ssize_t jshift = j * jstride;
            for(ssize_t k=kmin; k<kmax; k++){
                const ssize_t ijk = ishift + jshift + k;

                const double qs_tmp = fmax(qs[ijk], 0.0);
                CLIMA_snow_melt(qs_tmp, density[k], temperature[ijk], &qs_tendency[ijk]);
            }
        }
    }
    return;
}
