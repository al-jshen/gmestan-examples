#include /shared_functions.stan

data {

    int<lower=1> N; // total length of data

    int<lower=0, upper=1> pos_obs[N]; // boolean array indicating whether position (ra, dec) is observed
    int<lower=0, upper=1> dist_obs[N]; // boolean array indicating whether distance is observed
    int<lower=0, upper=1> pm_obs[N]; // boolean array indicating whether pm is observed
    int<lower=0, upper=1> vlos_obs[N]; // boolean array indicating whether vlos is observed

    // position measurements
    // ra, dec in degrees
    // dist in kpc
    row_vector[N] ra_measured;
    row_vector[N] dec_measured;
    row_vector[N] dist_measured;

    // position measurement uncertainties
    row_vector[N] ra_err;
    row_vector[N] dec_err;
    row_vector[N] dist_err;

    // velocity measurements
    // pmra, pmdec in mas/yr
    // vlos in km/s
    row_vector[N] pmra_measured;
    row_vector[N] pmdec_measured;
    row_vector[N] vlos_measured;

    // velocity measurement uncertainties */
    row_vector[N] pmra_err;
    row_vector[N] pmdec_err;
    row_vector[N] vlos_err;

    vector[N] pos_corr; // ra with dec
    vector[N] pm_corr; // pmra with pmdec
    vector[N] ra_pmra_corr;
    vector[N] dec_pmdec_corr;
    vector[N] ra_pmdec_corr;
    vector[N] dec_pmra_corr;

    real alpha_mean;
    real alpha_sigma;
    real beta_mean;
    real beta_sigma;

}

transformed data {

    /* int grainsize = 1; */

    row_vector[2] pg_mean = [64.6317304 ,  0.44004421];
    cholesky_factor_cov[2] pg_sigma = cholesky_decompose([
      [4.76112184e+01, 1.65061266e-01],
      [1.65061266e-01, 2.09403334e-03]
    ]);

    real deg2rad = pi() / 180.;

    row_vector[N] ra_rad = ra_measured * deg2rad;
    row_vector[N] dec_rad = dec_measured * deg2rad;


    // icrs to gctc matrices and vectors
    matrix[3, 3] R = [
        [-0.05487395617553902, -0.8734371822248346, -0.48383503143198114],
        [0.4941107627040048, -0.4448286178025452, 0.7469819642829028],
        [-0.8676654903323697, -0.1980782408317943, 0.4559842183620723]
    ];
    matrix[3, 3] H = [
        [0.9999967207734917, 0.0, 0.002560945579906427],
        [0.0, 1.0, 0.0],
        [-0.002560945579906427, 0.0, 0.9999967207734917]
    ];
    matrix[3, 3] A = H * R;
    vector[3] offset = to_vector([8.122, 0., 0.]);
    vector[3] solarmotion = to_vector([12.9, 245.6, 7.78]);

    // array of N 2x2 covariance matrices
    cholesky_factor_cov[2] pos_cov_mats[N];
    vector[2] pos_measured[N];

    cholesky_factor_cov[2] pm_cov_mats[N];
    vector[2] pm_measured[N];

    cholesky_factor_cov[4] pos_pm_cov_mats[N];
    vector[4] pos_pm_measured[N];

    real dummy[N] = rep_array(1, N);

    for (i in 1:N) {

        // use biggest covariance matrix possible

        if (pos_obs[i] && pm_obs[i]) {
            /* diag(S) * R * diag(S) */
            /* where S is a vector of standard deviations (errs) */
            /* and R is the correlation matrix */

            // make correlation matrix
            // diagonals are 1s
            for (j in 1:4) {
                pos_pm_cov_mats[i, j, j] = 1.;
            }
            // fill in with correlation data
            pos_pm_cov_mats[i, 1, 2] = pos_corr[i];
            pos_pm_cov_mats[i, 2, 1] = pos_corr[i];
            pos_pm_cov_mats[i, 1, 3] = ra_pmra_corr[i];
            pos_pm_cov_mats[i, 3, 1] = ra_pmra_corr[i];
            pos_pm_cov_mats[i, 1, 4] = ra_pmdec_corr[i];
            pos_pm_cov_mats[i, 4, 1] = ra_pmdec_corr[i];
            pos_pm_cov_mats[i, 2, 3] = dec_pmra_corr[i];
            pos_pm_cov_mats[i, 3, 2] = dec_pmra_corr[i];
            pos_pm_cov_mats[i, 2, 4] = dec_pmdec_corr[i];
            pos_pm_cov_mats[i, 4, 2] = dec_pmdec_corr[i];
            pos_pm_cov_mats[i, 3, 4] = pm_corr[i];
            pos_pm_cov_mats[i, 4, 3] = pm_corr[i];

            // then convert it to a covariance matrix
            pos_pm_cov_mats[i] = quad_form_diag(pos_pm_cov_mats[i], [ra_err[i], dec_err[i], pmra_err[i], pmdec_err[i]]);

            // and apply cholesky decomposition
            pos_pm_cov_mats[i] = cholesky_decompose(pos_pm_cov_mats[i]);

        } else {
            // these won't be used
            pos_pm_cov_mats[i] = diag_matrix(to_vector([1, 1, 1, 1]));
        }

        if (pos_obs[i]) {
            // first construct a correlation matrix
            pos_cov_mats[i, 1, 1] = 1.; // ra correlation with itself is 1
            pos_cov_mats[i, 2, 2] = 1.; // dec correlation with itself is 1
            pos_cov_mats[i, 1, 2] = pos_corr[i]; // correlation between ra and dec
            pos_cov_mats[i, 2, 1] = pos_corr[i]; // correlation between ra and dec

            // then convert it to a covariance matrix
            pos_cov_mats[i] = quad_form_diag(pos_cov_mats[i], [ra_err[i], dec_err[i]]);

            // and apply cholesky decomposition
            pos_cov_mats[i] = cholesky_decompose(pos_cov_mats[i]);
        }  else {
            // these won't be used
            pos_cov_mats[i] = diag_matrix(to_vector([1, 1]));
        }

        if (pm_obs[i]) {
            // first construct a correlation matrix
            pm_cov_mats[i, 1, 1] = 1.; // pmra correlation with itself is 1
            pm_cov_mats[i, 2, 2] = 1.; // pmdec correlation with itself is 1
            pm_cov_mats[i, 1, 2] = pm_corr[i]; // correlation between pmra and pmdec
            pm_cov_mats[i, 2, 1] = pm_corr[i]; // correlation between pmra and pmdec

            // then convert it to a covariance matrix
            pm_cov_mats[i] = quad_form_diag(pm_cov_mats[i], [pmra_err[i], pmdec_err[i]]);

            // and apply cholesky decomposition
            pm_cov_mats[i] = cholesky_decompose(pm_cov_mats[i]);
        }  else {
            // these won't be used
            pm_cov_mats[i] = diag_matrix(to_vector([1, 1]));
        }
    }

    for (i in 1:N) {
        pos_measured[i] = [ra_measured[i], dec_measured[i]]';
        pm_measured[i] = [pmra_measured[i], pmdec_measured[i]]';
        pos_pm_measured[i] = [ra_measured[i], dec_measured[i], pmra_measured[i], pmdec_measured[i]]';
    }

}

parameters {

    // parameters for df
    real<lower=0, upper=1> p_gamma; // power-law slope of gravitational potential
    real<lower=0> p_phi0;
    real<upper=1> p_beta;
    real<lower=max([3., p_beta * (2. - p_gamma) + p_gamma / 2.])> p_alpha; // power-law slope of satellite population

    // position parameters for each tracer
    /* row_vector<lower=0, upper=360.>[N] ra; // degrees */
    /* row_vector<lower=-90., upper=90.>[N] dec; // degrees */
    row_vector<lower=0, upper=max(dist_measured + 5 * dist_err)>[N] dist; // kpc

    // velocity parameters for each tracer
    row_vector<lower=min(pmra_measured - 5 * pmra_err), upper=max(pmra_measured + 5 * pmra_err)>[N] pmra; // mas/yr
    row_vector<lower=min(pmdec_measured - 5 * pmdec_err), upper=max(pmdec_measured + 5 * pmdec_err)>[N] pmdec; // mas/yr
    row_vector<lower=min(vlos_measured - 5 * vlos_err), upper=max(vlos_measured + 5 * vlos_err)>[N] vlos; // km/s

}

transformed parameters {

    real beta_std = (p_beta - beta_mean) / beta_sigma;
    real alpha_std = (p_alpha - alpha_mean) / alpha_sigma;

    /* real<lower=0> p_phi0 = phi0_raw + 18.21597964588313; */

    /* vector[2] pos[N]; */
    /* vector[2] pm[N]; // mas/yr */
    matrix[4, N] pos_pm = [ra_measured, dec_measured, pmra, pmdec]; // mas/yr for the pm part

    matrix[3, N] vels_sph;

    row_vector[N] dist_std = (dist_measured - dist) ./ dist_err;
    row_vector[N] vlos_std = (vlos_measured - vlos) ./ vlos_err;

    matrix[3, N] y;

    matrix[3, N] pos_gc;

    profile("transform_vels") {
        vels_sph = transform_vels_vec(ra_rad, dec_rad, dist, pmra, pmdec, vlos, R, H, offset, solarmotion) ./ 100.; // units of 100km/s
    }

    profile("y and position transformation") {
        row_vector[N] vt_sq = square(vels_sph[2]) + square(vels_sph[3]);
        y[2] = sqrt(vt_sq);
        y[1] = sqrt(square(vels_sph[1]) + vt_sq);
        pos_gc = transform_pos_vec(ra_rad, dec_rad, dist, R, H);;
        y[3] = sqrt(columns_dot_self(pos_gc));;
    }

}


model {

    // hyperpriors ----------------

    profile("hyperpriors") {
        [p_phi0, p_gamma] ~ multi_normal_cholesky(pg_mean, pg_sigma);
        /* phi0_raw ~ gamma(48.884347, 1./0.9755005176256963); */
        /* p_beta ~ normal(beta_mean, beta_sigma); */
        /* p_alpha ~ normal(alpha_mean, alpha_sigma); */
        beta_std ~ std_normal();
        alpha_std ~ std_normal();
    }

    // no explicit priors on "true" parameters because the DF is the prior

    profile("observation") {

        // full process with support for missing data
        for (i in 1:N) {
            if (pos_obs[i] && pm_obs[i]) {
                pos_pm_measured[i] ~ multi_normal_cholesky(pos_pm[1:4,i], pos_pm_cov_mats[i]);
            } else {
                if (pos_obs[i]) {
                    pos_measured[i] ~ multi_normal_cholesky(pos_pm[1:2,i], pos_cov_mats[i]);
                }
                if (pm_obs[i]) {
                    pm_measured[i] ~ multi_normal_cholesky(pos_pm[3:4,i], pm_cov_mats[i]);
                }
            }
            if (dist_obs[i]) {
                dist_std[i] ~ std_normal();
            }
            if (vlos_obs[i]) {
                vlos_std[i] ~ std_normal();
            }
        }

        // using complete data here. uncomment this as necessary (might be a little faster).

        /* for (i in 1:N) { */
        /*     pos_pm_measured[i] ~ multi_normal_cholesky(pos_pm[1:4,i], pos_pm_cov_mats[i]); */
        /* } */
        /* dist_std ~ std_normal(); */
        /* vlos_std ~ std_normal(); */

    }

    // likelihood ---------------------

    profile("df") {
        // single core non-vectorized
        /* target += partial_df_lupdf(dummy | 1, N, y, p_phi0, p_gamma, p_alpha, p_beta); */

        // multicore non-vectorized
        /* target += reduce_sum(partial_df_lupdf, dummy, grainsize, y, p_phi0, p_gamma, p_alpha, p_beta); */

        // single core vectorized. this is fastest.
        y ~ df_vec(p_phi0, p_gamma, p_alpha, p_beta);

        // multicore vectorize
        /* target += reduce_sum(partial_df_vec_lupdf, dummy, grainsize, y, p_phi0, p_gamma, p_alpha, p_beta); */
    }
}

generated quantities {
}

