functions {
  real df_lpdf(array[] real y, real phi0, real g, real b, real a) {
    // y[1] = v, y[2] = v_t, y[3] = r;
    real E = -square(y[1]) / 2. + phi0 / pow(y[3], g);
    real L = y[3] * y[2];
    
    // terms in numerator
    real num_t1 = -2 * b * log(L);
    real num_t2 = (b * (g - 2.) / g + a / g - 3. / 2.) * log(E);
    real num_t3 = lgamma(a / g - 2. * b / g + 1.);
    
    // terms in denominator
    real denom_t1 = log(pi() * sqrt(pi()) * pow(2., -b + 3. / 2.))
                    + lgamma(1. - b);
    real denom_t2 = (-2. * b / g + a / g) * log(phi0);
    real denom_t3 = lgamma(b * (g - 2.) / g + a / g - 1. / 2.);
    
    real numerator = num_t1 + num_t2 + num_t3;
    real denominator = denom_t1 + denom_t2 + denom_t3;
    
    // this ~should~ never happen
    if (E < 0.) {
      return negative_infinity();
    }
    
    return numerator - denominator;
  }
  
  real shifted_gamma_lpdf(real y, real a, real b, real shift) {
    real shifted_y = y - shift;
    return a * log(b) - lgamma(a) + (a - 1) * log(shifted_y) - b * shifted_y;
  }
  
  vector get_angles(real x, real y, real z) {
    real projR = sqrt(square(x) + square(y));
    real R = sqrt(square(x) + square(y) + square(z));
    
    real sintheta = y / projR;
    real costheta = x / projR;
    
    real sinphi = projR / R;
    real cosphi = z / R;
    
    return to_vector([sintheta, costheta, sinphi, cosphi]);
  }
  
  vector transform_vels(real ra, real dec, real pmra, real pmdec, real vlos,
                        real plx, real sintheta, real costheta, real sinphi,
                        real cosphi) {
    matrix[3, 3] T = [[-0.06699, -0.87276, -0.48354],
                      [0.49273, -0.45035, 0.74458],
                      [-0.86760, -0.18837, 0.46020]];
    real deg2rad = pi() / 180.;
    real sinra = sin(ra * deg2rad);
    real cosra = cos(ra * deg2rad);
    real sindec = sin(dec * deg2rad);
    real cosdec = cos(dec * deg2rad);
    matrix[3, 3] A = [[cosra * cosdec, -sinra, -cosra * sindec],
                      [sinra * cosdec, cosra, -sinra * sindec],
                      [sindec, 0, cosdec]];
    matrix[3, 3] B = T * A;
    real k = 4.74057;
    vector[3] solarmotion = to_vector([11.1, 232.24, 7.25]); // including rotation
    
    vector[3] dat = to_vector([vlos, k * pmra / plx, k * pmdec / plx]);
    vector[3] uvw = B * dat + solarmotion;
    
    matrix[3, 3] ptz_mat = [[costheta, sintheta, 0],
                            [-sintheta, costheta, 0], [0, 0, 1]];
    vector[3] ptz = ptz_mat * uvw;
    
    matrix[3, 3] rtp_mat = [[cosphi, 0, sinphi], [0, 1, 0],
                            [-sinphi, 0, cosphi]];
    
    vector[3] rtp = rtp_mat * ptz;
    
    return rtp / 100; // in units of 100km/s
  }
}
data {
  int<lower=1> N; // total length of data
  
  vector[N] ra;
  vector[N] dec;
  vector[N] plx;
  
  vector[N] Xgc;
  vector[N] Ygc;
  vector[N] Zgc;
  
  // measurements
  vector[N] pmra_measured;
  vector[N] pmdec_measured;
  vector[N] vlos_measured;
  vector[N] r_measured;
  
  // measurement uncertainties
  vector[N] pmra_err;
  vector[N] pmdec_err;
  vector[N] vlos_err;
  vector[N] r_err;
}
transformed data {
  vector[N] sintheta;
  vector[N] costheta;
  vector[N] sinphi;
  vector[N] cosphi;
  
  for (i in 1 : N) {
    vector[4] angles = get_angles(Xgc[i], Ygc[i], Zgc[i]);
    sintheta[i] = angles[1];
    costheta[i] = angles[2];
    sinphi[i] = angles[3];
    cosphi[i] = angles[4];
  }
}
parameters {
  real<lower=0., upper=1.> p_gamma; // power-law slope of gravitational potential
  real<lower=1., upper=200.> p_phi0; // scale factor for gravitational potential */
  real<lower=-3., upper=1.> p_beta; // velocity anisotropy parameter
  real<lower=max(to_vector([3., p_beta * (2. - p_gamma) + p_gamma / 2.])),
       upper=10.> p_alpha; // power-law slope of satellite population
  
  vector<lower=min(pmra_measured - 5 * pmra_err),
         upper=max(pmra_measured + 5 * pmra_err)>[N] pmra;
  
  vector<lower=min(pmdec_measured - 5 * pmdec_err),
         upper=max(pmdec_measured + 5 * pmdec_err)>[N] pmdec;
  
  vector<lower=min(vlos_measured - 5 * vlos_err),
         upper=max(vlos_measured + 5 * vlos_err)>[N] vlos;
  vector<lower=0, upper=max(r_measured + 5 * r_err)>[N] r;
}
transformed parameters {
  
}
model {
  // hyperpriors ----------------
  
  // for the GC dataset
  p_gamma ~ normal(0.5, 0.06);
  p_phi0 ~ uniform(1., 200);
  p_beta ~ uniform(-0.5, 1.);
  p_alpha ~ shifted_gamma(2.99321604, 2.82409927, 3.);
  
  // observation process ----------------------
  
  pmra_measured ~ normal(pmra, pmra_err);
  pmdec_measured ~ normal(pmdec, pmdec_err);
  vlos_measured ~ normal(vlos, vlos_err);
  r_measured ~ normal(r, r_err);
  
  // likelihood ---------------------
  
  for (i in 1 : N) {
    vector[3] vels_sph = transform_vels(ra[i], dec[i], pmra[i], pmdec[i],
                                        vlos[i], plx[i], sintheta[i],
                                        costheta[i], sinphi[i], cosphi[i]);
    array[3] real y;
    y[2] = sqrt(square(vels_sph[2]) + square(vels_sph[3]));
    y[1] = sqrt(square(vels_sph[1]) + square(y[2]));
    y[3] = r[i];
    y ~ df(p_phi0, p_gamma, p_beta, p_alpha);
  }
}
generated quantities {
  
}

