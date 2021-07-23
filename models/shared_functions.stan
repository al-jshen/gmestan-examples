functions {

    real df_lpdf(real[] y, real phi0, real g, real a, real b) {

        // y[1] = v, y[2] = v_t, y[3] = r;
        real E = -square(y[1]) / 2. + phi0 / pow(y[3], g);
        real L = y[3] * y[2];

        // terms in numerator
        real num_t1 = -2 * b * log(L);
        real num_t2 = (b * (g - 2.) / g + a / g - 3. / 2.) * log(E);
        real num_t3 = lgamma(a / g - 2. * b / g + 1.);

        // terms in denominator
        real denom_t1 = log(pi() * sqrt(pi()) * pow(2., -b + 3. / 2.)) + lgamma(1. - b);
        real denom_t2 = (-2. * b / g + a / g) * log(phi0);
        real denom_t3 = lgamma(b * (g - 2.) / g + a / g - 1. / 2.);

        real numerator = num_t1 + num_t2 + num_t3;
        real denominator = denom_t1 + denom_t2 + denom_t3;

        /* print("E ", E); */

        // this ~should~ never happen
        if (E < 0.) {
            return negative_infinity();
        }

        return numerator - denominator;
    }

    real df_vec_lpdf(matrix y, real phi0, real g, real a, real b) {

        int N = size(y[1]);
        // y has shape [3, N], each column contains v, v_t, r, in that order, N observations (columns)
        row_vector[N] E = -square(y[1]) ./ 2. + phi0 ./ pow(y[3], g);
        row_vector[N] L = y[3] .* y[2];

        // save some computations
        real ag = a / g;
        real bg = b / g;
        real bg2gag = b * (g - 2) / g + a / g;

        // terms in numerator
        row_vector[N] num_t1 = -2 * b * log(L);
        row_vector[N] num_t2 = (bg2gag - 3. / 2.) * log(E);
        real num_t3 = lgamma(ag - 2. * bg  + 1.);

        // terms in denominator
        real denom_t1 = log(pi() * sqrt(pi()) * pow(2., -b + 3. / 2.)) + lgamma(1. - b);
        real denom_t2 = (-2. * bg + ag) * log(phi0);
        real denom_t3 = lgamma(bg2gag - 1. / 2.);

        row_vector[N] numerator = num_t1 + num_t2 + num_t3;
        real denominator = denom_t1 + denom_t2 + denom_t3;

        return sum(numerator - denominator);
    }

    real partial_df_lpdf(real[] dummy, int start, int end, real[,] slice_y, real phi0, real g, real a, real b) {
        // for N = end - start, slice_y is an [N, 3] array of reals
        // y[i] is a length 3 array of reals, containing
        // y[1] = v, y[2] = v_t, y[3] = r
        real sum = 0.;
        for (i in start:end) {
            sum += df_lupdf(slice_y[i] | phi0, g, a, b);
        }
        return sum;
    }

    real partial_df_vec_lpdf(real[] dummy, int start, int end, matrix slice_y, real phi0, real g, real a, real b) {
        return df_vec_lupdf(slice_y[,start:end] | phi0, g, a, b);
    }

    real shifted_gamma_lpdf(real y, real a, real b, real shift) {
        real shifted_y = y - shift;
        return a * log(b) - lgamma(a) + (a - 1) * log(shifted_y) - b * shifted_y;
    }


    // ra, dec in radians
    // dist in kpc
    // R is icrs -> gctc matrix
    // H is icrs -> gctc matrix
    // offset is icrs -> gctc offset vector (galcen_dist, 0, 0)
    vector transform_pos(real ra, real dec, real dist, matrix R, matrix H, vector offset) {

        real x_icrs = dist * cos(ra) * cos(dec);
        real y_icrs = dist * sin(ra) * cos(dec);
        real z_icrs = dist * sin(dec);
        vector[3] r_icrs = to_vector([x_icrs, y_icrs, z_icrs]);

        vector[3] r_gal = R * r_icrs;
        r_gal -= offset;
        r_gal = H * r_gal;

        /* print("position transformation"); */
        /* print("ra ", ra, "dec ", dec, "dist ", dist); */
        /* print("r_gal ", r_gal); */
        return r_gal;
    }

    // same as transform_pos, but vectorized. will return length-N array of length-3 vectors
    matrix transform_pos_vec(row_vector ra, row_vector dec, row_vector dist, matrix R, matrix H) {
        // number of tracers
        int N = size(ra);

        row_vector[N] dist_cosdec = dist .* cos(dec);
        row_vector[N] x_icrs = dist_cosdec .* cos(ra);
        row_vector[N] y_icrs = dist_cosdec .* sin(ra);
        row_vector[N] z_icrs = dist .* sin(dec);

        // append_row is like np.stack, just using [] is even faster
        matrix[3, N] r_icrs = [x_icrs, y_icrs, z_icrs]; // = append_row(x_icrs, append_row(y_icrs, z_icrs));

        matrix[3, N] r_gal = R * r_icrs;
        // get the first row (x_gal) and apply offset
        r_gal[1] -= 8.122;
        r_gal = H * r_gal;
        return r_gal;
    }

    // converts cartesian position to spherical position
    vector c2s_pos(real x, real y, real z) {
        real dist = sqrt(dot_self([x, y, z]));
        real theta = asin(z / dist);
        real phi = atan2(y, x);

        return to_vector([dist, theta, phi]);
    }

    matrix c2s_pos_vec(row_vector x, row_vector y, row_vector z) {
        int N = size(x);

        matrix[3, N] cart_pos = [x, y, z];

        row_vector[N] dist = sqrt(columns_dot_self(cart_pos));
        row_vector[N] theta = asin(z ./ dist);
        row_vector[N] phi;

        // atan2 not vectorized
        for (i in 1:N) {
            phi[i] = atan2(y[i], x[i]);
        }

        return [dist, theta, phi];

    }

    // converts spherical position to cartesian position
    vector s2c_pos(real r, real theta, real phi) {
        real x = r * cos(phi) * cos(theta);
        real y = r * sin(phi) * cos(theta);
        real z = r * sin(theta);

        return to_vector([x, y, z]);
    }

    matrix s2c_pos_vec(row_vector r, row_vector theta, row_vector phi) {
        int N = size(r);
        row_vector[N] r_costheta = r .* cos(theta);
        row_vector[N] x = r_costheta .* cos(phi);
        row_vector[N] y = r_costheta .* sin(phi);
        row_vector[N] z = r .* sin(theta);

        return [x, y, z];

    }

    // converts cartesian velocity to spherical velocity
    vector c2s_vel(real x, real y, real z, real vx, real vy, real vz) {

        vector[3] sph_pos = c2s_pos(x, y, z);
        real dist = sph_pos[1];
        real lat = sph_pos[2];
        real lon = sph_pos[3];
        real proj_dist = sqrt(square(x) + square(y));

        real vr = dot_product([x, y, z], [vx, vy, vz]) / dist;

        real mu_theta = (z * (x * vx + y * vy) - square(proj_dist) * vz) / square(dist) / proj_dist;
        real vtheta = -mu_theta * dist;

        real mu_phi = (x * vy - y * vx) / square(proj_dist);
        real vphi = mu_phi * dist * cos(lat);

        return to_vector([vr, vtheta, vphi]);
    }

    matrix c2s_vel_vec(row_vector x, row_vector y, row_vector z, row_vector vx, row_vector vy, row_vector vz) {
        int N = size(x);

        matrix[3, N] sph_pos = c2s_pos_vec(x, y, z);

        row_vector[N] dist = sph_pos[1];
        row_vector[N] lat = sph_pos[2];
        row_vector[N] lon = sph_pos[3];
        row_vector[N] proj_dist_sq = square(x) + square(y);

        matrix[3, N] cart_pos = [x, y, z]; // = append_row(x, append_row(y, z));
        matrix[3, N] cart_vel = [vx, vy, vz]; // = append_row(vx, append_row(vy, vz));

        row_vector[N] vr = columns_dot_product(cart_pos, cart_vel) ./ dist;

        row_vector[N] mu_theta = (z .* (x .* vx + y .* vy) - proj_dist_sq .* vz) ./ square(dist) ./ sqrt(proj_dist_sq);
        row_vector[N] vtheta = -mu_theta .* dist;

        row_vector[N] mu_phi = (x .* vy - y .* vx) ./ proj_dist_sq;
        row_vector[N] vphi = mu_phi .* dist .* cos(lat);

        return [vr, vtheta, vphi];

    }

    // converts spherical velocity to cartesian velocity
    vector s2c_vel(real r, real theta, real phi, real vr, real vtheta, real vphi) {

        real vx = vr * cos(phi) * cos(theta) - vphi * sin(phi) - vtheta * cos(phi) * sin(theta);
        real vy = vr * sin(phi) * cos(theta) + vphi * cos(phi) - vtheta * sin(phi) * sin(theta);
        real vz = vr * sin(theta) + vtheta * cos(theta);

        return to_vector([vx, vy, vz]);
    }

    matrix s2c_vel_vec(row_vector r, row_vector theta, row_vector phi, row_vector vr, row_vector vtheta, row_vector vphi) {
        int N = size(r);

        // compute once and reuse
        row_vector[N] sintheta = sin(theta);
        row_vector[N] costheta = cos(theta);
        row_vector[N] sinphi = sin(phi);
        row_vector[N] cosphi = cos(phi);

        row_vector[N] vx = vr .* cosphi .* costheta - vphi .* sinphi - vtheta .* cosphi .* sintheta;
        row_vector[N] vy = vr .* sinphi .* costheta + vphi .* cosphi - vtheta .* sinphi .* sintheta;
        row_vector[N] vz = vr .* sintheta + vtheta .* costheta;

        return [vx, vy, vz];
    }


    // converts heliocentric velocity to galactocentric velocity
    // ra, dec in radians
    // dist in kpc
    // pmra, pmdec in mas/yr
    // vlos in km/s
    vector transform_vels(real ra, real dec, real dist, real pmra, real pmdec, real vlos, matrix R, matrix H, vector offset, vector solarmotion) {
        real km_per_kpc = 3.085677581491367e+16;
        real mas_to_unitless = 4.84813681109536e-09;
        real s_per_yr = 31557600.0;

        real vra = dist * pmra * km_per_kpc * mas_to_unitless / s_per_yr; // in km/s
        real vdec = dist * pmdec * km_per_kpc * mas_to_unitless / s_per_yr; // in km/s

        vector[3] r_gal = transform_pos(ra, dec, dist, R, H, offset);
        vector[3] v_icrs = s2c_vel(dist, dec, ra, vlos, vdec, vra);

        matrix[3, 3] A = H * R;
        vector[3] v_gal = A * v_icrs + solarmotion;

        vector[3] v_gal_sph = c2s_vel(r_gal[1], r_gal[2], r_gal[3], v_gal[1], v_gal[2], v_gal[3]);
        /* print("velocity transformation"); */
        /* print("ra ", ra, "dec ", dec, "dist ", dist, "pmra ", pmra, "pmdec ", pmdec, "vlos ", vlos); */
        /* print("v_gal_sph ", v_gal_sph); */
        return v_gal_sph;
    }

    // r_gal is [3, N] matrix with rows xgc, ygc, zgc
    matrix transform_vels_vec(row_vector ra, row_vector dec, row_vector dist, row_vector pmra, row_vector pmdec, row_vector vlos, matrix R, matrix H, vector offset, vector solarmotion) {

        int N = size(ra);

        /* this is km_per_kpc * mas_to_unitless / s_per_yr */
        real conversion_factor = 4.740470463533349;
        row_vector[N] dist_with_conversion = dist * conversion_factor;

        row_vector[N] vra = dist_with_conversion .* pmra; // km/s
        row_vector[N] vdec = dist_with_conversion .* pmdec; // km/s

        matrix[3, N] r_gal = transform_pos_vec(ra, dec, dist, R, H);
        matrix[3, N] v_icrs = s2c_vel_vec(dist, dec, ra, vlos, vdec, vra);

        matrix[3, 3] A = H * R;
        matrix[3, N] v_gal = A * v_icrs;
        v_gal[1] += solarmotion[1];
        v_gal[2] += solarmotion[2];
        v_gal[3] += solarmotion[3];

        matrix[3, N] v_gal_sph = c2s_vel_vec(r_gal[1], r_gal[2], r_gal[3], v_gal[1], v_gal[2], v_gal[3]);
        return v_gal_sph;
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


    vector transform_vels_old(real ra, real dec, vector pm, real vlos, real plx, real sintheta, real costheta, real sinphi, real cosphi) {
        matrix[3,3] T = [
            [-0.06699, -0.87276, -0.48354],
            [ 0.49273, -0.45035,  0.74458],
            [-0.86760, -0.18837,  0.46020]
        ];
        real deg2rad = pi() / 180.;
        real sinra = sin(ra * deg2rad);
        real cosra = cos(ra * deg2rad);
        real sindec = sin(dec * deg2rad);
        real cosdec = cos(dec * deg2rad);
        real pmra = pm[1];
        real pmdec = pm[2];
        matrix[3,3] A = [
            [cosra * cosdec, -sinra, -cosra * sindec],
            [sinra * cosdec, cosra , -sinra * sindec],
            [sindec,              0,          cosdec]
        ];
        matrix[3, 3] B = T * A;
        real k = 4.74057;
        vector[3] solarmotion = to_vector([11.1, 232.24, 7.25]); // including rotation

        vector[3] dat = to_vector([vlos, k * pmra / plx, k * pmdec / plx]);
        vector[3] uvw = B * dat + solarmotion;

        matrix[3, 3] ptz_mat = [
            [costheta, sintheta, 0],
            [-sintheta, costheta, 0],
            [0, 0, 1]
        ];
        vector[3] ptz = ptz_mat * uvw;

        matrix[3, 3] rtp_mat = [
            [cosphi, 0, sinphi],
            [0, 1, 0],
            [-sinphi, 0, cosphi]
        ];

        vector[3] rtp = rtp_mat * ptz;

        return rtp / 100; // in units of 100km/s
    }

}
