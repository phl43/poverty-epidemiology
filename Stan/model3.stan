data {
    int<lower=1> N; // length of data
    int<lower=0> t[N]; // times of data
    int<lower=1> P; // number of covariates
    int<lower=1> seed_length; // number of days for seed
    int<lower=1> prediction_length; // length of prediction
    matrix[seed_length + prediction_length, P] X; // features matrix
    int<lower=0> C_h[N]; // incidence for all variants except B.1.1.7
    int<lower=0> C_v[N]; // incidence for B.1.1.7
    real<lower=0> mu; // mean of generation time
    real<lower=0> sigma; // coefficient of variation of generation time
}

transformed data {
    int total_length = seed_length + prediction_length;
    vector[total_length] times;
    real cv = sigma / mu;
    real alpha = 1 / cv^2;
    real beta = 1 / (cv^2 * mu);
    vector[total_length] gi;
    vector[total_length] gi_rev;

    // compute generation time distribution
    for (i in 1:total_length) {
        times[i] = i;
        gi[i] = gamma_cdf(i, alpha, beta) - gamma_cdf(i - 1, alpha, beta);
    }

    // reverse generation time for convolution
    for(i in 1:total_length) {
        gi_rev[i] = gi[total_length - i + 1];
    }
}

parameters {
    real<lower=0> seed_h;
    real<lower=0> seed_v;
    real<lower=0> R0_h;
    real<lower=0> R0_v;
    vector[P] gamma;
    real<lower=0> phi_h;
    real<lower=0> phi_v;
}

transformed parameters {
    vector[total_length] Rt_h;
    vector[total_length] Rt_v;
    real<lower=0> Rt[total_length] = rep_array(.0, total_length);
    vector[total_length] expected_C_h = append_row(rep_vector(seed_h, seed_length), rep_vector(0, prediction_length));
    vector[total_length] expected_C_v= append_row(rep_vector(seed_v, seed_length), rep_vector(0, prediction_length));

    Rt_h = R0_h * exp(-X * gamma);
    Rt_v = R0_v * exp(-X * gamma);

    for (i in (seed_length + 1):total_length) {
        real convolution_h = dot_product(head(expected_C_h, i - 1), tail(gi_rev, i - 1));
        real convolution_v = dot_product(head(expected_C_v, i - 1), tail(gi_rev, i - 1));

        expected_C_h[i] = Rt_h[i] * convolution_h;
        expected_C_v[i] = Rt_v[i] * convolution_v;

        Rt[i] = Rt_h[i] * convolution_h / (convolution_h + convolution_v) + Rt_v[i] * convolution_v / (convolution_h + convolution_v);
    }
}

model {
    seed_h ~ normal(10000, 5000);
    seed_v ~ normal(100, 100);
    R0_h ~ uniform(0, 5);
    R0_v ~ uniform(0, 5);
    gamma ~ normal(0, 1);
    phi_h ~ normal(0,5);
    phi_v ~ normal(0,5);
    for (i in 1:N) {
        C_h[i] ~ neg_binomial_2(expected_C_h[t[i]], phi_h);
    }
    for (i in 1:N) {
        C_v[i] ~ neg_binomial_2(expected_C_v[t[i]], phi_v);
    }
}
