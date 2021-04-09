data {
    int<lower=1> N; // length of data
    int<lower=0> t[N]; // times of data
    int<lower=0> N_t[N]; // number of cases for each data point
    int<lower=0> C_v[N]; // number of B.1.1.7 cases for each data point
    real<lower=0> mu; // mean of generation time
    real<lower=0> sigma; // standard deviation of generation time
    real<lower=0> R_h; // reproduction number of the historical lineage
    int<lower=1> prediction_length; // length of prediction
}

transformed data {
    real cv;
    real a;
    real b;
    real r_h;

    cv = sigma / mu;
    a = 1 / cv^2;
    b = a / mu;
    r_h = b * (pow(R_h, 1 / a) - 1);
}

parameters {
    real<lower=0> p_0;
    real alpha;
}

transformed parameters {
    real<lower=0, upper=1> p_t[N];
    real<lower=0> R_v = (1 + alpha) * R_h;
    real r_v = b * (pow(R_v, 1 / a) - 1);

    for (i in 1:N) {
        p_t[i] = 1 / (1 + (1 - p_0) / p_0 * exp((r_h - r_v) * t[i]));
    }
}

model {
    p_0 ~ uniform(0, 1);
    alpha ~ uniform(-1, 1);
    C_v ~ binomial(N_t, p_t);
}

generated quantities {
    real<lower=0> expected_p_t[prediction_length];
    real<lower=0> expected_C_h[prediction_length];
    real<lower=0> expected_C_v[prediction_length];

    for (i in 1:prediction_length) {
        expected_p_t[i] = 1 / (1 + (1 - p_0) / p_0 * exp((r_h - r_v) * i));

        expected_C_h[i] = (N_t[1] - C_v[1]) * exp(r_h * (i - t[1] + 1));
        expected_C_v[i] = C_v[1] * exp(r_v * (i - t[1] + 1));
    }
}
