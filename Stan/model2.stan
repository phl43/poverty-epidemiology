data {
    int<lower=1> N; // length of data
    int<lower=0> t[N]; // times of data
    int<lower=0> C_h[N]; // number of non-B.1.1.7 cases for each data point
    int<lower=0> C_v[N]; // number of B.1.1.7 cases for each data point
    real<lower=0> mu; // mean of generation time
    real<lower=0> sigma; // standard deviation of generation time
    int<lower=1> prediction_length; // length of prediction
}

transformed data {
    real cv;
    real a;
    real b;

    cv = sigma / mu;
    a = 1 / cv^2;
    b = a / mu;
}

parameters {
    real<lower=0> seed_h;
    real<lower=0> seed_v;
    real<lower=0> R_h;
    real<lower=0> R_v;
    real<lower=0> phi_h;
    real<lower=0> phi_v;
}

transformed parameters {
    real r_h = b * (pow(R_h, 1 / a) - 1);
    real r_v = b * (pow(R_v, 1 / a) - 1);
    real<lower=0> expected_C_h[N];
    real<lower=0> expected_C_v[N];

    for (i in 1:N) {
        expected_C_h[i] = seed_h * exp(r_h * t[i]);
        expected_C_v[i] = seed_v * exp(r_v * t[i]);
    }
}

model {
    seed_h ~ normal(15000, 5000);
    seed_v ~ normal(100, 500);
    R_h ~ normal(1, 5);
    R_v ~ normal(1, 5);
    phi_h ~ normal(0,5);
    phi_v ~ normal(0,5);
    for (i in 1:N) {
        C_h[i] ~ neg_binomial_2(expected_C_h[i], phi_h);
    }
    for (i in 1:N) {
        C_v[i] ~ neg_binomial_2(expected_C_v[i], phi_v);
    }
}

generated quantities {
    real<lower=0> predicted_C_h[prediction_length];
    real<lower=0> predicted_C_v[prediction_length];

    for (i in 1:prediction_length) {
        predicted_C_h[i] = seed_h * exp(r_h * i);
        predicted_C_v[i] = seed_v * exp(r_v * i);
    }
}
