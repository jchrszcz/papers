data {
  int<lower=1> N;
  int<lower=1> J;
  int<lower=1> K;
  int condition[N];
  int idx[N];
  vector[N] prob;
  vector[N] prompt;
  vector[3] a;
}
parameters {
  simplex[3] phi;
  real<lower=0> mu[J];
  real<lower=0> sigma[J];
  matrix[K,J] mu_idx;
  matrix<lower=0>[K,J] sigma_idx;
  real<lower=0> sigma_mu[J];
  real<lower=0> sigma_sigma[J];
}
transformed parameters {
  vector[N] y;
  vector[N] alpha;
  vector[N] beta;
  vector[K] e_mu;
  vector[K] e_sigma;
  
  for (n in 1:N) {
    y[n] = Phi((prompt[n] - mu_idx[idx[n], condition[n]]) / sigma_idx[idx[n], condition[n]]);
    alpha[n] = sigma_idx[idx[n], condition[n]] * y[n];
    beta[n] = sigma_idx[idx[n], condition[n]] * (1.0 - y[n]);
  }

  for (k in 1:K) {
    e_mu[k] = mu_idx[k,1] - mu[1];
    e_sigma[k] = sigma_idx[k,1] - sigma[1];
  }
}
model {
  phi ~ dirichlet(a);
  mu ~ normal(0,1);
  sigma ~ normal(0,1);
  sigma_mu ~ normal(0,1);
  sigma_sigma ~ normal(0,1);
  for (k in 1:K) {
    for (j in 1:J) {
      mu_idx[k,j] ~ normal(mu[j], sigma_mu[j]);
      sigma_idx[k,j] ~ normal(sigma[j], sigma_sigma[j]);
    }
  }
  for (n in 1:N) {
    if (prob[n] == 0) {
      target += bernoulli_lpmf(1 | phi[1]); 
    }
    else if (prob[n] == 1) {
      target += bernoulli_lpmf(1 | phi[2]); 
    }
    else {
      target += bernoulli_lpmf(1 | phi[3]) + beta_lpdf(prob[n] | alpha[n], beta[n]); 
    }
  }
}