data {
  int<lower=0> N;
  vector[N] d_meas;
  vector[N] d_std;
  int nvoc;
  int<lower=1, upper=nvoc> voc[N];
  int<lower=0> trans[N];
  int<lower=1> HH[N];
  int<lower=1> n_ind_age;
  int<lower=1, upper=n_ind_age> index_age[N];
  real<lower=0> crowding[N]; 
  int<lower=0, upper=1> comorb[N];
}

parameters {
  real beta;
  real gamma_ind[n_ind_age];
  real kappa;
  real epsilon;
  real<lower=0> phi;
  real d[N];
  real zeta_mu[nvoc];
  real<lower=0.3> hp_sd[nvoc];
  real zeta_sd[nvoc];
  real hp_mu_mu;
  real<lower=0.3> hp_mu_sd;
  real hp_sd_mu;
  real<lower=0.3> hp_sd_sd;
}

transformed parameters {
  real<lower=0> theta = phi;
  real zeta[nvoc];
  real pbar[N];
  for (n in 1:nvoc){
     zeta[n] = zeta_mu[n] + hp_sd[n]*zeta_sd[n];
  }
  for (n in 1:N){
    pbar[n] = inv_logit(zeta[voc[n]]+gamma_ind[index_age[n]]+beta*d[n]+kappa*crowding[n]+epsilon*comorb[n]);
    
  }
}

model {
  beta ~ normal(0, 2); // decline rate param
  gamma_ind ~ normal(0, 2);
  zeta_sd ~ std_normal();
  zeta_mu ~ normal(hp_mu_mu, hp_mu_sd);
  hp_mu_mu ~ normal(0, 2);
  hp_mu_sd ~ normal(0, 2);
  hp_sd ~ normal(hp_sd_mu, hp_sd_sd);
  hp_sd_mu ~ normal(0, 2);
  hp_sd_sd ~ normal(0, 2);
  kappa ~ normal(0, 2);
  epsilon ~ normal(0, 2);
  d_meas ~ normal(d, d_std);
  d ~ normal(0, 4);
  phi ~ exponential(1);
  
  for (n in 1:N){
    target += beta_binomial_lpmf(trans[n] |  HH[n], theta*pbar[n]/(1-pbar[n]), theta);
  }
}

generated quantities {
  vector[51] x_gen_pre_alpha;
  vector[51] x_gen_alpha;
  vector[51] x_gen_delta;

  int c;
  
  c = 0;
   for (tr in -20:30){
     c = c+1;
     x_gen_pre_alpha[c] = inv_logit(zeta[1] + gamma_ind[2] + beta*tr/10 + kappa);
     x_gen_alpha[c] = inv_logit(zeta[2] + gamma_ind[2] + beta*tr/10 + kappa);
     x_gen_delta[c] = inv_logit(zeta[3] + gamma_ind[2] + beta*tr/10 + kappa);
   }
   
   real log_lik[N];
   for (n in 1:N){
   log_lik[n] = beta_binomial_lpmf(trans[n] |  HH[n], theta*pbar[n]/(1-pbar[n]), theta);
   }
}
