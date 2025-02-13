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
}

parameters {
  real beta;
  real<lower=0> phi;
  real offset;
  real d[N];
}

transformed parameters {
  real<lower=0> theta = phi;
  real pbar[N];
  for (n in 1:N){
    pbar[n] = inv_logit(offset + beta*d[n]);
  }
}

model {
  beta ~ normal(0, 2); // decline rate param
  d_meas ~ normal(d, d_std);
  d ~ normal(0, 4);
  phi ~ exponential(1);
  
  for (n in 1:N){
    target += beta_binomial_lpmf(trans[n] |  HH[n], theta*pbar[n]/(1-pbar[n]), theta);
  }
}

generated quantities {
  vector[51] x_gen;

  int c;
  
  c = 0;
   for (tr in -20:30){
     c = c+1;
     x_gen[c] = inv_logit(offset + beta*tr/10);
   }
   
   real log_lik[N];

   for (n in 1:N){
   log_lik[n] = beta_binomial_lpmf(trans[n] |  HH[n], theta*pbar[n]/(1-pbar[n]), theta);
   }
}
