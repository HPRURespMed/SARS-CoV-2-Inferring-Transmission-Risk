data {
  int<lower=0> N;
  vector[N] d_meas;
  vector[N] d_std;
  int nvoc;
  int<lower=0, upper=nvoc> voc[N];
  int<lower=0> trans_vacc[N];
  int<lower=0> trans_unvacc[N];
  int<lower=0> HH_vacc[N];
  int<lower=0> HH_unvacc[N];
  int<lower=1> n_ind_age;
  int<lower=1, upper=n_ind_age> index_age[N];
  //real<lower=0> crowding[N];
  int<lower=-1, upper=1> ind_vacc[N];
  //int con_vacc[N];
  int<lower=0, upper=1> antiviral[N];
  int<lower=0, upper=1> chronic[N];
  //int<lower=1> age[N];
}

// The parameters accepted by the model.
parameters {
  real alpha;
  real beta;
  //real kappa;
  // real zeta[nvoc];
  real d[N];
  // real zeta_mu[nvoc];
  // real<lower=0.1> hp_sd[nvoc];
  // real zeta_sd[nvoc];
  // real hp_mu_mu;
  // real<lower=0.1> hp_mu_sd;
  // real hp_sd_mu;
  // real<lower=0.1> hp_sd_sd;
  // real gamma_ind[n_ind_age];
  // real delta;
  real<lower=1> theta;
  //real<lower=1> theta_unvacc;

  //real gamma;
  //real sigma;
}

transformed parameters {
  //real zeta[nvoc];
  //real<lower=0, upper=1> pbar_vacc[N];
  real<lower=0, upper=1> pbar[N];
  // for (n in 1:nvoc){
  //   zeta[n] = zeta_mu[n] + hp_sd[n]*zeta_sd[n];
  // }
  
  for (n in 1:N){
    //pbar_vacc[n] = inv_logit(beta*d[n] + gamma_ind[index_age[n]] + zeta[voc[n]] + kappa*antiviral[n]+delta*chronic[n]);
    pbar[n] = inv_logit(alpha + beta*d[n]);

  }
}

// The model to be estimated.
model {
  alpha ~ normal(0, 2);
  beta ~ normal(0, 2); // decline rate param
  // gamma_ind ~ normal(0, 2);
  // kappa ~ normal(0, 2);
  // delta ~ normal(0, 2);
  // theta_vacc ~ exponential(1);
  theta ~ exponential(1);

  // zeta ~ normal(0, 2); // voc param
  
   // Hierarchical non-centered voc param
  // zeta_sd ~ std_normal();
  // zeta_mu ~ normal(hp_mu_mu, hp_mu_sd);
  // hp_mu_mu ~ normal(0, 2);
  // hp_mu_sd ~ normal(0, 2);
  // hp_sd ~ normal(hp_sd_mu, hp_sd_sd);
  // hp_sd_mu ~ normal(0, 2);
  // hp_sd_sd ~ normal(0, 2);
  // kappa ~ normal(0, 2);
  d_meas ~ normal(d, d_std);
  d ~ normal(0, 2);
  //sigma ~ normal(0, 2);
  for (n in 1:N){
    //trans[n] ~ bernoulli_logit(alpha + beta*d[n] + zeta[voc[n]]);
    //trans[n] ~ binomial_logit(HH[n], beta*d[n] + gamma*age[n] + zeta[voc[n]]);
    //trans_unvacc[n] ~ binomial_logit(HH_unvacc[n], alpha + beta*d[n] + gamma_ind[index_age[n]] + zeta[voc[n]] + kappa*antiviral[n]+delta*chronic[n]);
    //trans_vacc[n] ~ binomial_logit(HH_vacc[n], beta*d[n] + gamma_ind[index_age[n]] + zeta[voc[n]] + kappa*antiviral[n]+delta*chronic[n]);
    //trans_vacc[n] ~ binomial_logit(HH_vacc[n], beta*d[n] + gamma_ind[index_age[n]] + zeta[voc[n]] + kappa*antiviral[n]+delta*chronic[n]);
    //trans_unvacc[n] ~ binomial_logit(HH_unvacc[n], alpha + beta*d[n] + gamma_ind[index_age[n]] + zeta[voc[n]] + kappa*antiviral[n]+delta*chronic[n]);
    //target += beta_binomial_lpmf(trans_vacc[n] |  HH_vacc[n], theta_vacc*pbar_vacc[n]/(1-pbar_vacc[n]), theta_vacc);
    target += beta_binomial_lpmf(trans_unvacc[n] + trans_vacc[n] |  HH_vacc[n]+ HH_unvacc[n], theta*pbar[n]/(1-pbar[n]), theta);

  }
}

generated quantities {
  //vector[51] x_gen_pre_alpha;
  //vector[51] x_gen_alpha;
  //vector[51] x_gen_delta;
  vector[56] x_gen;

  int c;
  
  c = 0;
   for (tr in -20:35){
     c = c+1;
     // Vaccinated Delta
     x_gen[c] = inv_logit( alpha + beta*tr/10);
     //x_gen_alpha[c] = inv_logit(beta*tr/10+gamma_ind[2]+zeta[2]+kappa);
     //x_gen_pre_alpha[c] = inv_logit(beta*tr/10+gamma_ind[2]+zeta[1]+kappa);
   }
   
   real log_lik[N];
   for (n in 1:N){
   log_lik[n] = beta_binomial_lpmf(trans_unvacc[n] + trans_vacc[n] |  HH_vacc[n]+ HH_unvacc[n], theta*pbar[n]/(1-pbar[n]), theta);
   }
}
