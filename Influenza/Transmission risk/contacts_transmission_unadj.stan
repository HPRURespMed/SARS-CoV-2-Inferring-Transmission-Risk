data {
  int<lower=0> N;
  vector[N] d_meas;
  vector[N] d_std;
  //int nvoc;
  //int<lower=1, upper=nvoc> voc[N];
  int<lower=0, upper=1> trans[N];
  //int n_con_age;
  //int n_ind_age;
  //int<lower=1, upper=n_con_age> contact_age[N];
  //int<lower=1, upper=n_ind_age> index_age[N];
  //real<lower=0> crowding[N];
  //int<lower=1> HH[N];
}

// The parameters accepted by the model.
parameters {
  real alpha;
  real beta;
  //real kappa;
  // real zeta[nvoc];
  real d[N];
  //real zeta_mu[nvoc];
  //real<lower=0.3> hp_sd[nvoc];
  //real zeta_sd[nvoc];
  //real hp_mu_mu;
  //real<lower=0.3> hp_mu_sd;
  //real hp_sd_mu;
  //real<lower=0.3> hp_sd_sd;
  // real gamma;
  
  //real gamma_con[n_con_age];
  //real gamma_ind[n_ind_age];
  //real sigma;
}

// transformed parameters {
//   real zeta[nvoc];
//   for (n in 1:nvoc){
//     zeta[n] = zeta_mu[n] + hp_sd[n]*zeta_sd[n];
//   }
// }

// The model to be estimated.
model {
  alpha ~ normal(0, 2);
  beta ~ normal(0, 2); // decline rate param
  // gamma_con ~ normal(0, 2);
  // gamma_ind ~ normal(0, 2);

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
    //trans[n] ~ bernoulli_logit(alpha+beta*d[n]+gamma_con[contact_age[n]]+gamma_ind[index_age[n]]+zeta[voc[n]]+kappa*crowding[n]);
    // trans[n] ~ bernoulli_logit(HH[n], alpha + beta*d[n] + zeta[voc[n]]);
    trans[n] ~ bernoulli_logit(alpha+beta*d[n]);
  }
}

generated quantities {
  vector[81] x_gen;


  int c;
  
  c = 0;
   for (tr in -20:60){
     c = c+1;
     x_gen[c] = inv_logit(alpha + beta*tr/10);
     // Vaccinated Delta
     // x_gen_delta[c] = inv_logit(alpha + beta*tr/10+zeta[3]+gamma_con[2]+gamma_ind[2]+kappa);
     // x_gen_alpha[c] = inv_logit(alpha + beta*tr/10+zeta[2]+gamma_con[2]+gamma_ind[2]+kappa);
     // x_gen_pre_alpha[c] = inv_logit(alpha + beta*tr/10+zeta[1]+gamma_con[2]+gamma_ind[2]+kappa);
   }
}
