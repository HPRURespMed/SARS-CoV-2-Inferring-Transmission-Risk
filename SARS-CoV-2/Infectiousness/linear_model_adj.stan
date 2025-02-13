data {
  int<lower=0> N;
  vector[N] x_std;
  vector[N] x_meas;
  int<lower=0, upper=1> vacc[N];
  int nvoc;
  int<lower=1, upper=nvoc> voc[N];
  vector[N] d_meas;
  vector[N] d_std;
  int<upper=120> age[N];
}

parameters {
  vector[N] x;
  real beta;
  real<lower=0.5> sigma_x;
  real gamma;
  vector[N] d;
  real zeta_mu[nvoc];
  real<lower=0.5> hp_sd[nvoc];
  real zeta_sd[nvoc];
  real hp_mu_mu;
  real<lower=0.5> hp_mu_sd;
  real hp_sd_mu;
  real<lower=0.5> hp_sd_sd;
}

transformed parameters {
  real zeta[nvoc];
  for (n in 1:nvoc){
    zeta[n] = zeta_mu[n] + hp_sd[n]*zeta_sd[n];
  }
}

model {
  beta ~ normal(0, sd(x_meas)/sd(d_meas));

  zeta_sd ~ std_normal();
  zeta_mu ~ normal(hp_mu_mu, hp_mu_sd);
  hp_mu_mu ~ normal(mean(x_meas)*0.5, 0.5*sd(x_meas));
  hp_mu_sd ~ normal(0, 4*sd(x_meas));
  hp_sd ~ normal(hp_sd_mu, hp_sd_sd);
  hp_sd_mu ~ normal(0, 4*sd(x_meas));
  hp_sd_sd ~ normal(0, 4*sd(x_meas));

  sigma_x ~ normal(0, 4*sd(x_meas));
  gamma ~ normal(0, 2*sd(x_meas)); 
  x_meas ~ normal(x, x_std);  
  d_meas ~ normal(d, d_std);
  d ~ normal(mean(d_meas), 4*sd(d_meas));
  x ~ normal(mean(x_meas), 4*sd(x_meas));
  
  for (n in 1:N){
    x[n] ~ normal(beta*d[n] + zeta[voc[n]]+gamma*vacc[n], sigma_x);
  }
}

generated quantities {
  vector[71] x_gen_pre_alpha;
  vector[71] x_gen_alpha;
  vector[71] x_gen_delta_unvacc;
  vector[71] x_gen_delta_vacc;

  int c;
  c = 0;
   for (tr in -20:50){
     c = c+1;
     x_gen_delta_vacc[c] = beta*tr/10+zeta[2]+gamma;
     x_gen_delta_unvacc[c] = beta*tr/10+zeta[2];
     x_gen_alpha[c] = beta*tr/10+zeta[1];
     x_gen_pre_alpha[c] = beta*tr/10+zeta[3];
   }
}

