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
  real alpha;
  real beta;
  real<lower=0.5> sigma_x;
  vector[N] d;
}

model {
  beta ~ normal(0, sd(x_meas)/sd(d_meas));
  alpha ~ normal(0.5*mean(x_meas), 0.5*sd(x_meas));
  sigma_x ~ normal(0, 4*sd(x_meas));
  x_meas ~ normal(x, x_std);
  d_meas ~ normal(d, d_std);
  d ~ normal(mean(d_meas), 4*sd(d_meas));
  x ~ normal(mean(x_meas), 4*sd(x_meas));
  
  for (n in 1:N){
    x[n] ~ normal(beta*d[n] + alpha, sigma_x);
  }
}

generated quantities {
  vector[71] x_gen;

  int c;
  c = 0;
   for (tr in -20:50){
     c = c+1;
     x_gen[c] = beta*tr/10+alpha;
   }
}

