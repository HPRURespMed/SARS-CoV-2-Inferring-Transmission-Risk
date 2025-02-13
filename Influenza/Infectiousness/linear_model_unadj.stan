data {
  int<lower=0> N;
  vector[N] x_std;
  vector[N] x_meas;
  vector[N] d_meas;
  vector[N] d_std;
}

parameters {
  vector[N] x;
  vector[N] d;
  real<lower=0> alpha;
  real beta;
  real<lower=5> sigma_x;
}

model {
  beta ~ normal(0, sd(x_meas)/sd(d_meas));
  alpha ~ normal(mean(x_meas)*0.5, sd(x_meas)*0.5);
  sigma_x ~ normal(0, sd(x_meas));
  x_meas ~ normal(x, x_std);    // measurement model
  d_meas ~ normal(d, d_std);
  d ~ normal(mean(d_meas), sd(d_meas));
  x ~ normal(mean(x_meas), sd(x_meas));
  
  for (n in 1:N){
    x[n] ~ normal(alpha+beta*d[n], sigma_x);
  }
}

generated quantities {
  vector[71] x_gen;
  int c;
  
  c = 0;
   for (tr in -20:50){
     c = c+1;
     x_gen[c] = alpha+beta*tr/10;
   }
}

