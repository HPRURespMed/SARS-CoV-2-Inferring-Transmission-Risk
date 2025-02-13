data {
  int<lower=0> N;
  vector[N] d_meas;
  vector[N] d_std;
  int<lower=0, upper=1> trans[N];
}

parameters {
  real alpha;
  real beta;
  real d[N];
}

model {
  alpha ~ normal(0, 2);
  beta ~ normal(0, 2); // decline rate param
  
  d_meas ~ normal(d, d_std);
  d ~ normal(0, 2);
  for (n in 1:N){
    trans[n] ~ bernoulli_logit(alpha+beta*d[n]);
  }
}

generated quantities {
  vector[51] x_gen;

  int c;
  
  c = 0;
   for (tr in -20:30){
     c = c+1;
     x_gen[c] = inv_logit(alpha + beta*tr/10);
   }
}
