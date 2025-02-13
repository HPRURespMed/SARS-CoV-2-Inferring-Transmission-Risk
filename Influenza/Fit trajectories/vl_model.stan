functions {
  real logVL(real t, row_vector vp) {
    real vl;
    vl=log(vp[1])+log(1/vp[2]+1/vp[3])-log_sum_exp(-log(vp[3])-t/vp[2],-log(vp[2])+t/vp[3]); 
    return vl;
  }
}

// The input data is a vector 'y' of length 'N'.
data {
  int<lower=0> N; // n.o. participants
  int<lower=0> M; // total n.o. observations
  vector[M] t;
  vector[M] x;
  int<lower=0> obs_id[M]; // participant identifier
  int<lower=1> NG; // number of groups to be fitted (1 for only 1 group)
  int<lower=1> group_index[N];
  real<lower=0> lod; // limit of VL detection
  
  // priors that don't vary across subjects
  real<lower=0> pr_v_s;        // prior mean of SD of PCR measurement of CT
  real<lower=0> pr_v_s_sd;     // SD of prior for SD of PCR measurement of CT
  real<lower=0> pr_fp;         // prior mean of prob of false positive
  real<lower=0> pr_fp_sd;      // SD of prior for mean of prob of false positive
  real<lower=0> pr_fp_v_s;     // prior for SD of "error" CT distribution 
  real<lower=0> pr_fp_v;       // width of prior for "error" CT distribution mean and SD

  // hyperpriors
  vector[3] hp_v; // hyperprior on kinetic parameter group means
  vector<lower=0>[3] hp_v_sd; // hyperprior on kinetic parameter group SDs
  vector<lower=0>[3] hp_v_sd_sd; // hyperprior on within-group SD
  real t_mu; //mean of t_max
  real<lower=0> eta; // prior on correl matrix
  real<lower=0> t_sd;    // standard deviation of t_max
}

transformed data {
  vector[3] zeros = rep_vector(0, 3);
}

// The parameters accepted by the model
parameters {
  // hyperparameters
  matrix[3, NG] vgrow;        // viral growth rate
  vector<lower=0>[3] v_sd;     // coeff of var
  cholesky_factor_corr[3] Lc;
  
  // parameters  
  real<lower=1> v_s;          // sigma for VL
  vector[N] p_t_max;          // infection time for each subject
  matrix[N, 3] n_v;            

  real l_fp;           // false negative probability
  real fp_ct_mean;
  real<lower=1> fp_ct_sd;
}

transformed parameters {
  real ll_total;              // log likelihood for all subjects
  vector[M] log_lik;          // log likelihood for each observation
  matrix<lower=0>[N, 3] p_v;
  real<lower=0> fp;
  real mu;
  {
    int id;
    real l_fpm;
    // Non-centered parametrisation
    for(n in 1:N) {
      // p_v[n,] = exp(v_sd' .* n_v[n,]);
      p_v[n,] = exp(vgrow[,group_index[n]]' + v_sd' .* n_v[n,]);
      // p_v[n,] = exp(v[,group_index[n]]' + v_sd' .* n_v[n,]);
    }
    fp=exp(-l_fp);
    l_fpm=log(1.0-fp);
  
		  
  // likelihood
  for (m in 1:M){
    id = obs_id[m];
    mu = logVL(t[m]-p_t_max[id], p_v[id,]);
    if(x[m] <= 0.0) {
      log_lik[m] = log_sum_exp(l_fpm+normal_lcdf(0.0 | mu, v_s), -l_fp+normal_lcdf(0.0 | fp_ct_mean, fp_ct_sd));
    } else {
      log_lik[m] = log_sum_exp(l_fpm+normal_lpdf(x[m] | mu, v_s), -l_fp+normal_lpdf(x[m] | fp_ct_mean, fp_ct_sd));
    }
  }
  ll_total=sum(log_lik);
  }
}

  model {
  
  // likelihood of observed ct values
  target += ll_total;
  
  //hyperpriors 
  for(i in 1:NG){
    vgrow[, i] ~ normal(hp_v, hp_v_sd);
  }
  
  v_sd ~ normal(0, hp_v_sd_sd);
  Lc ~ lkj_corr_cholesky(eta);
  for(n in 1:N)
    n_v[n,]~multi_normal_cholesky(zeros, Lc);

  // priors
  p_t_max~normal(t_mu, t_sd); // time of peak VL (max observed CT in data at day=0)

  // params which don't vary across subjects
  v_s~normal(pr_v_s,pr_v_s_sd); // sd of CT measurements 
  l_fp~normal(pr_fp,pr_fp_sd); //  prior for "error" probability
  fp_ct_mean~normal(0,pr_fp_v);
  fp_ct_sd~normal(pr_fp_v_s,pr_fp_v);
  
  target+=ll_total;
}

generated quantities {

  vector[N] p_ln_v_max; // peak ln(VL) for each subject
  vector[N] v_a; // growth rate
  vector[N] v_b; // decay rate
  vector<lower=0>[42*N] pred_vl; // VL trajectory for each subject
  int i;
  real tpi;
  real vl;
  vector<lower=0>[N] auc;
  real<lower=0> tot_vl;

  i = 0;
  for (n in 1:N){
    v_a[n]=p_v[n, 2];
    v_b[n]=p_v[n, 3];
    p_ln_v_max[n]=log(p_v[n, 1]);
    tot_vl = 0.0;

    for(j in 1:168) {
      tpi=(j-1)*0.25-14;
      vl = logVL(tpi-p_t_max[n], p_v[n, ]);
      if(vl < 0 || is_inf(vl)) vl=0;  // report VL below detection limit as 0
      tot_vl += vl*0.25;
      if (floor(tpi)==tpi){
        i = i+1;
        pred_vl[i] = vl;
      }
    }
    auc[n] = tot_vl;
  }
}
