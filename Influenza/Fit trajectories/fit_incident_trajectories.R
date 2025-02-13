library(readxl, plyr)
library(dplyr)
library(rstan)

df_merged <- readRDS(file="../Data/incident_trajectories.RDS") 

df_merged$group <- 1

# Copy number
options(mc.cores = 4)
rstan_options(auto_write = TRUE)
model <- stan_model("vl_model.stan")

x <- df_merged$vl
x <- x - 2.954243

t <- df_merged$t
obs_id <- df_merged$obs_id[!is.na(x)]
x <- x[!is.na(x)]
x[x<0] <- 0 
N <- length(unique(obs_id))
group_index <- df_merged[!is.na(x),] %>%
  group_by(obs_id) %>%
  summarise(group_index= unique(group))
group_index <- group_index$group_index
NG <- length(unique(group_index))
M <- length(x)

# Pass data to model
data.stan<-list(M=M,
                N=N,
                NG=NG,
                group_index=group_index, 
                t=t,
                x=x,
                obs_id=obs_id,
                pr_fp=3,         # mean of prior for error proportion
                pr_fp_sd=3,      # SD of prior of error proportion
                pr_fp_v=3,       # prior for mean of error CT distribution
                pr_fp_v_s=3,
                hp_v = c(5.75, 0.8, 2),
                hp_v_sd = c(5, 0.48, 1),
                hp_v_sd_sd = c(2.5, 0.64, 1.0),
                pr_v_s=3,
                pr_v_s_sd=3,
                t_sd=5,
                t_mu=0,
                lod=0,
                eta=1
)

initf <- function(chain_id) {
  list(
    v_s = 3 + (runif(1) - 0.5),
    vgrow = array(c(5, 1, 2), dim = c(3, NG)),
    v_sd = c(5, 0.5, 0.5), fp_ct_sd=1.51+(runif(1)-0.5)
  )
}

fit = sampling(model, chains=1, cores=1, data=data.stan,
               iter=1000, init=initf, init_r=0.1, warmup=300, seed=497101, thin=2,
               control=list(max_treedepth=12, adapt_delta = 0.994, stepsize=0.01))

save(fit, file = "fit_vl.Rdata")
