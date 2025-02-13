library(readxl, plyr)
library(dplyr)
library(rstan)

df_merged <- readRDS(file="incident_trajectories.RDS") 
df_merged <- df_merged %>%
  mutate(group = case_when(
    (vaccinated == FALSE) & (WGS == "Pre-Alpha") ~ 1,
    (vaccinated == FALSE) & (WGS == "Alpha") ~ 1,
    (vaccinated == FALSE) & (WGS == "Delta") ~ 1,
    (vaccinated == TRUE) & (WGS == "Delta") ~ 2
  ))  %>%
  filter(!is.na(group)) # Remove participants which don't fit into any group

df_merged$obs_id <- as.integer(factor(df_merged$id_sub, levels=unique(df_merged$id_sub)))

# PFU
options(mc.cores = 8)
rstan_options(auto_write = TRUE)
model <- stan_model("vl_model.stan")

exclude_pfu <- c(12,18,23,25,41,56)

# Filter to select only those with pfu measurements
df_pfu <- df_merged[df_merged$pfu_exist, ]
df_pfu$obs_id_old <-
  as.integer(factor(df_pfu$id_sub, levels = unique(df_pfu$id_sub)))
df_pfu <- df_pfu[!df_pfu$obs_id_old%in%exclude_pfu, ]
df_pfu$obs_id <- as.integer(factor(df_pfu$id_sub, levels = unique(df_pfu$id_sub)))

x <- df_pfu$pfu
x[x==10] <- 10.01
x[x==1] <- 10
x <- log(x) - log(10)
t <- df_pfu$day[!is.na(x)]
obs_id <- df_pfu$obs_id[!is.na(x)]
x <- x[!is.na(x)]
N <- length(unique(obs_id))
group_index <- df_pfu[!is.na(x),] %>%
  group_by(obs_id) %>%
  summarise(group_index= unique(group))
group_index <- group_index$group_index
NG <- length(unique(group_index))
M <- length(x)

# Pass data to model
data.stan=list(M=M,
               N=N,
               NG=NG,
               group_index=group_index, 
               t=t,
               x=x,
               obs_id=obs_id,
               pr_fp=3.0,         # mean of prior for error proportion
               pr_fp_sd=2.0,      # SD of prior of error proportion
               pr_fp_v=2.0,       # prior for mean of error CT distribution
               pr_fp_v_s=2.2,  
               hp_v = c(12.0, 0.8, 2),
               hp_v_sd = c(12.0, 0.48, 1.8),
               hp_v_sd_sd = c(10.0, 0.64, 1.8),
               pr_v_s=1.2,
               pr_v_s_sd=0.2,
               t_sd=5,
               t_mu=5,
               lod=2.3,
               eta=1
)

initf <- function(chain_id) {
  list(
    v_s = 3 + (runif(1) - 0.5),
    vgrow = array(c(10, 0.8, 2), dim = c(3, NG)),
    v_sd = c(2, 0.5, 0.5), fp_ct_sd = 2.0 + 1.0*(runif(1)-0.5),
    fp_ct_mean=0.2+(runif(1)-0.5), l_fp = 3+(runif(1)-0.5)
  )
}

fit_pfu = sampling(model, chains=8, cores=8, data=data.stan,
           iter=8000, init=initf, init_r=0.1, warmup=3000, seed=497101, thin=2,
           control=list(max_treedepth=14, adapt_delta = 0.9994, stepsize=0.01))

save(fit_pfu, file = "fit_pfu.Rdata")
