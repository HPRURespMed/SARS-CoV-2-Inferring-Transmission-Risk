library(fitdistrplus)
library(readxl, plyr)
library(dplyr)
library(gridExtra)
library(grid)
library(ggplot2)
library(lattice)
library(rstan)
library(ggplot2)
library(Rlab)
library(distr)
library(pracma)
library(sfsmisc)
library(magrittr)
library(ggpubr)
library(Hmisc)
library(tidyr)
library(stringr)
library(npreg)
library(tableone)
library(boot)

options(mc.cores = 4)
rstan_options(auto_write = TRUE)

model <- stan_model("contacts_transmission_unadj.stan") # Adjusted model
#model <- stan_model("contacts_transmission_adj.stan") # Unadjusted model

# Functions
extract_data <- function(part){
  df_part <- df_merged[df_merged$obs_id==part, ]
  t <- df_part$study_day
  vl <- log(df_part$copy)
  return(data.frame("t"=t, "vl"=vl))
}

posterior_copy_i <- function(ii) {
  return(
    data.frame(
      "v_b" = posterior_copy$v_b[, ii],
      "p_ln_v_max" = exp(posterior_copy$p_ln_v_max[, ii]),
      "v_s" = posterior_copy$v_s,
      "fp" = posterior_copy$fp
    )
  )
}

df_filt <- readRDS("sars_contacts.RDS")
df_index <- readRDS("sars_index.RDS")
#df_filt <- readRDS("sars_contacts_sensitivity.RDS") # Sero-negative cohort

d_meas <- df_filt$decline_time
d_meas_i <- df_index$decline_time
sd_d <- sd(d_meas_i)
mean_d <- mean(d_meas_i)
d_std <- df_filt$std_time/sd_d
d_meas <- (d_meas - mean_d)/sd_d

N <- length(d_meas)

voc <- as.integer(factor(df_filt$variants, levels=c("Pre-alpha", "Alpha", "Delta")))
nvoc <- length(unique(voc))

trans <- as.integer(df_filt$sarspos)

contact_age <- as.integer(cut(as.integer(df_filt$age), breaks=c(0, 19, 35, 69)))
index_age <- as.integer(cut(as.integer(df_filt$index_age), breaks=c(0, 19, 35, 69)))
n_contact_age <- length(unique(contact_age)) 
n_index_age <- length(unique(index_age)) 

comorb_ind <- as.integer(factor(df_filt$comorb_preg.x, levels=c("No", "Yes")))-1
comorb_con <- as.integer(factor(df_filt$comorb_preg.y, levels=c("No", "Yes")))-1

data.stan = list(
  N = N,
  d_std = d_std,
  d_meas = d_meas,
  nvoc = nvoc,
  voc = voc,
  trans = trans,
  n_con_age = n_contact_age,
  n_ind_age = n_index_age,
  contact_age = contact_age,
  index_age = index_age,
  crowding = df_filt$people_beds,
  comorb_ind = comorb_ind,
  comorb_con = comorb_con
)

fit <- sampling(model, chains=4, cores=4, data=data.stan,
                iter=5000, warmup=2000, seed=497101)

posterior_contact <- rstan::extract(fit)

mean(posterior_contact$beta>0)

prob_1 <- inv.logit(posterior_contact$alpha)
prob_2 <- inv.logit(posterior_contact$alpha + posterior_contact$beta)

diff_prob <- prob_2 - prob_1

mean(diff_prob)*100
quantile(diff_prob, prob=0.025)*100
quantile(diff_prob, prob=0.975)*100

mean(diff_prob/prob_1)*100
quantile(diff_prob/prob_1, prob=0.025)*100
quantile(diff_prob/prob_1, prob=0.975)*100

