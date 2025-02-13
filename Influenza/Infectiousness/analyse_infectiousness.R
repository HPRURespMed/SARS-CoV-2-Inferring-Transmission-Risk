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
library(svglite)
library(matrixStats)

load("../Fit trajectories/fit_vl.Rdata")

posterior <- rstan::extract(fit)

df_merged <- readRDS(file="../Data/incident_trajectories.RDS") 
df_merged$group <- 1

posterior_copy_i <- function(ii) {
  return(
    data.frame(
      "v_a" = posterior$v_a[, ii],
      "v_b" = posterior$v_b[, ii],
      "p_ln_v_max" = exp(posterior$p_ln_v_max[, ii]),
      "v_s" = posterior$v_s,
      "fp" = posterior$fp,
      "auc" = posterior$auc[, ii]
    )
  )
}

vl_function <- function(row) {
  vl<-log(row[2])+log(1/row[3]+1/row[4])-logSumExp(c(-log(row[4])-row[1]/row[3],-log(row[3])+row[1]/row[4])); 
  return(vl)
}

calc_vl <- function(group) {
  cc <- 0
  for (part in group) {
    cc <- cc+1
    for (j in (1:50)) {
      tpi <- (j - 1) * 0.5 - 8
      p_t_max <- posterior$p_t_max[, part]
      vmax <- exp(posterior$p_ln_v_max[, part])
      va <- posterior$v_a[, part]
      vb <- posterior$v_b[, part]
      df_vl <-
        data.frame(
          "t" = tpi - p_t_max,
          "vmax" = vmax,
          "va" = va,
          "vb" = vb
        )
      vl <- apply(df_vl, 1, vl_function)
      vl <- ifelse(vl < 0 | is.infinite(vl) | is.na(vl), 0, vl)
      vl <- exp(vl)
      
      if ((cc == 1)&(j==1)) {
        df_vl_all <- data.frame("x" = vl)
        colnames(df_vl_all) <- paste0("cc", cc, "j",1)
      } else {
        df_vl_all[, ncol(df_vl_all) + 1] <- vl
        colnames(df_vl_all)[ncol(df_vl_all)] <- paste0("cc", cc, "j", j)
      }
    }
  }
  return(df_vl_all)
}

df_vl_tot <- calc_vl(unique(df_merged$obs_id))

calc_power_auc_v2 <- function(part, power_coeff) {
  s_len <- length(df_vl_tot[, 1])
  pred_len <- 50
  pred_vl_df <- data.frame()
  pred_vl_df[1:s_len, (1:pred_len)] <-
    df_vl_tot[, (1:pred_len) + (part - 1) * pred_len]^power_coeff
  return(rowSums(pred_vl_df)*0.5-pred_len*0.5)
}

df_extract <- data.frame("obs_id"=unique(df_merged$obs_id), "ID"=unique(df_merged$ID))
df_extract <- df_extract %>%
  rowwise() %>%
  mutate("decline_time"=median(posterior_copy_i(obs_id)$v_b)) %>%
  mutate("std_time"=sd(posterior_copy_i(obs_id)$v_b)) %>%
  mutate("growth_time"=median(posterior_copy_i(obs_id)$v_a)) %>%
  mutate("sd_growth_time"=sd(posterior_copy_i(obs_id)$v_a)) %>%
  mutate("auc_copy"=mean(posterior_copy_i(obs_id)$auc)) %>%
  mutate("sd_auc"=sd(posterior_copy_i(obs_id)$auc)) %>%
  mutate("copy_max"=mean(posterior_copy_i(obs_id)$p_ln_v_max)) %>%
  mutate("copy_peak_sd"=sd(posterior_copy_i(obs_id)$p_ln_v_max)) %>%
  mutate("low_CrI"=quantile(posterior_copy_i(obs_id)$v_b, prob=0.025)) %>%
  mutate("up_CrI"=quantile(posterior_copy_i(obs_id)$v_b, prob=0.975)) %>%
  mutate("CrI_width"=up_CrI-low_CrI) %>%
  mutate("std"=sd(posterior_copy_i(obs_id)$v_b)) %>%
  mutate("power_auc"=mean(calc_power_auc_v2(obs_id, 0.1))) %>%
  mutate("power_auc_sd"=sd(calc_power_auc_v2(obs_id, 0.1))) %>%
  mutate("weight"=1/CrI_width^2)


df_cov <- readRDS(file="../Data/influenza_data.RDS")
df_comb <- merge(df_extract, df_cov, by="ID", all.x=TRUE)

options(mc.cores = 4)
rstan_options(auto_write = TRUE)
model <- stan_model("linear_model_unadj.stan")

df_filt <- df_comb

x_meas <- df_filt$power_auc
x_std <- df_filt$power_auc_sd

vacc <- df_filt$antiviral
vacc[is.na(vacc)] <- 0
N <- length(x_meas)

voc <- as.vector(ones(1, nrow(df_filt)))
nvoc <- length(unique(voc))

d_meas <- df_filt$decline_time
d_std <- df_filt$std_time

age <- as.integer(df_filt$age)
age[is.na(age)] <- as.integer(median(age, na.rm = TRUE))

data.stan = list(
  N = N,
  x_std = x_std,
  x_meas = x_meas,
  d_meas = d_meas,
  d_std = d_std
)

fit <- sampling(
  model,
  chains = 4,
  cores = 4,
  data = data.stan,
  iter = 2000,
  warmup = 1000,
  seed = 497101,
  control = list(
    max_treedepth = 15,
    adapt_delta = 0.9994,
    stepsize = 0.001
  )
)

posterior_dr <- rstan::extract(fit)

sum(posterior_dr$beta>0)/length(posterior_dr$beta)

quantile(posterior_dr$beta, prob=0.5)
quantile(posterior_dr$beta, prob=0.025)
quantile(posterior_dr$beta, prob=0.975)

extract_pred_auc <- function(pred_auc) {
  
  pred_auc <- as.data.frame(pred_auc)
  
  pred_vl_0.5 <- pred_auc %>%
    summarise(across(.cols = everything(),  ~ quantile(.x, probs = 0.5))) %>% t(.)
  
  pred_vl_0.05 <- pred_auc %>%
    summarise(across(.cols = everything(),  ~ quantile(.x, probs = 0.05))) %>% t(.)
  
  pred_vl_0.95 <- pred_auc %>%
    summarise(across(.cols = everything(),  ~ quantile(.x, probs = 0.95))) %>% t(.)
  
  mu <-
    data.frame(
      "t" = (-20:50) / 10,
      "medx" = pred_vl_0.5,
      "CI_l" = pred_vl_0.05,
      "CI_u" = pred_vl_0.95
    )
  return(mu)
}

df_scatter <- data.frame(x=d_meas, y=x_meas, ex=d_std, antiviral=factor(vacc))

sd_d <- 1
mean_d <- 0
sd_x <- 1
mean_x <- 0

g_plot <- ggplot() +
  geom_ribbon(
    data = extract_pred_auc(posterior_dr$x_gen),
    aes(x = t*sd_d+mean_d,
        ymin = CI_l*sd_x+mean_x,
        ymax = CI_u*sd_x+mean_x),
    fill = "grey",
    alpha = 0.5
  ) +
  
  geom_point(data = df_scatter, aes(
    x = x,
    y = y
  ),
  col="black",
  size=2.0) +

  geom_line(data = extract_pred_auc(posterior_dr$x_gen),
            aes(x = t*sd_d+mean_d,
                y = medx*sd_x+mean_x
            ), 
            col = "black",
            size=2) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.title = element_blank(),
    legend.position = "bottom"
  ) +
  xlim(min(df_scatter$x) - 0.1, max(df_scatter$x) + 1) +
  xlab("RNA VL decline time [days]") +
  ylab("AUC RNA VL/ml") +
  scale_x_continuous(limits = c(min(df_scatter$x) - 0.1, max(df_scatter$x) +
                                  0.1),
                     expand = c(0, 0)) +
  scale_y_continuous(limits = c(min(df_scatter$y) - 1, max(df_scatter$y)+1),
                     expand = c(0, 0)) +
  theme(axis.text = element_text(size = 15)) +
  theme(axis.title = element_text(size = 15)) +
  theme(axis.line = element_line(size = 0.75)) +
  scale_color_brewer(palette = "Set1")# +

svglite("Figure3C.svg", width = 4*1.68, height = 4)
print(g_plot)
dev.off()
