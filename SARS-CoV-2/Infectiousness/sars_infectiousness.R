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
library(packrat)

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

options(mc.cores = 4)
rstan_options(auto_write = TRUE)

#model <- stan_model("linear_model_adj.stan") # Adjusted
model <- stan_model("linear_model_unadj.stan") # Unadjusted

df_filt <- readRDS("incident_cases.RDS")

x_meas <- df_filt$power_auc
x_std <- df_filt$sd_power_auc

vacc <-
  as.integer(factor(df_filt$vax_status,
                    levels = c("unvaccinated", "two doses"))) - 1
N <- length(x_meas)

voc <- as.integer(as.factor(df_filt$WGS))
nvoc <- length(unique(voc))

d_meas <- df_filt$decline_time
d_std <- df_filt$decline_time_std

age <- as.integer(df_filt$age)
age <- age

data.stan = list(
  N = N,
  x_std = x_std,
  x_meas = x_meas,
  vacc = vacc,
  nvoc = nvoc,
  voc = voc,
  d_meas = d_meas,
  d_std = d_std,
  age=age
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
    max_treedepth = 13,
    adapt_delta = 0.994,
    stepsize = 0.01
  )
)

posterior_dr <- rstan::extract(fit)

mean(posterior_dr$beta>0)

d_e <- posterior_dr$beta*sd(d_meas)/sd(x_meas)

mean_d_e <- format(mean(d_e), digits=2)
low_ci <- format(quantile(d_e, prob=0.025)[[1]], digits=2)
high_ci <- format(quantile(d_e, prob=0.975)[[1]], digits=2)

mean_d_e
c_i <- paste0("'", mean_d_e, " (",low_ci,", ", high_ci, ")")
c_i

# Figure 3A
df_scatter <- data.frame(x=d_meas, y=x_meas, ex=d_std, WGS=df_filt$WGS, vax_status=df_filt$vax_status)

df_scatter <- df_scatter %>%
  mutate(vax_var = case_when(
    WGS == "Pre-Alpha" ~ "Pre-alpha",
    WGS == "Alpha" ~ "Alpha",
    WGS == "Delta" & vax_status == "unvaccinated" ~ "Delta unvaccinated",
    WGS == "Delta" & vax_status == "two doses" ~ "Delta vaccinated"
  )) %>%
  mutate(vax_var = factor(factor(vax_var, levels=c("Pre-alpha", "Alpha", "Delta unvaccinated", "Delta vaccinated"))))

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
  size=2.0) +
  
  geom_line(data = extract_pred_auc(posterior_dr$x_gen),
            aes(x = t*sd_d+mean_d,
                y = medx*sd_x+mean_x
            ), 
            col = "black", size=2) +
  theme_classic() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.title = element_blank(),
    legend.position = "bottom"
  ) +
  theme(axis.text = element_text(size = 15)) +
  theme(axis.title = element_text(size = 15)) +
  theme(axis.line = element_line(size = 0.75)) +
  xlim(min(df_scatter$x) - 0.1, max(df_scatter$x) + 1) +
  xlab("RNA VL decline time [days]") +
  ylab("AUC RNA VL/ml") +
  scale_x_continuous(limits = c(min(df_scatter$x) - 0.1, max(df_scatter$x) +
                                  0.1),
                     expand = c(0, 0)) +
  scale_y_continuous(limits = c(min(df_scatter$y) - 1.5, max(df_scatter$y) +
                                  1),
                     expand = c(0, 0)) +
  scale_color_brewer(palette = "Set1")

svglite("Figure3A.svg", width = 4*1.68, height = 4)
print(g_plot)
dev.off()

df_filt <- readRDS("incident_cases_pfu.RDS")

# AUC PFU
x_meas <- df_filt$mean_auc
x_std <- df_filt$sd_auc

vacc <-
  as.integer(factor(df_filt$vax_status,
                    levels = c("unvaccinated", "two doses"))) - 1
N <- length(x_meas)

voc <- as.integer(as.factor(df_filt$WGS))
nvoc <- length(unique(voc))

d_meas <- df_filt$pfu_time
d_std <- df_filt$pfu_time_std

age <- as.integer(df_filt$age)
age <- age#-median(age)

data.stan = list(
  N = N,
  x_std = x_std,
  x_meas = x_meas,
  vacc = vacc,
  nvoc = nvoc,
  voc = voc,
  d_meas = d_meas,
  d_std = d_std,
  age=age
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
    max_treedepth = 13,
    adapt_delta = 0.994,
    stepsize = 0.01
  )
)

posterior_dr_pfu <- rstan::extract(fit)

mean(posterior_dr_pfu$beta>0)

d_e <- posterior_dr_pfu$beta*sd(d_meas)/sd(x_meas)

mean_d_e <- format(mean(d_e), digits=2)
low_ci <- format(quantile(d_e, prob=0.025)[[1]], digits=2)
high_ci <- format(quantile(d_e, prob=0.975)[[1]], digits=2)

mean_d_e
c_i <- paste0("'", mean_d_e, " (",low_ci,", ", high_ci, ")")
c_i

# Figure 3B
df_scatter <- data.frame(x=d_meas, y=x_meas, ex=d_std, WGS=df_filt$WGS, vax_status=df_filt$vax_status)

df_scatter <- df_scatter %>%
  mutate(vax_var = case_when(
    WGS == "Pre-Alpha" ~ "Pre-alpha",
    WGS == "Alpha" ~ "Alpha",
    WGS == "Delta" & vax_status == "unvaccinated" ~ "Delta unvaccinated",
    WGS == "Delta" & vax_status == "two doses" ~ "Delta vaccinated"
  )) %>%
  mutate(vax_var = factor(factor(vax_var, levels=c("Pre-alpha", "Alpha", "Delta unvaccinated", "Delta vaccinated"))))

sd_d <- 1
mean_d <- 0
sd_x <- 1
mean_x <- 0

g_plot <- ggplot() +
  geom_ribbon(
    data = extract_pred_auc(posterior_dr_pfu$x_gen),
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
  size=2.0) +
  
  geom_line(data = extract_pred_auc(posterior_dr_pfu$x_gen),
            aes(x = t*sd_d+mean_d,
                y = medx*sd_x+mean_x
            ), 
            col = "black", size=2) +
  theme_classic() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.title = element_blank(),
    legend.position = "bottom"
  ) +
  theme(axis.text = element_text(size = 15)) +
  theme(axis.title = element_text(size = 15)) +
  theme(axis.line = element_line(size = 0.75)) +
  xlim(min(df_scatter$x) - 0.1, max(df_scatter$x) + 1) +
  xlab("RNA VL decline time [days]") +
  ylab("AUC RNA VL/ml") +
  scale_x_continuous(limits = c(min(df_scatter$x) - 0.1, max(df_scatter$x) +
                                  0.1),
                     expand = c(0, 0)) +
  scale_y_continuous(limits = c(min(df_scatter$y) - 1.5, max(df_scatter$y) +
                                  1),
                     expand = c(0, 0)) +
  scale_color_brewer(palette = "Set1")

svglite("Figure3B.svg", width = 4*1.68, height = 4)
print(g_plot)
dev.off()