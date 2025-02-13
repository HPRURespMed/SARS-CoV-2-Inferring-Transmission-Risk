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
library(svglite)
library(boot)
library(cowplot)
library(Hmisc)

load("../Fit trajectories/fit.Rdata") 
posterior <- rstan::extract(fit)
remove(fit)

df_index <- readRDS(file = "../Data/influenza_sar.RDS")
df_merged <- readRDS(file="../Data/index_trajectories.RDS")
df_merged$group <- 1

posterior_copy_i <- function(ii) {
  return(
    data.frame(
      "v_b" = posterior$v_b[, ii],
      "p_ln_v_max" = exp(posterior$p_ln_v_max[, ii]),
      "v_s" = posterior$v_s,
      "fp" = posterior$fp
    )
  )
}

extract_pred_logit <- function(pred) {
  
  pred <- as.data.frame(pred)
  
  pred_vl_0.5 <- pred %>%
    summarise(across(.cols = everything(),  ~ quantile(.x, probs = 0.5))) %>% t(.)
  
  pred_vl_0.05 <- pred %>%
    summarise(across(.cols = everything(),  ~ quantile(.x, probs = 0.05))) %>% t(.)
  
  pred_vl_0.95 <- pred %>%
    summarise(across(.cols = everything(),  ~ quantile(.x, probs = 0.95))) %>% t(.)
  
  mu <-
    data.frame(
      "t" = (-20:35) / 10,
      "medx" = pred_vl_0.5,
      "CI_l" = pred_vl_0.05,
      "CI_u" = pred_vl_0.95
    )
  return(mu)
}

extract_pred_conf <- function(pred) {
  
  pred <- as.data.frame(pred)
  
  pred_vl_0.5 <- pred %>%
    summarise(across(.cols = everything(),  ~ quantile(.x, probs = 0.5))) %>% t(.)
  
  pred_vl_0.05 <- pred %>%
    summarise(across(.cols = everything(),  ~ quantile(.x, probs = 0.05))) %>% t(.)
  
  pred_vl_0.95 <- pred %>%
    summarise(across(.cols = everything(),  ~ quantile(.x, probs = 0.95))) %>% t(.)
  
  mu <-
    data.frame(
      "t" = (0:99),
      "medx" = pred_vl_0.5,
      "CI_l" = pred_vl_0.05,
      "CI_u" = pred_vl_0.95
    )
  return(mu)
}

# Extract results
df_extract <- data.frame("obs_id"=unique(df_merged$obs_id), "ID"=unique(df_merged$ID))
df_extract <- df_extract %>%
  rowwise() %>%
  mutate("decline_time"=median(posterior_copy_i(obs_id)$v_b)) %>%
  mutate("std_time"=sd(posterior_copy_i(obs_id)$v_b)) %>%
  mutate("peak_vl"=median(log10(posterior_copy_i(obs_id)$p_ln_v_max))) %>%
  mutate("std_peak"=sd(log10(posterior_copy_i(obs_id)$p_ln_v_max)))

df_comb <- merge(df_extract, df_index, by="ID")

cowling_contacts <- readRDS(file="../Data/influenza_contacts.RDS")

df_extract_index <- df_extract %>%
  #filter(ID %in% df_traj_v2$ID) %>%
  filter(ID %in% df_comb$ID) %>%
  mutate(ID_index = ID) %>%
  #mutate(obs_id_index = obs_id) %>%
  dplyr::select(-obs_id, -ID)
df_comb_contacts <- merge(cowling_contacts, df_extract_index, by="ID_index")

df_cov <- readRDS(file="../Data/influenza_data.RDS")

df_cov_index <- df_cov %>%
  filter(is_index==TRUE) %>%
  mutate(index_vacc = vacc) %>%
  mutate(index_age = age) %>%
  mutate(index_antiviral = antiviral) %>%
  mutate(index_sex = sex) %>%
  mutate(index_subtype = subtype) %>%
  mutate(index_chronic = chronic_disease) %>%
  mutate(ID_index = ID) %>%
  dplyr::select(ID_index, household_id, index_vacc, index_antiviral, index_age, index_sex, index_subtype, index_chronic)

df_cov_contacts <- merge(df_comb_contacts, df_cov_index, by="ID_index")

df_filt <- df_cov_contacts[df_cov_contacts$index_vacc>-1&df_cov_contacts$vacc>-1&df_cov_contacts$index_subtype>0, ]

df_cov_index_cohort <- filter(df_cov_index, df_cov_index$ID_index%in%df_filt$ID_index)

options(mc.cores = 4)
rstan_options(auto_write = TRUE)
model <- stan_model("contacts_transmission_unadj.stan")

df_filt <- df_cov_contacts[df_cov_contacts$index_vacc>-1&df_cov_contacts$vacc>-1&df_cov_contacts$index_subtype>0, ]
d_meas <- df_filt$decline_time
d_meas_i <- df_extract_index$decline_time
sd_d <- sd(d_meas_i)
mean_d <- mean(d_meas_i)
d_std <- df_filt$std_time/sd_d
d_meas <- (d_meas - mean_d)/sd_d
N <- length(d_meas)
voc <- df_filt$index_subtype
nvoc <- length(unique(voc))

trans <- as.integer(df_filt$pcr_pos)

contact_age <- as.integer(cut(as.integer(df_filt$age), breaks=c(-1, 19, 35, 150)))
index_age <- as.integer(cut(as.integer(df_filt$index_age), breaks=c(-1, 19, 35, 150)))

n_contact_age <- length(unique(1:max(contact_age))) 
n_index_age <- length(unique(1:max(index_age))) 

ind_vacc <- df_filt$index_vacc
con_vacc <- df_filt$vacc

antiviral <- df_filt$index_antiviral

chronic_ind <- df_filt$index_chronic
chronic_con <- df_filt$chronic_disease

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
  antiviral = antiviral,
  ind_vacc = ind_vacc,
  con_vacc = con_vacc,
  chronic_ind = chronic_ind,
  chronic_con = chronic_con
)

fit <- sampling(model, chains=4, cores=4, data=data.stan,
                iter=5000, warmup=2000, seed=497101, control=list(max_treedepth=12, adapt_delta = 0.994, stepsize=0.01))

posterior_contact <- rstan::extract(fit)
sum(posterior_contact$beta>0)/length(posterior_contact$beta)

prob_1 <- inv.logit(posterior_contact$alpha)
prob_2 <- inv.logit(posterior_contact$alpha + posterior_contact$beta)

hist(prob_2)

diff_prob <- prob_2 - prob_1

mean(diff_prob)*100
quantile(diff_prob, prob=0.025)*100
quantile(diff_prob, prob=0.975)*100

mean(diff_prob/prob_1)*100
quantile(diff_prob/prob_1, prob=0.025)*100
quantile(diff_prob/prob_1, prob=0.975)*100

# Secondary attack rate model
df_filt <- df_cov_contacts[df_cov_contacts$index_vacc>-1&df_cov_contacts$vacc>-1&df_cov_contacts$index_subtype>0, ]

df_sar <- data.frame("ID_index" = unique(df_filt$ID_index))

trans_vacc <- c()
trans_unvacc <- c()
HH_vacc <- c()
HH_unvacc <- c()
for (id in df_sar$ID_index) {
  trans_vacc <-
    c(trans_vacc, sum(df_filt$pcr_pos[df_filt$vacc == 1 &
                                        df_filt$ID_index %in% id]))
  trans_unvacc <-
    c(trans_unvacc, sum(df_filt$pcr_pos[df_filt$vacc == 0 &
                                          df_filt$ID_index %in% id]))
  HH_vacc <-
    c(HH_vacc, length(df_filt$pcr_pos[df_filt$vacc == 1 &
                                        df_filt$ID_index %in% id]))
  HH_unvacc <-
    c(HH_unvacc, length(df_filt$pcr_pos[df_filt$vacc == 0 &
                                          df_filt$ID_index %in% id]))
}
df_sar$trans_vacc <- trans_vacc
df_sar$trans_unvacc <- trans_unvacc
df_sar$HH_vacc <- HH_vacc
df_sar$HH_unvacc <- HH_unvacc

options(mc.cores = 4)

rstan_options(auto_write = TRUE)
model <- stan_model("index_transmission_unadj.stan")

df_cov_index_cohort <- filter(df_cov_index, df_cov_index$ID_index%in%df_filt$ID_index)

df_filt_index <- merge(df_cov_index_cohort, df_sar, by="ID_index")
df_filt_index <- merge(df_filt_index, df_extract_index, by="ID_index")

d_meas <- df_filt_index$decline_time
d_std <- df_filt_index$std_time/sd(d_meas)

sd_d <- sd(d_meas)
mean_d <- mean(d_meas)

d_meas <- (d_meas - mean(d_meas))/sd(d_meas)
N <- length(d_meas)

ind_vacc <- df_filt_index$index_vacc
antiviral <- df_filt_index$index_antiviral
chronic <- df_filt_index$index_chronic

voc <- df_filt_index$index_subtype
nvoc <- length(unique(voc))

HH_vacc <- df_filt_index$HH_vacc
trans_vacc <- as.integer(df_filt_index$trans_vacc)
HH_unvacc <- df_filt_index$HH_unvacc
trans_unvacc <- as.integer(df_filt_index$trans_unvacc)

trans <- trans_vacc + trans_unvacc
HH <- HH_vacc + HH_unvacc

sar <- trans/HH

index_age <- as.integer(cut(as.integer(df_filt_index$index_age), breaks=c(0, 19, 35, 150)))
n_ind_age <- length(unique(index_age)) 

data.stan = list(
  N = N,
  d_std = d_std,
  d_meas = d_meas,
  nvoc = nvoc,
  voc = voc,
  trans_unvacc = trans_unvacc,
  HH_unvacc = HH_unvacc,
  trans_vacc = trans_vacc,
  HH_vacc = HH_vacc,
  index_age=index_age,
  n_ind_age=n_ind_age,
  ind_vacc=ind_vacc,
  antiviral=antiviral,
  chronic=chronic
)

fit <- sampling(model, chains=4, cores=4, data=data.stan,
                iter=6000, warmup=3000, seed=497101, control=list(max_treedepth=12, adapt_delta = 0.994, stepsize=0.01))

posterior_index <- rstan::extract(fit)

sum(posterior_index$beta>0)/length(posterior_index$beta)

prob_1 <- inv.logit(posterior_index$alpha)
prob_2 <- inv.logit(posterior_index$alpha + posterior_index$beta)

diff_prob <- (prob_2 - prob_1)

quantile(diff_prob/prob_1, prob=0.5)*100
mean(diff_prob/prob_1)*100
quantile(diff_prob/prob_1, prob=0.025)*100
quantile(diff_prob/prob_1, prob=0.975)*100

df_scatter <- data.frame(x=d_meas, y=trans_vacc+trans_unvacc, z=HH_unvacc+HH_vacc, ex=d_std)
df_scatter <- df_scatter[order(df_scatter$x),]

# Bin data
df_bin <- df_scatter

n_cuts <- 4
breaks <- df_scatter$x[c(1, 38*c(1,2,3)+1, 153)]

df_bin$interval <- findInterval(df_bin$x, breaks[1:n_cuts])
df_bin$mid <- breaks[df_bin$interval] + diff(breaks)[df_bin$interval]/2

df_bin_sar <- data.frame(interval = sort(unique(df_bin$interval)), mid = sort(unique(df_bin$mid)))
df_bin_sar$bin_sar <- df_bin_sar$mid
df_bin_sar$CI_lower <- df_bin_sar$mid
df_bin_sar$CI_upper <- df_bin_sar$mid
df_bin_sar$n_inf <- df_bin_sar$mid
df_bin_sar$n_tot <- df_bin_sar$mid
for (i in df_bin_sar$interval) {
  df_bin_sar$bin_sar[df_bin_sar$interval == i] <-
    sum(df_bin$y[df_bin$interval == i])/sum(df_bin$z[df_bin$interval == i])
  df_bin_sar$CI_lower[df_bin_sar$interval == i] <-
    binconf(x = sum(df_bin$y[df_bin$interval == i]), n = sum(df_bin$z[df_bin$interval ==
                                                                        i]))[2]
  df_bin_sar$CI_upper[df_bin_sar$interval == i] <-
    binconf(x = sum(df_bin$y[df_bin$interval == i]), n = sum(df_bin$z[df_bin$interval ==
                                                                        i]))[3]
  df_bin_sar$n_inf[df_bin_sar$interval == i] <- sum(df_bin$y[df_bin$interval == i])
  df_bin_sar$n_tot[df_bin_sar$interval == i] <- sum(df_bin$z[df_bin$interval == i])
  
}

df_scatter <-
  data.frame(
    x = d_meas,
    y = (trans_unvacc + trans_vacc) / (HH_unvacc + HH_vacc),
    ex = d_std,
    WGS = voc
  )

df_pred <- extract_pred_logit(posterior_index$x_gen)
df_pred <- df_pred[df_pred$t>min(df_scatter$x), ]

g_plot <- ggplot() +
  geom_ribbon(
    data = df_pred,
    aes(
      x = t * sd_d + mean_d,
      ymin = CI_l,
      ymax = CI_u
    ),
    fill = "gray",
    alpha = 0.75
  ) +
  geom_line(data = df_pred,
            aes(x = t * sd_d + mean_d,
                y = medx)) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.title = element_blank(),
    legend.position = "bottom"
  ) +
  xlim(0.5, 3.2) +
  xlab("Index RNA VL decline time [days]") +
  ylab("Secondary attack rate") +
  geom_vline(xintercept = mean_d, linetype = "dashed", size=0.75) +
  
  geom_point(
    data = df_bin_sar,
    aes(x = mid * sd_d + mean_d,
        y = bin_sar),
    col = "black",
    size = 3
  ) +
  
  geom_errorbar(
    data = df_bin_sar,
    aes(
      x = mid * sd_d + mean_d,
      ymin = CI_lower,
      ymax = CI_upper
    ),
    width = 0.2
  ) +
  geom_vline(
    xintercept = breaks[1] * sd_d + mean_d,
    linetype = "solid",
    size = 0.75
  ) +
  geom_vline(
    xintercept = breaks[2] * sd_d + mean_d,
    linetype = "solid",
    size = 0.75
  ) +
  geom_vline(
    xintercept = breaks[3] * sd_d + mean_d,
    linetype = "solid",
    size = 0.75
  ) +
  geom_vline(
    xintercept = breaks[4] * sd_d + mean_d,
    linetype = "solid",
    size = 0.75
  ) +
  geom_vline(
    xintercept = breaks[5] * sd_d + mean_d,
    linetype = "solid",
    size = 0.75
  ) +
  geom_label(aes(x = df_bin_sar$mid * sd_d + mean_d + 0.18, y = df_bin_sar$bin_sar + 0.01, label = paste0(df_bin_sar$n_inf, "/", df_bin_sar$n_tot)), fill = "white", size=5) +
  
  theme(axis.text = element_text(size = 15)) +
  theme(axis.title = element_text(size = 15)) +
  theme(axis.line = element_line(size = 0.75))

x_hist <- ggplot() +
  geom_histogram(data = df_scatter,
                 aes(x = x * sd_d + mean_d),
                 fill = "gray",
                 col = "black") +
  xlab("") +
  ylab("Frequency") +
  ylim(c(0, 16)) +
  xlim(0.5, 3.2) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  theme(axis.text = element_text(size = 15)) +
  theme(axis.title = element_text(size = 15)) +
  theme(axis.line.y = element_line(size = 0.75)) +
  theme(axis.line.x = element_line(size = 0.0)) +
  
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank()
  )

aligned_x_hist <- align_plots(x_hist, g_plot, align = "v")[[1]]

svglite("Figure5B.svg", width = 4*1.68, height = 4)

plot_grid(
  aligned_x_hist
  , g_plot
  , ncol = 1
  , nrow = 2
  , rel_heights = c(0.6, 1)
  , rel_widths = c(1, 0.2)
)

dev.off()
