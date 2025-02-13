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
library(loo)
library(svglite)
library(boot)
library(cowplot)

# Functions
extract_data <- function(part){
  df_part <- df_merged[df_merged$obs_id==part, ]
  t <- df_part$study_day
  vl <- log(df_part$copy)
  return(data.frame("t"=t, "vl"=vl))
}

# pfu_to_copy_index(ii)
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
      "t" = (-20:30) / 10,
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


df_index <- readRDS("sars_index.RDS")
#df_index <- readRDS("sars_index_seronegative.RDS") # Seronegative cohort

options(mc.cores = 4)
rstan_options(auto_write = TRUE)

model <- stan_model("index_transmission_unadj.stan") # Unadjusted model
#model <- stan_model("index_transmission_adj.stan") # Adjusted model

d_meas <- df_index$decline_time
d_std <- df_index$std_time/sd(d_meas)

sd_d <- sd(d_meas)
mean_d <- mean(d_meas)
d_meas <- (d_meas - mean(d_meas))/sd(d_meas)

N <- length(d_meas)
voc <- as.integer(factor(df_index$variants, levels=c("Pre-alpha", "Alpha", "Delta")))
nvoc <- length(unique(voc))

HH <- df_index$household_size
trans <- as.integer(df_index$sar*HH)

comorb <- as.integer(factor(df_index$comorb_preg, levels=c("No", "Yes")))-1

index_age <- as.integer(cut(as.integer(df_index$age), breaks=c(0, 19, 35, 69)))
n_ind_age <- length(unique(index_age)) 
crowding <- df_index$people_beds

data.stan = list(
  N = N,
  d_std = d_std,
  d_meas = d_meas,
  nvoc = nvoc,
  voc = voc,
  trans = trans,
  HH = HH,
  index_age=index_age,
  n_ind_age=n_ind_age,
  crowding=crowding,
  comorb=comorb
)

fit <- sampling(model, chains=4, cores=4, data=data.stan,
                iter=6000, warmup=3000, seed=497101,
                control = list(
                  max_treedepth = 13,
                  adapt_delta = 0.99994,
                  stepsize = 0.0001
                )
)

posterior_index <- rstan::extract(fit)

mean(posterior_index$beta>0)

prob_1 <- inv.logit(posterior_index$offset)
prob_2 <- inv.logit(posterior_index$offset + posterior_index$beta)

odds_1 <- (prob_1)/(1-prob_1)
odds_2 <- (prob_2)/(1-prob_2)

diff_prob <- (prob_2 - prob_1)/prob_1

mean(diff_prob)*100
quantile(diff_prob, prob=0.025)*100
quantile(diff_prob, prob=0.975)*100


library(Hmisc)
df_scatter <- data.frame(x=d_meas, y=trans, z=HH, ex=d_std)
df_scatter <- df_scatter[order(df_scatter$x),]

# Bin data
df_bin <- df_scatter

n_cuts <- 3

breaks <- df_scatter$x[c(1, 13*c(1,2)+1, 38)]

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

table(df_bin$interval)

# Figure 5A
df_scatter <- data.frame(x=d_meas, y=df_index$sar, ex=d_std, WGS=factor(factor(df_index$variants, levels=c("Pre-alpha", "Alpha", "Delta"))))

df_pred <- extract_pred_logit(posterior_index$x_gen)
df_pred <- df_pred[df_pred$t>=min(df_scatter$x), ]

g_plot <- ggplot() +
  geom_point(data=df_scatter, aes(x=x * sd_d + mean_d, y=y), alpha=0) +
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
  
  geom_label(aes(x = df_bin_sar$mid * sd_d + mean_d + 0.15, y = df_bin_sar$bin_sar + 0.01, label = paste0(df_bin_sar$n_inf, "/", df_bin_sar$n_tot)), fill = "white", size=5) +
  
  
  theme(axis.text = element_text(size = 15)) +
  theme(axis.title = element_text(size = 15)) +
  theme(axis.line = element_line(size = 0.75)) 

library(cowplot)
x_hist <- ggplot() +
  geom_histogram(data = df_scatter,
                 aes(x = x * sd_d + mean_d),
                 fill = "gray",
                 col = "black") +
  xlab("") +
  ylab("Frequency") +
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
svglite("Figure5A.svg", width = 4*1.68, height = 4)

plot_grid(
  aligned_x_hist
  , g_plot
  , ncol = 1
  , nrow = 2
  , rel_heights = c(0.6, 1)
  , rel_widths = c(1, 0.2)
)

dev.off()


