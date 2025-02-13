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

# The data can be downloaded from: https://datadryad.org/stash/dataset/doi:10.5061/dryad.1p3kn
data <- read.csv(file="data.csv", col.names=seq(1:40), header = 0)

colnames(data)[1:4] <- c("household_id", "member_id", "household_size", "pcr_pos")
colnames(data)[c(8,9,10,14)] <- c("age", "sex", "vacc", "antiviral")
colnames(data)[16:28] <- paste0("D", seq(from=0, to=12))
colnames(data)[29] <- "chronic_disease"
colnames(data)[32] <- "subtype"
colnames(data)[33] <- "co-index"

data <- data %>%
  mutate(ID = paste0(paste0("H", data$household_id), paste0("M", data$member_id))) %>%
  mutate(is_index = member_id == 0)
data$obs_id <- seq(1:nrow(data))
saveRDS(data, file="influenza_data.RDS")

df_vl <-
  data.frame("obs_id" = character(),
             "t" = integer(),
             "vl" = double())
for (id in data$obs_id) {
  line_data <- data[data$obs_id == id, ]
  vl <- line_data[16:28]
  colnames(vl) <- NULL
  vl <- t(vl)[, 1]
  
  time <- seq(from = 0, to = 12)
  
  loc_max <- which(vl == max(vl))[1]
  after_peak <- loc_max:length(vl)
  vl <- vl[after_peak]
  time <- time[after_peak]
  
  if (sum(vl > 2.954243) > 1) {
    time <- time[!vl == -1]
    time <- time - min(time)
    vl <- vl[!vl == -1]
    df_vl <-
      rbind(
        df_vl,
        data.frame(
          "obs_id" = id,
          "ID" = line_data$ID,
          "t" = time,
          "vl" = vl,
          "index" = line_data$is_index
        )
      )
  }
}

df_vl <- df_vl[!df_vl$vl == -1 & df_vl$index,]
df_vl$obs_id <- as.integer(factor(df_vl$ID, levels=unique(df_vl$ID)))
saveRDS(df_vl, file="index_trajectories.RDS")

# Bootstrap: select two samples per index case
df_vl <-
  data.frame("obs_id" = character(),
             "t" = integer(),
             "vl" = double())
for (id in data$obs_id) {
  line_data <- data[data$obs_id == id,]
  vl <- line_data[16:28]
  colnames(vl) <- NULL
  vl <- t(vl)[, 1]
  select1 <- vl > -1
  time <- seq(from = 0, to = 12)
  vl <- vl[select1]
  time <- time[select1]
  
  if (sum(vl > 2.954243) > 1) {
    select2 <- vl > 2.954243
    vl <- vl[select2]
    time <- time[select2]
    loc_max <- which(vl == max(vl))[1]
    after_peak <- loc_max:length(vl)
    vl <- vl[after_peak]
    time <- time[after_peak]
    
    select3 <- sort(sample(1:length(vl), min(c(length(
      vl
    ), 2))))
    
    vl <- vl[select3]
    time <- time[select3]
    
    loc_max <- which(vl == max(vl))[1]
    after_peak2 <- loc_max:length(vl)
    vl <- vl[after_peak2]
    time <- time[after_peak2]
    
    if (sum(vl > 2.954243) > 1) {
      time <- time[!vl == -1]
      time <- time - min(time)
      vl <- vl[!vl == -1]
      
      df_vl <-
        rbind(
          df_vl,
          data.frame(
            "obs_id" = id,
            "ID" = line_data$ID,
            "t" = time,
            "vl" = vl,
            "index" = line_data$is_index,
            "n_pos" = sum(vl > 2.954243)
          )
        )
    }
  }
}

df_vl <- df_vl[!df_vl$vl == -1 & df_vl$index,]
df_vl$obs_id <- as.integer(factor(df_vl$ID, levels=unique(df_vl$ID)))
df_vl <- df_vl[!df_vl$vl==-1 & df_vl$n_pos>1 & df_vl$index & df_vl$vl>2.954243, ]
df_vl$obs_id <- as.integer(factor(df_vl$ID, levels=unique(df_vl$ID)))

df_boot <- data.frame()
for (id in unique(df_vl$obs_id)){
  df_row <- df_vl[df_vl$obs_id==id, ]
  df_row <- df_row[sort(sample(1:nrow(df_row), 2, replace=FALSE)), ]
  df_boot <- rbind(df_boot, df_row)
}
saveRDS(df_vl, file="bootstrap_index_trajectories.RDS")

# Incident trajectories
df_vl <-
  data.frame("obs_id" = character(),
             "t" = integer(),
             "vl" = double())
c_id <- 0
for (id in data$obs_id) {
  line_data <- data[data$obs_id == id, ]
  vl <- line_data[16:28]
  colnames(vl) <- NULL
  vl <- t(vl)[, 1]
  
  time <- seq(from = 0, to = 12)
  
  loc_max <- which(vl == max(vl))[1]
  loc_first <- which(vl > 2.954243)[1]
  loc_last <- which(vl > 2.954243)
  loc_last <- loc_last[length(loc_last)]
  
  if (sum(vl > 2.954243) > 1) {
    if (!loc_max %in% c(loc_first, loc_last)) {
      c_id <- c_id + 1
      
      time <- time[!vl == -1]
      time <- time - min(time)
      vl <- vl[!vl == -1]
      
      df_vl <-
        rbind(
          df_vl,
          data.frame(
            "obs_id" = c_id,
            "ID" = line_data$ID,
            "t" = time,
            "vl" = vl,
            "index" = line_data$is_index
          )
        )
    }
  }
}

df_vl <- df_vl[!df_vl$vl == -1,]
df_vl$vl <- df_vl$vl - 2.954243

saveRDS(df_vl, file="incident_trajectories.RDS")
