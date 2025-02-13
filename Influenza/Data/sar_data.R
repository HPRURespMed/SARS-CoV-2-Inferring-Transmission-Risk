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

df_data <- readRDS(file="influenza_data.RDS")

df_index <- df_data %>%
  filter(is_index==TRUE) 

df_contacts <- df_data %>%
  filter(is_index==FALSE) %>%
  group_by(household_id) %>%
  summarise(secondary_transmissions = sum(pcr_pos))

df_household <- df_data %>%
  filter(is_index==FALSE) %>%
  group_by(household_id) %>%
  summarise(hh_size = length(pcr_pos))

df_sar <- merge(df_contacts, df_household, by="household_id") 

df_sar <- df_sar %>%
  mutate(sar = secondary_transmissions/(hh_size)) %>%
  dplyr::select(secondary_transmissions, sar, household_id, hh_size)

df_index_sar <- merge(df_index, df_sar, by="household_id")

saveRDS(df_index_sar, file="influenza_sar.RDS")

df_contacts <- df_data %>%
  filter(is_index==FALSE)

index_id <- df_data %>%
  filter(is_index==TRUE) %>%
  mutate(ID_index = ID) %>%
  mutate(obs_id_index = obs_id) %>%
  dplyr::select(ID_index, obs_id_index, household_id)

df_contacts <- merge(df_contacts, index_id, by="household_id")

saveRDS(df_contacts, file="influenza_contacts.RDS")



