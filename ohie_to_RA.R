library(haven)
library(dplyr)
library(ggplot2)
library(splines)
library(mgcv)
library(readr)
library(lfe)

descriptives <- read_dta("~/Downloads/OHIE_Public_Use_Files/OHIE_Data/oregonhie_descriptive_vars.dta")
hh_lvl <- descriptives %>% group_by(household_id) %>%
  summarise(n_obs = n())
descriptives <- descriptives %>% inner_join(hh_lvl) %>%
  filter(n_obs == 1)
descriptives$n_obs <- NULL

ed <- read_dta("~/Downloads/OHIE_Public_Use_Files/OHIE_Data/oregonhie_ed_vars.dta")
ctrl_cols <- c(ed %>% colnames %>% (function(.) .[grepl('pre',.)]) %>%
  (function(.) .[grepl('num',.)]), 'charg_tot_pre_ed')


inp <- read_dta("~/Downloads/OHIE_Public_Use_Files/OHIE_Data/oregonhie_inperson_vars.dta") %>%
  mutate(gender_inp = coalesce(gender_inp,-1), age_inp = coalesce(age_inp, -1),
         health_last12_inp = coalesce(health_last12_inp, -1),
         edu_inp = coalesce(edu_inp, -1)) %>%
  select(person_id, gender_inp, age_inp, health_last12_inp, edu_inp)

state <- read_dta("~/Downloads/OHIE_Public_Use_Files/OHIE_Data/oregonhie_stateprograms_vars.dta") %>%
  mutate(D = ohp_all_ever_inperson) %>%
  select(person_id, D)

descriptives %>% inner_join(ed, by = 'person_id') %>% inner_join(inp, by = 'person_id') %>%
  inner_join(state, by = 'person_id')


dat <- descriptives %>% inner_join(ed, by = 'person_id') %>% 
  inner_join(inp, by = 'person_id') %>%
  inner_join(state, by = 'person_id') %>%
  mutate(Y = ed_charg_tot_ed>0, W = treatment) %>%
  (function(.) .[,c('Y', 'D', 'W',ctrl_cols, 'gender_inp', 'age_inp', 'health_last12_inp', 'edu_inp')])

dat %>% felm(formula(paste('Y~W+',paste(c(ctrl_cols,'gender_inp', 'age_inp', 'health_last12_inp', 'edu_inp'),collapse='+'),sep='')), data = .) %>% summary()
dat %>% felm(D ~ W, data = .) %>% summary
dat %>% felm(formula(paste('D~W+',paste(c(ctrl_cols,'gender_inp', 'age_inp', 'health_last12_inp', 'edu_inp'),collapse='+'),sep='')), data = .) %>% summary()


setwd('~/Documents/ECON/RegressionAdjustment')
dat %>% write_csv('dat_for_RA.csv')
