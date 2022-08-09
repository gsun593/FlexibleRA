library(haven)
library(dplyr)
library(ggplot2)
library(splines)
library(mgcv)
library(readr)
library(lfe)

# OHIE
#####
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

descriptives %>% inner_join(ed, by = 'person_id') %>% left_join(inp, by = 'person_id') %>%
  inner_join(state, by = 'person_id')


dat <- descriptives %>% inner_join(ed, by = 'person_id') %>% 
  left_join(inp, by = 'person_id') %>%
  inner_join(state, by = 'person_id') %>%
  mutate(Y = ed_charg_tot_ed>0, W = treatment) %>%
  (function(.) .[,c('Y', 'D', 'W',ctrl_cols, 'gender_inp', 'age_inp', 'health_last12_inp', 'edu_inp')]) %>%
  mutate(gender_inp = coalesce(gender_inp, -2), age_inp = coalesce(age_inp, -2),
         health_last12_inp = coalesce(health_last12_inp, -2), edu_inp = coalesce(edu_inp, -2))

dat %>% felm(formula(paste('Y~W+',paste(c(ctrl_cols,'gender_inp', 'age_inp', 'health_last12_inp', 'edu_inp'),collapse='+'),sep='')), data = .) %>% summary()
dat %>% felm(D ~ W, data = .) %>% summary
dat %>% felm(formula(paste('D~W+',paste(c(ctrl_cols,'gender_inp', 'age_inp', 'health_last12_inp', 'edu_inp'),collapse='+'),sep='')), data = .) %>% summary()


dat %>% write_csv('dat_for_RA_OHIE.csv')
#####

# Ferraro Price
#####
dat <- read_dta("~/Downloads/ferraroprice/090113_TotWatDat_cor_merge_Price.dta")
dat %>% mutate(Y = jun07_3x+jul07_3x+aug07_3x+sep07_3x, W = treatment) %>%
  select(Y,W,jun06,jul06,aug06,sep06,oct06,nov06,dec06,jan07,feb07,mar07,apr07_3x,may07_3x) %>%
  write_csv('dat_for_RA_ferraroprice.csv')

dat <- read_dta("~/Downloads/ferraroprice/090113_TotWatDat_cor_merge_Price.dta")
dat %>% mutate(Y = jun07_3x+jul07_3x+aug07_3x+sep07_3x, W = treatment,
               jun06=log(jun06+1),
               jul06=log(jul06+1),aug06=log(aug06+1),sep06=log(sep06+1),oct06=log(oct06+1),
               nov06=log(nov06+1),dec06=log(dec06+1),jan07=log(jan07+1),feb07=log(feb07+1),
               mar07=log(mar07+1),apr07_3x=log(apr07_3x+1),may07_3x=log(may07_3x+1)) %>%
  select(Y,W,jun06,jul06,aug06,sep06,oct06,nov06,dec06,jan07,feb07,mar07,apr07_3x,may07_3x) %>%
  write_csv('dat_for_RA_ferraroprice_logged.csv')
#####

# CogX
#####
all_dat <- read_dta("~/Downloads/CHECC_data/2185_DATA.dta")
demog <- read_dta('~/Downloads/CHECC_data/demog_covariates.dta')
mothers <- read_dta('~/Downloads/CHECC_data/mothers_edu.dta')
grades <- read_dta("~/Downloads/CHECC_data/parcc_disc_grades_wnewadditions.dta")
cog_noncog <- read_dta("~/Downloads/CHECC_data/pre_mid_post_sl_cog_ncog.dta")


outcome <- all_dat %>% filter(!is.na(attendance_cogx_dummy), ismatched==1) %>%
  mutate(Y = cog_sl) %>% select(child,Y)

demog_covariates <- demog %>% filter(!is.na(child)) %>%
  mutate(W = as.numeric(treatment!='control'), has_bw = as.numeric(!is.na(birthweight)),
         bw = coalesce(birthweight,0),
         hinc = coalesce(household_income, 0),
         has_mother_age_cb = as.numeric(!is.na(mother_age_childbirth)),
         mother_age_cb = coalesce(mother_age_childbirth,0),
         mother_education_pre = coalesce(mother_education_pre, 0),
         mother_employed_pre = coalesce(mother_employed_pre, 0),
         hl = case_when(home_language == '' ~ 0,
                        home_language == 'english' ~ 1,
                        home_language == 'spanish' ~ 2,
                        home_language == 'other' ~ 3,
                        home_language == 'english/spanish' ~ 4,
                        home_language == 'english/other' ~ 2)) %>%
  select(child, W, has_bw, bw, has_mother_age_cb, mother_age_cb, mother_education_pre, mother_employed_pre,hinc,hl)

adl_covariates <- all_dat %>% filter(!is.na(attendance_cogx_dummy), ismatched==1, age_pre < 5*12, treatment!='prek') %>%
  mutate(female = gender == 'female') %>%
  select(child, female, age_pre, cog_pre, ncog_pre, has_cog_pre, has_ncog_pre, race_b, race_h,
  )

outcome %>% inner_join(demog_covariates) %>% inner_join(adl_covariates) %>%
  group_by(W) %>% summarise(mean(is.na(Y)), sqrt(var(is.na(Y)) / n()),
                            n())


setwd('~/Documents/ECON/RegressionAdjustment')
dat <- outcome %>% filter(!is.na(Y)) %>% inner_join(demog_covariates) %>% inner_join(adl_covariates)
dat$child <- NULL
dat %>% lm(Y ~ W, data = .) %>% summary

dat%>% ggplot + geom_point(aes(x=cog_pre, y=Y)) +
  geom_smooth(aes(x=cog_pre,y=Y),method='lm')

dat %>% write_csv('dat_for_RA_CHECC.csv')
#####
