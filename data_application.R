library(fst)
library(parallel)
library(foreign)
library(forcats)
library(CausalGPS)
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggh4x)
library(truncnorm)
library(caret)
library(rpart)
library(rpart.plot)
library(scales)
library(data.table)
source('/nfs/home/N/nik864/shared_space/ci3_analysis/nkc_subgroupId/HERF/Code/CCIT/evaluate.sequence.R')
source('/nfs/home/N/nik864/shared_space/ci3_analysis/nkc_subgroupId/HERF/Code/CCIT/run_simu.R')
source('/nfs/home/N/nik864/shared_space/ci3_analysis/nkc_subgroupId/HERF/Code/CCIT/create.sequence.R')
source('/nfs/home/N/nik864/shared_space/ci3_analysis/nkc_subgroupId/HERF/Code/CCIT/rpart_funcs.R')
source('/nfs/home/N/nik864/shared_space/ci3_analysis/nkc_subgroupId/HERF/Code/CCIT/stratified_GPS_matching.R')
source('/nfs/home/N/nik864/shared_space/ci3_analysis/nkc_subgroupId/HERF/Code/CCIT/split_dataset.R')
source('/nfs/home/N/nik864/shared_space/ci3_analysis/nkc_subgroupId/HERF/Code/CCIT/generate_synthetic_data_covs.R')
source('/nfs/home/N/nik864/shared_space/ci3_analysis/nkc_subgroupId/HERF/Code/CCIT/generate_synthetic_data_outcome.R')

dir_out = '/nfs/home/N/nik864/shared_space/ci3_analysis/nkc_subgroupId/HERF/'

# QD Outcomes, Exposures, Covariates, and Offsets
load("~/shared_space/ci3_analysis/josey_erc_strata/Data/national_merged2016_qd.RData")
national_merged2016_qd$time_count<-national_merged2016_qd$followup_year_plus_one-national_merged2016_qd$followup_year

# convert to categorical variables
national_merged2016_qd <- national_merged2016_qd %>%
  mutate_at(vars(race, sex, dual, entry_age_break), as.factor)
# dichotomize age variable -  65 to 74, 75 and above 
national_merged2016_qd <- national_merged2016_qd %>%
  mutate(age_bin = fct_collapse(entry_age_break, low = c('1', '2'),
                                high = c('3', '4', '5', '6', '7', '8')))
# dichotomize race variable
national_merged2016_qd <- national_merged2016_qd %>%
  mutate(race_bin = fct_collapse(as.factor(race), white = c('1'),
                                 non_white = c('0', '2', '3', '4', '5', '6')))

save(national_merged2016_qd, file = paste0(dir_out, 'national_merged2016_qd.RData'))
load(paste0(dir_out, 'national_merged2016_qd.RData'))

# four potential effect modifiers (all binary)
effect_modifiers <- c('sex', 
                      'dual', 'race_bin', 'age_bin')

confounders_names <- c("year", "mean_bmi","smoke_rate","hispanic","pct_blk","medhouseholdincome","medianhousevalue"  ,
                       "poverty","education","popdensity","pct_owner_occ","summer_tmmx","winter_tmmx","summer_rmax" ,      
                       "winter_rmax","region")

dead_personyear <- national_merged2016_qd %>%
  group_by(zip, year, sex, dual, race_bin, age_bin) %>%
  summarise_at(vars(dead, time_count), sum)

save(dead_personyear, file = paste0(dir_out, 'dead_personyear2.RDats'))
load(file = paste0(dir_out, 'dead_personyear.RDats'))

confounders_names <- names(national_merged2016_qd)[13:27]
confounders <- national_merged2016_qd %>%
  group_by(zip, year, sex, dual, race_bin, age_bin) %>%
  summarise_at(vars(confounders_names), min)
save(confounders, file = paste0(dir_out, 'confounders2.RData'))

rm(national_merged2016_qd)
aggregate_data_qd <- merge(dead_personyear, confounders, 
                           by = c('zip', 'year', 'sex', 'dual', 'race_bin', 'age_bin'))


# discretize treatment values for gps matching
a.vals <- seq(min(aggregate_data_qd$pm25_ensemble), max(aggregate_data_qd$pm25_ensemble), length.out = 20)
delta_n <- a.vals[2] - a.vals[1]

aggregate_data_qd <- 
  aggregate_data_qd %>%
  mutate(treat_level = cut(pm25_ensemble, breaks = a.vals)) %>%
  mutate_at(vars(effect_modifiers), as.factor)

aggregate_data_qd <- aggregate_data_qd %>% mutate(mortality = dead/time_count) %>%
  mutate_at(vars(region, year), as.factor)

# split data into subsamples, stratifying on exposure level bins
aggregate_data_qd <- 
  split_dataset(data = aggregate_data_qd, exposure_bins = aggregate_data_qd$treat_level)

save(aggregate_data_qd, file = paste0(dir_out, 'aggregate_data_qd.RData'))
load(paste0(dir_out, 'aggregate_data_qd.RData'))

exploration.aggregate_data_qd <-
  aggregate_data_qd %>% filter(subsample == 'exploration')

val.aggregate_data_qd <-
  aggregate_data_qd %>% filter(subsample == 'validation')

inference.aggregate_data_qd <-
  aggregate_data_qd %>% filter(subsample == 'inference')

exploration.aggregate_data_qd <-
  exploration.aggregate_data_qd %>% tibble::rowid_to_column('orig_id')

val.aggregate_data_qd <-
  val.aggregate_data_qd %>% tibble::rowid_to_column('orig_id')

inference.aggregate_data_qd <-
  inference.aggregate_data_qd %>% tibble::rowid_to_column('orig_id')

exploration.aggregate_data_qd_sub <- exploration.aggregate_data_qd %>%
  filter(sex == levels(sex)[1], race_bin == levels(race_bin)[1], age_bin == levels(age_bin)[1], dual == levels(dual)[1])

pseudo_pop_fit <- generate_pseudo_pop(exploration.aggregate_data_qd[,'mortality'],  # why is the outcome required?
                                      exploration.aggregate_data_qd[,'pm25_ensemble'],
                                      #1:nrow(sub_pop),
                                      #sub_pop[,]
                                      select(exploration.aggregate_data_qd, confounders_names),
                                      ci_appr = "matching",
                                      pred_model = "sl",
                                      gps_model = "parametric",
                                      use_cov_transform = TRUE,
                                      transformers = list("pow2", "pow3"),
                                      bin_seq = a.vals,
                                      sl_lib = c("m_xgboost"),
                                      params = list('xgb_max_depth' = c(3,4,5),
                                                    'xgb_nrounds'=c(10,20,30,40,50,60),
                                                    #"xgb_nrounds"=50,
                                                    #"xgb_max_depth"=6,
                                                    "xgb_eta"=0.3,
                                                    "xgb_min_child_weight"=1),
                                      nthread=5, # number of cores, you can change,
                                      covar_bl_method = "absolute",
                                      covar_bl_trs = 0.3,
                                      covar_bl_trs_type = "mean",
                                      trim_quantiles = c(0,1), # trimed, you can change
                                      optimized_compile = FALSE, #created a column counter for how many times matched,
                                      max_attempt = 1,
                                      matching_fun = "matching_l1",
                                      delta_n = delta_n, # you can change this to the one you used in previous analysis,
                                      scale = 1.0)
exploration.aggregate_data_qd_sub.matched <- pseudo_pop_fit$pseudo_pop
# plot overall covariate balance
unmatched.explor_sub <- exploration.aggregate_data_qd_sub %>% 
  as.data.table
matched.explor_sub <- exploration.aggregate_data_qd_sub.matched %>% 
  as.data.table

unmatched_cor_sub <- 
  absolute_corr_fun(unmatched.explor_sub[,'pm25_ensemble'], unmatched.explor_sub[,c("mean_bmi","smoke_rate","hispanic","pct_blk","medhouseholdincome","medianhousevalue"  ,
                                                                            "poverty","education","popdensity","pct_owner_occ","summer_tmmx","winter_tmmx","summer_rmax" ,      
                                                                            "winter_rmax","region")])
matched_cor_sub <- 
  absolute_corr_fun(matched.explor_sub[,'w'], matched.explor_sub[,c("mean_bmi","smoke_rate","hispanic","pct_blk","medhouseholdincome","medianhousevalue"  ,
                                                                        "poverty","education","popdensity","pct_owner_occ","summer_tmmx","winter_tmmx","summer_rmax" ,      
                                                                        "winter_rmax","region")])

abs_cor_sub <- 
  rbind(data.frame(cov = c("mean_bmi","smoke_rate","hispanic","pct_blk","medhouseholdincome","medianhousevalue"  ,
                           "poverty","education","popdensity","pct_owner_occ","summer_tmmx","winter_tmmx","summer_rmax" ,      
                           "winter_rmax","region"),
                   cor = unmatched_cor_sub$absolute_corr, matched = FALSE), 
        data.frame(cov = c("mean_bmi","smoke_rate","hispanic","pct_blk","medhouseholdincome","medianhousevalue"  ,
                           "poverty","education","popdensity","pct_owner_occ","summer_tmmx","winter_tmmx","summer_rmax" ,      
                           "winter_rmax","region"),
                   #names(synth_data)[3:length(names(synth_data))],
                   cor = matched_cor_sub$absolute_corr, matched = TRUE))
ggplot(abs_cor_sub, aes(x = cov, y = cor, color = matched, group = matched)) + geom_point() + geom_line() + 
  theme(axis.text.x = element_text(angle = 45))
ggsave(paste0(dir_out, 'data_application_subset_balance.png'))

exploration.aggregate_data_qd.matched <- stratified_GPS_matching(exploration.aggregate_data_qd, delta_n = delta_n, bin_seq = a.vals,
                                                                 exposure_name = 'pm25_ensemble', confounders_names = confounders_names, 
                                                                 em_names = effect_modifiers, outcome_name = 'mortality')
exploration.aggregate_data_qd.matched <- 
  exploration.aggregate_data_qd.matched %>%
  mutate(id = row_number())

# plot overall covariate balance
unmatched.explor <- exploration.aggregate_data_qd %>% 
  as.data.table
matched.explor <- exploration.aggregate_data_qd.matched %>% 
  as.data.table

unmatched_cor <- 
  absolute_corr_fun(unmatched.explor[,'pm25_ensemble'], unmatched.explor[,c("mean_bmi","smoke_rate","hispanic","pct_blk","medhouseholdincome","medianhousevalue"  ,
                                                                            "poverty","education","popdensity","pct_owner_occ","summer_tmmx","winter_tmmx","summer_rmax" ,      
                                                                            "winter_rmax","region")])
matched_cor <- 
  absolute_corr_fun(matched.explor[,'pm25_ensemble'], matched.explor[,c("mean_bmi","smoke_rate","hispanic","pct_blk","medhouseholdincome","medianhousevalue"  ,
                                                                            "poverty","education","popdensity","pct_owner_occ","summer_tmmx","winter_tmmx","summer_rmax" ,      
                                                                            "winter_rmax","region")])

abs_cor <- 
  rbind(data.frame(cov = confounders_names,
                   cor = unmatched_cor$absolute_corr, matched = FALSE), 
        data.frame(cov = confounders_names,
                   #names(synth_data)[3:length(names(synth_data))],
                   cor = matched_cor$absolute_corr, matched = TRUE))
ggplot(abs_cor, aes(x = cov, y = cor, color = matched, group = matched)) + geom_point() + geom_line() + 
  theme(axis.text.x = element_text(angle = 45))


overall_effect1 <- lm(Y ~ w + 1, data = exploration.aggregate_data_qd.matched)$coefficients['w']

parms.used1 <- list(w  = exploration.aggregate_data_qd.matched$w,
                    # covariates = data.used.cont.cont[, 3:dim(data.used.cont.cont)[2]],
                    Y   = exploration.aggregate_data_qd.matched$Y,
                    overall_effect = overall_effect1)   
ulist.used <- define_u.list(stopping.rule = TRUE)
f <- as.formula('id ~ sex + dual + race_bin + age_bin')
lists1 <- create.sequence(exploration.aggregate_data_qd.matched, ulist.used, parms.used1, f)
