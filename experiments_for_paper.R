source('/nfs/home/N/nik864/shared_space/ci3_analysis/nkc_subgroupId/HERF/Code/CCIT/evaluate.sequence.R')
source('/nfs/home/N/nik864/shared_space/ci3_analysis/nkc_subgroupId/HERF/Code/CCIT/run_simu.R')
source('/nfs/home/N/nik864/shared_space/ci3_analysis/nkc_subgroupId/HERF/Code/CCIT/create.sequence.R')
source('/nfs/home/N/nik864/shared_space/ci3_analysis/nkc_subgroupId/HERF/Code/CCIT/rpart_funcs.R')
source('/nfs/home/N/nik864/shared_space/ci3_analysis/nkc_subgroupId/HERF/Code/CCIT/stratified_GPS_matching.R')
source('/nfs/home/N/nik864/shared_space/ci3_analysis/nkc_subgroupId/HERF/Code/CCIT/split_dataset.R')
source('/nfs/home/N/nik864/shared_space/ci3_analysis/nkc_subgroupId/HERF/Code/CCIT/generate_synthetic_data_covs.R')
source('/nfs/home/N/nik864/shared_space/ci3_analysis/nkc_subgroupId/HERF/Code/CCIT/generate_synthetic_data_outcome.R')
library("devtools")
install_github("fasrc/CausalGPS", ref="master")
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

lambdas <- c(3000, 4000, 5000)

######################## Setting 1: 4 groups, same intercept
correct_splits <- list(c('em2', 'em1', 'em1'), c('em1', 'em2', 'em2'))
noise.var <- c('em3', 'em4')
beta <- 1

true_trt_effect_func <- function(df) {
  df %>%
    mutate_at(vars(em1, em2, em3, em4), function(e) {as.numeric(as.character(e))}) %>%
    mutate(eff = 10*(1 + beta * (0.5 * em1 - 0.8 * em2))) %>%
    mutate_at(vars(em1, em2, em3, em4), as.factor)
}

results <- run_simu(lambdas = lambdas, correct_splits = correct_splits, num_exposure_cats = 10, noise.var = noise.var, gps_spec = 1, em_spec = 1,
                    true_trt_effect_func = true_trt_effect_func, beta = beta,
                    sample_size = 2000, heterogenous_intercept = FALSE, stopping.rule = TRUE, n_trials = 0)#, 
                    exploration.sample_covs = exploration.sample_covs, val.sample_covs = val.sample_covs, inference.sample_covs = inference.sample_covs,
                    matched.exploration.sample = matched.exploration.sample, matched.validation.sample = matched.validation.sample,
                    matched.inference.sample = matched.inference.sample)
tree.data <- results$results

save(tree.data, file = paste0(dir_out,'selected.tree.data.with.stopping.setting1.lambda3000-5000.100trials.resampled.outcome.RData'))

################## Setting 2: 3 groups, same intercept
correct_splits <- list(c('em2', 'em1'))
noise.var <- c('em3', 'em4')

true_trt_effect_func <- function(df) {
  df %>%
    mutate_at(vars(em1, em2, em3, em4), function(e) {as.numeric(as.character(e))}) %>%
    mutate(eff = 10*(1 + 0.5 * em1 * em2 + 0.2 * em2)) %>%
    mutate_at(vars(em1, em2, em3, em4), as.factor)
  
}

results <- run_simu(lambdas = lambdas, correct_splits = correct_splits, num_exposure_cats = 10, noise.var = noise.var, gps_spec = 1, em_spec = 2,
                    true_trt_effect_func = true_trt_effect_func, beta = beta,
                    sample_size = 2000, heterogenous_intercept = FALSE, stopping.rule = TRUE, n_trials = 100, 
                    exploration.sample_covs = results$exploration.sample_covs, val.sample_covs = results$val.sample_covs, 
                    inference.sample_covs = results$inference.sample_covs,
                    matched.exploration.sample = results$matched.exploration.sample, matched.validation.sample = results$matched.validation.sample,
                    matched.inference.sample = results$matched.inference.sample)

tree.data <- results$results

save(tree.data, file = paste0(dir_out,'selected.tree.data.with.stopping.setting2.lambda3000-5000.100trials.resampled.outcome.RData'))

######################## Setting 3: 4 groups, different intercept
correct_splits <- list(c('em2', 'em1', 'em1'), c('em1', 'em2', 'em2'))
noise.var <- c('em3', 'em4')

true_trt_effect_func <- function(df) {
  df %>%
    mutate_at(vars(em1, em2, em3, em4), function(e) {as.numeric(as.character(e))}) %>%
    mutate(eff = 10*(1 + 0.5 * em1 - 0.8 * em2)) %>%
    mutate_at(vars(em1, em2, em3, em4), as.factor)
  
}

results <- run_simu(lambdas = lambdas, correct_splits = correct_splits, num_exposure_cats = 10, noise.var = noise.var, gps_spec = 1, em_spec = 1,
                    true_trt_effect_func = true_trt_effect_func, beta = beta,
                    sample_size = 2000, heterogenous_intercept = TRUE, stopping.rule = TRUE, n_trials = 100, 
                    exploration.sample_covs = results$exploration.sample_covs, val.sample_covs = results$val.sample_covs, 
                    inference.sample_covs = results$inference.sample_covs,
                    matched.exploration.sample = results$matched.exploration.sample, matched.validation.sample = results$matched.validation.sample,
                    matched.inference.sample = results$matched.inference.sample)

tree.data <- results$results

save(tree.data, file = paste0(dir_out,'selected.tree.data.with.stopping.setting3.lambda3000-5000.100trials.resampled.outcome.RData'))

################## Setting 4: 3 groups, different intercept
correct_splits <- list(c('em2', 'em1'))
noise.var <- c('em3', 'em4')

true_trt_effect_func <- function(df) {
  df %>%
    mutate_at(vars(em1, em2, em3, em4), function(e) {as.numeric(as.character(e))}) %>%
    mutate(eff = 10*(1 + 0.5 * em1 * em2 + 0.2 * em2)) %>%
    mutate_at(vars(em1, em2, em3, em4), as.factor)
  
}

results <- run_simu(lambdas = lambdas, correct_splits = correct_splits, num_exposure_cats = 10, noise.var = noise.var, gps_spec = 1, em_spec = 2,
                    true_trt_effect_func = true_trt_effect_func, beta = beta,
                    sample_size = 2000, heterogenous_intercept = TRUE, stopping.rule = TRUE, n_trials = 100, 
                    exploration.sample_covs = results$exploration.sample_covs, val.sample_covs = results$val.sample_covs, 
                    inference.sample_covs = results$inference.sample_covs,
                    matched.exploration.sample = results$matched.exploration.sample, matched.validation.sample = results$matched.validation.sample,
                    matched.inference.sample = results$matched.inference.sample)

tree.data <- results$results

save(tree.data, file = paste0(dir_out,'selected.tree.data.with.stopping.setting4.lambda3000-5000.100trials.resampled.outcome.RData'))


########################## Setting 5: no effect modification
correct_splits <- list(c())
noise.var <- c('em1', 'em2', 'em3', 'em4')
true_trt_effect_func <- function(df) {
  df %>%
    mutate(eff = 10)
}

results <- run_simu(lambdas = lambdas, correct_splits = correct_splits, num_exposure_cats = 10, noise.var = noise.var, gps_spec = 1, em_spec = 0,
                    true_trt_effect_func = true_trt_effect_func, beta = beta,
                    sample_size = 2000, heterogenous_intercept = FALSE, stopping.rule = TRUE, n_trials = 100, 
                    exploration.sample_covs = results$exploration.sample_covs, val.sample_covs = results$val.sample_covs, 
                    inference.sample_covs = results$inference.sample_covs,
                    matched.exploration.sample = results$matched.exploration.sample, matched.validation.sample = results$matched.validation.sample,
                    matched.inference.sample = results$matched.inference.sample)

tree.data <- results$results

save(tree.data, file = paste0(dir_out,'selected.tree.data.with.stopping.setting5.lambda3000-5000.100trials.resampled.outcome.RData'))

###################### Combining results for different settings
selected.tree.data.combined <- data.frame()
for (i in 1:5){
  load(file = paste0(dir_out,'selected.tree.data.with.stopping.setting', i, '.lambda3000-5000.100trials.resampled.outcome.RData'))
  selected.tree.data.combined <- selected.tree.data.combined %>% rbind(tree.data %>% cbind(setting = i))
}

setting.labs <- c('Setting 1', 'Setting 2', 'Setting 3', 'Setting 4', 'Setting 5')
names(setting.labs) <- c(1,2,3,4,5)

selected.tree.data.combined %>%
  filter(lambda == 3000) %>%
  tidyr::gather("metric", "value", c(bias, mse)) %>%
  ggplot(aes(y = value, group = lambda)) + geom_boxplot() + 
  facet_grid(rows = vars(metric), cols = vars(setting), labeller = labeller(setting = setting.labs), scales = 'free_y') + 
  labs(y=NULL) + 
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
ggsave(file = paste0(dir_out, 'mse.biase.with.stopping.png'))


############################ Setting 1 with Different Beta Values
betas <- c(0.15,0.2,0.3,1)
for (beta in betas) {
  true_trt_effect_func <- function(df) {
    df %>%
      mutate_at(vars(em1, em2, em3, em4), function(e) {as.numeric(as.character(e))}) %>%
      mutate(eff = 10*(1 + beta * (0.5 * em1 - 0.8 * em2))) %>%
      mutate_at(vars(em1, em2, em3, em4), as.factor)
  }
  
  results <- run_simu(lambdas = lambdas, correct_splits = correct_splits, num_exposure_cats = 10, noise.var = noise.var, gps_spec = 1, em_spec = 1,
                      true_trt_effect_func = true_trt_effect_func, beta = beta,
                      sample_size = 2000, heterogenous_intercept = FALSE, stopping.rule = TRUE, n_trials = 100, 
                      exploration.sample_covs = results$exploration.sample_covs, val.sample_covs = results$val.sample_covs, 
                      inference.sample_covs = results$inference.sample_covs,
                      matched.exploration.sample = results$matched.exploration.sample, matched.validation.sample = results$matched.validation.sample,
                      matched.inference.sample = results$matched.inference.sample)
  tree.data <- results$results
  save(tree.data, file = paste0(dir_out,'selected.tree.data.with.stopping.setting1.lambda3000-5000.100trials.resampled.outcome.beta', beta, '.RData'))
}


## results for setting 1 with different beta
selected.tree.data.combined <- data.frame()
for (beta in betas){
  load(file = paste0(dir_out,'selected.tree.data.with.stopping.setting1.lambda3000-5000.100trials.resampled.outcome.beta', beta, '.RData'))
  selected.tree.data.combined <- selected.tree.data.combined %>% rbind(selected.tree.data %>% cbind(beta = beta))
}

beta.labs <- sapply(betas, function(x) {paste0('Beta = ', x)})
names(beta.labs) <- as.character(betas)

selected.tree.data.combined %>%
  filter(lambda == 3000) %>%
  tidyr::gather("metric", "value", c(bias, mse)) %>%
  ggplot(aes(y = value, group = lambda)) + geom_boxplot() + 
  facet_grid(rows = vars(metric), cols = vars(beta), labeller = labeller(beta = beta.labs), scales = 'free_y') + 
  labs(y=NULL) + 
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
ggsave(file = paste0(dir_out, 'mse.bias.with.stopping.beta.png'))

selected.tree.data.combined %>%
  filter(lambda == 3000) %>%
  group_by(beta) %>%
  summarise(percent.correct.tree = mean(selected.correct.splits))

# selected.tree.data.combined %>%
#   filter(lambda == 3000) %>%
#   #filter(selected.correct.splits) %>%
#   ggplot(aes(y = mse, group = lambda)) + geom_boxplot() + 
#   facet_wrap(vars(beta), nrow = 1) + 
#   scale_x_continuous(breaks = lambdas)
# ggsave(file = paste0(dir_out, 'mse.with.stopping.beta.png'))
# 
# selected.tree.data.combined %>%
#   filter(lambda == 3000) %>%
#   ggplot(aes(y = bias, group = lambda)) + geom_boxplot() + 
#   facet_wrap(vars(beta), nrow = 1) + 
#   scale_x_continuous(breaks = lambdas)
# ggsave(file = paste0(dir_out, 'bias.with.stopping.beta.png'))
# 
# 
# 
# selected.tree.data.combined %>%
#   filter(lambda == 3000) %>%
#   ggplot(aes(y = mse, group = lambda)) + geom_boxplot() + 
#   facet_wrap(vars(setting), nrow = 1, labeller = labeller(setting = setting.labs)) + 
#   scale_x_continuous(breaks = lambdas)
# ggsave(file = paste0(dir_out, 'mse.with.stopping2.png'))
# 
# selected.tree.data.combined %>%
#   filter(lambda == 3000) %>%
#   ggplot(aes(y = bias, group = lambda)) + geom_boxplot() + 
#   facet_wrap(vars(setting), nrow = 1, labeller = labeller(setting = setting.labs)) + 
#   scale_x_continuous(breaks = lambdas)
# ggsave(file = paste0(dir_out, 'bias.with.stopping2.png'))
# 
# selected.tree.data.combined %>% group_by(setting) %>%
#   summarise_at(vars(numb.noise, selected.correct.splits, bias, mse), mean)
# 
# 
# 
# 
# 
# 
# 
# load(paste0(dir_out,'no.stopping.RData'))
# 
# 
# selected.tree.data %>%
#   ggplot(aes(x = lambda, y = mse, group = lambda)) + geom_boxplot() + 
#   scale_x_continuous(breaks = lambdas)
# ggsave(paste0(dir_out, 'MSE.no.stopping100trials.png'))
# 
# 
# 
# 
# 
# 
# 
# selected.tree.data %>%
#   ggplot(aes(x = lambda, fill = selected.correct.splits, y = ..count../100)) + geom_histogram(position = 'dodge') + 
#   scale_x_continuous(breaks = lambdas) + 
#   labs(y = 'Proportion of Simulation trials') 
# ggsave(paste0(dir_out, 'Prop.correct.trees.no.stopping.png'))
# 
# selected.tree.data %>%
#   ggplot(aes(x = lambda, y = numb.noise, group = lambda)) + geom_boxplot() + 
#   scale_x_continuous(breaks = lambdas) + 
#   labs(y = 'Proportion of Simulation trials')
# ggsave(paste0(dir_out, 'Prop.correct.trees.no.stopping.png'))
# 
# selected.tree.data %>%
#   mutate('Pruned to Correct Size' = (floor(tree.size) == 3)) %>%
#   ggplot(aes(x = lambda, fill = `Pruned to Correct Size`, y = ..count../100)) + geom_histogram(position = 'dodge') + 
#   scale_x_continuous(breaks = lambdas) + 
#   labs(y = 'Proportion of Simulation trials')
# ggsave(paste0(dir_out, 'Pruning_accuracy.no.stopping.png'))
# ggsave(paste0(dir_out, 'Pruning_accuracy_bootstrap_n100000.png'))
# 
# 
# selected.tree.data %>%
#   ggplot(aes(x = lambda, y = floor(tree.size), group = lambda)) + geom_boxplot() + 
#   geom_hline(yintercept = 3) + 
#   scale_x_continuous(breaks = lambdas)
# ggsave(paste0(dir_out, 'Pruned_tree_size_crossfit.png'))
# ggsave(paste0(dir_out, 'Pruned_tree_size_bootstrapn100000.png'))
# 
# 
# 
# selected.tree.data %>%
#   ggplot(aes(x = lambda, y = selected.tree.size, group = lambda)) + geom_boxplot() + 
#   scale_x_continuous(breaks = lambdas)
# ggsave(paste0(dir_out, 'Pruned_tree_size.png'))
# 
# 
# selected.tree.data %>%
#   mutate('Pruned to Correct Size' = as.logical(selected.correct.splits)) %>%
#   ggplot(aes(x = lambda, fill = `Pruned to Correct Size`)) + geom_histogram(position = 'dodge') + 
#   scale_x_continuous(breaks = lambdas)
# ggsave(paste0(dir_out, 'Pruning_accuracy.png'))
# 
# 
# selected.tree.data %>%
#   group_by(lambda) %>%
#   summarise(prop_correct = sum(selected.correct.splits)/n(), 
#             avg_tree_size = mean(selected.tree.size), 
#             prop_correct_in_sequence = sum(correct.tree.in.sequence)/n(),
#             prop_correct_size = sum(selected.tree.size == 3)/n())
# 
# ret <- evaluate_simulated_continuous_outcome(lambda = 1, correct_splits = correct_splits, num_exposure_cats = 20)
