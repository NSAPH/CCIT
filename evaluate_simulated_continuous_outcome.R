library(dplyr)
library(truncnorm)
library(CausalGPS)
library(caret)
library(rpart)

dir_out = '/nfs/home/N/nik864/shared_space/ci3_analysis/nkc_subgroupId/HERF/'

# wrapper function to simulate data, fit tree, and evaluate outcomes
evaluate_simulated_continuous_outcome <- function(lambdas, correct_splits, noise.var, num_exposure_cats, gps_spec = 1, em_spec = 1, 
                                                  sample_size = 20000, heterogenous_intercept = FALSE, stopping.rule) {
  synth_data <- generate_syn_data_het(sample_size = 20000, outcome_type = 'continuous',
                                       gps_spec = 1, em_spec = em_spec, cova_spec = 1, heterogenous_intercept = heterogenous_intercept, 
                                       em_as_confounder = FALSE, outcome_sd = 1, beta = 0.3)
  
  # discretize treatment values for gps matching
  a.vals <- seq(min(synth_data$treat), max(synth_data$treat), length.out = num_exposure_cats)
  delta_n <- a.vals[2] - a.vals[1]
  synth_data <- 
    synth_data %>%
    mutate(treat_level = cut(treat, breaks = a.vals))

  synth_data <- 
    split_dataset(data = synth_data, exposure_bins = synth_data$treat_level)
  
  exploration.sample <- 
    synth_data %>% filter(subsample == 'exploration')
  
  val.sample <- 
    synth_data %>% filter(subsample == 'validation')
  
  inference.sample <- 
    synth_data %>% filter(subsample == 'inference')
  
  matched.exploration.sample <- 
    stratified_GPS_matching(exploration.sample, delta_n) %>%
    mutate(id = row_number()) %>%
    #rename(w = treat) %>%
    mutate_at(vars(em1, em2, em3, em4), as.factor)
  
  matched.validation.sample <- 
    stratified_GPS_matching(val.sample, delta_n) %>%#
    mutate(id = row_number()) %>%
    #rename(w = treat) %>%
    mutate_at(vars(em1, em2, em3, em4), as.factor)
    
  overall_effect1 <- lm(Y ~ w + 1, data = matched.exploration.sample)$coefficients['w']
    
  parms.used1 <- list(w  = matched.exploration.sample$w,
                     # covariates = data.used.cont.cont[, 3:dim(data.used.cont.cont)[2]],
                     Y   = matched.exploration.sample$Y,
                     overall_effect = overall_effect1)   
  ulist.used <- define_u.list(stopping.rule = stopping.rule)
  lists1 <- create.sequence(matched.exploration.sample, ulist.used, parms.used1)

  tree.list1 <- lists1[[1]]
  g.h.list1 <- lists1[[2]]
  tree.sizes1 <- sapply(tree.list1, function(t) {ifelse(is.null(t$splits), 0,nrow(t$splits))})
  
  # tree size selected via the validation sample (no cross-fitting for now)
  selected.tree.sizes <- data.frame()
  for (b in 1:1) {
    boot_sample_index <- sample(1:nrow(matched.validation.sample), replace = TRUE)
    boot.matched.validation.sample <- matched.validation.sample[boot_sample_index, ]
    
    complex.vals1 <- evaluate.sequence(tree.list1, boot.matched.validation.sample, matched.exploration.sample, lambdas)
    
    selected.tree.size1 <- complex.vals1 %>% group_by(lambda) %>%
      filter(complex.val == max(complex.val)) %>%
      select(tree.size) %>%
      ungroup()
    selected.tree.sizes <- rbind(selected.tree.sizes, selected.tree.size1)
    
  }
  
  selected.tree.size <- selected.tree.sizes %>%
    group_by(lambda) %>%
    summarise_at(vars(tree.size), function(s) {round(mean(s))}) %>%
    ungroup()
  
  print(selected.tree.size$tree.size)
  print(tree.sizes1)
  selected.trees <- lapply(selected.tree.size$tree.size, function(s) {tree.list1[[which(tree.sizes1 == s)]]})
  names(selected.trees) <- lambdas
  est.treatment.effects <-  lapply(selected.trees, 
                                   function(t) {mutate(inference.sample, pred = predict(t, newdata = inference.sample))})
    
  
  true_trt_effects <- true_trt_effect_func(inference.sample)
  
  mse <- sapply(est.treatment.effects, function(est.treatment.effect) { mean((est.treatment.effect$pred - true_trt_effects$eff)^2)})
  
  bias <- sapply(est.treatment.effects, function(est.treatment.effect) { mean((est.treatment.effect$pred - true_trt_effects$eff))})
  
  # TRUE if tree is exactly correct - CIT paper uses the correct number of splits, but shouldn't order effect the correct number??
  selected.correct.splits <- sapply(selected.trees, function(t) {list(row.names(t$splits)) %in% correct_splits})
    
  
  # TRUE if one of the trees in the generated sequence is correct
  correct.tree.in.sequence <- any(sapply(tree.list1, function(t) {list(row.names(t$splits)) %in% correct_splits}))
  
  # Number of noise variables selected
  numb.noise <- sapply(selected.trees, function(t) {sum(t$frame$var %in% noise.var)})
  
  return(selected.tree.size %>%
           mutate(selected.correct.splits = selected.correct.splits,
                  correct.tree.in.sequence = correct.tree.in.sequence,
                  mse = mse, 
                  bias = bias,
                  numb.noise = numb.noise))
}

# complex.vals1 <- data.frame(complex.val = evaluate.sequence(tree.list1, matched.validation.sample, matched.exploration.sample, lambda), 
#            pruning_step = 1:length(tree.list1),
#            tree.size = tree.sizes1)
# 
# selected.tree.size1 <- complex.vals1[which(complex.vals1$complex.val == max(complex.vals1$complex.val)), 'tree.size']
# 
# overall_effect2 <- lm(Y ~ w + 1, data = matched.validation.sample)$coefficients['w']
# 
# parms.used2 <- list(w  = matched.validation.sample$w,
#                     # covariates = data.used.cont.cont[, 3:dim(data.used.cont.cont)[2]],
#                     Y   = matched.validation.sample$Y,
#                     overall_effect = overall_effect2)   
# 
# lists2 <- create.sequence(matched.validation.sample, ulist.used, parms.used2)
# 
# tree.list2 <- lists2[[1]]
# g.h.list2 <- lists2[[2]]
# 
# tree.sizes2 <- sapply(tree.list2, function(t) {nrow(t$splits)})
# complex.vals2 <- data.frame(complex.val = evaluate.sequence(tree.list2, matched.exploration.sample, matched.validation.sample, lambda), 
#                             pruning_step = 1:length(tree.list2),
#                             tree.size = tree.sizes2)
# 
# selected.tree.size2 <- complex.vals2[which(complex.vals2$complex.val == max(complex.vals2$complex.val)), 'tree.size']

#selected.tree.size <- mean(c(selected.tree.size1, selected.tree.size2))

# pps calculation doesn't complete 
# pps = 1
# for(i in 1:(nrow(inference.sample)-1)){
#   for(j in (i+1):nrow(inference.sample)){
#     a=b=0
#     if(true_trt_effects$eff[i] == true_trt_effects$eff[j]){a = 1}
#     
#     # Treat observations in the 1-observation terminal node be in the same node
#     if(is.na(est.treatment.effects$pred[i]) | is.na(est.treatment.effects$pred[j])) {
#       b = 1
#     } else if (est.treatment.effects$pred[i] == est.treatment.effects$pred[j]) {
#       b = 1
#     }
#     pps = pps - abs(a-b)/choose(nrow(inference.sample), 2)
#   }
# }

