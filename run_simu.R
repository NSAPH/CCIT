
# wrapper function to simulate data, fit tree, and evaluate outcomes
run_simu <- function(lambdas, correct_splits, true_trt_effect_func, num_exposure_cats, noise.var, gps_spec = 1,
                     sample_size = 20000, em_spec = 1, beta, outcome_sd = 1, heterogenous_intercept = FALSE, stopping.rule, n_trials,
                     exploration.sample_covs = NULL, inference.sample_covs = NULL, val.sample_covs = NULL, 
                     matched.exploration.sample = NULL, matched.validation.sample = NULL, matched.inference.sample = NULL) {
  
  # if the matched covariate dataset (everything but the outcome) is not passed in, generate one
  if (is.null(matched.exploration.sample) | is.null(matched.validation.sample) | is.null(matched.inference.sample) | 
      is.null(exploration.sample_covs) | is.null(val.sample_covs) | is.null(inference.sample_covs)) {
    # generate covariate data
    synth_data_covs <- generate_syn_data_covs(sample_size = sample_size, gps_spec = gps_spec)
    
    # discretize treatment values for gps matching
    a.vals <- seq(min(synth_data_covs$treat), max(synth_data_covs$treat), length.out = num_exposure_cats)
    delta_n <- a.vals[2] - a.vals[1]
    
    synth_data_covs <- 
      synth_data_covs %>%
      mutate(treat_level = cut(treat, breaks = a.vals)) %>%
      mutate_at(vars(em1, em2, em3, em4), as.factor)
    
    # split data into subsamples, stratifying on exposure level bins
    synth_data_covs <- 
      split_dataset(data = synth_data_covs, exposure_bins = synth_data_covs$treat_level)
    
    exploration.sample_covs <-
      synth_data_covs %>% filter(subsample == 'exploration')
    
    val.sample_covs <-
      synth_data_covs %>% filter(subsample == 'validation')
    
    inference.sample_covs <-
      synth_data_covs %>% filter(subsample == 'inference')
    
    # GPS matching within subsamples and effect modifier groups
    matched.exploration.sample <- 
      stratified_GPS_matching(exploration.sample_covs, delta_n, bin_seq = a.vals, exposure_name = 'treat', 
                              confounders_names = c('cf1','cf2', 'cf3', 'cf4', 'cf5', 'cf6'), 
                              em_names = c('em1', 'em2', 'em3', 'em4'), 
                              outcome_name = NA) %>%
      mutate(id = row_number()) %>%
      #rename(w = treat) %>%
      mutate_at(vars(em1, em2, em3, em4), as.factor)
    
    matched.validation.sample <- 
      stratified_GPS_matching(val.sample_covs, delta_n, bin_seq = a.vals, exposure_name = 'treat', 
                              confounders_names = c('cf1','cf2', 'cf3', 'cf4', 'cf5', 'cf6'), 
                              em_names = c('em1', 'em2', 'em3', 'em4'), 
                              outcome_name = NA) %>%
      mutate(id = row_number()) %>%
      #rename(w = treat) %>%
      mutate_at(vars(em1, em2, em3, em4), as.factor)
    
    matched.inference.sample <- 
      stratified_GPS_matching(inference.sample_covs, delta_n, bin_seq = a.vals, exposure_name = 'treat', 
                              confounders_names = c('cf1','cf2', 'cf3', 'cf4', 'cf5', 'cf6'), 
                              em_names = c('em1', 'em2', 'em3', 'em4'), 
                              outcome_name = NA) %>%
      mutate(id = row_number()) %>%
      #rename(w = treat) %>%
      mutate_at(vars(em1, em2, em3, em4), as.factor)
  }
  
  results <- data.frame()
  if (n_trials > 0) {
    for (i in 1:n_trials){
      # generate outcome data
      exploration.sample_outcome <- 
        do.call(generate_syn_data_outcome, 
                c(as.list(exploration.sample_covs %>% select(-treat_level, -subsample)), beta = beta, em_spec = em_spec, outcome_sd = outcome_sd)) %>%
        tibble::rowid_to_column()
      
      validation.sample_outcome <- 
        do.call(generate_syn_data_outcome, 
                c(as.list(val.sample_covs %>% select(-treat_level, -subsample)), beta = beta, em_spec = em_spec, outcome_sd = outcome_sd)) %>%
        tibble::rowid_to_column()
      
      inference.sample_outcome <- 
        do.call(generate_syn_data_outcome, 
                c(as.list(inference.sample_covs %>% select(-treat_level, -subsample)), beta = beta, em_spec = em_spec, outcome_sd = outcome_sd)) %>%
        tibble::rowid_to_column()
      
      # add outcome data into matched samples
      matched.exploration.sample.outcomes <- 
        matched.exploration.sample %>% 
        select(-Y, -row_index) %>%
        left_join(exploration.sample_outcome, by = c('orig_id' = 'rowid', 'cf1', 'cf2', 'cf3', 'cf4', 'cf5', 'cf6', 'em1', 'em2', 'em3', 'em4',
                                                     'treat'))
      # p <- matched.exploration.sample.outcomes %>%
      #   ggplot(aes(x = w, y = Y, color = interaction(em1, em2))) + geom_point()
      
      matched.validation.sample.outcomes <- 
        matched.validation.sample %>% 
        select(-Y, -row_index) %>%
        left_join(validation.sample_outcome, by = c('orig_id' = 'rowid', 'cf1', 'cf2', 'cf3', 'cf4', 'cf5', 'cf6', 'em1', 'em2', 'em3', 'em4',
                                                    'treat'))
      matched.inference.sample.outcomes <- 
        matched.inference.sample %>% 
        select(-Y, -row_index) %>%
        left_join(inference.sample_outcome, by = c('orig_id' = 'rowid', 'cf1', 'cf2', 'cf3', 'cf4', 'cf5', 'cf6', 'em1', 'em2', 'em3', 'em4',
                                                   'treat'))
      
      overall_effect1 <- lm(Y ~ w + 1, data = matched.exploration.sample.outcomes)$coefficients['w']
      
      parms.used1 <- list(w  = matched.exploration.sample.outcomes$w,
                          # covariates = data.used.cont.cont[, 3:dim(data.used.cont.cont)[2]],
                          Y   = matched.exploration.sample.outcomes$Y,
                          overall_effect = overall_effect1)   
      ulist.used <- define_u.list(stopping.rule = stopping.rule)
      lists1 <- create.sequence(matched.exploration.sample, ulist.used, parms.used1)
      
      tree.list1 <- lists1[[1]]
      g.h.list1 <- lists1[[2]]
      tree.sizes1 <- sapply(tree.list1, function(t) {ifelse(is.null(t$splits), 0,nrow(t$splits))})
      
      # tree size selected via the validation sample (no cross-fitting for now)
      selected.tree.sizes <- data.frame()
      for (b in 1:1) {
        boot_sample_index <- sample(1:nrow(matched.validation.sample.outcomes), replace = TRUE)
        boot.matched.validation.sample <- matched.validation.sample.outcomes[boot_sample_index, ]
        
        complex.vals1 <- evaluate.sequence(tree.list1, boot.matched.validation.sample, matched.exploration.sample.outcomes, lambdas)
        
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
      
      #print(selected.tree.size$tree.size)
      #print(tree.sizes1)
      selected.trees <- lapply(selected.tree.size$tree.size, function(s) {tree.list1[[which(tree.sizes1 == s)]]})
      names(selected.trees) <- lambdas
      
      # determine which selected subgroup each inference observation is in. THere are a number of packages that do this, but the most recent versions
      # are not installed
      matched.inference.sample.subgroups <-  lapply(selected.trees, 
                                                    function(t) {
                                                      matched.inference.sample.new <- matched.inference.sample.outcomes
                                                      matched.inference.sample.new$subgroup <- NA
                                                      leaf_nodes <- t$frame %>% tibble::rownames_to_column() %>% filter(var == '<leaf>') %>% 
                                                        .$rowname %>% as.numeric()
                                                      for(n in leaf_nodes) {
                                                        rule <- path.rpart(t, n)
                                                        if(nrow(t$frame) == 1) {
                                                          ind <- 1:nrow(matched.inference.sample)
                                                        }
                                                        else {
                                                          rule_2 <- sapply(rule[[1]][-1], function(x) strsplit(x, '(?<=[><=])(?=[^><=])|(?<=[^><=])(?=[><=])', perl = TRUE))
                                                          ind <- apply(do.call(cbind, 
                                                                               lapply(rule_2, function(x) eval(call(ifelse(x[2] == '=', '==', x[2]), 
                                                                                                                    matched.inference.sample.outcomes[,x[1]], as.numeric(x[3]))))), 1, all)
                                                        }
                                                        matched.inference.sample.new[ind, 'subgroup'] <- n
                                                      }
                                                      matched.inference.sample.new <- matched.inference.sample.new %>%
                                                        mutate_at(vars(subgroup), as.factor)
                                                      return(matched.inference.sample.new)
                                                    })
      est.treatment.effects <- lapply(matched.inference.sample.subgroups, 
                                      function(df) {
                                        df %>%
                                          group_by(subgroup) %>%
                                          do(model = lm(Y ~ w + 1, data = .)) %>% mutate(est = summary(model)$coefficients['w', 1]) %>%
                                        right_join(df, by = 'subgroup')
                                        #print(df %>% filter(subgroup == 4))
                                        })
      # est.treatment.effects.dishonest <-  lapply(selected.trees,
      #                                  function(t) {mutate(matched.inference.sample, est = predict(t, newdata = matched.inference.sample))})
      # print(count(est.treatment.effects.dishonest[[1]], est))
      
      true_trt_effects <- true_trt_effect_func(matched.inference.sample.outcomes)
      
      # print(table(est.treatment.effects[[1]]$est - true_trt_effects$eff))
      # print(count(est.treatment.effects[[1]],est))
      # print(count(true_trt_effects,eff))
      # 
      # print(head(est.treatment.effects[[1]]$est))
      # print(head(true_trt_effects$eff))
      
      mse <- sapply(est.treatment.effects, function(est.treatment.effect) { mean((est.treatment.effect$est - true_trt_effects$eff)^2)})
      
      bias <- sapply(est.treatment.effects, function(est.treatment.effect) { mean((est.treatment.effect$est - true_trt_effects$eff))})
      
      # TRUE if tree is exactly correct - CIT paper uses the correct number of splits, but shouldn't order effect the correct number??
      selected.correct.splits <- sapply(selected.trees, function(t) {list(row.names(t$splits)) %in% correct_splits})
      
      
      # TRUE if one of the trees in the generated sequence is correct
      correct.tree.in.sequence <- any(sapply(tree.list1, function(t) {list(row.names(t$splits)) %in% correct_splits}))
      
      #print(sapply(tree.list1, function(t) {list(row.names(t$splits))}))
      # Number of noise variables selected
      numb.noise <- sapply(selected.trees, function(t) {sum(t$frame$var %in% noise.var)})
      
      iter_results <- 
        selected.tree.size %>%
        mutate(selected.correct.splits = selected.correct.splits,
               correct.tree.in.sequence = correct.tree.in.sequence,
               mse = mse, 
               bias = bias,
               numb.noise = numb.noise, 
               iter = i)
      results <- rbind(results, iter_results)
    }
  }
  return(list(results = results, matched.exploration.sample = matched.exploration.sample, 
              matched.inference.sample = matched.inference.sample,
              matched.validation.sample = matched.validation.sample,
              inference.sample_covs = inference.sample_covs,
              exploration.sample_covs = exploration.sample_covs, 
              val.sample_covs = val.sample_covs))
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

