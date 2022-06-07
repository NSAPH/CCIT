# stratified_GPS_matching <- function(data, delta_n, bin_seq = bin_seq) {
#   
#   data <- data %>% tibble::rowid_to_column('orig_id')
#   
#   matched <- data.frame()
#   
#   for(i in levels(data$em1)) {
#     for (j in levels(data$em2)) {
#       for (k in levels(data$em3)) {
#         for (l in levels(data$em4)) {
#           sub_pop <- data %>%
#             filter(em1 == i, em2 == j, em3 == k, em4 == l) %>%
#             tibble::rowid_to_column('new_id') #row index of sub-pop
#           pseudo_pop_fit <- generate_pseudo_pop(1:nrow(sub_pop), # why is the outcome required?
#                                                 sub_pop$treat, 
#                                                 sub_pop[, sapply(1:6, function(x) {paste0('cf',x)})],
#                                                 ci_appr = "matching",
#                                                 pred_model = "sl",
#                                                 gps_model = "parametric",
#                                                 use_cov_transform = TRUE,
#                                                 transformers = list("pow2", "pow3"),
#                                                 bin_seq = bin_seq,
#                                                 sl_lib = c("m_xgboost"),
#                                                 params = list('xgb_max_depth' = c(3,4,5),
#                                                               'xgb_nrounds'=c(10,20,30,40,50,60),
#                                                   #"xgb_nrounds"=50,
#                                                               #"xgb_max_depth"=6,
#                                                               "xgb_eta"=0.3,
#                                                               "xgb_min_child_weight"=1),
#                                                 nthread=5, # number of cores, you can change,
#                                                 covar_bl_method = "absolute",
#                                                 covar_bl_trs = 0.1,
#                                                 covar_bl_trs_type = "mean",
#                                                 trim_quantiles = c(0.5,0.95), # trimed, you can change
#                                                 optimized_compile = FALSE, #created a column counter for how many times matched,
#                                                 max_attempt = 1,
#                                                 matching_fun = "matching_l1",
#                                                 delta_n = delta_n, # you can change this to the one you used in previous analysis,
#                                                 scale = 1.0)
#           sub_pop_pseudo <- pseudo_pop_fit$pseudo_pop %>%
#             left_join(sub_pop, by = c('row_index' = 'new_id', 'cf1','cf2', 'cf3', 'cf4', 'cf5', 'cf6'))
#           # need to convert row_index to original row_index in data
#           matched <- rbind(matched, sub_pop_pseudo)
#           print(absolute_corr_fun(as.data.table(sub_pop)[,'treat'], as.data.table(sub_pop)[,c('cf1', 'cf2', 'cf3', 'cf4', 'cf5', 'cf6')]))
#           print(absolute_corr_fun(as.data.table(sub_pop_pseudo)[,'treat'], as.data.table(sub_pop_pseudo)[,c('cf1', 'cf2', 'cf3', 'cf4', 'cf5', 'cf6')]))
#         }
#       }
#     }
#   }
#   return(matched)
# }


stratified_GPS_matching <- function(data, delta_n, bin_seq, exposure_name, confounders_names, em_names, outcome_name) {
  
  if (length(em_names) > 0) {
    em <- em_names[1]
    em_lvls <- levels(data[,em])
    print(em)
    rbind(stratified_GPS_matching(data %>% filter(get(em) == em_lvls[1]), delta_n = delta_n, bin_seq = bin_seq, exposure_name = exposure_name, 
                            confounders_names = confounders_names, em_names = em_names[em != em_names], outcome_name = outcome_name),
    stratified_GPS_matching(data %>% filter(get(em) == em_lvls[2]), delta_n = delta_n, bin_seq = bin_seq, exposure_name = exposure_name, 
                            confounders_names = confounders_names, em_names = em_names[em != em_names], outcome_name = outcome_name))
    
  }
  else {
    if(nrow(data) > 0) {
      print(nrow(data))
      sub_pop <- data %>%
        tibble::rowid_to_column('new_id') #row index of sub-pop
      print(head(select(sub_pop, confounders_names)[,'cf1']))
      pseudo_pop_fit <- generate_pseudo_pop(ifelse(is.na(outcome_name), 1:nrow(sub_pop), sub_pop[,outcome_name]),  # why is the outcome required?
                                            sub_pop[,exposure_name],
                                            #1:nrow(sub_pop),
                                            #sub_pop[,]
                                            select(sub_pop, confounders_names),
                                            ci_appr = "matching",
                                            pred_model = "sl",
                                            gps_model = "parametric",
                                            use_cov_transform = TRUE,
                                            transformers = list("pow2", "pow3"),
                                            bin_seq = bin_seq,
                                            sl_lib = c("m_xgboost"),
                                            params = list('xgb_max_depth' = c(3,4,5),
                                                          'xgb_nrounds'=c(10,20,30,40,50,60),
                                                          #"xgb_nrounds"=50,
                                                          #"xgb_max_depth"=6,
                                                          "xgb_eta"=0.3,
                                                          "xgb_min_child_weight"=1),
                                            nthread=5, # number of cores, you can change,
                                            covar_bl_method = "absolute",
                                            covar_bl_trs = 0.1,
                                            covar_bl_trs_type = "mean",
                                            trim_quantiles = c(0.05,0.95), # trimed, you can change
                                            optimized_compile = FALSE, #created a column counter for how many times matched,
                                            max_attempt = 1,
                                            matching_fun = "matching_l1",
                                            delta_n = delta_n, # you can change this to the one you used in previous analysis,
                                            scale = 1.0)
      sub_pop_pseudo <- pseudo_pop_fit$pseudo_pop %>%
        left_join(sub_pop, by = c('row_index' = 'new_id', confounders_names))
      return(sub_pop_pseudo)
    }
    else {
      print('here')
      return(data.frame())}
  }
}
