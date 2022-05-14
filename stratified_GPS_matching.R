stratified_GPS_matching <- function(data, delta_n, bin_seq = bin_seq) {
  
  data <- data %>% tibble::rowid_to_column('orig_id')
  
  matched <- data.frame()
  
  for(i in levels(data$em1)) {
    for (j in levels(data$em2)) {
      for (k in levels(data$em3)) {
        for (l in levels(data$em4)) {
          sub_pop <- data %>%
            filter(em1 == i, em2 == j, em3 == k, em4 == l) %>%
            tibble::rowid_to_column('new_id')
          pseudo_pop_fit <- generate_pseudo_pop(1:nrow(sub_pop), # why is the outcome required?
                                                sub_pop$treat, 
                                                sub_pop[, sapply(1:6, function(x) {paste0('cf',x)})],
                                                ci_appr = "matching",
                                                pred_model = "sl",
                                                gps_model = "parametric",
                                                use_cov_transform = FALSE,
                                                #transformers = list("pow2", "pow3"),
                                                bin_seq = bin_seq,
                                                sl_lib = c("m_xgboost"),
                                                params = list("xgb_nrounds"=50,
                                                              "xgb_max_depth"=6,
                                                              "xgb_eta"=0.3,
                                                              "xgb_min_child_weight"=1),
                                                nthread=5, # number of cores, you can change,
                                                covar_bl_method = "absolute",
                                                covar_bl_trs = 0.5,
                                                covar_bl_trs_type = "mean",
                                                trim_quantiles = c(0,1), # trimed, you can change
                                                optimized_compile = FALSE, #created a column counter for how many times matched,
                                                max_attempt = 1,
                                                matching_fun = "matching_l1",
                                                delta_n = delta_n, # you can change this to the one you used in previous analysis,
                                                scale = 1.0)
          sub_pop_pseudo <- pseudo_pop_fit$pseudo_pop %>%
            left_join(sub_pop, by = c('row_index' = 'new_id', 'cf1','cf2', 'cf3', 'cf4', 'cf5', 'cf6'))
          # need to convert row_index to original row_index in data
          matched <- rbind(matched, sub_pop_pseudo)
        }
      }
    }
  }
  return(matched)
}
