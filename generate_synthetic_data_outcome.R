#' @title
#' Generate Synthetic Data for CausalGPS Package - Modified by Nina
#'
#' @description
#' Generates synthetic data set based on different GPS models and covariates.
#'
#'Required Libraries: truncnorm, dplyr
generate_syn_data_outcome <- function(cf1, cf2, cf3, cf4, cf5, cf6, em1, em2, em3, em4, treat, 
                                      outcome_sd = 1,
                              em_spec = 1, heterogenous_intercept = FALSE, 
                              beta = 1) {
  cf <- as.matrix(cbind(cf1, cf2, cf3, cf4))
  # produce outcome Y
  #mu <- as.numeric()
  Y <- as.numeric()
  em1 <- as.numeric(as.character(em1))
  em2 <- as.numeric(as.character(em2))
  em3 <- as.numeric(as.character(em3))
  em4 <- as.numeric(as.character(em4))
  # em2 = 0, em1 = 0 -> 10
  # em2 = 0, em1 = 1 -> 15
  # em2 = 1, em1 = 1 -> 7
  # em2 = 1, em2 = 0 -> 2
  if (em_spec == 0) {
    mu <- -10 - sum(c(2,2,3,-1) * cf) - 2 * (cf5 + cf6) + 10 * treat
  }
  if (em_spec == 1) {
    mu <- -10 - sum(c(2, 2, 3, -1) * cf) - 2 * (cf5 + cf6) +
      10 * (1 + beta * (0.5 * em1- 0.8 * em2)) * treat
  }
  else if (em_spec == 2) {
    mu <- -10 - sum(c(2, 2, 3, -1) * cf) - 2 * (cf5 + cf6) +
      10 * (1 + beta * (0.5 * em1 * em2 + 0.2 * em2)) * treat
  }
  if (heterogenous_intercept) {
    mu <- mu + 20 * (em1 + 2*em2)
  }
  Y <- mu + stats::rnorm(length(treat), mean=0, sd=outcome_sd)

  # for (i in 1:length(treat)) {
  #   if (em_spec == 0) {
  #     mu[i] <- -10 - sum(c(2,2,3,-1) * cf[i, ]) - 2 * (cf5[i] + cf6[i]) + 10 * treat[i]
  #   }
  #   if (em_spec == 1) {
  #     mu[i] <- -10 - sum(c(2, 2, 3, -1) * cf[i, ]) - 2 * (cf5[i] + cf6[i]) + 
  #       10 * (1 + beta * (0.5 * em1[i] - 0.8 * em2[i])) * treat[i]
  #   }
  #   else if (em_spec == 2) {
  #     mu[i] <- -10 - sum(c(2, 2, 3, -1) * cf[i, ]) - 2 * (cf5[i] + cf6[i]) + 
  #       10 * (1 + beta * (0.5 * em1[i] * em2[i] + 0.2 * em2[i])) * treat[i]
  #   }
  #   if (heterogenous_intercept) {
  #     mu[i] <- mu[i] + 20 * (em1[i] + 2*em2[i])
  #   }
  #   Y[i] <- mu[i] + stats::rnorm(1, mean=0, sd=outcome_sd)
  # }
  # if (outcome_type == 'binary') {
  #   simulated_data <- 
  #     as.data.frame(cbind(treat_level = treat_df$treat_level, avg_treat = treat_df$avg_treat, 
  #                         mu, em1, em2, em3, em4)) %>%
  #     mutate(treat_level = cut(treat, breaks = 20)) %>%
  #     group_by(em1, em2, em3, em4, treat_level, avg_treat, mu) %>%
  #     summarise(Y = rpois(1, n()*mu), 
  #                         n = n()) %>%
  #     mutate_at(vars(em1, em2, em3, em4), as.factor)
  #   #colnames(simulated_data)[2:12]<-c("cf1","cf2","cf3","cf4","cf5","cf6","em1","em2","em3","em4")
  # }

    simulated_data<-data.frame(cbind(Y,treat,cf, cf5, cf6, em1, em2, em3, em4)) %>%
      mutate_at(vars(em1, em2, em3, em4), as.factor)
    colnames(simulated_data)[3:12]<-c("cf1","cf2","cf3","cf4","cf5","cf6","em1","em2","em3","em4")

  return(simulated_data)
}
 

