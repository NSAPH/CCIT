#' @title
#' Generate Synthetic Data for CausalGPS Package - Modified by Nina
#'
#' @description
#' Generates synthetic data set based on different GPS models and covariates.
#'
#'Required Libraries: truncnorm, dplyr
generate_syn_data_covs <- function(sample_size=10000, 
                                  gps_spec = 1) {
  
  if (sample_size < 0 || !is.numeric(sample_size)){
    stop("'sample_size' should be a positive ineteger numer.")
  }
  
  size <- sample_size
  
  #pre-treatment variables (confounders)
  cf  <- MASS::mvrnorm(n = size,
                       mu = c(0,0,0,0),
                       Sigma = matrix(c(1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1),
                                      ncol=4))
  colnames(cf) <- c('cf1', 'cf2', 'cf3', 'cf4')
  
  cf5 <- sample(c((-2):2), size, replace = TRUE)
  cf6 <- stats::runif(size, min=-3, max=3)
  
  # true effect modifiers
  em1 <- rbinom(size, 1, 0.5)
  em2 <- rbinom(size, 1, 0.4)
  
  # null effect modifiers
  em3 <- rbinom(size, 1, 0.3)
  em4 <- rbinom(size, 1, 0.2)
  if (gps_spec == 1) {
    # modified Xiao's propensity function to ensure treatment >= 0. Note that mu is not the expected propensity.
    mu <- 3
    treat <- rtruncnorm(size,a=0,b=Inf,mean=mu,sd=5)
    # currently not using other gps_spec values
  } else if (gps_spec == 2) {
    mu <- (- 0.8 + 0.1 * cf[ ,1] + 0.1 * cf[ ,2] - 0.1 * cf[ ,3]
           + 0.2 * cf[ ,4] + 0.1 * cf5 + 0.1 * cf6) * 9 + 17 
    #treat <- rnorm(size,mean=mu,sd=5)
    
    treat <- rtruncnorm(size,a=0,b=Inf,mean=mu,sd=5)
    # treat <- ((- 0.8 + 0.1 * cf[ ,1] + 0.1 * cf[ ,2] - 0.1 * cf[ ,3]
    #            + 0.2 * cf[ ,4] + 0.1 * cf5 + 0.1 * cf6) * 15 + 22 + stats::rt(size,2))
    # 
    # treat[which(treat < (-5))] <- (-5)
    # treat[which(treat > (25))] <- (25)
    
  } else if (gps_spec == 3) {
    
    treat <- ((- 0.8 + 0.1 * cf[ , 1] + 0.1 * cf[ , 2]- 0.1 *cf[ ,3] + 0.2 * cf [ , 4]
               + 0.1 * cf5 + 0.1 * cf6) * 9
              + 1.5 * cf[ , 3] ^ 2 + stats::rnorm(size, mean = 0, 5) + 15)
    
  } else if (gps_spec == 4) {
    
    treat <- (49 * exp((-0.8 + 0.1 * cf[ ,1] + 0.1 * cf[ , 2] - 0.1 * cf[ , 3]
                        + 0.2 * cf[ , 4] + 0.1 * cf5 + 0.1 * cf6))
              / (1 + exp((-0.8 + 0.1 * cf[,1] + 0.1 * cf[ , 2] - 0.1 * cf[ , 3]
                          + 0.2 * cf[ , 4] + 0.1 * cf5 + 0.1 * cf6))) - 6 + stats::rnorm(size, sd=5))
    
  } else if (gps_spec == 5) {
    
    treat <- (42 / (1 + exp((-0.8 + 0.1 * cf[ , 1] + 0.1 * cf[ , 2]- 0.1 * cf[ , 3]
                             + 0.2 * cf[,4] + 0.1 * cf5 + 0.1 * cf6))) - 18 + stats::rnorm(size,sd=5))
    
  } else if (gps_spec == 6) {
    
    treat <- (log(abs(-0.8 + 0.1 * cf[ , 1] + 0.1 * cf[ , 2] - 0.1 * cf[ , 3]
                      + 0.2 * cf[ , 4] + 0.1 * cf5 + 0.1 * cf6)) * 7 + 13 + stats::rnorm(size,sd=4))
    
  } else if (gps_spec == 7) {
    
    treat <- ((-0.8 + 0.1 * cf[,1] + 0.1 * cf[,2] - 0.1 * cf[,3] + 0.2 * cf[,4]
               + 0.1 * cf5 + 0.1 * cf6) * 15 + 22 + stats::rt(size,2)) #+ rcauchy(size)
  } else {
    
    stop(paste("gps_spec: ", gps_spec, ", is not a valid value."))
    
  }
  
  # not currently using this
  # if (em_as_confounder) {
  #   treat <- treat + em1*3
  # }
  
  simulated_covariates <- as.data.frame(cf) %>% 
    mutate(cf5 = cf5, cf6 = cf6, em1 = em1, em2 = em2, em3 = em3, em4 = em4, treat = treat)
    
    # list(cf = cf, cf5 = cf5, cf6 = cf6, em1 = em1, em2 = em2, em3 = em3, em4 = em4, treat = treat)
    # colnames(simulated_covariates)[2:11]<-c("cf1","cf2","cf3","cf4","cf5","cf6","em1","em2","em3","em4")
  
  return(simulated_covariates)
}


