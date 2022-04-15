# define functions for custom rpart function
# required libraries: rpart, rpart.plot
define_u.list <- function(stopping.rule = TRUE) {
  criterion <- function(Y,w) {
    lmod <- lm(Y ~ w + 1)
    lmod$coefficients['w']
  }
  # initialization function for rpart
  itemp <- function(y, offset, parms, wt) {
    
    sfun <- function(yval, dev, wt, ylevel, digits ) {
      paste(" mean=", format(signif(yval, digits)),
            ", MSE=" , format(signif(dev/wt, digits)),
            sep = '')
    }
    
    environment(sfun) <- .GlobalEnv
    list(y = c(y), parms = parms, numresp = 1, numy = 1, summary = sfun)
  }
  
  # evaluation function for rpart - I don't think I can compute the complexity measure as defined by Su without x. I could pass as a parm
  etemp <- function(y, wt, parms) {
    Y = parms$Y[y]
    w = parms$w[y]
    list(label = criterion(Y,w), deviance = sd(Y))
  }
  
  # splitting function for rpart
  stemp <- function(y, wt, x, parms, continuous){
    # only 1 possible split per covariate
    Y <- parms$Y[y]
    w <- parms$w[y]
    EM <- x
    ux <- sort(unique(x))
    means <- data.frame(Y = Y, EM = EM) %>% group_by(EM) %>% summarise_at(vars(Y), mean) %>% .$Y
    ord <- order(means)
    lmod <- lm(Y ~ EM*w)
    t <- abs(coef(summary(lmod))["EM:w", "t value"])
    df <- df.residual(lmod)
    inter <- abs(coef(summary(lmod))["EM:w", "Estimate"])
    # if there is only a single observation in any EM category, set goodness to 0
    if (min(table(EM)) <= 1) {
      goodness = c(0)
    } 
    else if (stopping.rule && inter <= parms$overall_effect/10 ) {
      goodness = c(0)
    }
    # if the estimated interaction effect is less than 1/10 of the overall effect, set goodness to 0
    # else if (t <= 500) {
    # else if (inter <= parms$overall_effect/10) {
    #   goodness = c(0)
    # }
    # otherwise, set goodness to the t-statistic for the interaction term
    else {
      goodness = c(t**2)
      #goodness = c(p)
    }
    list(goodness = goodness, direction = ux[ord])
  }
  
  ulist.used <- list(eval = etemp, split = stemp, init = itemp)
  return(ulist.used)
}


# criterion <- function(Y,w) {
#   lmod <- glm(Y ~ w + 1, family = 'poisson')
#   lmod$coefficients['w']
# }

    
