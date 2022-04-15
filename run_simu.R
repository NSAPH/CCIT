
correct_splits <- list()
noise.var <- c('em1', 'em2', 'em3', 'em4')
true_trt_effect_func <- function(df) {
  df %>%
    mutate(eff = 10)
}


######################## Setting 1: 4 groups, same intercept
correct_splits <- list(c('em2', 'em1', 'em1'), c('em1', 'em2', 'em2'))
noise.var <- c('em3', 'em4')
beta <- 0.15

true_trt_effect_func <- function(df) {
  df %>%
    mutate_at(vars(em1, em2, em3, em4), function(e) {as.numeric(as.character(e))}) %>%
    mutate(eff = 10*(1 + beta * (0.5 * em1 - 0.8 * em2))) %>%
    mutate_at(vars(em1, em2, em3, em4), as.factor)
  
}

lambdas <- c(1.5, 1.75, 2, 3, 4)
lambdas <- c(3000, 4000, 5000)

selected.tree.data <- data.frame()
for (i in 1:50) {
  
  ret <- evaluate_simulated_continuous_outcome(lambda = lambdas, correct_splits = correct_splits, noise.var = noise.var, 
                                               num_exposure_cats = 20, em_spec = 1, stopping.rule = TRUE)
  selected.tree.data <- rbind(selected.tree.data, as.data.frame(ret))
}

all(selected.tree.data$selected.correct.splits )

save(selected.tree.data, file = paste0(dir_out,'selected.tree.data.with.stopping.setting1.lambda3000-5000.50trials.beta0.15.RData'))

load(paste0(dir_out,'selected.tree.data.with.stopping.setting1.beta0.3.lambda3000-5000.50trials.RData'))

selected.tree.data %>%
  ggplot(aes(x = lambda, y = bias, group = lambda)) + geom_boxplot() + 
  scale_x_continuous(breaks = lambdas)
ggsave(paste0(dir_out, 'MSE.no.stopping100trials.png'))


################## Setting 2: 3 groups, same intercept
correct_splits <- list(c('em2', 'em1'))
noise.var <- c('em3', 'em4')

true_trt_effect_func <- function(df) {
  df %>%
    mutate_at(vars(em1, em2, em3, em4), function(e) {as.numeric(as.character(e))}) %>%
    mutate(eff = 10*(1 + 0.5 * em1 * em2 + 0.2 * em2)) %>%
    mutate_at(vars(em1, em2, em3, em4), as.factor)
  
}

#lambdas <- c(1.5, 1.75, 2, 3, 4)

lambdas <- c(3000, 4000, 5000)

selected.tree.data <- data.frame()
for (i in 1:50) {
  
  ret <- evaluate_simulated_continuous_outcome(lambda = lambdas, correct_splits = correct_splits, noise.var = noise.var, 
                                               num_exposure_cats = 20, em_spec = 2)
  selected.tree.data <- rbind(selected.tree.data, as.data.frame(ret))
}

save(selected.tree.data, file = paste0(dir_out,'selected.tree.data.with.stopping.setting2.sd10.lambda3000-5000.50trials.RData'))

load(paste0(dir_out,'no.stopping.RData'))


selected.tree.data %>%
  ggplot(aes(x = lambda, y = mse, group = lambda)) + geom_boxplot() + 
  scale_x_continuous(breaks = lambdas)
ggsave(paste0(dir_out, 'MSE.no.stopping100trials.png'))



######################## Setting 3: 4 groups, different intercept
correct_splits <- list(c('em2', 'em1', 'em1'), c('em1', 'em2', 'em2'))
noise.var <- c('em3', 'em4')

true_trt_effect_func <- function(df) {
  df %>%
    mutate_at(vars(em1, em2, em3, em4), function(e) {as.numeric(as.character(e))}) %>%
    mutate(eff = 10*(1 + 0.5 * em1 - 0.8 * em2)) %>%
    mutate_at(vars(em1, em2, em3, em4), as.factor)
  
}

#lambdas <- c(1.5, 1.75, 2, 3, 4)

lambdas <- c(3000, 4000, 5000)

selected.tree.data <- data.frame()
for (i in 1:50) {
  
  ret <- evaluate_simulated_continuous_outcome(lambda = lambdas, correct_splits = correct_splits, noise.var = noise.var, 
                                               num_exposure_cats = 20, heterogenous_intercept = TRUE)
  selected.tree.data <- rbind(selected.tree.data, as.data.frame(ret))
}

save(selected.tree.data, file = paste0(dir_out,'selected.tree.data.with.stopping.setting3.lambda3000-5000.50trials.RData'))

load(paste0(dir_out,'no.stopping.RData'))

selected.tree.data %>%
  ggplot(aes(x = lambda, y = mse, group = lambda)) + geom_boxplot() + 
  scale_x_continuous(breaks = lambdas)
ggsave(paste0(dir_out, 'MSE.with.stopping50trials.png'))


selected.tree.data %>%
  ggplot(aes(x = lambda, y = bias, group = lambda)) + geom_boxplot() + 
  scale_x_continuous(breaks = lambdas)

selected.tree.data %>% summarise_at(vars(numb.noise, selected.correct.splits, bias, mse), mean)

################## Setting 4: 3 groups, different intercept
correct_splits <- list(c('em2', 'em1'))
noise.var <- c('em3', 'em4')

true_trt_effect_func <- function(df) {
  df %>%
    mutate_at(vars(em1, em2, em3, em4), function(e) {as.numeric(as.character(e))}) %>%
    mutate(eff = 10*(1 + 0.5 * em1 * em2 + 0.2 * em2)) %>%
    mutate_at(vars(em1, em2, em3, em4), as.factor)
  
}

#lambdas <- c(1.5, 1.75, 2, 3, 4)

lambdas <- c(3000, 4000, 5000)

selected.tree.data <- data.frame()
for (i in 1:50) {
  
  ret <- evaluate_simulated_continuous_outcome(lambda = lambdas, correct_splits = correct_splits, noise.var = noise.var, 
                                               num_exposure_cats = 20, em_spec = 2, heterogenous_intercept = TRUE)
  selected.tree.data <- rbind(selected.tree.data, as.data.frame(ret))
}

save(selected.tree.data, file = paste0(dir_out,'selected.tree.data.with.stopping.setting4.lambda3000-5000.50trials.RData'))

## results for setting 1 with different beta
selected.tree.data.combined <- data.frame()
for (beta in c('0.1', '0.15', '0.2', '0.3', '1')){
  load(file = paste0(dir_out,'selected.tree.data.with.stopping.setting', 1, '.lambda3000-5000.50trials.beta', beta, '.RData'))
  selected.tree.data.combined <- selected.tree.data.combined %>% rbind(selected.tree.data %>% cbind(beta = beta))
}

selected.tree.data.combined %>%
  filter(lambda == 3000) %>%
  ggplot(aes(y = mse, group = lambda)) + geom_boxplot() + 
  facet_wrap(vars(beta), nrow = 1) + 
  scale_x_continuous(breaks = lambdas)
ggsave(file = paste0(dir_out, 'mse.with.stopping.beta.png'))

selected.tree.data.combined %>%
  filter(lambda == 3000) %>%
  ggplot(aes(y = bias, group = lambda)) + geom_boxplot() + 
  facet_wrap(vars(beta), nrow = 1) + 
  scale_x_continuous(breaks = lambdas)
ggsave(file = paste0(dir_out, 'bias.with.stopping.beta.png'))


### Combining results for different settings
selected.tree.data.combined <- data.frame()
for (i in 1:5){
  load(file = paste0(dir_out,'selected.tree.data.with.stopping.setting', i, '.lambda3000-5000.50trials.RData'))
  selected.tree.data.combined <- selected.tree.data.combined %>% rbind(selected.tree.data %>% cbind(setting = i))
}

setting.labs <- c('Setting 1', 'Setting 2', 'Setting 3', 'Setting 4', 'Setting 5')
names(setting.labs) <- c(1,2,3,4,5)

selected.tree.data.combined %>%
  filter(lambda == 3000) %>%
  ggplot(aes(y = mse, group = lambda)) + geom_boxplot() + 
  facet_wrap(vars(setting), nrow = 1, labeller = labeller(setting = setting.labs)) + 
  scale_x_continuous(breaks = lambdas)
ggsave(file = paste0(dir_out, 'mse.with.stopping2.png'))

selected.tree.data.combined %>%
  filter(lambda == 3000) %>%
  ggplot(aes(y = bias, group = lambda)) + geom_boxplot() + 
  facet_wrap(vars(setting), nrow = 1, labeller = labeller(setting = setting.labs)) + 
  scale_x_continuous(breaks = lambdas)
ggsave(file = paste0(dir_out, 'bias.with.stopping2.png'))

selected.tree.data.combined %>% group_by(setting) %>%
  summarise_at(vars(numb.noise, selected.correct.splits, bias, mse), mean)







load(paste0(dir_out,'no.stopping.RData'))


selected.tree.data %>%
  ggplot(aes(x = lambda, y = mse, group = lambda)) + geom_boxplot() + 
  scale_x_continuous(breaks = lambdas)
ggsave(paste0(dir_out, 'MSE.no.stopping100trials.png'))







selected.tree.data %>%
  ggplot(aes(x = lambda, fill = selected.correct.splits, y = ..count../100)) + geom_histogram(position = 'dodge') + 
  scale_x_continuous(breaks = lambdas) + 
  labs(y = 'Proportion of Simulation trials') 
ggsave(paste0(dir_out, 'Prop.correct.trees.no.stopping.png'))

selected.tree.data %>%
  ggplot(aes(x = lambda, y = numb.noise, group = lambda)) + geom_boxplot() + 
  scale_x_continuous(breaks = lambdas) + 
  labs(y = 'Proportion of Simulation trials')
ggsave(paste0(dir_out, 'Prop.correct.trees.no.stopping.png'))

selected.tree.data %>%
  mutate('Pruned to Correct Size' = (floor(tree.size) == 3)) %>%
  ggplot(aes(x = lambda, fill = `Pruned to Correct Size`, y = ..count../100)) + geom_histogram(position = 'dodge') + 
  scale_x_continuous(breaks = lambdas) + 
  labs(y = 'Proportion of Simulation trials')
ggsave(paste0(dir_out, 'Pruning_accuracy.no.stopping.png'))
ggsave(paste0(dir_out, 'Pruning_accuracy_bootstrap_n100000.png'))


selected.tree.data %>%
  ggplot(aes(x = lambda, y = floor(tree.size), group = lambda)) + geom_boxplot() + 
  geom_hline(yintercept = 3) + 
  scale_x_continuous(breaks = lambdas)
ggsave(paste0(dir_out, 'Pruned_tree_size_crossfit.png'))
ggsave(paste0(dir_out, 'Pruned_tree_size_bootstrapn100000.png'))



selected.tree.data %>%
  ggplot(aes(x = lambda, y = selected.tree.size, group = lambda)) + geom_boxplot() + 
  scale_x_continuous(breaks = lambdas)
ggsave(paste0(dir_out, 'Pruned_tree_size.png'))


selected.tree.data %>%
  mutate('Pruned to Correct Size' = as.logical(selected.correct.splits)) %>%
  ggplot(aes(x = lambda, fill = `Pruned to Correct Size`)) + geom_histogram(position = 'dodge') + 
  scale_x_continuous(breaks = lambdas)
ggsave(paste0(dir_out, 'Pruning_accuracy.png'))


selected.tree.data %>%
  group_by(lambda) %>%
  summarise(prop_correct = sum(selected.correct.splits)/n(), 
            avg_tree_size = mean(selected.tree.size), 
            prop_correct_in_sequence = sum(correct.tree.in.sequence)/n(),
            prop_correct_size = sum(selected.tree.size == 3)/n())

ret <- evaluate_simulated_continuous_outcome(lambda = 1, correct_splits = correct_splits, num_exposure_cats = 20)
