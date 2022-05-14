

################ generate synthetic data with effect modifiers ###################################################

# continuous outcome. No EM effect on intercept
synth_data_covs <- generate_syn_data_covs(sample_size = 20000, 
                                    gps_spec = 1)

synth_data_outcome <- do.call(generate_syn_data_outcome, mutate(synth_data_covs, beta = 0.5))

# # continuous outcome. EM effect on intercept
# synth_data <- generate_syn_data_het(sample_size = 20000, outcome_type = 'continuous',
#                                     gps_spec = 1, em_spec = 2, cova_spec = 1, heterogenous_intercept = TRUE, 
#                                     em_as_confounder = FALSE)
# 
# # count/binary outcome. 
# synth_data3 <- generate_syn_data_het(sample_size = 200000, outcome_type = 'binary',
#                                     gps_spec = 1, em_spec = 0, cova_spec = 1, heterogenous_intercept = FALSE, 
#                                     em_as_confounder = FALSE)
synth_data3 %>% 
  mutate(treat_level = cut(treat, breaks = 20)) %>%
  mutate(`Effect Modifiers` = interaction(em1, em2)) %>%
  mutate(`Effect Modifiers` = recode_factor(`Effect Modifiers`, 
                                            `0.0` = "E1 = 0, E2 = 0",
                                            `1.0` = "E1 = 1, E2 = 0",
                                            `0.1` = "E1 = 0, E2 = 1",
                                            `1.1` = "E1 = 1, E2 = 1")) %>%
  
  group_by(treat_level, `Effect Modifiers`, mu) %>%
  #summarise(Y = sum(Y), n = n())
  summarise(outcome_rate = sum(Y)/n()) %>%
  ggplot(aes(x = treat_level, color = `Effect Modifiers`, group = `Effect Modifiers`)) + geom_point(aes(y = outcome_rate)) + geom_line(aes(y = mu)) + 
  theme(axis.text.x = element_text(angle = 90))
ggsave(filename = 'effect_modification_binary.png')

synth_data <- synth_data2
  ################################## Optionally, plot data to visualize #####################################
# plot exposure vs each covariate
synth_data %>% 
  pivot_longer(c(cf1,cf2,cf3,cf4,cf5,cf6), names_to = 'cf', values_to = 'cov_value') %>%
  ggplot(aes(x = cov_value, y = treat)) + 
  geom_point() + 
  facet_wrap(vars(cf))

# plot outcome vs exposure by effect modifier group - continuous outcome
matched.exploration.sample %>% 
  mutate_at(vars(em1, em2), as.factor) %>%
  mutate(`Effect Modifiers` = interaction(em1, em2)) %>%
  mutate(`Effect Modifiers` = recode_factor(`Effect Modifiers`, 
                                            `0.0` = "E1 = 0, E2 = 0",
                                            `1.0` = "E1 = 1, E2 = 0",
                                            `0.1` = "E1 = 0, E2 = 1",
                                            `1.1` = "E1 = 1, E2 = 1")) %>%
  rename(exposure = treat) %>%
  rename(outcome = Y) %>%
  ggplot(aes(x = exposure, color = `Effect Modifiers`)) + geom_point(aes(y = outcome))# %>%
ggsave(filename = paste0(dir_out, 'simulation_setting1_beta0.5.png'))

synth_data <- synth_data3 %>% select(-mu)
####################### Split Data, stratified by exposure level #####################################
#Using 20 bins, Priyanka uses 100
a.vals <- seq(min(synth_data$treat), max(synth_data$treat), length.out = 20)
delta_n <- a.vals[2] - a.vals[1]

synth_data <- 
  synth_data %>%
  mutate(treat_level = cut(treat, breaks = a.vals))

# try 1/3, 1/3, 1/3 split
synth_data <- 
  split_dataset(data = synth_data, props = c(0.33, 0.33, 0.33), exposure_bins = synth_data$treat_level)

exploration.sample <- 
  synth_data %>% filter(subsample == 'exploration')

val.sample <- 
  synth_data %>% filter(subsample == 'validation')

################## GPS Matching ################################################
matched.exploration.sample <- stratified_GPS_matching(exploration.sample, delta_n)

matched.validation.sample <- stratified_GPS_matching(exploration.sample, delta_n)

# absolute correlation unmatched
data.table::setDT(synth_data)
cor_val_unmatched <- absolute_corr_fun(synth_data[,2], synth_data[, 3:(length(synth_data))])
cor_val_unmatched <- absolute_corr_fun(synth_data[,2], synth_data[, 3:(length(synth_data)-2)])

print(cor_val_unmatched)

matched.exploration.sample <- matched.exploration.sample %>%
  mutate_at(vars(em1, em2, em3, em4), as.factor)

# absolute correlation matched
data.table::setDT(matched.exploration.sample)
cor_val_matched <- absolute_corr_fun(matched.exploration.sample[,2], 
                                     matched.exploration.sample[, 6:length(matched.exploration.sample)])
print(cor_val_matched)

abs_cor = data.frame(cov = c('em1', 'em2', 'em3', 'em4', 'cf1', 'cf2', 'cf3', 'cf4', 'cf5', 'cf6'),
                     #names(synth_data)[3:length(names(synth_data))],
                     unmatched = cor_val_unmatched$absolute_corr, matched = cor_val_matched$absolute_corr) %>%
  pivot_longer(cols = c(unmatched, matched), names_to = 'dataset', values_to = 'absolute correlation')

ggplot(abs_cor, aes(x = cov, y = `absolute correlation`, color = dataset, group = dataset)) + geom_point() + geom_line()
ggsave('abs_cor.png')

################## Tree Splitting ################################################

# format data for tree splitting
matched.exploration.sample <- matched.exploration.sample %>%
  mutate(id = row_number()) %>%
  mutate_at(vars(em1, em2, em3, em4), as.factor)

# format data for tree splitting
matched.validation.sample <- matched.validation.sample %>%
  mutate(id = row_number()) %>%
  mutate_at(vars(em1, em2, em3, em4), as.factor)

overall_effect <- lm(Y ~ w + 1, data = matched.exploration.sample)$coefficients['w']

#overall_effect <- glm(Y ~ w + 1, data = matched.exploration.sample, family = 'poisson')$coefficients['w']

parms.used <- list(w  = matched.exploration.sample$w,
                   # covariates = data.used.cont.cont[, 3:dim(data.used.cont.cont)[2]],
                   Y   = matched.exploration.sample$Y,
                   overall_effect = overall_effect)   


lists <- create.sequence(matched.exploration.sample, ulist.used, parms.used)

tree.list <- lists[[1]]
g.h.list <- lists[[2]]

png(file = 'final_tree.png')
rpart.plot(tree.list[[1]])
dev.off()

library(grid)
library(ggplotify)
pdf(file = 'tree_list.pdf')
par(mfrow = c(3,1))
lapply(tree.list, function(t) {rpart.plot(t)})
dev.off()

tree.sizes <- sapply(tree.list, function(t) {nrow(t$splits)})

lambdas <- c(2,3,4,6,10,20,30,50)

complex.vals.all <- 
  lapply(lambdas, function(l) {data.frame(complex.val = evaluate.sequence(tree.list, matched.validation.sample, l), 
                                          lambda = l,
                                          pruning_step = 1:length(tree.list),
                                          tree.size = tree.sizes)})

complex.vals.all <- bind_rows(complex.vals.all) %>%
  mutate_at(vars(lambda), as.factor)

max_complexity <- complex.vals.all %>%
  group_by(lambda) %>%
  filter(complex.val == max(complex.val))

complex.vals.all %>%
  filter(tree.size != 0) %>%
  ggplot(aes(x = pruning_step, y = complex.val, color = lambda, group = lambda)) + geom_line() + geom_point(data = max_complexity)

complex.vals.all %>%
  filter(tree.size != 0) %>%
  ggplot(aes(x = tree.size, y = complex.val, color = lambda, group = lambda)) + geom_line() + geom_point(data = max_complexity) + 
  scale_x_continuous(breaks = pretty_breaks()) + 
  labs(x = 'Tree size (Number of Splits)', y = 'Complexity Value')
  ggsave('VaryingLambda.png')

png(file = 'pruned_tree.png')
rpart.plot(tree.list[[which(tree.sizes == 5)]])
dev.off()

selected.tree.size <- as.numeric(max_complexity[which(max_complexity$lambda == 50), 'tree.size'])
selected.tree <- tree.list[[which(tree.sizes == selected.tree.size)]]
rpart.plot(tree.list[[which(tree.sizes == as.numeric(selected.tree.size))]])

list(row.names(selected.tree$splits)) %in% correct_splits

correct_splits <- list(c('em2', 'em1', 'em1'), c('em1', 'em2', 'em2'))
