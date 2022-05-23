# load functions and libraries
source("predict_dpmm_fit.R")
source("dpmm_fit.R")
source("posterior_dpmm.R")
source("plot_ggpairs.R")
source("plot_alpha.R")

################################################################################
############################### Dataset 1 ######################################
################################################################################

# Showcase:
### - Fitting DPMM to continuous values
### - Standardising continuous values
### - Conditional predictions for continuous values
### - Conditional samples during model fit

# load dataset
dataset <- readRDS("datasets/dataset_test1.rds")

### standardise data
posteriors <- runModel(dataset, mcmc_iterations = 2500, L = 6, standardise = TRUE)

# Dataset 1 has two clusters of values
# We introduce missingness in covariate X1 for a specific cluster.
# This tests the capability of making conditional predictions 
#   for patients with missing variables.
rows <- 501:550
dataset_missing <- dataset
dataset_missing_predict <- dataset_missing[rows,]
dataset_missing_predict[,1] <- as.numeric(NA)

# making predictions
posteriors.dpmmfit <- predict.dpmm_fit(posteriors, dataset_missing_predict, samples = seq(1500,2500, 25))

# collect credible intervals for predictions
low_q <- NULL
median_q <- NULL
high_q <- NULL
for (i in 1:length(posteriors.dpmmfit)) {
  low_q <- c(low_q, quantile(unlist(posteriors.dpmmfit[[i]][]), probs = c(0.05)))
  median_q <- c(median_q, quantile(unlist(posteriors.dpmmfit[[i]][]), probs = c(0.5)))
  high_q <- c(high_q, quantile(unlist(posteriors.dpmmfit[[i]][]), probs = c(0.95)))
}

# Because the model used standardisation we need to
#   standardise original values using standardisation values from the model
standardised_values <- (dataset_missing[rows,1]-posteriors$mean_values[1])/posteriors$sd_values[1]

# Plot predicted values vs observed values
# Note: Dataset 1 contains two separate clusters 
#   with little explanation of structure within cluster,
#   therefore we only expect predictions within the correct cluster
cbind(standardised_values,low_q, median_q,high_q) %>%
  as.data.frame() %>%
  set_names(c("observed", "low","median","high")) %>%
  ggplot() +
  geom_point(aes(x = observed, y = median)) +
  geom_errorbar(aes(x = observed, ymin = low, ymax = high)) +
  geom_abline(aes(intercept = 0, slope = 1))

# We can draw random samples from the DPMM in order to understand
#   the distribution of values
plot.ggpairs(posteriors, iterations = seq(1500,2500,1))

# We can compare random samples from the DPMM with the original data
plot.ggpairs(posteriors, dataset, nburn = 1500)

################
##### Here is the same analysis but without standardising continuous values

posteriors <- runModel(dataset, mcmc_iterations = 2500, L = 6, standardise = FALSE)


# Dataset 1 has two clusters of values
# We introduce missingness in covariate X1 for a specific cluster.
# This tests the capability of making conditional predictions 
#   for patients with missing variables.
rows <- 1:50
dataset_missing <- dataset
dataset_missing_predict <- dataset_missing[rows,]
dataset_missing_predict[,1] <- as.numeric(NA)

# making predictions
posteriors.dpmmfit <- predict.dpmm_fit(posteriors, dataset_missing_predict, samples = seq(1500,2500, 25))

# collect credible intervals for predictions
low_q <- NULL
median_q <- NULL
high_q <- NULL
for (i in 1:length(posteriors.dpmmfit)) {
  low_q <- c(low_q, quantile(unlist(posteriors.dpmmfit[[i]][]), probs = c(0.05)))
  median_q <- c(median_q, quantile(unlist(posteriors.dpmmfit[[i]][]), probs = c(0.5)))
  high_q <- c(high_q, quantile(unlist(posteriors.dpmmfit[[i]][]), probs = c(0.95)))
}

# Plot predicted values vs observed values
# Note: Dataset 1 contains two separate clusters 
#   with little explanation of structure within cluster,
#   therefore we only expect predictions within the correct cluster
cbind(dataset_missing[1:50,1],low_q, median_q,high_q) %>%
  as.data.frame() %>%
  set_names(c("observed", "low","median","high")) %>%
  ggplot() +
  geom_point(aes(x = observed, y = median)) +
  geom_errorbar(aes(x = observed, ymin = low, ymax = high)) +
  geom_abline(aes(intercept = 0, slope = 1))

# We can draw random samples from the DPMM in order to understand
#   the distribution of values
plot.ggpairs(posteriors, iterations = seq(1500,2500,1))

# We can compare random samples from the DPMM with the original data
plot.ggpairs(posteriors, dataset, nburn = 1500)


#########################################
##### Fitting the model in the presence of missing values
#########################################


# Dataset 1 has two clusters of values
# We introduce missingness in covariate X1 for both clusters.
# This is done by choosing an individual from each cluster
#   and removing the X1 value
dataset_missing <- dataset
dataset_missing[1,1] <- NA
dataset_missing[501,1] <- NA

# Run the DPMM model with incomplete data
posteriors <- runModel(dataset_missing, mcmc_iterations = 2500, L = 6, standardise = TRUE)

# Plot predicted values vs observed values
# Note: Dataset 1 contains two separate clusters 
#   with little explanation of structure within cluster,
#   therefore we only expect predictions within the correct cluster
ggplot() +
  geom_density(aes(x = posteriors$samples$`x_cont_miss[1, 1]`[1000:2500])) +
  geom_vline(aes(xintercept = (dataset[1,1]-posteriors$mean_values[1])/posteriors$sd_values[1]))

ggplot() +
  geom_density(aes(x = posteriors$samples$`x_cont_miss[2, 1]`[1000:2500])) +
  geom_vline(aes(xintercept = (dataset[501,1]-posteriors$mean_values[1])/posteriors$sd_values[1]))

################
##### Here is the same analysis but without standardising continuous values

# Run the DPMM model with incomplete data
posteriors <- runModel(dataset_missing, mcmc_iterations = 2500, L = 6, standardise = FALSE)

# Plot predicted values vs observed values
# Note: Dataset 1 contains two separate clusters 
#   with little explanation of structure within cluster,
#   therefore we only expect predictions within the correct cluster
ggplot() +
  geom_density(aes(x = posteriors$samples$`x_cont_miss[1, 1]`[1000:2500])) +
  geom_vline(aes(xintercept = dataset[1,1]))

ggplot() +
  geom_density(aes(x = posteriors$samples$`x_cont_miss[2, 1]`[1000:2500])) +
  geom_vline(aes(xintercept = dataset[501,1]))



################################################################################
############################### Dataset 2 ######################################
################################################################################

# Showcase:
### - Fitting DPMM to continuous + categorical values
### - Conditional predictions for continuous values from categorical values
### - Conditional predictions for categorical values from continuous values
### - Conditional samples during model fit from categorical values

# load dataset
dataset <- readRDS("datasets/dataset_test2.rds")

posteriors <- runModel(dataset, mcmc_iterations = 2500, L = 6, standardise = FALSE)

# Dataset 2 has two clusters of values
# We introduce missingness in covariate X2 for a specific cluster.
# This tests the capability of making conditional predictions 
#   for patients with missing variables where only X3 (categorical variable)
#   is able to tell which cluster the values belongs to
rows <- 950:1000
dataset_missing <- posteriors$dataset
dataset_missing <- dataset_missing[rows,]
dataset_missing[,2] <- as.numeric(NA)

# making predictions
posteriors.dpmmfit <- predict.dpmm_fit(posteriors, dataset_missing, samples = seq(1500,2500, 25))

# collect credible intervals for predictions
low_q <- NULL
median_q <- NULL
high_q <- NULL
for (i in 1:length(posteriors.dpmmfit)) {
  low_q <- c(low_q, quantile(unlist(posteriors.dpmmfit[[i]][]), probs = c(0.05)))
  median_q <- c(median_q, quantile(unlist(posteriors.dpmmfit[[i]][]), probs = c(0.5)))
  high_q <- c(high_q, quantile(unlist(posteriors.dpmmfit[[i]][]), probs = c(0.95)))
}

# Plot predicted values vs observed values
# Note: Dataset 2 contains two separate clusters 
#   with little explanation of structure within cluster,
#   therefore we only expect predictions within the correct cluster
cbind(posteriors$dataset[rows,2],low_q, median_q,high_q) %>%
  as.data.frame() %>%
  set_names(c("observed", "low","median","high")) %>%
  ggplot() +
  geom_point(aes(x = observed, y = median)) +
  geom_errorbar(aes(x = observed, ymin = low, ymax = high)) +
  geom_abline(aes(intercept = 0, slope = 1))

##########
##########

# Dataset 2 has two clusters of values
# We introduce missingness in covariate X3 for a specific cluster.
# This tests the capability of making conditional predictions 
#   for patients with missing variables where only X2 (continuous variable)
#   is able to tell which cluster the values belongs to
rows <- 950:1000
dataset_missing <- posteriors$dataset
dataset_missing <- dataset_missing[rows,]
dataset_missing[,3] <- factor(NA, levels = levels(posteriors$dataset[,3]))

# making predictions
posteriors.dpmmfit <- predict.dpmm_fit(posteriors, dataset_missing, samples = seq(1500,2500, 25))

# collect credible intervals for predictions
predicted <- NULL
for (i in 1:length(posteriors.dpmmfit)) {
  var <- table(unlist(posteriors.dpmmfit[i]))
  predicted <- c(predicted, names(var[var==max(var)]))
}
predicted <- as.numeric(predicted)

# Table predicted values vs observed values
# Note: Dataset 2 contains two separate clusters 
#   with the continuous values (X1, X2) perfectly explaining
#   which cluster the individual is part of
cbind(posteriors$dataset[rows,3],predicted) %>%
  as.data.frame() %>%
  set_names(c("observed", "predicted")) %>%
  mutate(observed = factor(observed),
         predicted = factor(predicted)) %>%
  table()

# We can draw random samples from the DPMM in order to understand
#   the distribution of values
plot.ggpairs(posteriors, iterations = seq(1500,2500,1))

# We can compare random samples from the DPMM with the original data
plot.ggpairs(posteriors, dataset, nburn = 1500)


#########################################
##### Fitting the model in the presence of missing values
#########################################


# Dataset 2 has two clusters of values
# We introduce missingness in covariate X3 for both clusters.
# This is done by choosing an individual from each cluster
#   and removing the X3 value
dataset_missing <- dataset
dataset_missing[2,3] <- NA
dataset_missing[501,3] <- NA

# Run the DPMM model with incomplete data
posteriors <- runModel(dataset_missing, mcmc_iterations = 2500, L = 6, standardise = FALSE)

# Table predicted values vs observed values
# Note: Dataset 2 contains two separate clusters
#   with the continuous values (X1, X2) perfectly explaining
#   which cluster the individual is part of
cbind(posteriors$samples$`x_disc_miss[1, 1]`[1000:2500], rep(as.numeric(dataset[2,3]), 1500)) %>%
  as.data.frame() %>%
  set_names(c("Predicted","True")) %>%
  table()

cbind(posteriors$samples$`x_disc_miss[2, 1]`[1000:2500], rep(as.numeric(dataset[501,3]), 1500)) %>%
  as.data.frame() %>%
  set_names(c("Predicted","True")) %>%
  table()

##########
##########

# Dataset 2 has two clusters of values
# We introduce missingness in covariate X2 for both clusters.
# This is done by choosing an individual from each cluster
#   and removing the X2 value
dataset_missing <- dataset
dataset_missing[1,2] <- NA
dataset_missing[501,2] <- NA

# Run the DPMM model with incomplete data
posteriors <- runModel(dataset_missing, mcmc_iterations = 2500, L = 6, standardise = FALSE)

# Plot predicted values vs observed values
# Note: Dataset 2 contains two separate clusters 
#   with little explanation of structure within cluster,
#   therefore we only expect predictions within the correct cluster
ggplot() +
  geom_density(aes(x = posteriors$samples$`x_cont_miss[1, 2]`[1000:2500])) +
  geom_vline(aes(xintercept = dataset[1,2]))

ggplot() +
  geom_density(aes(x = posteriors$samples$`x_cont_miss[2, 2]`[1000:2500])) +
  geom_vline(aes(xintercept = dataset[501,2]))

##########
##########

# Dataset 2 has two clusters of values
# We introduce missingness in all continuous variables (X1, X2) for both clusters.
# This is done by choosing an individual from each cluster
#   and removing the X1 and X2 values
dataset_missing <- dataset
dataset_missing[c(1,501),c(1,2)] <- NA

# Run the DPMM model with incomplete data
posteriors <- runModel(dataset_missing, mcmc_iterations = 2500, L = 6, standardise = FALSE)

# Plot predicted values vs observed values
# Note: Dataset 2 contains two separate clusters
#   with X3 providing enough information regarding the correct cluster
#   but not enough information for the structure within the cluster
ggplot() +
  geom_point(aes(x = posteriors$samples$`x_cont_miss[1, 1]`[1000:2500], y = posteriors$samples$`x_cont_miss[1, 2]`[1000:2500]), colour = "darkorange") +
  geom_point(aes(x = dataset[1,1], y = dataset[1,2]), colour = "darkorange3", fill = "black", size = 3) +
  geom_point(aes(x = posteriors$samples$`x_cont_miss[2, 1]`[1000:2500], y = posteriors$samples$`x_cont_miss[2, 2]`[1000:2500]), colour = "royalblue") +
  geom_point(aes(x = dataset[501,1], y = dataset[501,2]), colour = "royalblue4", size = 3)



################################################################################
############################### Dataset 3 ######################################
################################################################################

# Showcase:
### - Fitting DPMM to categorical values
### - Conditional predictions for categorical values
### - Conditional samples during model fit from categorical values

# load dataset
dataset <- readRDS("datasets/dataset_test3.rds")

posteriors <- runModel(dataset, mcmc_iterations = 2500, L = 6, standardise = FALSE)

# Dataset 3 has two clusters of values
# We introduce missingness in covariate X1 for a specific cluster.
# There is clear separation between clusters therefore 
#   X2 explains X1 well
rows <- 950:1000
dataset_missing <- dataset
dataset_missing <- dataset_missing[rows,]
dataset_missing[,1] <- factor(NA, levels = levels(posteriors$dataset[,1]))

# making predictions
posteriors.dpmmfit <- predict.dpmm_fit(posteriors, dataset_missing, samples = seq(1500,2500, 25))

# collect credible intervals for predictions
predicted <- NULL
for (i in 1:length(posteriors.dpmmfit)) {
  var <- table(unlist(posteriors.dpmmfit[i]))
  predicted <- c(predicted, names(var[var==max(var)]))
}
predicted <- as.numeric(predicted)

# Table predicted values vs observed values
# Note: Dataset 2 contains two separate clusters
#   with X2 perfectly explaining X1
cbind(dataset[rows,1],predicted) %>%
  as.data.frame() %>%
  set_names(c("observed", "predicted")) %>%
  mutate(observed = factor(observed),
         predicted = factor(predicted)) %>%
  table()

# We can draw random samples from the DPMM in order to understand
#   the distribution of values
plot.ggpairs(posteriors, iterations = seq(1500,2500,1))

# We can compare random samples from the DPMM with the original data
plot.ggpairs(posteriors, dataset, nburn = 1500)


#########################################
##### Fitting the model in the presence of missing values
#########################################


# Dataset 3 has two clusters of values
# We introduce missingness in covariate X1 for both clusters.
# This is done by choosing an individual from each cluster
#   and removing the X3 value
dataset_missing <- dataset
dataset_missing[1,1] <- factor(NA, levels = levels(dataset[,1]))
dataset_missing[501,1] <- factor(NA, levels = levels(dataset[,1]))

# Run the DPMM model with incomplete data
posteriors <- runModel(dataset_missing, mcmc_iterations = 2500, L = 6, standardise = FALSE)

# Table predicted values vs observed values
# Note: Dataset 3 contains two separate clusters
#   with X2 perfectly explaining X1
cbind(posteriors$samples$`x_disc_miss[1, 1]`[1000:2500], dataset[1,1]) %>%
  as.data.frame() %>%
  set_names(c("Prediction","True Value")) %>%
  table()

cbind(posteriors$samples$`x_disc_miss[2, 1]`[1000:2500], dataset[501,1]) %>%
  as.data.frame() %>%
  set_names(c("Prediction","True Value")) %>%
  table()



################################################################################
############################### Dataset 4 ######################################
################################################################################

# Showcase:
### - Fitting DPMM to continuous + categorical values
### - Conditional predictions for continuous values
### - Conditional samples during model fit from continuous + categorical values

# load dataset
dataset <- readRDS("datasets/dataset_test4.rds")

posteriors <- runModel(dataset, mcmc_iterations = 2500, L = 6, standardise = FALSE)

# Dataset 4 has a specific structure
# We introduce missingness in covariate X1
#   and use all other variables to explain the missing variable.
rows <- sample(1:nrow(posteriors$dataset), 100)
dataset_missing <- posteriors$dataset
dataset_missing <- dataset_missing[rows,]
dataset_missing[,1] <- as.numeric(NA)

# making predictions
posteriors.dpmmfit <- predict.dpmm_fit(posteriors, dataset_missing, samples = seq(1500,2500, 25))

# collect credible intervals for predictions
low_q <- NULL
median_q <- NULL
high_q <- NULL
for (i in 1:length(posteriors.dpmmfit)) {
  low_q <- c(low_q, quantile(unlist(posteriors.dpmmfit[[i]][]), probs = c(0.05)))
  median_q <- c(median_q, quantile(unlist(posteriors.dpmmfit[[i]][]), probs = c(0.5)))
  high_q <- c(high_q, quantile(unlist(posteriors.dpmmfit[[i]][]), probs = c(0.95)))
}

# Plot predicted values vs observed values
# Note: Dataset 4 has a strong specific structure
#   and we use all other variables to explaing the missingness
cbind(posteriors$dataset[rows,1],low_q, median_q,high_q) %>%
  as.data.frame() %>%
  set_names(c("observed", "low","median","high")) %>%
  ggplot() +
  geom_point(aes(x = observed, y = median)) +
  geom_errorbar(aes(x = observed, ymin = low, ymax = high)) +
  geom_abline(aes(intercept = 0, slope = 1))

# We can draw random samples from the DPMM in order to understand
#   the distribution of values
plot.ggpairs(posteriors, iterations = seq(1500,2500,1))

# We can compare random samples from the DPMM with the original data
plot.ggpairs(posteriors, dataset, nburn = 1500)



################################################################################
############################### Dataset 5 ######################################
################################################################################

# Showcase:
### - Fitting DPMM to continuous + categorical values
### - Conditional predictions for continuous values
### - Conditional samples during model fit from continuous + categorical values

# load dataset
dataset <- readRDS("datasets/dataset_test5.rds")

posteriors <- runModel(dataset, mcmc_iterations = 2500, L = 6, standardise = FALSE)

# We can draw random samples from the DPMM in order to understand
#   the distribution of values
plot.ggpairs(posteriors, iterations = seq(1500,2500,1))

# We can compare random samples from the DPMM with the original data
plot.ggpairs(posteriors, dataset, nburn = 1500)


#########################################
##### Fitting the model in the presence of missing values
#########################################


# Dataset 5 has two clusters of values
# We introduce missingness in covariate X1 to explore 
#   the prediction of continuous variables and
#   we introduce missingness in the covariate X3 where X4
#   is able to explain the missing variable
dataset_missing <- dataset
dataset_missing[1,1] <- NA
dataset_missing[2,3] <- NA

# Run the DPMM model with incomplete data
posteriors <- runModel(dataset_missing, mcmc_iterations = 2500, L = 6, standardise = FALSE)

# Plot predicted values vs observed values
ggplot() +
  geom_density(aes(x = posteriors$samples$`x_cont_miss[1, 1]`[1000:2500])) +
  geom_vline(aes(xintercept = dataset[1,1]))

# Table predicted values vs observed values
# Note: Variable X3 is well explained by X4
cbind(posteriors$samples$`x_disc_miss[2, 1]`[1000:2500], rep(dataset[2,3], rep = 1500)) %>%
  as.data.frame() %>%
  set_names(c("Predicted","True")) %>%
  table()

##########
##########


# Dataset 5 has two clusters of values
# We introduce missingness in covariates X1,X2,X3 to explore 
#   the prediction of continuous + categorical variables and
#   these are conditional on X4 which splits the data 
#   into a specific cluster
dataset_missing <- dataset
dataset_missing[1,1:3] <- NA

# Run the DPMM model with incomplete data
posteriors <- runModel(dataset_missing, mcmc_iterations = 2500, L = 6, standardise = FALSE)

# Plot predicted values vs observed values
# Note: we expect only one cluster represented
ggplot() +
  geom_point(aes(x = posteriors$samples$`x_cont_miss[1, 1]`[1000:2500], y = posteriors$samples$`x_cont_miss[1, 2]`[1000:2500]))

# Table predicted values vs observed values
# Note: Variable X3 is well explained by X4
cbind(posteriors$samples$`x_disc_miss[1, 1]`[1000:2500], rep(dataset[1,3], rep = 1500)) %>%
  as.data.frame() %>%
  set_names(c("Predicted","True")) %>%
  table()

### if both categorical are gone and X1, with X2 being one of the clusters. This should give the right X3 = 1, X4 = 1 or 2 


##########
##########


# Dataset 5 has two clusters of values
# We introduce missingness in covariates X1,X3,X4 to explore 
#   the prediction of continuous + categorical variables,
#   with X2 representing one of the clusters, X3 should only 
#   have one category predicted, X4 should only have two categories predicted
dataset_missing <- dataset
dataset_missing[1,c(1,3,4)] <- NA

# Run the DPMM model with incomplete data
posteriors <- runModel(dataset_missing, mcmc_iterations = 2500, L = 6, standardise = FALSE)

# Table predicted values vs observed values
# Note: Variable X3 is well explained by X2
cbind(posteriors$samples$`x_disc_miss[1, 1]`[1000:2500], rep(dataset[1,3], rep = 1500)) %>%
  as.data.frame() %>%
  set_names(c("Predicted","True")) %>%
  table()

# Table predicted values vs observed values
# Note: Variable X4 is well explained by X2 but should have two categories
cbind(posteriors$samples$`x_disc_miss[1, 2]`[1000:2500], rep(dataset[1,4], rep = 1500)) %>%
  as.data.frame() %>%
  set_names(c("Predicted","True")) %>%
  table()

# Plot predicted values vs observed values
# Note: we expect only one cluster represented
ggplot() +
  geom_density(aes(x = posteriors$samples$`x_cont_miss[1, 1]`[1000:2500])) +
  geom_vline(aes(xintercept = dataset[1,1]))



################################################################################
############################### Dataset 6 ######################################
################################################################################

# Showcase:
### - Fitting DPMM to continuous values with high correlation
### - Conditional predictions for continuous values
### - Conditional samples during model fit from continuous values

# load dataset
dataset <- readRDS("datasets/dataset_test6.rds")

posteriors <- runModel(dataset, mcmc_iterations = 2500, L = 12, standardise = FALSE)

# Dataset 6 has specific clusters of data
# We introduce missingness in covariate X1 and
#   X2 explains the variable well
dataset_missing <- dataset
dataset_missing_predict <- dataset_missing[1:50,]
dataset_missing_predict[,1] <- as.numeric(NA)
posteriors.dpmmfit <- predict.dpmm_fit(posteriors, dataset_missing_predict, samples = seq(1500,2500, 25))

# collect credible intervals for predictions
low_q <- NULL
median_q <- NULL
high_q <- NULL
for (i in 1:length(posteriors.dpmmfit)) {
  low_q <- c(low_q, quantile(unlist(posteriors.dpmmfit[[i]][]), probs = c(0.05)))
  median_q <- c(median_q, quantile(unlist(posteriors.dpmmfit[[i]][]), probs = c(0.5)))
  high_q <- c(high_q, quantile(unlist(posteriors.dpmmfit[[i]][]), probs = c(0.95)))
}

# Plot predicted values vs observed values
# Note: Dataset 6 has specific clusters of data
#   and X1 is well explained by X2
cbind(dataset_missing[1:50,1],low_q, median_q,high_q) %>%
  as.data.frame() %>%
  set_names(c("observed", "low","median","high")) %>%
  ggplot() +
  geom_point(aes(x = observed, y = median)) +
  geom_errorbar(aes(x = observed, ymin = low, ymax = high)) +
  geom_abline(aes(intercept = 0, slope = 1))

# We can draw random samples from the DPMM in order to understand
#   the distribution of values
plot.ggpairs(posteriors, iterations = seq(1500,2500,1))

# We can compare random samples from the DPMM with the original data
plot.ggpairs(posteriors, dataset, nburn = 1500)


#########################################
##### Fitting the model in the presence of missing values
#########################################


# Dataset 6 has two clusters of values
# We introduce missingness in covariates X1 to explore 
#   the prediction of continuous values from two
#   different clustersdataset_missing <- dataset
dataset_missing[1,1] <- NA
dataset_missing[501,1] <- NA

# Run the DPMM model with incomplete data
posteriors <- runModel(dataset_missing, mcmc_iterations = 2500, L = 6, standardise = FALSE)

# Plot predicted values vs observed values
# Note: Dataset 6 has specific clusters of data
#   and X1 is well explained by X2
ggplot() +
  geom_density(aes(x = posteriors$samples$`x_cont_miss[1, 1]`[1000:2500])) +
  geom_vline(aes(xintercept = dataset[1,1]))

ggplot() +
  geom_density(aes(x = posteriors$samples$`x_cont_miss[2, 1]`[1000:2500])) +
  geom_vline(aes(xintercept = dataset[501,1]))



################################################################################
############## Demonstrating uses of alpha function ############################
################################################################################

# load dataset
dataset <- readRDS("datasets/dataset_test4.rds")

# Run DPMM model
posteriors.1 <- runModel(dataset, mcmc_iterations = 2500,  L = 12, standardise = FALSE)

# Plot of: 
#   - number of components used
#   - average number of individuals in ranked components
#   - trace plot for alpha values
plot.alpha(posteriors.1, nburn = 1500, thinning = 1000)

# Run multiple DPMM models
posteriors.2 <- runModel(dataset, mcmc_iterations = 2500,  L = 12, standardise = FALSE)

posteriors.3 <- runModel(dataset, mcmc_iterations = 2500,  L = 12, standardise = FALSE)

# Combine all DPMM models into a single list
posteriors <- list(posteriors.1, posteriors.2, posteriors.3)

# Plots values
plot.alpha(posteriors, nburn = 1500, thinning = 1000)

# We can draw random samples from the DPMM in order to understand
#   the distribution of values. When using a list of Model posteriors
#   'iterations' are the samples used from each Model posteriors
plot.ggpairs(posteriors, iterations = seq(2000,2500,1))

# We can compare random samples from the DPMM with the original data
#   When using a list of Model posteriors, 'nburn' corresponds to
#   the initial samples removed from each Model posteriors
plot.ggpairs(posteriors, dataset, nburn = 2000)



