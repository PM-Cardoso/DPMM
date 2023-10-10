## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  fig.dim=c(6,6),
  comment = "#>"
)

## ----message = FALSE, warning = FALSE-----------------------------------------
library(DPMM)
library(tidyverse)

## ----message = FALSE, warning = FALSE-----------------------------------------
data(dataset_1)

GGally::ggpairs(dataset_1)

## ----message = FALSE, warning = FALSE-----------------------------------------
posteriors <- runModel(dataset_1, 
                       mcmc_iterations = 2500, 
                       L = 6, 
                       mcmc_chains = 2, 
                       standardise = TRUE)

## ----message = FALSE, warning = FALSE-----------------------------------------
plot_ggpairs(posteriors, 
             newdata = dataset_1, 
             nburn = 1500)

## ----message = FALSE, warning = FALSE-----------------------------------------
rows <- 1:50
dataset_missing <- dataset_1
dataset_missing_predict <- dataset_missing[rows,]
dataset_missing_predict[,"X1"] <- as.numeric(NA)

## ----message = FALSE, warning = FALSE-----------------------------------------
posteriors.prediction <- predict_dpmm_fit(posteriors, 
                                          newdata = dataset_missing_predict,
                                          samples = seq(1500,2500, 25))

## ----message = FALSE, warning = FALSE-----------------------------------------
low_q <- NULL
median_q <- NULL
high_q <- NULL
for (i in 1:length(posteriors.prediction)) {
  low_q <- c(low_q, quantile(unlist(posteriors.prediction[[i]][]), probs = c(0.05)))
  median_q <- c(median_q, quantile(unlist(posteriors.prediction[[i]][]), probs = c(0.5)))
  high_q <- c(high_q, quantile(unlist(posteriors.prediction[[i]][]), probs = c(0.95)))
}

## ----message = FALSE, warning = FALSE-----------------------------------------
standardised_values <- (dataset_missing[rows,"X1"]-posteriors$mean_values[1])/posteriors$sd_values[1]

## ----message = FALSE, warning = FALSE-----------------------------------------
cbind(standardised_values,low_q, median_q,high_q) %>%
  as.data.frame() %>%
   rlang::set_names(c("observed", "low","median","high")) %>%
  ggplot() +
  geom_point(aes(x = observed, y = median)) +
  geom_errorbar(aes(x = observed, ymin = low, ymax = high)) +
  geom_abline(aes(intercept = 0, slope = 1))

## ----message = FALSE, warning = FALSE-----------------------------------------
posteriors <- runModel(dataset_1, 
                       mcmc_iterations = 2500, 
                       L = 6, 
                       mcmc_chains = 2, 
                       standardise = FALSE)

## ----message = FALSE, warning = FALSE-----------------------------------------
plot_ggpairs(posteriors, dataset_1, nburn = 1500)

## ----message = FALSE, warning = FALSE-----------------------------------------
posteriors.prediction <- predict_dpmm_fit(posteriors, 
                                          newdata = dataset_missing_predict,
                                          samples = seq(1500,2500, 25))

## ----message = FALSE, warning = FALSE-----------------------------------------
low_q <- NULL
median_q <- NULL
high_q <- NULL
for (i in 1:length(posteriors.prediction)) {
  low_q <- c(low_q, quantile(unlist(posteriors.prediction[[i]][]), probs = c(0.05)))
  median_q <- c(median_q, quantile(unlist(posteriors.prediction[[i]][]), probs = c(0.5)))
  high_q <- c(high_q, quantile(unlist(posteriors.prediction[[i]][]), probs = c(0.95)))
}

cbind(dataset_missing[1:50,"X1"],low_q, median_q,high_q) %>%
  as.data.frame() %>%
   rlang::set_names(c("observed", "low","median","high")) %>%
  ggplot() +
  geom_point(aes(x = observed, y = median)) +
  geom_errorbar(aes(x = observed, ymin = low, ymax = high)) +
  geom_abline(aes(intercept = 0, slope = 1))

## ----message = FALSE, warning = FALSE-----------------------------------------
dataset_missing <- dataset_1
dataset_missing[1,"X1"] <- NA
dataset_missing[501,"X1"] <- NA

## ----message = FALSE, warning = FALSE-----------------------------------------
posteriors <- runModel(dataset_missing, 
                       mcmc_iterations = 2500, 
                       L = 6, 
                       mcmc_chains = 2, 
                       standardise = TRUE)

## ----message = FALSE, warning = FALSE-----------------------------------------
intercept_1 = (dataset_1[1,"X1"]-posteriors$mean_values[1])/posteriors$sd_values[1]
ggplot() +
  geom_density(aes(x = as_tibble(posteriors$samples$chain1)%>%
                     select(`x_cont_miss[1, 1]`)%>%
                     dplyr::slice(1000:2500)%>%
                     unlist())) +
  geom_vline(aes(xintercept = intercept_1))

intercept_2 = (dataset_1[501,"X1"]-posteriors$mean_values[1])/posteriors$sd_values[1]
ggplot() +
  geom_density(aes(x = as_tibble(posteriors$samples$chain1)%>%
                     select(`x_cont_miss[2, 1]`)%>%
                     dplyr::slice(1000:2500)%>%
                     unlist())) +
  geom_vline(aes(xintercept = intercept_2))


## ----message = FALSE, warning = FALSE-----------------------------------------
posteriors <- runModel(dataset_missing, 
                       mcmc_iterations = 2500, 
                       L = 6, 
                       mcmc_chains = 2, 
                       standardise = FALSE)

## ----message = FALSE, warning = FALSE-----------------------------------------
ggplot() +
  geom_density(aes(x = as_tibble(posteriors$samples$chain1)%>%
                     select(`x_cont_miss[1, 1]`)%>%
                     dplyr::slice(1000:2500)%>%
                     unlist())) +
  geom_vline(aes(xintercept = dataset_1[1,"X1"]))

ggplot() +
  geom_density(aes(x = as_tibble(posteriors$samples$chain1)%>%
                     select(`x_cont_miss[2, 1]`)%>%
                     dplyr::slice(1000:2500)%>%
                     unlist())) +
  geom_vline(aes(xintercept = dataset_1[501,"X1"]))

## ----message = FALSE, warning = FALSE-----------------------------------------
data(dataset_2)

GGally::ggpairs(dataset_2)

## ----message = FALSE, warning = FALSE-----------------------------------------
posteriors <- runModel(dataset_2, 
                       mcmc_iterations = 2500, 
                       L = 6, 
                       mcmc_chains = 2)

## ----message = FALSE, warning = FALSE-----------------------------------------
plot_ggpairs(posteriors, 
             newdata = dataset_2, 
             nburn = 1500)

## ----message = FALSE, warning = FALSE-----------------------------------------
rows <- 950:1000
dataset_missing <- dataset_2
dataset_missing_predict <- dataset_missing[rows,]
dataset_missing_predict[,"X2"] <- as.numeric(NA)

## ----message = FALSE, warning = FALSE-----------------------------------------
posteriors.prediction <- predict_dpmm_fit(posteriors, 
                                          newdata = dataset_missing_predict, 
                                          samples = seq(1500,2500, 25))

## ----message = FALSE, warning = FALSE-----------------------------------------
low_q <- NULL
median_q <- NULL
high_q <- NULL
for (i in 1:length(posteriors.prediction)) {
  low_q <- c(low_q, quantile(unlist(posteriors.prediction[[i]][]), probs = c(0.05)))
  median_q <- c(median_q, quantile(unlist(posteriors.prediction[[i]][]), probs = c(0.5)))
  high_q <- c(high_q, quantile(unlist(posteriors.prediction[[i]][]), probs = c(0.95)))
}

## ----message = FALSE, warning = FALSE-----------------------------------------
standardised_values <- (dataset_missing[rows,"X2"]-posteriors$mean_values[2])/posteriors$sd_values[2]

## ----message = FALSE, warning = FALSE-----------------------------------------
cbind(standardised_values,low_q, median_q,high_q) %>%
  as.data.frame() %>%
   rlang::set_names(c("observed", "low","median","high")) %>%
  ggplot() +
  geom_point(aes(x = observed, y = median)) +
  geom_errorbar(aes(x = observed, ymin = low, ymax = high)) +
  geom_abline(aes(intercept = 0, slope = 1))

## ----message = FALSE, warning = FALSE-----------------------------------------
rows <- 950:1000
dataset_missing <- dataset_2
dataset_missing_predict <- dataset_missing[rows,]
dataset_missing_predict[,"X3"] <- factor(NA, levels = levels(dataset_2[,"X3"]))

## ----message = FALSE, warning = FALSE-----------------------------------------
posteriors.prediction <- predict_dpmm_fit(posteriors, 
                                          newdata = dataset_missing_predict, 
                                          samples = seq(1500,2500, 25))

## ----message = FALSE, warning = FALSE-----------------------------------------
predicted <- NULL
for (i in 1:length(posteriors.prediction)) {
  var <- table(unlist(posteriors.prediction[i]))
  predicted <- c(predicted, names(var[var==max(var)]))
}
predicted <- as.numeric(predicted)

cbind(dataset_2[rows,"X3"],predicted) %>%
  as.data.frame() %>%
   rlang::set_names(c("observed", "predicted")) %>%
  mutate(observed = factor(observed),
         predicted = factor(predicted)) %>%
  table()

## ----message = FALSE, warning = FALSE-----------------------------------------
dataset_missing <- dataset_2
dataset_missing[2,"X3"] <- NA
dataset_missing[501,"X3"] <- NA

## ----message = FALSE, warning = FALSE-----------------------------------------
posteriors <- runModel(dataset_missing, 
                       mcmc_iterations = 2500,
                       L = 6, 
                       standardise = FALSE)

## ----message = FALSE, warning = FALSE-----------------------------------------
imputed_values <- as_tibble(posteriors$samples$chain1)%>%
  select(`x_disc_miss[1, 1]`)%>%
  dplyr::slice(1000:2500)%>%
  unlist()
cbind(imputed_values, rep(as.numeric(dataset_2[2,"X3"]), 1501)) %>%
  as.data.frame() %>%
   rlang::set_names(c("Predicted","True")) %>%
  table()

imputed_values <- as_tibble(posteriors$samples$chain1)%>%
  select(`x_disc_miss[2, 1]`)%>%
  dplyr::slice(1000:2500)%>%
  unlist()
cbind(imputed_values, rep(as.numeric(dataset_2[501,"X3"]), 1501)) %>%
  as.data.frame() %>%
   rlang::set_names(c("Predicted","True")) %>%
  table()


## ----message = FALSE, warning = FALSE-----------------------------------------
dataset_missing <- dataset_2
dataset_missing[1,"X2"] <- NA
dataset_missing[501,"X2"] <- NA

## ----message = FALSE, warning = FALSE-----------------------------------------
posteriors <- runModel(dataset_missing, 
                       mcmc_iterations = 2500,
                       L = 6, 
                       standardise = FALSE)

## ----message = FALSE, warning = FALSE-----------------------------------------
imputed_values <- as_tibble(posteriors$samples$chain1)%>%
  select(`x_cont_miss[1, 2]`)%>%
  dplyr::slice(1000:2500)%>%
  unlist()
ggplot() +
  geom_density(aes(x = imputed_values)) +
  geom_vline(aes(xintercept = dataset_2[1,"X2"]))

imputed_values <- as_tibble(posteriors$samples$chain1)%>%
  select(`x_cont_miss[2, 2]`)%>%
  dplyr::slice(1000:2500)%>%
  unlist()
ggplot() +
  geom_density(aes(x = imputed_values)) +
  geom_vline(aes(xintercept = dataset_2[501,"X2"]))


## ----message = FALSE, warning = FALSE-----------------------------------------
dataset_missing <- dataset_2
dataset_missing[c(1,501),c("X1", "X2")] <- NA

## ----message = FALSE, warning = FALSE-----------------------------------------
posteriors <- runModel(dataset_missing, 
                       mcmc_iterations = 2500,
                       L = 6, 
                       standardise = FALSE)

## ----message = FALSE, warning = FALSE-----------------------------------------
ggplot() +
  geom_point(aes(x = as_tibble(posteriors$samples$chain1)%>%
                   select(`x_cont_miss[1, 1]`)%>%
                   dplyr::slice(1000:2500)%>%
                   unlist(), 
                 y = as_tibble(posteriors$samples$chain1)%>%
                   select(`x_cont_miss[1, 2]`)%>%
                   dplyr::slice(1000:2500)%>%
                   unlist()), colour = "darkorange") +
  geom_point(aes(x = dataset_2[1,"X1"], y = dataset_2[1,"X2"]), 
             colour = "darkorange3", 
             fill = "black", 
             size = 3) +
  geom_point(aes(x = as_tibble(posteriors$samples$chain1)%>%
                   select(`x_cont_miss[2, 1]`)%>%
                   dplyr::slice(1000:2500)%>%
                   unlist(), 
                 y = as_tibble(posteriors$samples$chain1)%>%
                   select(`x_cont_miss[2, 2]`)%>%
                   dplyr::slice(1000:2500)%>%
                   unlist()), colour = "royalblue") +
  geom_point(aes(x = dataset_2[501,"X1"], y = dataset_2[501,"X2"]), 
             colour = "royalblue4", 
             fill = "black", 
             size = 3)

## ----message = FALSE, warning = FALSE-----------------------------------------
data(dataset_3)

GGally::ggpairs(dataset_3)

## ----message = FALSE, warning = FALSE-----------------------------------------
posteriors <- runModel(dataset_3, 
                       mcmc_iterations = 2500, 
                       L = 6, 
                       mcmc_chains = 2)

## ----message = FALSE, warning = FALSE-----------------------------------------
## PROBLEM HERE
plot_ggpairs(posteriors, 
             newdata = dataset_3, 
             nburn = 1500)

## ----message = FALSE, warning = FALSE-----------------------------------------
rows <- 950:1000
dataset_missing <- dataset_3
dataset_missing_predict <- dataset_missing[rows,]
dataset_missing_predict[,"X1"] <- factor(NA, levels = levels(dataset_3[,"X1"]))

## ----message = FALSE, warning = FALSE-----------------------------------------
posteriors.prediction <- predict_dpmm_fit(posteriors, 
                                          newdata = dataset_missing_predict, 
                                          samples = seq(1500,2500, 25))

## ----message = FALSE, warning = FALSE-----------------------------------------
predicted <- NULL
for (i in 1:length(posteriors.prediction)) {
  var <- table(unlist(posteriors.prediction[i]))
  predicted <- c(predicted, names(var[var==max(var)]))
}
predicted <- as.numeric(predicted)

cbind(dataset_3[rows,"X1"],predicted) %>%
  as.data.frame() %>%
   rlang::set_names(c("observed", "predicted")) %>%
  mutate(observed = factor(observed),
         predicted = factor(predicted)) %>%
  table()

## ----message = FALSE, warning = FALSE-----------------------------------------
dataset_missing <- dataset_3
dataset_missing[1,"X1"] <- factor(NA, levels = levels(dataset_3[,"X1"]))
dataset_missing[501,"X1"] <- factor(NA, levels = levels(dataset_3[,"X1"]))

## ----message = FALSE, warning = FALSE-----------------------------------------
posteriors <- runModel(dataset_missing, 
                       mcmc_iterations = 2500,
                       L = 6, 
                       standardise = FALSE)

## ----message = FALSE, warning = FALSE-----------------------------------------
imputed_values <- as_tibble(posteriors$samples$chain1)%>%
  select(`x_disc_miss[1, 1]`)%>%
  dplyr::slice(1000:2500)%>%
  unlist()
cbind(imputed_values, dataset_3[1,"X1"]) %>%
  as.data.frame() %>%
   rlang::set_names(c("Prediction","True Value")) %>%
  table()

imputed_values <- as_tibble(posteriors$samples$chain1)%>%
  select(`x_disc_miss[2, 1]`)%>%
  dplyr::slice(1000:2500)%>%
  unlist()
cbind(imputed_values, dataset_3[501,"X1"]) %>%
  as.data.frame() %>%
   rlang::set_names(c("Prediction","True Value")) %>%
  table()

## ----message = FALSE, warning = FALSE-----------------------------------------
data(dataset_4)

GGally::ggpairs(dataset_4)

## ----message = FALSE, warning = FALSE-----------------------------------------
posteriors <- runModel(dataset_4, 
                       mcmc_iterations = 2500, 
                       L = 6, 
                       standardise = FALSE)

## ----message = FALSE, warning = FALSE-----------------------------------------
plot_ggpairs(posteriors, 
             newdata = dataset_4, 
             nburn = 1500)

## ----message = FALSE, warning = FALSE-----------------------------------------
rows <- sample(1:nrow(dataset_4), 100)
dataset_missing <- dataset_4
dataset_missing_predict <- dataset_missing[rows,]
dataset_missing_predict[,"Y"] <- as.numeric(NA)

## ----message = FALSE, warning = FALSE-----------------------------------------
posteriors.prediction <- predict_dpmm_fit(posteriors, 
                                          dataset_missing_predict, 
                                          samples = seq(1500,2500, 25))

## ----message = FALSE, warning = FALSE-----------------------------------------
low_q <- NULL
median_q <- NULL
high_q <- NULL
for (i in 1:length(posteriors.prediction)) {
  low_q <- c(low_q, quantile(unlist(posteriors.prediction[[i]][]), probs = c(0.05)))
  median_q <- c(median_q, quantile(unlist(posteriors.prediction[[i]][]), probs = c(0.5)))
  high_q <- c(high_q, quantile(unlist(posteriors.prediction[[i]][]), probs = c(0.95)))
}

cbind(dataset_4[rows,"Y"],low_q, median_q,high_q) %>%
  as.data.frame() %>%
   rlang::set_names(c("observed", "low","median","high")) %>%
  ggplot() +
  geom_point(aes(x = observed, y = median)) +
  geom_errorbar(aes(x = observed, ymin = low, ymax = high)) +
  geom_abline(aes(intercept = 0, slope = 1))

## ----message = FALSE, warning = FALSE-----------------------------------------
data(dataset_5)

GGally::ggpairs(dataset_5)

## ----message = FALSE, warning = FALSE-----------------------------------------
dataset_missing <- dataset_5
dataset_missing[1,"X1"] <- NA
dataset_missing[2,"X3"] <- factor(NA, levels = levels(dataset_5[,"X3"]))

## ----message = FALSE, warning = FALSE-----------------------------------------
posteriors <- runModel(dataset_missing, 
                       mcmc_iterations = 2500, 
                       L = 6, 
                       mcmc_chains = 2, 
                       standardise = FALSE)

## ----message = FALSE, warning = FALSE-----------------------------------------
imputed_values <- as_tibble(posteriors$samples$chain1)%>%
  select(`x_cont_miss[1, 1]`)%>%
  dplyr::slice(1000:2500)%>%
  unlist()
ggplot() +
  geom_density(aes(x = imputed_values)) +
  geom_vline(aes(xintercept = dataset_5[1,"X1"]))

imputed_values <- as_tibble(posteriors$samples$chain1)%>%
  select(`x_disc_miss[2, 1]`)%>%
  dplyr::slice(1000:2500)%>%
  unlist()
cbind(imputed_values, rep(dataset_5[2,"X3"], rep = 1500)) %>%
  as.data.frame() %>%
   rlang::set_names(c("Predicted","True")) %>%
  table()

## ----message = FALSE, warning = FALSE-----------------------------------------
dataset_missing <- dataset_5
dataset_missing[1,c("X1", "X2")] <- NA
dataset_missing[1,"X3"] <- factor(NA, levels = levels(dataset_5[,"X3"]))

## ----message = FALSE, warning = FALSE-----------------------------------------
posteriors <- runModel(dataset_missing, 
                       mcmc_iterations = 2500, 
                       L = 6, 
                       standardise = FALSE)

## ----message = FALSE, warning = FALSE-----------------------------------------
imputed_values_1 <- as_tibble(posteriors$samples$chain1)%>%
  select(`x_cont_miss[1, 1]`)%>%
  dplyr::slice(1000:2500)%>%
  unlist()
imputed_values_2 <- as_tibble(posteriors$samples$chain1)%>%
  select(`x_cont_miss[1, 2]`)%>%
  dplyr::slice(1000:2500)%>%
  unlist()
ggplot() +
  geom_point(aes(x = imputed_values_1, y = imputed_values_2), alpha = 0.6) +
  geom_point(aes(x = dataset_5[1,"X1"], y = dataset_5[1,"X2"]))

# Table predicted values vs observed values
# Note: Variable X3 is well explained by X4
imputed_values <- as_tibble(posteriors$samples$chain1)%>%
        select(`x_disc_miss[1, 1]`)%>%
        dplyr::slice(1000:2500)%>%
        unlist()
cbind(imputed_values, rep(dataset_5[1,"X3"], rep = 1500)) %>%
  as.data.frame() %>%
   rlang::set_names(c("Predicted","True")) %>%
  table()

## ----message = FALSE, warning = FALSE-----------------------------------------
dataset_missing <- dataset_5
dataset_missing[1,"X1"] <- NA
dataset_missing[1,"X3"] <- factor(NA, levels = levels(dataset_5[,"X3"]))
dataset_missing[1,"X4"] <- factor(NA, levels = levels(dataset_5[,"X4"]))

## ----message = FALSE, warning = FALSE-----------------------------------------
posteriors <- runModel(dataset_missing, 
                       mcmc_iterations = 2500, 
                       L = 6, 
                       standardise = FALSE)

## ----message = FALSE, warning = FALSE-----------------------------------------
imputed_values <- as_tibble(posteriors$samples$chain1)%>%
  select(`x_disc_miss[1, 1]`)%>%
  dplyr::slice(1000:2500)%>%
  unlist()
cbind(imputed_values, rep(dataset_5[1,"X3"], rep = 1500)) %>%
  as.data.frame() %>%
   rlang::set_names(c("Predicted","True")) %>%
  table()

## ----message = FALSE, warning = FALSE-----------------------------------------
imputed_values <- as_tibble(posteriors$samples$chain1)%>%
  select(`x_disc_miss[1, 2]`)%>%
  dplyr::slice(1000:2500)%>%
  unlist()
cbind(imputed_values, rep(dataset_5[1,"X4"], rep = 1500)) %>%
  as.data.frame() %>%
   rlang::set_names(c("Predicted","True")) %>%
  table()

## ----message = FALSE, warning = FALSE-----------------------------------------
imputed_values <- as_tibble(posteriors$samples$chain1)%>%
  select(`x_cont_miss[1, 1]`)%>%
  dplyr::slice(1000:2500)%>%
  unlist()
ggplot() +
  geom_density(aes(x = imputed_values)) +
  geom_vline(aes(xintercept = dataset_5[1,"X1"]))

## ----message = FALSE, warning = FALSE-----------------------------------------
data(dataset_6)

GGally::ggpairs(dataset_6)

## ----message = FALSE, warning = FALSE-----------------------------------------
posteriors <- runModel(dataset_6, 
                       mcmc_iterations = 2500, 
                       L = 6, 
                       standardise = FALSE)

## ----message = FALSE, warning = FALSE-----------------------------------------
plot_ggpairs(posteriors, 
             newdata = dataset_6, 
             nburn = 1500)

## ----message = FALSE, warning = FALSE-----------------------------------------
rows <- 1:50
dataset_missing <- dataset_6
dataset_missing_predict <- dataset_missing[rows,]
dataset_missing_predict[,"X1"] <- as.numeric(NA)

## ----message = FALSE, warning = FALSE-----------------------------------------
posteriors.prediction <- predict_dpmm_fit(posteriors, 
                                          dataset_missing_predict, 
                                          samples = seq(1500,2500, 25))

## ----message = FALSE, warning = FALSE-----------------------------------------
low_q <- NULL
median_q <- NULL
high_q <- NULL
for (i in 1:length(posteriors.prediction)) {
  low_q <- c(low_q, quantile(unlist(posteriors.prediction[[i]][]), probs = c(0.05)))
  median_q <- c(median_q, quantile(unlist(posteriors.prediction[[i]][]), probs = c(0.5)))
  high_q <- c(high_q, quantile(unlist(posteriors.prediction[[i]][]), probs = c(0.95)))
}

cbind(dataset_6[rows,"X1"],low_q, median_q,high_q) %>%
  as.data.frame() %>%
   rlang::set_names(c("observed", "low","median","high")) %>%
  ggplot() +
  geom_point(aes(x = observed, y = median)) +
  geom_errorbar(aes(x = observed, ymin = low, ymax = high)) +
  geom_abline(aes(intercept = 0, slope = 1))

## ----message = FALSE, warning = FALSE-----------------------------------------
dataset_missing <- dataset_6
dataset_missing[1,"X1"] <- NA
dataset_missing[501,"X1"] <- NA

## ----message = FALSE, warning = FALSE-----------------------------------------
posteriors <- runModel(dataset_missing, 
                       mcmc_iterations = 2500, 
                       L = 6, 
                       standardise = FALSE)

## ----message = FALSE, warning = FALSE-----------------------------------------
imputed_values <- as_tibble(posteriors$samples$chain1)%>%
  select(`x_cont_miss[1, 1]`)%>%
  dplyr::slice(1000:2500)%>%
  unlist()
ggplot() +
  geom_density(aes(x = imputed_values)) +
  geom_vline(aes(xintercept = dataset_6[1,"X1"]))

imputed_values <- as_tibble(posteriors$samples$chain1)%>%
  select(`x_cont_miss[2, 1]`)%>%
  dplyr::slice(1000:2500)%>%
  unlist()
ggplot() +
  geom_density(aes(x = imputed_values)) +
  geom_vline(aes(xintercept = dataset_6[501,"X1"]))


