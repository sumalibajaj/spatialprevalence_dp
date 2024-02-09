library(dplyr)
library(tidyr)
source("src/R/functions/function_dp_gibbs_univariate.R")

# ------------------------------------------------------------------------------
# GETTING PREVALENCE DATA READY FOR CLUSTERING
# ------------------------------------------------------------------------------
dat <- read.csv("data/raw/ltla_in_region_debiased_prev.csv")
dat <- dat %>%
  filter(week == 1)

# ------------------------------------------------------------------------------
# PARAMETERS AND INITIAL VALUES FOR DIRICHLET PROCESS GIBBS SAMPLER
# ------------------------------------------------------------------------------
# the base distribution here is Normal
# set the parameter values and initial values to be used in the DP Gibbs sampler
alpha <- 2 # concentration parameter >0 (smaller = few new clusters)
mu0 <- 0 # mean of base distribution
sigma0 <- 0.5 # sd of base distribution
sigma <- 0.5 # sd of likelihood
c_init <- rep(1, nrow(dat)) # initially all points belong to the same cluster
# c_init <- seq(1:nrow(dat)) # initially all points belong to different clusters
maxIters <- 100


# ------------------------------------------------------------------------------
# INFERENCE (RUNNING THE GIBBS SAMPLER)
# ------------------------------------------------------------------------------

results <- dp_gibbs(
  data = dat$mean_prev %>% as.data.frame(), alpha = alpha, mu0 = mu0,
  sigma0 = sigma0, sigma = sigma, c_init = c_init,
  maxIters = maxIters
)

results_plot <- dat %>%
  mutate(s.no = 1:n())
results_plot <- cbind(results_plot, c_init)
results_plot <- cbind(results_plot, results) %>%
  rename(`0` = c_init)

saveRDS(results_plot, file = "data/processed/cluster_assignment.rds")
