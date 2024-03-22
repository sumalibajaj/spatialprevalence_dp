library(dplyr)
library(ggplot2)
library(kernlab)
library(rstan)
library(RColorBrewer)
source("src/R/functions/function_dp_gibbs_multivariate.R")
source("src/R/functions/function_summarise_clusterassignment.R") # summary over iters
source("src/R/functions/function_utility.R")
source("src/R/functions/function_plotting.R")

# To select the same number of initial clusters when running multiple chains
set.seed(12)
# ------------------------------------------------------------------------------
# GETTING PREVALENCE DATA READY FOR CLUSTERING
# ------------------------------------------------------------------------------
dat_og <- read.csv("data/raw/ltla_in_region_debiased_prev.csv")

# Code to select weeks for which you want to do clustering etc
weeks_all <- unique(dat_og$week_date)

# Delta invasion and establishment
weeks <- weeks_all[20:24]
  
data <- dat_og %>% 
  filter(week_date %in% weeks) %>% 
  dplyr::select(location_fine, week_date, mean_prev)

data <- reshape(data, idvar = "location_fine", timevar = "week_date", direction = "wide")
data_location_fine <- data[1]
data <- data[,-1] %>% as.matrix()

# ------------------------------------------------------------------------------
# PARAMETERS AND INITIAL VALUES FOR DIRICHLET PROCESS GIBBS SAMPLER
# ------------------------------------------------------------------------------
# the base distribution here is Normal
# set the parameter values and initial values to be used in the DP Gibbs sampler
alpha <- 10 # concentration parameter >0 (smaller = few new clusters)
mu0 <- matrix(rep(0, length(weeks)), ncol = length(weeks), byrow = TRUE) # mean of base distribution
sigma0 <- diag(length(weeks)) * 1 # sd of base distribution
sigma <- diag(length(weeks))/100000 # sd of likelihood
maxIters <- 150

test_clusters <- alpha*log(1+(nrow(data)/alpha))
# ------------------------------------------------------------------------------
# INFERENCE (RUNNING THE GIBBS SAMPLER)
# ------------------------------------------------------------------------------

# Chain 1 - initially all points belong to different clusters
results_c1 <- dp_gibbs(
  data = data, alpha = alpha, mu0 = mu0,
  sigma0 = sigma0, sigma = sigma, c_init = seq(1:nrow(data)),
  maxIters = maxIters
)

# Chain 2 - initially all points belong to the same cluster
results_c2 <- dp_gibbs(
  data = data, alpha = alpha, mu0 = mu0,
  sigma0 = sigma0, sigma = sigma, c_init = rep(1, nrow(data)),
  maxIters = maxIters
)

# Chain 3 and more - number of clusters fixed and points assigned accordingly
n_more_chains <- 3
# final list with all chains including the two above
results_cn <- vector(mode='list', n_more_chains + 2)
# these results have an extra column with number of initial clusters
results_cn <- run_multiple_chains(n_more_chains, 
                        data = data, alpha = alpha, mu0 = mu0,
                        sigma0 = sigma0, sigma = sigma, c_init = c_init_r,
                        maxIters = maxIters,
                        results_list = results_cn)

# add the results from chain 1 and 2
results_cn[[length(results_cn) - 1]] <- results_c1
results_cn[[length(results_cn)]] <- results_c2


# Calculate how many clusters is each chain create at the last iteration
final_iter_nclust_vec <- sapply(results_cn, length_unique_last_column)

# FINAL CHAIN - select and plot 
# Choose one chain proportional the frequency of final #chains
mode_nclust <- names(sort(-table(final_iter_nclust_vec)))[1] %>% as.numeric() # MODE!
# Get indices where value is equal to mode
indices <- which(final_iter_nclust_vec == mode_nclust)
# Randomly pick an index i.e. chain
random_index <- sample(indices, 1)
results_final <- results_cn[[random_index]]


plot_and_save_last_iter(data = data, results_final = results_final)
 
 
 
 
################################################################################ 
 # RELATED TO CONVERGENCE
  # Check if chains have converged
  n_clusters_c1 <- apply(results_c1, MARGIN = 2, FUN = length_unique)  
table(n_clusters_c1)
  n_clusters_c2 <- apply(results_c2, MARGIN = 2, FUN = length_unique)
table(n_clusters_c2)
  n_clusters_c3 <- apply(results_c3, MARGIN = 2, FUN = length_unique)
  
  # discarding first half of the samples
  n_clusters <- cbind(n_clusters_c1[round(maxIters/2):maxIters],
                      n_clusters_c2[round(maxIters/2):maxIters])
  Rhat(sims = n_clusters)

  
  # saveRDS(results_plot, file = paste0("data/processed/clustertrend_assignment_", i, ".rds"))
  

  names(results_plot)
  
  # RELATED TO SUMMARIZING WILL NEED FOR STATS ANALYSIS BUT IN DIFFERENT RSCRIPT
  
  result_final_cluster <- summary_cluster_assign(results = results_c1, 
                                                 maxIters = maxIters,
                                                 data = data)
  
  result_final_cluster1 <- result_final_cluster %>% 
    as.data.frame() %>%
    pivot_longer(!final_clusters, names_to = "week_date", values_to = "mean_prev") %>%
    group_by(week_date, final_clusters) %>%
    summarise(mean_prev = mean(mean_prev))
  
  ggplot(data = result_final_cluster1, aes(x = week_date, y = mean_prev, group = final_clusters)) +
    geom_line(aes(color = final_clusters), show.legend = FALSE)
  






