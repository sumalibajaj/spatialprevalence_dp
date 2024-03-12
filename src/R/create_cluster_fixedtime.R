library(dplyr)
library(tidyr)
library(kernlab)
source("src/R/functions/function_dp_gibbs_univariate.R")

# ------------------------------------------------------------------------------
# GETTING PREVALENCE DATA READY FOR CLUSTERING
# ------------------------------------------------------------------------------
dat_og <- read.csv("data/raw/ltla_in_region_debiased_prev.csv")

# Code to select weeks for which you want to do clustering etc
weeks <- unique(dat_og$week_date)
weeks <- sample(weeks, 3)

for(i in weeks){
  print(paste0("Week ", i))
  # running clustering for each week separately and store results
  dat <- dat_og %>%
    filter(week_date == i)
  
  # ------------------------------------------------------------------------------
  # PARAMETERS AND INITIAL VALUES FOR DIRICHLET PROCESS GIBBS SAMPLER
  # ------------------------------------------------------------------------------
  # the base distribution here is Normal
  # set the parameter values and initial values to be used in the DP Gibbs sampler
  alpha <- 10 # concentration parameter >0 (smaller = few new clusters)
  mu0 <- 0 # mean of base distribution
  sigma0 <- 0.5 # sd of base distribution
  sigma <- 0.005 # sd of likelihood
  # c_init <- rep(1, nrow(dat)) # initially all points belong to the same cluster
  c_init <- seq(1:nrow(dat)) # initially all points belong to different clusters
  maxIters <- 50
  
  test_clusters <- alpha*log(1+(nrow(dat)/alpha))
  
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
  
  saveRDS(results_plot, file = paste0("data/processed/cluster_assignment_", i, ".rds"))
  
  
  # Cluster summary
  # Keep only columns of cluster assignments
  assignments <- results
  
  # Fake assignment matrix to see if code is doing the right thing
  # assignments <- matrix(c(10, 10, 21,
  #                         23, 5, 23), nrow = 3, ncol = 2)
  # maxIters <- 2
  
  # Create an empty matrix
  matrix_size <- nrow(assignments)
  matrix <- matrix(0, nrow = matrix_size, ncol = matrix_size)
  
  # 
  # Iterate over assignments in the second half
  for (col_i in ((maxIters %/% 2) + 1):maxIters) {
    assignment <- assignments[, col_i]
    
    # Iterate over characters and indices in the assignment
    for (ii in seq_along(assignment)) {
      for (jj in seq_along(assignment)) {
        # Check if indices are equal
        if (ii == jj) {
          matrix[ii, jj] <- matrix[ii, jj] + 1
        } else if (ii != jj && assignment[ii] == assignment[jj]) {
          matrix[ii, jj] <- matrix[ii, jj] + 1
        }
      }
    }
  }
  
  # Similarity matrix between 0 and 1
  matrix <- matrix/(maxIters-(maxIters %/% 2))
  
  # number of clusters at each iteration
  len_unique <- function(x) length(unique(x))
  num_clusters <- apply(results, MARGIN = 2, FUN = len_unique)
  mode_num_clusters <- names(sort(-table(num_clusters)))[1] %>% as.numeric() # MODE!
  mode_num_clusters
  
  # Perform spectral clustering on this similarity or affinity matrix
  final_clusters <- specc(matrix, centers = mode_num_clusters)
  final_clusters <- as.numeric(final_clusters)
  
  length(unique(final_clusters)) == mode_num_clusters # check final no. of clusters is as expected
  
  # assign final cluster assignments to prevalence data
  dat <- cbind(dat, final_clusters)
  
  saveRDS(dat, file = paste0("data/processed/final_cluster_assignment_", i, ".rds"))
  
}





