# ------------------------------------------------------------------------------
# FUNCTIONS FOR GIBBS SAMPLING
# ------------------------------------------------------------------------------

# Code modified from here:
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6583910/

library(dplyr)
library(tidyr)

# function for normalising a vector
normalized.vector <- function(vector) {
  return(vector / sum(vector))
}

# function for Gibbs sampling for the DP
dp_gibbs <- function(data, alpha, mu0, sigma0, sigma_y, c_init, maxIters) {
  data_n <- nrow(data)
  tau0 <- 1 / sigma0
  tau <- 1 / sigma
  
  # initialize the DP Gibbs sampler
  y_c <- c_init
  n_k <- as.vector(table(y_c)) # number of points in clusters
  c_n <- length(n_k) # total number of clusters
  
  # save cluster assignments for each point and iteration
  result <- matrix(NA, nrow = data_n, ncol = maxIters)
  
  # Stating DP Gibbs sampler
  for (iter in 1:maxIters) {
    if (iter %% 25 == 0){
      print(paste0("Iteration: ", iter)) 
    }
    for (i in 1:data_n) {
      # STEP1: remove data point from its current cluster
      # remove the data point under consideration from its cluster
      c_i <- y_c[i]
      n_k[c_i] <- n_k[c_i] - 1
      
      # if there are 0 datapoints in the cluster now, then do the following
      # so there is no empty cluster
      if (n_k[c_i] == 0) {
        n_k[c_i] <- n_k[c_n] # last cluster to replace this empty cluster
        y_c[y_c == c_n] <- c_i # move points from last cluster to empty cluster
        n_k <- n_k[-c_n] # throw away the last empty cluster now
        c_n <- c_n - 1
      }
      # makes sure this doesn't get counted as a cluster in sum_data step below
      y_c[i] <- -1
      
      # STEP2: Assign data point to a cluster
      # each potential cluster has associated probability of being joined
      c_join_lp <- rep(NA, c_n + 1)
      
      # 2a: join an existing cluster
      # how far is the point to all the exiting points in the cluster
      for (c_exist_i in 1:c_n) {
        # sum all the points in this cluster
        y_exist_i <- which(y_c == c_exist_i)
        sum_data <- sum(data[y_c == c_exist_i, ])
        mean_p <- (sum_data * tau + mu0 * tau0) / (n_k[c_exist_i] * tau + tau0)
        sd_p <- 1 / (n_k[c_exist_i] * tau + tau0) + sigma
        c_join_lp[c_exist_i] <- log(n_k[c_exist_i]) +
          dnorm(data[i, ], mean = mean_p, sd = sd_p, log = TRUE) - log(data_n - 1 + alpha)
      }
      
      # join a new cluster
      c_join_lp[c_n + 1] <- log(alpha) +
        dnorm(data[i, ], mean = mu0, sd = sigma0 + sigma, log = TRUE) - log(data_n - 1 + alpha)
      
      # which cluster will the point join
      cluster_join <- sample(1:(c_n + 1), 1, replace = TRUE, prob = normalized.vector(exp(c_join_lp)))
      
      # new cluster
      if (cluster_join == c_n + 1) {
        n_k <- c(n_k, 0)
        c_n <- c_n + 1
      }
      y_c[i] <- cluster_join
      n_k[cluster_join] <- n_k[cluster_join] + 1 # increase no. of points in cluster
    }
    result[, iter] <- y_c # cluster membership of N observations
  }
  return(result)
}
