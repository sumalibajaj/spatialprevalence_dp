# Code modified from here:
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6583910/

library(dplyr)
library(ggplot2)
library(tidyr)
library(gganimate)

# ------------------------------------------------------------------------------
# DATA SIMULATION
# ------------------------------------------------------------------------------

# create some data
set.seed(11)
n <- 50
m1 <- 10
sd1 <- 1
y1 <- rnorm(n = n, mean = m1, sd = sd1)
m2 <- 11
sd2 <- 1
y2 <- rnorm(n = n, mean = m2, sd = sd2)
data <- data.frame(c(y1, y2)) %>% rename(values = 1)
hist(data$values, breaks = nrow(data)) # visualise the data


# ------------------------------------------------------------------------------
# PARAMETERS AND INITIAL VALUES FOR DIRICHLET PROCESS GIBBS SAMPLER
# ------------------------------------------------------------------------------

# set the parameter values and initial values to be used in the DP Gibbs sampler
alpha <- 1 # concentration parameter >0 (smaller = few new clusters)
mu0 <- 0 # mean of base distribution
sigma0 <- 2 # sd of base distribution
sigma <- 0.2 # sd of likelihood
c_init <- rep(1, nrow(data)) # initially all points belong to the same cluster
c_init <- seq(1:nrow(data)) # initially all points belong to different clusters
maxIters <- 50


# ------------------------------------------------------------------------------
# FUNCTIONS FOR GIBBS SAMPLING
# ------------------------------------------------------------------------------

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


# ------------------------------------------------------------------------------
# INFERENCE (RUNNING THE GIBBS SAMPLER)
# ------------------------------------------------------------------------------

results <- dp_gibbs(
  data = data, alpha = alpha, mu0 = mu0,
  sigma0 = sigma0, sigma = sigma, c_init = c_init,
  maxIters = maxIters
)


# ------------------------------------------------------------------------------
# PLOTTING THE CLUSTER ASSIGNMENTS
# ------------------------------------------------------------------------------

# Plot the cluster assignment
results_plot <- data %>%
  mutate(s.no = 1:n())
results_plot <- cbind(results_plot, c_init)
results_plot <- cbind(results_plot, results) %>%
  rename(`0` = c_init)

# Plot initial and later iteration
ggplot(data = results_plot, aes(x = s.no, y = values)) +
  geom_point(aes(color = as.factor(`0`)))
ggplot(data = results_plot, aes(x = s.no, y = values)) +
  geom_point(aes(color = as.factor(`50`)))

# Plot number of clusters by iterations
results_plot_iter <- pivot_longer(data = results_plot,
                                  cols = c(3:ncol(results_plot)),
                                  names_to = "iter",
                                  values_to = "cluster") %>%
  group_by(iter) %>%
  summarize(clusters = length(unique(cluster))) %>%
  mutate(iter = as.numeric(iter)) %>%
  arrange(iter)

ggplot(data = results_plot_iter, aes(x = iter, y = clusters)) +
  geom_line()
