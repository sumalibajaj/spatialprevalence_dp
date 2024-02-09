library(ggplot2)
library(gganimate)
# ------------------------------------------------------------------------------
# PLOTTING THE CLUSTER ASSIGNMENTS
# ------------------------------------------------------------------------------

results_plot <- readRDS("data/processed/cluster_assignment.rds")

# Plot initial and later iteration
ggplot(data = results_plot, aes(x = s.no, y = mean_prev)) +
  geom_point(aes(color = as.factor(`0`))) +
  theme(legend.position = "none")

ggplot(data = results_plot, aes(x = s.no, y = mean_prev)) +
  geom_point(aes(color = as.factor(`100`))) +
  theme(legend.position = "none")


# Plot number of clusters by iterations
results_plot_iter <- pivot_longer(data = results_plot,
                                  cols = c(14:ncol(results_plot)),
                                  names_to = "iter",
                                  values_to = "cluster") %>%
  group_by(iter) %>%
  summarize(clusters = length(unique(cluster))) %>%
  mutate(iter = as.numeric(iter)) %>%
  arrange(iter)

ggplot(data = results_plot_iter, aes(x = iter, y = clusters)) +
  geom_line()

