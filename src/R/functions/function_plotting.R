plot_and_save_last_iter <- function(data, results_final){
  # Keeping prevalence data and final cluster assignment
  results_plot <- data %>%
    as.data.frame() %>%
    cbind(results_final) %>%
    select(head(names(.), length(weeks)), tail(names(.), 1)) %>%
    rename(final_clusters = names(.)[ncol(.)]) %>%
    pivot_longer(!final_clusters, names_to = "week_date", values_to = "mean_prev") %>%
    mutate(final_clusters = as.factor(final_clusters))
  
  results_plot_sum <- results_plot %>%
    group_by(week_date, final_clusters) %>%
    summarise(mean_prev = mean(mean_prev))

  
  p_cl_obs <- ggplot(data = results_plot, aes(x = week_date, y = mean_prev, group = final_clusters)) +
    geom_jitter(aes(color = final_clusters), show.legend = FALSE, alpha = 0.3) +
    geom_line(data = results_plot_sum, aes(color = final_clusters), show.legend = FALSE) +
    theme_bw() +
    scale_color_brewer(palette = "Set2")
  
  print(p_cl_obs)
  
  data_location_fine <- data[1]
  
  # assign final cluster assignments to prevalence data
  dat_map <- cbind(data_location_fine, data) %>%
    cbind(results_final) %>%
    as.data.frame() %>%
    select(head(names(.), length(weeks)), tail(names(.), 1)) %>%
    rename(final_clusters = names(.)[ncol(.)])
  
  saveRDS(dat_map, file = "data/processed/final_clustertrend_assignment.rds")
}
