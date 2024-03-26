plot_and_save_last_iter <- function(data, data_location_fine, results_final, mmyy){
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
  
  # assign final cluster assignments to prevalence data
  dat_map <- cbind(data_location_fine, data) %>%
    cbind(results_final) %>%
    as.data.frame() %>%
    select(head(names(.), length(weeks)), tail(names(.), 1)) %>%
    rename(final_clusters = names(.)[ncol(.)])
  
  saveRDS(dat_map, file = paste0("data/processed/final_clustertrend_assignment_", 
                                 mmyy, ".rds"))
  
  
  # England map
  Eng_map <- st_read(dsn = "data/raw/Local_Authority_Districts_(December_2022)_Boundaries_UK_BFC/",
                     layer = "LAD_DEC_2022_UK_BFC")
  
  Eng_map <- subset(Eng_map, startsWith(LAD22CD, "E"))
  
  # final cluster assignment
  ltla_cluster <- dat_map
  
  a_temp <- ltla_cluster %>%
    rename(LAD22CD = location_fine)
  map_and_data_temp <- sp::merge(Eng_map, a_temp, by = "LAD22CD")
  a_temp_plot <- tm_shape(map_and_data_temp) +
    tm_fill("final_clusters", border.alpha = 0,
            title = mmyy, style = "cat", palette = "Set2") +
    tm_layout(legend.outside = TRUE, frame = FALSE)
  a_temp_plot
  
  print(a_temp_plot)
  
  tmap_save(a_temp_plot, paste0("outputs/maps/cluster_", mmyy,".png"))
  save(a_temp_plot, file = paste0("data/processed/cluster_map_", mmyy,".RData"))
  ggsave(p_cl_obs, file = paste0("outputs/trends/cluster_", mmyy,".pdf"))
  saveRDS(dat_map, file = paste0("data/processed/final_clustertrend_assignment_", 
                                 mmyy, ".rds"))
  
}



map_and_save_last_iter <- function(mmyy){
} 
