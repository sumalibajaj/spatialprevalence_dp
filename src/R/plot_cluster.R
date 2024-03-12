library(ggplot2)
library(rgdal) # to read in the shapefile
library(tmap) # to plot the map
library(sf)

# ------------------------------------------------------------------------------
# PLOTTING THE CLUSTER ASSIGNMENTS
# ------------------------------------------------------------------------------
# all cluster assignments over iterations
results_plot <- readRDS("data/processed/cluster_assignment.rds")

# Plot initial and later iteration
ggplot(data = results_plot, aes(x = s.no, y = mean_prev)) +
  geom_point(aes(color = as.factor(`0`))) +
  theme(legend.position = "none")

ggplot(data = results_plot, aes(x = s.no, y = mean_prev)) +
  geom_point(aes(color = as.factor(`50`))) +
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

ggplot(data = results_plot_iter %>% filter(clusters<10), aes(x = iter, y = clusters)) +
  geom_line()


# ------------------------------------------------------------------------------
# PLOTTING THE CLUSTER ASSIGNMENTS IN MAPS
# ------------------------------------------------------------------------------
weeks

# England map
Eng_map <- st_read(dsn = "data/raw/Local_Authority_Districts_(December_2022)_Boundaries_UK_BFC/",
                   layer = "LAD_DEC_2022_UK_BFC")

Eng_map <- subset(Eng_map, startsWith(LAD22CD, "E"))

for (i in weeks){
  filename_temp <- paste0("data/processed/final_cluster_assignment_", i, ".rds")
  print(filename_temp)
  
  # final cluster assignment
  ltla_cluster <- readRDS(filename_temp)

  a_temp <- ltla_cluster %>%
    rename(LAD22CD = location_fine)
  map_and_data_temp <- sp::merge(Eng_map, a_temp, by = "LAD22CD")
  a_temp_plot <- tm_shape(map_and_data_temp) +
    tm_fill("final_clusters", border.alpha = 0,
            title = i, style = "cat") +
    tm_layout(legend.outside = TRUE, frame = FALSE)
  
  tmap_save(a_temp_plot, paste0("outputs/cluster_",i,".png"))
  save(a_temp_plot, file = paste0("data/processed/cluster_map_",i,".RData"))
}


# final cluster assignment
ltla_cluster <- readRDS("data/processed/final_cluster_assignment.rds")

# England map
Eng_map <- st_read(dsn = "data/raw/Local_Authority_Districts_(December_2022)_Boundaries_UK_BFC/",
                   layer = "LAD_DEC_2022_UK_BFC")

Eng_map <- subset(Eng_map, startsWith(LAD22CD, "E"))

a_temp <- ltla_cluster %>%
  rename(LAD22CD = location_fine)
map_and_data_temp <- sp::merge(Eng_map, a_temp, by = "LAD22CD")
a_temp_plot <- tm_shape(map_and_data_temp) +
  tm_fill("final_clusters", border.alpha = 0,
          style = "cat") +
  tm_layout(legend.outside = TRUE, frame = FALSE)
a_temp_plot


for(i in 1:length(dates_imp)){
  print(i)
  # print(weeks_imp[i])
  
  # temp dataset of specific week
  a_temp <- ltla_cluster %>%
    rename(LAD23CD = location_fine) %>%
    filter(week == weeks_imp[i])
  
  # merging this temp dataset with map data
  map_and_data_temp <- sp::merge(Eng_map, a_temp, by = "LAD22CD")
  
  a_temp_plot <- tm_shape(map_and_data_temp) +
    tm_fill("change", border.alpha = 0, title = dates_imp[i],
            style = "cont", palette = "-inferno") +
    tm_layout(legend.outside = TRUE, frame = FALSE)
  
  tmap_save(a_temp_plot, paste0("outputs/prev_variant_ltla_ratio_t",i,".png"))
  save(a_temp_plot, file = paste0("data/processed/prev_variant_ltla_ratio_t",i,".RData"))
}


tmap_grob(a_temp_plot)

ggsave("outputs/Figure1.pdf", AA, width = 10, height = 12)
