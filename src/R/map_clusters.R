library(ggplot2)
library(tmap) # to plot the map
library(sf)
# ------------------------------------------------------------------------------
# PLOTTING THE CLUSTER ASSIGNMENTS IN MAPS
# ------------------------------------------------------------------------------
weeks

# England map
Eng_map <- st_read(dsn = "data/raw/Local_Authority_Districts_(December_2022)_Boundaries_UK_BFC/",
                   layer = "LAD_DEC_2022_UK_BFC")

Eng_map <- subset(Eng_map, startsWith(LAD22CD, "E"))

filename_temp <- "data/processed/final_clustertrend_assignment.rds"
print(filename_temp)

# final cluster assignment
ltla_cluster <- readRDS(filename_temp)

a_temp <- ltla_cluster %>%
  rename(LAD22CD = location_fine)
map_and_data_temp <- sp::merge(Eng_map, a_temp, by = "LAD22CD")
a_temp_plot <- tm_shape(map_and_data_temp) +
  tm_fill("final_clusters", border.alpha = 0,
          title = i, style = "cat", palette = "Set2") +
  tm_layout(legend.outside = TRUE, frame = FALSE)
a_temp_plot

# tmap_save(a_temp_plot, paste0("outputs/cluster_",i,".png"))
# save(a_temp_plot, file = paste0("data/processed/cluster_map_",i,".RData"))


