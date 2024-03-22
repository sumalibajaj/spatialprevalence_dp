library(dplyr)

# ------------------------------------------------------------------------------
# GETTING PREVALENCE DATA READY FOR CLUSTERING
# ------------------------------------------------------------------------------
dat_og <- read.csv("data/raw/ltla_in_region_debiased_prev.csv")



# Code to select weeks for which you want to do clustering etc
ltla_all <- unique(dat_og$location_fine)
ltla <- sample(ltla_all, 1)
dat <- dat_og %>% filter(location_fine == ltla)

ggplot(data = dat, aes(x = week_date, y = mean_prev)) +
  geom_point() +
  geom_line()
