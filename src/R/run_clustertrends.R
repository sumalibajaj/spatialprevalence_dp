# Run the clustering with outputs:
# cluster trend lines with observed data
# maps of clusters
source("src/R/create_clustertrends.R")
dat_og <- read.csv("data/raw/ltla_in_region_debiased_prev.csv")

# Code to select weeks for which you want to do clustering etc
weeks_monthyear <- tibble(weeks = unique(dat_og$week_date) ) %>%
  mutate(monthyear = format(as.Date(weeks), "%Y-%m"))

for(temp_monthyear in "2022-03"){
  print(temp_monthyear)
  temp_weeks <- weeks_monthyear %>% 
    filter(monthyear == temp_monthyear) %>%
    select(weeks) %>%
    pull() %>%
    as.Date()
  
  # Run for a given month
  if(length(temp_weeks) <= 1)
    print("only one week in this month, so moving on to next month")
  else
    create_clustertrends(dat_og = dat_og, weeks = temp_weeks, n_more_chains = 1)
}


