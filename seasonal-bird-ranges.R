# Exploration of bird data in USGS analysis
# ER Zylstra
# 17 March 2025

library(dplyr)
library(stringr)
library(rnpn) # Will need to check version
library(ebirdst)
library(terra)
library(tidyterra)
library(ggplot2)
library(cowplot)

# Species-phenophase-location combinations used in analysis -------------------#

df <- read.csv("data/birds_mammals_phenophases_locations.csv") %>%
  select(-kingdom)

# Which phenophases are included?
count(df, phenophase_description)

# Summary of species-php combinations
spp_ph <- df %>%
  group_by(common_name, phenophase_description) %>%
  summarize(n_sites = n_distinct(site_id),
            min_length = min(n),
            max_length = max(n),
            .groups = "keep") %>%
  data.frame()

# Summary by species
spp <- df %>%
  group_by(group, common_name) %>%
  summarize(n_php = n_distinct(phenophase_description),
            n_sites = n_distinct(site_id),
            n_series = n(),
            min_length = min(n),
            max_length = max(n),
            .groups = "keep") %>%
  data.frame()
count(spp, group)
# 71 bird spp, 15 mammal spp

# Exploring eBird products ----------------------------------------------------#

ebird_spp <- ebirdst_runs %>% data.frame()
str(ebird_spp)

# Using AMGO as an example

filter(ebird_spp, str_detect(str_to_lower(common_name), "american goldfinch"))
# Want to extract date ranges for each season (species specific) ###############

# There are 9-km or 27-km resolution maps, and smoothed/raw maps
ebirdst_download_status(species = "amegfi",
                        path = ebirdst_data_dir(), # Default. Might want to change
                        download_abundance = FALSE,
                        download_ranges = TRUE,
                        dry_run = FALSE,
                        pattern = "raw_9km",  # Only download raw to save space?
                        show_progress = TRUE)

amgo_range <- load_ranges(species = "amegfi", 
                          resolution = "9km", 
                          smoothed = FALSE,
                          path = ebirdst_data_dir())
amgo_range <- vect(amgo_range)

plot_b <- ggplot(data = subset(amgo_range, amgo_range$season == "breeding")) +
  geom_spatvector(aes(fill = season), fill = "#F8766D", show.legend = FALSE) +
  labs(title = element_text("Breeding"))
plot_nb <- ggplot(data = subset(amgo_range, amgo_range$season == "nonbreeding")) +
  geom_spatvector(aes(fill = season), fill = "#abd9e9", show.legend = FALSE) +
  labs(title = element_text("Nonbreeding"))
plot_m1 <- ggplot(data = subset(amgo_range, amgo_range$season == "prebreeding_migration")) +
  geom_spatvector(aes(fill = season), fill = "#ffffbf", show.legend = FALSE) +
  labs(title = element_text("Pre-breeding migration"))
plot_m2 <- ggplot(data = subset(amgo_range, amgo_range$season == "postbreeding_migration")) +
  geom_spatvector(aes(fill = season), fill = "#abdda4", show.legend = FALSE) +
  labs(title = element_text("Post-breeding migration"))

# Extract locations with AMGO observations
amgo_locs <- df %>%
  filter(common_name == "American goldfinch") %>%
  select(-group, -phenophase_description, -n) %>%
  distinct()
amgo_locs$loc <- as.numeric(rownames(amgo_locs))

# Add locations to maps
amgo_locs_v <- vect(amgo_locs, geom = c("longitude", "latitude"), 
                    crs = "epsg:4326", keepgeom = TRUE)
plot_b <- plot_b + geom_spatvector(data = amgo_locs_v)
plot_nb <- plot_nb + geom_spatvector(data = amgo_locs_v)
plot_m1 <- plot_m1 + geom_spatvector(data = amgo_locs_v)
plot_m2 <- plot_m2 + geom_spatvector(data = amgo_locs_v)

plot_grid(plot_m1, plot_b, plot_m2, plot_nb, nrow = 2)

# For each location, extract seasons bird is present
ext_ranges <- terra::extract(amgo_range, amgo_locs_v)
ranges <- left_join(ext_ranges, select(data.frame(amgo_locs_v), -common_name),
                    by = c("id.y" = "loc")) %>%
  group_by(species_code, common_name, site_id, latitude, longitude) %>%
  summarize(b = ifelse("breeding" %in% season, 1, 0),
            nb = ifelse("nonbreeding" %in% season, 1, 0),
            m1 = ifelse("prebreeding_migration" %in% season, 1, 0),
            m2 = ifelse("postbreeding_migration" %in% season, 1, 0),
            .groups = "keep") %>%
  data.frame()
ranges

