# Exploration of bird data in USGS analysis
# ER Zylstra
# 19 March 2025

library(dplyr)
library(stringr)
library(lubridate)
library(rnpn) # Will need to check version
library(ebirdst)
library(terra)
library(tidyterra)
library(ggplot2)
library(cowplot)

# Species-phenophase-location combinations used in analysis -------------------#

df <- read.csv("data/birds_mammals_phenophases_locations.csv") %>%
  select(-kingdom)

# Remove spring peeper (amphibian) from dataframe
df <- df %>%
  filter(!grepl("peeper", common_name))

# Correct spelling for catbird
df <- df %>%
  mutate(common_name = str_replace(common_name, "grey catbird", "gray catbird"))

# Exclude mammals
df_birds <- df %>%
  filter(group == "Bird")

# Which phenophases are included?
count(df_birds, phenophase_description) %>% arrange(desc(n))
  # Live individuals       
  # Calls or song (birds)
  # Individuals at a feeding station
  # Fruit/seed consumption
  # Singing individuals (birds)
  # Insect consumption
  # Flower visitation
  # Nest building (birds)

# Live individuals & Individuals at a feeding station just indicate PRESENCE, so
# we should probably exclude series for resident species or migratory species
# that are supposed to be present at that location at the beginning of the year

# Calls or song, Singing individuals, and Nest builing indicate a BEHAVIOR, so
# it seems like all series are worth including

# Not sure about the other phenophases: Fruit/seed consumption, 
# Insect consumption, and Flower visitation

# Summary of species-php combinations
spp_ph <- df_birds %>%
  group_by(common_name, phenophase_description) %>%
  summarize(n_sites = n_distinct(site_id),
            min_length = min(n),
            max_length = max(n),
            .groups = "keep") %>%
  data.frame()

# Summary by species
spp <- df_birds %>%
  group_by(group, common_name) %>%
  summarize(n_php = n_distinct(phenophase_description),
            n_sites = n_distinct(site_id),
            n_series = n(),
            min_length = min(n),
            max_length = max(n),
            .groups = "keep") %>%
  data.frame()
count(spp, group)
# 70 bird spp, 15 mammal spp

# Exploring eBird products ----------------------------------------------------#

ebird_spp <- ebirdst_runs %>% data.frame()
str(ebird_spp)

# From https://science.ebird.org/en/status-and-trends/faq#seasons:

# Resident (i.e., non-migratory) species are identified by having TRUE in the 
# is_resident column of ebirdst_runs, and these species are assessed across the 
# whole year rather than seasonally. 

# The seasonal dates define the weeks that fall within each season. Breeding and 
# non-breeding season dates are defined for each species as the weeks during 
# those seasons when the species’ population does not move. For this reason, 
# these seasons are also described as stationary periods. Migration periods are 
# defined as the periods of movement between the stationary non-breeding and 
# breeding seasons. Note that for many species these migratory periods include 
# not only movement from breeding grounds to non-breeding grounds, but also 
# post-breeding dispersal, molt migration, and other movements.

# A rating of 0 implies this season failed review and model results should not 
# be used at all for this period. Ratings of 1-3 correspond to a gradient of 
# more to less extrapolation and/or omission, and we often use a traffic light 
# analogy when referring to them:
# 1: low quality, extensive extrapolation and/or omission, but at least some 
# regions have estimates that are accurate; can be used with caution in certain 
# regions.
# 2: medium quality, some extrapolation and/or omission; use with caution.
# 3: high quality, very little or no extrapolation and/or omission.

# Want to extract date ranges for each season (species specific) ###############

# Extract range information for all non-resident species ----------------------#

# Don't need to do anything for resident species, as they're year-round in
# locations where they're present

# Logical indicating whether to display maps with seaonal ranges with NPN
# locations overlaid
map <- FALSE

# Dataframe with NPN bird species
birds <- spp %>%
  filter(group == "Bird") %>%
  select(common_name, n_php, n_sites, n_series) %>%
  # Make common names all lowercase to match up easily with eBird common_names
  mutate(common_name_l = str_to_lower(common_name)) %>%
  mutate(ebird = ifelse(common_name_l %in% str_to_lower(ebird_spp$common_name), 
                        1, 0))

# Check that all birds in NPN dataset have ebird range maps
filter(birds, ebird == 0)

# Attach ebird info to NPN bird list
ebirds <- ebird_spp %>%
  mutate(common_name = str_to_lower(common_name)) %>%
  select(-contains("trends"), -rsquared, -beta0, -scientific_name)
birds <- birds %>%
  left_join(ebirds, by = c("common_name_l" = "common_name"))

# Download range maps for non-resident NPN bird species
# For now, just going to download high-resolution raw maps. Won't download or
# overwrite file unless we set force = TRUE.
  # for (spp6 in birds$species_code[birds$is_resident == FALSE]) {
  #   ebirdst_download_status(species = spp6,
  #                           path = ebirdst_data_dir(), # Default. Might want to change
  #                           download_abundance = FALSE,
  #                           download_ranges = TRUE,
  #                           dry_run = FALSE,
  #                           pattern = "raw_9km",  # Only download raw to save space?
  #                           show_progress = TRUE)
  # }

# Loop through species and extract range information
for (i in 1:nrow(birds)) {
  
  cn <- birds$common_name[i]
  code <- birds$species_code[i]
  
  # Extract locations and convert to SpatVector
  locs <- df_birds %>% 
    filter(common_name == cn) %>%
    select(-group, -phenophase_description, -n) %>%
    distinct()
  
  # If species is migratory, extract information about seasonal ranges
  if (birds$is_resident[i] == FALSE) {
  
    # Convert NPN locations to SpatVector
    locs$loc <- as.numeric(rownames(locs))
    locs <- vect(locs, geom = c("longitude", "latitude"), 
                 crs = "epsg:4326", keepgeom = TRUE)

    # Load species range
    range <- load_ranges(species = code, 
                         resolution = "9km", 
                         smoothed = FALSE,
                         path = ebirdst_data_dir())
    range <- vect(range)
    
    if (map) {
      plot_b <- ggplot(data = subset(range, range$season == "breeding")) +
        geom_spatvector(aes(fill = season), fill = "#F8766D", show.legend = FALSE) +
        geom_spatvector(data = locs) +
        labs(title = element_text(paste0(cn, ": Breeding")))
      plot_nb <- ggplot(data = subset(range, range$season == "nonbreeding")) +
        geom_spatvector(aes(fill = season), fill = "#abd9e9", show.legend = FALSE) +
        geom_spatvector(data = locs) +
        labs(title = element_text(paste0(cn, ": Nonbreeding")))
      plot_m1 <- ggplot(data = subset(range, range$season == "prebreeding_migration")) +
        geom_spatvector(aes(fill = season), fill = "#ffffbf", show.legend = FALSE) +
        geom_spatvector(data = locs) +
        labs(title = element_text(paste0(cn, ": Pre-breeding migration")))
      plot_m2 <- ggplot(data = subset(range, range$season == "postbreeding_migration")) +
        geom_spatvector(aes(fill = season), fill = "#abdda4", show.legend = FALSE) +
        geom_spatvector(data = locs) +
        labs(title = element_text(paste0(cn, ": Post-breeding migration")))
      print(plot_grid(plot_m1, plot_b, plot_m2, plot_nb, nrow = 2))
    }
    
    # For each location, extract seasons bird is present
    ext_ranges <- terra::extract(range, locs)
    ranges <- left_join(select(ext_ranges, - common_name), 
                        data.frame(locs),
                        by = c("id.y" = "loc")) %>%
      group_by(species_code, common_name, site_id, latitude, longitude) %>%
      summarize(b = ifelse("breeding" %in% season, 1, 0),
                nb = ifelse("nonbreeding" %in% season, 1, 0),
                m1 = ifelse("prebreeding_migration" %in% season, 1, 0),
                m2 = ifelse("postbreeding_migration" %in% season, 1, 0),
                .groups = "keep") %>%
      mutate(resident = FALSE) %>%
      data.frame()
  
  } else {
    ranges <- cbind(species_code = code, locs) %>%
      mutate(b = 1,
             nb = 1,
             m1 = 1,
             m2 = 1,
             resident = TRUE)
  }
  
  # Merge info for all species
  if (i == 1) {
    spp_ranges <- ranges
  } else {
    spp_ranges <- rbind(spp_ranges, ranges)
  }
  message("Appended range information for ", cn)
}

# For each migratory species, identify season that encompasses start of calendar 
# year (Jan 1), 
  # Could also do this for water or summer year, but not doing now...
  # oct1_doy <- yday(as.Date("2022-10-01"))
  # jul1_doy <- yday(as.Date("2022-07-01"))

mig_birds_nb <- birds %>%
  filter(is_resident == FALSE) %>%
  select(common_name, contains("nonbreeding")) %>%
  mutate(season = "nb")
colnames(mig_birds_nb) <- str_remove_all(colnames(mig_birds_nb),
                                         "nonbreeding_")
mig_birds_m1 <- birds %>%
  filter(is_resident == FALSE) %>%
  select(common_name, contains("prebreeding_migration")) %>%
  mutate(season = "m1")
colnames(mig_birds_m1) <- str_remove_all(colnames(mig_birds_m1),
                                         "prebreeding_migration_") 
mig_birds_b <- birds %>%
  filter(is_resident == FALSE) %>%
  select(common_name, breeding_quality, breeding_start, breeding_end) %>%
  mutate(season = "b")
colnames(mig_birds_b) <- str_remove_all(colnames(mig_birds_b), "breeding_")
mig_birds_m2 <- birds %>%
  filter(is_resident == FALSE) %>%
  select(common_name, contains("postbreeding_migration")) %>%
  mutate(season = "m2")
colnames(mig_birds_m2) <- str_remove_all(colnames(mig_birds_m2),
                                         "postbreeding_migration_") 
  
mig_birds <- rbind(mig_birds_nb, mig_birds_m1, mig_birds_b, mig_birds_m2) %>%
  filter(!is.na(start)) %>%
  mutate(doy1 = yday(start),
         doy2 = yday(end)) %>%
  # When season overlaps Jan 1, make start date negative
  mutate(doy1 = ifelse(doy1 > doy2, -1 *(366 - doy1), doy1))

# Duplicate season that overlaps Jan 1 with start = doy and end = 365 + doy
mig_birds_add <- mig_birds %>%
  filter(doy1 < 0) %>%
  rename(doy1_old = doy1,
         doy2_old = doy2) %>%
  mutate(doy1 = yday(start),
         doy2 = 365 + yday(end)) %>%
  select(-c(doy1_old, doy2_old))

# Find which season overlaps Jan 1
mig_birds <- rbind(mig_birds, mig_birds_add) %>%
  rowwise() %>%
  mutate(jan1 = ifelse(1 %in% doy1:doy2, 1, 0)) %>%
  ungroup() %>%
  data.frame()
mig_birdspp <- mig_birds %>%
  group_by(common_name) %>%
  summarize(jan1 = season[jan1 == 1]) %>%
  data.frame()

# Add season to spp_ranges dataframe
spp_ranges <- spp_ranges %>%
  left_join(mig_birdspp, by = "common_name")

# For each species and location, indicate whether species is supposed to be 
# present on Jan 1
spp_ranges$present_jan1 <- NA
for (i in 1:nrow(spp_ranges)) {
  if (spp_ranges$resident[i] == TRUE) {
    spp_ranges$present_jan1[i] <- 1
  } else {
    season <- spp_ranges$jan1[i]
    spp_ranges$present_jan1[i] <- ifelse(spp_ranges[i, season] == 1, 1, 0)
  }
}

# Append information about presence of species on Jan 1 to original dataframe
# (with a row for each series = species-site-phenophase)
df_birds <- df_birds %>%
  left_join(select(spp_ranges, common_name, site_id, 
                   resident, jan1, present_jan1), 
            by = c("common_name", "site_id"))

count(df_birds, resident, present_jan1)
# 135 series with resident species (always present)
# 268 series with migratory species that were likely present as of Jan 1
# 113 series with migratory species that were likely absent as of Jan 1

# NEXT STEP:
# make suggestions for include/exclude that take phenophase into account
