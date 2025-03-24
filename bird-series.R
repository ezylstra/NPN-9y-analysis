# Exploration of bird data in USGS analysis
# ER Zylstra
# 20 March 2025

library(dplyr)
library(stringr)
library(lubridate)
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
  filter(!str_detect(common_name, "peeper"))

# Correct classification of titmice (birds, but listed as mammals)
df <- df %>%
  mutate(group = ifelse(str_detect(common_name, "titmouse"), "Bird", group))

# Correct spelling for catbird
df <- df %>%
  mutate(common_name = str_replace(common_name, "grey catbird", "gray catbird"))

# Exclude mammals
df_birds <- df %>%
  filter(group == "Bird")

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
# 72 bird spp

# Exploring eBird products ----------------------------------------------------#
# (Note: needed to get eBird key before downloading any products)

ebird_spp <- ebirdst_runs %>% data.frame()
# str(ebird_spp)

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
      p_b <- ggplot(data = subset(range, range$season == "breeding")) +
        geom_spatvector(aes(fill = season), 
                        fill = "#F8766D", show.legend = FALSE) +
        geom_spatvector(data = locs) +
        labs(title = element_text(paste0(cn, ": Breeding")))
      p_nb <- ggplot(data = subset(range, range$season == "nonbreeding")) +
        geom_spatvector(aes(fill = season), 
                        fill = "#abd9e9", show.legend = FALSE) +
        geom_spatvector(data = locs) +
        labs(title = element_text(paste0(cn, ": Nonbreeding")))
      p_m1 <- ggplot(data = subset(range, range$season == "prebreeding_migration")) +
        geom_spatvector(aes(fill = season), 
                        fill = "#ffffbf", show.legend = FALSE) +
        geom_spatvector(data = locs) +
        labs(title = element_text(paste0(cn, ": Pre-breeding migration")))
      p_m2 <- ggplot(data = subset(range, range$season == "postbreeding_migration")) +
        geom_spatvector(aes(fill = season), 
                        fill = "#abdda4", show.legend = FALSE) +
        geom_spatvector(data = locs) +
        labs(title = element_text(paste0(cn, ": Post-breeding migration")))
      
      print(plot_grid(p_m1, p_b, p_m2, p_nb, nrow = 2))
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
  # Could also do this for water or summer year, but not doing so now...
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
# (row for each observation series = species-site-phenophase)
df_birds <- df_birds %>%
  left_join(select(spp_ranges, common_name, site_id, 
                   resident, jan1, present_jan1), 
            by = c("common_name", "site_id"))

# Categorizing series by phenophase type --------------------------------------#

count(df_birds, phenophase_description) %>% arrange(desc(n))
# Live individuals       
# Calls or song (birds)
# Individuals at a feeding station
# Fruit/seed consumption
# Singing individuals (birds)
# Insect consumption
# Flower visitation
# Nest building (birds)

# Add phenophase_type label (making these explicit for easy changes)
df_birds <- df_birds %>%
  mutate(phenophase_type = case_when(
    phenophase_description == "Live individuals" ~ "presence",
    phenophase_description == "Individuals at a feeding station" ~ "presence",
    phenophase_description == "Calls or song (birds)" ~ "presence",
    phenophase_description == "Singing individuals (birds)" ~ "behavior",
    phenophase_description == "Nest building (birds)" ~ "behavior",
    phenophase_description == "Fruit/seed consumption" ~ "feeding",
    phenophase_description == "Insect consumption" ~ "feeding",
    phenophase_description == "Flower visitation" ~ "feeding",
    .default = NA
  ))

phpt <- df_birds %>%
  count(phenophase_description, phenophase_type) %>%
  arrange(phenophase_type) %>%
  mutate(rule_word = case_when(
    phenophase_type == "presence" ~ "Include if absent on 1 Jan",
    .default = "Include"
  ))



# Suggestions for including/excluding series in analyses ----------------------#

# If phenophase is a behavior: include
# If phenophase is related to feeding: include (see emails from Ellen)
# If phenophase indicates presence: include migratory species that should not be
  # at that location at the start of the year

df_birds <- df_birds %>%
  mutate(include = case_when(
    phenophase_type %in% c("behavior", "feeding") ~ 1,
    present_jan1 == 0 ~ 1,
    .default = 0
  ))
include <- filter(df_birds, include == 1)

# Looking for series with redundant information -------------------------------#

# Indications of species presence at the same site:
presence <- include %>%
  filter(phenophase_type == "presence") %>%
  group_by(common_name, site_id) %>%
  summarize(live = ifelse("Live individuals" %in% phenophase_description, 1, 0),
            station = ifelse("Individuals at a feeding station" %in% phenophase_description, 1, 0),
            calls = ifelse("Calls or song (birds)" %in% phenophase_description, 1, 0),
            .groups = "keep") %>%
  data.frame()

presence_count <- count(presence, live, station, calls)
# No species-sites where there is a feeding station series but no live individual series
# 10 species-sites where there are both series
# All combinations of live individuals and calls/song series

# Remove all "Individuals at feeding station" series since they're always
# accompanied by "Live individuals" series
include <- include %>%
  filter(phenophase_description != "Individuals at a feeding station")

# If a species-site combination has "Calls or song" and "Live individuals",
# then remove the "Calls or song" series
presence_dups <- presence %>%
  filter(calls == 1 & live == 1) %>%
  select(common_name, site_id) %>%
  mutate(phenophase_description = "Calls or song (birds)",
         remove = 1)
include <- include %>%
  left_join(presence_dups, 
            by = c("common_name", "site_id", "phenophase_description")) %>%
  filter(is.na(remove)) %>%
  select(-remove)

# Summarizing series that are left --------------------------------------------#

# Create "include2" variable to indicate whether a series should be included in 
# analyses after accounting for all factors
include <- include %>%
  rename(include2 = include)
df_birds <- df_birds %>%
  left_join(select(include, common_name, site_id, phenophase_description, include2),
            by = c("common_name", "site_id", "phenophase_description")) %>%
  mutate(include2 = replace_na(include2, replace = 0))

# Summarize by species
spp_include <- include %>%
  group_by(common_name, resident) %>%
  summarize(n_series = n(),
            n_sites = n_distinct(site_id),
            n_php_types = n_distinct(phenophase_type),
            n_php = n_distinct(phenophase_description),
            .groups = "keep") %>%
  arrange(common_name, .locale = "en") %>%
  data.frame()

# Summarize by phenophase
php_include <- include %>%
  group_by(phenophase_type, phenophase_description) %>%
  summarize(n_series = n(),
            n_spp = n_distinct(common_name),
            n_sites = n_distinct(site_id),
            .groups = "keep") %>%
  arrange(factor(phenophase_type, levels = c("presence", "behavior", "feeding")), desc(n_series)) %>%
  data.frame()

