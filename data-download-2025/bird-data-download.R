# Download, process bird data
# 11 Feb 2026

library(rnpn)
library(ggplot2)
library(dplyr)
library(tidyr)
library(lubridate)
library(stringr)
library(ebirdst)
library(terra)

# Download/Load individual phenometric data -----------------------------------#

# dl_orig <- npn_download_individual_phenometrics(
#   request_source = "erinz",
#   functional_types = c("Bird"),
#   years = 2011:2025
# ) %>% data.frame()
# 
# # Write to file
# write.csv(dl_orig,
#           "data/orig-downloads/bird-data-download.csv",
#           row.names = FALSE)

dat_orig <- read.csv("data/orig-downloads/bird-data-download.csv")

# Alter one common name that's spelled wrong in NPN database
dat_orig <- dat_orig %>%
  mutate(common_name = case_when(
    common_name == "grey catbird" ~ "gray catbird",
    .default = common_name
  ))

# Data processing steps -------------------------------------------------------#

# 1. Remove unnecessary columns to reduce size of dataset
# 2. Add broad phenophase groups/categories
# 3. Filter by phenophase (remove dead animal series)
# 4. Filter data to exclude sites outside of continental US (48 states)
# 5. Figure out whether any species-state-phenophases should be evaluated on 
#    something other than calendar year
# 6. Remove series with fewer than 9 years of data
# 7. Remove presence-related phenophase series when the species is expected 
#    to be present at the beginning of each calendar year
# 8. Remove redundant series 

# Then, create three datasets (prior no within 7 days, 14 days, or no limits)
# For each:
  # Remove yeses that don't meet prior no criteria
  # Filter data so we just keep the first yes in each year
  # Limit series to those that have first yeses in at least 9 years

# 1. Simplify dataset ---------------------------------------------------------#

dat <- dat_orig %>%
  select(-c(elevation_in_meters, 
            last_yes_year, last_yes_month, last_yes_day, last_yes_doy, 
            last_yes_julian_date, numdays_until_next_no))

# 2. Add broad phenophase categories ------------------------------------------#

# Create new individualID-phenophase column
dat$ind_phen <- paste0(dat$individual_id, "_", dat$phenophase_description)

# Load csv with broader phenophase categories
lookup_table <- read.csv("NPN_synthesis_phenophase_key_new.csv")

##### Edit phenophase groups for Calls or song (birds) #######
# Per Ellen G: 
# We want to classify “Calls or song” as a presence-type phenophase, as 
# observers were instructed to record a yes for this phenophase when they heard
# a species vocalize but did not see it. This phenophase is distinct from 
# “Singing individuals”, which is used to record observations of bird 
# vocalizations associated with breeding, mating, or territorial behaviors.
lookup_table <- lookup_table %>%
  mutate(phenophase = ifelse(phenophase_description == "Calls or song (birds)",
                             "Adults observed", 
                             phenophase)) %>%
  mutate(phenophase = ifelse(phenophase_description == "Singing individuals (birds)", 
                             "Singing individuals",
                             phenophase))

# Strip any leading/trailing whitespaces (which caused things like "Courtship "
# to go unmatched)
lookup_table <- lookup_table %>%
  select(phenophase_description, phenophase) %>%
  mutate(across(phenophase_description:phenophase, str_trim))

# Add broader phenophase categories to data
dat <- dat %>%
  left_join(lookup_table, by = "phenophase_description")
# Check:
# count(dat, phenophase, phenophase_description)

# 3. Filter by phenophase -----------------------------------------------------#

# Remove series for observations of dead animals
dat <- dat %>%
  filter(phenophase != "Dead individuals")

# 4. Filter data by location --------------------------------------------------#

# Can only keep sites in the continental US (lower 48 states) due to climate
# data availability. Not all sites have the state listed, so we'll first need to
# add state and then filter. 

# Postal codes for lower 48 states
states48 <- state.abb[! state.abb %in% c("AK", "HI")]

# Load shapefile with US state boundaries
states <- vect("states/cb_2017_us_state_500k.shp")
# Reproject to WGS84, which is datum that NPN uses
states <- terra::project(states, "epsg:4326")
# Subset
states <- terra::subset(states, states$STUSPS %in% states48)

# Get state codes for sites with missing entries
dat <- dat %>%
  filter(is.na(state) | state %in% states48)
state_fill <- dat %>%
  select(site_id, longitude, latitude, state) %>%
  distinct()
state_fillv <- vect(state_fill, 
                    geom = c("longitude", "latitude"), 
                    crs = "epsg:4326")
state_new <- terra::extract(states, state_fillv)
state_fill <- cbind(state_fill, state_new = state_new$STUSPS)
# Attach to data
dat <- dat %>%
  left_join(select(state_fill, site_id, state_new), by = "site_id")
  # Look at differences:
  # dat %>%
  #   mutate(same = ifelse(state == state_new, 1, 0)) %>%
  #   count(same, state, state_new)
# Select which state code to use and remove any sites outside lower 48
dat <- dat %>%
  mutate(state_new = case_when(
    # Select new state code when present
    !is.na(state_new) ~ state_new,
    # Select old state code if present but new state code wasn't (likely
    # because location falls just outside state boundary in shapefile)
    !is.na(state) ~ state,
    # Otherwise, leave as NA (and will remove from dataset)
    .default = NA
  )) %>%
  select(-state) %>%
  rename(state = state_new) %>%
  filter(!is.na(state))

# 5. Remove series with <9 years of data --------------------------------------#

dat <- dat %>%
  group_by(ind_phen) %>%
  mutate(series_yrs = n_distinct(first_yes_year)) %>%
  filter(series_yrs > 8) %>%
  ungroup() %>%
  data.frame()
# Series years includes all yeses, regardless of prior nos

# 6. Use anything other than calendar year for birds? -------------------------#

# Looked at what we classified as summer year previously and it's not clear this
# is necessary, especially if we filter out series where the species is
# present during winter or year round. Will use calendar year for all.

# Create columns to match up with datasets for other functional groups
# that use summer or water year. Create day-of-period (DOP) variable that here,
# will be the same as first_yes_doy. Create year variable that here, will be the 
# same as first_yes_year.
dat <- dat %>%
  mutate(first_yes_date = parse_date_time(x = paste(first_yes_year, first_yes_doy),
                                          orders = "yj")) %>%
  mutate(yeartype = "calendar",
         dop = first_yes_doy,
         year = first_yes_year)

# 7. Use eBird data to identify series where species is present on Jan 1 ------#

# See this report for more information about investigation of bird data series:
# https://erinzylstra.quarto.pub/exploration-of-bird-data-for-usgs-analysis/

# The idea: If a species is a year-round resident, or if it migrates but spends 
# the first part of the year in a region that encompasses the site of interest, 
# then assessing trends in the date on which live individuals were first 
# observed at that site does not provide phenological insights. In these 
# instances, the first date a phenophase was observed probably tells us more 
# about when observers visited the site than about species ecology. For 
# migratory species, we'll use seasonal range maps from eBird to identify 
# whether a species was likely to be present at a site at the beginning of the 
# year. To do this, we'll download seasonal range maps (at a 9-km resolution) 
# for all species classified as migratory. For each species, we'll determine 
# which season (breeding, nonbreeding, pre-breeding migration, post-breeding 
# migration) overlapped January 1st. Then, by inspecting whether a site was 
# located within the expected range for that season, we can determine whether 
# the species was likely to be present at that site at the beginning of the 
# year.

# So for presence-related phenophases (Live individuals, Individuals at a 
# feeding station, Calls or song) we'll exclude series for species that are 
# likely present at that site on January 1st. Note that we are classifying 
# “Calls or song” as a presence-type phenophase, as observers were instructed to 
# record a yes for this phenophase when they heard a species vocalize but did 
# not see it. This phenophase is distinct from “Singing individuals”, which is 
# used to record observations of bird vocalizations associated with breeding, 
# mating, or territorial behaviors.

# Extract presence-related phenophases
datp <- dat %>%
  filter(phenophase == "Adults observed")

# Get list of eBird species
ebird_spp <- ebirdst_runs %>% data.frame()
# Change one common name so we can match things up easily
ebird_spp <- ebird_spp %>%
  mutate(common_name = case_when(
    common_name == "Northern/Southern House Wren" ~ "House Wren",
    .default = common_name))

# Dataframe with NPN bird species
birds <- datp %>%
  distinct(common_name) %>%
  # Make common names all lowercase to match up easily with eBird common_names
  mutate(common_name_l = str_to_lower(common_name)) %>%
  mutate(ebird = ifelse(common_name_l %in% str_to_lower(ebird_spp$common_name), 
                        1, 0))
# Check that all birds in NPN dataset have ebird range maps
# filter(birds, ebird == 0) 

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
#                           path = ebirdst_data_dir(),
#                           download_abundance = FALSE,
#                           download_ranges = TRUE,
#                           dry_run = FALSE,
#                           pattern = "raw_9km",
#                           show_progress = TRUE)
# }

# Loop through species and extract range information
for (i in 1:nrow(birds)) {
  
  cn <- birds$common_name[i]
  code <- birds$species_code[i]
  
  # Extract locations and convert to SpatVector
  locs <- datp %>% 
    filter(common_name == cn) %>%
    select(site_id, latitude, longitude) %>%
    distinct()
  # A few site locations fell in "holes" in polygon layers for some species,
  # so adjusting longitude just slightly to extract values
  locs <- locs %>%
    mutate(longitude = case_when(
      site_id == 27407 ~ -122.0,
      site_id == 10145 ~ -83.35,
      site_id == 6091 ~ -118.85,
      site_id == 2803 ~ -105.5,
      .default = longitude
    ))
  
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

    # For each location, extract seasons bird is present
    ext_ranges <- terra::extract(range, locs)
    ranges <- left_join(select(ext_ranges, -common_name), 
                        data.frame(locs),
                        by = c("id.y" = "loc")) %>%
      group_by(species_code, site_id, latitude, longitude) %>%
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
# Check if at least one season encompasses Jan 1
mig_birdspp <- mig_birds %>%
  group_by(common_name) %>%
  summarize(n_jan1 = sum(jan1)) %>%
  data.frame()
filter(mig_birdspp, n_jan1 == 0)
# For 3 species, Jan 1 falls between fall migration and non-breeding season
# according to eBird dates. We'll assume Jan 1 is in non-breeding season since 
# that's true for majority of species
mig_birdspp <- mig_birds %>%
  group_by(common_name) %>%
  mutate(n_jan1 = sum(jan1)) %>%
  ungroup() %>%
  mutate(jan1 = case_when(
    n_jan1 == 0 & season == "nb" ~ 1,
    .default = jan1
  )) %>%
  group_by(common_name) %>%
  summarize(jan1 = season[jan1 == 1]) %>%
  data.frame()

# Add season to spp_ranges dataframe
spp_ranges <- spp_ranges %>%
  left_join(select(birds, common_name, species_code), by = "species_code") %>%
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
# with presence-related phenophases and exclude series when the species is 
# expected to be presesnt on Jan 1.
datp <- datp %>%
  left_join(select(spp_ranges, common_name, site_id, resident, present_jan1),
            by = c("common_name", "site_id")) %>%
  filter(present_jan1 == 0)

# Merge all series data back together
datnp <- dat %>%
  filter(phenophase != "Adults observed") %>%
  mutate(resident = NA,
         present_jan1 = NA)
dat <- rbind(datp, datnp)

# 8. Remove redundant series --------------------------------------------------#
# For some species, there may be instances where a positive observation of one 
# phenophase is usually or always associated with a positive observation of 
# another phenophase. For instance, if an observer reported that there were 
# individuals at a feeding station, they likely also reported that they observed 
# live individuals on the same day. Including data series for both of these 
# phenophases for the same species at the same site is likely to be redundant, 
# and thus, one of the data series should probably be excluded from analyses. 
# That said, all series should be included if there are reasons to believe that 
# information from related series are independent and provide unique insights 
# into species phenology. 

# See whether same site_ids have individuals at feeding station and live
# individuals, since that would be redundant
station_redundant <- dat %>%
  group_by(common_name, site_id) %>%
  summarize(feeding = ifelse("Individuals at a feeding station" %in% phenophase_description, 1, 0),
            live = ifelse("Live individuals" %in% phenophase_description, 1, 0),
            .groups = "keep") %>%
  data.frame()
count(station_redundant, feeding, live)
# For every species-site combination that had a data series for “Individuals at 
# a feeding station”, there was also a data series for “Live individuals”. We'll
# remove feeding station series.
station_redundant_r <- station_redundant %>% 
  filter(feeding == 1 & live == 1) %>%
  select(common_name, site_id) %>%
  mutate(phenophase_description = "Individuals at a feeding station", 
         remove_stationr = 1)
# Remove redundancies 
dat <- dat %>%
  left_join(station_redundant_r, 
            by = c("common_name", "site_id", "phenophase_description")) %>%
  mutate(remove_stationr = replace_na(remove_stationr, 0)) %>%
  filter(remove_stationr == 0)

# See whether same site_ids have calls or song (different than singing 
# individuals!) and live individuals, since that would be redundant
calls_redundant <- dat %>%
  group_by(common_name, site_id) %>%
  summarize(calls = ifelse("Calls or song (birds)" %in% phenophase_description, 1, 0),
            live = ifelse("Live individuals" %in% phenophase_description, 1, 0),
            .groups = "keep") %>%
  data.frame()
count(calls_redundant, calls, live)
# Many species-site combinations that had both calls or song and live individual
# series. Will remove calls or song series
calls_redundant_r <- calls_redundant %>% 
  filter(calls == 1 & live == 1) %>%
  select(common_name, site_id) %>%
  mutate(phenophase_description = "Calls or song (birds)", 
         remove_callsr = 1)
# Remove redundancies 
dat <- dat %>%
  left_join(calls_redundant_r, 
            by = c("common_name", "site_id", "phenophase_description")) %>%
  mutate(remove_callsr = replace_na(remove_callsr, 0)) %>%
  filter(remove_callsr == 0)

# Clean up dataframe
dat <- dat %>%
  select(-c(remove_stationr, remove_callsr))

# Create dataset with no restrictions on prior nos ----------------------------#

# Create new dataset with ID-phenophase-year column
dat_all <- dat %>%
  mutate(ind_phen_year = paste0(ind_phen, "_", year))

# Filter data to keep just the first yes in each year
dat_all <- dat_all %>% 
  group_by(ind_phen_year) %>%
  filter(dop == min(dop))

# Calculate the number of years for each series, and remove any with fewer than
# 9 years
dat_all <- dat_all %>%
  group_by(ind_phen) %>%
  mutate(series_yrs = n()) %>%
  filter(series_yrs > 8) %>% 
  ungroup() %>%
  data.frame()

# Write to file (keep commented out so we don't accidentally overwrite)
# write.csv(dat_all,
#           "data/out/bird-series-allyeses-thru2025.csv",
#           row.names = FALSE)

# Create dataset with yeses preceded by a no within 14 days -------------------#

# Create new dataset with ID-phenophase-year column
dat_14 <- dat %>%
  mutate(ind_phen_year = paste0(ind_phen, "_", year))

# Remove observations that did not have a prior "no" within 14 days
dat_14 <- dat_14 %>%
  filter(!is.na(numdays_since_prior_no) & numdays_since_prior_no <= 14)

# Filter data to keep just the first yes in each year
dat_14 <- dat_14 %>% 
  group_by(ind_phen_year) %>%
  filter(dop == min(dop))

# Calculate the number of years for each series, and remove any with fewer than
# 9 years
dat_14 <- dat_14 %>%
  group_by(ind_phen) %>%
  mutate(series_yrs = n()) %>%
  filter(series_yrs > 8) %>% 
  ungroup() %>%
  data.frame()

# Write to file (keep commented out so we don't accidentally overwrite)
# write.csv(dat_14,
#           "data/out/bird-series-prior14-thru2025.csv",
#           row.names = FALSE)

# Create dataset with yeses preceded by a no within 7 days --------------------#

# Create new dataset with ID-phenophase-year column
dat_7 <- dat %>%
  mutate(ind_phen_year = paste0(ind_phen, "_", year))

# Remove observations that did not have a prior "no" within 7 days
dat_7 <- dat_7 %>%
  filter(!is.na(numdays_since_prior_no) & numdays_since_prior_no <= 7)

# Filter data to keep just the first yes in each year
dat_7 <- dat_7 %>% 
  group_by(ind_phen_year) %>%
  filter(dop == min(dop))

# Calculate the number of years for each series, and remove any with fewer than
# 9 years
dat_7 <- dat_7 %>%
  group_by(ind_phen) %>%
  mutate(series_yrs = n()) %>%
  filter(series_yrs > 8) %>% 
  ungroup() %>%
  data.frame()

# Write to file (keep commented out so we don't accidentally overwrite)
# write.csv(dat_7,
#           "data/out/bird-series-prior7-thru2025.csv",
#           row.names = FALSE)
