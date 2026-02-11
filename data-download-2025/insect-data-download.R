# Download, process insect data
# 11 Feb 2026

library(rnpn)
library(ggplot2)
library(dplyr)
library(tidyr)
library(lubridate)
library(stringr)
library(terra)

# Download/Load individual phenometric data -----------------------------------#

# dl_orig <- npn_download_individual_phenometrics(
#   request_source = "erinz",
#   functional_types = c("Insect"),
#   years = 2011:2025
# ) %>% data.frame()
# 
# # Write to file
# write.csv(dl_orig,
#           "data-download-2025/data-original/insect-data-download.csv",
#           row.names = FALSE)

dat_orig <- read.csv("data-download-2025/data-original/insect-data-download.csv")

# Data processing steps -------------------------------------------------------#

# 1. Remove unnecessary columns to reduce size of dataset
# 2. Add broad phenophase groups/categories
# 3. Filter by phenophase (remove dead animal series)
# 4. Filter data to exclude sites outside of continental US (48 states)
# 5. Determine whether any species-state-phenophases should be evaluated on 
#    something other than calendar year

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
lookup_table <- read.csv("data-download-2025/NPN_synthesis_phenophase_key_new.csv")

# Strip any leading/trailing whitespaces (which might cause phenophases to
# go unmatched)
lookup_table <- lookup_table %>%
  select(phenophase_description, phenophase) %>%
  mutate(across(phenophase_description:phenophase, str_trim))

# Couple small edits to lookup table:
lookup_table <- lookup_table %>%
  mutate(phenophase = ifelse(phenophase == "#N/A", NA, phenophase))
lt_add <- data.frame(
  phenophase_description = "Mating (one on top or end to end)",
  phenophase = "Mating activity"
)
lookup_table <- rbind(lookup_table, lt_add)

# Add broader phenophase categories to data
dat <- dat %>%
  left_join(lookup_table, by = "phenophase_description")
# Check:
# count(dat, phenophase, phenophase_description)

# 3. Filter by phenophase -----------------------------------------------------#

# Remove series for observations of dead animals and drone cells 
dat <- dat %>%
  filter(!is.na(phenophase) & phenophase != "Dead individuals")

# 4. Filter data by location --------------------------------------------------#

# Can only keep sites in the continental US (lower 48 states) due to climate
# data availability. Not all sites have the state listed, so we'll first need to
# add state and then filter. 

# Postal codes for lower 48 states
states48 <- state.abb[! state.abb %in% c("AK", "HI")]

# Load shapefile with US state boundaries
states <- vect("data-download-2025/states/cb_2017_us_state_500k.shp")
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

# 5. Use anything other than calendar year for insects? -----------------------#

# No prior knowledge or expectation that summer or water year would be more 
# appropriate than calendar year. Will use calendar year for all.

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
#           "data-download-2025/data-out/insect-series-allyeses-thru2025.csv",
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
#           "data-download-2025/data-out/insect-series-prior14-thru2025.csv",
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
#           "data-download-2025/data-out/insect-series-prior7-thru2025.csv",
#           row.names = FALSE)
