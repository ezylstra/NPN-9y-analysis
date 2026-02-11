# Download, process plant data
# 10 Feb 2026

library(rnpn)
library(ggplot2)
library(dplyr)
library(tidyr)
library(lubridate)
library(stringr)
library(terra)

# Download/Load individual phenometric data -----------------------------------#

# # Get taxonomic information from NPN database
# spp <- npn_species() %>% data.frame()
# count(spp, kingdom, class_id, class_common_name)
# # Class IDs = 13:16
# 
# # Download individual phenometric data
# dl_orig <- npn_download_individual_phenometrics(
#   request_source = "erinz",
#   class_ids = 13:16,
#   years = 2009:2025) %>%
#   data.frame()
#   
# write.csv(dl_orig,
#           "data/orig-downloads/plant-data-download.csv",
#           row.names = FALSE)

dat_orig <- read.csv("data/orig-downloads/plant-data-download.csv")

# Data processing steps -------------------------------------------------------#

# 1. Remove unnecessary columns to reduce size of dataset
# 2. Add broad phenophase groups/categories
# 3. Limit leaf budburst/senescence data to first/second half of year
# 4. Filter data to exclude sites outside of continental US (48 states)
# 5. Filter by phenophase
# 6. Remove series with fewer than 9 years of data
# 7. Identify species-state-phenophases that should be evaluated on 
#    something other than a calendar year

# Then, create three datasets (prior no within 7 days, 14 days, or no limits)
# For each:
  # Remove yeses that don't meet prior no criteria
  # Filter data so we just keep the first yes in calendar/water/summer year
  # Limit series to those that have first yeses in at least 9 years

# 1. Simplify dataset ---------------------------------------------------------#

dat <- dat_orig %>%
  select(-c(elevation_in_meters, kingdom, class_id, class_name,
            class_common_name, last_yes_year, last_yes_month, last_yes_day, 
            last_yes_doy, last_yes_julian_date, numdays_until_next_no))

# 2. Add broad phenophase categories ------------------------------------------#

# Create new individualID-phenophase column
dat$ind_phen <- paste0(dat$individual_id, "_", dat$phenophase_description)

# Load csv with broader phenophase categories
lookup_table <- read.csv("NPN_synthesis_phenophase_key_new.csv")

# Strip any leading/trailing whitespaces (which might cause phenophases to
# go unmatched)
lookup_table <- lookup_table %>%
  select(phenophase_description, phenophase) %>%
  mutate(across(phenophase_description:phenophase, str_trim))

# Add broader phenophase categories to data
dat <- dat %>%
  left_join(lookup_table, by = "phenophase_description")
# Check:
# count(dat, phenophase, phenophase_description)

# 3. Filter leaf data ---------------------------------------------------------#

# Remove leaf budburst first yeses that occurred after DOY 200
dat <- dat %>%
  filter(!(phenophase == "Leaf budburst" & first_yes_doy > 200))

# Remove colored or falling leaf first yeses that occurred before DOY 201
dat <- dat %>%
  filter(!(phenophase == "Colored leaves" & first_yes_doy < 201)) %>%
  filter(!(phenophase == "Falling leaves" & first_yes_doy < 201))

# 4. Filter data by location and species --------------------------------------#

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
  # Look at differences
    # dat %>%
    #   mutate(same = ifelse(state == state_new, 1, 0)) %>%
    #   count(same, state, state_new)
  # A few odd ones, but mostly all fine
# Select which state code to use and remove any sites outside lower 48
dat <- dat %>%
  mutate(state_new = case_when(
    !is.na(state_new) ~ state_new,
    !is.na(state) ~ state,
    .default = NA
  )) %>%
  select(-state) %>%
  rename(state = state_new) %>%
  filter(!is.na(state))

# Remove series for general "citrus" (species unknown; could be hybrids)
dat <- dat %>%
  filter(common_name != "citrus")

# 5. Filter by phenophase -----------------------------------------------------#

# We're only using 5 phenophases for analyses. Removing all others
dat <- dat %>%
  filter(phenophase %in% c("Leaf budburst",
                           "Open flowers", 
                           "Ripe fruits",
                           "Colored leaves",
                           "Falling leaves"))

# 6. Remove series with fewer than 9 years ------------------------------------#
# Doing this now to reduce the number of series we need to evalaute and make
# rules about yeartype for

dat <- dat %>%
  group_by(ind_phen) %>%
  mutate(series_yrs = n_distinct(first_yes_year)) %>%
  filter(series_yrs > 8) %>%
  ungroup() %>%
  data.frame()
# Series years includes all yeses, regardless of prior nos

# 6. Evaluating whether we need something other than calendar year ------------#

# Create functions to calculate day of wateryr and summeryr
wateryr_calc = function(x, start.month = 10){
  x = as.Date(x)
  start.yr = year(x) - 1*(month(x) < start.month)
  start.date = make_date(start.yr, start.month, 1)
  as.integer(x - start.date + 1)
}
summeryr_calc = function(x, start.month = 7){
  x = as.Date(x)
  start.yr = year(x) - 1*(month(x) < start.month)
  start.date = make_date(start.yr, start.month, 1)
  as.integer(x - start.date + 1)
}

# First, add a first_yes_date column to dataframe, along with summer/water years
# and dosy/dowy to make everything easier
  # Summer year = first_yes_year for Jan-Jun observations and first_yes_year + 1 for Jul-Dec
  # Water year = first_yes_year for Jan-Sep observations and first_yes_year + 1 for Oct-Dec
dat <- dat %>%
  mutate(first_yes_date = parse_date_time(x = paste(first_yes_year, first_yes_doy),
                                          orders = "yj")) %>%
  mutate(summer_year = case_when(
    first_yes_month %in% 7:12 ~ first_yes_year + 1,
    .default = first_yes_year
  )) %>%
  mutate(water_year = case_when(
    first_yes_month %in% 10:12 ~ first_yes_year + 1,
    .default = first_yes_year
  )) %>%
  mutate(dosy = summeryr_calc(first_yes_date),
         dowy = wateryr_calc(first_yes_date))

# Calendar year seems appropriate in most, but not all, cases. We'll create
# rules for species-phenophase combinations by state. Already filtered 
# leaf-related phenophases by date, so only need to evaluate this for open 
# flowers and ripe fruit phenophases

# Code below that's commented out was used to identify what 
# species-phenophase-state combinations should be changed from calendar year. 
# Used that to create a separate csv file  (flower-fruit-yeartype.csv) that 
# summarized findings.

# # Open flowers
#   # Identify earliest first yes in calendar year for each series
#   flowers_cal <- dat %>%
#     filter(phenophase == "Open flowers") %>%
#     arrange(ind_phen, first_yes_date) %>%
#     group_by(ind_phen, first_yes_year) %>%
#     mutate(earliest = ifelse(first_yes_julian_date == min(first_yes_julian_date),
#                              1, 0)) %>%
#     ungroup() %>%
#     data.frame() %>%
#     filter(earliest == 1)
#   # Summarize data by series
#     # Identify those that have first yeses in Dec and Jan
#     # Count number of unique months with first yeses
#     # Calculate the span between earliest and latest first yes
#   flowers_cal_s <- flowers_cal %>%
#     group_by(ind_phen, common_name, individual_id, site_id, state, 
#              phenophase_description) %>%
#     summarize(nyrs = n_distinct(first_yes_year),
#               decjan = ifelse(12 %in% first_yes_month & 1 %in% first_yes_month, 
#                               1, 0),
#               nmonths = n_distinct(first_yes_month),
#               gap = max(first_yes_doy) - min(first_yes_doy), 
#               .groups = "keep") %>% 
#     data.frame()
#   # Summarize by species and state to see which combinations have the highest 
#   # number of proportion of problematic series
#   flowers_prob <- flowers_cal_s %>%
#     group_by(state, common_name) %>%
#     summarize(nseries = n(),
#               n_decjan = sum(decjan == 1),
#               n_months6 = sum(nmonths >= 6),
#               n_gap180 = sum(gap >= 180),
#               .groups = "keep") %>%
#     data.frame() %>%
#     mutate(prop_decjan = round(n_decjan / nseries, 2),
#            prop_months6 = round(n_months6 / nseries, 2),
#            prop_gap180 = round(n_gap180 / nseries, 2)) %>%
#     filter(n_months6 + n_gap180 + n_decjan > 0)
#   # Look at summaries by state. Focus especially on species-states that have 
#   # high proportions of series with first yeses in Dec & Jan, then those that 
#   # have a high proportion of series with first yeses in 6 or more months. 
#   # Looked at iNat seasonality curves to see if there's evidence that flowering
#   # period might span more than one calendar year.
#   arrange(flowers_prob, common_name, state) 
# 
# # Ripe fruit
#   # Identify earliest first yes in calendar year for each series
#   fruits_cal <- dat %>%
#     filter(phenophase == "Ripe fruits") %>%
#     arrange(ind_phen, first_yes_date) %>%
#     group_by(ind_phen, first_yes_year) %>%
#     mutate(earliest = ifelse(first_yes_julian_date == min(first_yes_julian_date),
#                              1, 0)) %>%
#     ungroup() %>%
#     data.frame() %>%
#     filter(earliest == 1)
#   # Summarize data by series
#     # Identify those that have first yeses in Dec and Jan
#     # Count number of unique months with first yeses
#     # Calculate the span between earliest and latest first yes
#   fruits_cal_s <- fruits_cal %>%
#     group_by(ind_phen, common_name, individual_id, site_id, state, 
#              phenophase_description) %>%
#     summarize(nyrs = n_distinct(first_yes_year),
#               decjan = ifelse(12 %in% first_yes_month & 1 %in% first_yes_month, 
#                               1, 0),
#               nmonths = n_distinct(first_yes_month),
#               gap = max(first_yes_doy) - min(first_yes_doy), 
#               .groups = "keep") %>% 
#     data.frame()
#   # Summarize by species and state to see which combinations have the highest 
#   # number of proportion of problematic series
#   fruits_prob <- fruits_cal_s %>%
#     group_by(state, common_name) %>%
#     summarize(nseries = n(),
#               n_decjan = sum(decjan == 1),
#               n_months6 = sum(nmonths >= 6),
#               n_gap180 = sum(gap >= 180),
#               .groups = "keep") %>%
#     data.frame() %>%
#     mutate(prop_decjan = round(n_decjan / nseries, 2),
#            prop_months6 = round(n_months6 / nseries, 2),
#            prop_gap180 = round(n_gap180 / nseries, 2)) %>%
#     filter(n_months6 + n_gap180 + n_decjan > 0)
#   # Look at summaries. Focus on species-states that have high proportions of 
#   # series with first yeses in Dec & Jan, then those that have a high proportion 
#   # of series with first yeses in 6 or more months. Not too worried about
#   # combinations with very few series that have a large span of dates (gaps).
#   # Looked at iNat seasonality curves to see if there's evidence that fruiting
#   # period might span more than one calendar year.
#   arrange(fruits_prob, common_name, state)

# Load csv with yeartype for subset of flowering, fruiting phenophase series
yt <- read.csv("flower-fruit-yeartype.csv")

# Attach yeartype classification to data (and assuming calendar year for all
# series not included in the csv file)
dat <- dat %>%
  left_join(select(yt, phenophase, state, common_name, yeartype),
            by = c("phenophase", "state", "common_name")) %>%
  mutate(yeartype = ifelse(is.na(yeartype), "calendar", yeartype))
  
# Create day-of-period (DOP) variable that selects first_yes_doy for calendar
# year series, dosy for summer year series, and dowy for water year series. 
# Create year variable that selects first_yes_year for calendar year series 
# and summer_year/water_year for summer or water year series, respectively.
dat <- dat %>%
  mutate(dop = case_when(
    yeartype == "calendar" ~ first_yes_doy,
    yeartype == "summer" ~ dosy,
    yeartype == "water" ~ dowy,
    .default = NA
  )) %>%
  mutate(year = case_when(
    yeartype == "calendar" ~ first_yes_year,
    yeartype == "summer" ~ summer_year,
    yeartype == "water" ~ water_year,
    .default = NA
  ))

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
#           "data/out/plant-series-allyeses-thru2025.csv",
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
#           "data/out/plant-series-prior14-thru2025.csv",
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
#           "data/out/plant-series-prior7-thru2025.csv",
#           row.names = FALSE)
