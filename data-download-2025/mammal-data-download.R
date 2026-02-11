# Download, process mammal data
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
#   functional_types = c("Mammal"),
#   years = 2011:2025
# ) %>% data.frame()
# 
# # Write to file
# write.csv(dl_orig,
#           "data-download-2025/data-original/mammal-data-download.csv",
#           row.names = FALSE)

dat_orig <- read.csv("data-download-2025/data-original/mammal-data-download.csv")

# Data processing steps -------------------------------------------------------#

# 1. Remove unnecessary columns to reduce size of dataset
# 2. Add broad phenophase groups/categories
# 3. Filter by phenophase (remove dead animal series)
# 4. Filter data to exclude sites outside of continental US (48 states)
# 5. Identify species-state-phenophases that should be evaluated on 
#    something other than a calendar year
# 6. Remove series with fewer than 9 years of data
# 7. Remove presence-related phenophase series when the species is expected 
#    to be active at the beginning of each calendar year
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
lookup_table <- read.csv("data-download-2025/NPN_synthesis_phenophase_key_new.csv")

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

# 5. Evaluate whether we need something other than calendar year --------------#

# Calendar year seems appropriate in almost all cases. The exception that was
# identified previously was mating-related phenophases for elk in CO. Elk 
# mating season in CO is in Sept-Oct, but can be observed after the new
# year. Using summer year (Jul-Jun) is more appropriate than water year.

# # See distribution of first yeses for elk mating phenophases
# dat %>%
#   filter(common_name == "elk") %>%
#   filter(state == "CO") %>%
#   filter(phenophase == "Mating activity") %>%
#   ggplot() +
#   geom_point(aes(x = first_yes_year, y = first_yes_doy)) +
#   facet_grid(~ind_phen)

# Create function to calculate day of summeryr
summeryr_calc = function(x, start.month = 7){
  x = as.Date(x)
  start.yr = year(x) - 1*(month(x) < start.month)
  start.date = make_date(start.yr, start.month, 1)
  as.integer(x - start.date + 1)
}

# Add a first_yes_date column to dataframe, along with summer year and dosy to
# make everything easier
  # Summer year = first_yes_year for Jan-Jun observations and first_yes_year + 1 for Jul-Dec
dat <- dat %>%
  mutate(first_yes_date = parse_date_time(x = paste(first_yes_year, first_yes_doy),
                                          orders = "yj")) %>%
  mutate(summer_year = case_when(
    first_yes_month %in% 7:12 ~ first_yes_year + 1,
    .default = first_yes_year
  )) %>%
  mutate(dosy = summeryr_calc(first_yes_date))

# Attach yeartype classification to data 
dat <- dat %>%
  mutate(yeartype = case_when(
    common_name == "elk" & state == "CO" & phenophase == "Mating activity" ~ "summer",
    .default = "calendar"
  ))

# Create day-of-period (DOP) variable that is equal to first_yes_doy for
# calendar year series and dosy for summer year series. Create year variable 
# that is equal to first_yes_year for calendar year series and summer_year for 
# summer year series.
dat <- dat %>%
  mutate(dop = case_when(
    yeartype == "calendar" ~ first_yes_doy,
    yeartype == "summer" ~ dosy,
    .default = NA
  )) %>%
  mutate(year = case_when(
    yeartype == "calendar" ~ first_yes_year,
    yeartype == "summer" ~ summer_year,
    .default = NA
  ))

# 6. Remove series with <9 years of data --------------------------------------#

dat <- dat %>%
  group_by(ind_phen) %>%
  mutate(series_yrs = n_distinct(year)) %>%
  filter(series_yrs > 8) %>%
  ungroup() %>%
  data.frame()
# Series years includes all yeses, regardless of prior nos

# 7. Identify series where species is active on Jan 1 -------------------------#

# See this report for more information about investigation of mammal data series:
# https://erinzylstra.quarto.pub/exploration-of-mammal-data-for-usgs-analysis/

# The idea: If a species is active year-round, then assessing trends in the date 
# on which live individuals were first observed at that site does not provide 
# phenological insights. In these instances, the first date a phenophase was 
# observed probably tells us more about when observers visited the site than 
# about species ecology. We consulted outside sources to determine whether a 
# species was likely to be active year-round, or whether it was likely to be 
# inactive at the beginning of the year (because it was in hibernation or 
# entered torpor). We extracted information about seasonal activity from the 
# Peterson Field Guide to the Mammals of North America (2006 edition) and from 
# species accounts on University of Michigan’s Animal Diversity Web 
# animaldiversity.org, and noted whether each species at each site was likely to 
# be active or inactive at the beginning of the year.

# Understanding when species exhibit behaviors associated with breeding or 
# mating or when young animals are first observed and how the timing of these 
# phenophases may shift over time or in response to climate is of interest, 
# regardless of whether the species is active year-round or not. Thus, it seems 
# reasonable to include all of these phenophases in analyses. Similarly, it 
# would make sense to include all feeding phenophases that relate to seasonal 
# resources (e.g., nuts, fruits, seeds) in analyses.

# In contrast, we'll exclude presence-related phenophases (Live individuals) and 
# general "feeding" series from analyses for species at locations where they're 
# expected to be present and active on Jan 1.

# Extract presence-related phenophases
datp <- dat %>%
  filter(phenophase_description %in% c("Live individuals", "Feeding"))

# Load file with information about species' seasonal activity
mammal_traits <- read.csv("data-download-2025/mammal_species_traits.csv")

# Append information about seasonal activity to original dataframe with
# presence-related phenophases and exclude series when the species is expected
# to be active on Jan 1.
datp <- datp %>%
  left_join(select(mammal_traits, common_name, site_id, seasonally_inactive),
            by = c("common_name", "site_id")) %>%
  mutate(active_jan1 = ifelse(seasonally_inactive == "no", 1, 0)) %>%
  filter(active_jan1 == 0)

# Merge all series data back together
datnp <- dat %>%
  filter(!phenophase_description %in% c("Live individuals", "Feeding")) %>%
  mutate(seasonally_inactive = NA,
         active_jan1 = NA)
dat <- rbind(datp, datnp)

# 8. Remove redundant series --------------------------------------------------#
# For some species, there may be instances where a positive observation of one 
# phenophase is usually or always associated with a positive observation of 
# another phenophase. In particular, every time an observer reported seeing an 
# individual feeding (without specifying the food resource), they also reported 
# that they observed live individuals on the same day. Including data series for 
# both of these phenophases for the same species at the same site is likely to 
# be redundant, and thus, one of the data series should probably be excluded 
# from analyses.

# See whether same sites have feeding and live individuals series, since that 
# would be redundant
feeding_redundant <- dat %>%
  group_by(common_name, site_id) %>%
  summarize(feeding = ifelse("Feeding" %in% phenophase_description, 1, 0),
            live = ifelse("Live individuals" %in% phenophase_description, 1, 0),
            .groups = "keep") %>%
  data.frame()
count(feeding_redundant, feeding, live)
# For every species-site combination that had a "Feeding" data series, there was 
# also a data series for “Live individuals”. We'll remove feeding series.
feeding_redundant_r <- feeding_redundant %>% 
  filter(feeding == 1 & live == 1) %>%
  select(common_name, site_id) %>%
  mutate(phenophase_description = "Feeding", 
         remove_feedingr = 1)
# Remove redundancies 
dat <- dat %>%
  left_join(feeding_redundant_r, 
            by = c("common_name", "site_id", "phenophase_description")) %>%
  mutate(remove_feedingr = replace_na(remove_feedingr, 0)) %>%
  filter(remove_feedingr == 0)

# Clean up dataframe
dat <- dat %>%
  select(-remove_feedingr)

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
#           "data-download-2025/data-out/mammal-series-allyeses-thru2025.csv",
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
#           "data-download-2025/data-out/mammal-series-prior14-thru2025.csv",
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
#           "data-download-2025/data-out/mammal-series-prior7-thru2025.csv",
#           row.names = FALSE)
