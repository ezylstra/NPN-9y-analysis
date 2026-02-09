# Download, process insect data
# 3 Feb 2026

library(rnpn)
library(ggplot2)
library(dplyr)
library(tidyr)
library(lubridate)
library(stringr)

# Download/Load individual phenometric data -----------------------------------#

# dl_orig <- npn_download_individual_phenometrics(
#   request_source = "erinz",
#   functional_types = c("Insect"),
#   years = 2011:2025
# ) %>% data.frame()
# 
# # Write to file
# write.csv(dl_orig,
#           "data/orig-downloads/insect-data-download.csv",
#           row.names = FALSE)

dat <- read.csv("data/orig-downloads/insect-data-download.csv")

# Data processing steps -------------------------------------------------------#
# Note: evaluating first yeses in each calendar year for all series

# 1. Add broad phenophase groups/categories
# 2. Remove any dead animal series

# Then, create three datasets (prior no within 7 days, 14 days, or no limits)
# For each:
  # Remove yeses that don't meet prior no criteria
  # Filter data so we just keep the first yes in each year
  # Limit series to those that have first yeses in at least 9 years

# 1. Add broad phenophase categories ------------------------------------------#

# Create new individualID-phenophase column
dat$ind_phen <- paste0(dat$individual_id, "_", dat$phenophase_description)

# Load csv with broader phenophase categories
lookup_table <- read.csv("NPN_synthesis_phenophase_key_new.csv")

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

# 2. Remove dead animal series ------------------------------------------------#

# Remove series for observations of dead animals and drone cells 
dat <- dat %>%
  filter(!is.na(phenophase) & phenophase != "Dead individuals")

# Create dataset with no restrictions on prior nos ----------------------------#

# Create new dataset with ID-phenophase-year column
dat_all <- dat %>%
  mutate(ind_phen_year = paste0(ind_phen, "_", first_yes_year))
  
# Filter data to keep just the first yes in each year
dat_all <- dat_all %>% 
  group_by(ind_phen_year) %>%
  filter(first_yes_doy == min(first_yes_doy))

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
#           "data/out/insect-series-allyeses-thru2025.csv",
#           row.names = FALSE)

# Create dataset with yeses preceded by a no within 14 days -------------------#

# Create new dataset with ID-phenophase-year column
dat_14 <- dat %>%
  mutate(ind_phen_year = paste0(ind_phen, "_", first_yes_year))

# Remove observations that did not have a prior "no" within 14 days
dat_14 <- dat_14 %>%
  filter(!is.na(numdays_since_prior_no) & numdays_since_prior_no <= 14)

# Filter data to keep just the first yes in each year
dat_14 <- dat_14 %>% 
  group_by(ind_phen_year) %>%
  filter(first_yes_doy == min(first_yes_doy))

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
#           "data/out/insect-series-prior14-thru2025.csv",
#           row.names = FALSE)

# Create dataset with yeses preceded by a no within 7 days --------------------#

# Create new dataset with ID-phenophase-year column
dat_7 <- dat %>%
  mutate(ind_phen_year = paste0(ind_phen, "_", first_yes_year))

# Remove observations that did not have a prior "no" within 7 days
dat_7 <- dat_7 %>%
  filter(!is.na(numdays_since_prior_no) & numdays_since_prior_no <= 7)

# Filter data to keep just the first yes in each year
dat_7 <- dat_7 %>% 
  group_by(ind_phen_year) %>%
  filter(first_yes_doy == min(first_yes_doy))

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
#           "data/out/insect-series-prior7-thru2025.csv",
#           row.names = FALSE)
