# Looking at individual phenometrics data that were used for USGS analysis
# 4 Sep 2025
# ER Zylstra

library(tidyr)
library(dplyr)
library(lubridate)
library(stringr)
library(ggplot2)
library(rnpn)

# Create functions to calculate day of wateryr and summeryr
wateryr_calc = function(x, start.month = 10){
  x = as.Date(x)
  start.yr = year(x) - (month(x) < start.month)
  start.date = make_date(start.yr, start.month, 1)
  as.integer(x - start.date + 1)
}
summeryr_calc = function(x, start.month = 7){
  x = as.Date(x)
  start.yr = year(x) - (month(x) < start.month)
  start.date = make_date(start.yr, start.month, 1)
  as.integer(x - start.date + 1)
}

# Load csv file with data used for 9+ yr analysis
df_orig <- read.csv("first-yes-years/new_9year_alldata_mar2025.csv")

# Remove unnecssary columns and format things to make it simpler
df <- df_orig %>%
  select(-c(first_yes_julian_date, last_yes_julian_date, 
            gdd, gddf, elevation_in_meters,
            tmax_winter, tmax_spring, tmax_summer, tmax_fall, tmax, tmaxf,
            tmin_winter, tmin_spring, tmin_summer, tmin_fall, tmin, tminf,
            prcp_winter, prcp_spring, prcp_summer, prcp_fall, prcp,
            acc_prcp, daylength, ind_phen, ind_phen_year, MAT)) %>%
  rename(site = site_id,
         lat = latitude,
         lon = longitude, 
         id = individual_id)

# Are state labels always there and reliable?
count(df, state) # Note that there are some sites in Puerto Rico
# 405 blank, but they're all for one site in the UP in MI
count(filter(df, state == ""), site, lat, lon)
df <- df %>%
  mutate(state = ifelse(state == "", "MI", state))

# Summarize info for each "series" --------------------------------------------#
df <- df %>%
  mutate(seriesid = paste(common_name, id, phenophase_id, sep = "_"))

series <- df %>%
  mutate(yes_date = parse_date_time(x = paste(first_yes_year, first_yes_doy),
                                    orders = "yj")) %>%
  # First, calculate day of water year and day of summer year
  mutate(dowy = wateryr_calc(yes_date)) %>%
  mutate(dosy = summeryr_calc(yes_date)) %>%
  group_by(seriesid, common_name, genus, species, kingdom, id, site, state,
           phenophase_description) %>%
  summarize(nyrs = n(),
            # Number of unique months in series
            nmonths = n_distinct(first_yes_month),
            # Calculate length of time between min/max doy for each yeartype
            calyr_gap = max(first_yes_doy) - min(first_yes_doy),
            wateryr_gap = max(dowy) - min(dowy),
            summeryr_gap = max(dosy) - min(dosy),
            # Identify if there are any first yeses in each month
            fy_jul = ifelse(7 %in% first_yes_month, 1, 0),
            fy_aug = ifelse(8 %in% first_yes_month, 1, 0),
            fy_sep = ifelse(9 %in% first_yes_month, 1, 0),
            fy_oct = ifelse(10 %in% first_yes_month, 1, 0),
            fy_nov = ifelse(11 %in% first_yes_month, 1, 0),
            fy_dec = ifelse(12 %in% first_yes_month, 1, 0),
            fy_jan = ifelse(1 %in% first_yes_month, 1, 0),
            fy_feb = ifelse(2 %in% first_yes_month, 1, 0),
            fy_mar = ifelse(3 %in% first_yes_month, 1, 0),
            fy_apr = ifelse(4 %in% first_yes_month, 1, 0),
            fy_may = ifelse(5 %in% first_yes_month, 1, 0),
            fy_jun = ifelse(6 %in% first_yes_month, 1, 0),
            .groups = "keep") %>%
  mutate(fy_1 = ifelse(fy_jan + fy_feb + fy_mar > 0, 1, 0),
         fy_2 = ifelse(fy_apr + fy_may + fy_jun > 0, 1, 0),
         fy_3 = ifelse(fy_jul + fy_aug + fy_sep > 0, 1, 0),
         fy_4 = ifelse(fy_oct + fy_nov + fy_dec > 0, 1, 0)) %>%
  # Calculate how many quarters have first yeses
  mutate(nquarters = fy_1 + fy_2 + fy_3 + fy_4) %>%
  # Calculate the shortest length of time between min/max doy across yeartypes
  mutate(min_gap = min(calyr_gap, wateryr_gap, summeryr_gap)) %>%
  select(-c(fy_jun, fy_jul, fy_aug, fy_mar, fy_apr, fy_may)) %>%
  data.frame()

# Are all valid species? ------------------------------------------------------#

# count(series, species)
series %>%
  filter(species == "spp.") %>%
  count(kingdom, genus, species, common_name, state, site)

# Some taxa are spp., but we think they're probably ok to include. The exception
# is Citrus spp. in Tucson, AZ (since they'd be heavily reliant on supplemental 
# water and phenology could vary among species too)
series <- series %>%
  filter(genus != "Citrus") %>%
  select(-c(genus, species))
df <- df %>%
  filter(genus != "Citrus") %>%
  select(-c(genus, species))

# Check that there's only one first yes date for each individual, year --------#
df %>%
  distinct(id, phenophase_description, first_yes_year) %>%
  nrow() ==
  nrow(df)
# This is the case.

# Add in Janet's phenophase and taxonomic groups ------------------------------#
categories <- df_orig %>%
  distinct(common_name, group, phenophase_description, phenophase)
series <- series %>%
  left_join(categories, by = c("common_name", "phenophase_description")) %>%
  select(-seriesid)

# Filter bird series ----------------------------------------------------------#
# See this report for more information about investigation of bird data series:
# https://erinzylstra.quarto.pub/exploration-of-bird-data-for-usgs-analysis/

# Load csv with information from ebird 
ebird <- read.csv("first-yes-years/birds_resident_presentJan1.csv")

# A little formatting:
ebird <- ebird %>%
  mutate(common_name = ifelse(common_name == "gray catbird", "grey catbird",
                              common_name)) %>%
  rename(site = site_id)

# Want to exclude series that indicate species presence when the species is
# expected to be present as of Jan 1
series_b <- series %>%
  left_join(select(ebird, common_name, site, phenophase_description, 
                   present_jan1),
            by = c("common_name", "site", "phenophase_description")) %>%
  mutate(present_jan1 = replace_na(present_jan1, 0)) %>%
  filter(!(present_jan1 == 1 & 
             phenophase_description %in% c("Calls or song (birds)",
                                           "Individuals at a feeding station",
                                           "Live individuals")))

# See whether same sites have individuals at feeding station and live
# individuals, since that would be redundant
station_redundant <- series_b %>%
  filter(group == "Bird") %>%
  group_by(common_name, site) %>%
  summarize(feeding = ifelse("Individuals at a feeding station" %in% phenophase_description, 1, 0),
            live = ifelse("Live individuals" %in% phenophase_description, 1, 0),
            .groups = "keep") %>%
  data.frame()
station_redundant_r <- station_redundant %>% 
  filter(feeding == 1 & live == 1) %>%
  select(common_name, site) %>%
  mutate(phenophase_description = "Individuals at a feeding station", 
         remove_stationr = 1)
# Remove redundancies 
series_b <- series_b %>%
  left_join(station_redundant_r, 
            by = c("common_name", "site", "phenophase_description")) %>%
  mutate(remove_stationr = replace_na(remove_stationr, 0)) %>%
  filter(remove_stationr == 0)

# See whether same sites have calls or song (different than singing individuals!)
# and live individuals, since that would also be redundant
calls_redundant <- series_b %>%
  filter(group == "Bird") %>%
  group_by(common_name, site) %>%
  summarize(calls = ifelse("Calls or song (birds)" %in% phenophase_description, 1, 0),
            live = ifelse("Live individuals" %in% phenophase_description, 1, 0),
            .groups = "keep") %>%
  data.frame()
calls_redundant_r <- calls_redundant %>% 
  filter(calls == 1 & live == 1) %>%
  select(common_name, site) %>%
  mutate(phenophase_description = "Calls or song (birds)", 
         remove_callsr = 1)
# Remove redundancies 
series_b <- series_b %>%
  left_join(calls_redundant_r, 
            by = c("common_name", "site", "phenophase_description")) %>%
  mutate(remove_callsr = replace_na(remove_callsr, 0)) %>%
  filter(remove_callsr == 0)

# Filter mammal series ----------------------------------------------------------#
# See this report for more information about investigation of mammal data series:
# https://erinzylstra.quarto.pub/exploration-of-mammal-data-for-usgs-analysis/

# Load csv with information about seasonal activity for each species-
# site combination.
mammal_traits <- read.csv("first-yes-years/mammal_species_traits.csv")

# A little formatting:
mammal_traits <- mammal_traits %>% rename(site = site_id)

# Append information about seasonal activity to series
series_mb <- series_b %>%
  left_join(select(mammal_traits, common_name, site, seasonally_inactive),
            by = c("common_name", "site")) %>%
  mutate(active_jan1 = ifelse(is.na(seasonally_inactive) | 
                                seasonally_inactive %in% c("hibernate", "torpor"),
                              0, 1)) %>%
  filter(!(active_jan1 == 1 &
             phenophase_description %in% c("Feeding", "Live individuals")))

# See whether same sites have feeding and live individuals, since that would be 
# redundant
feeding_redundant <- series_mb %>%
  filter(group == "Mammal") %>%
  group_by(common_name, site) %>%
  summarize(feeding = ifelse("Feeding" %in% phenophase_description, 1, 0),
            live = ifelse("Live individuals" %in% phenophase_description, 1, 0),
            .groups = "keep") %>%
  data.frame()
feeding_redundant_r <- feeding_redundant %>% 
  filter(feeding == 1 & live == 1) %>%
  select(common_name, site) %>%
  mutate(phenophase_description = "Feeding", 
         remove_feedingr = 1)
# Remove redundancies 
series_mb <- series_mb %>%
  left_join(feeding_redundant_r, 
            by = c("common_name", "site", "phenophase_description")) %>%
  mutate(remove_feedingr = replace_na(remove_feedingr, 0)) %>%
  filter(remove_feedingr == 0)

# Remove unnecessary columns
series_mb <- series_mb %>%
  select(-c(present_jan1, remove_stationr, remove_callsr,
            seasonally_inactive, active_jan1, remove_feedingr))

# Look at potential data filters (by series) ----------------------------------#

# Two potential data filters:
# Number of months in the series >= 6
# Length of time between first and last yes date for all yeartypes (calendar, 
# water, summer) >= 180 days

series_mb %>%
  filter(kingdom == "Plantae") %>%
  group_by(phenophase) %>%
  summarize(nseries = n(),
            nseries_good = sum(nmonths < 6 & min_gap < 180)) %>%
  mutate(prop_good = round(nseries_good/nseries, 2)) %>%
  data.frame()
# Only 68% of ripe fruit series meet both criteria; 94% of open flowers

series_mb %>%
  filter(kingdom == "Animalia") %>%
  group_by(group) %>%
  summarize(nseries = n(),
            nseries_good = sum(nmonths < 6 & min_gap < 180)) %>%
  mutate(prop_good = round(nseries_good/nseries, 2)) %>%
  data.frame()
# 83% of birds; 78% herps; 90% insects; 76% mammals

# Look at potential data filters by spp-php-states ----------------------------#
combos <- series_mb %>%
  group_by(kingdom, group, common_name, 
           phenophase_description, phenophase, state) %>%
  summarize(nseries = n(),
            nmonths_max = max(nmonths),
            gap_max = max(min_gap),
            .groups = "keep") %>%
  mutate(good = ifelse(nmonths_max < 6 & gap_max < 180, 1, 0)) %>%
  mutate(ngood = nseries * good) %>%
  data.frame()

combos %>%
  filter(kingdom == "Plantae") %>%
  group_by(phenophase) %>%
  summarize(nseries = sum(nseries),
            nseries_good = sum(ngood)) %>%
  mutate(prop_good = round(nseries_good/nseries, 2)) %>%
  data.frame()
# Then only keep 51% of ripe fruit series! and 84% of open flowers

combos %>%
  filter(kingdom == "Animalia") %>%
  group_by(group) %>%
  summarize(nseries = sum(nseries),
            nseries_good = sum(ngood)) %>%
  mutate(prop_good = round(nseries_good/nseries, 2)) %>%
  data.frame()
# Then only keep 78% of birds, 75% herps, 90% insects, 76% mammals

# Look at a few species with known issues -------------------------------------#

# American witchhazel
filter(series_mb, common_name == "American witchhazel") %>%
  select(common_name, id, site, state, phenophase_description, nyrs, nmonths,
         calyr_gap, wateryr_gap, summeryr_gap, min_gap, nquarters) %>%
  filter(nmonths > 5 | min_gap >= 180) %>%
  arrange(phenophase_description)

# # Look at first yes months for each phenophase:
# filter(df, common_name == "American witchhazel") %>%
#   mutate(first_yes_month = factor(first_yes_month)) %>%
#   ggplot(aes(x = first_yes_month)) +
#   geom_bar() +
#   facet_grid(phenophase_description ~ state, scales = "free_y")
# # Things look fine for BLB, colored leaves, falling leaves
# 
# # Open flowers
# df %>%
#   filter(common_name == "American witchhazel",
#          phenophase_description == "Open flowers") %>%
#   mutate(first_yes_month = factor(first_yes_month,
#                                   levels = 1:12)) %>%
#   ggplot(aes(x = first_yes_month)) +
#   geom_bar() +
#   facet_grid(id ~ state) +
#   scale_x_discrete(drop = FALSE) +
#   labs(x = "First yes month", y = "Count",
#        title = "American witchhazel, Open flowers")
# # Open flowers in ME and NY spread out a lot
# 
# # Fruit phenopohases
# df %>%
#   mutate(idstate = paste(id, state, sep = "_")) %>%
#   filter(common_name == "American witchhazel",
#          phenophase_description %in% c("Recent fruit or seed drop",
#                                        "Ripe fruits")) %>%
#   mutate(first_yes_month = factor(first_yes_month,
#                                   levels = 1:12)) %>%
#   ggplot(aes(x = first_yes_month)) +
#   geom_bar() +
#   geom_vline(xintercept = 9.5, linetype = "dotted", color = "steelblue3") +
#   facet_grid(idstate ~ phenophase_description, scales = "free_y") +
#   scale_x_discrete(drop = FALSE) +
#   labs(x = "First yes month", y = "Count",
#        title = "American witchhazel, fruit phenophases")
# # Fruit phenophases in ME look ok with calendar year
# # Fruit phenophases in NY spread out a lot

# Tuliptree
filter(series_mb, common_name == "tuliptree") %>%
  select(common_name, id, site, state, phenophase_description, nyrs, nmonths,
         calyr_gap, wateryr_gap, summeryr_gap, min_gap, nquarters) %>%
  filter(nmonths > 5 | min_gap >= 180) %>%
  arrange(phenophase_description)

# Sweetgum 
filter(series_mb, common_name == "sweetgum") %>%
  select(common_name, id, site, state, phenophase_description, nyrs, nmonths,
         calyr_gap, wateryr_gap, summeryr_gap, min_gap, nquarters) %>%
  filter(nmonths > 5 | min_gap >= 180) %>%
  arrange(phenophase_description)

# Coyotebrush
filter(series_mb, common_name == "coyotebrush") %>%
  select(common_name, id, site, state, phenophase_description, nyrs, nmonths,
         calyr_gap, wateryr_gap, summeryr_gap, min_gap, nquarters) %>%
  filter(nmonths > 5 | min_gap >= 180) %>%
  arrange(phenophase_description)

# Jojoba
filter(series_mb, common_name == "jojoba") %>%
  select(common_name, id, site, state, phenophase_description, nyrs, nmonths,
         calyr_gap, wateryr_gap, summeryr_gap, min_gap, nquarters) %>%
  filter(nmonths > 5 | min_gap >= 180) %>%
  arrange(phenophase_description)

# Oysterwood
filter(series_mb, common_name == "oysterwood") %>%
  select(common_name, id, site, state, phenophase_description, nyrs, nmonths,
         calyr_gap, wateryr_gap, summeryr_gap, min_gap, nquarters) %>%
  filter(nmonths > 5 | min_gap >= 180) %>%
  arrange(phenophase_description)

# Use months/gap criteria to remove individual series -------------------------#
seriesf <- series_mb %>%
  filter(nmonths < 6 & min_gap < 180)

# Group by species-phenoophase-state, pick yeartype ---------------------------#
sps <- seriesf %>%
  group_by(kingdom, common_name, group, phenophase_description, phenophase,
           state) %>%
  summarize(nseries = n(),
            calyr_180 = sum(calyr_gap < 180),
            wateryr_180 = sum(wateryr_gap < 180),
            summeryr_180 = sum(summeryr_gap < 180),
            calyr_md = median(calyr_gap),
            wateryr_md = median(wateryr_gap),
            summeryr_md = median(summeryr_gap),
            .groups = "keep") %>%
  mutate(calyr_p = calyr_180/nseries,
         wateryr_p = wateryr_180/nseries,
         summeryr_p = summeryr_180/nseries) %>%
  data.frame()

count(sps, calyr_p == 1, wateryr_p == 1, summeryr_p ==1)
# There are only 4 spp-php-state combinations where one of the yeartypes doesn't 
# work for all series

yearissues <- sps %>%
  filter(calyr_p < 1 & wateryr_p < 1 & summeryr_p < 1)
yearissues
# They are all in CA

# Checked, but there weren't any clear latitudinal patterns
# seriesf %>%
#   filter(common_name == yearissues$common_name[4] &
#            phenophase_description == yearissues$phenophase_description[4] &
#            state == yearissues$state[4]) %>%
#   select(common_name, id, site, state, phenophase_description, nyrs, nmonths,
#          calyr_gap, wateryr_gap, summeryr_gap) %>%
#   left_join(distinct(df, site, lat, lon), by = "site")

# For these, will pick the yeartype that has the highest proportion of series 
# with gaps of <180 days. If there's a tie, then select the yeartypethat has the 
# lowest median gap length.
yeartypes <- c("calendar", "water", "summer")
sps <- sps %>%
  rowwise() %>%
  mutate(md = yeartypes[order(c_across(calyr_md:summeryr_md), decreasing = FALSE)[1]]) %>%
  ungroup() %>%
  mutate(yeartype = case_when(
    calyr_p == 1 ~ "calendar",
    wateryr_p == 1 ~ "water",
    summeryr_p == 1 ~ "summer",
    calyr_p > wateryr_p & calyr_p > summeryr_p ~ "calendar",
    wateryr_p > calyr_p & wateryr_p > summeryr_p ~ "water",
    summeryr_p > calyr_p & summeryr_p > wateryr_p ~ "summer",
    .default = md
  )) %>%
  data.frame()

# Create list of series to include --------------------------------------------#
series_dl <- seriesf %>%
  select(kingdom, common_name, id, site, state, phenophase_description, 
         nyrs, nmonths) %>%
  # Add in yeartype
  left_join(select(sps, common_name, phenophase_description, state, yeartype),
            by = c("common_name", "phenophase_description", "state")) %>%
  # Add in phenophase_id, species ID, and lat/lon
  left_join(distinct(df, common_name, phenophase_description, phenophase_id, 
                     species_id, site, lat, lon),
            by = c("common_name", "phenophase_description", "site"))

# Check calendar year data ----------------------------------------------------#
# jojoba <- df %>%
#   filter(id == 117061, phenophase_description == "Open flowers") %>%
#   select(site, species_id, phenophase_id, first_yes_year, first_yes_doy, 
#          numdays_since_prior_no)
# jojoba
# 
# jojoba_dl <- npn_download_individual_phenometrics(
#   request_source = "erinz",
#   years = 2009:2024,
#   individual_ids = 117061,
#   phenophase_ids = 501
# ) %>% data.frame()
# 
# jojoba_dl2 <- jojoba_dl %>%
#   select(first_yes_year, first_yes_doy, numdays_since_prior_no) %>%
#   arrange(first_yes_year, first_yes_doy) %>%
#   distinct(first_yes_year, .keep_all = TRUE)
# 
# # Check if they're the same:
# all.equal(select(jojoba, first_yes_year, first_yes_doy, numdays_since_prior_no),
#           jojoba_dl2)
# # Yes

# So, we shouldn't have to re-download any of the calendar year series

# Download water year data ----------------------------------------------------#
wy <- filter(series_dl, yeartype == "water")

# Could grab data for each series based on individual_id and phenophase_id. 
# However, the number of years in the time series could increase or decrease 
# relative to a time series based on the calendar year. For that reason, it 
# might be worth downloading all data available for that species (and 
# phenophase) in the state in case we pick up new series.

wy_sppst <- wy %>%
  distinct(species_id, state, phenophase_id)

# Can use years = 2009:2024 for all (first water yr = Oct 2009 - Sep 2010)
# Note that water year 2024-2025 isn't complete yet, but that won't matter since
# we're very close to the end and only use water year if most observations are 
# many months before that.

# Can skip climate data (using PRISM now)

# for (i in 1:nrow(wy_sppst)) {
#   dl_temp <- npn_download_individual_phenometrics(
#     request_source = "erinz",
#     years = 2009:2024,
#     species_ids = wy_sppst$species_id[i],
#     states = wy_sppst$state[i],
#     phenophase_ids = wy_sppst$phenophase_id[i],
#     period_start = "10-01",
#     period_end = "09-30",
#     additional_fields = "site_name"
#   )
#   if (i == 1) {
#     wy_dl <- dl_temp
#   } else {
#     wy_dl <- rbind(wy_dl, dl_temp)
#   }
# }
# 
# wy_dl <- data.frame(wy_dl)
# write.csv(wy_dl, "first-yes-years/wateryr_download.csv", row.names = FALSE)

wy_dl <- read.csv("first-yes-years/wateryr_download.csv")

# Important note:
# First yes year, doy are related to the CALENDAR year, not the water year
# If we want to select the first yes in a water year, need to create that 
# variable first.
# water year = first_yes_year if month %in% 10:12, first_yes_year - 1 otherwise
# This is different than how it is typically done 
# (see https://water.usgs.gov/nwc/explain_data.html), but using this because
# then both water and summer year will be associated with the calendar year
# in which they start.
# Create dowy variable using function created at start of script (after creating 
# observation date)

wy_dl2 <- wy_dl %>%
  # Calculate water year
  mutate(water_year = ifelse(first_yes_month %in% 10:12,
                             first_yes_year,
                             first_yes_year - 1)) %>%
  # Calculate yes date
  mutate(yes_date = parse_date_time(x = paste(first_yes_year, first_yes_doy),
                                    orders = "yj")) %>%
  # Calculate day of water year (dowy)
  mutate(dowy = wateryr_calc(yes_date)) %>%
  # Remove all but the first yes in each water year
  arrange(individual_id, phenophase_id, water_year, dowy) %>%
  distinct(individual_id, phenophase_id, water_year, .keep_all = TRUE) %>%
  data.frame()

# Summarize information for each series
wy_dl_series <- wy_dl2 %>%
  rename(id = individual_id, 
         site = site_id) %>%
  group_by(common_name, genus, species, kingdom, id, site, site_name, state,
           phenophase_description) %>%
  summarize(nwyrs = n(),
            # Number of unique months in series
            wy_nmonths = n_distinct(first_yes_month),
            # Calculate length of time between min/max dowy
            wateryr_gap = max(dowy) - min(dowy),
            .groups = "keep") %>%
  data.frame()

# Want to see whether series we were looking for are in this new dataset:
wy <- wy %>%
  left_join(select(wy_dl_series, id, phenophase_description, nwyrs, wy_nmonths, 
                   wateryr_gap),
            by = c("id", "phenophase_description"))
# Plants at NEON sites are problematic because the site and plant IDs change
# annually (so we can't match them up). Will assess these separately.

# For non-NEON sites (series that we could match up with IDs):
wy %>% 
  filter(!is.na(nwyrs)) %>%
  count(wy_nmonths > 5, wateryr_gap >= 180)
# 47 of 67 met criteria (most failed because there were > 180 days between
# earliest and latest date)

# For NEON sites:
filter(wy, is.na(nwyrs)) %>% nrow() # 10 series in Janet's dataset
wy_dl_series %>%
  filter(grepl(".phenology.", site_name)) %>%
  filter(nwyrs > 8 & wy_nmonths < 6 & wateryr_gap < 180) %>%
  nrow() 
# 10 series in new dataset

# What series are they?
wy %>%
  filter(is.na(nwyrs)) %>%
  select(common_name, state, phenophase_description, nyrs, nmonths)
wy_dl_series %>%
  filter(grepl(".phenology.", site_name)) %>%
  filter(nwyrs > 8 & wy_nmonths < 6 & wateryr_gap < 180) %>%
  select(common_name, site_name, state, phenophase_description, nwyrs,
         wy_nmonths, wateryr_gap)
# More longleaf pine in GA and monkeypod in PR

wy_dl_series %>%
  filter(grepl(".phenology.", site_name)) %>%
  filter(common_name == "creosote bush" & nwyrs > 8) %>%
  select(common_name, site_name, state, id, phenophase_description, nwyrs,
         wy_nmonths, wateryr_gap)
# No series with creosote in AZ (all have 3-5 months of first yeses that
# span > 270 water days)

# Are we gaining some new (non-NEON) series?
wy_dl_series %>%
  filter(!grepl(".phenology.", site_name)) %>%
  filter(nwyrs > 8 & wy_nmonths < 6 & wateryr_gap < 180) %>%
  select(common_name, site_name, state, id, phenophase_description, nwyrs,
         wy_nmonths, wateryr_gap) %>%
  nrow() # 48
# So just gaining one new series.

# Identify wateryr series that meet the filtering criteria (nmonths < 6,
# wateryr_gap < 180)
wy_dl_series <- wy_dl_series %>%
  mutate(keep = ifelse(nwyrs > 8 & wy_nmonths < 6 & wateryr_gap < 180, 1, 0))
# Total of 58

# Extract data for those series from wy_dl2
wy_dl_filtered <- wy_dl2 %>%
  left_join(select(wy_dl_series, id, phenophase_description, keep),
            by = c("individual_id" = "id", 
                   "phenophase_description" = "phenophase_description")) %>%
  filter(keep == 1) %>%
  select(-keep)

# Download summer year data ---------------------------------------------------#
sy <- filter(series_dl, yeartype == "summer")

sy_sppst <- sy %>%
  distinct(species_id, state, phenophase_id)

# for (i in 1:nrow(sy_sppst)) {
#   dl_temp <- npn_download_individual_phenometrics(
#     request_source = "erinz",
#     years = 2009:2024,
#     species_ids = sy_sppst$species_id[i],
#     states = sy_sppst$state[i],
#     phenophase_ids = sy_sppst$phenophase_id[i],
#     period_start = "07-01",
#     period_end = "06-30",
#     additional_fields = "site_name"
#   )
#   if (i == 1) {
#     sy_dl <- dl_temp
#   } else {
#     sy_dl <- rbind(sy_dl, dl_temp)
#   }
# }
# 
# sy_dl <- data.frame(sy_dl)
# write.csv(sy_dl, "first-yes-years/summeryr_download.csv", row.names = FALSE)

sy_dl <- read.csv("first-yes-years/summeryr_download.csv")

# Important note:
# First yes year, doy are related to the CALENDAR year, not the summer year
# If we want to select the first yes in a summer year, need to create that 
# variable first.
# summer year = first_yes_year if month %in% 7:12, first_yes_year - 1 otherwise
# Create dosy variable using function created at start of script (after creating 
# observation date)

sy_dl2 <- sy_dl %>%
  # Calculate water year
  mutate(summer_year = ifelse(first_yes_month %in% 7:12,
                              first_yes_year,
                              first_yes_year - 1)) %>%
  # Calculate yes date
  mutate(yes_date = parse_date_time(x = paste(first_yes_year, first_yes_doy),
                                    orders = "yj")) %>%
  # Calculate day of water year (dowy)
  mutate(dosy = summeryr_calc(yes_date)) %>%
  # Remove all but the first yes in each water year
  arrange(individual_id, phenophase_id, summer_year, dosy) %>%
  distinct(individual_id, phenophase_id, summer_year, .keep_all = TRUE) %>%
  data.frame()

# Summarize information for each series
sy_dl_series <- sy_dl2 %>%
  rename(id = individual_id, 
         site = site_id) %>%
  group_by(common_name, genus, species, kingdom, id, site, site_name, state,
           phenophase_description) %>%
  summarize(nsyrs = n(),
            # Number of unique months in series
            sy_nmonths = n_distinct(first_yes_month),
            # Calculate length of time between min/max dowy
            summeryr_gap = max(dosy) - min(dosy),
            .groups = "keep") %>%
  data.frame()

# Want to see whether series we were looking for are in this new dataset:
sy <- sy %>%
  left_join(select(sy_dl_series, id, phenophase_description, nsyrs, sy_nmonths, 
                   summeryr_gap),
            by = c("id", "phenophase_description"))
sum(is.na(sy$summeryr_gap))
# Doesn't look like we have any plants at NEON sites 

# For non-NEON sites (series that we could match up with IDs):
sy %>% 
  count(sy_nmonths > 5, summeryr_gap >= 180)
# 167 of 182 series meet filtering criteria

# Are we gaining some new series?
sy_dl_series %>%
  filter(nsyrs > 8 & sy_nmonths < 6 & summeryr_gap < 180) %>%
  select(common_name, site_name, state, id, phenophase_description, nsyrs,
         sy_nmonths, summeryr_gap) %>%
  nrow() # 225
# So just gaining 58 series (225 - 167)

# Identify summeryr series that meet the filtering criteria (nmonths < 6,
# wateryr_gap < 180)
sy_dl_series <- sy_dl_series %>%
  mutate(keep = ifelse(nsyrs > 8 & sy_nmonths < 6 & summeryr_gap < 180, 1, 0))
# Total of 225

# Extract data for those series from sy_dl2
sy_dl_filtered <- sy_dl2 %>%
  left_join(select(sy_dl_series, id, phenophase_description, keep),
            by = c("individual_id" = "id", 
                   "phenophase_description" = "phenophase_description")) %>%
  filter(keep == 1) %>%
  select(-keep)

# Extract calendar year series data we're keeping -----------------------------#
cy <- filter(series_dl, yeartype == "calendar")

# Formatting, so we can extract data from original file
cy <- cy %>%
  mutate(keep = 1) %>%
  rename(individual_id = id)

# Extract data from original dataframe
cy_filtered <- df_orig %>%
  left_join(select(cy, individual_id, phenophase_description, keep),
            by = c("individual_id", "phenophase_description")) %>%
  mutate(keep = replace_na(keep, 0)) %>%
  filter(keep == 1) %>%
  select(-c(keep, gdd, gddf, 
            tmax_winter, tmax_spring, tmax_summer, tmax_fall, tmax, tmaxf,
            tmin_winter, tmin_spring, tmin_summer, tmin_fall, tmin, tminf,
            prcp_winter, prcp_spring, prcp_summer, prcp_fall, prcp,
            acc_prcp, daylength, 
            ind_phen, ind_phen_year, MAT, phenophase, group)) %>%
  mutate(yes_date = parse_date_time(x = paste(first_yes_year, first_yes_doy),
                                    orders = "yj")) %>%
  mutate(yeartype = "calendar",
         water_year = NA,
         dowy = NA,
         summer_year = NA,
         dosy = NA)

# Merge all the datasets together ---------------------------------------------#
wy_filtered <- wy_dl_filtered %>%
  mutate(yeartype = "water", .after = yes_date) %>%
  relocate(water_year, .after = yeartype) %>%
  mutate(summer_year = NA,
         dosy = NA) %>%
  select(-site_name)

sy_filtered <- sy_dl_filtered %>%
  mutate(yeartype = "summer", .after = yes_date) %>%
  relocate(summer_year, .after = yeartype) %>%
  mutate(water_year = NA, .after = yeartype) %>%
  mutate(dowy = NA, .after = water_year) %>%
  select(-site_name)

alldata <- rbind(cy_filtered, wy_filtered, sy_filtered)

# Add in Janet's phenophase categories
alldata <- alldata %>%
  left_join(select(categories, common_name, phenophase_description, phenophase),
            by = c("common_name", "phenophase_description"))

# Add in functional groups, but get data from rnpn database
rnpn_spp <- npn_species() %>%
  select(kingdom, common_name, functional_type) %>%
  distinct() %>%
  data.frame()

alldata <- alldata %>%
  left_join(rnpn_spp, by = c("kingdom", "common_name"))

# Create columns with year, day of year appropriate for each series:
alldata <- alldata %>%
  mutate(first_yes = case_when(
    yeartype == "calendar" ~ first_yes_doy,
    yeartype == "water" ~ dowy,
    yeartype == "summer" ~ dosy
  )) %>%
  mutate(year = case_when(
    yeartype == "calendar" ~ first_yes_year,
    yeartype == "water" ~ water_year,
    yeartype == "summer" ~ summer_year
  ))

# Check that all series meet filtering criteria
alldatas <- alldata %>%
  group_by(kingdom, functional_type, common_name, individual_id, 
           phenophase_description, phenophase, yeartype) %>%
  summarize(nyrs = n_distinct(year),
            nmonths = n_distinct(first_yes_month),
            gap = max(first_yes) - min(first_yes),
            .groups = "keep") %>%
  data.frame()

filter(alldatas, nmonths > 5 | gap >= 180) 
# 4 series -- they're all species-state-php combinations where one yeartype 
# didn't fit all series well. Will remove these 4 series from final dataset

alldatas <- alldatas %>%
  mutate(remove = ifelse(nmonths > 5 | gap >= 180, 1, 0))
alldata <- alldata %>%
  left_join(select(alldatas, individual_id, phenophase_description, remove),
            by = c("individual_id", "phenophase_description")) %>%
  filter(remove == 0) %>%
  select(-remove)
alldatas <- alldatas %>%
  filter(remove == 0) %>%
  select(-remove)

# Write final dataset to file:
# write.csv(alldata, "first-yes-years/final-dataset.csv", row.names = FALSE)

# Summarize info on summer/water year series and save to file:
# alldata %>%
#   filter(yeartype != "calendar") %>%
#   group_by(yeartype, functional_type, common_name, phenophase_id, 
#            phenophase_description, state, site_id, latitude, longitude,
#            elevation_in_meters) %>% 
#   summarize(earliest_yes_DOY = min(first_yes_doy),
#             lastest_yes_DOY = max(first_yes_doy),
#             .groups = "keep") %>%
#   data.frame() %>%
#   arrange(functional_type, common_name, state, phenophase_description,
#           .locale = "en") %>%
#   write.csv("first-yes-years/final-dataset-water-summer-series.csv",
#             row.names = FALSE)

# Look at a few things with the new dataset.....................................

# Original number of plant series, by phenophase:
df_orig %>%
  filter(kingdom == "Plantae") %>%
  distinct(individual_id, phenophase_description, phenophase) %>%
  group_by(phenophase) %>%
  summarize(nseries = n())
# New number of plant series, by phenophase:
alldata %>%
  filter(kingdom == "Plantae") %>%
  distinct(individual_id, phenophase_description, phenophase) %>%
  group_by(phenophase) %>%
  summarize(nseries = n())

# How do things look for one of the problematic series:
filter(alldata, common_name == "tuliptree", 
       phenophase_description == "Ripe fruits") %>%
  ggplot(aes(x = year, y = first_yes)) + 
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  facet_wrap(~individual_id)
filter(alldata, common_name == "tuliptree", 
       phenophase_description == "Recent fruit or seed drop") %>%
  ggplot(aes(x = year, y = first_yes)) + 
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  facet_wrap(~individual_id)

filter(alldata, common_name == "sweetgum", 
       phenophase_description == "Open flowers") %>%
  ggplot(aes(x = year, y = first_yes)) + 
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  facet_wrap(~individual_id)

filter(alldata, common_name == "coyotebrush", 
       phenophase_description == "Open flowers") %>%
  ggplot(aes(x = year, y = first_yes)) + 
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  facet_wrap(~individual_id)
filter(alldata, common_name == "coyotebrush", 
       phenophase_description == "Ripe fruits") %>%
  ggplot(aes(x = year, y = first_yes)) + 
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  facet_wrap(~individual_id)
filter(alldata, common_name == "coyotebrush", 
       phenophase_description == "Recent fruit or seed drop") %>%
  ggplot(aes(x = year, y = first_yes)) + 
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  facet_wrap(~individual_id)
