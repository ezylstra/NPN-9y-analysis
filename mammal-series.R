# Exploration of mammal data in USGS analysis
# ER Zylstra
# 24 March 2025

library(rnpn)
library(dplyr)
library(stringr)
library(tidyr)

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

# Remove "Dead individuals" phenophase
df <- df %>%
  filter(phenophase_description != "Dead individuals")

# Exclude birds
df_mammals <- df %>%
  filter(group == "Mammal")

# Summarize by species
df_mammals %>%
  group_by(common_name) %>%
  summarize(n_series = n(),
            n_phenophases = n_distinct(phenophase_description),
            n_sites = n_distinct(site_id)) %>%
  arrange(desc(n_series), desc(n_sites)) %>%
  data.frame()

# Obtain information about seasonal activity ----------------------------------#

# Unfortunately, there doesn't seem to be a database for mammal species that has 
# information about seasonal activity and/or hibernation. (Looked at PanTHERIA,
# MammalBase, and IUCN). Extracting information about seasonal activity from:
  # Mammals of North America - Peterson Field Guide 2006 edition
  # Animal Diversity Web (University of Michigan; animaldiversity.org)
# Additional input from E. Posthumus

# Write species-site combinations to file:
mammals_file <- "data/mammal_species.csv"
if (!file.exists(mammals_file)) {
  spp_sites <- df_mammals %>%
    select(common_name, site_id, latitude, longitude) %>%
    distinct() %>%
    arrange(common_name, site_id)
  npn_spp <- npn_species() %>%
    filter(common_name %in% unique(spp_sites$common_name)) %>%
    data.frame()
  spp_sites <- spp_sites %>%
    left_join(select(npn_spp, common_name, genus, species), by = "common_name")
  write.csv(spp_sites, mammals_file, row.names = FALSE)
}

# Created a new csv with information about seasonal activity for each species-
# site combination. Load file:
mammal_traits <- read.csv("data/mammal_species_traits.csv")

# Append information about hiberation/winter torpor
df_mammals <- df_mammals %>%
  left_join(select(mammal_traits, common_name, site_id, seasonally_inactive),
            by = c("common_name", "site_id")) %>%
  mutate(active_jan1 = ifelse(seasonally_inactive == "no", 1, 0))

# Categorizing series by phenophase type --------------------------------------#

count(df_mammals, phenophase_description, common_name)
# Live individuals       
# Feeding
# Young individuals
# Fruit/seed consumption
# Nut gathering
# Male combat (horn/antler)
# Males vocalizing

# Live individuals indicates the species is present and active: ACTIVE

# Young individuals also indicates PRESENCE, but it's likely to be seasonal for
# most species, so maybe best to distinguish by labeling as YOUNG

# Male combat & Males vocalizing indicate a BEHAVIOR

# Fruit/seed consumption and Nut gathering indicate FEEDING on resources that 
# may be seasonal (and thus interesting to think about phenological changes).
# Feeding is more general, and it's not clear that there would be a seasonal 
# signal for most mammals. If they're present and active, they'll be feeding 
# all the time...

# Add phenophase_type label (making these explicit for easy changes)
df_mammals <- df_mammals %>%
  mutate(phenophase_type = case_when(
    phenophase_description == "Live individuals" ~ "activity",
    phenophase_description == "Young individuals" ~ "young",
    phenophase_description == "Male combat (horn/antler)" ~ "behavior",
    phenophase_description == "Males vocalizing" ~ "behavior",
    phenophase_description == "Feeding" ~ "feeding",
    phenophase_description == "Fruit/seed consumption" ~ "feeding",
    phenophase_description == "Nut gathering" ~ "feeding",
    .default = NA
  ))

# Suggestions for including/excluding series in analyses ----------------------#

# If phenophase is a behavior: include
# If phenophase is related to feeding on some particular resource: include
# If phenophase indicates presence of young: include
# If phenophase indicates presence/activity or feeding (on anything): include 
  # species if it is not usually active at the start of the year

phpt <- df_mammals %>%
  count(phenophase_description, phenophase_type) %>%
  mutate(rule_word = case_when(
    phenophase_type == "activity" ~ "Include if inactive on 1 Jan",
    phenophase_description == "Feeding" ~ "Include if inactive on 1 Jan",
    .default = "Include"
  )) %>%
  arrange(factor(phenophase_type, 
                 levels = c("young", "active", "feeding", "behavior")))

df_mammals <- df_mammals %>%
  mutate(include = case_when(
    phenophase_type %in% c("behavior", "young") ~ 1,
    phenophase_description %in% c("Fruit/seed consumption", "Nut gathering") ~ 1,
    active_jan1 == 0 ~ 1,
    .default = 0
  ))

# checks:
# count(df_mammals, phenophase_description, active_jan1, include)
# count(df_mammals, include)

include <- filter(df_mammals, include == 1)

# Looking for series with redundant information -------------------------------#

# Indications of species activity and general feeding observations:
activity <- include %>%
  filter(phenophase_description %in% c("Feeding", "Live individuals")) %>%
  group_by(common_name, site_id) %>%
  summarize(live = ifelse("Live individuals" %in% phenophase_description, 1, 0),
            feeding = ifelse("Feeding" %in% phenophase_description, 1, 0),
            .groups = "keep") %>%
  data.frame()

activity_count <- count(activity, live, feeding)
# No species-sites where there is a feeding series but no live individual series

# Remove all "Feeding" series since they're always accompanied by 
# "Live individuals" series
include <- include %>%
  filter(phenophase_description != "Feeding")

# Summarizing series that are left --------------------------------------------#

# Create "include2" variable to indicate whether a series should be included in 
# analyses after accounting for all factors
include <- include %>%
  rename(include2 = include)
df_mammals <- df_mammals %>%
  left_join(select(include, common_name, site_id, phenophase_description, include2),
            by = c("common_name", "site_id", "phenophase_description")) %>%
  mutate(include2 = replace_na(include2, replace = 0))

# Summarize by species
spp_include <- include %>%
  group_by(common_name, seasonally_inactive) %>%
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
  arrange(factor(phenophase_type, 
                 levels = c("young", "active", "feeding", "behavior")), 
          desc(n_series)) %>%
  data.frame()
