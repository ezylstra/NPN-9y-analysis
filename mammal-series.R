# Exploration of mammal data in USGS analysis
# ER Zylstra
# 20 March 2025

library(dplyr)
library(stringr)
library(lubridate)
library(terra)
library(tidyterra)
library(ggplot2)
library(cowplot)

rm(list = ls())

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

# Categorizing series by phenophase type --------------------------------------#

count(df_mammals, phenophase_description) %>% arrange(desc(n))
# Live individuals       
# Feeding
# Young individuals
# Fruit/seed consumption
# Nut gathering
# Male combat (horn/antler)
# Males vocalizing

# Live individuals & Young individuals indicate PRESENCE (with young provding
# more specific information that is likely to be seasonal)
# Male combat & Males vocalizing indicate a BEHAVIOR
# Feeding, Fruit/seed consumption, and Nut gathering indicate FEEDING

# Add phenophase_type label (making these explicit for easy changes)
df_mammals <- df_mammals %>%
  mutate(phenophase_type = case_when(
    phenophase_description == "Live individuals" ~ "presence",
    phenophase_description == "Young individuals" ~ "presence",
    phenophase_description == "Male combat (horn/antler)" ~ "behavior",
    phenophase_description == "Males vocalizing" ~ "behavior",
    phenophase_description == "Feeding" ~ "feeding",
    phenophase_description == "Fruit/seed consumption" ~ "feeding",
    phenophase_description == "Nut gathering" ~ "feeding",
    .default = NA
  ))

# Suggestions for including/excluding series in analyses ----------------------#

# If phenophase is any behavior: include
# If phenophase indicates feeding: include after removing any redundancies
# If phenophase is "Young individuals": include
# If phenophase is "Live individuals": include if species is inactive at start
# of year (or is migratory, though I'm not sure any of these mammals are)


# Need to figure out:
# hibernates (ie, inactive at start of year at that location)?

filter(df_mammals, phenophase_description == "Young individuals") %>% 
  arrange(common_name, latitude)
# when young are present

# get seasonal range maps from MOL? (not sure if this is necessary)



