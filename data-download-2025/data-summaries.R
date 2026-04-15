# Summarize data (using 14-day prior no subset)
# 15 Apr 2026

library(dplyr)
library(tidyr)
library(lubridate)
library(stringr)

# Load data
birds <- read.csv("data-download-2025/data-out/bird-series-prior14-thru2025.csv")
herps <- read.csv("data-download-2025/data-out/herp-series-prior14-thru2025.csv")
insects <- read.csv("data-download-2025/data-out/insect-series-prior14-thru2025.csv")
mammals <- read.csv("data-download-2025/data-out/mammal-series-prior14-thru2025.csv")
plants <- read.csv("data-download-2025/data-out/plant-series-prior14-thru2025.csv")

# List of phenophases ---------------------------------------------------------#
php_birds <- birds %>%
  group_by(phenophase_description) %>%
  summarize(n_series = n_distinct(ind_phen)) %>%
  data.frame() %>%
  mutate(kingdom = "Animalia", .before = phenophase_description) %>%
  mutate(lifeform_phenophase = "Bird", .after = "kingdom")
php_herps <- herps %>%
  group_by(phenophase_description) %>%
  summarize(n_series = n_distinct(ind_phen)) %>%
  data.frame() %>%
  mutate(kingdom = "Animalia", .before = phenophase_description) %>%
  mutate(lifeform_phenophase = "Herpetofauna", .after = "kingdom")
php_insects <- insects %>%
  group_by(phenophase_description) %>%
  summarize(n_series = n_distinct(ind_phen)) %>%
  data.frame() %>%
  mutate(kingdom = "Animalia", .before = phenophase_description) %>%
  mutate(lifeform_phenophase = "Insect", .after = "kingdom")
php_mammals <- mammals %>%
  group_by(phenophase_description) %>%
  summarize(n_series = n_distinct(ind_phen)) %>%
  data.frame() %>%
  mutate(kingdom = "Animalia", .before = phenophase_description) %>%
  mutate(lifeform_phenophase = "Mammal", .after = "kingdom")
php_plants <- plants %>%
  group_by(phenophase, phenophase_description) %>%
  summarize(n_series = n_distinct(ind_phen), .groups = "drop") %>%
  data.frame() %>%
  rename(lifeform_phenophase = phenophase) %>%
  mutate(kingdom = "Plantae", .before = lifeform_phenophase)

# Merge all for Table S1:
phenophases <- rbind(php_birds,
                     php_herps,
                     php_insects,
                     php_mammals,
                     php_plants)
# Write to file:
write.csv(phenophases, 
          "data-download-2025/data-out/tableS1-phenophases.csv",
          row.names = FALSE)

# Year type designations ------------------------------------------------------#

yt_mammals <- mammals %>%
  filter(yeartype != "calendar") %>%
  group_by(yeartype, common_name, phenophase_description, state) %>%
  summarize(n_series = n_distinct(ind_phen), .groups = "drop") %>%
  data.frame() %>%
  mutate(group = "Mammal", .before = yeartype)
yt_plants <- plants %>%
  filter(yeartype != "calendar") %>%
  group_by(yeartype, common_name, phenophase_description, state) %>%
  summarize(n_series = n_distinct(ind_phen), .groups = "drop") %>%
  data.frame() %>%
  mutate(group = "Plant", .before = yeartype) %>%
  arrange(common_name, phenophase_description, .locale = "en")

# Merge all for Table S2:
yeartypes <- rbind(yt_mammals, yt_plants)
# Write to file
write.csv(yeartypes, 
          "data-download-2025/data-out/tableS2-yeartypes.csv",
          row.names = FALSE)
