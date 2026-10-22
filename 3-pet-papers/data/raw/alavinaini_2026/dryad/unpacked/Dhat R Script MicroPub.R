# Bobcat and Domestic Cat Temporal Overlap Analysis
#
# This script calculates temporal activity overlap between bobcats and domestic cats
# using solar-time transformed detection timestamps and the coefficient of overlap (Dhat).
#
# Required input files:
#   1. Bobcat_Data.csv
#   2. Domestic_Cat_Data.csv
#
# Instructions:
#   Place this script in the same folder as the two CSV files.
#   In RStudio, set the working directory to the source file location.
#   Run the script from top to bottom.
#
# Outputs:
#   - Dhat estimates for solo-site and shared-site comparisons
#   - overlap plots
#   - bootstrap confidence intervals

# Install these packages once, if needed:
# install.packages(c("overlap", "lubridate", "dplyr", "readr", "suntools"))
library(overlap)
library(lubridate)
library(dplyr)
library(readr)
library(suntools)

# Before running this script, place it in the same folder as Bobcat_Data.csv and Domestic_Cat_Data.csv.
# In RStudio, set the working directory to the script folder:
# Session > Set Working Directory > To Source File Location
bobcats <- read_csv("Bobcat_Data.csv", col_types = cols(.default = col_character()))
domestic_cats <- read_csv("Domestic_Cat_Data.csv", col_types = cols(.default = col_character()))
bobcats <- bobcats %>%
  mutate(species_group = "bobcat")

domestic_cats <- domestic_cats %>%
  mutate(species_group = "domestic_cat")
all_detections <- bind_rows(bobcats, domestic_cats)
all_detections <- all_detections %>%
  rename(
    timestamp_utc = `Timestamp (UTC)`,
    site_status = Status,
    site = SiteAbbr,
    latitude = Latitude,
    longitude = Longitude
  )
all_detections <- all_detections %>%
  mutate(
    site_status = tolower(site_status),
    timestamp_utc = ymd_hms(timestamp_utc, tz = "UTC"),
    timestamp_local = with_tz(timestamp_utc, tzone = "America/Chicago"),
    local_date = as.Date(timestamp_local),
    hour_decimal = hour(timestamp_local) +
      minute(timestamp_local) / 60 +
      second(timestamp_local) / 3600,
    clock_rad = hour_decimal * 2 * pi / 24,
    # Use one non-site-specific location for every detection: 1910 Harris
    # County Courthouse (latitude 29.760890, longitude -95.359704).
    latitude = 29.760890,
    longitude = -95.359704
  )
# Coordinates must be longitude first, latitude second
coords <- cbind(
  all_detections$longitude,
  all_detections$latitude
)

# Use POSIXct local timestamp, NOT local_date
all_detections$suntime_rad <- sunTime(
  clockTime = all_detections$clock_rad,
  Dates = all_detections$timestamp_local,
  Coords = coords
)
summary(all_detections$suntime_rad)
range(all_detections$suntime_rad, na.rm = TRUE)
sum(is.na(all_detections$suntime_rad))
class(all_detections$timestamp_local)
table(all_detections$site_status, all_detections$species_group)
unique(all_detections$site_status)
site_summary <- all_detections %>%
  filter(!is.na(site), !is.na(species_group)) %>%
  distinct(site, species_group) %>%
  group_by(site) %>%
  summarise(
    has_bobcat = any(species_group == "bobcat"),
    has_domestic_cat = any(species_group == "domestic_cat"),
    .groups = "drop"
  ) %>%
  mutate(
    site_category = case_when(
      has_bobcat & has_domestic_cat ~ "both",
      has_bobcat & !has_domestic_cat ~ "bobcat_solo",
      !has_bobcat & has_domestic_cat ~ "domestic_cat_solo",
      TRUE ~ NA_character_
    )
  )

site_summary
table(site_summary$site_category)
all_detections <- all_detections %>%
  left_join(
    site_summary %>% select(site, site_category),
    by = "site"
  )
table(all_detections$site_category, all_detections$species_group)
all_detections_clean <- all_detections %>%
  filter(
    !is.na(suntime_rad),
    !is.na(site_category),
    !is.na(species_group)
  )
table(all_detections_clean$site_category, all_detections_clean$species_group)
solo_bobcat_times <- all_detections_clean %>%
  filter(
    site_category == "bobcat_solo",
    species_group == "bobcat"
  ) %>%
  pull(suntime_rad)

solo_domestic_cat_times <- all_detections_clean %>%
  filter(
    site_category == "domestic_cat_solo",
    species_group == "domestic_cat"
  ) %>%
  pull(suntime_rad)

solo_bobcat_times <- solo_bobcat_times[!is.na(solo_bobcat_times)]
solo_domestic_cat_times <- solo_domestic_cat_times[!is.na(solo_domestic_cat_times)]

length(solo_bobcat_times)
length(solo_domestic_cat_times)
both_bobcat_times <- all_detections_clean %>%
  filter(
    site_category == "both",
    species_group == "bobcat"
  ) %>%
  pull(suntime_rad)

both_domestic_cat_times <- all_detections_clean %>%
  filter(
    site_category == "both",
    species_group == "domestic_cat"
  ) %>%
  pull(suntime_rad)

both_bobcat_times <- both_bobcat_times[!is.na(both_bobcat_times)]
both_domestic_cat_times <- both_domestic_cat_times[!is.na(both_domestic_cat_times)]

length(both_bobcat_times)
length(both_domestic_cat_times)
solo_estimator <- ifelse(
  min(length(solo_bobcat_times), length(solo_domestic_cat_times)) < 50,
  "Dhat1",
  "Dhat4"
)

both_estimator <- ifelse(
  min(length(both_bobcat_times), length(both_domestic_cat_times)) < 50,
  "Dhat1",
  "Dhat4"
)

solo_estimator
both_estimator
solo_dhat <- overlapEst(
  solo_bobcat_times,
  solo_domestic_cat_times,
  type = solo_estimator
)

both_dhat <- overlapEst(
  both_bobcat_times,
  both_domestic_cat_times,
  type = both_estimator
)

solo_dhat
both_dhat
dhat_results <- data.frame(
  comparison = c("solo_sites", "both_sites"),
  n_bobcat = c(length(solo_bobcat_times), length(both_bobcat_times)),
  n_domestic_cat = c(length(solo_domestic_cat_times), length(both_domestic_cat_times)),
  estimator = c(solo_estimator, both_estimator),
  dhat = c(solo_dhat, both_dhat)
)

dhat_results
overlapPlot(
  solo_bobcat_times,
  solo_domestic_cat_times,
  main = paste0(
    "SunTime Activity Overlap: Solo Sites\n",
    "Dhat = ", round(solo_dhat, 3),
    " | Estimator = ", solo_estimator
  )
)

legend(
  "topright",
  legend = c("Bobcats at bobcat-only sites", "Domestic cats at cat-only sites"),
  lty = c(1, 2),
  bty = "n"
)
overlapPlot(
  both_bobcat_times,
  both_domestic_cat_times,
  main = paste0(
    "SunTime Activity Overlap: Shared Sites\n",
    "Dhat = ", round(both_dhat, 3),
    " | Estimator = ", both_estimator
  )
)

legend(
  "topright",
  legend = c("Bobcats at shared sites", "Domestic cats at shared sites"),
  lty = c(1, 2),
  bty = "n"
)
set.seed(123)

solo_boot <- bootstrap(
  solo_bobcat_times,
  solo_domestic_cat_times,
  nb = 1000,
  type = solo_estimator
)

solo_ci <- bootCI(solo_dhat, solo_boot)

solo_ci
set.seed(123)

both_boot <- bootstrap(
  both_bobcat_times,
  both_domestic_cat_times,
  nb = 1000,
  type = both_estimator
)

both_ci <- bootCI(both_dhat, both_boot)

both_ci
png("suntime_overlap_solo_sites.png", width = 1000, height = 800)

overlapPlot(
  solo_bobcat_times,
  solo_domestic_cat_times,
  main = paste0(
    "SunTime Activity Overlap: Solo Sites\n",
    "Dhat = ", round(solo_dhat, 3),
    " | Estimator = ", solo_estimator
  )
)

legend(
  "topright",
  legend = c("Bobcats at bobcat-only sites", "Domestic cats at cat-only sites"),
  lty = c(1, 2),
  bty = "n"
)

dev.off()
png("suntime_overlap_both_sites.png", width = 1000, height = 800)

overlapPlot(
  both_bobcat_times,
  both_domestic_cat_times,
  main = paste0(
    "SunTime Activity Overlap: Shared Sites\n",
    "Dhat = ", round(both_dhat, 3),
    " | Estimator = ", both_estimator
  )
)

legend(
  "topright",
  legend = c("Bobcats at shared sites", "Domestic cats at shared sites"),
  lty = c(1, 2),
  bty = "n"
)

dev.off()
