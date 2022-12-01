source(here::here("R", "00_setup.R"))

observed_data <- read_rds(here("data", "observed_data", "descriptive_data.rds"))

sites <- observed_data$sites_grids$select_site %>%
  bind_rows() %>%
  select(date_set, site_id = unique_site, site, landuse, site_easting, site_northing, village, visit, grid_number, trap_id = trap_uid, trap_easting, trap_northing, elevation)

detections <- observed_data$detections %>%
  rename(trap_id = trap_uid) %>%
  left_join(sites %>%
              select(site_id, trap_id, date_set))

tn <- sites %>%
  group_by(village, visit) %>%
  summarise(tn = n())

season <- tibble(month = 1:12, season = c(rep("Dry", 4), rep("Rainy", 6), rep("Dry", 2)))
  
season_detection <- detections %>%
  mutate(month = month(date_set)) %>%
  left_join(tn) %>%
  left_join(season) %>%
  group_by(clean_names, tn, visit, season) %>%
  summarise(n = n(),
            tn = mean(tn)) %>%
  ungroup() %>%
  rowwise() %>%
  mutate(detection_rate = (n/tn)*100) %>%
  ggplot() +
  geom_col(aes(x = visit, y = detection_rate, fill = season)) +
  facet_wrap(~ clean_names) +
  theme_bw() +
  labs(fill = "Season",
       y = "Detection rate (/100 TN)",
       x = "Visit")
