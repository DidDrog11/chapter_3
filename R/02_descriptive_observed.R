source(here::here("R", "00_setup.R"))

observed_data <- read_rds(here("data", "observed_data", "descriptive_data.rds"))

# This dataframe includes a row for each trap placed within a grid. The four nights of trapping observations are grouped as a single replicate. 
sites <- observed_data$sites_grids$select_site %>%
  bind_rows() %>%
  select(date_set, site_id = unique_site, site, landuse, site_easting, site_northing, village, visit, grid_number, trap_id = trap_uid, trap_easting, trap_northing, elevation)

# This dataframe includes a row for each rodent individual trapped during the study. 
detections <- observed_data$detections %>%
  rename(trap_id = trap_uid) %>%
  left_join(sites %>%
              select(site_id, trap_id, date_set))

# This calculates the number of trap nights performed at each site
tn <- sites %>%
  group_by(village, visit) %>%
  summarise(tn = n() * 4)

tn_village_landuse <- sites %>%
  group_by(village, landuse, visit) %>%
  mutate(tn = n() * 4) %>%
  distinct(trap_id, village, visit, landuse, tn)

# Allocate seasons to months
season <- tibble(month = 1:12, season = c(rep("Dry", 4), rep("Rainy", 6), rep("Dry", 2)))

grids <- list()

for(i in 1:length(observed_data$sites_grids$grids_for_plotting)) {
  
  grid_site <- names(observed_data$sites_grids$grids_for_plotting[i])
  
  grids[[i]] <- tibble(observed_data$sites_grids$grids_for_plotting[[i]]) %>%
    mutate(grid_id = grid_site)
  
}

grids <- grids %>%
  bind_rows() %>%
  st_as_sf() %>%
  rename(geometry = 1)

# Description trapping effort -------------------------------------------

number_trap_nights <- sum(tn$tn)

number_visits <- length(unique(sites$visit))

median_number_trap_visit <- tn %>%
  group_by(village) %>%
  summarise(tn_median = median(tn),
            tn_min = min(tn),
            tn_max = max(tn))

# Raw area of traps -------------------------------------------------------
trap_locations <- sites %>%
  select(village, grid_number, visit, trap_id, trap_easting, trap_northing) %>%
  st_as_sf(coords = c("trap_easting", "trap_northing"), crs = SL_UTM) %>%
  group_by(village, grid_number, visit) %>%
  summarise()

trap_areas <- trap_locations %>%
  st_convex_hull()

trap_areas$area <- st_area(trap_areas$geometry)

tibble(trap_areas %>%
  filter(grid_number != 7)) %>% 
  summarise(mean(area))

# Description trapping locations ------------------------------------------

land_use_type <- sites %>%
  tibble() %>%
  group_by(village, landuse) %>%
  summarise(tn = n() * 4)

# plot of grids by village
traps_in_grid <- mapply(X = sites %>%
                          st_as_sf(coords = c("trap_easting", "trap_northing"), crs = SL_UTM) %>%
                          mutate(trap_number = str_split(trap_id, pattern = "_", simplify = TRUE)[, 4]) %>%
                          group_by(village, grid_number) %>%
                          group_split(),
                        Y = grids %>%
                          mutate(village = str_split(grid_id, "_", simplify = TRUE)[, 1],
                                 grid = as.numeric(str_split(grid_id, "_", simplify = TRUE)[, 2])) %>%
                          group_by(village, grid) %>%
                          group_split(),
                        function(X, Y) {
                          list(st_join(X, Y, st_within))
                        }) %>%
  do.call(rbind, .) %>%
  select(village = village.x, visit, grid_number, trap_number, landuse, grid_id, trap_point = geometry)

grids_with_traps <- mapply(X = sites %>%
                             st_as_sf(coords = c("trap_easting", "trap_northing"), crs = SL_UTM) %>%
                             mutate(trap_number = str_split(trap_id, pattern = "_", simplify = TRUE)[, 4]) %>%
                             group_by(village, grid_number) %>%
                             group_split(),
                           Y = grids %>%
                             mutate(village = str_split(grid_id, "_", simplify = TRUE)[, 1],
                                    grid = as.numeric(str_split(grid_id, "_", simplify = TRUE)[, 2])) %>%
                             group_by(village, grid) %>%
                             group_split(),
                           function(X, Y) {
                             list(st_join(Y, X, st_contains, left = FALSE))
                           }) %>%
  do.call(rbind, .) %>%
  select(village = village.x, grid_number, landuse, grid_id, grid_polygon = geometry)

f <- function(i) paste(i, collapse="_")

grids_with_traps$grp <- factor(sapply(st_equals(grids_with_traps), f))

grids_with_traps <- grids_with_traps %>%
  group_by(grp) %>%
  mutate(tn = n() * 4) %>%
  distinct() %>%
  group_by(village) %>%
  group_split()

osm_bbox <- lapply(X = grids_with_traps, function(X) {X %>% group_by(village) %>% mutate(geometry = st_union(.)) %>% distinct(village, geometry)}) %>%
  do.call(rbind, .) %>%
  st_as_sf(crs = SL_UTM) %>%
  st_transform(crs = default_CRS) %>%
  group_split()

names(osm_bbox) <- c("baiama", "lalehun", "lambayama", "seilama")

if(!file.exists(here("data", "observed_data", "baiama_raster.tif"))) {
  
  bg = lapply(X = osm_bbox, function(X) {osm.raster(extract_bbox(st_bbox(X)), zoomin = + 2)})
  bg = lapply(X = bg, function(X) {rast(X) %>%
      project(y = SL_UTM)})
  
  writeRaster(bg$baiama, here("data", "observed_data", "baiama_raster.tif"), overwrite = TRUE)
  writeRaster(bg$lalehun, here("data", "observed_data", "lalehun_raster.tif"), overwrite = TRUE)
  writeRaster(bg$lambayama, here("data", "observed_data", "lambayama_raster.tif"), overwrite = TRUE)
  writeRaster(bg$seilama, here("data", "observed_data", "seilama_raster.tif"), overwrite = TRUE)
}

fig_1_df <- grids_with_traps

write_rds(fig_1_df, here("data", "observed_data", "fig_1_df.rds"))

# Description rodents trapped ---------------------------------------------

number_rodents <- nrow(detections)

trap_success_rate <- round(number_rodents/number_trap_nights * 100, 1)

trap_success_rate_df <- detections %>%
  select(rodent_id, trap_id) %>%
  left_join(sites,
            by = "trap_id") %>%
  distinct() %>%
  drop_na() %>%
  group_by(village, landuse) %>%
  summarise(n_rodents = n()) %>%
  left_join(land_use_type, by = c("village", "landuse")) %>%
  mutate(trap_success = round(n_rodents/tn * 100, 1))


# Description trap rate by season ----------------------------------------------

landuse_order <- c("village", "agriculture", "forest")
names(landuse_order) <- c("Village", "Agriculture", "Forest")

season_detection <- detections %>%
  mutate(month = month(date_set)) %>%
  left_join(tn) %>%
  left_join(season) %>%
  group_by(clean_names, village, tn, season) %>%
  summarise(n = n(),
            tn = unique(tn)) %>%
  ungroup() %>%
  rowwise() %>%
  mutate(trap_rate = (n/tn)*100,
         clean_names = str_to_sentence(str_replace(clean_names, "_", " "))) %>%
  ggplot() +
  geom_boxplot(aes(x = season, y = trap_rate, fill = season)) +
  facet_wrap(~ clean_names) +
  theme_bw() +
  labs(fill = element_blank(),
       y = "Trap rate (/100 TN)",
       x = "Season")

season_detection_landuse <- detections %>%
  mutate(month = month(date_set)) %>%
  left_join(tn_village_landuse, by = c("trap_id", "village", "visit")) %>%
  left_join(season) %>%
  group_by(clean_names, village, landuse, tn, season) %>%
  summarise(n = n(),
            tn = unique(tn)) %>%
  ungroup() %>%
  rowwise() %>%
  mutate(trap_rate = (n/tn)*100,
         clean_names = str_to_sentence(str_replace(clean_names, "_", " ")),
         landuse = factor(landuse, levels = landuse_order, labels = names(landuse_order))) %>%
  ggplot() +
  geom_boxplot(aes(x = landuse, y = trap_rate, fill = season), position = position_dodge(preserve = "single")) +
  facet_wrap(~ clean_names, scales = "free_y") +
  theme_bw() +
  labs(fill = element_blank(),
       y = "Trap rate (/100 TN)",
       x = "Landuse")

mastomys_season_detection <- detections %>%
  filter(clean_names == "mastomys_spp") %>%
  mutate(month = month(date_set)) %>%
  left_join(tn_village_landuse, by = c("trap_id", "village", "visit")) %>%
  left_join(season) %>%
  group_by(clean_names, village, landuse, tn, season) %>%
  summarise(n = n(),
            tn = unique(tn)) %>%
  ungroup() %>%
  rowwise() %>%
  mutate(trap_rate = (n/tn)*100,
         clean_names = str_to_sentence(str_replace(clean_names, "_", " ")),
         landuse = factor(landuse, levels = landuse_order, labels = names(landuse_order))) %>%
  ggplot() +
  geom_boxplot(aes(x = landuse, y = trap_rate, fill = season), position = position_dodge(preserve = "single")) +
  facet_wrap(~ village, scales = "free_y") +
  theme_bw() +
  labs(fill = element_blank(),
       y = "Trap rate (/100 TN)",
       x = "Landuse")

# Description species trapped ---------------------------------------------

species_trapped <- detections %>%
  group_by(clean_names, village) %>%
  summarise(n = n()) %>%
  drop_na(clean_names) %>%
  group_by(clean_names) %>%
  mutate(N = sum(n)) %>%
  arrange(-N)

species_order <- unique(str_to_sentence(str_replace_all(species_trapped$clean_names, "_", " ")))
village_order <- c("baiama", "lalehun", "lambayama", "seilama")
names(village_order) <- c("Baiama", "Lalehun", "Lambayama", "Seilama")

species_trap_success_rate <- detections %>%
  select(site_id, clean_names) %>%
  group_by(site_id, clean_names) %>%
  summarise(n = n()) %>%
  left_join(sites %>%
              tibble() %>%
              select(site_id, village, landuse),
            by = "site_id") %>%
  distinct() %>%
  group_by(clean_names, village, landuse) %>%
  summarise(n = sum(n)) %>%
  drop_na(clean_names) %>%
  left_join(., trap_success_rate_df %>%
              group_by(village, landuse) %>%
              summarise(tn = sum(tn))) %>%
  rowwise() %>%
  mutate(trap_success = round(n/tn * 100, 1),
         value = paste0(n, " (", if(trap_success <= 0.1) "<0.1" else trap_success, "%)"),
         clean_names = factor(str_to_sentence(str_replace_all(clean_names, "_", " ")), levels = species_order),
         landuse = factor(landuse, levels = landuse_order, labels = names(landuse_order)),
         village = factor(village, levels = village_order, labels = names(village_order))) %>%
  select(clean_names, village, landuse, value) %>%
  arrange(landuse) %>%
  pivot_wider(names_from = landuse, values_from = value, values_fill = "-") %>%
  arrange(clean_names, village)

species_trap_success_village <- detections %>%
  select(rodent_id, site_id, clean_names) %>%
  left_join(sites %>%
              tibble() %>%
              select(site_id, village, landuse),
            by = "site_id") %>%
  distinct()  %>%
  group_by(clean_names, village, landuse) %>%
  summarise(n = n()) %>%
  drop_na(clean_names) %>%
  left_join(., trap_success_rate_df %>%
              group_by(village) %>%
              summarise(tn = sum(tn))) %>%
  group_by(clean_names, village) %>%
  summarise(n = sum(n),
            tn = median(tn)) %>%
  ungroup() %>%
  rowwise() %>%
  mutate(trap_success = round(n/tn * 100, 1),
         value = paste0(n, " (", if(trap_success <= 0.1) "<0.1" else trap_success, "%)"),
         clean_names = factor(str_to_sentence(str_replace_all(clean_names, "_", " ")), levels = species_order),
         village = factor(village, levels = village_order, labels = names(village_order)),
         Combined = value) %>%
  select(clean_names, village, Combined)

table_1b <- left_join(species_trap_success_village, species_trap_success_rate, by = c("clean_names", "village"))

write_rds(table_1b, here("output", "table_1b.rds"))

# Table 1b as plot --------------------------------------------------------

fig_2_df <- detections %>%
  left_join(sites,
            by = c("village", "visit", "grid_number", "trap_id", "site_id", "date_set")) %>%
  group_by(landuse, village, clean_names) %>%
  summarise(n_detected = n()) %>%
  left_join(trap_success_rate_df %>%
              select(-trap_success), by = c("village", "landuse")) %>%
  rowwise() %>%
  mutate(detection_rate = n_detected/tn * 1000,
         clean_names = str_to_sentence(str_replace(clean_names, "_", " ")),
         village = factor(str_to_sentence(village), levels = c("Baiama", "Lalehun", "Lambayama", "Seilama")),
         landuse = factor(str_to_sentence(landuse), levels = c("Forest", "Agriculture", "Village")))

fig_2 <- fig_2_df %>%
  ggplot(aes(x = village, y = detection_rate, fill = landuse, label = n_detected)) +
  geom_col(position = position_dodge2(preserve = "single", width = 0.7)) +
  geom_text(size = 3, position = position_dodge2(preserve = "single", width = 0.8)) +
  scale_y_continuous(breaks = scales::pretty_breaks()) +
  scale_fill_manual(values = landuse_palette) +
  facet_wrap(~ clean_names, scales = "free_y") +
  labs(x = element_blank(),
       y = "Detection rate per 1,000 trap nights",
       fill = "Landuse") +
  theme_bw()

save_plot(plot = fig_2, filename = here("output", "Figure_2.png"), base_width = 13, base_height = 7)

# Species diversity -------------------------------------------------------

richness <- sites %>%
  tibble() %>%
  distinct(village, grid_number, site_id, landuse) %>%
  left_join(detections, by = c("village", "grid_number", "site_id")) %>%
  group_by(village, grid_number, site_id, landuse, clean_names) %>%
  summarise(n = n()) %>%
  ungroup()

richness$n[is.na(richness$clean_names)] <- 0
richness$n[!is.na(richness$clean_names)] <- 1

richness_village <- richness %>%
  select(village, clean_names, n) %>%
  filter(n != 0) %>%
  distinct() %>%
  group_by(village) %>%
  summarise(species_richness = n()) %>%
  mutate(village = factor(village, levels = village_order, labels = names(village_order)))

richness_landuse <- richness %>%
  select(landuse, clean_names, n) %>%
  filter(n != 0) %>%
  distinct() %>%
  group_by(landuse) %>%
  summarise(species_richness = n()) %>%
  mutate(landuse = factor(landuse, levels = landuse_order, labels = names(landuse_order)),
         village = "Combined")

richness_landuse_village <- richness %>%
  select(village, landuse, clean_names, n) %>%
  filter(n != 0) %>%
  distinct() %>%
  group_by(village, landuse) %>%
  summarise(species_richness = n()) %>%
  mutate(landuse = factor(landuse, levels = landuse_order, labels = names(landuse_order)),
         village = factor(village, levels = village_order, labels = names(village_order)))

richness_combined <- bind_rows(richness_village, richness_landuse, richness_landuse_village)

diversity <- sites %>%
  tibble() %>%
  distinct(village, grid_number, site_id, trap_id, landuse) %>%
  left_join(detections, by = c("village", "grid_number", "trap_id", "site_id")) %>%
  group_by(village, grid_number, site_id, landuse, clean_names) %>%
  summarise(n = sum(n())) %>%
  ungroup()

diversity$n[is.na(diversity$clean_names)] <- 0

diversity_village <- diversity %>%
  filter(!is.na(clean_names)) %>%
  group_by(village, clean_names) %>%
  summarise(n = sum(n)) %>%
  group_by(village) %>%
  summarise(N = sum(n),
            shannon_diversity = diversity(n, index = "shannon", MARGIN = 2)) %>%
  mutate(village = factor(village, levels = village_order, labels = names(village_order))) %>%
  arrange(village)

diversity_landuse <- diversity %>%
  filter(!is.na(clean_names)) %>%
  group_by(landuse, clean_names) %>%
  summarise(n = sum(n)) %>%
  group_by(landuse) %>%
  summarise(N = sum(n),
            shannon_diversity = diversity(n, index = "shannon", MARGIN = 2)) %>%
  mutate(landuse = factor(landuse, levels = landuse_order, labels = names(landuse_order)),
         village = "Combined") %>%
  arrange(landuse)

diversity_landuse_village <- diversity %>%
  filter(!is.na(clean_names)) %>%
  group_by(village, landuse, clean_names) %>%
  summarise(n = sum(n)) %>%
  group_by(village, landuse) %>%
  summarise(N = sum(n),
            shannon_diversity = diversity(n, index = "shannon", MARGIN = 2)) %>%
  arrange(-shannon_diversity) %>%
  mutate(landuse = factor(landuse, levels = landuse_order, labels = names(landuse_order)),
         village = factor(village, levels = village_order, labels = names(village_order))) %>%
  arrange(village, landuse)

diversity_combined <- bind_rows(diversity_village, diversity_landuse, diversity_landuse_village)

trap_denominator <- sites %>%
  tibble() %>%
  group_by(village, visit, site_id, landuse) %>%
  summarise(traps = n() * 4) %>%
  group_by(village, landuse) %>%
  summarise(TN = sum(traps)) %>%
  group_by(landuse) %>%
  mutate(TN_landuse = sum(TN)) %>%
  group_by(village) %>%
  mutate(TN_village = sum(TN)) %>%
  ungroup()

trap_village <- trap_denominator %>%
  select(village, TN = TN_village) %>%
  distinct() %>%
  mutate(village = factor(village, levels = village_order, labels = names(village_order)))

trap_landuse <- trap_denominator %>%
  select(landuse, TN  = TN_landuse) %>%
  distinct() %>%
  mutate(landuse = factor(landuse, levels = landuse_order, labels = names(landuse_order)),
         village = "Combined")

trap_landuse_village <- trap_denominator %>%
  select(village, landuse, TN) %>%
  distinct() %>%
  mutate(landuse = factor(landuse, levels = landuse_order, labels = names(landuse_order)),
         village = factor(village, levels = village_order, labels = names(village_order)))

trap_combined <- bind_rows(trap_village, trap_landuse, trap_landuse_village)

table_1a <- richness_combined %>%
  left_join(diversity_combined, by = c("village", "landuse")) %>%
  left_join(trap_combined, by = c("village", "landuse")) %>%
  select(village, landuse, N, TN, species_richness, shannon_diversity) %>%
  mutate(village = factor(village, levels = c("Combined", str_to_sentence(village_order)), labels = c("All villages", names(village_order))),
         landuse = case_when(is.na(landuse) ~ "Combined",
                             TRUE ~ as.character(landuse)),
         landuse = factor(landuse, levels = c("Village",
                                              "Agriculture",
                                              "Forest",
                                              "Combined")),
         shannon_diversity = round(shannon_diversity, 2),
         TN = paste0(TN, " (", round(N/TN * 100, 1), "%)")) %>%
  arrange(village, landuse)

write_rds(table_1a, here("output", "table_1a.rds"))

# NEEDS FURTHER READING IF GOING TO USE THIS #

dissimilarity_df <- richness %>%
  filter(!is.na(clean_names)) %>%
  group_by(village, landuse, clean_names) %>%
  summarise(n = sum(n)) %>%
  pivot_wider(names_from = clean_names, values_from = n, values_fill = 0) %>%
  ungroup() %>%
  select(-village, -landuse)

jaccard_dis <- vegdist(dissimilarity_df, method = "jaccard", binary = FALSE)
jaccard_hclus <- hclust(jaccard_dis)

jaccard_labs <- richness %>%
  filter(!is.na(clean_names)) %>%
  group_by(village, landuse, clean_names) %>%
  summarise(n = sum(n)) %>%
  pivot_wider(names_from = clean_names, values_from = n, values_fill = 0) %>%
  ungroup() %>%
  mutate(labs = paste0(str_sub(village, end = 3), "_", landuse)) %>%
  pull(labs)

jaccard_hclus[["labels"]] = jaccard_labs[jaccard_hclus[["order"]]]

plot(jaccard_hclus)

# Species accumulation graphs for supplementary
baiama_accum <- richness %>%
  filter(village == "baiama") %>%
  group_by(site_id, clean_names) %>%
  summarise(n = sum(n)) %>%
  ungroup() %>%
  pivot_wider(names_from = clean_names, values_from = n, values_fill = 0) %>%
  select(-site_id, -`NA`) %>%
  specaccum(comm = ., method = "exact")

lalehun_accum <- richness %>%
  filter(village == "lalehun") %>%
  group_by(site_id, clean_names) %>%
  summarise(n = sum(n)) %>%
  ungroup() %>%
  pivot_wider(names_from = clean_names, values_from = n, values_fill = 0) %>%
  select(-site_id, -`NA`) %>%
  specaccum(comm = ., method = "exact")

lambayama_accum <- richness %>%
  filter(village == "lambayama") %>%
  group_by(site_id, clean_names) %>%
  summarise(n = sum(n)) %>%
  ungroup() %>%
  pivot_wider(names_from = clean_names, values_from = n, values_fill = 0) %>%
  select(-site_id, -`NA`) %>%
  specaccum(comm = ., method = "exact")

seilama_accum <- richness %>%
  filter(village == "seilama") %>%
  group_by(site_id, clean_names) %>%
  summarise(n = sum(n)) %>%
  ungroup() %>%
  pivot_wider(names_from = clean_names, values_from = n, values_fill = 0) %>%
  select(-site_id, -`NA`) %>%
  specaccum(comm = ., method = "exact")

combined_accumulation <- tibble(Village = c(rep("Baiama", length(baiama_accum$sites)),
                                            rep("Lalehun", length(lalehun_accum$sites)),
                                            rep("Lambayama", length(lambayama_accum$sites)),
                                            rep("Seilama", length(seilama_accum$sites))),
                                Sites = c(baiama_accum$sites, lalehun_accum$sites, lambayama_accum$sites, seilama_accum$sites),
                                Richness = c(baiama_accum$richness, lalehun_accum$richness, lambayama_accum$richness, seilama_accum$richness),
                                sd = c(baiama_accum$sd, lalehun_accum$sd, lambayama_accum$sd, seilama_accum$sd))

accumulation_plot <- ggplot(combined_accumulation) +
  geom_line(aes(x = Sites, y = Richness, colour = Village)) +
  geom_ribbon(aes(x = Sites, ymin = Richness - sd, ymax = Richness + sd, colour = Village, fill = Village), alpha = 0.2) +
  scale_colour_manual(values = village_palette) +
  scale_fill_manual(values = village_palette) +
  theme_bw()

save_plot(plot = accumulation_plot, base_width = 10, filename = here("output", "species_accumulation.pdf"))


# Discussion description --------------------------------------------------
# Trap success in buildings, village and other for comparison
indoor_traps <- observed_data$sites_grids$select_site %>%
  bind_rows() %>%
  filter(site_habitat == "village_inside") %>%
  select(site_id = unique_site, village, visit, trap_uid) %>%
  mutate(tn = 4,
         landuse = "village_inside") %>%
  group_by(village, visit) %>%
  mutate(tn = sum(tn))


detections_visit_indoors <- detections %>%
  group_by(site_id, visit, village, clean_names) %>%
  summarise(n = n()) %>%
  filter(site_id %in% indoor_traps$site_id) %>%
  left_join(indoor_traps %>%
              tibble() %>%
              select(landuse, visit, village, tn) %>%
              distinct(),
            by = c("visit", "village")) %>%
  drop_na(landuse) %>%
  group_by(clean_names, village, visit, landuse) %>%
  summarise(n = sum(n),
            tn = unique(tn)) %>%
  group_by(village, visit) %>%
  summarise(n_all = sum(n),
            tn = unique(tn))

sum(detections_visit_indoors$n_all)/sum(detections_visit_indoors$tn)
