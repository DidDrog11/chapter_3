
# Read data exported from the Rodent trapping repository (https://github.com/DidDrog11/rodent_trapping) #
# which cleans the data entered by the field team on the location of traps and rodents trapped.         #

combined_data <- readRDS(gzcon(url("https://github.com/DidDrog11/rodent_trapping/raw/main/data/data_for_export/combined_data.rds")))

write_rds(combined_data, here("data", "input", "combined_data.rds"))

# Prepare data into sites and detections ----------------------------------

detections <- combined_data$rodent_data %>%
  filter(trap_uid %in% combined_data$trap_data$trap_uid) %>% # only keep rodents for sites with coordinates
  filter(!str_detect(rodent_uid, "BAM")) %>% # remove Bambawo
  filter(!str_detect(rodent_uid, "4_BAI_018")) %>% # the trap number for this rodent of 4_BAI_020 is wrong, remove this one for now
  drop_na(clean_names) %>% # remove individuals without a genus identification
  mutate(village = str_split(as.character(trap_uid), "_", simplify = TRUE)[, 1],
         visit = as.numeric(str_split(as.character(trap_uid), "_", simplify = TRUE)[, 2]),
         trap_night = as.numeric(str_split(as.character(trap_uid), "_", simplify = TRUE)[, 3]),
         grid_number = as.numeric(str_split(as.character(trap_uid), "_", simplify = TRUE)[, 4]),
         trap_number = as.numeric(str_split(as.character(trap_uid), "_", simplify = TRUE)[, 5])) %>%
  select(village, visit, grid_number, trap_number, trap_uid, clean_names) %>%
  mutate(clean_names = case_when(clean_names == "mus_spp" & grid_number %in% c(6, 7) ~ "mus_musculus",
                                 clean_names == "mus_spp" ~ "mus_minutoides",
                                 TRUE ~ clean_names)) %>% # for now we will assign all village trapped mus to mus_musculus and all others to mus_minutoides
  distinct(village, visit, grid_number, trap_number, clean_names, .keep_all = TRUE) # keep all distinct detections of each species at each site for each replicate

# Each four night trapping activity will be considered as a single replicate  #
# The exact location of a trap varied between replicates. To incorporate      #
# replicates we will re-number the trap sites based on the closes located     #
# trap in a subsequent visit. To do this we convert coordinates to projected  #
# UTM 29N for Sierra Leone this is EPSG:32629.                                #

sites <- combined_data$trap_data %>%
  select(date_set, village, trap_uid, visit, grid_number, trap_number, habitat, proximity, site_use, elevation, geometry) %>%
  filter(village != "bambawo") %>% # remove Bambawo as only used for one replicate
  mutate(visit = as.numeric(as.character(visit)),
         grid_number = as.numeric(as.character(grid_number)),
         grid_number = case_when(grid_number == 6 ~ 7,
                                 TRUE ~ grid_number), # combining 6 and 7 as overlap spatially
         longitude = st_coordinates(geometry)[, 1],
         latitude = st_coordinates(geometry)[, 2]) %>%
  tibble(.) %>%
  select(-geometry) %>%
  distinct(village, visit, grid_number, trap_number, longitude, latitude, .keep_all = TRUE) %>%
  select(village, visit, grid_number, trap_number, longitude, latitude) %>%
  st_as_sf(coords = c("longitude", "latitude")) %>%
  st_set_crs(value = "EPSG:4326") %>%
  st_transform(crs = "EPSG:32629") %>%
  mutate(easting = st_coordinates(geometry)[, 1],
         northing = st_coordinates(geometry)[, 2])

assign_traps_to_cells <- function(all_sites = sites) {
  
  select_site <- all_sites %>%
    group_by(village, grid_number) %>%
    group_split()
  
  # make a 25 m^2 grid for each site, with a 1 meter buffer from the traps  #
  
  grids <- lapply(select_site, function(x) {
    st_make_grid(st_buffer(x, 1), cellsize = 5, square = TRUE)
  })
  
  sites_in_grid <- mapply(function(X, Y) {
                                 containing_grid <- st_within(X, Y)
                               },
                          X = select_site,
                          Y = grids)
  
}

a <- which.min(st_distance(t[[2]][1, ], t[[1]]))
b <- which.min(st_distance(t[[2]][-a], t[[1]][2, ]))

a <- st_buffer(t[[1]], 2)

b <- st_within(t[[2]], a) %>%
  lengths > 0

t[[2]][st_within(t[[2]], a) %>%
    lengths > 0, ]
t[[1]][st_contains(a, t[[2]]) %>%
         lengths > 0, ]

t[[2]][st_within(t[[2]], a) %>%
    lengths > 0, ]

closest_trap <- c(as.numeric(NA))

for(i in seq(unique(t[[2]]$trap_number))) {
  
  closest_trap[i] <- which.min(st_distance(t[[2]][i, ], t[[1]][!t[[1]]$trap_number %in% closest_trap, ]))
  
}

t[[1]]$nearest_2 <- as.numeric(NA)

for(i in unique(t[[1]]$trap_number)) {
  
  nearest <- st_nearest_feature(t[[1]][i, ], t[[2]])
  
  t[[1]] <- t[[1]] %>%
    mutate(nearest_2 = case_when(trap_number == i ~ as.numeric(nearest),
                                 TRUE ~ as.numeric(nearest_2)))
  
  t[[2]] <- t[[2]][-nearest, ]
  
}

nearest <- tibble(visit_2_tn = c(t[[2]]$trap_number))

for(i in seq(unique(t[[1]]$trap_number))) {
  
  distance = c(st_distance(t[[1]][i, ], t[[2]]))
  
  nearest[, paste0(t[[1]][i, ]$trap_number)] <- distance
  
}

a <- nearest %>%
  pivot_longer(cols = 2:50, names_to = "visit_1", values_to = "distance") %>%
  mutate(rank_distance = rank(distance, ties.method = "first")) %>%
  group_by(visit_2_tn) %>%
  mutate(rank_closest = rank(distance, ties.method = "first")) %>%
  arrange(rank_closest, rank_distance)
  
p <- a %>%
  group_by(visit_2_tn) %>%
  slice_min(rank_closest)
first <- a %>%
  group_by(visit_1) %>%
  slice_min(rank_distance) %>%
  ungroup() %>%
  group_by(visit_2_tn) %>%
  slice_min(rank_closest)

second <- a %>%
  filter(!visit_2_tn %in% first$visit_2_tn,
         !visit_1 %in% first$visit_1) %>%
  group_by(visit_1) %>%
  slice_min(rank_distance) %>%
  ungroup() %>%
  group_by(visit_2_tn) %>%
  slice_min(rank_closest)

third <- a %>%
  filter(!visit_2_tn %in% c(first$visit_2_tn, second$visit_2_tn),
         !visit_1 %in% c(first$visit_1, second$visit_1)) %>%
  group_by(visit_1) %>%
  slice_min(rank_distance) %>%
  ungroup() %>%
  group_by(visit_2_tn) %>%
  slice_min(rank_closest)

fourth <- a %>%
  filter(!visit_2_tn %in% c(first$visit_2_tn, second$visit_2_tn, third$visit_2_tn),
         !visit_1 %in% c(first$visit_1, second$visit_1, third$visit_1)) %>%
  group_by(visit_1) %>%
  slice_min(rank_distance) %>%
  ungroup() %>%
  group_by(visit_2_tn) %>%
  slice_min(rank_closest)

fifth <- a %>%
  filter(!visit_2_tn %in% c(first$visit_2_tn, second$visit_2_tn, third$visit_2_tn, fourth$visit_2_tn),
         !visit_1 %in% c(first$visit_1, second$visit_1, third$visit_1, fourth$visit_1)) %>%
  group_by(visit_1) %>%
  slice_min(rank_distance) %>%
  ungroup() %>%
  group_by(visit_2_tn) %>%
  slice_min(rank_closest)

sixth <- a %>%
  filter(!visit_2_tn %in% c(first$visit_2_tn, second$visit_2_tn, third$visit_2_tn, fourth$visit_2_tn, fifth$visit_2_tn),
         !visit_1 %in% c(first$visit_1, second$visit_1, third$visit_1, fourth$visit_1, fifth$visit_2_tn)) %>%
  group_by(visit_1) %>%
  slice_min(rank_distance) %>%
  ungroup() %>%
  group_by(visit_2_tn) %>%
  slice_min(rank_closest)

a <- tibble(trap_number = c(t[[2]]$trap_number),
              distance = c(st_distance(t[[1]][1, ], t[[2]])))

st_nearest_feature(t[[1]], t[[2]])

trap_buffer <- st_buffer(sites, 5)

sf_int <- st_intersects(trap_buffer, sites)

ggplot() +
  geom_sf(data = sites %>% filter(village == "lambayama")) +
  geom_sf(data = sites[c(sf_int[[130]]), ], aes(colour = trap_number))
