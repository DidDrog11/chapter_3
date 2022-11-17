source(here::here("R", "00_setup.R"))

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
         grid_number = case_when(grid_number == 6 ~ 7,
                                 TRUE ~ grid_number),
         trap_number = as.numeric(str_split(as.character(trap_uid), "_", simplify = TRUE)[, 5])) %>%
  select(village, visit, grid_number, trap_number, trap_uid, clean_names) %>%
  mutate(clean_names = case_when(clean_names == "mus_spp" & grid_number %in% c(6, 7) ~ "mus_musculus",
                                 clean_names == "mus_spp" ~ "mus_minutoides",
                                 TRUE ~ clean_names)) %>% # for now we will assign all village trapped mus to mus_musculus and all others to mus_minutoides
  mutate(visit = case_when(str_detect(village, "baiama|lambayama") & visit %in% c(1:4) ~ visit + 2,
                           TRUE ~ visit), # change the visit for baiama and lambayama to keep consistent numberings for dates
         trap_uid = factor(paste0(village, "_", visit, "_", grid_number, "_", trap_number))) %>% # change the trap_uid to reflect this
  distinct(village, visit, grid_number, trap_number, clean_names, .keep_all = TRUE) # keep all distinct detections of each species at each site for each replicate

# Each four night trapping activity will be considered as a single replicate  #
# The exact location of a trap varied between replicates. To incorporate      #
# replicates we will re-number the trap sites based on the closes located     #
# trap in a subsequent visit. To do this we convert coordinates to projected  #
# UTM 29N for Sierra Leone this is EPSG:32629.                                #

sites <- combined_data$trap_data %>%
  select(date_set, village, trap_uid, visit, grid_number, trap_number, habitat_group, site_habitat, site_use, elevation, geometry) %>%
  filter(village != "bambawo") %>% # remove Bambawo as only used for one replicate
  mutate(visit = as.numeric(as.character(visit)),
         grid_number = as.numeric(as.character(grid_number)),
         grid_number = case_when(grid_number == 6 ~ 7,
                                 TRUE ~ grid_number), # combining 6 and 7 as overlap spatially
         longitude = st_coordinates(geometry)[, 1],
         latitude = st_coordinates(geometry)[, 2],
         landuse = case_when(site_habitat == "forest" ~ "forest",
                             str_detect(site_habitat, "village") ~ "village",
                             str_detect(site_habitat, "banana|cassava|fallow|agriculture|palm|rice") ~ "agriculture",
                             TRUE ~ habitat_group),
         visit = case_when(str_detect(village, "baiama|lambayama") & visit %in% c(1:4) ~ visit + 2,
                           TRUE ~ visit), # change the visit for baiama and lambayama to keep consistent numberings for dates
         trap_uid = factor(paste0(village, "_", visit, "_", grid_number, "_", trap_number))) %>% # change the trap_uid to reflect this
  tibble(.) %>%
  select(-geometry) %>%
  distinct(village, visit, grid_number, trap_number, landuse, longitude, latitude, .keep_all = TRUE) %>%
  select(village, visit, grid_number, trap_number, landuse, longitude, latitude) %>%
  st_as_sf(coords = c("longitude", "latitude")) %>%
  st_set_crs(value = default_CRS) %>%
  st_transform(crs = SL_UTM) %>%
  mutate(trap_easting = st_coordinates(geometry)[, 1],
         trap_northing = st_coordinates(geometry)[, 2])

assign_traps_to_cells <- function(all_sites = sites) {
  
  # create list where each element is a grid from a village containing all  #
  # trap locations                                                          #
  
  select_site <- all_sites %>%
    group_by(village, grid_number) %>%
    group_split()
  
  # make a 49 m^2 grid for each site, with a 1 meter buffer from the traps  #
  
  grids <- lapply(select_site, function(x) {
    st_make_grid(st_buffer(x, 1), cellsize = 7, square = TRUE)
  })
  
  # identify the cells of the grids that contain traps                      #
  
  sites_in_grid <- mapply(function(X, Y) {
    containing_grid <- st_within(X, Y)
  },
  X = select_site,
  Y = grids)
  
  # allocate these cells to the trap locations                              #
  
  select_site <- mapply(function(X, Y) {
    list(X %>%
           mutate(site = c(unlist(Y))) %>%
           tibble() %>%
           select(-geometry))
  },
  X = select_site,
  Y = sites_in_grid)
  
  # get the centroid of each of these grid cells and append to the trap     #
  # locations. Name each grid cell uniquely based on village_grid_site      #
  
  grid_coords <- lapply(grids, function(x) {
    site_easting <- st_coordinates(st_centroid(x))[, 1]
    site_northing <- st_coordinates(st_centroid(x))[, 2]
    site <- c(1:length(x))
    
    return(tibble(site_easting = site_easting,
                  site_northing = site_northing,
                  site = site))
  })
  
  select_site <- mapply(function(X, Y) {
    
    list(X %>%
           left_join(Y, by = "site") %>%
           mutate(unique_site = paste0(village, "_", grid_number, "_", site),
                  trap_uid = paste0(village, "_", visit, "_", grid_number, "_", trap_number)))
    
  },
  X = select_site,
  Y = grid_coords)
  
  # name the element in the list based on the village and grid_number     #
  
  names(select_site) <- lapply(select_site, function(x) {
    
    name <- paste0(unique(x$village), "_", unique(x$grid_number))
    
    return(name)
  })
  
  names(grid_coords) <- lapply(select_site, function(x) {
    
    name <- paste0(unique(x$village), "_", unique(x$grid_number))
    
    return(name)
  })
  
  names(grids) <- names(grid_coords)
  
  return(list(select_site = select_site,
              grid_coords = grid_coords,
              grids_for_plotting = grids))
  
}

sites_grids <- assign_traps_to_cells(sites)

write_rds(sites_grids, here("data", "synthetic_data", "sites_grids.rds"))

sites_in_grid <- sites_grids$select_site
grid_coords <- sites_grids$grid_coords

for(i in 1:length(grid_coords)) {
  
  village_grid <- names(grid_coords[i])
  
  grid_coords[[i]] <- grid_coords[[i]] %>%
    mutate(unique_site = paste0(village_grid, "_", site))
  
}

grid_coords <- bind_rows(grid_coords)

# we can visualise the location of the newly produced sites below       #

visualise_sites_in_grid <- lapply(sites_in_grid, function(x) {
  
  x %>% 
    group_by(visit, site, site_easting, site_northing) %>%
    summarise(TN = n() * 4) %>%
    st_as_sf(coords = c("site_easting", "site_northing")) %>%
    st_set_crs(value = SL_UTM) %>%
    ggplot() +
    geom_sf(aes(colour = TN)) +
    facet_wrap(~ visit) +
    theme_bw()
  
})

# Produce synthetic data --------------------------------------------------
# As data is still be collected we will produce a synthetic dataset     #
# data will effectively be duplicated, with the max visit number = 12   #

duplicate_detections_lal_sei <- detections %>%
  filter(str_detect(village,  "lalehun|seilama")) %>%
  mutate(visit = visit + max(visit),
         trap_uid = factor(paste0(village, "_", visit, "_", grid_number, "_", trap_number))) %>%
  filter(visit <= 12)

duplicate_detections_lam_bai <- detections %>%
  filter(str_detect(village,  "lambayama|baiama")) %>%
  mutate(visit = (visit - 2) + max(visit),
         trap_uid = factor(paste0(village, "_", visit, "_", grid_number, "_", trap_number))) %>%
  filter(visit <= 12)

duplicate_detections <- bind_rows(duplicate_detections_lal_sei, duplicate_detections_lam_bai)
  
duplicate_sites <- lapply(sites_in_grid, function(x) {
  
  if(str_detect(unique(x$village), "lalehun|seilama")) {
    x %>%
      mutate(visit = visit + max(detections$visit),
             trap_uid = paste0(village, "_", visit, "_", grid_number, "_", trap_number),
             unique_site = paste0(village, "_", grid_number, "_", site))
  } else {
    
    x %>%
      mutate(visit = (visit - 2) + max(detections$visit),
             trap_uid = paste0(village, "_", visit, "_", grid_number, "_", trap_number),
             unique_site = paste0(village, "_", grid_number, "_", site))
    
  }
  
})

# Give visits by village consecutive numbers

synthetic_detections <- bind_rows(detections, duplicate_detections)

synthetic_sites <- bind_rows(bind_rows(sites_in_grid), bind_rows(duplicate_sites) %>%
                               filter(visit <= 12)) %>%
  arrange(unique_site)

synthetic_data <- list(synthetic_detections = synthetic_detections,
                       synthetic_sites = synthetic_sites)

write_rds(synthetic_data, here("data", "synthetic_data", "synthetic_data.rds"))

# Chapter 4 will use the same data format
chapter_4_data <- list(detections = detections,
                       sites = bind_rows(sites_in_grid))

write_rds(chapter_4_data, here("data", "data_for_export", "chapter_4_extract.rds"))

# Detection covariates ----------------------------------------------------

# Trap nights -------------------------------------------------------------
# Multiple traps are allocated to a single grid cell and 4 trap nights were conducted per trap
# We will use this as a measure of effort for detection

trap_nights <- lapply(sites_in_grid, function(X) {
  
  X <- X %>%
    group_by(unique_site, visit) %>%
    summarise(n_traps = n()) %>%
    mutate(trap_nights = n_traps * 4)
  
  return(X)
}) %>%
  bind_rows() %>%
  ungroup()

duplicate_trap_nights <- trap_nights %>%
  mutate(visit = visit + 6) %>%
  filter(visit <= 12)

combined_trap_nights <- bind_rows(trap_nights, duplicate_trap_nights) %>%
  select(unique_site, visit, trap_nights)

# Add date_set to the sites dataset to calculate moon and rainfall    #

current_dates <- trap_nights %>%
  left_join(bind_rows(sites_in_grid) %>%
              select(trap_uid, unique_site)) %>%
  distinct(trap_uid) %>%
  left_join(tibble(combined_data$trap_data) %>%
              select(date_set, village, visit, grid_number, trap_number) %>%
              mutate(visit = case_when(str_detect(village, "lambayama|baiama") & as.numeric(as.character(visit)) < 7 ~ as.numeric(as.character(visit)) + 2,
                                       TRUE ~ as.numeric(as.character(visit))),
                     grid_number = case_when(grid_number == as.character(6) ~ as.numeric(7),
                                             TRUE ~ as.numeric(as.character(grid_number))),
                     trap_uid = factor(paste0(village, "_", visit, "_", grid_number, "_", trap_number)),
                     date_set = case_when(village == "lalehun" & date_set == "2020-11-30" ~ ymd("2020-12-01"), # Change date to prevent big difference between nights rainfall
                                          TRUE ~ ymd(date_set))))

# Create future dates for the simulated visits, assuming every 3 months #
last_visit <- current_dates %>%
  group_by(village) %>%
  summarise(last_visit = max(visit),
            date_set = max(date_set))

last_visit$date_set[last_visit$village == "seilama"] <- ymd("2022-04-12") + days(5)

future_visit <- list()

for(i in 1:6) {
  
  future_visit[[i]] <- last_visit %>%
    mutate(visit = last_visit + i,
           date_set = ymd(date_set) + months(3 * i))
  
}

future_visit <- bind_rows(future_visit) %>%
  select(village, visit, date_set)
  
# Date set for all sites                                                #
date_set <- combined_trap_nights %>%
  mutate(village = str_split(unique_site, "_", simplify = TRUE)[ , 1]) %>%
  left_join(synthetic_data$synthetic_sites %>%
              select(unique_site, site_easting, site_northing)) %>%
  left_join(current_dates %>%
              distinct(village, visit, date_set),
            by = c("village", "visit")) %>%
  left_join(future_visit %>%
              distinct(village, visit, date_set), 
            by = c("village", "visit")) %>%
  mutate(date_set = coalesce(date_set.x, date_set.y)) %>%
  distinct(date_set, visit, unique_site, site_easting, site_northing)

# Monthly rainfall --------------------------------------------------------
# For the precipitation data we need to provide lon and lat points.     #
# We extract the centre of the trapping sites for this.                 #

tile_coords <- st_coordinates(st_centroid(st_as_sfc(st_bbox(combined_data$trap_data), crs = default_CRS)))

worldclim_tile("worldclim", var = "prec", res = 0.5, lon = tile_coords[1], lat = tile_coords[2], path = here("data", "geodata"))

precip_rast <- rast(here("data", "geodata", "wc2.1_tiles", "tile_30_wc2.1_30s_prec.tif"))

month_split <- date_set %>%
  group_by(unique_site, visit, site_easting, site_northing) %>%
  filter(date_set == min(date_set)) %>%
  mutate(month = month(date_set)) %>%
  st_as_sf(coords = c("site_easting", "site_northing"), crs = SL_UTM) %>%
  st_transform(crs = default_CRS) %>%
  group_by(month) %>%
  group_split()

month_split <- lapply(function(X) {
  
  # Extract month of interest
  month <- unique(X$month)
  # Subset monthly raster to month of interest
  precip_month <- precip_rast[[month]]
  # Convert spatial DF to vect for speed
  vect_X <- vect(X)
  # Append precipitation to month_split DF
  X$precipitation <- terra::extract(precip_month, vect_X)[, 2]
  
  return(X)
}, X = month_split)

# Moon phase --------------------------------------------------------------

month_split <- lapply(function(X) {
  
  date_trap <- ymd(unique(X$date_set))
  
  moon_fraction <- getMoonIllumination(date = date_trap) %>%
    select(date, fraction)
  
  X <- X %>%
    left_join(moon_fraction, by = c("date_set" = "date")) %>%
    rename(moon_fraction = fraction)
  
  return(X)
  
}, X = month_split)

rain_moon <- bind_rows(month_split) %>%
  tibble() %>%
  select(unique_site, visit, precipitation, moon_fraction)

# Detection covariates combined -------------------------------------------

detection_covariates <- left_join(combined_trap_nights, rain_moon, by = c("unique_site", "visit")) %>%
  arrange(unique_site)

write_rds(detection_covariates, here("data", "synthetic_data", "detection_covariates.rds"))


# Occurrence covariates ---------------------------------------------------
# The primary outcome is the effect of habitat type on occurrence we extract this from the site data

if(!file.exists(here("data", "synthetic_data", "occurrence_covariates.rds"))) {
  
  land_use <- synthetic_data$synthetic_sites %>%
    select(unique_site, landuse)
  
  duplicated_land_use <- land_use %>%
    group_by(unique_site) %>%
    mutate(n = n()) %>%
    filter(n >= 2)  %>%
    distinct(unique_site, landuse) %>%
    ungroup()
  
  land_use <- land_use %>% 
    filter(!unique_site %in% duplicated_land_use$unique_site) %>%
    bind_rows(duplicated_land_use) %>%
    arrange(unique_site)
  
  # Get the bounding box of the village and buffer it by 100m before downloading from OSM
  
  distance_from_building <- function(data = synthetic_data$synthetic_sites, village_name) {
    
    osm <- opq(st_as_sfc(st_bbox(data %>%
                                 filter(village == village_name) %>%
                                 st_as_sf(coords = c("site_easting", "site_northing"), crs = SL_UTM) %>%
                                 st_transform(crs = default_CRS)), crs = default_CRS) %>%
               st_buffer(dist = 100) %>%
               st_bbox()) %>%
      add_osm_feature(key = "building") %>%
      osmdata_sf()
    
    buildings <- osm$osm_polygons %>%
      st_transform(crs = SL_UTM) %>%
      st_union()
    
    sites <- synthetic_data$synthetic_sites %>%
      filter(village == village_name) %>%
      distinct(unique_site, site_easting, site_northing) %>%
      st_as_sf(coords = c("site_easting", "site_northing"), crs = SL_UTM) %>%
      mutate(distance_building = as.numeric(st_distance(., buildings)))
  }
  
  distance_building <- lapply(X = c("baiama", "lalehun", "lambayama", "seilama"), function(X) {
    distance_from_building(village_name = X)
  }) %>%
    bind_rows() %>%
    tibble() %>%
    distinct(unique_site, distance_building)
  
  # We also use the distance from the centre of the village site these are stored as coordinates
  village_coords <- tibble(village = c("baiama", "lalehun", "lambayama", "seilama"),
                           X = c(-11.268454, -11.0803, -11.198249, -11.193628469657279),
                           Y = c(7.83708, 8.197533, 7.854131, 8.122285428353395)) %>%
    st_as_sf(coords = c("X", "Y"), crs = default_CRS) %>%
    st_transform(crs = SL_UTM)
  
  distance_from_centre <- synthetic_data$synthetic_sites %>%
    distinct(unique_site, site_easting, site_northing) %>%
    st_as_sf(coords = c("site_easting", "site_northing"), crs = SL_UTM) %>%
    mutate(distance_centre = case_when(str_detect(unique_site, "baiama") ~ as.numeric(st_distance(., village_coords %>%
                                                                                                  filter(village == "baiama"))),
                                       str_detect(unique_site, "lalehun") ~ as.numeric(st_distance(., village_coords %>%
                                                                                                   filter(village == "lalehun"))),
                                       str_detect(unique_site, "lambayama") ~ as.numeric(st_distance(., village_coords %>%
                                                                                                     filter(village == "lambayama"))),
                                       str_detect(unique_site, "seilama") ~ as.numeric(st_distance(., village_coords %>%
                                                                                                   filter(village == "seilama"))),
                                       TRUE ~ as.numeric(NA))) %>%
    tibble() %>%
    distinct(unique_site, distance_centre)
  
  # Elevation will also be used as an occurrence covariate
  # This method is currently failing due to an invalid/expired certificate
  # elevation_3s(lon = tile_coords[1], lat = tile_coords[2], path = here("data", "geodata"))
  # Will use the elevatr package instead for this we need a data.frame with the centre of each village
  
  elevation_rast <- get_elev_raster(locations = village_coords %>%
                                    st_transform(crs = default_CRS), prj = default_CRS, z = 12) %>%
    rast()
  
  elevation_vect <- vect(synthetic_data$synthetic_sites %>%
                         distinct(unique_site, site_easting, site_northing) %>%
                         st_as_sf(coords = c("site_easting", "site_northing"), crs = SL_UTM) %>%
                         st_transform(crs = default_CRS))
  
  elevation_vect$elevation <- terra::extract(elevation_rast, elevation_vect)[, 2]
  
  elevation <- data.frame(elevation_vect)
  
  
  # Occurrence covariates combined ------------------------------------------
  
  occurrence_covariates <- synthetic_data$synthetic_sites %>%
    distinct(village, unique_site) %>%
    left_join(land_use, by = c("unique_site")) %>%
    left_join(distance_building, by = c("unique_site")) %>%
    left_join(distance_from_centre, by = c("unique_site")) %>%
    left_join(elevation, by = c("unique_site")) %>%
    arrange(unique_site)
  
  write_rds(occurrence_covariates, here("data", "synthetic_data", "occurrence_covariates.rds"))
  
} else {
  
  occurrence_covariates <- read_rds(here("data", "synthetic_data", "occurrence_covariates.rds"))
  
}

# Site coordinates --------------------------------------------------------

coords <- synthetic_data$synthetic_sites %>%
  distinct(unique_site, site_easting, site_northing) %>%
  arrange(unique_site)

coord_array <- array(data = c(coords$site_easting, coords$site_northing), dim = c(nrow(coords), 2), dimnames = list(c(coords$unique_site), c("X", "Y")))

write_rds(coord_array, here("data", "synthetic_data", "site_coords.rds"))


