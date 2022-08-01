
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
         trap_number = as.numeric(str_split(as.character(trap_uid), "_", simplify = TRUE)[, 5]),
         trap_uid = factor(paste0(village, "_", visit, "_", grid_number, "_", trap_number))) %>%
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
  
  return(select_site)
  
}

sites_in_grid <- assign_traps_to_cells(sites)

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
# data will effectively be duplicated                                   #

duplicate_detections <- detections %>%
  mutate(visit = visit + 6,
         trap_uid = factor(paste0(village, "_", visit, "_", grid_number, "_", trap_number)))
  
duplicate_sites <- lapply(sites_in_grid, function(x) {
  
  x %>%
    mutate(visit = visit + 6,
           trap_uid = paste0(village, "_", visit, "_", grid_number, "_", trap_number))
  
})

synthetic_detections <- bind_rows(detections, duplicate_detections)

synthetic_sites <- bind_rows(bind_rows(sites_in_grid), bind_rows(duplicate_sites))

synthetic_data <- list(synthetic_detections = synthetic_detections,
                       synthetic_sites = synthetic_sites)

write_rds(synthetic_data, here("data", "synthetic_data", "synthetic_data.rds"))
