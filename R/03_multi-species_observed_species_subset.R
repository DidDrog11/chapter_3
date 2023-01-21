source(here::here("R", "00_setup.R"))

# First step is to make the data into the required format. The load_data.R  #
# script has read in the data and performed the initial reshaping of the    #
# data. The next steps are to do this for the multi-species model           #

observed_data <- read_rds(here("data", "observed_data", "descriptive_data.rds"))

detections <- observed_data$detections %>%
  rename(trap_id = trap_uid,
         species = clean_names)

sites <- observed_data$sites_grids$select_site %>%
  bind_rows() %>%
  select(date_set, site_id = unique_site, site, landuse, site_easting, site_northing, village, visit, grid_number, trap_id = trap_uid, trap_easting, trap_northing, elevation)

trap_nights <- read_rds(here("data", "observed_data", "trap_nights.rds"))

data_for_4 <- write_rds(x = list(detections = detections, sites = sites), here("data", "data_for_export", "chapter_4_extract.rds"))

# Produce y --------------------------------------

# produce a long format of detections with a single record per site, visit  #
# and species                                                               #

y_long <- detections %>%
  left_join(sites %>%
              select(trap_id, site_id),
            by = "trap_id") %>%
  group_by(site_id, visit, species) %>%
  summarise(count = n())

included_species <- y_long %>%
  group_by(species) %>%
  summarise(count = sum(count)) %>%
  filter(count > 12) %>%
  pull(species)

y_long <- y_long %>%
  filter(species %in% included_species)

# names of the "species" trapped
sp_codes <- sort(unique(y_long$species))
# names of the sites trapping occurred at
site_match <- tibble(site_code = 1:length(unique(sites$site_id)),
                     site_id = unique(sites$site_id))
site_codes <- site_match$site_code

# define which sites were surveyed at each replicate
trap_nights_df <- site_match %>%
  left_join(trap_nights) %>%
  arrange(visit) %>%
  mutate(trap_nights = case_when(trap_nights > 1 ~ 1,
                                 TRUE ~ as.numeric(NA))) %>%
  pivot_wider(names_from = visit, values_from = trap_nights) %>%
  arrange(site_code) %>%
  select(-site_id)

# associate sites with number in site_match
y_long <- y_long %>%
  left_join(site_match, by = "site_id") %>%
  group_by(site_code) %>%
  mutate(non_0_site_code = cur_group_id())

write_rds(y_long, here("data", "observed_data", "y_long.rds"))

non_0_site_code <- unique(y_long$site_code)

# number of species
N <- length(sp_codes)

# number of replicates
K <- max(sites$visit)

# number of sites
J <- length(site_codes)

if(!file.exists(here("data", "observed_data", "y_sp_subset.rds"))) {
  
  y = array(NA, dim = c(N, J, K), dimnames = list(sp_codes, site_codes[1:J], 1:K))
  
  
  for(j in 1:J) { # loop through sites
    for(k in 1:K) { # loop through replicates
      # extract data for current site/replicate combination
      curr_df <- y_long %>%
        filter(site_code == site_codes[j],
               visit == k)
      # if plot j had a detection during replicate k curr_df will have at least 1 row, if not no rodent was observed
      if(nrow(curr_df) > 0) {
        # extract the species that were observed during this site/replicate
        curr_sp <- which(sp_codes %in% curr_df$species)
        # set value to 1 for species that were observed
        y[curr_sp, j, k] <- 1
        # and to 0 otherwise
        y[-curr_sp, j, k] <- 0
      } else 
        # if plot j was not sampled at replicate k (+1 as first column is the site code) set it to NA
        if(is.na(trap_nights_df[j, k+1])) {
          y[1:N, j, k] <- NA
        } else 
          # otherwise no rodent was trapped
        { 
          y[1:N, j, k] <- 0 
        }
    }
  }
  
  str(y)
  
  # this produces our y array which can be saved to save time on repeat runs
  
  write_rds(y, here("data", "observed_data", "y_sp_subset.rds"))
  
} else {
  
  y <- read_rds(here("data", "observed_data", "y_sp_subset.rds"))
  
}

# summarise the total number of observations for each species at distinct sites
observed_species <- tibble(species = names(apply(y, 1, sum, na.rm = TRUE)),
                           observed = apply(y, 1, sum, na.rm = TRUE),
                           check_number = y_long %>%
                             distinct(site_code, visit, species) %>%
                             group_by(species) %>%
                             summarise(n = n()) %>%
                             arrange(species) %>%
                             pull(n))

# Produce detection covariates --------------------------------------------
# here we add covariates that can impact the probability of detecting a rodent if it is present

if(!file.exists(here("data", "observed_data", "det_covs_sp_subset.rds"))) {
  
  raw_det <- read_rds(here("data", "observed_data", "detection_covariates.rds")) %>%
    left_join(site_match, by = c("site_id"))
  
  precip_mat <- matrix(NA, nrow = J, ncol = K)
  moon_mat <- matrix(NA, nrow = J, ncol = K)
  tn_mat <- matrix(NA, nrow = J, ncol = K)
  
  for (j in 1:J) { # Loop through sites
    for (k in 1:K) { # Loop through replicate surveys
      curr_vals <- raw_det %>%
        filter(site_code == site_codes[j], visit == k)
      # If the site was surveyed for the given replicate, 
      # extract the first date and time value. 
      if (nrow(curr_vals) > 0) {
        precip_mat[j, k] <- curr_vals$precipitation[1]
        moon_mat[j, k] <- curr_vals$moon_fraction[1] 
        tn_mat[j, k] <- curr_vals$trap_nights[1] 
      }
    } # k (replicates)
  } # j (sites) 
  
  det_covs <- list(precipitation = precip_mat,
                   moon_fraction = moon_mat,
                   trap_nights = tn_mat)
  
  write_rds(det_covs, here("data", "observed_data", "det_covs_sp_subset.rds"))
  
} else {
  
  det_covs <- read_rds(here("data", "observed_data", "det_covs_sp_subset.rds"))
  
}

# Produce occurrence covariates -------------------------------------------
raw_occ <- read_rds(here("data", "observed_data", "occurrence_covariates.rds")) %>%
  left_join(site_match, by = c("site_id")) %>%
  mutate(pop_quartile = case_when(str_detect(pop_quartile, "first") ~ 1,
                                  str_detect(pop_quartile, "second") ~ 2,
                                  str_detect(pop_quartile, "third") ~ 3,
                                  str_detect(pop_quartile, "fourth") ~ 4),
         setting = case_when(village == "lambayama" ~ "peri-urban",
                             TRUE ~ "rural"),
         group_landuse = case_when(village == "lambayama" & landuse == "agriculture" ~ "agriculture - peri-urban",
                                   village == "lambayama" & landuse == "village" ~ "village - peri-urban",
                                   village != "lambayama" & landuse == "agriculture" ~ "agriculture - rural",
                                   village != "lambayama" & landuse == "village" ~ "village - rural",
                                   TRUE ~ landuse))


landuse_mat <- matrix(NA, nrow = J, ncol = 1)
village_mat <- matrix(NA, nrow = J, ncol = 1)
setting_mat <- matrix(NA, nrow = J, ncol = 1)
group_landuse_mat <- matrix(NA, nrow = J, ncol = 1)
building_mat <- matrix(NA, nrow = J, ncol = 1)
dist_village_mat <- matrix(NA, nrow = J, ncol = 1)
elevation_mat <- matrix(NA, nrow = J, ncol = 1)
population_mat <- matrix(NA, nrow = J, ncol = 1)
population_q_mat <- matrix(NA, nrow = J, ncol = 1)

for(j in 1:J) {
  landuse_mat[[j]] <- raw_occ$landuse[[j]]
  village_mat[[j]] <- raw_occ$village[[j]]
  setting_mat[[j]] <- raw_occ$setting[[j]]
  group_landuse_mat[[j]] <- raw_occ$group_landuse[[j]]
  building_mat[[j]] <- raw_occ$distance_building[[j]]
  dist_village_mat[[j]] <- raw_occ$distance_centre[[j]]
  elevation_mat[[j]] <- raw_occ$elevation[[j]]
  population_mat[[j]] <- raw_occ$population[[j]]
  population_q_mat[[j]] <- raw_occ$pop_quartile[[j]]
}

occ_covs <- list(landuse = factor(landuse_mat, levels = c("forest", "agriculture", "village")),
                 village = factor(village_mat),
                 setting = factor(setting_mat, levels = c("rural", "peri-urban")),
                 group_landuse = factor(group_landuse_mat, levels = c("forest", "agriculture - rural", "agriculture - peri-urban",
                                                                      "village - rural", "village - peri-urban")),
                 distance_building = building_mat,
                 distance_village = dist_village_mat,
                 elevation = elevation_mat,
                 population = population_mat,
                 population_q = population_q_mat)

write_rds(occ_covs, here("data", "observed_data", "occ_covs.rds"))
write_rds(raw_occ, here("data", "observed_data", "occ_covariates.rds"))

# Format site coordinates -------------------------------------------------
coords <- read_rds(here("data", "observed_data", "site_coords.rds"))


# Spatial model data ------------------------------------------------------
# Distances between sites
dist_sites <- dist(coords)
# Exponential covariance model
cov_model <- "exponential"

n_factors = 1

lambda_inits <- matrix(0, N, n_factors)

diag(lambda_inits) <- 1

lambda_inits[lower.tri(lambda_inits)] <- rnorm(sum(lower.tri(lambda_inits)))

# Create list object ------------------------------------------------------

data_msom <- list(y = y,
                  occ.covs = occ_covs,
                  det.covs = det_covs,
                  coords = coords)

data_msom_spatial <- list(y = y,
                          occ.covs = as.data.frame(occ_covs),
                          det.covs = det_covs,
                          coords = coords)

# Multi-species occupancy model -------------------------------------------

# Model structure intercept only
# Occurrence
occ_ms_formula_int <- ~ 1

det_ms_formula_int <- ~ scale(precipitation) + moon_fraction + scale(trap_nights)

# Model structure 1
# Occurrence
occ_ms_formula_1 <- ~ landuse + village + scale(distance_building) + scale(elevation)
# Detection
det_ms_formula_1 <- ~ scale(precipitation) + moon_fraction + scale(trap_nights)

# Model structure 2 random intercept, fixed predictor in individual level
# Occurrence
occ_ms_formula_2 <- ~ landuse + (1|population_q) + scale(distance_building) + scale(elevation)
# Detection
det_ms_formula_2 <- ~ scale(precipitation) + moon_fraction + scale(trap_nights)

# Model structure 3 random intercept, random slope
# Occurrence
occ_ms_formula_3 <- ~ landuse + scale(distance_building) + scale(elevation)
# Detection
det_ms_formula_3 <- ~ scale(precipitation) + moon_fraction + scale(trap_nights)

# Model structure 3b village is a random effect
# Occurrence
occ_ms_formula_3b <- ~ landuse + scale(distance_building)
# Detection
det_ms_formula_3b <- ~ scale(precipitation) + moon_fraction + scale(trap_nights)

# Model structure 3c quartile of population density is a random effect
# Occurrence
occ_ms_formula_3c <- ~ landuse + scale(elevation)
# Detection
det_ms_formula_3c <- ~ scale(precipitation) + moon_fraction + scale(trap_nights)

# Model structure 3d population density is added as a scaled effect
# Occurrence
occ_ms_formula_3d <- ~ landuse + scale(distance_building) + scale(elevation) + scale(population)
# Detection
det_ms_formula_3d <- ~ scale(precipitation) + moon_fraction + scale(trap_nights)

# Model structure 3e landuse and village only
# Occurrence
occ_ms_formula_3e <- ~ landuse + village
# Detection
det_ms_formula_3e <- ~ scale(precipitation) + moon_fraction + scale(trap_nights)

# Model structure 3f landuse, village and elevation
# Occurrence
occ_ms_formula_3f <- ~ landuse + village + scale(elevation)
# Detection
det_ms_formula_3f <- ~ scale(precipitation) + moon_fraction + scale(trap_nights)

# Model structure 3g landuse, village and elevation
# Occurrence
occ_ms_formula_3g <- ~ landuse + village + scale(distance_building)
# Detection
det_ms_formula_3g <- ~ scale(precipitation) + moon_fraction + scale(trap_nights)

# Model structure 4 (spatial)
# Occurrence
occ_ms_formula_4 <- ~ landuse  + village + scale(distance_building) + scale(elevation)
# Detection
det_ms_formula_4 <- ~ scale(precipitation) + moon_fraction + scale(trap_nights)

# Initial values non spatial
ms_inits <- list(alpha.comm = 0,
                 beta.comm = 0,
                 beta = 0,
                 alpha = 0,
                 tau.sq.beta = 1,
                 tau.sq.alpha = 1,
                 z = apply(y, c(1, 2), max, na.rm = TRUE))

# Initial values spatial
ms_inits_spatial <- list(alpha.comm = 0,
                         beta.comm = 0,
                         beta = 0,
                         alpha = 0,
                         tau.sq.beta = 1,
                         tau.sq.alpha = 1,
                         lambda = lambda_inits, 
                         phi = 3 / mean(dist_sites),
                         z = apply(data_msom$y, c(1, 2), max, na.rm = TRUE))

# Prior values
ms_priors <- list(beta.comm.normal = list(mean = 0, var = 2.72),
                  alpha.comm.normal = list(mean = 0, var = 2.72),
                  tau.sq.beta.ig = list(a = 0.1, b = 0.1),
                  alpha.sq.beta.ig = list(a = 0.1, b = 0.1))

# Prior values spatial
min_dist <- min(dist_sites)
max_dist <- max(dist_sites)

ms_priors_spatial <- list(beta.comm.normal = list(mean = 0, var = 2.72),
                          alpha.comm.normal = list(mean = 0, var = 2.72), 
                          tau.sq.beta.ig = list(a = 0.1, b = 0.1), 
                          tau.sq.alpha.ig = list(a = 0.1, b = 0.1),
                          sigma.sq.ig = list(a = 2, b = 2), 
                          phi.unif = list(a = 3 / max_dist, b = 3 / min_dist))

## Intercept only model ----------------------------------------------------
# Run intercept only model takes ~ 25 mins
dir.create(here("data", "observed_model_output"))

if(!file.exists(here("data", "observed_model_output", "intercept_only_sp_subset.rds"))) {
  
  out_ms_int <- msPGOcc(occ.formula = occ_ms_formula_int, 
                        det.formula = det_ms_formula_int, 
                        data = data_msom, 
                        inits = ms_inits, 
                        n.samples = 30000, 
                        priors = ms_priors, 
                        verbose = TRUE, 
                        n.report = 6000, 
                        n.burn = 10000,
                        n.thin = 50, 
                        n.chains = 3)
  
  write_rds(out_ms_int, here("data", "observed_model_output", "intercept_only_sp_subset.rds"))
  
} else {
  
  out_ms_int <- read_rds(here("data", "observed_model_output", "intercept_only_sp_subset.rds"))
  
}

summary(out_ms_int, level = "both")

waic_int <- waicOcc(out_ms_int)

if(!file.exists(here("data", "observed_model_output", "ppc_ms_out_int_sp_subset.rds"))) {
  
  ppc_ms_out_int <- ppcOcc(out_ms_int, 'chi-squared', group = 1)
  
  write_rds(ppc_ms_out_int, here("data", "observed_model_output", "ppc_ms_out_int_sp_subset.rds"))
  
} else {
  
  ppc_ms_out_int <- read_rds(here("data", "observed_model_output", "ppc_ms_out_int_sp_subset.rds"))
  
}

summary(ppc_ms_out_int)

## Model 1 --------------------------------------------------------------
# Run model takes ~30 mins

if(!file.exists(here("data", "observed_model_output", "model_1_sp_subset.rds"))) {
  
  out_ms_1 <- msPGOcc(occ.formula = occ_ms_formula_1, 
                      det.formula = det_ms_formula_1, 
                      data = data_msom, 
                      inits = ms_inits, 
                      n.samples = 30000, 
                      priors = ms_priors, 
                      n.omp.threads = 1, 
                      verbose = TRUE, 
                      n.report = 6000, 
                      n.burn = 10000,
                      n.thin = 50, 
                      n.chains = 3)
  
  write_rds(out_ms_1, here("data", "observed_model_output", "model_1_sp_subset.rds"))
  
} else {
  
  out_ms_1 <- read_rds(here("data", "observed_model_output", "model_1_sp_subset.rds"))
  
}

summary(out_ms_1, level = "both")

if(!file.exists(here("data", "observed_model_output", "ppc_ms_out_1_sp_subset.rds"))) {
  
  ppc_ms_out_1 <- ppcOcc(out_ms_1, 'chi-squared', group = 1)
  
  write_rds(ppc_ms_out_1, here("data", "observed_model_output", "ppc_ms_out_1_sp_subset.rds"))
  
} else {
  
  ppc_ms_out_1 <- read_rds(here("data", "observed_model_output", "ppc_ms_out_1_sp_subset.rds"))
  
}

summary(ppc_ms_out_1)

waic_1 <- waicOcc(out_ms_1)


## Model 2 -----------------------------------------------------
# Run model 2 takes ~30 mins
if(!file.exists(here("data", "observed_model_output", "model_2_sp_subset.rds"))) {
  
  out_ms_2 <- msPGOcc(occ.formula = occ_ms_formula_2, 
                      det.formula = det_ms_formula_2, 
                      data = data_msom, 
                      inits = ms_inits, 
                      n.samples = 30000, 
                      priors = ms_priors, 
                      n.omp.threads = 1, 
                      verbose = TRUE, 
                      n.report = 6000, 
                      n.burn = 10000,
                      n.thin = 50, 
                      n.chains = 3)
  
  write_rds(out_ms_2, here("data", "observed_model_output", "model_2_sp_subset.rds"))
  
} else {
  
  out_ms_2 <- read_rds(here("data", "observed_model_output", "model_2_sp_subset.rds"))
  
}

summary(out_ms_2, level = "both")

waic_2 <- waicOcc(out_ms_2)

if(!file.exists(here("data", "observed_model_output", "ppc_ms_out_2_sp_subset.rds"))) {
  
  ppc_ms_out_2 <- ppcOcc(out_ms_2, 'chi-squared', group = 1)
  
  write_rds(ppc_ms_out_2, here("data", "observed_model_output", "ppc_ms_out_2_sp_subset.rds"))
  
} else {
  
  ppc_ms_out_2 <- read_rds(here("data", "observed_model_output", "ppc_ms_out_2_sp_subset.rds"))
  
}

summary(ppc_ms_out_2)



## Model 3 -----------------------------------------------------------------

if(!file.exists(here("data", "observed_model_output", "model_3_sp_subset.rds"))) {
  
  
  out_ms_3 <- msPGOcc(occ.formula = occ_ms_formula_3, 
                      det.formula = det_ms_formula_3, 
                      data = data_msom, 
                      inits = ms_inits,
                      n.samples = 30000, 
                      priors = ms_priors, 
                      n.omp.threads = 1, 
                      verbose = TRUE, 
                      n.report = 6000, 
                      n.burn = 10000,
                      n.thin = 50, 
                      n.chains = 3)
  
  write_rds(out_ms_3, here("data", "observed_model_output", "model_3_sp_subset.rds"))
  
} else {
  
  out_ms_3 <- read_rds(here("data", "observed_model_output", "model_3_sp_subset.rds"))
  
}

summary(out_ms_3, level = "both")

waic_3 <- waicOcc(out_ms_3)

if(!file.exists(here("data", "observed_model_output", "ppc_ms_out_3_sp_subset.rds"))) {
  
  ppc_ms_out_3 <- ppcOcc(out_ms_3, 'chi-squared', group = 1)
  
  write_rds(ppc_ms_out_3, here("data", "observed_model_output", "ppc_ms_out_3_sp_subset.rds"))
  
} else {
  
  ppc_ms_out_3 <- read_rds(here("data", "observed_model_output", "ppc_ms_out_3_sp_subset.rds"))
  
}

summary(ppc_ms_out_3)



## Model 3b ----------------------------------------------------------------

# Run model 3b

if(!file.exists(here("data", "observed_model_output", "model_3b_sp_subset.rds"))) {
  
  out_ms_3b <- msPGOcc(occ.formula = occ_ms_formula_3b, 
                       det.formula = det_ms_formula_3b, 
                       data = data_msom, 
                       inits = ms_inits, 
                       n.samples = 30000, 
                       priors = ms_priors, 
                       n.omp.threads = 1, 
                       verbose = TRUE, 
                       n.report = 6000, 
                       n.burn = 10000,
                       n.thin = 50, 
                       n.chains = 3)
  
  write_rds(out_ms_3b, here("data", "observed_model_output", "model_3b_sp_subset.rds"))
  
} else {
  
  out_ms_3b <- read_rds(here("data", "observed_model_output", "model_3b_sp_subset.rds"))
  
}

summary(out_ms_3b, level = "both")

waic_3b <- waicOcc(out_ms_3b)


if(!file.exists(here("data", "observed_model_output", "ppc_ms_out_3b_sp_subset.rds"))) {
  
  ppc_ms_out_3b <- ppcOcc(out_ms_3b, 'chi-squared', group = 1)
  
  write_rds(ppc_ms_out_3b, here("data", "observed_model_output", "ppc_ms_out_3b_sp_subset.rds"))
  
} else {
  
  ppc_ms_out_3b <- read_rds(here("data", "observed_model_output", "ppc_ms_out_3b_sp_subset.rds"))
  
}

summary(ppc_ms_out_3b)

## Model 3c ----------------------------------------------------------------

# Run model 3c

if(!file.exists(here("data", "observed_model_output", "model_3c_sp_subset.rds"))) {
  
  out_ms_3c <- msPGOcc(occ.formula = occ_ms_formula_3c, 
                       det.formula = det_ms_formula_3c, 
                       data = data_msom, 
                       inits = ms_inits, 
                       n.samples = 30000, 
                       priors = ms_priors, 
                       n.omp.threads = 1, 
                       verbose = TRUE, 
                       n.report = 6000, 
                       n.burn = 10000,
                       n.thin = 50, 
                       n.chains = 3)
  
  write_rds(out_ms_3c, here("data", "observed_model_output", "model_3c_sp_subset.rds"))
  
} else {
  
  out_ms_3c <- read_rds(here("data", "observed_model_output", "model_3c_sp_subset.rds"))
  
}

summary(out_ms_3c, level = "both")

waic_3c <- waicOcc(out_ms_3c)


if(!file.exists(here("data", "observed_model_output", "ppc_ms_out_3c_sp_subset.rds"))) {
  
  ppc_ms_out_3c <- ppcOcc(out_ms_3c, 'chi-squared', group = 1)
  
  write_rds(ppc_ms_out_3c, here("data", "observed_model_output", "ppc_ms_out_3c_sp_subset.rds"))
  
} else {
  
  ppc_ms_out_3c <- read_rds(here("data", "observed_model_output", "ppc_ms_out_3c_sp_subset.rds"))
  
}

summary(ppc_ms_out_3c)

## Model 3d ----------------------------------------------------------------

# Run model 3d

if(!file.exists(here("data", "observed_model_output", "model_3d_sp_subset.rds"))) {
  
  out_ms_3d <- msPGOcc(occ.formula = occ_ms_formula_3d, 
                       det.formula = det_ms_formula_3d, 
                       data = data_msom, 
                       inits = ms_inits, 
                       n.samples = 30000, 
                       priors = ms_priors, 
                       n.omp.threads = 1, 
                       verbose = TRUE, 
                       n.report = 6000, 
                       n.burn = 10000,
                       n.thin = 50, 
                       n.chains = 3)
  
  write_rds(out_ms_3d, here("data", "observed_model_output", "model_3d_sp_subset.rds"))
  
} else {
  
  out_ms_3d <- read_rds(here("data", "observed_model_output", "model_3d_sp_subset.rds"))
  
}

summary(out_ms_3d, level = "both")

waic_3d <- waicOcc(out_ms_3d)


if(!file.exists(here("data", "observed_model_output", "ppc_ms_out_3d_sp_subset.rds"))) {
  
  ppc_ms_out_3d <- ppcOcc(out_ms_3d, 'chi-squared', group = 1)
  
  write_rds(ppc_ms_out_3d, here("data", "observed_model_output", "ppc_ms_out_3d_sp_subset.rds"))
  
} else {
  
  ppc_ms_out_3d <- read_rds(here("data", "observed_model_output", "ppc_ms_out_3d_sp_subset.rds"))
  
}

summary(ppc_ms_out_3d)

## Model 3e ----------------------------------------------------------------

# Run model 3e

if(!file.exists(here("data", "observed_model_output", "model_3e_sp_subset.rds"))) {
  
  out_ms_3e <- msPGOcc(occ.formula = occ_ms_formula_3e, 
                       det.formula = det_ms_formula_3e, 
                       data = data_msom, 
                       inits = ms_inits, 
                       n.samples = 30000, 
                       priors = ms_priors, 
                       n.omp.threads = 1, 
                       verbose = TRUE, 
                       n.report = 6000, 
                       n.burn = 10000,
                       n.thin = 50, 
                       n.chains = 3)
  
  write_rds(out_ms_3e, here("data", "observed_model_output", "model_3e_sp_subset.rds"))
  
} else {
  
  out_ms_3e <- read_rds(here("data", "observed_model_output", "model_3e_sp_subset.rds"))
  
}

summary(out_ms_3e, level = "both")

waic_3e <- waicOcc(out_ms_3e)


if(!file.exists(here("data", "observed_model_output", "ppc_ms_out_3e_sp_subset.rds"))) {
  
  ppc_ms_out_3e <- ppcOcc(out_ms_3e, 'chi-squared', group = 1)
  
  write_rds(ppc_ms_out_3e, here("data", "observed_model_output", "ppc_ms_out_3e_sp_subset.rds"))
  
} else {
  
  ppc_ms_out_3e <- read_rds(here("data", "observed_model_output", "ppc_ms_out_3e_sp_subset.rds"))
  
}

summary(ppc_ms_out_3e)

## Model 3f ----------------------------------------------------------------

# Run model 3f

if(!file.exists(here("data", "observed_model_output", "model_3f_sp_subset.rds"))) {
  
  out_ms_3f <- msPGOcc(occ.formula = occ_ms_formula_3f, 
                       det.formula = det_ms_formula_3f, 
                       data = data_msom, 
                       inits = ms_inits, 
                       n.samples = 30000, 
                       priors = ms_priors, 
                       n.omp.threads = 1, 
                       verbose = TRUE, 
                       n.report = 6000, 
                       n.burn = 10000,
                       n.thin = 50, 
                       n.chains = 3)
  
  write_rds(out_ms_3f, here("data", "observed_model_output", "model_3f_sp_subset.rds"))
  
} else {
  
  out_ms_3f <- read_rds(here("data", "observed_model_output", "model_3f_sp_subset.rds"))
  
}

summary(out_ms_3f, level = "both")

waic_3f <- waicOcc(out_ms_3f)


if(!file.exists(here("data", "observed_model_output", "ppc_ms_out_3f_sp_subset.rds"))) {
  
  ppc_ms_out_3f <- ppcOcc(out_ms_3f, 'chi-squared', group = 1)
  
  write_rds(ppc_ms_out_3f, here("data", "observed_model_output", "ppc_ms_out_3f_sp_subset.rds"))
  
} else {
  
  ppc_ms_out_3f <- read_rds(here("data", "observed_model_output", "ppc_ms_out_3f_sp_subset.rds"))
  
}

summary(ppc_ms_out_3f)

## Model 3g ----------------------------------------------------------------

# Run model 3g

if(!file.exists(here("data", "observed_model_output", "model_3g_sp_subset.rds"))) {
  
  out_ms_3g <- msPGOcc(occ.formula = occ_ms_formula_3g, 
                       det.formula = det_ms_formula_3g, 
                       data = data_msom, 
                       inits = ms_inits, 
                       n.samples = 30000, 
                       priors = ms_priors, 
                       n.omp.threads = 1, 
                       verbose = TRUE, 
                       n.report = 6000, 
                       n.burn = 10000,
                       n.thin = 50, 
                       n.chains = 3)
  
  write_rds(out_ms_3g, here("data", "observed_model_output", "model_3g_sp_subset.rds"))
  
} else {
  
  out_ms_3g <- read_rds(here("data", "observed_model_output", "model_3g_sp_subset.rds"))
  
}

summary(out_ms_3g, level = "both")

waic_3g <- waicOcc(out_ms_3g)


if(!file.exists(here("data", "observed_model_output", "ppc_ms_out_3g_sp_subset.rds"))) {
  
  ppc_ms_out_3g <- ppcOcc(out_ms_3g, 'chi-squared', group = 1)
  
  write_rds(ppc_ms_out_3g, here("data", "observed_model_output", "ppc_ms_out_3g_sp_subset.rds"))
  
} else {
  
  ppc_ms_out_3g <- read_rds(here("data", "observed_model_output", "ppc_ms_out_3g_sp_subset.rds"))
  
}

summary(ppc_ms_out_3g)

## Model 4 (Spatial) -------------------------------------------------------
# Time to run ~ 36 minutes

if(!file.exists(here("data", "observed_model_output", "model_4_sp_subset.rds"))) {
  
  out_ms_4 <- sfMsPGOcc(occ.formula = occ_ms_formula_4, 
                        det.formula = det_ms_formula_4, 
                        data = data_msom_spatial, 
                        inits = ms_inits_spatial, 
                        n.batch = 500,
                        n.burn = 3000,
                        batch.length = 25,
                        n.thin = 20,
                        accept.rate = 0.43,
                        priors = ms_priors_spatial,
                        n.factors = n_factors,
                        cov.model = cov_model, 
                        n.omp.threads = 1, 
                        verbose = TRUE, 
                        NNGP = TRUE,
                        n.report = 50, 
                        n.chains = 4,
                        tuning = list(phi = 2))
  
  write_rds(out_ms_4, here("data", "observed_model_output", "model_4_sp_subset.rds"))
  
  gc()
  
} else {
  
  out_ms_4 <- read_rds(here("data", "observed_model_output", "model_4_sp_subset.rds"))
  
}

summary(out_ms_4, level = "both")

waic_4 <- waicOcc(out_ms_4)

if(!file.exists(here("data", "observed_model_output", "ppc_ms_out_4_sp_subset.rds"))) {
  
  gc()
  
  ppc_ms_out_4 <- ppcOcc(out_ms_4, 'chi-squared', group = 1)

  write_rds(ppc_ms_out_4, here("data", "observed_model_output", "ppc_ms_out_4_sp_subset.rds"))

} else {
  
  ppc_ms_out_4 <- read_rds(here("data", "observed_model_output", "ppc_ms_out_4_sp_subset.rds"))
  
}

summary(ppc_ms_out_4)

all_models <- tibble(model = c("out_ms_int", "out_ms_1", "out_ms_2", "out_ms_3",
                                "out_ms_3b", "out_ms_3c", "out_ms_3d", "out_ms_3e",
                               "out_ms_4"),
                      terms = as.character(c(occ_ms_formula_int, occ_ms_formula_1, occ_ms_formula_2, occ_ms_formula_3,
                                             occ_ms_formula_3b, occ_ms_formula_3c, occ_ms_formula_3d, occ_ms_formula_3e,
                                             occ_ms_formula_4)),
                      waic = c(waic_int[3], waic_1[3], waic_2[3], waic_3[3], waic_3b[3], waic_3c[3], waic_3d[3], waic_3e[3], waic_4[3]),
                      com_bpc = c(0.56, .61, 0.48, 0.51, 0.53, 0.52, 0.51, 0.51, 0.48),
                      max_bpc = c(0.8, 0.92, 0.72, 0.77, 0.78, 0.77, 0.74, 0.73, 0.7),
                      min_bpc = c(0.18, 0.15, 0.13, 0.11, 0.1, 0.11, 0.13, 0.13, 0.13))

