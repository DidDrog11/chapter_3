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
                                  str_detect(pop_quartile, "fourth") ~ 4))

landuse_mat <- matrix(NA, nrow = J, ncol = 1)
village_mat <- matrix(NA, nrow = J, ncol = 1)
building_mat <- matrix(NA, nrow = J, ncol = 1)
dist_village_mat <- matrix(NA, nrow = J, ncol = 1)
elevation_mat <- matrix(NA, nrow = J, ncol = 1)
population_mat <- matrix(NA, nrow = J, ncol = 1)
population_q_mat <- matrix(NA, nrow = J, ncol = 1)

for(j in 1:J) {
  landuse_mat[[j]] <- raw_occ$landuse[[j]]
  village_mat[[j]] <- raw_occ$village[[j]]
  building_mat[[j]] <- raw_occ$distance_building[[j]]
  dist_village_mat[[j]] <- raw_occ$distance_centre[[j]]
  elevation_mat[[j]] <- raw_occ$elevation[[j]]
  population_mat[[j]] <- raw_occ$population[[j]]
  population_q_mat[[j]] <- raw_occ$pop_quartile[[j]]
}

occ_covs <- list(landuse = factor(landuse_mat),
                 village = factor(village_mat),
                 distance_building = building_mat,
                 distance_village = dist_village_mat,
                 elevation = elevation_mat,
                 population = population_mat,
                 population_q = population_q_mat)

write_rds(occ_covs, here("data", "observed_data", "occ_covs.rds"))

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
  
  ppc_ms_out_3e <- ppcOcc(out_ms_3d, 'chi-squared', group = 1)
  
  write_rds(ppc_ms_out_3e, here("data", "observed_model_output", "ppc_ms_out_3e_sp_subset.rds"))
  
} else {
  
  ppc_ms_out_3e <- read_rds(here("data", "observed_model_output", "ppc_ms_out_3e_sp_subset.rds"))
  
}

summary(ppc_ms_out_3e)

## Model 4 (Spatial) -------------------------------------------------------
# Time to run ~ 36 minutes

if(!file.exists(here("data", "observed_model_output", "model_4_sp_subset.rds"))) {
  
  out_ms_4 <- sfMsPGOcc(occ.formula = occ_ms_formula_4, 
                        det.formula = det_ms_formula_4, 
                        data = data_msom_spatial, 
                        inits = ms_inits_spatial, 
                        n.batch = 400,
                        batch.length = 25,
                        accept.rate = 0.43,
                        priors = ms_priors_spatial,
                        n.factors = n_factors,
                        cov.model = cov_model, 
                        n.omp.threads = 1, 
                        verbose = TRUE, 
                        NNGP = TRUE,
                        n.report = 50, 
                        n.chains = 3,
                        tuning = list(phi = 0.5))
  
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



# Interpretation ----------------------------------------------------------

# Check mixing
plot(out_ms_4$beta.samples, density = FALSE)

plot(out_ms_4$alpha.samples, density = FALSE)

# Posterior predictive checks
# Bayesian p-value around 0.5 indicates adequate model fit, with values less than 0.1 or greater than 0.9 indicating poor fit.
summary(ppc_ms_out_4)

ppc.df <- tibble(ppv = c(rep(1:1200, times = 2)),
                 "crocidura" = c(ppc_ms_out_4$fit.y[ ,1], ppc_ms_out_4$fit.y.rep[, 1]),
                 "lophuromys" = c(ppc_ms_out_4$fit.y[ ,2], ppc_ms_out_4$fit.y.rep[, 2]),
                 "mastomys" = c(ppc_ms_out_4$fit.y[ ,3], ppc_ms_out_4$fit.y.rep[, 3]),
                 "minutoides" = c(ppc_ms_out_4$fit.y[ ,4], ppc_ms_out_4$fit.y.rep[, 4]),
                 "musculus" = c(ppc_ms_out_4$fit.y[ ,5], ppc_ms_out_4$fit.y.rep[, 5]),
                 "praomys" = c(ppc_ms_out_4$fit.y[ ,6], ppc_ms_out_4$fit.y.rep[, 6]),
                 "rattus" = c(ppc_ms_out_4$fit.y[ ,7], ppc_ms_out_4$fit.y.rep[, 7]),
                 fit = c(rep("True", times = 1200), rep("Fitted", times = 1200))) %>%
  arrange(ppv) %>%
  pivot_longer(cols = c("crocidura", "lophuromys", "mastomys", "minutoides", "musculus", "praomys", "rattus"), names_to = "Species",
               values_to = "Fitted_value") %>%
  pivot_wider(names_from = fit, values_from = "Fitted_value") %>%
  mutate(discrepancy = case_when(True > Fitted ~ "True > Fit",
                                 TRUE ~ "Fit > True"))

ggplot(ppc.df, aes(x = True, y = Fitted, colour = discrepancy)) +
  geom_point() +
  geom_abline() +
  facet_wrap(~ Species, scales = "free") +
  theme_bw()


# Probability of occurrence -----------------------------------------------

d0 <- as.data.frame.table(out_ms_4$psi.samples)
d1 <- d0 %>%
  mutate(Site = as.integer(Var3),
         Species = factor(Var2, labels = sp_codes)) %>%
  group_by(Site, Species) %>%
  summarise(Mean_psi = mean(Freq),
            SD_psi = SD(Freq))

# Model checks

check_model <- left_join(d1, y_long %>%
                           rename(Site = site_code,
                                  Species = species) %>%
                           group_by(Site, Species) %>%
                           summarise(observed = sum(count)),
                         by = c("Site", "Species")) %>%
  mutate(observed_bin = case_when(observed > 0 ~ 1,
                                  is.na(observed) ~ 0))

visualise_check <- ggplot(check_model) +
  geom_point(aes(x = Mean_psi, y = observed_bin)) +
  facet_wrap(~ Species, ncol = 1)

# Interpretation of species occurrence by landuse
species_order_plots <- c("Mastomys spp", "Rattus spp", "Mus musculus", "Crocidura spp", "Praomys spp", "Lophuromys spp", "Mus minutoides")

## Species occurrence by landuse -------------------------------------------

plot_m4 <- d1 %>%
  left_join(raw_occ, by = c("Site" = "site_code")) %>%
  mutate(landuse = factor(landuse, levels = c("forest", "agriculture", "village"), 
                          labels = c("Forest", "Agriculture", "Village")),
         Species = factor(str_to_sentence(str_replace_all(Species, "_", " ")),
                          levels = species_order_plots),
         village = str_to_sentence(village),
         peri_urban = case_when(village == "Lambayama" ~ "Peri-Urban",
                                TRUE ~ "Rural"))

landuse_plot <- plot_m4 %>%
  ggplot() +
  geom_jitter(aes(y = Mean_psi, x = landuse, colour = landuse), alpha = 0.2) + 
  geom_violin(aes(y = Mean_psi, x = landuse, fill = landuse)) + 
  facet_wrap(~ Species, nrow = 2) +
  scale_fill_manual(values = landuse_palette) +
  scale_colour_manual(values = landuse_palette) +
  guides(colour = "none") +
  theme_bw() +
  labs(y = "Probability of occurrence (ψ)",
       x = element_blank(),
       fill = element_blank())

save_plot(plot = landuse_plot, filename = here("output", "Figure_3.png"), base_width = 11, base_height = 9)

summaries <- plot_m4 %>%
  group_by(Species, village, landuse) %>%
  summarise(median_psi = median(Mean_psi),
            IQR_lower = IQR(Mean_psi))


## Species occurrence by landuse split on peri-urban/rural -----------------

urbanisation_by_landuse <- plot_m4 %>%
  ggplot() +
  geom_point(aes(y = Mean_psi, x = peri_urban, colour = landuse, fill = landuse), alpha = 0.2, position = position_jitterdodge(dodge.width = 0.9)) + 
  geom_violin(aes(y = Mean_psi, x = peri_urban, fill = landuse)) + 
  facet_wrap(~ Species, nrow = 2) +
  scale_fill_manual(values = landuse_palette) +
  scale_colour_manual(values = landuse_palette) +
  theme_bw() +
  guides(colour = "none") +
  labs(y = "Probability of occurrence (ψ)",
       x = element_blank(),
       fill = "Landuse",
       colour = element_blank())

save_plot(plot = urbanisation_by_landuse, filename = here("output", "Figure_4.png"), base_width = 8, base_height = 8)


## Species occurrence by distance from building and elevation ----------------------------

scaling_pred <- raw_occ %>%
  mutate(scaled_distance_building = scale(distance_building)[,1],
         scaled_elevation = scale(elevation)[,1],
         bin_forest = case_when(landuse == "forest" ~ 1,
                                TRUE ~ 0),
         bin_village = case_when(landuse == "village" ~ 1,
                                 TRUE ~ 0),
         bin_lal = case_when(village == "lalehun" ~ 1,
                             TRUE ~ 0),
         bin_lam = case_when(village == "lambayama" ~ 1,
                             TRUE ~ 0),
         bin_sei = case_when(village == "seilama" ~ 1,
                             TRUE ~ 0),
         intercept = 1)

# Marginal effect of distance_building
distance_plot <- plot_m4 %>%
  ungroup() %>%
  mutate(scaled_distance_building = scale(distance_building)[, 1]) %>%
  ggplot(aes(x = scaled_distance_building, y = Mean_psi, colour = landuse, group = landuse)) +
  geom_point() +
  geom_smooth() +
  facet_wrap(~ Species + village) +
  theme_bw()

save_plot(plot = distance_plot, filename = here("output", "distance_marginal.png"), base_width = 8, base_height = 8)

# Marginal effect of elevation
elevation_plot <- plot_m4 %>%
  ungroup() %>%
  mutate(scaled_elevation = scale(elevation)[, 1]) %>%
  ggplot(aes(x = scaled_elevation, y = Mean_psi, colour = landuse, group = landuse)) +
  geom_point() +
  geom_smooth() +
  facet_wrap(~ Species + village) +
  theme_bw()

save_plot(plot = elevation_plot, filename = here("output", "elevation_marginal.png"), base_width = 8, base_height = 8)

# Probability of detection ------------------------------------------------

d2 <- as.data.frame.table(out_ms_4$z.samples)
d3 <- d2 %>%
  mutate(Site = as.integer(Var3),
         Species = factor(Var2, labels = sp_codes)) %>%
  group_by(Site, Species) %>%
  summarise(Mean_z = mean(Freq),
            SD_z = SD(Freq))

precipitation <- scale(c(det_covs$precipitation))
min_precipitation <- min(precipitation, na.rm = TRUE)
max_precipitation <- max(precipitation, na.rm = TRUE)
pred_precipitation <- seq(from = min_precipitation, to = max_precipitation, length.out = 100)
mean_precipitation <- mean(pred_precipitation, na.rm = TRUE)

moon <- c(det_covs$moon_fraction)
min_moon <- min(moon, na.rm = TRUE)
max_moon <- max(moon, na.rm = TRUE)
pred_moon <- seq(from = min_moon, to = max_moon, length.out = 100)
mean_moon <- mean(moon, na.rm = TRUE)

tn <- scale(c(det_covs$trap_nights))
min_trap_nights <- min(tn, na.rm = TRUE)
max_trap_nights <- max(tn, na.rm = TRUE)
pred_trap_nights <- seq(from = min_trap_nights, to = max_trap_nights, length.out = 100)
mean_trap_nights <- mean(pred_trap_nights, na.rm = TRUE)


## Probability of detection - precipitation --------------------------------

precipitation_prediction <- cbind(1, pred_precipitation, mean_moon, mean_trap_nights)
out_detection_precipitation <- predict(out_ms_4, precipitation_prediction, type = "detection")

precipitation_estimates <- c(apply(out_detection_precipitation$p.0.samples, c(2, 3), mean))
plot_precipitation <- data.frame(detection_probability = precipitation_estimates,
                                 Species = rep(sp_codes, 100),
                                 scaled_precipitation = rep(pred_precipitation, each = N)) %>%
  mutate(Species = factor(str_to_sentence(str_replace_all(Species, "_", " ")),
                          levels = species_order_plots),
         Precipitation = scaled_precipitation * attr(precipitation, "scaled:scale") + attr(precipitation, 'scaled:center'))

precipitation_plot <- plot_precipitation %>%
  ggplot() +
  geom_line(aes(x = Precipitation, y = detection_probability)) +
  theme_bw() + 
  scale_y_continuous(limits = c(0, 1)) + 
  facet_wrap(~ Species) + 
  labs(x = 'Mean monthly rainfall (mm)', y = 'Detection Probability') 


save_plot(plot = precipitation_plot, filename = here("output", "precipitation_marginal.png"), base_width = 8, base_height = 8)

## Probability of detection - moon --------------------------------

moon_prediction <- cbind(1, mean_precipitation, pred_moon, mean_trap_nights)
out_detection_moon <- predict(out_ms_4, moon_prediction, type = "detection")

moon_estimates <- c(apply(out_detection_moon$p.0.samples, c(2, 3), mean))
plot_moon <- data.frame(detection_probability = moon_estimates,
                                 Species = rep(sp_codes, 100),
                                 Moon = rep(pred_moon, each = N)) %>%
  mutate(Species = factor(str_to_sentence(str_replace_all(Species, "_", " ")),
                          levels = species_order_plots))

moon_plot <- plot_moon %>%
  ggplot() +
  geom_line(aes(x = Moon, y = detection_probability)) +
  theme_bw() + 
  scale_y_continuous(limits = c(0, 1)) + 
  facet_wrap(~ Species) + 
  labs(x = 'Moon fraction', y = 'Detection Probability')

save_plot(plot = moon_plot, filename = here("output", "moon_marginal.png"), base_width = 8, base_height = 8)

## Probability of detection - trap night --------------------------------

trap_night_prediction <- cbind(1, mean_precipitation, mean_moon, pred_trap_nights)
out_detection_trap_night <- predict(out_ms_4, trap_night_prediction, type = "detection")

trap_night_estimates <- c(apply(out_detection_trap_night$p.0.samples, c(2, 3), mean))
plot_trap_night <- data.frame(detection_probability = trap_night_estimates,
                                 Species = rep(sp_codes, 100),
                                 scaled_trap_night = rep(pred_trap_nights, each = N)) %>%
  mutate(Species = factor(str_to_sentence(str_replace_all(Species, "_", " ")),
                          levels = species_order_plots),
         trap_night = scaled_trap_night * attr(tn, "scaled:scale") + attr(tn, 'scaled:center'))

trap_night_plot <- plot_trap_night %>%
  ggplot() +
  geom_line(aes(x = trap_night, y = detection_probability)) +
  theme_bw() + 
  scale_y_continuous(limits = c(0, 1)) + 
  facet_wrap(~ Species) + 
  labs(x = 'Trap nights', y = 'Detection Probability') 

save_plot(plot = trap_night_plot, filename = here("output", "tn_marginal.png"), base_width = 8, base_height = 8)