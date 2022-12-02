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
  left_join(site_match, by = c("site_id"))

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

occ_covs <- list(landuse = landuse_mat,
                 village = village_mat,
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
occ_ms_formula_2 <- ~ landuse + (1|village) + scale(elevation)
# Detection
det_ms_formula_2 <- ~ scale(precipitation) + moon_fraction + scale(trap_nights)

# Model structure 3 random intercept, random slope
# Occurrence
occ_ms_formula_3 <- ~ landuse + (landuse|village) + scale(elevation)
# Detection
det_ms_formula_3 <- ~ scale(precipitation) + moon_fraction + scale(trap_nights)

# Model structure 3b village is a random effect
# Occurrence
occ_ms_formula_3b <- ~ landuse + (landuse|village) + scale(distance_building) + scale(elevation)
# Detection
det_ms_formula_3b <- ~ scale(precipitation) + moon_fraction + scale(trap_nights)

# Model structure 3c quartile of population density is a random effect
# Occurrence
occ_ms_formula_3c <- ~ landuse + (landuse|population_q) + scale(distance_building) + scale(elevation)
# Detection
det_ms_formula_3c <- ~ scale(precipitation) + moon_fraction + scale(trap_nights)

# Model structure 3d population density is added as a scaled effect
# Occurrence
occ_ms_formula_3d <- ~ landuse + scale(distance_building) + scale(elevation) + scale(population)
# Detection
det_ms_formula_3d <- ~ scale(precipitation) + moon_fraction + scale(trap_nights)

# Model structure 4 (spatial)
# Occurrence
occ_ms_formula_4 <- ~ landuse + landuse*village + scale(elevation)
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

if(!file.exists(here("data", "observed_observed_model_output", "intercept_only_sp_subset.rds"))) {
  
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

waicOcc(out_ms_int)

if(!file.exists(here("data", "observed_model_output", "ppc_ms_out_int_sp_subset.rds"))) {
  
  ppc_ms_out_int <- ppcOcc(out_ms_int, 'chi-squared', group = 1)
  
  write_rds(ppc_ms_out_int, here("data", "observed_model_output", "ppc_ms_out_int_sp_subset.rds"))
  
} else {
  
  ppc_ms_out_int <- read_rds(here("data", "observed_model_output", "ppc_ms_out_int_sp_subset.rds"))
  
}

summary(ppc_ms_out_int)

## Model 1 --------------------------------------------------------------
# Run model takes ~25 mins

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

waicOcc(out_ms_1)


## Model 2 -----------------------------------------------------
# Run model 2 takes ~27 mins
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

waicOcc(out_ms_2)

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
                      n.samples = 5000, 
                      priors = ms_priors, 
                      n.omp.threads = 1, 
                      verbose = TRUE, 
                      n.report = 1000, 
                      n.burn = 1500,
                      n.thin = 50, 
                      n.chains = 3)
  
  write_rds(out_ms_3, here("data", "observed_model_output", "model_3_sp_subset.rds"))
  
} else {
  
  out_ms_3 <- read_rds(here("data", "observed_model_output", "model_3_sp_subset.rds"))
  
}

summary(out_ms_3, level = "both")

waicOcc(out_ms_3)

if(!file.exists(here("data", "observed_model_output", "ppc_ms_out_3_sp_subset.rds"))) {
  
  ppc_ms_out_3 <- ppcOcc(out_ms_3, 'chi-squared', group = 1)
  
  write_rds(ppc_ms_out_3, here("data", "observed_model_output", "ppc_ms_out_3_sp_subset.rds"))
  
} else {
  
  ppc_ms_out_3 <- read_rds(here("data", "observed_model_output", "ppc_ms_out_3_sp_subset.rds"))
  
}

ppc_ms_out_3 <- ppcOcc(out_ms_3, 'chi-squared', group = 1)

summary(ppc_ms_out_3)



## Model 3b ----------------------------------------------------------------

# Run model 3b
data_msom_f <- data_msom
data_msom_f$occ.covs$village[data_msom_f$occ.covs$village == "baiama"] <- 1
data_msom_f$occ.covs$village[data_msom_f$occ.covs$village == "lalehun"] <- 2
data_msom_f$occ.covs$village[data_msom_f$occ.covs$village == "lambayama"] <- 3
data_msom_f$occ.covs$village[data_msom_f$occ.covs$village == "seilama"] <- 4
data_msom_f$occ.covs$village <- as.numeric(data_msom_f$occ.covs$village)

if(!file.exists(here("data", "observed_model_output", "model_3b_sp_subset.rds"))) {
  
  out_ms_3b <- msPGOcc(occ.formula = occ_ms_formula_3b, 
                       det.formula = det_ms_formula_3b, 
                       data = data_msom_f, 
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

waicOcc(out_ms_3b)


if(!file.exists(here("data", "observed_model_output", "ppc_ms_out_3b_sp_subset.rds"))) {
  
  ppc_ms_out_3b <- ppcOcc(out_ms_3b, 'chi-squared', group = 1)
  
  write_rds(ppc_ms_out_3b, here("data", "observed_model_output", "ppc_ms_out_3b_sp_subset.rds"))
  
} else {
  
  ppc_ms_out_3b <- read_rds(here("data", "observed_model_output", "ppc_ms_out_3b_sp_subset.rds"))
  
}

summary(ppc_ms_out_3b)




## Model 4 (Spatial) -------------------------------------------------------

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
                        n.report = 100, 
                        n.burn = 2000,
                        n.thin = 20, 
                        n.chains = 3,
                        tuning = list(phi = 0.5))
  
  write_rds(out_ms_4, here("data", "observed_model_output", "model_4_sp_subset.rds"))
  
} else {
  
  out_ms_4 <- read_rds(here("data", "observed_model_output", "model_4_sp_subset.rds"))
  
}

summary(out_ms_4, level = "both")

waicOcc(out_ms_4)

if(!file.exists(here("data", "observed_model_output", "ppc_ms_out_4_sp_subset.rds"))) {
  
  ppc_ms_out_4 <- ppcOcc(out_ms_4, 'chi-squared', group = 1)
  
  write_rds(ppc_ms_out_4, here("data", "observed_model_output", "ppc_ms_out_4_sp_subset.rds"))
  
} else {
  
  ppc_ms_out_4 <- read_rds(here("data", "observed_model_output", "ppc_ms_out_4_sp_subset.rds"))
  
}

summary(ppc_ms_out_4)

predict(out_ms_4)

all_species <- tibble(model = c("out_ms_int", "out_ms_1", "out_ms_2", "out_ms_3", "out_ms_3b", "out_ms_4"),
                      waic = c(5666, 5297, 5291, 5289, 5056, 5215),
                      com_bpc = c(0.63, .61, 0.61, 0.59, 0.27, 0.60),
                      max_bpc = c(0.87, 0.92, 0.93, 0.91, 0.55, 0.93),
                      min_bpc = c(0.19, 0.15, 0.16, 0.12, 0, 0.17))


# Predictions -------------------------------------------------------------


## Intercept only ----------------------------------------------------------
results_int <- list()

for(i in 1:nrow(raw_occ)) {
  
  results_int[[i]] <- as.data.frame(out_ms_int$psi.samples[, , i])
  names(results_int[[i]]) <- sp_codes
  results_int[[i]] <- results_int[[i]] %>%
    pivot_longer(cols = everything(), names_to = "species", values_to = "occurrence") %>%
    group_by(species) %>%
    summarise(mean_occurrence = mean(occurrence),
              sd_occurrence = sd(occurrence))
  
  results_int[[i]] <- bind_cols(raw_occ[i, ], results_int[[i]])
  
}

results_int <- do.call(rbind, results_int)

ggplot(data = results_int, aes(x = mean_occurrence, y = species, fill = species)) +
  geom_density_ridges() +
  theme_bw() +
  labs(title = "Probability of occurrence - Intercept only",
       y = "Species",
       x = "Psi",
       colour = element_blank())

number_individuals_int <- list()

for(i in 1:nrow(raw_occ)) {
  
  number_individuals_int[[i]] <- as.data.frame(out_ms_int$z.samples[, , i])
  names(number_individuals_int[[i]]) <- sp_codes
  number_individuals_int[[i]] <- number_individuals_int[[i]] %>%
    pivot_longer(cols = everything(), names_to = "species", values_to = "presence") %>%
    group_by(species) %>%
    summarise(modelled_presence = sum(presence))
  
  number_individuals_int[[i]] <- bind_cols(raw_occ[i, ], number_individuals_int[[i]])
  
}

number_individuals_int <- do.call(rbind, number_individuals_int)

modelled_individuals_int <- number_individuals_int %>%
  group_by(species) %>%
  summarise(`Intercept only` = sum(modelled_presence)/1200)

compare_predicted <- left_join(observed_species, modelled_individuals_int, by = "species")

## Model 1 -----------------------------------------------------------------
results_m1 <- list()

for(i in 1:nrow(raw_occ)) {
  
  results_m1[[i]] <- as.data.frame(out_ms_1$psi.samples[, , i])
  names(results_m1[[i]]) <- sp_codes
  results_m1[[i]] <- results_m1[[i]] %>%
    pivot_longer(cols = everything(), names_to = "species", values_to = "occurrence")
  
  results_m1[[i]] <- bind_cols(raw_occ[i, ], results_m1[[i]]) %>%
    mutate(landuse = factor(landuse, levels = c("forest", "agriculture", "village"), labels = c("Forest", "Agriculture", "Village")),
           village = factor(village, levels = c("baiama", "lalehun", "seilama", "lambayama"), labels = c("Baiama", "Lalehun", "Seilama", "Lambayama")))
  
}

results_m1 <- do.call(rbind, results_m1)

summary_m1 <- results_m1 %>%
  group_by(species, landuse, village) %>%
  summarise(Median = round(median(occurrence), 4),
            Lower_quartile = round(quantile(occurrence, 0.25), 4),
            Upper_quartile = round(quantile(occurrence, 0.75), 4))

m1_figure <- ggplot(data = results_m1, aes(x = occurrence, y = species, fill = village)) +
  geom_boxplot(outlier.shape = NA) +
  theme_bw() +
  labs(title = "Probability of occurrence - Model 1",
       y = "Probability of occurrence",
       x = "Species",
       fill = "Village",
       colour = element_blank()) +
  facet_wrap(~ landuse, ncol = 1)

save_plot(plot = m1_figure, filename = here("output", "model_1_plot_sp_subset.png"),
          base_width = 8,
          base_height = 10)

number_individuals_1 <- list()

for(i in 1:nrow(raw_occ)) {
  
  number_individuals_1[[i]] <- as.data.frame(out_ms_1$z.samples[, , i])
  names(number_individuals_1[[i]]) <- sp_codes
  number_individuals_1[[i]] <- number_individuals_1[[i]] %>%
    pivot_longer(cols = everything(), names_to = "species", values_to = "presence") %>%
    group_by(species) %>%
    summarise(modelled_presence = sum(presence))
  
  number_individuals_1[[i]] <- bind_cols(raw_occ[i, ], number_individuals_1[[i]])
  
}

number_individuals_1 <- do.call(rbind, number_individuals_1)

modelled_individuals_1 <- number_individuals_1 %>%
  group_by(species) %>%
  summarise(`Model 1` = sum(modelled_presence)/1200)

compare_predicted <- left_join(compare_predicted, modelled_individuals_1, by = "species")

## Model 2 --------------------------------------------

X_2 <- matrix(c(1, 0, 0, 0, 0,
                1, 1, 0, 0, 0,
                1, 0, 1, 0, 0,
                1, 0, 0, 1, 0,
                1, 0, 0, 0, 1),
              nrow = 5, ncol =  5, byrow = TRUE)

pred_2 <- predict(out_ms_2, X_2)

interp_2 <- list()

for(i in 1:nrow(X_2)) {
  
  interp_2[[i]] <- pred_2$psi.0.samples[, , i]
  colnames(interp_2[[i]]) <- sp_codes
  
}

names(interp_2) <- c("Agriculture", "Fallow", "Forest", "Village_inside", "Village_outside")

interp_2 <- bind_rows(as.data.frame(interp_2[[1]]) %>%
                        mutate(landuse = "Agriculture") %>%
                        tibble(),
                      as.data.frame(interp_2[[2]]) %>%
                        mutate(landuse = "Fallow") %>%
                        tibble(),
                      as.data.frame(interp_2[[3]]) %>%
                        mutate(landuse = "Forest") %>%
                        tibble(),
                      as.data.frame(interp_2[[4]]) %>%
                        mutate(landuse = "Village_inside") %>%
                        tibble(),
                      as.data.frame(interp_2[[5]]) %>%
                        mutate(landuse = "Village_outside") %>%
                        tibble()) %>% 
  pivot_longer(cols = c(-landuse), values_to = "occupancy", names_to = "species")

ggplot(interp_2) +
  geom_point(aes(x = occupancy, y = species, colour = landuse), position = position_jitter()) +
  facet_wrap(~ landuse) +
  theme_bw()


## Model 3 --------------------------------------------

X_3 <- matrix(c(1, 0, 0, 0, 0, 0, 0, 0,
                1, 1, 0, 0, 0, 0, 0, 0,
                1, 0, 1, 0, 0, 0, 0, 0,
                1, 0, 0, 1, 0, 0, 0, 0,
                1, 0, 0, 0, 1, 0, 0, 0,
                1, 0, 0, 0, 0, 1, 0, 0,
                1, 1, 0, 0, 0, 1, 0, 0,
                1, 0, 1, 0, 0, 1, 0, 0,
                1, 0, 0, 1, 0, 1, 0, 0,
                1, 0, 0, 0, 1, 1, 0, 0,
                1, 0, 0, 0, 0, 0, 1, 0,
                1, 1, 0, 0, 0, 0, 1, 0,
                1, 0, 1, 0, 0, 0, 1, 0,
                1, 0, 0, 1, 0, 0, 1, 0,
                1, 0, 0, 0, 1, 0, 1, 0,
                1, 0, 0, 0, 0, 0, 0, 1,
                1, 1, 0, 0, 0, 0, 0, 1,
                1, 0, 1, 0, 0, 0, 0, 1,
                1, 0, 0, 1, 0, 0, 0, 1,
                1, 0, 0, 0, 1, 0, 0, 1),
              nrow = 20, ncol =  8, byrow = TRUE)


pred_3 <- predict(out_ms_3, X_3)

interp_3 <- list()

for(i in 1:nrow(X_3)) {
  
  interp_3[[i]] <- as.data.frame(pred_3$psi.0.samples[, , i])
  colnames(interp_3[[i]]) <- sp_codes
  
}

landuse_names <- c("Agriculture", "Fallow", "Forest", "Village_inside", "Village_outside")
village_names <- c("Baiama", "Lalehun", "Lambayama", "Seilama")

assign_names <- as.data.frame(X_3) %>%
  mutate(village = case_when(V1 == 1 & V6 == 0 & V7 == 0 & V8 == 0 ~ "Baiama",
                             V6 == 1 ~ "Lalehun",
                             V7 == 1 ~ "Lambayama",
                             V8 == 1 ~ "Seilama"),
         landuse = case_when(V1 == 1 & V2 == 0 & V3 == 0 & V4 == 0  & V5 == 0 ~ "Agriculture",
                             V2 == 1 ~ "Fallow",
                             V3 == 1 ~ "Forest",
                             V4 == 1 ~ "Village_inside",
                             V5 == 1 ~ "Village_outside"))

for(i in 1:length(interp_3)) {
  
  interp_3[[i]] <- interp_3[[i]] %>%
    mutate(village = assign_names$village[i],
           landuse = assign_names$landuse[i])
  
}

interp_3 <- do.call(rbind, interp_3) %>% 
  pivot_longer(cols = c(-landuse, -village), values_to = "occupancy", names_to = "species")

ggplot(interp_3) +
  geom_boxplot(aes(x = occupancy, y = species, colour = village)) +
  facet_wrap(~ landuse) +
  theme_bw()

## Model 3b --------------------------------------------

X_3b <- matrix(c(1, 0, 0, 0, 0, 1,
                 1, 1, 0, 0, 0, 1,
                 1, 0, 1, 0, 0, 1,
                 1, 0, 0, 1, 0, 1,
                 1, 0, 0, 0, 1, 1,
                 1, 0, 0, 0, 0, 2,
                 1, 1, 0, 0, 0, 2,
                 1, 0, 1, 0, 0, 2,
                 1, 0, 0, 1, 0, 2,
                 1, 0, 0, 0, 1, 2,
                 1, 0, 0, 0, 0, 3,
                 1, 1, 0, 0, 0, 3,
                 1, 0, 1, 0, 0, 3,
                 1, 0, 0, 1, 0, 3,
                 1, 0, 0, 0, 1, 3,
                 1, 0, 0, 0, 0, 4,
                 1, 1, 0, 0, 0, 4,
                 1, 0, 1, 0, 0, 4,
                 1, 0, 0, 1, 0, 4,
                 1, 0, 0, 0, 1, 4),
               nrow = 20, ncol =  6, byrow = TRUE)


# pred_3b <- predict(out_ms_3b, X_3b)
# 
# interp_3b <- list()
# 
# for(i in 1:nrow(X_3b)) {
#   
#   interp_3b[[i]] <- as.data.frame(pred_3b$psi.0.samples[, , i])
#   colnames(interp_3b[[i]]) <- sp_codes
#   
# }
# 
# landuse_names <- c("Agriculture", "Fallow", "Forest", "Village_inside", "Village_outside")
# village_names <- c("Baiama", "Lalehun", "Lambayama", "Seilama")
# 
# assign_names <- as.data.frame(X_3b) %>%
#   mutate(village = case_when(V1 == 1 & V6 == 0 & V7 == 0 & V8 == 0 ~ "Baiama",
#                              V6 == 1 ~ "Lalehun",
#                              V7 == 1 ~ "Lambayama",
#                              V8 == 1 ~ "Seilama"),
#          landuse = case_when(V1 == 1 & V2 == 0 & V3 == 0 & V4 == 0  & V5 == 0 ~ "Agriculture",
#                              V2 == 1 ~ "Fallow",
#                              V3 == 1 ~ "Forest",
#                              V4 == 1 ~ "Village_inside",
#                              V5 == 1 ~ "Village_outside"))
# 
# for(i in 1:length(interp_3b)) {
#   
#   interp_3b[[i]] <- interp_3b[[i]] %>%
#     mutate(village = assign_names$village[i],
#            landuse = assign_names$landuse[i])
#   
# }
# 
# interp_3b <- do.call(rbind, interp_3b) %>% 
#   pivot_longer(cols = c(-landuse, -village), values_to = "occupancy", names_to = "species")
# 
# ggplot(interp_3b) +
#   geom_boxplot(aes(x = occupancy, y = species, colour = village)) +
#   facet_wrap(~ landuse) +
#   theme_bw()


## Model 4 -----------------------------------------------------------------


