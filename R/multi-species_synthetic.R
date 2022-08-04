# First step is to make the data into the required format. The load_data.R  #
# script has read in the data and performed the initial reshaping of the    #
# data. The next steps are to do this for the multi-species model           #

synthetic_data <- read_rds(here("data", "synthetic_data", "synthetic_data.rds"))

# Produce y --------------------------------------

# produce a long format of detections with a single record per site, visit  #
# and species                                                               #

y_long <- synthetic_data$synthetic_detections %>%
  rename(species = clean_names) %>%
  left_join(synthetic_data$synthetic_sites %>%
              select(trap_uid, unique_site),
            by = "trap_uid") %>%
  group_by(unique_site, visit, species) %>%
  summarise(count = n())
  

# names of the "species" trapped
sp_codes <- sort(unique(y_long$species))
# names of the sites trapping occurred at
site_match <- tibble(site_code = 1:length(unique(synthetic_data$synthetic_sites$unique_site)),
                     unique_site = unique(synthetic_data$synthetic_sites$unique_site))
site_codes <- site_match$site_code

# associate sites with number in site_match
y_long <- y_long %>%
  left_join(site_match, by = "unique_site")

# number of species
N <- length(sp_codes)

# number of replicates
K <- max(synthetic_data$synthetic_sites$visit)

# number of sites
J <- length(site_codes)

if(!file.exists(here("data", "synthetic_data", "y.rds"))) {
  
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
      } else {
        y[1:N, j, k] <- 0
      }
    }
  }
  
  bai_lam_sites <- site_match$site_code[str_detect(site_match$unique_site, "baiama|lambayama")]
  # set the unsampled replicates in Baiama and Lambayama to NA
  y[1:N, bai_lam_sites, 9:10] <- NA
  
  
  str(y)
  
  # this produces our y array which can be saved to save time on repeat runs
  
  write_rds(y, here("data", "synthetic_data", "y.rds"))
  
} else {
  
  y <- read_rds(here("data", "synthetic_data", "y.rds"))
  
}

# summarise the total number of observations for each species
apply(y, 1, sum, na.rm = TRUE)


# Produce detection covariates --------------------------------------------
# here we add covariates that can impact the probability of detecting a rodent if it is present
raw_det <- read_rds(here("data", "synthetic_data", "detection_covariates.rds")) %>%
  left_join(site_match, by = c("unique_site"))

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

write_rds(det_covs, here("data", "synthetic_data", "det_covs.rds"))

# sites with NA trap_nights were not sampled during that replicate, we use this to recode y
for(n in 1:N) {
  for(k in 1:K) {
    y[n, 1:J, k][is.na(det_covs$trap_nights[, k])] <- NA
  }
}


# Produce occurrence covariates -------------------------------------------
raw_occ <- read_rds(here("data", "synthetic_data", "occurrence_covariates.rds")) %>%
  left_join(site_match, by = c("unique_site"))

landuse_mat <- matrix(NA, nrow = J, ncol = 1)
village_mat <- matrix(NA, nrow = J, ncol = 1)
building_mat <- matrix(NA, nrow = J, ncol = 1)
dist_village_mat <- matrix(NA, nrow = J, ncol = 1)
elevation_mat <- matrix(NA, nrow = J, ncol = 1)

for(j in 1:J) {
  landuse_mat[[j]] <- raw_occ$landuse[[j]]
  village_mat[[j]] <- raw_occ$village[[j]]
  building_mat[[j]] <- raw_occ$distance_building[[j]]
  dist_village_mat[[j]] <- raw_occ$distance_centre[[j]]
  elevation_mat[[j]] <- raw_occ$elevation[[j]]
}

occ_covs <- list(landuse = landuse_mat,
                 village = village_mat,
                 distance_building = building_mat,
                 distance_village = dist_village_mat,
                 elevation = elevation_mat)

write_rds(occ_covs, here("data", "synthetic_data", "occ_covs.rds"))

# Format site coordinates -------------------------------------------------
coords <- read_rds(here("data", "synthetic_data", "site_coords.rds"))

# Create list object ------------------------------------------------------

data_msom <- list(y = y,
                  occ.covs = occ_covs,
                  det.covs = det_covs,
                  coords = coords)


# Multi-species occupancy model -------------------------------------------

# Model structure intercept only
# Occurrence
occ_ms_formula_int <- ~ 1

# Model structure 1
# Occurrence
occ_ms_formula_1 <- ~ landuse + village + scale(distance_building) + scale(distance_village) + scale(elevation)
# Detection
det_ms_formula <- ~ scale(precipitation) + moon_fraction + scale(trap_nights)

# Initial values
ms_inits <- list(alpha.comm = 0,
                 beta.comm = 0,
                 beta = 0,
                 alpha = 0,
                 tau.sq.beta = 1,
                 tau.sq.alpha = 1,
                 z = apply(y, c(1, 2), max, na.rm = TRUE))

# Prior values
ms_priors <- list(beta.comm.normal = list(mean = 0, var = 2.72),
                  alpha.comm.normal = list(mean = 0, var = 2.72),
                  tau.sq.beta.ig = list(a = 0.1, b = 0.1),
                  alpha.sq.beta.ig = list(a = 0.1, b = 0.1))

# Intercept only model ----------------------------------------------------
# Run intercept only model
out_ms_int <- msPGOcc(occ.formula = occ_ms_formula_int, 
                    det.formula = det_ms_formula, 
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

write_rds(out_ms_int, here("data", "model_output", "intercept_only.rds"))

summary(out_ms_int, level = "both")

waicOcc(out_ms_int)

ppc_ms_out_int <- ppcOcc(out_ms_int, 'chi-squared', group = 1)

summary(ppc_ms_out_int)

X_0 <- array(1, dim = c(13, 1))
pred_int <- predict(out_ms_int, X.0)

pred_df_int <- data.frame(species = sp_codes,
                          mean_psi = apply(pred_int$psi.0.samples, 2, mean),
                          sd_psi = apply(pred_int$psi.0.samples, 2, sd))

pred_df_int %>%
  ggplot() +
  geom_point(aes(x = species, y = mean_psi)) +
  geom_errorbar(aes(x = species, ymin = mean_psi-sd_psi, ymax = mean_psi+sd_psi)) +
  theme_bw()

# Run model
out_ms_1 <- msPGOcc(occ.formula = occ_ms_formula_1, 
                    det.formula = det_ms_formula, 
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

summary(out_ms_1, level = "both")

ppc_ms_out_1 <- ppcOcc(out_ms_1, 'chi-squared', group = 1)

summary(ppc_ms_out_1)

waicOcc(out_ms_1)

