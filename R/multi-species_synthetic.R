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

if(!file.exists(here("data", "synthetic_data", "det_covs.rds"))) {
  
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
  
} else {
  
  det_covs <- read_rds(here("data", "synthetic_data", "det_covs.rds"))
  
}

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

# Model structure 2
# Occurrence
occ_ms_formula_2 <- ~ landuse
# Detection
det_ms_formula_2 <- ~ scale(precipitation)

# Model structure 3
# Occurrence
occ_ms_formula_3 <- ~ landuse + village
# Detection
det_ms_formula_3 <- ~ scale(precipitation)

# Model structure 3b village is a random effect
# Occurrence
occ_ms_formula_3b <- ~ landuse + (1|village)
# Detection
det_ms_formula_3b <- ~ scale(precipitation)

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

if(!file.exists(here("data", "model_output", "intercept_only.rds"))) {
  
  out_ms_int <- msPGOcc(occ.formula = occ_ms_formula_int, 
                        det.formula = det_ms_formula, 
                        data = data_msom, 
                        inits = ms_inits, 
                        n.samples = 30000, 
                        priors = ms_priors, 
                        verbose = TRUE, 
                        n.report = 6000, 
                        n.burn = 10000,
                        n.thin = 50, 
                        n.chains = 3)
  
  write_rds(out_ms_int, here("data", "model_output", "intercept_only.rds"))
  
} else {
  
  out_ms_int <- read_rds(here("data", "model_output", "intercept_only.rds"))
  
}

summary(out_ms_int, level = "both")

waicOcc(out_ms_int)

ppc_ms_out_int <- ppcOcc(out_ms_int, 'chi-squared', group = 1)

summary(ppc_ms_out_int)

X_0 <- array(1, dim = c(13, 1))
pred_int <- predict(out_ms_int, X_0)

pred_df_int <- data.frame(species = sp_codes,
                          mean_psi = apply(pred_int$psi.0.samples, 2, mean),
                          sd_psi = apply(pred_int$psi.0.samples, 2, sd))

pred_df_int %>%
  ggplot() +
  geom_point(aes(x = species, y = mean_psi, colour = species)) +
  geom_errorbar(aes(x = species, ymin = mean_psi-sd_psi, ymax = mean_psi+sd_psi, colour = species)) +
  theme_bw() +
  labs(title = "Probability of occurrence - Intercept only",
       y = "Mean Psi",
       x = "Species",
       colour = element_blank())

# Intermediate models -----------------------------------------------------
# Run model 2
if(!file.exists(here("data", "model_output", "model_2.rds"))) {
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
  
  write_rds(out_ms_2, here("data", "model_output", "model_2.rds"))
  
} else {
  
  out_ms_2 <- read_rds(here("data", "model_output", "model_2.rds"))
  
}

summary(out_ms_2, level = "both")

waicOcc(out_ms_2)

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

# Run model 3
if(!file.exists(here("data", "model_output", "model_3.rds"))) {
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
  
  write_rds(out_ms_3, here("data", "model_output", "model_3.rds"))
  
} else {
  
  out_ms_3 <- read_rds(here("data", "model_output", "model_3.rds"))
  
}

summary(out_ms_3, level = "both")

MCMCplot(out_ms_3$beta.samples)
MCMCplot(out_ms_3$alpha.samples)

waicOcc(out_ms_3)

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

# Run model 3b
data_msom_f <- data_msom
data_msom_f$occ.covs$village[data_msom_f$occ.covs$village == "baiama"] <- 1
data_msom_f$occ.covs$village[data_msom_f$occ.covs$village == "lalehun"] <- 2
data_msom_f$occ.covs$village[data_msom_f$occ.covs$village == "lambayama"] <- 3
data_msom_f$occ.covs$village[data_msom_f$occ.covs$village == "seilama"] <- 4
data_msom_f$occ.covs$village <- as.numeric(data_msom_f$occ.covs$village)

if(!file.exists(here("data", "model_output", "model_3b.rds"))) {
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
  
  write_rds(out_ms_3b, here("data", "model_output", "model_3b.rds"))
  
} else {
  
  out_ms_3b <- read_rds(here("data", "model_output", "model_3b.rds"))
  
}

summary(out_ms_3b, level = "both")

waicOcc(out_ms_3b)

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

# Full model --------------------------------------------------------------
# Run model

if(!file.exists(here("data", "model_output", "full_model.rds"))) {
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
  
  write_rds(out_ms_1, here("data", "model_output", "full_model.rds"))
  
} else {
  
  out_ms_1 <- read_rds(here("data", "model_output", "full_model.rds"))
  
}

summary(out_ms_1, level = "both")

ppc_ms_out_1 <- ppcOcc(out_ms_1, 'chi-squared', group = 1)

summary(ppc_ms_out_1)

waicOcc(out_ms_1)

# Produce prediction dataframe --------------------------------------------
X_1 <- raw_occ %>%
  select(village, landuse, distance_building, distance_centre, elevation) %>%
  mutate(village_value = 1,
         habitat_value = 1) %>%
  pivot_wider(names_from = village, values_from = village_value, values_fill = 0) %>%
  select(-baiama) %>%
  pivot_wider(names_from = landuse, values_from = habitat_value, values_fill = 0) %>%
  select(-agriculture) %>%
  mutate(intercept = 1,
         distance_building = (distance_building - mean(data_msom$occ.covs$distance_building)) / sd(data_msom$occ.covs$distance_building),
         distance_centre = (distance_centre - mean(data_msom$occ.covs$distance_village)) / sd(data_msom$occ.covs$distance_village),
         elevation = (elevation - mean(data_msom$occ.covs$elevation)) / sd(data_msom$occ.covs$elevation)) %>%
  select(intercept, fallow, forest, village_inside, village_outside, lalehun, lambayama, seilama, distance_building, distance_centre, elevation)

colnames_X_1 <- names(X_1)

X_1 <- array(data = unlist(X_1), dim = c(nrow(X_1), ncol(X_1)))

pred_1 <- predict(out_ms_1, X_1)

interp_1 <- as.data.frame(X_1)
names(interp_1) <- colnames_X_1

interp_1 <- interp_1 %>%
  mutate(village = case_when(intercept == 1 & lalehun == 0 & lambayama == 0 & seilama == 0 ~ "Baiama",
                             lalehun == 1 ~ "Lalehun",
                             lambayama == 1 ~ "Lambayama",
                             seilama == 1 ~ "Seilama"),
         landuse = case_when(intercept == 1 & fallow == 0 & forest == 0 & village_inside == 0 & village_outside == 0 ~ "Agriculture",
                             fallow == 1 ~ "Fallow",
                             forest == 1 ~ "Forest",
                             village_inside == 1 ~ "Village_inside",
                             village_outside == 1 ~ "Village_outside"))

pred_df <- list()

for(i in 1:nrow(interp_1)) {
  
  pred_df[[i]] <- as.data.frame(pred_1$psi.0.samples[ , , i]) %>%
    mutate(across(.cols = everything(), round, 5))
  names(pred_df[[i]]) <- sp_codes
  
  pred_df[[i]] <- pred_df[[i]] %>%
    mutate(village = interp_1$village[i],
           landuse = interp_1$landuse[i])
}

a <- do.call(rbind, pred_df) %>%
  pivot_longer(cols = c(-village, -landuse), names_to = "species", values_to = "psi") %>%
  filter(!str_detect(species, "gerbill|hybomy|hylom|lemniscom|malacomy"))

mus_lambayama <- a %>%
  filter(village == "Lambayama" & species == "mus_musculus") %>%
  group_by(landuse)

mean_mus_lambayama <- mus_lambayama %>%
  summarise(mean = mean(psi),
            median = median(psi))

ci_mus_lambayama <- lapply(group_split(mus_lambayama), function(X) {
  ecdf(X$psi)})

inv_ecdf <- function(f){
  x <- environment(f)$x
  y <- environment(f)$y
  approxfun(y, x)
}

mus_inside <- inv_ecdf(ci_mus_lambayama[[3]])
mus_inside(0.05)
mus_inside(0.95)

mus_outside <- inv_ecdf(ci_mus_lambayama[[4]])
mus_outside(0.05)
mus_outside(0.95)

mus_others <- a %>%
  filter(village != "Lambayama" & species == "mus_musculus") %>%
  group_by(landuse)

mean_mus_others <- mus_others %>%
  summarise(mean = mean(psi),
            median = median(psi))

ci_mus_others <- lapply(group_split(mus_others), function(X) {
  ecdf(X$psi)})

mus_inside <- inv_ecdf(ci_mus_others[[4]])
mus_inside(0.05)
mus_inside(0.95)

mus_outside <- inv_ecdf(ci_mus_others[[5]])
mus_outside(0.05)
mus_outside(0.95)

mastomys_all <- a %>%
  filter(species == "mastomys_spp") %>%
  group_by(landuse)

mean_mastomys_all <- mastomys_all %>%
  summarise(mean = mean(psi),
            median = median(psi))

ci_mastomys_all <- lapply(group_split(mastomys_all), function(X) {
  ecdf(X$psi)})

mastomys_inside <- inv_ecdf(ci_mastomys_all[[4]])
mastomys_inside(0.05)
mastomys_inside(0.95)

mastomys_outside <- inv_ecdf(ci_mastomys_all[[5]])
mastomys_outside(0.05)
mastomys_outside(0.95)

mastomys_ag <- inv_ecdf(ci_mastomys_all[[1]])
mastomys_ag(0.05)
mastomys_ag(0.95)

mastomys_forest <- inv_ecdf(ci_mastomys_all[[3]])
mastomys_forest(0.05)
mastomys_forest(0.95)

combined_plot <- ggplot(a) +
  geom_density_ridges(aes(x = psi, y = species, fill = village),
                      rel_min_height = 0.01) +
  facet_wrap(~ landuse, scales = "free") +
  theme_bw()

save_plot(plot = combined_plot, filename = here("output", "model_1_plot.pdf"),
          base_width = 16,
          base_height = 10)
