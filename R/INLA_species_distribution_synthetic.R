# extractFixedINLA
extractFixedINLA = function(model, model_name="mod", transform=FALSE){
  ff = model$summary.fixed
  ff$param = row.names(ff)
  ff$param[ ff$param == "(Intercept)" ] = "Intercept"
  names(ff)[3:5] = c("lower", "median", "upper")
  if(transform == TRUE){
    ff[ 1:5 ] = exp(ff[ 1:5 ])
  }
  ff
}

# --------- set up "co-occurrence study" ----------

set.seed(200)

synthetic_data <- read_rds(here("data", "synthetic_data", "synthetic_data.rds"))

# site
site_match <- tibble(site_code = 1:length(unique(synthetic_data$synthetic_sites$unique_site)),
                     unique_site = unique(synthetic_data$synthetic_sites$unique_site))
site_codes <- site_match$site_code

nsites = length(site_codes)

# locations
site_locations <- site_match %>%
  left_join(synthetic_data$synthetic_sites %>%
              select(unique_site, site_easting, site_northing) %>%
              distinct(),
            by = c("unique_site"))
plot(site_locations$site_easting, site_locations$site_northing)

# simulate so lower probability of rodent B in sites with rodent A present
# initial probability of rodent A = 0.5; original probability of rodent B = 0.55 but decreased prob when A present
sites = site_locations %>%
  select(site_code, site_easting, site_northing)

# trap_uid to site_code
trap_to_site <- synthetic_data$synthetic_sites %>%
  select(unique_site, trap_uid) %>%
  distinct() %>%
  left_join(site_match, by = c("unique_site")) %>%
  left_join(synthetic_data$synthetic_detections %>%
              mutate(detection = 1) %>%
              select(trap_uid, species = clean_names, detection) %>%
              filter(str_detect(species, "mastomys|mus_musculus")) %>%
              pivot_wider(names_from = species, values_from = detection, values_fill = 0),
            by = "trap_uid")

site_detections <- trap_to_site %>%
  filter(!is.na(mastomys_spp)) %>%
  select(-trap_uid) %>%
  distinct() %>%
  group_by(site_code) %>%
  summarise(across(.cols = where(is.numeric), sum))

site_detections[site_detections == 2] <- 1

site_non_detections <- trap_to_site %>%
  filter(!site_code %in% site_detections$site_code) %>%
  mutate(across(.cols = where(is.numeric), replace_na, 0)) %>%
  group_by(site_code) %>%
  summarise(across(.cols = where(is.numeric), sum)) %>%
  slice_sample(n = nrow(site_detections) * 2)

all_sites_detection <- bind_rows(site_detections, site_non_detections) %>%
  arrange(site_code) %>%
  left_join(sites, by = "site_code")

# map
all_sites_detection %>% 
  reshape2::melt(id.vars = c("site_code", "site_easting", "site_northing")) %>%
  ggplot() + 
  geom_point(aes(x = site_easting, y = site_northing, colour = value, size = value)) + 
  facet_wrap(~variable) + 
  theme_bw()

# create df for modelling
dd = all_sites_detection %>% 
  reshape2::melt(id.vars = c("site_code", "site_easting", "site_northing")) %>%
  rename("species" = variable, 
         "presence" = value)

# add "co-occurrence" variable
co_occ_vars <- unique(paste0("co_", dd$species))

co_value <- list()

for(n in 1:length(co_occ_vars)) {
  
  co_name <- co_occ_vars[n]
  species_name <- str_remove(co_name, "^co_")
  
  co_value[[n]] <- ifelse(dd$site_code %in% dd$site_code[dd$species == species_name & dd$presence == 1], 1, 0)
  
}

names(co_value) <- co_occ_vars
dd <- dd %>% 
  bind_cols(do.call(cbind, co_value))

# =============== fit joint likelihood model ===============

# multiple linear predictors, one for each rodent, each has a binomial likelihood
# i.e.
# Y_a = Mu_A + A.co.B # where Y_a is prob of species A, mu_A is the intercept for A, A.co.B is the co-occurrence with B
# Y_b = Mu_B + B.co.A
# both inferred jointly 

# create joint likelihood data: response is a matrix

nsites = length(unique(dd$site_code))

m.dat = list(Y = matrix(NA, length(unique(dd$species)) * nsites, length(unique(dd$species))))

for(n in 1:length(unique(dd$species))) {
  # each column is a rodent species
  m.dat$Y[nsites * (n-1) + 1:nsites * 1, n] = dd$presence[dd$species == unique(dd$species)[n]]
  
}

# intercept terms for each species (only "switched on" for that linear predictor)
m.dat$mu.A = rep(1:0, each = nsites)
m.dat$mu.B = rep(0:1, each = nsites)

# co-occurrence covariates (you would need to extend this to a matrix of n*n co-occurrences for n multi species)
# these are NA (i.e. excluded) from the linear predictor for the same species
m.dat$A.co.B = c(dd$co_mastomys_spp[dd$species == "mus_musculus"], rep(NA, nsites))
m.dat$B.co.A = c(rep(NA, nsites), dd$co_mus_musculus[dd$species == "mastomys_spp"])

# specify model
# priors
hyper1.iid = list(theta = list(prior = "pc.prec", param = c(1,0.01)))
control.fixed1 = list(mean.intercept = 0, # prior mean for intercept
                      prec.intercept = 1, # prior precision for intercept
                      mean = 0, # prior mean for fixed effects
                      prec = 1)  # prior precision for fixed effects

### fitINLAModel: fit and return INLA model with specified formula and family

#' @param formx inla formula object; i.e. created using formula(y ~ x + f())
#' @param family likelihood
#' @param config boolean (default FALSE): set config in compute to TRUE for inla.posterior.sample()
#' @param verbose verbose reporting on or off? default FALSE
#' @param return.marginals specify whether model should save and return marginals

fitINLAModel = function(formx, data, family, verbose=FALSE, config=FALSE, return.marginals=FALSE, inla.mode="experimental"){
  return(
    inla(formx,
         verbose = verbose,
         data = data,
         #data = inla.stack.data(stack1, spde=spde), # i have commented out the stack for spde as currently not fitting one
         family=family,
         control.fixed = control.fixed1,
         control.predictor=list(#A=inla.stack.A(stack1),  # i have commented out the stack for spde as currently not fitting one
           compute=TRUE,
           link=1),
         control.compute=list(cpo=TRUE, waic=TRUE, dic=TRUE,
                              config=config,
                              return.marginals=return.marginals),
         control.inla = list(strategy='adaptive', # adaptive gaussian
                             cmin=0), # fixing Q factorisation issue https://groups.google.com/g/r-inla-discussion-group/c/hDboQsJ1Mls)
         inla.mode = inla.mode)
  )
}


# formula including separate parameters for when A co-occurs with B and when B co-occurs with A
# with some work this could also include multiple SPDEs (1 per species)
formula = formula(Y ~ 0 + mu.A + mu.B + A.co.B + B.co.A)
m1 = fitINLAModel(formx = formula, data = m.dat, family = c("binomial", "binomial"))

# nicely shows the antagonism
summary(m1)

# summary model - colour by species (i.e. which linear predictor is the parameter inferred within)
# shows nicely the antagonism between the species
extractFixedINLA(m1) %>%
  dplyr::mutate(species = c("mus_musculus", "mastomys_natalensis", "mus_musculus", "mastomys_natalensis")) %>%
  ggplot() + 
  geom_point(aes(param, mean, col=species), size=2) + 
  geom_linerange(aes(param, ymin=lower, ymax=upper, col=species)) + 
  geom_hline(yintercept=0) + 
  xlab("Parameter") + ylab("Estimate (posterior mean + 95% cred int)") + 
  theme_bw() + coord_flip() 

# then calculate actual prob of occurrence

params = extractFixedINLA(m1)
A_noB = inla.link.invlogit(params$mean[1]) # 0.09
B_noA = inla.link.invlogit(params$mean[2]) # 0.29
A_withB = inla.link.invlogit(params$mean[1] + params$mean[3]) # 0.04
B_withA = inla.link.invlogit(params$mean[2] + params$mean[4]) # 0.16


# Compare to mastomys and rattus ------------------------------------------
# trap_uid to site_code
trap_to_site <- synthetic_data$synthetic_sites %>%
  select(unique_site, trap_uid) %>%
  distinct() %>%
  left_join(site_match, by = c("unique_site")) %>%
  left_join(synthetic_data$synthetic_detections %>%
              mutate(detection = 1) %>%
              select(trap_uid, species = clean_names, detection) %>%
              filter(str_detect(species, "mastomys|rattus")) %>%
              pivot_wider(names_from = species, values_from = detection, values_fill = 0),
            by = "trap_uid")

site_detections <- trap_to_site %>%
  filter(!is.na(mastomys_spp)) %>%
  select(-trap_uid) %>%
  distinct() %>%
  group_by(site_code) %>%
  summarise(across(.cols = where(is.numeric), sum))

site_detections[site_detections == 2] <- 1

site_non_detections <- trap_to_site %>%
  filter(!site_code %in% site_detections$site_code) %>%
  mutate(across(.cols = where(is.numeric), replace_na, 0)) %>%
  group_by(site_code) %>%
  summarise(across(.cols = where(is.numeric), sum)) %>%
  slice_sample(n = nrow(site_detections) * 2)

all_sites_detection <- bind_rows(site_detections, site_non_detections) %>%
  arrange(site_code) %>%
  left_join(sites, by = "site_code")

# map
all_sites_detection %>% 
  reshape2::melt(id.vars = c("site_code", "site_easting", "site_northing")) %>%
  ggplot() + 
  geom_point(aes(x = site_easting, y = site_northing, colour = value, size = value)) + 
  facet_wrap(~variable) + 
  theme_bw()

# create df for modelling
dd = all_sites_detection %>% 
  reshape2::melt(id.vars = c("site_code", "site_easting", "site_northing")) %>%
  rename("species" = variable, 
         "presence" = value)

# add "co-occurrence" variable
co_occ_vars <- unique(paste0("co_", dd$species))

co_value <- list()

for(n in 1:length(co_occ_vars)) {
  
  co_name <- co_occ_vars[n]
  species_name <- str_remove(co_name, "^co_")
  
  co_value[[n]] <- ifelse(dd$site_code %in% dd$site_code[dd$species == species_name & dd$presence == 1], 1, 0)
  
}

names(co_value) <- co_occ_vars
dd <- dd %>% 
  bind_cols(do.call(cbind, co_value))

nsites = length(unique(dd$site_code))

m.dat = list(Y = matrix(NA, length(unique(dd$species)) * nsites, length(unique(dd$species))))

for(n in 1:length(unique(dd$species))) {
  # each column is a rodent species
  m.dat$Y[nsites * (n-1) + 1:nsites * 1, n] = dd$presence[dd$species == unique(dd$species)[n]]
  
}

# intercept terms for each species (only "switched on" for that linear predictor)
m.dat$mu.B = rep(1:0, each = nsites)
m.dat$mu.C = rep(0:1, each = nsites)

# co-occurrence covariates (you would need to extend this to a matrix of n*n co-occurrences for n multi species)
# these are NA (i.e. excluded) from the linear predictor for the same species
m.dat$B.co.C = c(dd$co_rattus_spp[dd$species == "mastomys_spp"], rep(NA, nsites))
m.dat$C.co.B = c(rep(NA, nsites), dd$co_mastomys_spp[dd$species == "rattus_spp"])

# formula including separate parameters for when A co-occurs with B and when B co-occurs with A
# with some work this could also include multiple SPDEs (1 per species)
formula = formula(Y ~ 0 + mu.B + mu.C + B.co.C + C.co.B)
m2 = fitINLAModel(formx = formula, data = m.dat, family = c("binomial", "binomial"))

# nicely shows the antagonism
summary(m2)

# summary model - colour by species (i.e. which linear predictor is the parameter inferred within)
# shows nicely the antagonism between the species
extractFixedINLA(m2) %>%
  dplyr::mutate(species = c("mastomys_natalensis", "rattus_rattus", "mastomys_natalensis", "rattus_rattus")) %>%
  ggplot() + 
  geom_point(aes(param, mean, col=species), size=2) + 
  geom_linerange(aes(param, ymin=lower, ymax=upper, col=species)) + 
  geom_hline(yintercept=0) + 
  xlab("Parameter") + ylab("Estimate (posterior mean + 95% cred int)") + 
  theme_bw() + coord_flip() 

# then calculate actual prob of occurrence

params = extractFixedINLA(m2)
B_noC = inla.link.invlogit(params$mean[1]) # 0.24
C_noB = inla.link.invlogit(params$mean[2]) # 0.13
B_withC = inla.link.invlogit(params$mean[1] + params$mean[3]) # 0.18
C_withB = inla.link.invlogit(params$mean[2] + params$mean[4]) # 0.09
