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
site_codes <- sort(unique(synthetic_data$synthetic_sites$unique_site))

# number of species
N <- length(sp_codes)

# number of replicates
K <- max(synthetic_data$synthetic_sites$visit)

# number of sites
J <- length(site_codes)

if(!file.exists(here("data", "synthetic_data", "y.rds"))) {
  
  y = array(NA, dim = c(N, J, K), dimnames = list(sp_codes, site_codes, 1:K))
  
  
  for(j in 1:J) { # loop through sites
    for(k in 1:K) { # loop through replicates
      # extract data for current site/replicate combination
      curr_df <- y_long %>%
        filter(unique_site == site_codes[j],
               visit == k)
      # if plot j was sampled during replicate k curr.df will have at least 1 row, if not no rodent was observed
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
raw_det <- read_rds(here("data", "synthetic_data", "detection_covariates.rds"))


# Produce occurrence covariates -------------------------------------------


# Format site coordinates -------------------------------------------------
coords <- read_rds(here("data", "synthetic_data", "site_coords.rds"))

# Create list object ------------------------------------------------------




