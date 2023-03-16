# Load packages and project wide values
source(here::here("R", "00_setup.R"))

## Need to go through this as visit and sites have been updated in the initial cleaning process #
## Job for another day

# Load cleaned data
source(here("R", "01_load_data.R"))

# Produce descriptive tables
source(here("R", "02_descriptive_observed.R"))

# Produce multi-species occupancy models
source(here("R", "03_multi-observed_species_subset.R"))

# Produce region and locations map

