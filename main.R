# Load packages and project wide values
source(here::here("R", "00_setup.R"))

# Load cleaned data
source(here("R", "01_load_data.R"))

# Produce descriptive tables
source(here("R", "descriptive.R"))

# Produce network analysis
source(here("R", "network.R"))

# Produce multi-species occupancy models
source(here("R", "multi_species_models.R"))

# Produce species distributions
source(here("R", "species_distributions.R"))