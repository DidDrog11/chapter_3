# Load packages and project wide values
source(here::here("R", "setup.R"))

# Load cleaned data
source(here("R", "load_data.R"))

# Produce descriptive tables
source(here("R", "descriptive.R"))

# Produce network analysis
source(here("R", "network.R"))

# Produced multi-occupancy models
source(here("R", "models.R"))

# Produce species distributions
source(here("R", "species_distributions.R"))