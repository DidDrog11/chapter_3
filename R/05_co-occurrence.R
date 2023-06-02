source(here::here("R", "00_setup.R"))
install.packages("ggstatsplot")
library(ggstatsplot)

# Load data ---------------------------------------------------------------

y_long <- read_rds(here("data", "observed_data", "y_long.rds"))
occ_covs <- read_rds(here("data", "observed_data", "occ_covariates.rds"))
det_covs <- read_rds(here("data", "observed_data", "det_covs_sp_subset.rds"))

# included species
sp_codes <- sort(unique(y_long$species))
N <- length(sp_codes)

# Load models -------------------------------------------------------------

final_model <- read_rds(here("data", "observed_model_output", "model_dev", "final_model.rds"))
final_ppc <- read_rds(here("data", "observed_model_output", "model_dev", "final_ppc.rds"))

# Observed rodents --------------------------------------------------------
observed <- y_long %>%
  left_join(occ_covs) %>%
  janitor::tabyl(species, landuse)

observed_stratified <- y_long %>%
  left_join(occ_covs) %>%
  janitor::tabyl(species, setting, landuse)

psi_list <- as.data.frame.table(final_model$psi.samples) %>%
  mutate(Site = as.integer(Var3),
         Species = factor(Var2, labels = sp_codes)) %>%
  select(Site, Species, Freq) %>%
  group_by(Site, Species) %>%
  group_split()

psi_species <- lapply(psi_list, function(x) {
  
  samples <- 1
  
  sampled <- tibble(Site = rep(unique(x$Site), samples),
                    Species = rep(unique(x$Species), samples),
                    Psi = median(x$Freq),
                    Psi_mean = mean(x$Freq))
  
  return(sampled)
  
}) %>%
  bind_rows()


# Correlations for co-occurrence ------------------------------------------

cooccurrence_plot <- function(data = psi_species, species_1 = "mastomys_spp", species_2 = "rattus_spp") {
  
  paired_df <- data %>%
    filter(Species %in% c(species_1, species_2)) %>%
    select(Site, Species, Psi) %>%
    pivot_wider(names_from = Species, values_from = Psi) %>%
    left_join(occ_covs, by = c("Site" = "site_code")) %>%
    mutate(group_landuse = factor(group_landuse, levels = c("forest", "agriculture - rural", "agriculture - peri-urban", "village - rural", "village - peri-urban"),
                                  labels = c("Forest", "Agriculture - Rural", "Agriculture - Peri-urban", "Village - Rural", "Village - Peri-urban")))
  
  correlation <- paired_df %>%
    ggplot() +
    geom_point(aes(x = !! sym(species_1), y = !! sym(species_2), colour = group_landuse)) +
    scale_colour_manual(values = group_landuse_palette) +
    theme_bw() +
    labs(x = str_to_sentence(str_replace_all(species_1, "_", " ")),
         y = str_to_sentence(str_replace_all(species_2, "_", " ")),
         colour = "Stratified landuse")
  
  sp_1 <- paired_df %>%
    pull(species_1)
  sp_2 <- paired_df %>%
    pull(species_2)
  
  correlation_test_combined <- tibble(species_1 = species_1,
                                      species_2 = species_2,
                                      spearman_rho = cor.test(x = sp_1, y = sp_2, alternative = "two.sided", method = "spearman")$estimate,
                                      spearman_p = cor.test(x = sp_1, y = sp_2, alternative = "two.sided", method = "spearman")$p.value)
  
  if(observed %>% filter(species == species_1) %>% pull(forest) > 0 & 
     observed %>% filter(species == species_2) %>% pull(forest) > 0) {
    
    sp_1 <- paired_df %>%
      filter(str_detect(group_landuse, "Forest")) %>%
      pull(species_1)
    sp_2 <- paired_df %>%
      filter(str_detect(group_landuse, "Forest")) %>%
      pull(species_2)
    
    correlation_test_forest <- tibble(landuse = "Forest",
                                      species_1 = species_1,
                                      species_2 = species_2,
                                      spearman_rho = cor.test(x = sp_1, y = sp_2, alternative = "two.sided", method = "spearman")$estimate,
                                      spearman_p = cor.test(x = sp_1, y = sp_2, alternative = "two.sided", method = "spearman")$p.value)
  } else { 
    
    correlation_test_forest <- tibble(landuse = "Forest",
                                      species_1 = species_1,
                                      species_2 = species_2,
                                      spearman_rho = NA,
                                      spearman_p = NA)
    
  }
  
  if(observed %>% filter(species == species_1) %>% pull(agriculture) > 0 & 
     observed %>% filter(species == species_2) %>% pull(agriculture) > 0) {
    
    sp_1 <- paired_df %>%
      filter(str_detect(group_landuse, "Agriculture")) %>%
      pull(species_1)
    sp_2 <- paired_df %>%
      filter(str_detect(group_landuse, "Agriculture")) %>%
      pull(species_2)
    
    correlation_test_agriculture <- tibble(landuse = "Agriculture",
                                           species_1 = species_1,
                                           species_2 = species_2,
                                           spearman_rho = cor.test(x = sp_1, y = sp_2, alternative = "two.sided", method = "spearman")$estimate,
                                           spearman_p = cor.test(x = sp_1, y = sp_2, alternative = "two.sided", method = "spearman")$p.value)
    
    
  } else {
    
    correlation_test_agriculture <- tibble(landuse = "Agriculture",
                                           species_1 = species_1,
                                           species_2 = species_2,
                                           spearman_rho = NA,
                                           spearman_p = NA)
    
  }
  
  if(observed %>% filter(species == species_1) %>% pull(village) > 0 & 
     observed %>% filter(species == species_2) %>% pull(village) > 0) {
    
    sp_1 <- paired_df %>%
      filter(str_detect(group_landuse, "Village")) %>%
      pull(species_1)
    sp_2 <- paired_df %>%
      filter(str_detect(group_landuse, "Village")) %>%
      pull(species_2)
    
    correlation_test_village <- tibble(landuse = "Village",
                                       species_1 = species_1,
                                       species_2 = species_2,
                                       spearman_rho = cor.test(x = sp_1, y = sp_2, alternative = "two.sided", method = "spearman")$estimate,
                                       spearman_p = cor.test(x = sp_1, y = sp_2, alternative = "two.sided", method = "spearman")$p.value)
    
  } else {
    
    correlation_test_village <- tibble(landuse = "Village",
                                       species_1 = species_1,
                                       species_2 = species_2,
                                       spearman_rho = NA,
                                       spearman_p = NA)
    
  }
  
  correlation_tests <- bind_rows(correlation_test_forest,
                                 correlation_test_agriculture,
                                 correlation_test_village)
  
  return(list(correlation_plot = correlation,
              correlation_test_combined = correlation_test_combined,
              correlation_tests = correlation_tests))
  
}

# Correlations by landuse type ---------------------------------------------
species_list <- list()
sp_codes_2 <- c("mastomys_spp", "rattus_spp", "mus_musculus", "crocidura_spp", "praomys_spp", "lophuromys_spp", "mus_minutoides")

associations_to_test <- expand_grid(species_1 = sp_codes_2, species_2 = sp_codes_2) %>%
  mutate(key = paste0(pmin(species_1, species_2), pmax(species_1, species_2), sep = "")) %>%
  distinct(key, .keep_all = TRUE) %>%
  select(species_1, species_2) %>%
  group_by(species_1) %>%
  group_split()

correlation_list <- vector(mode = "list", length = length(sp_codes_2))

for(n in 1:length(sp_codes_2)) {
  
  correlation_list[[n]] <- vector(mode = "list", length = nrow(associations_to_test[[n]]))
  
  species_list[[n]] <- associations_to_test[[n]]
  
  for(i in 1:nrow(species_list[[n]])) {
    
    sp_1 <- species_list[[n]][i,1] %>%
      pull()
    sp_2 <- species_list[[n]][i,2] %>%
      pull()
    
    correlation_list[[n]][[i]] <- cooccurrence_plot(species_1 = sp_1, species_2 = sp_2)
    
    
  }
}

species_order_plots <- c("Mastomys spp", "Rattus spp", "Mus musculus", "Crocidura spp", "Praomys spp", "Lophuromys spp", "Mus minutoides")

correlation_df <- lapply(correlation_list, function(x) {
  lapply(x, function(y) {
    
    y$correlation_tests
    
  })
}) %>%
  bind_rows() %>%
  filter(species_1 != species_2) %>%
  mutate(sig = case_when(spearman_p <= 0.005 ~ TRUE,
                         is.na(spearman_p) ~ NA,
                         TRUE ~ FALSE),
         species_1 = factor(species_1, levels = c("mastomys_spp", "rattus_spp", "mus_musculus", "crocidura_spp", "praomys_spp", "lophuromys_spp", "mus_minutoides"),
                            labels = species_order_plots),
         species_2 = factor(species_2, levels = c("mastomys_spp", "rattus_spp", "mus_musculus", "crocidura_spp", "praomys_spp", "lophuromys_spp", "mus_minutoides"),
                            labels = species_order_plots),
         species_1 = fct_rev(species_1),
         landuse = factor(landuse, levels = c("Forest", "Agriculture", "Village")),
         strength = cut(spearman_rho, breaks = c(-1, -0.8, -0.6, -0.4, -0.2, -0.05, 0.05, 0.2, 0.4, 0.6, 0.8, 1),
                        labels = c("Very strong -ve", "Strong -ve", "Moderate -ve", "Weak -ve", "Very weak -ve", "No correlation", "Very weak +ve", "Weak +ve", "Moderate +ve", "Strong +ve", "Very strong +ve")),
         correlation_coef = round(spearman_rho, 2),
         pos_cor = case_when(correlation_coef < 0 ~ FALSE,
                             is.na(correlation_coef) ~ NA,
                             TRUE ~ TRUE))

correlation_plot <- correlation_df %>%
  ggplot() +
  geom_tile(aes(x = species_2, y = species_1, fill = correlation_coef)) +
  geom_label(data = correlation_df %>%
               filter(sig == TRUE) %>%
               mutate(correlation_coef = paste0(correlation_coef, "*")), 
             aes(x = species_2, y = species_1, label = correlation_coef), fontface = "bold") +
  geom_label(data = correlation_df %>%
               filter(sig == FALSE), 
             aes(x = species_2, y = species_1, label = correlation_coef)) +
  scale_fill_gradient2(low = "darkred", high = "darkblue", na.value = "grey", limits = c(-1, 1),
                       breaks = c(1, 0.5, 0, -0.5, -1), labels = c("Strong +ve", "", "None", "", "Strong -ve")) +
  scale_x_discrete(drop = FALSE, guide = guide_axis(n.dodge = 2)) +
  scale_y_discrete(drop = FALSE) +
  facet_wrap(~ landuse, ncol = 1) +
  labs(fill = "Strength of correlation",
       x = "Species",
       y = "Species") +
  theme_bw()

save_plot(correlation_plot, filename = here("output", "Figure_5.png"), base_height = 12, base_width = 8)


# Rearranging plot for poster ---------------------------------------------
# Function taken from https://stackoverflow.com/a/58734961
install.packages("lemon")
library(lemon)
shift_legend3 <- function(p) {
  pnls <- cowplot::plot_to_gtable(p) %>% gtable::gtable_filter("panel") %>%
    with(setNames(grobs, layout$name)) %>% purrr::keep(~identical(.x,zeroGrob()))
  
  if( length(pnls) == 0 ) stop( "No empty facets in the plot" )
  
  lemon::reposition_legend( p, "center", panel=names(pnls) )
}

poster_correlation_plot <- shift_legend3(correlation_df %>%
                                           ggplot() +
                                           geom_tile(aes(x = species_2, y = species_1, fill = correlation_coef)) +
                                           geom_label(data = correlation_df %>%
                                                        filter(sig == TRUE) %>%
                                                        mutate(correlation_coef = paste0(correlation_coef, "*")), 
                                                      aes(x = species_2, y = species_1, label = correlation_coef), fontface = "bold") +
                                           geom_label(data = correlation_df %>%
                                                        filter(sig == FALSE), 
                                                      aes(x = species_2, y = species_1, label = correlation_coef)) +
                                           scale_fill_gradient2(low = "darkred", high = "darkblue", na.value = "grey", limits = c(-1, 1),
                                                                breaks = c(1, 0.5, 0, -0.5, -1), labels = c("Strong +ve", "", "None", "", "Strong -ve")) +
                                           scale_x_discrete(drop = FALSE, guide = guide_axis(n.dodge = 2)) +
                                           scale_y_discrete(drop = FALSE) +
                                           facet_wrap(~ landuse, ncol = 2) +
                                           labs(fill = "Strength of correlation",
                                                x = element_blank(),
                                                y = element_blank()) +
                                           theme_bw())

save_plot(poster_correlation_plot, filename = here("output", "Figure_5.svg"), base_height = 8, base_width = 10)
