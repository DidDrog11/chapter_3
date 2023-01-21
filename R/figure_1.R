source(here::here("R", "00_setup.R"))

# A map of Sierra Leone in Africa -----------------------------------------

world <- ne_countries(scale = "medium", returnclass = "sf")

africa <- world %>%
  filter(str_detect(continent, "Africa")) %>%
  select(admin) %>%
  mutate(fill = case_when(str_detect(admin, "Sierra Leone") ~ "black",
                            TRUE ~ "grey")) %>%
  ms_filter_islands(min_area = 1e10)

sl_bbox <- africa %>%
  filter(str_detect(admin, "Sierra Leone")) %>%
  st_bbox() %>%
  st_as_sfc()

africa_map <- ggplot() +
  geom_sf(data = africa, aes(fill = fill)) +
  geom_sf(data = sl_bbox, fill = NA, colour = "black", size = 4) +
  scale_fill_manual(values = c("grey", "white")) +
  guides(fill = "none") +
  theme_void()

sle_sf <- geodata::gadm(country = "SLE", level = 2, path = here("data", "geodata")) %>%
  st_as_sf()

fig_1_palette <- c(village_palette, "#FFFFFF")
names(fig_1_palette) <- c(names(village_palette)[1:4], "poi")

poi <- tibble(name = c("Freetown", "Bo", "Kenema", "Baiama", "Lalehun", "Lambayama", "Seilama"),
              cat = c("poi", "poi", "poi", "Baiama", "Lalehun", "Lambayama", "Seilama"),
              lat = c(8.48708717953912, 7.966794221623807, 7.876161956810467, 7.837529372181356, 8.197392257077409, 7.850593096948891, 8.12230048178563),
              lon = c(-13.2356631741985, -11.740987026160457, -11.190811585001954, -11.268407665149846, -11.08032958100431, -11.196939025872055, -11.1935976318381)) %>%
  st_as_sf(coords = c("lon", "lat"), crs = default_CRS)

sl_map <- ggplot() + 
  geom_sf(data = sle_sf, fill = "grey") +
  geom_sf(data = sle_sf %>%
            filter(str_detect(NAME_2, "Kenema")), fill = "#FFD580") +
  geom_sf(data = poi, size = 0.8) +
  ggrepel::geom_label_repel(data = poi, aes(label = name, geometry = geometry, fill = cat), stat = "sf_coordinates", min.segment.length = 0) +
  labs(x = element_blank(),
       y = element_blank()) +
  scale_fill_manual(values = fig_1_palette) +
  theme_bw() +
  ggspatial::annotation_north_arrow(style = north_arrow_minimal()) +
  ggspatial::annotation_scale(location = "br") +
  guides(fill = "none") +
  labs(title = "A)")

sl_inset_map <- ggdraw() +
  draw_plot(sl_map) +
  draw_plot(africa_map, x = 0.7, y = 0.65, width = 0.3, height = 0.3)

ggsave2(plot = sl_inset_map, filename = here("output", "Figure_1a.png"), dpi = 300, width = 7, height = 6)


# Trap timeline -----------------------------------------------------------
trap_data <- read_rds(here("data", "input", "combined_data.rds"))
trap_data <- trap_data$trap_data

timeline <- trap_data %>%
  tibble() %>%
  filter(village != "bambawo") %>%
  group_by(visit, village)  %>%
  mutate(tn = n(),
         date_set = min(date_set)) %>%
  select(date_set, visit, village, tn) %>%
  distinct() %>%
  ungroup() %>%
  mutate(date_set = case_when(date_set == as.Date("2022-01-18") ~ as.Date("2022-01-21"),
                              date_set == as.Date("2022-04-12") & village == "lalehun" ~ as.Date("2022-04-17"),
                              date_set == as.Date("2022-04-28") & village == "lalehun" ~ as.Date("2022-08-08"),
                              date_set == as.Date("2022-10-18") ~ as.Date("2022-08-24"),
                              date_set == as.Date("2022-10-29") & village == "seilama" ~ as.Date("2022-11-06"),
                              TRUE ~ date_set),
         village = str_to_title(village))

timeline_plot <- ggplot(timeline) +
  geom_rect(aes(xmin = as.Date("2022-05-01"), xmax = as.Date("2022-11-01"), ymin = -Inf, ymax = Inf), fill = "lightblue", alpha = 0.2) +
  geom_rect(aes(xmin = as.Date("2021-05-01"), xmax = as.Date("2021-11-01"), ymin = -Inf, ymax = Inf), fill = "lightblue", alpha = 0.2) +
  geom_rect(aes(xmin = date_set, xmax = date_set + 4, ymin = 0, ymax = tn, fill = village)) +
  scale_fill_manual(values = village_palette) +
  theme_bw() +
  labs(x = "Visit date",
       y = "Trap nights",
       fill = "Village") +
  labs(title = "F)")

save_plot(plot = timeline_plot, filename = here("output", "Figure_1f.png"), base_height = 3, base_width = 6)
# Trap locations ----------------------------------------------------------

fig_1_df <- read_rds(here("data", "observed_data", "fig_1_df.rds"))

bg <- list()

bg$baiama <- rast(here("data", "observed_data", "baiama_raster.tif"))
bg$lalehun <- rast(here("data", "observed_data", "lalehun_raster.tif"))
bg$lambayama <- rast(here("data", "observed_data", "lambayama_raster.tif"))
bg$seilama <- rast(here("data", "observed_data", "seilama_raster.tif"))

grids_plot <- list()

breaks <- list(baiama = list(x = c(-11.265, -11.255),
                             y = c(7.838, 7.83, 7.824)),
               lalehun = list(x = c(-11.081, -11.079),
                              y = c(8.194, 8.196, 8.199)),
               lambayama = list(x = c(-11.198, -11.195, -11.192),
                                y = c(7.8515, 7.85, 7.8485)),
               seilama = list(x = c(-11.198, -11.195),
                              y = c(8.124, 8.122)))

for(i in 1:length(fig_1_df))  {
  
  grids_plot[[i]] <- ggplot() + 
    geom_spatraster_rgb(data = bg[[i]]) +
    geom_sf(data = fig_1_df[[i]] %>%
              mutate(landuse = factor(str_to_title(landuse), levels = c("Forest", "Agriculture", "Village"))),
            aes(fill = tn, colour = tn)) +
    scale_colour_viridis_c(limits = c(0, 100), direction = -1) +
    scale_fill_viridis_c(limits = c(0, 100), direction = -1) +
    guides(colour = "none") +
    facet_wrap(~ landuse) +
    labs(fill = "Number Trap-Nights",
         title = str_to_title(unique(fig_1_df[[i]]$village))) +
    coord_sf(expand = FALSE) +
    scale_x_continuous(breaks = breaks[[i]]$x) +
    scale_y_continuous(breaks = breaks[[i]]$y) +
    theme_bw() +
    annotation_scale()
  
}

save_plot(plot = grids_plot[[1]] +
            labs(title = "B)") +
            guides(fill = "none"), filename = here("output", "Figure_1b.png"), base_height = 3, base_width = 7)
save_plot(plot = grids_plot[[2]] +
            labs(title = "C)") +
            guides(fill = "none"), filename = here("output", "Figure_1c.png"), base_height = 4, base_width = 6)
save_plot(plot = grids_plot[[3]] +
            labs(title = "D)") +
            guides(fill = "none"), filename = here("output", "Figure_1d.png"), base_height = 3, base_width = 7)
save_plot(plot = grids_plot[[4]] +
            labs(title = "E)") +
            guides(fill = "none"), filename = here("output", "Figure_1e.png"), base_height = 3, base_width = 7)

combined_grids <- plot_grid(plotlist = list(grids_plot[[4]] +
                                              theme(legend.position = "none"),
                                            grids_plot[[2]] +
                                              theme(legend.position = "none"),
                                            grids_plot[[1]] +
                                              theme(legend.position = "none"),
                                            grids_plot[[3]] +
                                              theme(legend.position = "none")),
                            ncol = 2)

save_plot(plot = combined_grids, filename = here("output", "grid_locations.png"), base_height = 8)
