if (!require("pacman")) install.packages("pacman")

pkgs =
  c("coda",
    "cowplot",
    "DescTools",
    "elevatr",
    "fastDummies",
    "flextable",
    "geodata",
    "ggridges",
    "ggspatial",
    "ggtext",
    "googledrive",
    "here",
    "lubridate",
    "mapview",
    "osmdata",
    "RColorBrewer",
    "RhpcBLASctl",
    "rmapshaper",
    "rnaturalearth",
    "rnaturalearthdata",
    "rosm",
    "sf",
    "spOccupancy",
    "stars",
    "suncalc",
    "terra",
    "tidygraph",
    "tidyterra",
    "tidyverse",
    "vegan"
  )

pacman::p_load(pkgs, character.only = T)

default_CRS <- "EPSG:4326"
SL_UTM <- "EPSG:32629"

village_palette <- RColorBrewer::brewer.pal(n = 4, name = "Set1")
names(village_palette) <-  c("Lalehun", "Seilama", "Lambayama", "Baiama")

landuse_palette <- c("#00913A", "#FEC44F", "#A13B9E")
names(landuse_palette) <- c("Forest", "Agriculture", "Village")

all_species_order <- c("Mastomys spp", "Rattus spp", "Mus musculus", "Crocidura spp", "Praomys spp", "Lophuromys spp", "Mus minutoides", "Lemniscomys spp", "Malacomys spp", "Hylomyscus spp",
                       "Hybomys spp", "Dasymys spp", "Gerbillinae spp", "Gerbilliscus spp")

species_order_plots <- c("Mastomys spp", "Rattus spp", "Mus musculus", "Crocidura spp", "Praomys spp", "Lophuromys spp", "Mus minutoides")

group_landuse_palette <- c("#00913a", "#FEC44F", "#F7A820", "#A13B9E", "#5407A6")
names(group_landuse_palette) <- c("Forest", "Agriculture - Rural", "Agriculture - Peri-urban", "Village - Rural", "Village - Peri-urban")
season <- tibble(visit = 1:8, season = c(rep("Dry", 2), rep("Rainy", 2), rep("Dry", 2),rep("Rainy", 2)))
