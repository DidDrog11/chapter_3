if (!require("pacman")) install.packages("pacman")

pkgs =
  c("coda",
    "cowplot",
    "DescTools",
    "elevatr",
    "fastDummies",
    "flextable",
    "geodata",
    "ggmap",
    "ggnewscale",
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
    "tidyfast",
    "tidygraph",
    "tidyterra",
    "tidyverse",
    "vegan"
  )

pacman::p_load(pkgs, character.only = T)

default_CRS <- "EPSG:4326"
SL_UTM <- "EPSG:32629"

village_order <- c("baiama","lalehun", "lambayama", "seilama")

village_palette <- RColorBrewer::brewer.pal(n = 4, name = "Set1")
names(village_palette) <-  c("Lalehun", "Seilama", "Lambayama", "Baiama")

landuse_palette <- c("#00913A", "#FEC44F", "#A13B9E")
names(landuse_palette) <- c("Forest", "Agriculture", "Village")

all_species_order <- c("Mastomys natalensis", "Rattus rattus", "Mus musculus", "Crocidura olivieri", "Praomys rostratus", "Lophuromys sikapusi",
                       "Mus setulosus", "Crocidura buettikoferi", "Crocidura grandiceps", "Malacomys edwardsi", "Lemniscomys striatus",
                       "Hylomyscus simus", "Hybomys planifrons", "Mastomys erythroleucus", "Crocidura theresae", "Gerbilliscus guineae", "Dasymys rufulus")

species_order_plots <- c("Mastomys natalensis", "Rattus rattus", "Mus musculus", "Crocidura olivieri", "Praomys rostratus", "Lophuromys sikapusi",
                         "Mus setulosus", "Crocidura buettikoferi", "Crocidura grandiceps", "Malacomys edwardsi", "Lemniscomys striatus",
                         "Hylomyscus simus", "Hybomys planifrons", "Mastomys erythroleucus", "Crocidura theresae", "Gerbilliscus guineae", "Dasymys rufulus")

group_landuse_palette <- c("#00913a", "#FEC44F", "#F7A820", "#A13B9E", "#5407A6")
names(group_landuse_palette) <- c("Forest - Rural", "Agriculture - Rural", "Agriculture - Peri-urban", "Village - Rural", "Village - Peri-urban")
season <- tibble(visit = 1:10, season = c(rep("Dry", 2), rep("Rainy", 2), rep("Dry", 2), rep("Rainy", 2), rep("Dry", 2)))
