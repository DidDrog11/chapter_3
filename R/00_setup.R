if (!require("pacman")) install.packages("pacman")

pkgs =
  c("coda",
    "cowplot",
    "DescTools",
    "elevatr",
    "fastDummies",
    "flextable",
    "geodata",
    "ggraph",
    "ggridges",
    "ggspatial",
    "googledrive",
    "here",
    "igraph",
    "INLA",
    "lubridate",
    "mapview",
    "modelr",
    "osmdata",
    "RColorBrewer",
    "RhpcBLASctl",
    "rosm",
    "sf",
    "spOccupancy",
    "stars",
    "suncalc",
    "terra",
    "tidygraph",
    "tidyterra",
    "tidyverse",
    "tmap",
    "vegan"
  )

pacman::p_load(pkgs, character.only = T)

default_CRS <- "EPSG:4326"
SL_UTM <- "EPSG:32629"

village_palette <- c("#7a0177", "#fec44f", "#005120", "#253494")
names(village_palette) <-  c("Lalehun", "Seilama", "Lambayama", "Baiama")
