if (!require("pacman")) install.packages("pacman")

pkgs =
  c("coda",
    "cowplot",
    "elevatr",
    "geodata",
    "googledrive",
    "here",
    "INLA",
    "lubridate",
    "mapview",
    "osmdata",
    "RhpcBLASctl",
    "rosm",
    "sf",
    "spOccupancy",
    "stars",
    "suncalc",
    "terra",
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
