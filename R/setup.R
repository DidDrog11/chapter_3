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
    "sf",
    "spOccupancy",
    "stars",
    "suncalc",
    "terra",
    "tidyverse"
  )

pacman::p_load(pkgs, character.only = T)

default_CRS <- "EPSG:4326"
SL_UTM <- "EPSG:32629"
