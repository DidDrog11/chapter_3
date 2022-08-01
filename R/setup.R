if (!require("pacman")) install.packages("pacman")

pkgs =
  c("coda",
    "cowplot",
    "googledrive",
    "here",
    "INLA",
    "lubridate",
    "mapview",
    "RhpcBLASctl",
    "sf",
    "spOccupancy",
    "stars",
    "tidyverse"
  )

pacman::p_load(pkgs, character.only = T)

default_CRS <- "EPSG:4326"
SL_UTM <- "EPSG:32629"
