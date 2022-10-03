# Chapter 3

Planned analysis for chapter 3: Rodent trapping to explore rodent assemblage structure in Eastern Sierra Leone.

## Project structure

### Data

All data required to reproduce the analysis will be stored in `input`.

**Rodent data:** Cleaned rodent data will be obtained from `rodent_trapping.Rproj`. This code will therefore be separated into the cleaning project `rodent_trapping.Rproj` and this project for all analysis.

**Covariates:** Several covariates will be used, the data required for this will be imported and cleaned in this project.

### Code

Code can be run from the `main.R` script which will organise the functions and scripts used to conduct the analysis.

Scripts and functions will be contained in the `R` folder.

### Output

Plots and accompanying data will be produced and stored in an `output` folder at the point in which they are created. These will be subsequently imported into the report.

The manuscript or report will be stored in a `report` folder. This will be written as an `.Rmd` that can be produced as a`.docx` or `.pdf` file.

## Project aims

This chapter will describe the results of the rodent trapping study that was implemented in 4 villages in the Eastern Province of Sierra Leone. The project was designed to answer the following questions:

1.  What rodent species are prevalent at the study sites and how are species assemblages structured?
2.  Which species commonly co-occur and which species show evidence of competitive exclusion?
3.  Does rodent species diversity and richness differ between village site and habitat type?
4.  Do detection rates vary importantly over time?
5.  What is the potential impact of climate and land use change on *Lassa mammarenavirus* spillover risk?

### Future changes

  - [ ] Group agriculture and fallow, group village sites
  - [ ] Trial of rarely detected species as grouped by other
  - [ ] Run a spatial occupancy model
  - [ ] Remove network section which will be moved to chapter 4
