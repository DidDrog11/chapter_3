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

  - [x] Group agriculture and fallow, group village sites
  - [ ] Trial of rarely detected species as grouped by other - This doesn't make sense as some of the other species will be commensal/generalist/specialist
  - [x] Run a spatial occupancy model
  - [x] Remove network section which will be moved to chapter 4
  - [x] First component of analysis is species occurrence across landuse gradient
  - [x] Second component is species occurrence patterns across landuse gradient by urban scale
  - [ ] Should the model use human population density as a random effect rather than village as a fixed effect?
  - [ ] Do I want to predict probability of occupancy across Eastern Sierra Leone or limit it to the study areas? Prediction may be in effect replicating what I plan to do in chapter 5.
  - [ ] Does it make more sense to use the probability of occupancy data to explore co-occurrence as that way it is using the same modelling framework throughout?
  - [ ] Update introduction - add paragraph about species life histories generalist vs. specialist
  - [ ] Trapping protocol as supplementary, more specific detail into the methods for data collection  section
  - [x] Consistent definition of village site, grid site and landuse site
  - [x] Add species accumulation curves by village and by landuse setting
  - [x] Explore other options for co-occurrence
  
### Discussion points

  - [ ] 
