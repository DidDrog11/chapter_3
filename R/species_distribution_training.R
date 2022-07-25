

# n.b. installed the latest "testing" version of the INLA package from inla website
library(INLA)
library(tidyverse)

# extractFixedINLA
extractFixedINLA = function(model, model_name="mod", transform=FALSE){
  ff = model$summary.fixed
  ff$param = row.names(ff)
  ff$param[ ff$param == "(Intercept)" ] = "Intercept"
  names(ff)[3:5] = c("lower", "median", "upper")
  if(transform == TRUE){
    ff[ 1:5 ] = exp(ff[ 1:5 ])
  }
  ff
}

# --------- set up "co-occurrence study" ----------

set.seed(200)


# 100 sites
nsites = 100

# locations
sites_lat = rnorm(nsites, 0, 0.5)
sites_lon = rnorm(nsites, 0, 0.75)
plot(sites_lat, sites_lon)

# simulate so lower probability of rodent B in sites with rodent A present
# initial probability of rodent A = 0.5; original probability of rodent B = 0.55 but decreased prob when A present
sites = data.frame(site_id = 1:nsites, lon = sites_lon, lat = sites_lat)
sites$rodent_A = rbinom(nrow(sites), size=1, prob=0.5)
sites$rodent_B = rbinom(nrow(sites), size=1, prob=0.55 - (0.4 * sites$rodent_A))

# map
sites %>% 
  reshape2::melt(id.vars = 1:3) %>%
  ggplot() + 
  geom_point(aes(lon, lat, col=value, size=value)) + 
  facet_wrap(~variable) + 
  theme_bw()

# create df for modelling
dd = sites %>%
  reshape2::melt(id.vars = 1:3) %>%
  dplyr::rename("species" = variable, 
                "presence" = value)

# add "co-occurrence" variable
dd$co_rodentA = ifelse(dd$site_id %in% dd$site_id[dd$species == "rodent_A" & dd$presence == 1], 1, 0)
dd$co_rodentB = ifelse(dd$site_id %in% dd$site_id[dd$species == "rodent_B" & dd$presence == 1], 1, 0)

# =============== fit joint likelihood model ===============

# multiple linear predictors, one for each rodent, each has a binomial likelihood
# i.e.
# Y_a = Mu_A + A.co.B # where Y_a is prob of species A, mu_A is the intercept for A, A.co.B is the co-occurrence with B
# Y_b = Mu_B + B.co.A
# both inferred jointly 

# create joint likelihood data: response is a matrix
# each column is a rodent species
m.dat = list(Y = matrix(NA, 2 * nsites, 2))
m.dat$Y[1:nsites, 1] = dd$presence[dd$species == "rodent_A"]
m.dat$Y[nsites+1:nsites, 2] = dd$presence[dd$species=="rodent_B"]

# intercept terms for each species (only "switched on" for that linear predictor)
m.dat$mu.A = rep(1:0, each=nsites)
m.dat$mu.B = rep(0:1, each=nsites)

# co-occurrence covariates (you would need to extend this to a matrix of n*n co-occurrences for n multi species)
# these are NA (i.e. excluded) from the linear predictor for the same species
m.dat$A.co.B = c(dd$co_rodentB[dd$species == "rodent_A"], rep(NA, nsites))
m.dat$B.co.A = c(rep(NA, nsites), dd$co_rodentA[dd$species == "rodent_B"])

# specify model
# priors
hyper1.iid = list(theta = list(prior = "pc.prec", param = c(1,0.01)))
control.fixed1 = list(mean.intercept = 0, # prior mean for intercept
                      prec.intercept = 1, # prior precision for intercept
                      mean = 0, # prior mean for fixed effects
                      prec = 1)  # prior precision for fixed effects

### fitINLAModel: fit and return INLA model with specified formula and family

#' @param formx inla formula object; i.e. created using formula(y ~ x + f())
#' @param family likelihood
#' @param config boolean (default FALSE): set config in compute to TRUE for inla.posterior.sample()
#' @param verbose verbose reporting on or off? default FALSE
#' @param return.marginals specify whether model should save and return marginals

fitINLAModel = function(formx, data, family, verbose=FALSE, config=FALSE, return.marginals=FALSE, inla.mode="experimental"){
  return(
    inla(formx,
         verbose = verbose,
         data = data,
         #data = inla.stack.data(stack1, spde=spde), # i have commented out the stack for spde as currently not fitting one
         family=family,
         control.fixed = control.fixed1,
         control.predictor=list(#A=inla.stack.A(stack1),  # i have commented out the stack for spde as currently not fitting one
           compute=TRUE,
           link=1),
         control.compute=list(cpo=TRUE, waic=TRUE, dic=TRUE,
                              config=config,
                              return.marginals=return.marginals),
         control.inla = list(strategy='adaptive', # adaptive gaussian
                             cmin=0), # fixing Q factorisation issue https://groups.google.com/g/r-inla-discussion-group/c/hDboQsJ1Mls)
         inla.mode = inla.mode)
  )
}


# formula including separate parameters for when A co-occurs with B and when B co-occurs with A
# with some work this could also include multiple SPDEs (1 per species)
formula = formula(Y ~ 0 + mu.A + mu.B + A.co.B + B.co.A)
m1 = fitINLAModel(formx = formula, data = m.dat, family = c("binomial", "binomial"))

# nicely shows the antagonism
summary(m1)

# summary model - colour by species (i.e. which linear predictor is the parameter inferred within)
# shows nicely the antagonism between the species
extractFixedINLA(m1) %>%
  dplyr::mutate(species = c("A", "B", "A", "B")) %>%
  ggplot() + 
  geom_point(aes(param, mean, col=species), size=2) + 
  geom_linerange(aes(param, ymin=lower, ymax=upper, col=species)) + 
  geom_hline(yintercept=0) + 
  xlab("Parameter") + ylab("Estimate (posterior mean + 95% cred int)") + 
  theme_bw() + coord_flip() 

# then calculate actual prob of occurrence

# close to original fixed params
params = extractFixedINLA(m1)
A_noB = inla.link.invlogit(params$mean[1]) # A in absence of B = estimated 0.61
B_noA = inla.link.invlogit(params$mean[2]) # B in absence of A = estimated 0.39
A_withB = inla.link.invlogit(params$mean[1] + params$mean[3]) # A in presence of B = estimated 0.17
B_withA = inla.link.invlogit(params$mean[2] + params$mean[4]) # B in presence of A = estimated 0.07
