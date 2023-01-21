extract_coeff <- function(object,
                          quantiles = c(0.025, 0.5, 0.975),
                          digits = max(3L, getOption("digits") - 3L), ...) {
  
  tmp.1 <- t(apply(object$beta.comm.samples, 2, 
                   function(x) c(mean(x), sd(x))))
  colnames(tmp.1) <- c("Mean", "SD")
  tmp <- t(apply(object$beta.comm.samples, 2, 
                 function(x) quantile(x, prob = quantiles)))
  diags <- matrix(c(object$rhat$beta.comm, round(object$ESS$beta.comm, 0)), ncol = 2)
  colnames(diags) <- c('Rhat', 'ESS')
  
  occurrence_means = data.frame(round(cbind(tmp.1, tmp, diags), digits)) %>%
    tibble() %>%
    mutate(Coefficient = rownames(data.frame(round(cbind(tmp.1, tmp, diags), digits))),
           `Odds Mean` = exp(Mean),
           `Odds SD` = exp(SD),
           `Odds 2.5%` = exp(X2.5.),
           `Odds 50%` = exp(X50.),
           `Odds 97.5%` = exp(X97.5.)) %>%
    select(Coefficient,
           Mean, SD, "2.5%" = X2.5., "50%" = X50., "97.5%" = X97.5.,
           `Odds Mean`, `Odds SD`, `Odds 2.5%`, `Odds 50%`, `Odds 97.5%`,
           Rhat, ESS)
  
  occurrence_means_plot =  occurrence_means %>%
    mutate(Coefficient = fct_inorder(Coefficient),
           Coefficient = factor(Coefficient, labels = c("Intercept", "Landuse - Forest", "Landuse - Village", "Village - Lalehun",
                                                        "Village - Lambayama", "Village - Seilama", "Distance from building", "Elevation")),
           Coefficient = fct_rev(Coefficient),
           Probability = `Odds Mean`/ (1 + `Odds Mean`),
           Probability_25 = `Odds 2.5%`/ (1 + `Odds 2.5%`),
           Probability_975 = `Odds 97.5%`/ (1 + `Odds 97.5%`)) %>%
    ggplot() + 
    geom_point(aes(x = Probability, y  =  Coefficient)) +
    geom_segment(aes(x = Probability_25, xend = Probability_975, y = Coefficient, yend = Coefficient)) +
    scale_x_continuous(limits = c(0, 1)) +
    theme_bw() +
    labs(title = "Community level probability of occurrence")
  
  tmp.1 <- t(apply(object$tau.sq.beta.samples, 2, 
                   function(x) c(mean(x), sd(x))))
  colnames(tmp.1) <- c("Mean", "SD")
  tmp <- t(apply(object$tau.sq.beta.samples, 2, 
                 function(x) quantile(x, prob = quantiles)))
  diags <- matrix(c(object$rhat$tau.sq.beta, round(object$ESS$tau.sq.beta, 0)), ncol = 2)
  colnames(diags) <- c('Rhat', 'ESS')
  
  occurrence_variances = data.frame(round(cbind(tmp.1, tmp, diags), digits)) %>%
    tibble() %>%
    mutate(Coefficient = rownames(data.frame(round(cbind(tmp.1, tmp, diags), digits))),
           `Odds Mean` = exp(Mean),
           `Odds SD` = exp(SD),
           `Odds 2.5%` = exp(X2.5.),
           `Odds 50%` = exp(X50.),
           `Odds 97.5%` = exp(X97.5.)) %>%
    select(Coefficient,
           Mean, SD, "2.5%" = X2.5., "50%" = X50., "97.5%" = X97.5.,
           `Odds Mean`, `Odds SD`, `Odds 2.5%`, `Odds 50%`, `Odds 97.5%`,
           Rhat, ESS)
  
  # Detection
  tmp.1 <- t(apply(object$alpha.comm.samples, 2, 
                   function(x) c(mean(x), sd(x))))
  colnames(tmp.1) <- c("Mean", "SD")
  tmp <- t(apply(object$alpha.comm.samples, 2, 
                 function(x) quantile(x, prob = quantiles)))
  diags <- matrix(c(object$rhat$alpha.comm, round(object$ESS$alpha.comm, 0)), ncol = 2)
  colnames(diags) <- c('Rhat', 'ESS')
  
  detection_means = data.frame(round(cbind(tmp.1, tmp, diags), digits)) %>%
    tibble() %>%
    mutate(Coefficient = rownames(data.frame(round(cbind(tmp.1, tmp, diags), digits))),
           `Odds Mean` = exp(Mean),
           `Odds SD` = exp(SD),
           `Odds 2.5%` = exp(X2.5.),
           `Odds 50%` = exp(X50.),
           `Odds 97.5%` = exp(X97.5.)) %>%
    select(Coefficient,
           Mean, SD, "2.5%" = X2.5., "50%" = X50., "97.5%" = X97.5.,
           `Odds Mean`, `Odds SD`, `Odds 2.5%`, `Odds 50%`, `Odds 97.5%`,
           Rhat, ESS)
  
  detection_means_plot =  detection_means %>%
    mutate(Coefficient = fct_inorder(Coefficient),
           Coefficient = factor(Coefficient, labels = c("Intercept", "Precipitation", "Moon fraction", "Trapnights")),
           Coefficient = fct_rev(Coefficient),
           Probability = `Odds Mean`/ (1 + `Odds Mean`),
           Probability_25 = `Odds 2.5%`/ (1 + `Odds 2.5%`),
           Probability_975 = `Odds 97.5%`/ (1 + `Odds 97.5%`)) %>%
    ggplot() + 
    geom_point(aes(x = Probability, y  =  Coefficient)) +
    geom_segment(aes(x = Probability_25, xend = Probability_975, y = Coefficient, yend = Coefficient)) +
    scale_x_continuous(limits = c(0, 1)) +
    theme_bw() +
    labs(title = "Community level probability of detection")
  
  tmp.1 <- t(apply(object$tau.sq.alpha.samples, 2, 
                   function(x) c(mean(x), sd(x))))
  colnames(tmp.1) <- c("Mean", "SD")
  tmp <- t(apply(object$tau.sq.alpha.samples, 2, 
                 function(x) quantile(x, prob = quantiles)))
  diags <- matrix(c(object$rhat$tau.sq.alpha, round(object$ESS$tau.sq.alpha, 0)), ncol = 2)
  colnames(diags) <- c('Rhat', 'ESS')
  
  detection_variances = data.frame(round(cbind(tmp.1, tmp, diags), digits)) %>%
    tibble() %>%
    mutate(Coefficient = rownames(data.frame(round(cbind(tmp.1, tmp, diags), digits))),
           `Odds Mean` = exp(Mean),
           `Odds SD` = exp(SD),
           `Odds 2.5%` = exp(X2.5.),
           `Odds 50%` = exp(X50.),
           `Odds 97.5%` = exp(X97.5.)) %>%
    select(Coefficient,
           Mean, SD, "2.5%" = X2.5., "50%" = X50., "97.5%" = X97.5.,
           `Odds Mean`, `Odds SD`, `Odds 2.5%`, `Odds 50%`, `Odds 97.5%`,
           Rhat, ESS)
  
  # Species levels
  tmp.1 <- t(apply(object$beta.samples, 2, 
                   function(x) c(mean(x), sd(x))))
  colnames(tmp.1) <- c("Mean", "SD")
  tmp <- t(apply(object$beta.samples, 2, 
                 function(x) quantile(x, prob = quantiles)))
  diags <- matrix(c(object$rhat$beta, round(object$ESS$beta, 0)), ncol = 2)
  colnames(diags) <- c('Rhat', 'ESS')
  
  species_occurrence = data.frame(round(cbind(tmp.1, tmp, diags), digits)) %>%
    tibble() %>%
    mutate(Coefficient = rownames(data.frame(round(cbind(tmp.1, tmp, diags), digits))),
           `Odds Mean` = exp(Mean),
           `Odds SD` = exp(SD),
           `Odds 2.5%` = exp(X2.5.),
           `Odds 50%` = exp(X50.),
           `Odds 97.5%` = exp(X97.5.)) %>%
    select(Coefficient,
           Mean, SD, "2.5%" = X2.5., "50%" = X50., "97.5%" = X97.5.,
           `Odds Mean`, `Odds SD`, `Odds 2.5%`, `Odds 50%`, `Odds 97.5%`,
           Rhat, ESS)
  
  species_means_plot =  species_occurrence %>%
    mutate(Species = str_split(Coefficient, "-", simplify = TRUE)[, 2],
           Species = fct_inorder(Species),
           Species = factor(Species, labels = c("Crocidura", "Lophuromys", "Mastomys", "Mus minutoides",
                                                "Mus musculus", "Praomys", "Rattus")),
           Species = fct_rev(Species),
           Coefficient = str_split(Coefficient, "-", simplify = TRUE)[, 1],
           Coefficient = fct_inorder(Coefficient),
           Coefficient = factor(Coefficient, labels = c("Intercept", "Landuse - Forest", "Landuse - Village", "Village - Lalehun",
                                                        "Village - Lambayama", "Village - Seilama", "Distance from building", "Elevation")),
           Probability = `Odds Mean`/ (1 + `Odds Mean`),
           Probability_25 = `Odds 2.5%`/ (1 + `Odds 2.5%`),
           Probability_975 = `Odds 97.5%`/ (1 + `Odds 97.5%`)) %>%
    ggplot() + 
    geom_point(aes(x = Probability, y  =  Species, colour = Species)) +
    geom_segment(aes(x = Probability_25, xend = Probability_975, y = Species, yend = Species, colour = Species)) +
    facet_wrap(~ Coefficient, ncol = 1) +
    scale_x_continuous(limits = c(0, 1)) +
    theme_bw() +
    labs(title = "Species level probability of occurrence")
  
  tmp.1 <- t(apply(object$alpha.samples, 2, 
                   function(x) c(mean(x), sd(x))))
  colnames(tmp.1) <- c("Mean", "SD")
  tmp <- t(apply(object$alpha.samples, 2, 
                 function(x) quantile(x, prob = quantiles)))
  diags <- matrix(c(object$rhat$alpha, round(object$ESS$alpha, 0)), ncol = 2)
  colnames(diags) <- c('Rhat', 'ESS')
  
  species_detection = data.frame(round(cbind(tmp.1, tmp, diags), digits)) %>%
    tibble() %>%
    mutate(Coefficient = rownames(data.frame(round(cbind(tmp.1, tmp, diags), digits))),
           `Odds Mean` = exp(Mean),
           `Odds SD` = exp(SD),
           `Odds 2.5%` = exp(X2.5.),
           `Odds 50%` = exp(X50.),
           `Odds 97.5%` = exp(X97.5.)) %>%
    select(Coefficient,
           Mean, SD, "2.5%" = X2.5., "50%" = X50., "97.5%" = X97.5.,
           `Odds Mean`, `Odds SD`, `Odds 2.5%`, `Odds 50%`, `Odds 97.5%`,
           Rhat, ESS)
  
  return(list(occurrence_means = occurrence_means,
              occurrence_variances = occurrence_variances,
              detection_means = detection_means,
              detection_variances = detection_variances,
              species_occurrence = species_occurrence,
              species_detection = species_detection))  
}
