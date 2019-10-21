#' Fit BBS-style hierarchical model for estimating trends in point count survey data
#'
#' @param inputdata inputdata created by running \code{\link{setup_BBS_model}}
#' @param n.adapt number of iterations for adaptation, defaults to 500
#' @param n.update number of iterations for burn-in, defaults to 500
#' @param n.iter number of iterations for sampling, defauts to 1000
#' @param n.chains number of MCMC chains
#' @param overdispersion defaults to TRUE, to include a term in the model for
#' fitting overdispersion, as in BBS models
#' @param randomslope defaults to FALSE, because this is not included in BBS
#' models; change to TRUE to include a correlated random slope by transect, in
#' addition to a random intercept
#'
#' @return Returns a coda.samples object from rjags. Also displays summary
#' statistics from MCMCsummary and prints trace plots from MCMCplot.
#'
#' @export
#' @import rjags MCMCvis
#'

fit_BBS_model <- function(inputdata, n.adapt = 500, n.update = 500,
                          n.iter = 1000, n.chains = 3,
                          overdispersion = TRUE,
                          randomslope = FALSE,
                          ...) {
  if (randomslope == TRUE) {
    if (overdispersion == TRUE) {
      vars = c('mu.intercept', 'mu.slope', 'sigma.intercept', 'sigma.slope',
               'sigma.year', 'sigma.point', 'sigma.overdispersion',
               'slope', 'intercept',
               'abund', 'abund_transect', 'count_pred', 'count_pred_transect')

      modelstring = "
  model{
  ## ecological model
  for (i in 1:length(observed)){
    log(lambda[i]) = intercept[transect[i]] + slope[transect[i]] * zyear[i] +
      year_effect[year[i]] + point_effect[point[i]] + overdispersion[i]

    observed[i] ~ dpois(lambda[i])

    # overdisperson
    overdispersion[i] ~ dnorm(0, tau.overdispersion)
  }

  ## transect intercept and slope
  for (j in 1:ntransects) {
    intercept[j] <- B[j, 1]  # group level intercept
    slope[j]  <- B[j, 2]  # group level slope
    B[j, 1:2] ~ dmnorm(B.hat[j, 1:2], Tau.B)
    B.hat[j, 1] <- mu.intercept  # required by JAGS syntax
    B.hat[j, 2] <- mu.slope   # required by JAGS syntax
  }

  ## additional random intercepts for:
  ## - year
  for (t in 1:nyears){
    year_effect[t] ~ dnorm(0, tau.year)
  }
  ## - point (not nested within transect)
  for (k in 1:npoints){
    point_effect[k] ~ dnorm(0, tau.point)
  }

  ## PRIORS
  ## intercept, slope, and variance-covariance matrix:
  mu.intercept ~ dnorm(0, 1/10000)
  mu.slope ~ dnorm(0, 1/10000)

  sigma.intercept ~ dunif(0, 10)
  sigma.slope ~ dunif(0, 10)

  # correlation coefficient between alpha and beta per transect
  # random effects:
  # - correlation coefficient between transect-level slope and intercept
  rho ~ dunif(-1, 1)

  # - variance-covariance matrix:
  #     * first the diagonals = variance
  vcov[1,1] <- sigma.intercept^2
  vcov[2,2] <- sigma.slope^2

  #     * off-diagonals = covariance (rho x product of standard deviations)
  vcov[1,2] <- rho * sigma.intercept * sigma.slope
  vcov[2,1] <- vcov[1,2]

  # - inverse of covariance matrix
  Tau.B[1:2, 1:2] <- inverse(vcov[1:2, 1:2])

  ## Additional priors:
  sigma.year ~ dunif(0, 10)
  tau.year = 1/sigma.year^2

  sigma.point ~ dunif(0, 10)
  tau.point = 1/sigma.point^2

  sigma.overdispersion ~ dunif(0, 10)
  tau.overdispersion = 1/sigma.overdispersion^2

  ## ANNUAL ABUNDANCE INDICES:
  # - overall:
  B.hat_pred[1, 1] <- mu.intercept
  B.hat_pred[1, 2] <- mu.slope
  B_pred[1, 1:2] ~ dmnorm(B.hat_pred[1, 1:2], Tau.B)

  for (t in 1:nyears) {
    log.abund[t] <- B_pred[1, 1] + B_pred[1, 2] * (t-1) + year_effect[t] +
      0.5 * sigma.overdispersion^2
    abund[t] <- prop[t] * exp(log.abund[t])
  }

  # - by transect:
  for (t in 1:nyears) {
    for (j in 1:ntransects) {
      log(abund_transect[t, j]) <- intercept[j] + slope[j] * (t-1) +
        year_effect[t] + 0.5 * sigma.overdispersion^2
    }
  }

  ## PREDICTED VALUES:
  # - overall in each year
  for (i in 1:length(zyear.pred)) {
    log(count_pred[i]) <- B_pred[1, 1] + B_pred[1, 2] * zyear.pred[i]
  }


  # - by transect
  for (i in 1:length(zyear.pred)) {
    for (j in 1:ntransects) {
      log(count_pred_transect[i, j]) <- intercept[j] + slope[j] * zyear.pred[i]
    }
  }
}"} else {
  vars = c('mu.intercept', 'mu.slope', 'sigma.intercept', 'sigma.slope',
           'sigma.year', 'sigma.point', 'slope', 'intercept',
           'abund', 'abund_transect', 'count_pred', 'count_pred_transect')

  modelstring = "
  model{
  ## ecological model
  for (i in 1:length(observed)){
    log(lambda[i]) = intercept[transect[i]] + slope[transect[i]] * zyear[i] +
      year_effect[year[i]] + point_effect[point[i]]

    observed[i] ~ dpois(lambda[i])
  }

  ## transect intercept and slope
  for (j in 1:ntransects) {
    intercept[j] <- B[j, 1]  # group level intercept
    slope[j]  <- B[j, 2]  # group level slope
    B[j, 1:2] ~ dmnorm(B.hat[j, 1:2], Tau.B)
    B.hat[j, 1] <- mu.intercept  # required by JAGS syntax
    B.hat[j, 2] <- mu.slope   # required by JAGS syntax
  }

  ## additional random intercepts for:
  ## - year
  for (t in 1:nyears){
    year_effect[t] ~ dnorm(0, tau.year)
  }
  ## - point (not nested within transect)
  for (k in 1:npoints){
    point_effect[k] ~ dnorm(0, tau.point)
  }

  ## PRIORS
  ## intercept, slope, and variance-covariance matrix:
  mu.intercept ~ dnorm(0, 1/10000)
  mu.slope ~ dnorm(0, 1/10000)

  sigma.intercept ~ dunif(0, 10)
  sigma.slope ~ dunif(0, 10)

  # correlation coefficient between alpha and beta per transect
  # random effects:
  # - correlation coefficient between transect-level slope and intercept
  rho ~ dunif(-1, 1)

  # - variance-covariance matrix:
  #     * first the diagonals = variance
  vcov[1,1] <- sigma.intercept^2
  vcov[2,2] <- sigma.slope^2

  #     * off-diagonals = covariance (rho x product of standard deviations)
  vcov[1,2] <- rho * sigma.intercept * sigma.slope
  vcov[2,1] <- vcov[1,2]

  # - inverse of covariance matrix
  Tau.B[1:2, 1:2] <- inverse(vcov[1:2, 1:2])

  ## Additional priors:
  sigma.year ~ dunif(0, 10)
  tau.year = 1/sigma.year^2

  sigma.point ~ dunif(0, 10)
  tau.point = 1/sigma.point^2

  ## ANNUAL ABUNDANCE INDICES:
  # - overall:
  B.hat_pred[1, 1] <- mu.intercept
  B.hat_pred[1, 2] <- mu.slope
  B_pred[1, 1:2] ~ dmnorm(B.hat_pred[1, 1:2], Tau.B)

  for (t in 1:nyears) {
    log.abund[t] <- B_pred[1, 1] + B_pred[1, 2] * (t-1) + year_effect[t]
    abund[t] <- prop[t] * exp(log.abund[t])
  }

  # - by transect:
  for (t in 1:nyears) {
    for (j in 1:ntransects) {
      log(abund_transect[t, j]) <- intercept[j] + slope[j] * (t-1) +
        year_effect[t]
    }
  }

  ## PREDICTED VALUES:
  # - overall in each year
  for (i in 1:length(zyear.pred)) {
    log(count_pred[i]) <- B_pred[1, 1] + B_pred[1, 2] * zyear.pred[i]
  }


  # - by transect
  for (i in 1:length(zyear.pred)) {
    for (j in 1:ntransects) {
      log(count_pred_transect[i, j]) <- intercept[j] + slope[j] * zyear.pred[i]
    }
  }
}"
}
  } else if (randomslope == FALSE) {
    if (overdispersion == TRUE) {
      vars = c('intercept', 'slope', 'sigma.year', 'sigma.transect',
               'sigma.point', 'sigma.overdispersion',
               'abund', 'abund_transect', 'count_pred', 'count_pred_transect')

      modelstring = "
  model {
    ## ecological model
    for (i in 1:length(observed)){
      log(lambda[i]) = intercept + slope * zyear[i] +
        year_effect[year[i]] + transect_effect[transect[i]] +
        point_effect[point[i]] + overdispersion[i]

      observed[i] ~ dpois(lambda[i])

      # overdisperson by observation
      overdispersion[i] ~ dnorm(0, tau.overdispersion)
    }

    ## random year effect
    for (t in 1:nyears){
      year_effect[t] ~ dnorm(0, tau.year)
    }

    ## random transect effect
    for (j in 1:ntransects){
      transect_effect[j] ~ dnorm(0, tau.transect)
    }

    for (k in 1:npoints){
      point_effect[k] ~ dnorm(0, tau.point)
    }

    ## PRIORS
    intercept ~ dnorm(0, 1/10000)
    slope ~ dnorm(0, 1/10000)

    sigma.year ~ dunif(0,10)
    tau.year = 1/sigma.year^2

    sigma.transect ~ dunif(0, 10)
    tau.transect = 1/sigma.transect^2

    sigma.point ~ dunif(0, 10)
    tau.point = 1/sigma.point^2

    sigma.overdispersion ~ dunif(0, 10)
    tau.overdispersion = 1/sigma.overdispersion^2

    ## ANNUAL ABUNDANCE INDICES:
    # - overall:
    for (t in 1:nyears) {
      log.abund[t] <- intercept + year_effect[t] + slope * (t-1) +
        0.5 * sigma.transect^2 + 0.5 * sigma.overdispersion^2
      abund[t] <- prop[t] * exp(log.abund[t])
    }

    # - by transect:
    for (t in 1:nyears) {
      for (j in 1:ntransects) {
        log(abund_transect[t, j]) <- intercept + year_effect[t] + slope * (t-1) +
          transect_effect[j] + 0.5 * sigma.overdispersion^2
      }
    }

    ## PREDICTED VALUES:
    # - overall in each year
    for (i in 1:length(zyear.pred)) {
      log(count_pred[i]) <- intercept + slope * zyear.pred[i]
    }

    # - by transect
    for (i in 1:length(zyear.pred)) {
      for (j in 1:ntransects) {
        log(count_pred_transect[i, j]) <- intercept + slope * zyear.pred[i] +
          transect_effect[j]
      }
    }
  }"} else {
    vars = c('intercept', 'slope', 'sigma.year', 'sigma.transect', 'sigma.point',
             'abund', 'abund_transect', 'count_pred', 'count_pred_transect')

    modelstring = "
  model {
    ## ecological model
    for (i in 1:length(observed)){
      log(lambda[i]) = intercept + slope * zyear[i] +
        year_effect[year[i]] + transect_effect[transect[i]] +
        point_effect[point[i]]

      observed[i] ~ dpois(lambda[i])
    }

    ## random year effect
    for (t in 1:nyears){
      year_effect[t] ~ dnorm(0, tau.year)
    }

    ## random transect effect
    for (j in 1:ntransects){
      transect_effect[j] ~ dnorm(0, tau.transect)
    }

    for (k in 1:npoints){
      point_effect[k] ~ dnorm(0, tau.point)
    }

    ## PRIORS
    intercept ~ dnorm(0, 1/10000)
    slope ~ dnorm(0, 1/10000)

    sigma.year ~ dunif(0, 10)
    tau.year = 1/sigma.year^2

    sigma.transect ~ dunif(0, 10)
    tau.transect = 1/sigma.transect^2

    sigma.point ~ dunif(0, 10)
    tau.point = 1/sigma.point^2

    ## ANNUAL ABUNDANCE INDICES:
    # - overall:
    for (t in 1:nyears) {
      log.abund[t] <- intercept + year_effect[t] + slope * (t-1) +
        0.5 * sigma.transect^2
      abund[t] <- prop[t] * exp(log.abund[t])
    }

    # - by transect:
    for (t in 1:nyears) {
      for (j in 1:ntransects) {
        log(abund_transect[t, j]) <- intercept + year_effect[t] + slope * (t-1) +
          transect_effect[j]
      }
    }

    ## PREDICTED VALUES:
    # - overall in each year
    for (i in 1:length(zyear.pred)) {
      log(count_pred[i]) <- intercept + slope * zyear.pred[i]
    }

    # - by transect
    for (i in 1:length(zyear.pred)) {
      for (j in 1:ntransects) {
        log(count_pred_transect[i, j]) <- intercept + slope * zyear.pred[i] +
          transect_effect[j]
      }
    }
  }"
  }
  }

  jm = rjags::jags.model(file = textConnection(modelstring),
                         data = inputdata[-which(names(inputdata) %in% c('year.pred', 'dat'))],
                         n.adapt = n.adapt, n.chains = n.chains, ...)
  rjags::update.jags(jm, n.iter = n.update)
  results = rjags::coda.samples(jm, variable.names = vars, n.iter = n.iter)

  MCMCvis::MCMCsummary(results,
                       params = vars[-which(vars %in% c('abund',
                                                        'abund_transect',
                                                        'count_pred',
                                                        'count_pred_transect'))]) %>%
    print()

  MCMCvis::MCMCtrace(results,
                     params = vars[-which(vars %in% c('abund',
                                                      'abund_transect',
                                                      'count_pred',
                                                      'count_pred_transect'))],
                     pdf = FALSE)

  return(results)
}
