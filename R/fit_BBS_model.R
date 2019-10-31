#' Fit BBS-style hierarchical model for estimating trends in point count survey data
#'
#' @param inputdata Inputdata created by running \code{\link{setup_BBS_model}}
#' @param n.adapt Number of iterations for adaptation, defaults to 500
#' @param n.update Number of iterations for burn-in, defaults to 500
#' @param n.iter Number of iterations for sampling, defauts to 1000
#' @param n.chains Number of MCMC chains, defaults to 3
#' @param overdispersion Defaults to TRUE, to include a term in the model for
#' fitting overdispersion, as in BBS models
#' @param ... Additional arguments passed to \code{\link[rjags]{jags.model}}, e.g. initial values.
#'
#' @return Returns a coda.samples object from rjags. Also displays summary
#' statistics from MCMCsummary and prints trace plots from MCMCplot.
#'
#' @export
#' @import rjags
#' @importFrom stats update
#' @importFrom MCMCvis MCMCsummary MCMCtrace
#'

fit_BBS_model <- function(inputdata, n.adapt = 500, n.update = 500,
                          n.iter = 1000, n.chains = 3,
                          overdispersion = TRUE,
                          ...) {

  if (overdispersion == TRUE) {
      vars = c('intercept', 'slope', 'sigma.year', 'sigma.transect',
               'sigma.point', 'sigma.overdispersion', 'index', 'trend')

      modelstring = "
  model {
    ## ecological model
    for (i in 1:length(observed)){
      log(lambda[i]) = intercept[project[i]] + slope[project[i]] * zyear[i] +
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
    for (p in 1:nprojects) {
      intercept[p] ~ dnorm(0, 1/10000)
      slope[p] ~ dnorm(0, 1/10000)
    }

    sigma.year ~ dunif(0,10)
    tau.year = 1/sigma.year^2

    sigma.transect ~ dunif(0, 10)
    tau.transect = 1/sigma.transect^2

    sigma.point ~ dunif(0, 10)
    tau.point = 1/sigma.point^2

    sigma.overdispersion ~ dunif(0, 10)
    tau.overdispersion = 1/sigma.overdispersion^2

    ## ANNUAL ABUNDANCE INDICES:
    for (p in 1:nprojects) {
      for (t in 1:nyears) {
        log.index[t, p] <- intercept[p] + slope[p] * (t-1) + year_effect[t] +
          0.5 * sigma.transect^2 + 0.5 * sigma.overdispersion^2
        index[t, p] <- prop[t, p] * exp(log.index[t, p])
      }
    }

    ## PREDICTED VALUES:
    for (p in 1:nprojects) {
      for (i in 1:length(zyear.pred)) {
        log(trend[i, p]) <- intercept[p] + slope[p] * zyear.pred[i]
      }
    }
  }"
  } else if (overdispersion == FALSE) {
    vars = c('intercept', 'slope', 'sigma.year', 'sigma.transect',
             'sigma.point', 'index', 'trend')

    modelstring = "
  model {
    ## ecological model
    for (i in 1:length(observed)){
      log(lambda[i]) = intercept[project[i]] + slope[project[i]] * zyear[i] +
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
    for (p in 1:nprojects) {
      intercept[p] ~ dnorm(0, 1/10000)
      slope[p] ~ dnorm(0, 1/10000)
    }

    sigma.year ~ dunif(0, 10)
    tau.year = 1/sigma.year^2

    sigma.transect ~ dunif(0, 10)
    tau.transect = 1/sigma.transect^2

    sigma.point ~ dunif(0, 10)
    tau.point = 1/sigma.point^2

  ## ANNUAL ABUNDANCE INDICES:
  for (p in 1:nprojects) {
    for (t in 1:nyears) {
      log.index[t, p] <- intercept[p] + slope[p] * (t-1) + year_effect[t] +
        0.5 * sigma.transect^2
      index[t, p] <- prop[t, p] * exp(log.index[t, p])
    }
  }

  ## PREDICTED VALUES:
  for (p in 1:nprojects) {
    for (i in 1:length(zyear.pred)) {
      log(trend[i, p]) <- intercept[p] + slope[p] * zyear.pred[i]
    }
  }
  }"
  }

  jm = jags.model(file = textConnection(modelstring),
                         data = inputdata[-which(names(inputdata) %in% c('year.pred', 'dat'))],
                         n.adapt = n.adapt, n.chains = n.chains, ...)
  update(jm, n.iter = n.update)
  results = rjags::coda.samples(jm, variable.names = vars, n.iter = n.iter)

  MCMCvis::MCMCsummary(results,
                       params = vars[-which(vars %in% c('index', 'trend'))]) %>%
    print()

  MCMCvis::MCMCtrace(results,
                     params = vars[-which(vars %in% c('index', 'trend'))],
                     pdf = FALSE)

  return(results)
}
