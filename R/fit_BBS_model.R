#' Fit BBS-style hierarchical model for estimating trends in point count survey data
#' 
#' Uses \pkg{rjags} to fit a hierarchical model to point count data to estimate trends in a single species' abundance over time, with support for simultaneous estimation of separate overall trends for one or more "project" areas (e.g. POGO vs. PINN), with each project including many point count stations grouped into transects. Models include random intercepts for year within project, transects, and points. Optionally, models can also include a term for fitting overdispersion. Install \pkg{rjags} with: \code{install.packages("rjags")}
#' 
#' If plot = TRUE or print = TRUE, uses \pkg{MCMCvis} to plot traceplots and print a results summary table. Install it with: \code{install.packages("MCMCvis")}
#'
#' @param inputdata Inputdata created by running \code{\link{setup_BBS_model}}
#' @param n.adapt Number of iterations for adaptation, defaults to 1000
#' @param n.burnin Number of iterations for burn-in, defaults to 5000
#' @param n.sample Number of iterations for sampling, defauts to 10000
#' @param n.chains Number of MCMC chains, defaults to 3
#' @param overdispersion Defaults to TRUE, to include a term in the model for
#' fitting overdispersion, as in BBS models
#' @param print if TRUE, print table summarizing results
#' @param plot if TRUE, plot traceplots
#' @param ... Additional arguments passed to \code{\link[rjags]{jags.model}},
#' e.g. initial values.
#'
#' @return Returns a coda.samples object from rjags. Also displays summary
#' statistics from MCMCsummary and prints trace plots from MCMCplot.
#'
#' @export

fit_BBS_model <- function(inputdata, n.adapt = 1000, n.burnin = 5000,
                          n.sample = 10000, n.chains = 3,
                          overdispersion = TRUE, print = TRUE, plot = TRUE,
                          ...) {

  if (!requireNamespace("rjags", quietly = TRUE)) {
    stop(
      "Package 'rjags' is required for fit_BBS_model(). ",
      "Install it with install.packages('rjags').",
      "You may also need to install JAGS from: https://mcmc-jags.sourceforge.io/",
      call. = FALSE
    )
  }
  
  if (plot | print) {
    if (!requireNamespace("MCMCvis", quietly = TRUE)) {
      stop(
        "Package 'MCMCvis' is required for fit_BBS_model() with plot = TRUE or print = TRUE. ",
        "Install it with install.packages('MCMCvis').",
        call. = FALSE
      )
    }
  }
  
  
  if (overdispersion == TRUE) {
      vars = c('intercept', 'slope', 'sigma.year', 'sigma.transect',
               'sigma.point', 'sigma.overdispersion', 'index', 'trend')

      modelstring = "
  model {
    ## ecological model
    for (i in 1:length(observed)){
      log(lambda[i]) = intercept[project[i]] + slope[project[i]] * zyear[i] +
        year_effect[year[i], project[i]] + transect_effect[transect[i]] +
        point_effect[point[i]] + overdispersion[i]

      observed[i] ~ dpois(lambda[i])

      # overdisperson by observation
      overdispersion[i] ~ dnorm(0, tau.overdispersion)
    }

    ## random year effect - variance by project
    for (p in 1:nprojects) {
      for (t in 1:nyears){
        year_effect[t, p] ~ dnorm(0, tau.year[p])
      }
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

      sigma.year[p] ~ dunif(0,10)
      tau.year[p] = 1/sigma.year[p]^2
    }

    sigma.transect ~ dunif(0, 10)
    tau.transect = 1/sigma.transect^2

    sigma.point ~ dunif(0, 10)
    tau.point = 1/sigma.point^2

    sigma.overdispersion ~ dunif(0, 10)
    tau.overdispersion = 1/sigma.overdispersion^2

    ## ANNUAL ABUNDANCE INDICES:
    for (p in 1:nprojects) {
      for (t in 1:nyears) {
        log.index[t, p] <- intercept[p] + slope[p] * (t-1) + year_effect[t, p]
        # + 0.5 * sigma.transect^2 + 0.5 * sigma.overdispersion^2
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
        year_effect[year[i], project[i]] + transect_effect[transect[i]] +
        point_effect[point[i]]

      observed[i] ~ dpois(lambda[i])
    }

    ## random year effect
    for (p in 1:nprojects) {
      for (t in 1:nyears){
        year_effect[t, p] ~ dnorm(0, tau.year[p])
      }
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

      sigma.year[p] ~ dunif(0, 10)
      tau.year[p] = 1/sigma.year[p]^2
    }

    sigma.transect ~ dunif(0, 10)
    tau.transect = 1/sigma.transect^2

    sigma.point ~ dunif(0, 10)
    tau.point = 1/sigma.point^2

  ## ANNUAL ABUNDANCE INDICES:
  for (p in 1:nprojects) {
    for (t in 1:nyears) {
      log.index[t, p] <- intercept[p] + slope[p] * (t-1) + year_effect[t, p]
      # + 0.5 * sigma.transect^2
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

  rjags::load.module('glm')

  jm = rjags::jags.model(file = textConnection(modelstring),
                         data = inputdata[-which(names(inputdata) %in% c('year.pred', 'dat'))],
                         n.adapt = n.adapt,
                         n.chains = n.chains,
                         ...)
  
  stats::update(jm, n.iter = n.burnin)
  
  results = rjags::coda.samples(jm, variable.names = vars, n.iter = n.sample,
                                thin = 1)

  if (print) {
    MCMCvis::MCMCsummary(
      results,
      params = vars[-which(vars %in% c('index', 'trend'))]) |> 
      print()
  }
  
  if (plot) {
    MCMCvis::MCMCtrace(
      results,
      params = vars[-which(vars %in% c('index', 'trend'))], pdf = FALSE)
  }
  

  return(results)
}
